#! /usr/bin/env python

import sys
import os
import time
import glob
import multiprocessing
from operator import itemgetter
from subprocess import call

from epac.ete2 import Tree, SeqGroup
from epac.argparse import ArgumentParser,RawDescriptionHelpFormatter
from epac.config import SativaConfig,EpacConfig
from epac.raxml_util import RaxmlWrapper, FileUtils
from epac.json_util import RefJsonParser, RefJsonChecker, EpaJsonParser
from epac.taxonomy_util import TaxCode, Taxonomy
from epac.classify_util import TaxTreeHelper,TaxClassifyHelper
import epa_trainer

class LeaveOneTest:
    def __init__(self, config):
        self.cfg = config
        
        self.mis_fname = self.cfg.out_fname("%NAME%.mis")
        self.premis_fname = self.cfg.out_fname("%NAME%.premis")
        self.misrank_fname = self.cfg.out_fname("%NAME%.misrank")
        self.stats_fname = self.cfg.out_fname("%NAME%.stats")
        
        if os.path.isfile(self.mis_fname):
            print "\nERROR: Output file already exists: %s" % self.mis_fname
            print "Please specify a different job name using -n or remove old output files."
            self.cfg.exit_user_error()

        self.tmp_refaln = config.tmp_fname("%NAME%.refaln")
        self.reftree_lbl_fname = config.tmp_fname("%NAME%_lbl.tre")
        self.reftree_tax_fname = config.tmp_fname("%NAME%_tax.tre")
        self.optmod_fname = self.cfg.tmp_fname("%NAME%.opt")
        self.reftree_fname = self.cfg.tmp_fname("ref_%NAME%.tre")

        self.mislabels = []
        self.mislabels_cnt = []
        self.rank_mislabels = []
        self.rank_mislabels_cnt = []
        self.misrank_conf_map = {}
        
    def write_bid_tax_map(self, bid_tax_map, final):
        if self.cfg.debug:
            fname_suffix = "final" if final else "l1out"
            bid_fname = self.cfg.tmp_fname("%NAME%_" + "bid_tax_map_%s.txt" % fname_suffix)
            with open(bid_fname, "w") as outf:
              for bid, bid_rec in bid_tax_map.iteritems():
                outf.write("%s\t%s\t%d\t%f\n" % (bid, bid_rec[0], bid_rec[1], bid_rec[2]));    

    def write_assignments(self, assign_map, final):
        if self.cfg.debug:
            fname_suffix = "final" if final else "l1out"
            assign_fname = self.cfg.tmp_fname("%NAME%_" + "taxassign_%s.txt" % fname_suffix)
            with open(assign_fname, "w") as outf:
                for seq_name in assign_map.iterkeys():
                    ranks, lws = assign_map[seq_name]
                    outf.write("%s\t%s\t%s\n" % (seq_name, ";".join(ranks), ";".join(["%.3f" % l for l in lws])))

    def load_refjson(self, refjson_fname):
        try:
            self.refjson = RefJsonParser(refjson_fname)
        except ValueError:
            self.cfg.exit_user_error("ERROR: Invalid json file format!")
            
        #validate input json format 
        (valid, err) = self.refjson.validate()
        if not valid:
            self.cfg.log.error("ERROR: Parsing reference JSON file failed:\n%s", err)
            self.cfg.exit_user_error()
        
        self.rate = self.refjson.get_rate()
        self.node_height = self.refjson.get_node_height()
        self.origin_taxonomy = self.refjson.get_origin_taxonomy()
        self.tax_tree = self.refjson.get_tax_tree()
        self.cfg.compress_patterns = self.refjson.get_pattern_compression()

        self.bid_taxonomy_map = self.refjson.get_branch_tax_map()
        if not self.bid_taxonomy_map:
            # old file format (before 1.6), need to rebuild this map from scratch
            th = TaxTreeHelper(self.cfg, self.origin_taxonomy)
            th.set_mf_rooted_tree(self.tax_tree)
            th.set_bf_unrooted_tree(self.refjson.get_reftree())
            self.bid_taxonomy_map = th.get_bid_taxonomy_map()
            
        self.write_bid_tax_map(self.bid_taxonomy_map, final=False)

        reftree_str = self.refjson.get_raxml_readable_tree()
        self.reftree = Tree(reftree_str)
        self.reftree_size = len(self.reftree.get_leaves())

        # IMPORTANT: set EPA heuristic rate based on tree size!                
        self.cfg.resolve_auto_settings(self.reftree_size)
        # If we're loading the pre-optimized model, we MUST set the same rate het. mode as in the ref file        
        if self.cfg.epa_load_optmod:
            self.cfg.raxml_model = self.refjson.get_ratehet_model()

        self.classify_helper = TaxClassifyHelper(self.cfg, self.bid_taxonomy_map, self.rate, self.node_height)
        self.taxtree_helper = TaxTreeHelper(self.cfg, self.origin_taxonomy, self.tax_tree)
        
        tax_code_name = self.refjson.get_taxcode()
        self.tax_code = TaxCode(tax_code_name)
        
        self.taxonomy = Taxonomy(prefix=EpacConfig.REF_SEQ_PREFIX, tax_map=self.origin_taxonomy)
        self.tax_common_ranks = self.taxonomy.get_common_ranks()
#        print "Common ranks: ", self.tax_common_ranks

        self.mislabels_cnt = [0] * TaxCode.UNI_TAX_LEVELS
        self.rank_mislabels_cnt = [0] * TaxCode.UNI_TAX_LEVELS
        
    def run_epa_trainer(self):
        epa_trainer.run_trainer(self.cfg)

        if not os.path.isfile(self.cfg.refjson_fname):
            self.cfg.log.error("\nBuilding reference tree failed, see error messages above.")
            self.cfg.exit_fatal_error()
        
    def classify_seq(self, placement):
        edges = placement["p"]
        if len(edges) > 0:
            return self.classify_helper.classify_seq(edges)
        else:
            print "ERROR: no placements! something is definitely wrong!"

    def check_seq_tax_labels(self, seq_name, orig_ranks, ranks, lws):
        mis_rec = None
        
        num_common_ranks = len(self.tax_common_ranks)
        orig_rank_level = Taxonomy.lowest_assigned_rank_level(orig_ranks)
        new_rank_level = Taxonomy.lowest_assigned_rank_level(ranks)
        #if new_rank_level < 0 or (new_rank_level < num_common_ranks and orig_rank_level >= num_common_ranks):
#        if new_rank_level < 0:
        if len(ranks) == 0:
            mis_rec = {}
            mis_rec['name'] = seq_name
            mis_rec['orig_level'] = -1
            mis_rec['real_level'] = 0
            mis_rec['level_name'] = "[NotIngroup]"
            mis_rec['inv_level'] = -1 * mis_rec['real_level']  # just for sorting
            mis_rec['orig_ranks'] = orig_ranks
            mis_rec['ranks'] = []
            mis_rec['lws'] = [1.0]
            mis_rec['conf'] = mis_rec['lws'][0]
        else:
            mislabel_lvl = -1
            min_len = min(len(orig_ranks),len(ranks))
            for rank_lvl in range(min_len):
                if ranks[rank_lvl] != Taxonomy.EMPTY_RANK and ranks[rank_lvl] != orig_ranks[rank_lvl]:
                    mislabel_lvl = rank_lvl
                    break

            if mislabel_lvl >= 0:
                real_lvl = self.tax_code.guess_rank_level(orig_ranks, mislabel_lvl)
                mis_rec = {}
                mis_rec['name'] = seq_name
                mis_rec['orig_level'] = mislabel_lvl
                mis_rec['real_level'] = real_lvl
                mis_rec['level_name'] = self.tax_code.rank_level_name(real_lvl)[0]
                mis_rec['inv_level'] = -1 * mis_rec['real_level']  # just for sorting
                mis_rec['orig_ranks'] = orig_ranks
                mis_rec['ranks'] = ranks
                mis_rec['lws'] = lws
                mis_rec['conf'] = lws[mislabel_lvl]
    
        if mis_rec:
            self.mislabels.append(mis_rec)
            
        return mis_rec
        
    def filter_mislabels(self):
        filtered_mis = []
        for i in range(len(self.mislabels)):
            if self.mislabels[i]['conf'] >= self.cfg.conf_cutoff:
                filtered_mis.append(self.mislabels[i])
        
        self.mislabels = filtered_mis

    def check_rank_tax_labels(self, rank_name, orig_ranks, ranks, lws):
        mislabel_lvl = -1
        min_len = min(len(orig_ranks),len(ranks))
        for rank_lvl in range(min_len):
            if ranks[rank_lvl] != Taxonomy.EMPTY_RANK and ranks[rank_lvl] != orig_ranks[rank_lvl]:
                mislabel_lvl = rank_lvl
                break

        if mislabel_lvl >= 0:
            real_lvl = self.tax_code.guess_rank_level(orig_ranks, mislabel_lvl)
            mis_rec = {}
            mis_rec['name'] = rank_name
            mis_rec['orig_level'] = mislabel_lvl
            mis_rec['real_level'] = real_lvl
            mis_rec['level_name'] = self.tax_code.rank_level_name(real_lvl)[0]
            mis_rec['inv_level'] = -1 * real_lvl  # just for sorting
            mis_rec['orig_ranks'] = orig_ranks
            mis_rec['ranks'] = ranks
            mis_rec['lws'] = lws
            mis_rec['conf'] = lws[mislabel_lvl]
            self.rank_mislabels.append(mis_rec)
               
            return mis_rec
        else:
            return None                

    def mis_rec_to_string_old(self, mis_rec):
        lvl = mis_rec['orig_level']
        output = mis_rec['name'] + "\t"
        output += "%s\t%s\t%s\t%.3f\n" % (mis_rec['level_name'], 
            mis_rec['orig_ranks'][lvl], mis_rec['ranks'][lvl], mis_rec['lws'][lvl])
        output += ";".join(mis_rec['orig_ranks']) + "\n"
        output += ";".join(mis_rec['ranks']) + "\n"
        output += "\t".join(["%.3f" % conf for conf in mis_rec['lws']]) + "\n"
        return output

    def mis_rec_to_string(self, mis_rec):
        lvl = mis_rec['orig_level']
        uncorr_name = EpacConfig.strip_ref_prefix(self.refjson.get_uncorr_seqid(mis_rec['name']))
        uncorr_orig_ranks = self.refjson.get_uncorr_ranks(mis_rec['orig_ranks'])
        uncorr_ranks = self.refjson.get_uncorr_ranks(mis_rec['ranks'])
        output = uncorr_name + "\t"
      
        if lvl >= 0:
            output += "%s\t%s\t%s\t%.3f\t" % (mis_rec['level_name'], 
                uncorr_orig_ranks[lvl], uncorr_ranks[lvl], mis_rec['lws'][lvl])
        else:
            output += "%s\t%s\t%s\t%.3f\t" % (mis_rec['level_name'], 
                "NA", "NA", mis_rec['lws'][0])
        
        output += Taxonomy.lineage_str(uncorr_orig_ranks) + "\t"
        output += Taxonomy.lineage_str(uncorr_ranks) + "\t"
        output += ";".join(["%.3f" % conf for conf in mis_rec['lws']])
        if 'rank_conf' in mis_rec:
            output += "\t%.3f" % mis_rec['rank_conf']
        return output

    def sort_mislabels(self):
        self.mislabels = sorted(self.mislabels, key=itemgetter('inv_level', 'conf', 'name'), reverse=True)
        for mis_rec in self.mislabels:
            real_lvl = mis_rec["real_level"]
            self.mislabels_cnt[real_lvl] += 1
        
        if self.cfg.ranktest:
            self.rank_mislabels = sorted(self.rank_mislabels, key=itemgetter('inv_level', 'conf', 'name'), reverse=True)
            for mis_rec in self.rank_mislabels:
                real_lvl = mis_rec["real_level"]
                self.rank_mislabels_cnt[real_lvl] += 1
    
    def write_stats(self, toFile=False):
        self.cfg.log.info("Mislabeled sequences by rank:")
        seq_sum = 0
        rank_sum = 0
        stats = []
        for i in range(len(self.mislabels_cnt)):
            if i > 0:
                rname = self.tax_code.rank_level_name(i)[0].ljust(12)
            else:
                rname = "[NotIngroup]"
            if self.mislabels_cnt[i] > 0:
                seq_sum += self.mislabels_cnt[i]
#                    output = "%s:\t%d" % (rname, seq_sum)
                output = "%s:\t%d" % (rname, self.mislabels_cnt[i])
                if self.cfg.ranktest:
                    rank_sum += self.rank_mislabels_cnt[i]
                    output += "\t%d" % rank_sum
                self.cfg.log.info(output) 
                stats.append(output)

        if toFile:
            with open(self.stats_fname, "w") as fo_stat:
                for line in stats:
                    fo_stat.write(line + "\n")
    
    def write_mislabels(self, final=True):
        if final:
            out_fname = self.mis_fname
        else:
            out_fname = self.premis_fname
        
        with open(out_fname, "w") as fo_all:
            fields = ["SeqID", "MislabeledLevel", "OriginalLabel", "ProposedLabel", "Confidence", "OriginalTaxonomyPath", "ProposedTaxonomyPath", "PerRankConfidence"]
            if self.cfg.ranktest:
                fields += ["HigherRankMisplacedConfidence"]
            header = ";" + "\t".join(fields) + "\n"
            fo_all.write(header)
            if self.cfg.verbose and len(self.mislabels) > 0 and final:
                print "Mislabeled sequences:\n"
                print header 
            for mis_rec in self.mislabels:
                output = self.mis_rec_to_string(mis_rec)  + "\n"
                fo_all.write(output)
                if self.cfg.verbose and final:
                    print(output) 
                    
        if not final:
            return

        if self.cfg.ranktest:
            with open(self.misrank_fname, "w") as fo_all:
                fields = ["RankID", "MislabeledLevel", "OriginalLabel", "ProposedLabel", "Confidence", "OriginalTaxonomyPath", "ProposedTaxonomyPath", "PerRankConfidence"]
                header = ";" + "\t".join(fields)  + "\n"
                fo_all.write(header)
                if self.cfg.verbose  and len(self.rank_mislabels) > 0:
                    print "\nMislabeled higher ranks:\n"
                    print header 
                for mis_rec in self.rank_mislabels:
                    output = self.mis_rec_to_string(mis_rec) + "\n"
                    fo_all.write(output)
                    if self.cfg.verbose:
                        print(output) 
                        
        self.write_stats()
   
    def run_leave_subtree_out_test(self):
        job_name = self.cfg.subst_name("l1out_rank_%NAME%")
#        if self.jplace_fname:
#            jp = EpaJsonParser(self.jplace_fname)
#        else:        

        #create file with subtrees
        rank_tips = {}
        rank_parent = {}
        for node in self.tax_tree.traverse("postorder"):
            if node.is_leaf() or node.is_root():
                continue
            tax_path = node.name
            ranks = Taxonomy.split_rank_uid(tax_path)
            rank_lvl = Taxonomy.lowest_assigned_rank_level(ranks)
            if rank_lvl < 2:
                continue
                
            parent_ranks = Taxonomy.split_rank_uid(node.up.name)
            parent_lvl = Taxonomy.lowest_assigned_rank_level(parent_ranks)
            if parent_lvl < 1:
                continue
            
            rank_seqs = node.get_leaf_names()
            rank_size = len(rank_seqs)
            if rank_size < 2 or rank_size > self.reftree_size-4:
                continue

#            print rank_lvl, "\t", tax_path, "\t", rank_seqs, "\n"
            rank_tips[tax_path] = node.get_leaf_names()
            rank_parent[tax_path] = parent_ranks
                
        subtree_list = rank_tips.items()
        
        if len(subtree_list) == 0:
            return 0
            
        subtree_list_file = self.cfg.tmp_fname("treelist_%NAME%.txt")
        with open(subtree_list_file, "w") as fout:
            for rank_name, tips in subtree_list:
                fout.write("%s\n" % " ".join(tips))
        
        jp_list = self.raxml.run_epa(job_name, self.refalign_fname, self.reftree_fname, self.optmod_fname, 
            mode="l1o_subtree", subtree_fname=subtree_list_file)

        subtree_count = 0
        for jp in jp_list:
            placements = jp.get_placement()
            for place in placements:
                ranks, lws = self.classify_seq(place)
                tax_path = subtree_list[subtree_count][0]
                orig_ranks = Taxonomy.split_rank_uid(tax_path)
                rank_level = Taxonomy.lowest_assigned_rank_level(orig_ranks)
                rank_prefix = self.guess_rank_level_name(orig_ranks, rank_level)[0]
                rank_name = orig_ranks[rank_level]
                if not rank_name.startswith(rank_prefix):
                    rank_name = rank_prefix + rank_name
                parent_ranks = rank_parent[tax_path]
#                print orig_ranks, "\n", parent_ranks, "\n", ranks, "\n"
                mis_rec = self.check_rank_tax_labels(rank_name, parent_ranks, ranks, lws)
                if mis_rec:
                    self.misrank_conf_map[tax_path] = mis_rec['conf']
                subtree_count += 1

        return subtree_count    
        
    def run_leave_seq_out_test(self):
        job_name = self.cfg.subst_name("l1out_seq_%NAME%")
        placements = []
        if self.cfg.jplace_fname:
            if os.path.isdir(self.cfg.jplace_fname):
                jplace_fmask = os.path.join(self.cfg.jplace_fname, '*.jplace')
            else:
                jplace_fmask = self.cfg.jplace_fname

            jplace_fname_list = glob.glob(jplace_fmask)
            for jplace_fname in jplace_fname_list:
                jp = EpaJsonParser(jplace_fname)
                placements += jp.get_placement()
                
            config.log.debug("Loaded %d placements from %s\n", len(placements), jplace_fmask)
        else:        
            jp = self.raxml.run_epa(job_name, self.refalign_fname, self.reftree_fname, self.optmod_fname, mode="l1o_seq")
            placements = jp.get_placement()
            if self.cfg.output_interim_files:
                out_jplace_fname = self.cfg.out_fname("%NAME%.l1out_seq.jplace")
                self.raxml.copy_epa_jplace(job_name, out_jplace_fname, move=True, mode="l1o_seq")
        
        seq_count = 0
        l1out_ass = {}
        for place in placements:
            seq_name = place["n"][0]
            
            # get original taxonomic label
#            orig_ranks = self.get_orig_ranks(seq_name)
            orig_ranks =  self.taxtree_helper.get_seq_ranks_from_tree(seq_name)

            # get EPA tax label
            ranks, lws = self.classify_seq(place)
            l1out_ass[seq_name] = (ranks, lws)
            
            # check if they match
            mis_rec = self.check_seq_tax_labels(seq_name, orig_ranks, ranks, lws)
            # cross-check with higher rank mislabels
            if self.cfg.ranktest and mis_rec:
                rank_conf = 0
                for lvl in range(2,len(orig_ranks)):
                    tax_path = Taxonomy.get_rank_uid(orig_ranks, lvl)
                    if tax_path in self.misrank_conf_map:
                        rank_conf = max(rank_conf, self.misrank_conf_map[tax_path])
                mis_rec['rank_conf'] = rank_conf
            seq_count += 1

        self.write_assignments(l1out_ass, final=False)
            
        return seq_count    
        
    def run_final_epa_test(self):
        self.reftree_outgroup = self.refjson.get_outgroup()

        tmp_reftree = self.reftree.copy(method="newick") 
        name2refnode = {}
        for leaf in tmp_reftree.iter_leaves():
            name2refnode[leaf.name] = leaf        

        tmp_taxtree = self.tax_tree.copy(method="newick") 
        name2taxnode = {}
        for leaf in tmp_taxtree.iter_leaves():
            name2taxnode[leaf.name] = leaf        

        for mis_rec in self.mislabels:
            rname = mis_rec['name']
#            rname = EpacConfig.REF_SEQ_PREFIX + name

            if rname in name2refnode:
                name2refnode[rname].delete()
            else:
                print "Node not found in the reference tree: %s" % rname

            if rname in name2taxnode:
                name2taxnode[rname].delete()
            else:
                print "Node not found in the taxonomic tree: %s" % rname

        # remove unifurcation at the root
        if len(tmp_reftree.children) == 1:
            tmp_reftree = tmp_reftree.children[0]
            
        self.mislabels = []

        th = TaxTreeHelper(self.cfg, self.origin_taxonomy)
        th.set_mf_rooted_tree(tmp_taxtree)
            
        epa_result = self.run_epa_once(tmp_reftree)
        
        reftree_epalbl_str = epa_result.get_std_newick_tree()        
        placements = epa_result.get_placement()
        
        # update branchid-taxonomy mapping to account for possible changes in branch numbering
        reftree_tax = Tree(reftree_epalbl_str)
        th.set_bf_unrooted_tree(reftree_tax)
        bid_tax_map = th.get_bid_taxonomy_map()
        
        self.write_bid_tax_map(bid_tax_map, final=True)

        cl = TaxClassifyHelper(self.cfg, bid_tax_map, self.rate, self.node_height)
        
#        newtax_fname = self.cfg.subst_name("newtax_%NAME%.tre")
#        th.get_tax_tree().write(outfile=newtax_fname, format=3)

        final_ass = {}
        for place in placements:
            seq_name = place["n"][0]

            # get original taxonomic label
            orig_ranks = self.taxtree_helper.get_seq_ranks_from_tree(seq_name)

            # EXPERIMENTAL FEATURE - disabled for now!
            # It could happen that certain ranks were present in the "original" reference tree, but 
            # are completely missing in the pruned tree (e.g., all seqs of a species were considered "suspicious" 
            # after the leave-one-out test and thus pruned)
            # In this case, EPA has no chance to infer full original taxonomic annotation (=species) since the corresponding clade
            # is now missing. To account for this fact, we amend the original taxonomic annotation and set ranks missing from  
            # pruned tree to "Undefined".
#            orig_ranks = th.strip_missing_ranks(orig_ranks)
#            print orig_ranks

            # get EPA tax label
            ranks, lws = cl.classify_seq(place["p"])
            final_ass[seq_name] = (ranks, lws)

            #print seq_name, ": ", orig_ranks, "--->", ranks

            # check if they match
            mis_rec = self.check_seq_tax_labels(seq_name, orig_ranks, ranks, lws)

        self.write_assignments(final_ass, final=True)

    def run_epa_once(self, reftree):
        reftree_fname = self.cfg.tmp_fname("final_ref_%NAME%.tre")
        job_name = self.cfg.subst_name("final_epa_%NAME%")

        reftree.write(outfile=reftree_fname)

        # IMPORTANT: don't load the model, since it's invalid for the pruned true !!! 
        optmod_fname=""
        epa_result = self.raxml.run_epa(job_name, self.refalign_fname, reftree_fname, optmod_fname)

        if self.cfg.output_interim_files:
            out_jplace_fname = self.cfg.out_fname("%NAME%.final_epa.jplace")
            self.raxml.copy_epa_jplace(job_name, out_jplace_fname, move=True)

        return epa_result

    def run_test(self):
        self.raxml = RaxmlWrapper(self.cfg)

#        config.log.info("Number of sequences in the reference: %d\n", self.reftree_size)

        self.refjson.get_raxml_readable_tree(self.reftree_fname)
        self.refalign_fname = self.refjson.get_alignment(self.tmp_refaln)        
        self.refjson.get_binary_model(self.optmod_fname)
        
        if self.cfg.ranktest:
            config.log.info("Running the leave-one-rank-out test...\n")
            subtree_count = self.run_leave_subtree_out_test()
            
        config.log.info("Running the leave-one-sequence-out test...\n")
        self.run_leave_seq_out_test()

        if len(self.mislabels) > 0:
            config.log.info("Leave-one-out test identified %d suspicious sequences; running final EPA test to check them...\n", len(self.mislabels))
            if self.cfg.debug:
                self.write_mislabels(final=False)
            self.run_final_epa_test()

        self.filter_mislabels()
        self.sort_mislabels()
        self.write_mislabels()
        config.log.info("\nTotal mislabels: %d / %.2f %%", len(self.mislabels), (float(len(self.mislabels)) / self.reftree_size * 100))

def parse_args():
    parser = ArgumentParser(usage="%(prog)s -s ALIGNMENT -t TAXONOMY -x {BAC,BOT,ZOO,VIR} [options]",
    description=EpacConfig.SATIVA_INFO % "SATIVA",
    epilog="Example: sativa.py -s example/test.phy -t example/test.tax -x BAC",
    formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-s", dest="align_fname",
            help="""Reference alignment file (PHYLIP or FASTA). Sequences must be aligned, 
            their IDs must correspond to those in taxonomy file.""")
    parser.add_argument("-t", dest="taxonomy_fname",
            help="""Reference taxonomy file.""")
    parser.add_argument("-x", dest="taxcode_name", choices=["bac", "bot", "zoo", "vir"], type = str.lower,
            help="""Taxonomic code: BAC(teriological), BOT(anical), ZOO(logical), VIR(ological)""")
    parser.add_argument("-n", dest="output_name", default=None,
            help="""Job name, will be used as a prefix for output file names (default: taxonomy file name without extension)""")
    parser.add_argument("-o", dest="output_dir", default=".",
            help="""Output directory (default: current).""")
    parser.add_argument("-T", dest="num_threads", type=int, default=multiprocessing.cpu_count(),
            help="""Specify the number of CPUs (default: %d)""" % multiprocessing.cpu_count())
    parser.add_argument("-N", dest="rep_num", type=int, default=1, 
            help="""Number of RAxML tree searches (with distinct random seeds) to resolve multifurcation. Default: 1""")
    parser.add_argument("-v", dest="verbose", action="store_true",
            help="""Print additional info messages to the console.""")
    parser.add_argument("-R", dest="restart", action="store_true",
            help="""Resume execution after a premature termination (e.g., due to expired job time limit).
Run name of the previous (terminated) job must be specified via -n option.""")
    parser.add_argument("-c", dest="config_fname", default=None,
            help="Config file name.")
    parser.add_argument("-r", dest="ref_fname",
            help="""Specify the reference alignment and taxonomy in refjson format.""")
    parser.add_argument("-j", dest="jplace_fname", default=None,
            help="""Do not call RAxML EPA, use existing .jplace file as input instead. 
            This could be also a directory with *.jplace files.""")
    parser.add_argument("-p", dest="rand_seed", type=int, default=12345,
            help="""Random seed to be used with RAxML. Default: 12345""")
    parser.add_argument("-C", dest="conf_cutoff", type=float, default=0.,
            help="""Confidence cut-off between 0 and 1. Default: 0\n""")
    parser.add_argument("-P", dest="brlen_pv", type=float, default=0.,
            help="""P-value for branch length Erlang test. Default: 0=off\n""")
    parser.add_argument("-l", dest="min_lhw", type=float, default=0.,
            help="""A value between 0 and 1, the minimal sum of likelihood weight of
                    an assignment to a specific rank. This value represents a confidence 
                    measure of the assignment, assignments below this value will be discarded. 
                    Default: 0 to output all possbile assignments.""")
    parser.add_argument("-m", dest="mfresolv_method", choices=["thorough", "fast", "ultrafast"],
            default="thorough", help="""Method of multifurcation resolution: 
            thorough    use stardard constrainted RAxML tree search (default)
            fast        use RF distance as search convergence criterion (RAxML -D option)
            ultrafast   optimize model+branch lengths only (RAxML -f e option)""")
    parser.add_argument("-debug", dest="debug", action="store_true",
            help="""Debug mode, intermediate files will not be cleaned up.""")
    parser.add_argument("-ranktest", dest="ranktest", action="store_true",
            help="""Test for misplaced higher ranks.""")
    parser.add_argument("-tmpdir", dest="temp_dir", default=None,
            help="""Directory for temporary files.""")

    args = parser.parse_args()
    if len(sys.argv) == 1: 
        parser.print_help()
        sys.exit()
    check_args(args, parser)
    return args


def check_args(args, parser):    
    if args.ref_fname:
        if args.align_fname:
            print("WARNING: -r and -s options are mutually exclusive! Your alignment file will be ignored.\n")
        if args.taxonomy_fname:
            print("WARNING: -r and -t options are mutually exclusive! Your taxonomy file will be ignored.\n")
        if args.taxcode_name:
            print("WARNING: -r and -x options are mutually exclusive! The taxonomic code from reference file will be used.\n")
    elif not args.align_fname or not args.taxonomy_fname or not args.taxcode_name:
        print("ERROR: either reference in JSON format or taxonomy, alignment and taxonomic code name must be provided:\n")
        parser.print_help()
        sys.exit()
    
    if not os.path.exists(args.output_dir):
        print("Output directory does not exists: %s" % args.output_dir)
        sys.exit()

    #check if taxonomy file exists
    if args.taxonomy_fname and not os.path.isfile(args.taxonomy_fname):
        print "ERROR: Taxonomy file not found: %s" % args.taxonomy_fname
        sys.exit()

    #check if alignment file exists
    if args.align_fname and not os.path.isfile(args.align_fname):
        print "ERROR: Alignment file not found: %s" % args.align_fname
        sys.exit()

    if args.ref_fname and not os.path.isfile(args.ref_fname):
        print("Input reference json file does not exists: %s" % args.ref_fname)
        sys.exit()
    
    if args.jplace_fname and not (os.path.isfile(args.jplace_fname) or os.path.isdir(args.jplace_fname)):
        print("EPA placement file does not exists: %s" % args.jplace_fname)
        sys.exit()

    if args.min_lhw < 0 or args.min_lhw > 1.0:
         args.min_lhw = 0.0
    
    if args.conf_cutoff < 0 or args.conf_cutoff > 1.0:
         args.min_lhw = 0.0

    sativa_home = os.path.dirname(os.path.abspath(__file__))
    if not args.config_fname:
        args.config_fname = os.path.join(sativa_home, "sativa.cfg")
    if not args.temp_dir:
        args.temp_dir = os.path.join(sativa_home, "tmp")
    if not args.output_name:
        if args.taxonomy_fname:
            base_fname = args.taxonomy_fname
        else:
            base_fname = args.ref_fname
        args.output_name = os.path.splitext(base_fname)[0]
        
def print_run_info(config):
    print ""
    config.print_version("SATIVA")
    
    call_str = " ".join(sys.argv)
    config.log.info("SATIVA was called as follows:\n\n%s\n" % call_str)
    
    if config.verbose:
        config.log.info("Mislabels search is running with the following parameters:")
        if config.align_fname:
            config.log.info(" Alignment:                        %s", config.align_fname)
            config.log.info(" Taxonomy:                         %s", config.taxonomy_fname)
        if config.load_refjson:
            config.log.info(" Reference:                        %s", config.refjson_fname)
        if config.jplace_fname:
            config.log.info(" EPA jplace file:                  %s", config.jplace_fname)
        #config.log.info(" Min likelihood weight:            %f", args.min_lhw)
#        config.log.info(" Assignment method:                %s", args.method)
        config.log.info(" Output directory:                 %s", os.path.abspath(config.output_dir))
        config.log.info(" Job name / output files prefix:   %s", config.name)
        config.log.info(" Model of rate heterogeneity:      %s", config.raxml_model)
        config.log.info(" Confidence cut-off:               %f", config.conf_cutoff)
#        config.log.info(" P-value for branch length test:   %g", config.brlen_pv)
        config.log.info(" Number of threads:                %d", config.num_threads)
        config.log.info("")

    if config.debug:
        config.log.debug("Running in DEBUG mode, temp files will be saved to: %s\n", os.path.abspath(config.temp_dir))

if __name__ == "__main__":
    args = parse_args()
    config = SativaConfig(args)
    
    start_time = time.time()
    trainer_time = 0
    
    t = LeaveOneTest(config)
    print_run_info(config)

    if config.load_refjson:
        t.load_refjson(config.refjson_fname)
    else:
        config.log.info("*** STEP 1: Building the reference tree using provided alignment and taxonomic annotations ***\n")
        tr_start_time = time.time() 
        t.run_epa_trainer()
        trainer_time = time.time() - tr_start_time
        t.load_refjson(config.refjson_fname)
        config.log.info("*** STEP 2: Searching for mislabels ***\n")
    
    l1out_start_time = time.time()
    
    t.run_test()
    
    config.clean_tempdir()
        
    l1out_time = time.time() - l1out_start_time

    config.log.info("\nResults were saved to: %s", os.path.abspath(t.mis_fname))
    config.log.info("Execution log was saved to: %s\n", os.path.abspath(config.log_fname))

    elapsed_time = time.time() - start_time
    config.log.info("Analysis completed successfully, elapsed time: %.0f seconds (%.0fs reftree, %.0fs leave-one-out)\n", elapsed_time, trainer_time, l1out_time)
