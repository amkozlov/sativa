#! /usr/bin/env python

import sys
import os
import time
import glob
import multiprocessing
from operator import itemgetter
from subprocess import call

from epac.ete2 import Tree, SeqGroup
from epac.argparse import ArgumentParser
from epac.config import EpacConfig
from epac.raxml_util import RaxmlWrapper, FileUtils
from epac.json_util import RefJsonParser, RefJsonChecker, EpaJsonParser
from epac.taxonomy_util import TaxCode, Taxonomy
from epac.classify_util import TaxTreeHelper,TaxClassifyHelper
from epac.version import SATIVA_BUILD,SATIVA_RELEASE_DATE,SATIVA_RAXML_VER

class LeaveOneTest:
    def __init__(self, config, args):
        self.cfg = config
        self.method = args.method
        self.minlw = args.min_lhw
        self.jplace_fname = args.jplace_fname
        self.ranktest = args.ranktest
        
        self.mis_fname = self.cfg.out_fname("%NAME%.mis")
        self.premis_fname = self.cfg.out_fname("%NAME%.premis")
        self.misrank_fname = self.cfg.out_fname("%NAME%.misrank")
        self.stats_fname = self.cfg.out_fname("%NAME%.stats")
        
        if os.path.isfile(self.mis_fname):
            print "\nERROR: Output file already exists: %s" % self.mis_fname
            print "Please specify a different job name using -n or remove old output files."
            sys.exit()

        self.tmp_refaln = config.tmp_fname("%NAME%.refaln")
        self.reftree_lbl_fname = config.tmp_fname("%NAME%_lbl.tre")
        self.reftree_tax_fname = config.tmp_fname("%NAME%_tax.tre")
        self.optmod_fname = self.cfg.tmp_fname("%NAME%.opt")
        self.reftree_fname = self.cfg.tmp_fname("ref_%NAME%.tre")

        # switch off branch length filter
        self.brlen_pv = 0.

        self.mislabels = []
        self.mislabels_cnt = []
        self.rank_mislabels = []
        self.rank_mislabels_cnt = []
        self.misrank_conf_map = {}

    def load_refjson(self, refjson_fname):
        try:
            self.refjson = RefJsonParser(refjson_fname, ver="1.3")
        except ValueError:
            print("ERROR: Invalid json file format!")
            sys.exit()
            
        #validate input json format 
        self.refjson.validate()
        self.rate = self.refjson.get_rate()
        self.node_height = self.refjson.get_node_height()
        self.origin_taxonomy = self.refjson.get_origin_taxonomy()
        self.bid_taxonomy_map = self.refjson.get_bid_tanomomy_map()
        self.tax_tree = self.refjson.get_tax_tree()
        self.cfg.compress_patterns = self.refjson.get_pattern_compression()

        reftree_str = self.refjson.get_raxml_readable_tree()
        self.reftree = Tree(reftree_str)
        self.reftree_size = len(self.reftree.get_leaves())

        # IMPORTANT: set EPA heuristic rate based on tree size!                
        self.cfg.resolve_auto_settings(self.reftree_size)
        # If we're loading the pre-optimized model, we MUST set the same rate het. mode as in the ref file        
        if self.cfg.epa_load_optmod:
            self.cfg.raxml_model = self.refjson.get_ratehet_model()

        self.classify_helper = TaxClassifyHelper(self.cfg, self.bid_taxonomy_map, self.brlen_pv, self.rate, self.node_height)
        
        tax_code_name = self.refjson.get_taxcode()
        self.tax_code = TaxCode(tax_code_name)

        self.mislabels_cnt = [0] * TaxCode.UNI_TAX_LEVELS
        self.rank_mislabels_cnt = [0] * TaxCode.UNI_TAX_LEVELS
        
    def run_epa_trainer(self, refjson_fname, args):
        sativa_home = os.path.dirname(os.path.abspath(__file__))
        trainer_script = os.path.join(sativa_home, "epac", "epa_trainer.py")
        call_str = [trainer_script]
        call_str += ["-s", args.align_fname]
        call_str += ["-t", args.taxonomy_fname]
        call_str += ["-x", args.taxcode_name]
        call_str += ["-r", refjson_fname]
        call_str += ["-T", str(args.num_threads)]
        
        call_str += ["-c", args.config_fname]
        call_str += ["-o", args.output_dir]
        call_str += ["-n", self.cfg.name]
        call_str += ["-tmpdir", args.temp_dir]
        
        call_str += ["-C", "-no-hmmer"]
        
        if args.verbose:
            call_str += ["-v"]
        if args.debug:
            call_str += ["-debug"]
       
        #print call_str
        call(call_str)
        
        if not os.path.isfile(refjson_fname):
            self.cfg.log.error("\nBuilding reference tree failed, see error messages above.")
            sys.exit(-1)

    def cleanup(self):
        FileUtils.remove_if_exists(self.tmp_refaln)

    def classify_seq(self, placement):
        edges = placement["p"]
        if len(edges) > 0:
            return self.classify_helper.classify_seq(edges, self.method, self.minlw)
        else:
            print "ERROR: no placements! something is definitely wrong!"

    def check_seq_tax_labels(self, seq_name, orig_ranks, ranks, lws):
        mislabel_lvl = -1
        min_len = min(len(orig_ranks),len(ranks))
        for rank_lvl in range(min_len):
            if ranks[rank_lvl] != Taxonomy.EMPTY_RANK and ranks[rank_lvl] != orig_ranks[rank_lvl]:
                mislabel_lvl = rank_lvl
                break

        if mislabel_lvl >= 0:
            real_lvl = self.tax_code.guess_rank_level(orig_ranks, mislabel_lvl)
            mis_rec = {}
            mis_rec['name'] = EpacConfig.strip_ref_prefix(seq_name)
            mis_rec['orig_level'] = mislabel_lvl
            mis_rec['real_level'] = real_lvl
            mis_rec['level_name'] = self.tax_code.rank_level_name(real_lvl)[0]
            mis_rec['inv_level'] = -1 * real_lvl  # just for sorting
            mis_rec['orig_ranks'] = orig_ranks
            mis_rec['ranks'] = ranks
            mis_rec['lws'] = lws
            mis_rec['conf'] = lws[mislabel_lvl]
            self.mislabels.append(mis_rec)
            
            return mis_rec
        else:
            return None

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
        output = mis_rec['name'] + "\t"
        output += "%s\t%s\t%s\t%.3f\t" % (mis_rec['level_name'], 
            mis_rec['orig_ranks'][lvl], mis_rec['ranks'][lvl], mis_rec['lws'][lvl])
        output += Taxonomy.lineage_str(mis_rec['orig_ranks']) + "\t"
        output += Taxonomy.lineage_str(mis_rec['ranks']) + "\t"
        output += ";".join(["%.3f" % conf for conf in mis_rec['lws']])
        if 'rank_conf' in mis_rec:
            output += "\t%.3f" % mis_rec['rank_conf']
        return output

    def sort_mislabels(self):
        self.mislabels = sorted(self.mislabels, key=itemgetter('inv_level', 'conf'), reverse=True)
        for mis_rec in self.mislabels:
            real_lvl = mis_rec["real_level"]
            self.mislabels_cnt[real_lvl] += 1
        
        if self.ranktest:
            self.rank_mislabels = sorted(self.rank_mislabels, key=itemgetter('inv_level', 'conf'), reverse=True)
            for mis_rec in self.rank_mislabels:
                real_lvl = mis_rec["real_level"]
                self.rank_mislabels_cnt[real_lvl] += 1
    
    def write_mislabels(self, final=True):
        if final:
            out_fname = self.mis_fname
        else:
            out_fname = self.premis_fname
        
        with open(out_fname, "w") as fo_all:
            fields = ["SeqID", "MislabeledLevel", "OriginalLabel", "ProposedLabel", "Confidence", "OriginalTaxonomyPath", "ProposedTaxonomyPath", "PerRankConfidence"]
            if self.ranktest:
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

        if self.ranktest:
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

        self.cfg.log.info("Mislabels counts by ranks:")
        with open(self.stats_fname, "w") as fo_stat:
            seq_sum = 0
            rank_sum = 0
            for i in range(1, len(self.mislabels_cnt)):
                rname = self.tax_code.rank_level_name(i)[0].ljust(10)
                if self.mislabels_cnt[i] > 0:
                    seq_sum += self.mislabels_cnt[i]
#                    output = "%s:\t%d" % (rname, seq_sum)
                    output = "%s:\t%d" % (rname, self.mislabels_cnt[i])
                    if self.ranktest:
                        rank_sum += self.rank_mislabels_cnt[i]
                        output += "\t%d" % rank_sum
                    fo_stat.write(output + "\n")
                    self.cfg.log.info(output) 
       
    def get_orig_ranks(self, seq_name):
        nodes = self.tax_tree.get_leaves_by_name(seq_name)
        if len(nodes) != 1:
            print "FATAL ERROR: Sequence %s is not found in the taxonomic tree, or is present more than once!" % seq_name
            sys.exit()
        seq_node = nodes[0]
        orig_ranks = Taxonomy.split_rank_uid(seq_node.up.name)
        return orig_ranks
    
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
        if self.jplace_fname:
            jp = EpaJsonParser(self.jplace_fname)
        else:        
            jp = self.raxml.run_epa(job_name, self.refalign_fname, self.reftree_fname, self.optmod_fname, mode="l1o_seq")

        placements = jp.get_placement()
        seq_count = 0
        for place in placements:
            seq_name = place["n"][0]
            
            # get original taxonomic label
            orig_ranks = self.get_orig_ranks(seq_name)

            # get EPA tax label
            ranks, lws = self.classify_seq(place)
            # check if they match
            mis_rec = self.check_seq_tax_labels(seq_name, orig_ranks, ranks, lws)
            # cross-check with higher rank mislabels
            if self.ranktest and mis_rec:
                rank_conf = 0
                for lvl in range(2,len(orig_ranks)):
                    tax_path = Taxonomy.get_rank_uid(orig_ranks, lvl)
                    if tax_path in self.misrank_conf_map:
                        rank_conf = max(rank_conf, self.misrank_conf_map[tax_path])
                mis_rec['rank_conf'] = rank_conf
            seq_count += 1

        if not self.cfg.debug:
            self.raxml.cleanup(job_name)

        return seq_count    
        
    def run_final_epa_test(self):
        self.reftree_outgroup = self.refjson.get_outgroup()
        tmp_reftree = self.reftree.copy() 
        tmp_taxtree = self.tax_tree.copy() 
        for mis_rec in self.mislabels:
            name = mis_rec['name']
            rname = EpacConfig.REF_SEQ_PREFIX + name

            leaf_nodes = tmp_reftree.get_leaves_by_name(rname)
            if len(leaf_nodes) > 0:
                leaf_nodes[0].delete()
            else:
                print "Node not found in the reference tree: %s" % rname

            leaf_nodes = tmp_taxtree.get_leaves_by_name(rname)
            if len(leaf_nodes) > 0:
                leaf_nodes[0].delete()
            else:
                print "Node not found in the taxonomic tree: %s" % rname

        # remove unifurcation at the root
        if len(tmp_reftree.children) == 1:
            tmp_reftree = tmp_reftree.children[0]
            
        self.mislabels = []

        th = TaxTreeHelper(self.origin_taxonomy, self.cfg)
        th.set_mf_rooted_tree(tmp_taxtree)
            
        self.run_epa_once(tmp_reftree, th)
            

    def run_epa_once(self, reftree, th):
        reftree_fname = self.cfg.tmp_fname("final_ref_%NAME%.tre")
        job_name = self.cfg.subst_name("final_epa_%NAME%")

        reftree.write(outfile=reftree_fname)

        # IMPORTANT: don't load the model, since it's invalid for the pruned true !!! 
        optmod_fname=""
        epa_result = self.raxml.run_epa(job_name, self.refalign_fname, reftree_fname, optmod_fname)
        reftree_epalbl_str = epa_result.get_std_newick_tree()        
        placements = epa_result.get_placement()
        
        # update branchid-taxonomy mapping to account for possible changes in branch numbering
        reftree_tax = Tree(reftree_epalbl_str)
        th.set_bf_unrooted_tree(reftree_tax)
        bid_tax_map = th.get_bid_taxonomy_map()
        
        cl = TaxClassifyHelper(self.cfg, bid_tax_map, self.brlen_pv, self.rate, self.node_height)

        for place in placements:
            seq_name = place["n"][0]

            # get original taxonomic label
            orig_ranks = self.get_orig_ranks(seq_name)
            # get EPA tax label
            ranks, lws = cl.classify_seq(place["p"])
            # check if they match
            mis_rec = self.check_seq_tax_labels(seq_name, orig_ranks, ranks, lws)

        if not self.cfg.debug:
            self.raxml.cleanup(job_name)
            FileUtils.remove_if_exists(reftree_fname)
            
    def run_test(self):
        self.raxml = RaxmlWrapper(self.cfg)

#        config.log.info("Number of sequences in the reference: %d\n", self.reftree_size)

        self.refjson.get_raxml_readable_tree(self.reftree_fname)
        self.refalign_fname = self.refjson.get_alignment(self.tmp_refaln)        
        self.refjson.get_binary_model(self.optmod_fname)

        if self.ranktest:
            config.log.info("Running the leave-one-rank-out test...\n")
            subtree_count = self.run_leave_subtree_out_test()
            
        config.log.info("Running the leave-one-sequence-out test...\n")
        self.run_leave_seq_out_test()

        if len(self.mislabels) > 0:
            config.log.info("Leave-one-out test identified %d suspicious sequences; running final EPA test to check them...\n", len(self.mislabels))
            if self.cfg.debug:
                self.write_mislabels(final=False)
            self.run_final_epa_test()

        self.sort_mislabels()
        self.write_mislabels()
        config.log.info("\nPercentage of mislabeled sequences: %.2f %%", (float(len(self.mislabels)) / self.reftree_size * 100))

        if not self.cfg.debug:
            FileUtils.remove_if_exists(self.reftree_fname)
            FileUtils.remove_if_exists(self.optmod_fname)
            FileUtils.remove_if_exists(self.refalign_fname)

def parse_args():
    parser = ArgumentParser(description="Find putative mislabeled/misplaced sequences in a taxonomy.",
    epilog="Example: python sativa.py -s example/test.phy -t example/test.tax -x BAC")
    parser.add_argument("-s", dest="align_fname",
            help="""Reference alignment file. Sequences must be aligned, their IDs must correspond to those
in taxonomy file.""")
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
    parser.add_argument("-v", dest="verbose", action="store_true",
            help="""Print additional info messages to the console.""")

    parser.add_argument("-c", dest="config_fname", default=None,
            help="Config file name.")
    parser.add_argument("-r", dest="ref_fname",
            help="""Specify the reference alignment and taxonomy in json format.""")
    parser.add_argument("-j", dest="jplace_fname", default=None,
            help="""Do not call RAxML EPA, use existing .jplace file as input instead.""")
    parser.add_argument("-l", dest="min_lhw", type=float, default=0.,
            help="""A value between 0 and 1, the minimal sum of likelihood weight of
                    an assignment to a specific rank. This value represents a confidence 
                    measure of the assignment, assignments below this value will be discarded. 
                    Default: 0 to output all possbile assignments.""")
    parser.add_argument("-m", dest="method", default="1",
            help="""Assignment method 1 or 2
                    1: Max sum likelihood (default)
                    2: Max likelihood placement""")
#    parser.add_argument("-p", dest="p_value", type=float, default=0.001,
#            help="""P-value for branch length Erlang test. Default: 0.001\n""")

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
    
    if args.jplace_fname and not os.path.isfile(args.jplace_fname):
        print("EPA placement file does not exists: %s" % args.jplace_fname)
        sys.exit()

    if args.min_lhw < 0 or args.min_lhw > 1.0:
         args.min_lhw = 0.0
    
    if not (args.method == "1" or args.method == "2"):
        args.method == "1"

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
        
def print_run_info(config, args):
    print ""
    config.log.info("SATIVA %s, released on %s. Last version: https://github.com/amkozlov/sativa", SATIVA_BUILD, SATIVA_RELEASE_DATE)
    config.log.info("By A.Kozlov and J.Zhang, the Exelixis Lab. Based on RAxML %s by A.Stamatakis.\n", SATIVA_RAXML_VER)
    
    if config.verbose:
        config.log.info("Mislabels search is running with the following parameters:")
        if args.align_fname:
            config.log.info(" Alignment:                        %s", args.align_fname)
            config.log.info(" Taxonomy:                         %s", args.taxonomy_fname)
        if args.ref_fname:
            config.log.info(" Reference:                        %s", args.ref_fname)
        if args.jplace_fname:
            config.log.info(" EPA jplace file:                  %s", args.jplace_fname)
        #config.log.info(" Min likelihood weight:            %f", args.min_lhw)
#        config.log.info(" Assignment method:                %s", args.method)
    #    print(" P-value for branch length test:   %f" % args.p_value)
        config.log.info(" Output directory:                 %s", os.path.abspath(args.output_dir))
        config.log.info(" Job name / output files prefix:   %s", config.name)
        config.log.info(" Number of threads:                %d", config.num_threads)
        config.log.info("")

if __name__ == "__main__":
    args = parse_args()
    config = EpacConfig(args)
    
    start_time = time.time()
    trainer_time = 0
    
    t = LeaveOneTest(config, args)
    print_run_info(config, args)

    if config.refjson_fname:
        t.load_refjson(config.refjson_fname)
    else:
        refjson_fname = config.tmp_fname("%NAME%.refjson")
        config.log.info("*** STEP 1: Building the reference tree using provided alignment and taxonomic annotations ***\n")
        tr_start_time = time.time() 
        t.run_epa_trainer(refjson_fname, args)
        trainer_time = time.time() - tr_start_time
        t.load_refjson(refjson_fname)
        if not config.debug:
            FileUtils.remove_if_exists(refjson_fname)
        config.log.info("*** STEP 2: Searching for mislabels ***\n")
    
    l1out_start_time = time.time()
    
    t.run_test()
    if not config.debug:
        t.cleanup()
        
    l1out_time = time.time() - l1out_start_time

    config.log.info("\nResults were saved to: %s", os.path.abspath(t.mis_fname))
    config.log.info("Execution log was saved to: %s\n", os.path.abspath(config.log_fname))

    elapsed_time = time.time() - start_time
    config.log.info("Analysis completed successfully, elapsed time: %.0f seconds (%.0fs reftree, %.0fs leave-one-out)\n", elapsed_time, trainer_time, l1out_time)
