#!/usr/bin/env python

import sys
import os
import shutil
import datetime
import time
import logging
import multiprocessing
from string import maketrans

from epac.ete2 import Tree, SeqGroup
from epac.argparse import ArgumentParser,RawTextHelpFormatter
from epac.config import EpacConfig,EpacTrainerConfig
from epac.raxml_util import RaxmlWrapper, FileUtils
from epac.taxonomy_util import Taxonomy, TaxTreeBuilder
from epac.json_util import RefJsonBuilder
from epac.erlang import tree_param 
from epac.msa import hmmer
from epac.classify_util import TaxTreeHelper

class InputValidator:
    def __init__(self, config, input_tax, input_seqs, verbose=True): 
        self.cfg = config
        self.taxonomy = input_tax
        self.alignment = input_seqs
        self.verbose = verbose
        self.dupseq_sets = None
        self.merged_ranks = None
        self.corr_seqid = {}
        self.corr_ranks = {}
        self.gaps_trantab = maketrans("?N", "--")

    def validate(self):
        # following two checks are obsolete and disabled by default
        self.check_tax_disbalance()
        self.check_tax_duplicates()
        
        self.check_seq_ids()
        self.check_invalid_chars()
        
        self.check_identical_seqs()
        self.check_identical_ranks()
        
        self.taxonomy.close_taxonomy_gaps()
        
        return self.corr_ranks, self.corr_seqid, self.merged_ranks

    def normalize_gaps(self, seq):
        return seq.translate(self.gaps_trantab)
        
    def check_seq_ids(self):
        # check that seq IDs in taxonomy and alignment correspond
        self.mis_ids = []
        for sid in self.taxonomy.seq_ranks_map.iterkeys():
            unprefixed_sid = EpacConfig.strip_ref_prefix(sid)
            if not self.alignment.has_seq(unprefixed_sid):
                self.mis_ids.append(unprefixed_sid)
                
        if len(self.mis_ids) > 0 and self.verbose:
            errmsg = "ERROR: Following %d sequence(s) are missing in your alignment file:\n%s\n\n" % (len(self.mis_ids), "\n".join(self.mis_ids))
            errmsg += "Please make sure sequence IDs in taxonomic annotation file and in alignment are identical!\n"
            self.cfg.exit_user_error(errmsg)
            
        return self.mis_ids

    def check_invalid_chars(self):
        # check for invalid characters in rank names
        self.corr_ranks = self.taxonomy.normalize_rank_names()
        
        # check for invalid characters in sequence IDs
        self.corr_seqid = self.taxonomy.normalize_seq_ids()

        if self.verbose:
            for old_rank in sorted(self.corr_ranks.keys()):
                self.cfg.log.debug("NOTE: Following rank name contains illegal symbols and was renamed: %s --> %s", old_rank, self.corr_ranks[old_rank])
            if len(self.corr_ranks) > 0:
                self.cfg.log.debug("")
            for old_sid in sorted(self.corr_seqid.keys()):
                self.cfg.log.debug("NOTE: Following sequence ID contains illegal symbols and was renamed: %s --> %s" , old_sid, self.corr_seqid[old_sid])
            if len(self.corr_seqid) > 0:
                self.cfg.log.debug("")
            
        return self.corr_ranks, self.corr_seqid
        
    def check_identical_seqs(self):
        seq_hash_map = {}
        for name, seq, comment, sid in self.alignment.iter_entries():
            ref_seq_name = EpacConfig.REF_SEQ_PREFIX + name
            ref_seq_name = self.corr_seqid.get(ref_seq_name, ref_seq_name)
            if ref_seq_name in self.taxonomy.seq_ranks_map:
                seq_hash = hash(self.normalize_gaps(seq))
                if seq_hash in seq_hash_map:
                    seq_hash_map[seq_hash] += [name]
                else:
                    seq_hash_map[seq_hash] = [name]

        self.dupseq_count = 0
        self.dupseq_sets = []
        for seq_hash, seq_ids in seq_hash_map.iteritems():
            check_ids = seq_ids[:]
            while len(check_ids) > 1:
                # compare actual sequence strings, to account for a possible hash collision
                seq1 = self.normalize_gaps(self.alignment.get_seq(check_ids[0]))
                coll_ids = []
                dup_ids = [check_ids[0]]
                for i in range(1, len(check_ids)):
                    seq2 = self.normalize_gaps(self.alignment.get_seq(check_ids[i]))
                    if seq1 == seq2:
                        dup_ids += [check_ids[i]]
                    else:
                        # collision found, add put seq id on a list to be checked in the next iteration
                        coll_ids += [check_ids[i]]
 
                if len(dup_ids) > 1:
                    self.dupseq_sets += [dup_ids]
                    self.dupseq_count += len(dup_ids) - 1

                check_ids = coll_ids
                
        if self.verbose:
            for dup_ids in self.dupseq_sets:
                self.cfg.log.debug("NOTE: Following sequences are identical: %s", ", ".join(dup_ids))
            if self.dupseq_count > 0:
                self.cfg.log.debug("\nNOTE: Found %d sequence duplicates", self.dupseq_count)
                
        return self.dupseq_count, self.dupseq_sets

    def check_identical_ranks(self):
        if not self.dupseq_sets:
            self.check_identical_seqs()
        self.merged_ranks = {}
        for dup_ids in self.dupseq_sets:
            if len(dup_ids) > 1:
                duprank_map = {}
                for seq_name in dup_ids:
                    rank_id = self.taxonomy.seq_rank_id(seq_name)
                    duprank_map[rank_id] = duprank_map.get(rank_id, 0) + 1
                if len(duprank_map) > 1 and self.cfg.debug:
                    self.cfg.log.debug("Ranks sharing duplicates: %s\n", str(duprank_map))
                dup_ranks = []
                for rank_id, count in duprank_map.iteritems():
                    if count > self.cfg.taxa_ident_thres * self.taxonomy.get_rank_seq_count(rank_id):
                      dup_ranks += [rank_id]
                if len(dup_ranks) > 1:
                    prefix = "__TAXCLUSTER%d__" % (len(self.merged_ranks) + 1)
                    merged_rank_id = self.taxonomy.merge_ranks(dup_ranks, prefix)
                    self.merged_ranks[merged_rank_id] = dup_ranks

        if self.verbose:
            merged_count = 0
            for merged_rank_id, dup_ranks in self.merged_ranks.iteritems():
                dup_ranks_str = "\n".join([Taxonomy.rank_uid_to_lineage_str(rank_id) for rank_id in dup_ranks])
                self.cfg.log.warning("\nWARNING: Following taxa share >%.0f%% indentical sequences und thus considered indistinguishable:\n%s", self.cfg.taxa_ident_thres * 100, dup_ranks_str)
                merged_rank_str = Taxonomy.rank_uid_to_lineage_str(merged_rank_id)
                self.cfg.log.warning("For the purpose of mislabels identification, they were merged into one taxon:\n%s\n", merged_rank_str)
                merged_count += len(dup_ranks)
            
            if merged_count > 0:
                self.cfg.log.warning("WARNING: %d indistinguishable taxa have been merged into %d clusters.\n", merged_count, len(self.merged_ranks))

        return self.merged_ranks
            
    def check_tax_disbalance(self):
        # make sure we don't taxonomy "irregularities" (more than 7 ranks or missing ranks in the middle)
        action = self.cfg.wrong_rank_count
        if action != "ignore":
            autofix = action == "autofix"
            errs = self.taxonomy.check_for_disbalance(autofix)
            if len(errs) > 0:
                if action == "autofix":
                    print "WARNING: %d sequences with invalid annotation (missing/redundant ranks) found and were fixed as follows:\n" % len(errs)
                    for err in errs:
                        print "Original:   %s\t%s"   % (err[0], err[1])
                        print "Fixed as:   %s\t%s\n" % (err[0], err[2])
                elif action == "skip":
                    print "WARNING: Following %d sequences with invalid annotation (missing/redundant ranks) were skipped:\n" % len(errs)
                    for err in errs:
                        self.taxonomy.remove_seq(err[0])
                        print "%s\t%s" % err
                else:  # abort
                    print "ERROR: %d sequences with invalid annotation (missing/redundant ranks) found:\n" % len(errs)
                    for err in errs:
                        print "%s\t%s" % err
                    print "\nPlease fix them manually (add/remove ranks) and run the pipeline again (or use -wrong-rank-count autofix option)"
                    print "NOTE: Only standard 7-level taxonomies are supported at the moment. Although missing trailing ranks (e.g. species) are allowed,"
                    print "missing intermediate ranks (e.g. family) or sublevels (e.g. suborder) are not!\n"
                    self.cfg.exit_user_error()
        
    def check_tax_duplicates(self):
        # check for duplicate rank names
        action = self.cfg.dup_rank_names
        if action != "ignore":
            autofix = action == "autofix"
            dups = self.taxonomy.check_for_duplicates(autofix)
            if len(dups) > 0:
                if action == "autofix":
                    print "WARNING: %d sequences with duplicate rank names found and were renamed as follows:\n" % len(dups)
                    for dup in dups:
                        print "Original:    %s\t%s"   %  (dup[0], dup[1])
                        print "Duplicate:   %s\t%s"   %  (dup[2], dup[3])
                        print "Renamed to:  %s\t%s\n" %  (dup[2], dup[4])
                elif action == "skip":
                    print "WARNING: Following %d sequences with duplicate rank names were skipped:\n" % len(dups)
                    for dup in dups:
                        self.taxonomy.remove_seq(dup[2])
                        print "%s\t%s\n" % (dup[2], dup[3])
                else:  # abort
                    print "ERROR: %d sequences with duplicate rank names found:\n" % len(dups)
                    for dup in dups:
                        print "%s\t%s\n%s\t%s\n" % dup
                    print "Please fix (rename) them and run the pipeline again (or use -dup-rank-names autofix option)" 
                    self.cfg.exit_user_error()

class RefTreeBuilder:
    def __init__(self, config): 
        self.cfg = config
        self.mfresolv_job_name = self.cfg.subst_name("mfresolv_%NAME%")
        self.epalbl_job_name = self.cfg.subst_name("epalbl_%NAME%")
        self.optmod_job_name = self.cfg.subst_name("optmod_%NAME%")
        self.raxml_wrapper = RaxmlWrapper(config)
        
        self.outgr_fname = self.cfg.tmp_fname("%NAME%_outgr.tre")
        self.reftree_mfu_fname = self.cfg.tmp_fname("%NAME%_mfu.tre")
        self.reftree_bfu_fname = self.cfg.tmp_fname("%NAME%_bfu.tre")
        self.optmod_fname = self.cfg.tmp_fname("%NAME%.opt")
        self.lblalign_fname = self.cfg.tmp_fname("%NAME%_lblq.fa")
        self.reftree_lbl_fname = self.cfg.tmp_fname("%NAME%_lbl.tre")
        self.reftree_tax_fname = self.cfg.tmp_fname("%NAME%_tax.tre")
        self.brmap_fname = self.cfg.tmp_fname("%NAME%_map.txt")

    def load_alignment(self):
        in_file = self.cfg.align_fname
        self.input_seqs = None
        formats = ["fasta", "phylip_relaxed", "iphylip_relaxed", "phylip", "iphylip"]
        for fmt in formats:
            try:
                self.input_seqs = SeqGroup(sequences=in_file, format = fmt)
                break
            except:
                self.cfg.log.debug("Guessing input format: not " + fmt)
        if self.input_seqs == None:
            self.cfg.exit_user_error("Invalid input file format: %s\nThe supported input formats are fasta and phylip" % in_file)
            
    def validate_taxonomy(self):
        self.input_validator = InputValidator(self.cfg, self.taxonomy, self.input_seqs)
        self.input_validator.validate()
        
    def build_multif_tree(self):
        c = self.cfg
        
        tb = TaxTreeBuilder(c, self.taxonomy)
        (t, ids) = tb.build(c.reftree_min_rank, c.reftree_max_seqs_per_leaf, c.reftree_clades_to_include, c.reftree_clades_to_ignore)
        self.reftree_ids = frozenset(ids)
        self.reftree_size = len(ids)
        self.reftree_multif = t

        # IMPORTANT: select GAMMA or CAT model based on tree size!                
        self.cfg.resolve_auto_settings(self.reftree_size)

        if self.cfg.debug:
            refseq_fname = self.cfg.tmp_fname("%NAME%_seq_ids.txt")
            # list of sequence ids which comprise the reference tree
            with open(refseq_fname, "w") as f:
                for sid in ids:
                    f.write("%s\n" % sid)

            # original tree with taxonomic ranks as internal node labels
            reftax_fname = self.cfg.tmp_fname("%NAME%_mfu_tax.tre")
            t.write(outfile=reftax_fname, format=8)
        #    t.show()

    def export_ref_alignment(self):
        """This function transforms the input alignment in the following way:
           1. Filter out sequences which are not part of the reference tree
           2. Add sequence name prefix (r_)"""
        
        self.refalign_fname = self.cfg.tmp_fname("%NAME%_matrix.afa")
        with open(self.refalign_fname, "w") as fout:
            for name, seq, comment, sid in self.input_seqs.iter_entries():
                seq_name = EpacConfig.REF_SEQ_PREFIX + name
                if seq_name in self.input_validator.corr_seqid:
                  seq_name = self.input_validator.corr_seqid[seq_name]
                if seq_name in self.reftree_ids:
                    fout.write(">" + seq_name + "\n" + seq + "\n")

        # we do not need the original alignment anymore, so free its memory
        self.input_seqs = None

    def export_ref_taxonomy(self):
        self.taxonomy_map = {}
        
        for sid, ranks in self.taxonomy.iteritems():
            if sid in self.reftree_ids:
                self.taxonomy_map[sid] = ranks
            
        if self.cfg.debug:
            tax_fname = self.cfg.tmp_fname("%NAME%_tax.txt")
            with open(tax_fname, "w") as fout:
                for sid, ranks in self.taxonomy_map.iteritems():
                    ranks_str = self.taxonomy.seq_lineage_str(sid) 
                    fout.write(sid + "\t" + ranks_str + "\n")   

    def save_rooting(self):
        rt = self.reftree_multif

        tax_map = self.taxonomy.get_map()
        self.taxtree_helper = TaxTreeHelper(self.cfg, tax_map)
        self.taxtree_helper.set_mf_rooted_tree(rt)
        outgr = self.taxtree_helper.get_outgroup()
        outgr_size = len(outgr.get_leaves())
        outgr.write(outfile=self.outgr_fname, format=9)
        self.reftree_outgroup = outgr
        self.cfg.log.debug("Outgroup for rooting was saved to: %s, outgroup size: %d", self.outgr_fname, outgr_size)
            
        # remove unifurcation at the root
        if len(rt.children) == 1:
            rt = rt.children[0]
        
        # now we can safely unroot the tree and remove internal node labels to make it suitable for raxml
        rt.write(outfile=self.reftree_mfu_fname, format=9)

    # RAxML call to convert multifurcating tree to the strictly bifurcating one
    def resolve_multif(self):
        self.cfg.log.debug("\nReducing the alignment: \n")
        self.reduced_refalign_fname = self.raxml_wrapper.reduce_alignment(self.refalign_fname)
        
        self.cfg.log.debug("\nConstrained ML inference: \n")
        raxml_params = ["-s", self.reduced_refalign_fname, "-g", self.reftree_mfu_fname, "--no-seq-check", "-N", str(self.cfg.rep_num)] 
        if self.cfg.mfresolv_method  == "fast":
            raxml_params += ["-D"]
        elif self.cfg.mfresolv_method  == "ultrafast":
            raxml_params += ["-f", "e"]
        if self.cfg.restart and self.raxml_wrapper.result_exists(self.mfresolv_job_name):
            self.invocation_raxml_multif = self.raxml_wrapper.get_invocation_str(self.mfresolv_job_name)
            self.cfg.log.debug("\nUsing existing ML tree found in: %s\n", self.raxml_wrapper.result_fname(self.mfresolv_job_name))
        else:
            self.invocation_raxml_multif = self.raxml_wrapper.run(self.mfresolv_job_name, raxml_params)
#            self.invocation_raxml_multif = self.raxml_wrapper.run_multiple(self.mfresolv_job_name, raxml_params, self.cfg.rep_num)
            if self.cfg.mfresolv_method  == "ultrafast":
              self.raxml_wrapper.copy_result_tree(self.mfresolv_job_name, self.raxml_wrapper.besttree_fname(self.mfresolv_job_name))
              
        if self.raxml_wrapper.besttree_exists(self.mfresolv_job_name):        
            if not self.cfg.reopt_model:
                self.raxml_wrapper.copy_best_tree(self.mfresolv_job_name, self.reftree_bfu_fname)
                self.raxml_wrapper.copy_optmod_params(self.mfresolv_job_name, self.optmod_fname)
                self.invocation_raxml_optmod = ""
                job_name = self.mfresolv_job_name
            else:
                bfu_fname = self.raxml_wrapper.besttree_fname(self.mfresolv_job_name)
                job_name = self.optmod_job_name

                # RAxML call to optimize model parameters and write them down to the binary model file
                self.cfg.log.debug("\nOptimizing model parameters: \n")
                raxml_params = ["-f", "e", "-s", self.reduced_refalign_fname, "-t", bfu_fname, "--no-seq-check"]
                if self.cfg.raxml_model.startswith("GTRCAT") and not self.cfg.compress_patterns:
                    raxml_params +=  ["-H"]
                if self.cfg.restart and self.raxml_wrapper.result_exists(self.optmod_job_name):
                    self.invocation_raxml_optmod = self.raxml_wrapper.get_invocation_str(self.optmod_job_name)
                    self.cfg.log.debug("\nUsing existing optimized tree and parameters found in: %s\n", self.raxml_wrapper.result_fname(self.optmod_job_name))
                else:
                    self.invocation_raxml_optmod = self.raxml_wrapper.run(self.optmod_job_name, raxml_params)
                if self.raxml_wrapper.result_exists(self.optmod_job_name):
                    self.raxml_wrapper.copy_result_tree(self.optmod_job_name, self.reftree_bfu_fname)
                    self.raxml_wrapper.copy_optmod_params(self.optmod_job_name, self.optmod_fname)
                else:
                    errmsg = "RAxML run failed (model optimization), please examine the log for details: %s" \
                            % self.raxml_wrapper.make_raxml_fname("output", self.optmod_job_name)
                    self.cfg.exit_fatal_error(errmsg)
                    
            if self.cfg.raxml_model.startswith("GTRCAT"):
              mod_name = "CAT"
            else:
              mod_name = "GAMMA" 
            self.reftree_loglh = self.raxml_wrapper.get_tree_lh(job_name, mod_name)
            self.cfg.log.debug("\n%s-based logLH of the reference tree: %f\n" % (mod_name, self.reftree_loglh))

        else:
            errmsg = "RAxML run failed (mutlifurcation resolution), please examine the log for details: %s" \
                    % self.raxml_wrapper.make_raxml_fname("output", self.mfresolv_job_name)
            self.cfg.exit_fatal_error(errmsg)
            
    def load_reduced_refalign(self):
        formats = ["fasta", "phylip_relaxed"]
        for fmt in formats:
            try:
                self.reduced_refalign_seqs = SeqGroup(sequences=self.reduced_refalign_fname, format = fmt)
                break
            except:
                pass
        if self.reduced_refalign_seqs == None:
            errmsg = "FATAL ERROR: Invalid input file format in %s! (load_reduced_refalign)" % self.reduced_refalign_fname
            self.cfg.exit_fatal_error(errmsg)
    
    # dummy EPA run to label the branches of the reference tree, which we need to build a mapping to tax ranks    
    def epa_branch_labeling(self):
        # create alignment with dummy query seq
        self.refalign_width = len(self.reduced_refalign_seqs.get_seqbyid(0))
        self.reduced_refalign_seqs.write(format="fasta", outfile=self.lblalign_fname)
        
        with open(self.lblalign_fname, "a") as fout:
            fout.write(">" + "DUMMY131313" + "\n")        
            fout.write("A"*self.refalign_width + "\n")        
        
        # TODO always load model regardless of the config file settings?
        epa_result = self.raxml_wrapper.run_epa(self.epalbl_job_name, self.lblalign_fname, self.reftree_bfu_fname, self.optmod_fname, mode="epa_mp")
        self.reftree_lbl_str = epa_result.get_std_newick_tree()
        self.raxml_version = epa_result.get_raxml_version()
        self.invocation_raxml_epalbl = epa_result.get_raxml_invocation()

        if not self.raxml_wrapper.epa_result_exists(self.epalbl_job_name):        
            errmsg = "RAxML EPA run failed, please examine the log for details: %s" \
                    % self.raxml_wrapper.make_raxml_fname("output", self.epalbl_job_name)
            self.cfg.exit_fatal_error(errmsg)

    def epa_post_process(self):
        lbl_tree = Tree(self.reftree_lbl_str)
        self.taxtree_helper.set_bf_unrooted_tree(lbl_tree)
        self.reftree_tax = self.taxtree_helper.get_tax_tree()
        self.bid_ranks_map = self.taxtree_helper.get_bid_taxonomy_map()
        
        if self.cfg.debug:
            self.reftree_tax.write(outfile=self.reftree_tax_fname, format=3)
            with open(self.reftree_lbl_fname, "w") as outf:
                outf.write(self.reftree_lbl_str)
            with open(self.brmap_fname, "w") as outf:
                for bid, br_rec in self.bid_ranks_map.iteritems():
                    outf.write("%s\t%s\t%d\t%f\n" % (bid, br_rec[0], br_rec[1], br_rec[2]))

    def calc_node_heights(self):
        """Calculate node heights on the reference tree (used to define branch-length cutoff during classification step)
           Algorithm is as follows:
           Tip node or node resolved to Species level: height = 1 
           Inner node resolved to Genus or above:      height = min(left_height, right_height) + 1 
         """
        nh_map = {}
        dummy_added = False
        for node in self.reftree_tax.traverse("postorder"):
            if not node.is_root():
                if not hasattr(node, "B"):                
                    # In a rooted tree, there is always one more node/branch than in unrooted one
                    # That's why one branch will be always not EPA-labelled after the rooting
                    if not dummy_added: 
                        node.B = "DDD"
                        dummy_added = True
                        species_rank = Taxonomy.EMPTY_RANK
                    else:
                        errmsg = "FATAL ERROR: More than one tree branch without EPA label (calc_node_heights)"
                        self.cfg.exit_fatal_error(errmsg)
                else:
                    species_rank = self.bid_ranks_map[node.B][-1]
                bid = node.B
                if node.is_leaf() or species_rank != Taxonomy.EMPTY_RANK:
                    nh_map[bid] = 1
                else:
                    lchild = node.children[0]
                    rchild = node.children[1]
                    nh_map[bid] = min(nh_map[lchild.B], nh_map[rchild.B]) + 1

        # remove heights for dummy nodes, since there won't be any placements on them
        if dummy_added:
            del nh_map["DDD"]
            
        self.node_height_map = nh_map

    def __get_all_rank_names(self, root):
        rnames = set([])
        for node in root.traverse("postorder"):
            ranks = node.ranks
            for rk in ranks:
                rnames.add(rk)
        return rnames

    def mono_index(self):
        """This method will calculate monophyly index by looking at the left and right hand side of the tree"""
        children = self.reftree_tax.children
        if len(children) == 1:
            while len(children) == 1:
                children = children[0].children 
        if len(children) == 2:
            left = children[0]
            right =children[1]
            lset = self.__get_all_rank_names(left)
            rset = self.__get_all_rank_names(right)
            iset = lset & rset
            return iset
        else:
            print("Error: input tree not birfurcating")
            return set([])

    def build_hmm_profile(self, json_builder):
        print "Building the HMMER profile...\n"

        # this stupid workaround is needed because RAxML outputs the reduced
        # alignment in relaxed PHYLIP format, which is not supported by HMMER
        refalign_fasta = self.cfg.tmp_fname("%NAME%_ref_reduced.fa")
        self.reduced_refalign_seqs.write(outfile=refalign_fasta)

        hmm = hmmer(self.cfg, refalign_fasta)
        fprofile = hmm.build_hmm_profile()

        json_builder.set_hmm_profile(fprofile)
        
    def write_json(self):
        jw = RefJsonBuilder()

        jw.set_branch_tax_map(self.bid_ranks_map)
        jw.set_tree(self.reftree_lbl_str)
        jw.set_outgroup(self.reftree_outgroup)
        jw.set_ratehet_model(self.cfg.raxml_model)
        jw.set_tax_tree(self.reftree_multif)
        jw.set_pattern_compression(self.cfg.compress_patterns)
        jw.set_taxcode(self.cfg.taxcode_name)
        
        jw.set_merged_ranks_map(self.input_validator.merged_ranks)
        corr_ranks_reverse = dict((reversed(item) for item in self.input_validator.corr_ranks.items()))
        jw.set_corr_ranks_map(corr_ranks_reverse)
        corr_seqid_reverse = dict((reversed(item) for item in self.input_validator.corr_seqid.items()))
        jw.set_corr_seqid_map(corr_seqid_reverse)

        mdata = { "ref_tree_size": self.reftree_size, 
                  "ref_alignment_width": self.refalign_width,
                  "raxml_version": self.raxml_version,
                  "timestamp": str(datetime.datetime.now()),
                  "invocation_epac": self.invocation_epac,
                  "invocation_raxml_multif": self.invocation_raxml_multif,
                  "invocation_raxml_optmod": self.invocation_raxml_optmod,
                  "invocation_raxml_epalbl": self.invocation_raxml_epalbl,
                  "reftree_loglh": self.reftree_loglh
                }
        jw.set_metadata(mdata)

        seqs = self.reduced_refalign_seqs.get_entries()    
        jw.set_sequences(seqs)
        
        if not self.cfg.no_hmmer:
            self.build_hmm_profile(jw)

        orig_tax = self.taxonomy_map
        jw.set_origin_taxonomy(orig_tax)
        
        self.cfg.log.debug("Calculating the speciation rate...\n")
        tp = tree_param(tree = self.reftree_lbl_str, origin_taxonomy = orig_tax)
        jw.set_rate(tp.get_speciation_rate_fast())
        jw.set_nodes_height(self.node_height_map)
        
        jw.set_binary_model(self.optmod_fname)
        
        self.cfg.log.debug("Writing down the reference file...\n")
        jw.dump(self.cfg.refjson_fname)

    # top-level function to build a reference tree    
    def build_ref_tree(self):
        self.cfg.log.info("=> Loading taxonomy from file: %s ...\n" , self.cfg.taxonomy_fname)
        self.taxonomy = Taxonomy(prefix=EpacConfig.REF_SEQ_PREFIX, tax_fname=self.cfg.taxonomy_fname)
        self.cfg.log.info("==> Loading reference alignment from file: %s ...\n" , self.cfg.align_fname)
        self.load_alignment()
        self.cfg.log.info("===> Validating taxonomy and alignment ...\n")
        self.validate_taxonomy()
        self.cfg.log.info("====> Building a multifurcating tree from taxonomy with %d seqs ...\n" , self.taxonomy.seq_count())
        self.build_multif_tree()
        self.cfg.log.info("=====> Building the reference alignment ...\n")
        self.export_ref_alignment()
        self.export_ref_taxonomy()
        self.cfg.log.info("======> Saving the outgroup for later re-rooting ...\n")
        self.save_rooting()
        self.cfg.log.info("=======> Resolving multifurcation: choosing the best topology from %d independent RAxML runs ...\n" % self.cfg.rep_num)
        self.resolve_multif()
        self.load_reduced_refalign()
        self.cfg.log.info("========> Calling RAxML-EPA to obtain branch labels ...\n")
        self.epa_branch_labeling()
        self.cfg.log.info("=========> Post-processing the EPA tree (re-rooting, taxonomic labeling etc.) ...\n")
        self.epa_post_process()
        self.calc_node_heights()
        
        self.cfg.log.debug("\n==========> Checking branch labels ...")
        self.cfg.log.debug("shared rank names before training: %s", repr(self.taxonomy.get_common_ranks()))
        self.cfg.log.debug("shared rank names after  training: %s\n", repr(self.mono_index()))
        
        self.cfg.log.info("==========> Saving the reference JSON file: %s\n" % self.cfg.refjson_fname)
        self.write_json()

def parse_args():
    parser = ArgumentParser(description="Build a reference tree for EPA taxonomic placement.",
    epilog="Example: ./epa_trainer.py -t example/training_tax.txt -s example/training_seq.fa -n myref",
    formatter_class=RawTextHelpFormatter)
    parser.add_argument("-t", dest="taxonomy_fname", required=True,
            help="""Reference taxonomy file.""")
    parser.add_argument("-s", dest="align_fname", required=True,
            help="""Reference alignment file. Sequences must be aligned, their IDs must correspond to those
in taxonomy file.""")
    parser.add_argument("-r", dest="ref_fname",
            help="""Reference output file. It will contain reference alignment, phylogenetic tree and other
information needed for taxonomic placement of query sequences.""")
    parser.add_argument("-T", dest="num_threads", type=int, default=None,
            help="""Specify the number of CPUs.  Default: %d""" % multiprocessing.cpu_count())            
    parser.add_argument("-c", dest="config_fname", default=None,
            help="""Config file name.""")
    parser.add_argument("-o", dest="output_dir", default=".",
            help="""Output directory""")
    parser.add_argument("-n", dest="output_name", default=None,
            help="""Run name.""")
    parser.add_argument("-p", dest="rand_seed", type=int, default=None,
            help="""Random seed to be used with RAxML. Default: current system time.""")
    parser.add_argument("-m", dest="mfresolv_method", choices=["thorough", "fast", "ultrafast"],
            default="thorough", help="""Method of multifurcation resolution: 
            thorough    use stardard constrainted RAxML tree search (default)
            fast        use RF distance as search convergence criterion (RAxML -D option)
            ultrafast   optimize model+branch lengths only (RAxML -f e option)""")
    parser.add_argument("-N", dest="rep_num", type=int, default=1, 
            help="""Number of RAxML tree searches (with distinct random seeds). Default: 1""")
    parser.add_argument("-x", dest="taxcode_name", choices=["bac", "bot", "zoo", "vir"], type = str.lower,
            help="""Taxonomic code: BAC(teriological), BOT(anical), ZOO(logical), VIR(ological)""")
    parser.add_argument("-R", dest="restart", action="store_true",
            help="""Resume execution after a premature termination (e.g., due to expired job time limit).
Run name of the previous (terminated) job must be specified via -n option.""")
    parser.add_argument("-v", dest="verbose", action="store_true",
            help="""Print additional info messages to the console.""")
    parser.add_argument("-debug", dest="debug", action="store_true",
            help="""Debug mode, intermediate files will not be cleaned up.""")
    parser.add_argument("-no-hmmer", dest="no_hmmer", action="store_true",
            help="""Do not build HMMER profile.""")
    parser.add_argument("-dup-rank-names", dest="dup_rank_names", choices=["ignore", "abort", "skip", "autofix"],
            default="ignore", help="""Action to be performed if different ranks with same name are found: 
            ignore      do nothing
            abort       report duplicates and exit
            skip        skip the corresponding sequences (exlude from reference)
            autofix     make name unique by concatenating it with the parent rank's name""")
    parser.add_argument("-wrong-rank-count", dest="wrong_rank_count", choices=["ignore", "abort", "skip", "autofix"],
            default="ignore", help="""Action to be performed if lineage has less (more) than 7 ranks
            ignore      do nothing
            abort       report duplicates and exit
            skip        skip the corresponding sequences (exlude from reference)
            autofix     try to guess wich ranks should be added or removed (use with caution!)""")
    parser.add_argument("-tmpdir", dest="temp_dir", default=None,
            help="""Directory for temporary files.""")
    
    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()
    
    return args
 
def check_args(args):
    #check if taxonomy file exists
    if not os.path.isfile(args.taxonomy_fname):
        print "ERROR: Taxonomy file not found: %s" % args.taxonomy_fname
        sys.exit()

    #check if alignment file exists
    if not os.path.isfile(args.align_fname):
        print "ERROR: Alignment file not found: %s" % args.align_fname
        sys.exit()
        
    if not args.output_name:
        args.output_name = args.align_fname

    if not args.ref_fname:
        args.ref_fname = "%s.refjson" % args.output_name
    
    if args.output_dir and not os.path.dirname(args.ref_fname):
        args.ref_fname = os.path.join(args.output_dir, args.ref_fname)

    #check if reference json file already exists
    if os.path.isfile(args.ref_fname):
        print "ERROR: Reference tree file already exists: %s" % args.ref_fname
        print "Please delete it explicitely if you want to overwrite."
        sys.exit()
    
    #check if reference file can be created
    try:
        f = open(args.ref_fname, "w")
        f.close()
        os.remove(args.ref_fname)
    except:
        print "ERROR: cannot create output file: %s" % args.ref_fname
        print "Please check if directory %s exists and you have write permissions for it." % os.path.split(os.path.abspath(args.ref_fname))[0]
        sys.exit()
        
    if args.rep_num < 1 or args.rep_num > 1000:
        print "ERROR: Number of RAxML runs must be between 1 and 1000."
        sys.exit()

def which(program, custom_path=[]):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        path_list = custom_path
        path_list += os.environ["PATH"].split(os.pathsep)
        for path in path_list:
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None        
    
def check_dep(config):           
    if not config.no_hmmer:
        if not which("hmmalign", [config.hmmer_home]):
            print "ERROR: HMMER not found!"
            print "Please either specify path to HMMER executables in the config file" 
            print "or call this script with -no-hmmer option to skip building HMMER profile." 
            config.exit_user_error()
            
def run_trainer(config):
    check_dep(config)
    builder = RefTreeBuilder(config)
    builder.invocation_epac = " ".join(sys.argv)
    builder.build_ref_tree()
        
# -------
# MAIN
# -------
if __name__ == "__main__":
    args = parse_args()
    check_args(args)
    config = EpacTrainerConfig(args)

    print ""
    config.print_version("SATIVA-trainer")

    start_time = time.time()

    run_trainer(config)
    config.clean_tempdir()

    config.log.info("Reference JSON was saved to: %s", os.path.abspath(config.refjson_fname))
    config.log.info("Execution log was saved to: %s\n", os.path.abspath(config.log_fname))

    elapsed_time = time.time() - start_time
    config.log.info("Training completed successfully, elapsed time: %.0f seconds\n", elapsed_time)
