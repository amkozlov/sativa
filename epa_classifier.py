#! /usr/bin/env python
try:
    import sys
    import os
    import time
    import glob
    import multiprocessing
    import logging    
    
    from epac.ete2 import Tree, SeqGroup
    from epac.argparse import ArgumentParser
    from epac.config import EpacConfig
    from epac.raxml_util import RaxmlWrapper, FileUtils
    from epac.json_util import RefJsonParser, RefJsonChecker, EpaJsonParser
    from epac.msa import muscle, hmmer
    from epac.taxonomy_util import Taxonomy
    from epac.classify_util import TaxClassifyHelper
except ImportError, e:
    print("Some packages are missing, please re-downloand EPA-classifier")
    print e
    sys.exit()

class EpaClassifier:
    def __init__(self, config, args):
        self.cfg = config
        self.jplace_fname = args.jplace_fname
        self.ignore_refalign = args.ignore_refalign
        
        self.tmp_refaln = config.tmp_fname("%NAME%.refaln")
        #here is the final alignment file for running EPA
        self.epa_alignment = config.tmp_fname("%NAME%.afa")
        self.hmmprofile = config.tmp_fname("%NAME%.hmmprofile")
        self.tmpquery = config.tmp_fname("%NAME%.tmpquery")
        self.noalign = config.tmp_fname("%NAME%.noalign")
        self.seqs = None
        
        assign_fname = args.output_name + ".assignment.txt"
        self.out_assign_fname = os.path.join(args.output_dir, assign_fname)
        jplace_fname = args.output_name + ".jplace"
        self.out_jplace_fname = os.path.join(args.output_dir, jplace_fname)

        try:
            self.refjson = RefJsonParser(config.refjson_fname)
        except ValueError:
            self.cfg.exit_user_error("Invalid json file format: %s" % config.refjson_fname)
        #validate input json format 
        self.refjson.validate()
        self.bid_taxonomy_map = self.refjson.get_bid_tanomomy_map()
        self.reftree = self.refjson.get_reftree()
        self.rate = self.refjson.get_rate()
        self.node_height = self.refjson.get_node_height()
        self.cfg.compress_patterns = self.refjson.get_pattern_compression()

        self.classify_helper = TaxClassifyHelper(self.cfg, self.bid_taxonomy_map, args.p_value, self.rate, self.node_height)
        
    def require_muscle(self):
        basepath = os.path.dirname(os.path.abspath(__file__))
        if not os.path.exists(basepath + "/epac/bin/muscle"):
            errmsg = "The pipeline uses MUSCLE to merge alignments, please download the programm from:\n" + \
                     "http://www.drive5.com/muscle/downloads.htm\n" + \
                     "and specify path to your installation in the config file (sativa.cfg)\n"
            self.cfg.exit_user_error(errmsg)

    def require_hmmer(self):
        basepath = os.path.dirname(os.path.abspath(__file__))
        if not os.path.exists(basepath + "/epac/bin/hmmbuild") or not os.path.exists(basepath + "/epac/bin/hmmalign"):
            errmsg = "The pipeline uses HAMMER to align the query seqeunces, please download the programm from:\n" + \
                     "http://hmmer.janelia.org/\n" + \
                     "and specify path to your installation in the config file (sativa.cfg)\n"
            self.cfg.exit_user_error(errmsg)

    def align_to_refenence(self, noalign, minp = 0.9):
        refaln = self.refjson.get_alignment(fout = self.tmp_refaln)
        fprofile = self.refjson.get_hmm_profile(self.hmmprofile)
        
        # if there is no hmmer profile in json file, build it from scratch          
        if not fprofile:
            hmm = hmmer(self.cfg, refaln)
            fprofile = hmm.build_hmm_profile()
    
        hm = hmmer(config = self.cfg, refalign = refaln , query = self.tmpquery, refprofile = fprofile, discard = noalign, seqs = self.seqs, minp = minp)
        self.epa_alignment = hm.align()

    def merge_alignment(self, query_seqs):
        refaln = self.refjson.get_alignment_list()
        with open(self.epa_alignment, "w") as fout:
            for seq in refaln:
                fout.write(">" + seq[0] + "\n" + seq[1] + "\n")
            for name, seq, comment, sid in query_seqs.iter_entries():
                fout.write(">" + name + "\n" + seq + "\n")


    def checkinput(self, query_fname, minp = 0.9):
        formats = ["fasta", "phylip", "iphylip", "phylip_relaxed", "iphylip_relaxed"]
        for fmt in formats:
            try:
                self.seqs = SeqGroup(sequences=query_fname, format = fmt)
                break
            except:
                self.cfg.log.debug("Guessing input format: not " + fmt)
        if self.seqs == None:
            self.cfg.exit_user_error("Invalid input file format: %s\nThe supported input formats are fasta and phylip" % query_fname)

        if self.ignore_refalign:
            self.cfg.log.info("Assuming query file contains reference sequences, skipping the alignment step...\n")
            with open(self.epa_alignment, "w") as fout:
                for name, seq, comment, sid in self.seqs.iter_entries():
                    ref_name = EpacConfig.REF_SEQ_PREFIX + name
                    if ref_name in self.refjson.get_sequences_names():
                        seq_name = ref_name
                    else:
                        seq_name = EpacConfig.QUERY_SEQ_PREFIX + name
                    fout.write(">" + seq_name + "\n" + seq + "\n")
            return
            
        # add query seq name prefix to avoid confusion between reference and query sequences
        self.seqs.add_name_prefix(EpacConfig.QUERY_SEQ_PREFIX)
        
        self.seqs.write(format="fasta", outfile=self.tmpquery)
        self.cfg.log.info("Checking if query sequences are aligned ...")
        entries = self.seqs.get_entries()
        seql = len(entries[0][1])
        aligned = True
        for entri in entries[1:]:
            l = len(entri[1])
            if not seql == l:
                aligned = False
                break
        
        if aligned and len(self.seqs) > 1:
            self.cfg.log.info("Query sequences are aligned")
            refalnl = self.refjson.get_alignment_length()
            if refalnl == seql:
                self.cfg.log.info("Merging query alignment with reference alignment")
                self.merge_alignment(self.seqs)
            else:
                self.cfg.log.info("Merging query alignment with reference alignment using MUSCLE")
                self.require_muscle()
                refaln = self.refjson.get_alignment(fout = self.tmp_refaln)
                m = muscle(self.cfg)
                self.epa_alignment = m.merge(refaln, self.tmpquery)
        else:
            self.cfg.log.info("Query sequences are not aligned")
            self.cfg.log.info("Align query sequences to the reference alignment using HMMER")
            self.require_hmmer()
            self.align_to_refenence(self.noalign, minp = minp)

    def print_ranks(self, rks, confs, minlw = 0.0):
        ss = ""
        css = ""
        for i in range(len(rks)):
            conf = confs[i]
            if conf == confs[0] and confs[0] >=0.99:
                conf = 1.0
            if conf >= minlw:
                ss = ss + rks[i] + ";"
                css = css + "{0:.3f}".format(conf) + ";"
            else:
                break
        if ss == "":
            return None
        else:
            return ss[:-1] + "\t" + css[:-1]


    def classify(self, query_fname, method = "1", minlw = 0.0, pv = 0.02, minp = 0.9, ptp = False):
        if self.jplace_fname:
            jp = EpaJsonParser(self.jplace_fname)
        else:        
            self.checkinput(query_fname, minp)

            self.cfg.log.info("Running RAxML-EPA...\n")
            raxml = RaxmlWrapper(config)
            reftree_fname = self.cfg.tmp_fname("ref_%NAME%.tre")
            self.refjson.get_raxml_readable_tree(reftree_fname)
            optmod_fname = self.cfg.tmp_fname("%NAME%.opt")
            self.refjson.get_binary_model(optmod_fname)
            job_name = self.cfg.subst_name("epa_%NAME%")

            reftree_str = self.refjson.get_raxml_readable_tree()
            reftree = Tree(reftree_str)

            self.reftree_size = len(reftree.get_leaves())

            # IMPORTANT: set EPA heuristic rate based on tree size!                
            self.cfg.resolve_auto_settings(self.reftree_size)
            # If we're loading the pre-optimized model, we MUST set the same rate het. mode as in the ref file        
            if self.cfg.epa_load_optmod:
                self.cfg.raxml_model = self.refjson.get_ratehet_model()

            reduced_align_fname = raxml.reduce_alignment(self.epa_alignment)

            jp = raxml.run_epa(job_name, reduced_align_fname, reftree_fname, optmod_fname)
            
            raxml.copy_epa_jplace(job_name, self.out_jplace_fname, move=True)
        
        self.cfg.log.info("Assigning taxonomic labels based on EPA placements...\n")
 
        placements = jp.get_placement()
        
        if self.out_assign_fname:
            fo = open(self.out_assign_fname, "w")
        else:
            fo = None
        
        output2 = ""
        for place in placements:
            output = None
            taxon_name = place["n"][0]
            origin_taxon_name = EpacConfig.strip_query_prefix(taxon_name)
            edges = place["p"]
#            edges = self.erlang_filter(edges, p = pv)
            if len(edges) > 0:
                ranks, lws = self.classify_helper.classify_seq(edges, method, minlw)
                
                isnovo = self.novelty_check(place_edge = str(edges[0][0]), ranks =ranks, lws = lws, minlw = minlw)
                rankout = self.print_ranks(ranks, lws, minlw)
                
                if rankout == None:
                    output2 = output2 + origin_taxon_name+ "\t\t\t?\n"
                else:
                    output = "%s\t%s\t" % (origin_taxon_name, self.print_ranks(ranks, lws, minlw))
                    if isnovo: 
                        output += "*"
                    else:
                        output +="o"
                    if self.cfg.verbose:
                        print(output) 
                    if fo:
                        fo.write(output + "\n")
            else:
                output2 = output2 + origin_taxon_name+ "\t\t\t?\n"
        
        if os.path.exists(self.noalign):
            with open(self.noalign) as fnoa:
                lines = fnoa.readlines()
                for line in lines:
                    taxon_name = line.strip()[1:]
                    origin_taxon_name = EpacConfig.strip_query_prefix(taxon_name)
                    output = "%s\t\t\t?" % origin_taxon_name
                    if self.cfg.verbose:
                        print(output)
                    if fo:
                        fo.write(output + "\n")
        
        if self.cfg.verbose:
            print(output2)
        
        if fo:
            fo.write(output2)
            fo.close()

        #############################################
        #
        # EPA-PTP species delimitation
        #
        #############################################
        if ptp:
            full_aln = SeqGroup(self.epa_alignment)
            species_list = epa_2_ptp(epa_jp = jp, ref_jp = self.refjson, full_alignment = full_aln, min_lw = 0.5, debug = self.cfg.debug)
            
            self.cfg.log.debug("Species clusters:")
 
            if fout:
                fo2 = open(fout+".species", "w")
            else:
                fo2 = None

            for sp_cluster in species_list:
                translated_taxa = []
                for taxon in sp_cluster:
                    origin_taxon_name = EpacConfig.strip_query_prefix(taxon)
                    translated_taxa.append(origin_taxon_name)
                s = ",".join(translated_taxa)
                if fo2:
                    fo2.write(s + "\n")
                self.cfg.log.debug(s)

            if fo2:
                fo2.close()
        #############################################
        
    def novelty_check(self, place_edge, ranks, lws, minlw):
        """If the taxonomic assignment is not assigned to the genus level, 
        we need to check if it is due to the incomplete reference taxonomy or 
        it is likely to be something new:
        
        1. If the final ranks are assinged because of lw cut, that means with samller lw
        the ranks can be further assinged to lowers. This indicate the undetermined ranks 
        in the assignment is not due to the incomplete reference taxonomy, so the query 
        sequence is likely to be something new.
        
        2. Otherwise We check all leaf nodes' immediate lower rank below this ml placement point, 
        if they are not empty, output all ranks and indicate this could be novelty.
        """
        
        lowrank = 0
        for i in range(len(ranks)):
            if i < 6:
                """above genus level"""
                rk = ranks[i]
                lw = lws[i]
                if rk == "-":
                    break
                else:
                    lowrank = lowrank + 1
                    if lw >=0 and lw < minlw:
                        return True
        
        if lowrank >= 5 and not ranks[lowrank] == "-":
            return False
        else:
            placenode = self.reftree.search_nodes(B = place_edge)[0]
            if placenode.is_leaf():
                return False
            else:
                leafnodes = placenode.get_leaves()
                flag = True
                for leaf in leafnodes:
                    br_num = leaf.B
                    branks = self.bid_taxonomy_map[br_num]
                    if branks[lowrank] == "-":
                        flag = False
                        break
                        
                return flag

def print_options():
    print("usage: python epa_classifier.py -r example/reference.json -q example/query.fa -t 0.5 -v")
    print("Options:")
    print("    -r reference                   Specify the reference alignment and  taxonomy  in  json  format.\n")
    print("    -q query sequence              Specify the query seqeunces file.")
    print("                                   If the query sequences are aligned to  the  reference  alignment")
    print("                                   already, the epa_classifier will  classify the  queries  to  the")
    print("                                   lowest rank possible. If the  query sequences  are aligned,  but")
    print("                                   not to the reference, then a profile alignment will be perfermed")
    print("                                   to merge the two alignments. If  the  query  sequences  are  not")
    print("                                   aligned, then HMMER will be used to align  the  queries  to  the")
    print("                                   reference alignment. \n")
    print("    -t min likelihood weight       A value between 0 and 1, the minimal sum of likelihood weight of")
    print("                                   an assignment to a specific rank. This value represents a confi-")
    print("                                   dence measure of the assignment,  assignments  below  this value")
    print("                                   will be discarded. Default: 0 to output all possbile assignments.\n")
    print("    -o outputfile                  Specify the file name for output.\n")
    print("    -p p-value                     P-value for branch length Erlang test. Default: 0.02\n")
    print("    -minalign min-aligned%         Minimal percent of aligned sites.  Default: 0.9 (90%)\n")
    print("    -m method                      Assignment method 1 or 2")
    print("                                   1: Max sum likelihood (default)")
    print("                                   2: Max likelihood placement\n ")
    print("    -T numthread                   Specify the number of CPUs.\n")
    print("    -v                             Print the results on screen.\n")
    print("    --ptp                          Delimit species with PTP.\n")


def parse_args():
    parser = ArgumentParser(description="Assign taxonomy ranks to query sequences.")
    parser.add_argument("-r", dest="ref_fname",
            help="""Specify the reference alignment and taxonomy in json format.""")
    parser.add_argument("-q", dest="query_fname",
            help="""Specify the query seqeunces file.\n
                    If the query sequences are aligned to the reference alignment
                    already, the epa_classifier will classify the queries to the
                    lowest rank possible. If the query sequences are aligned, but
                    not to the reference, then a profile alignment will be performed
                    to merge the two alignments. If the query sequences are not
                    aligned, then HMMER will be used to align the queries to the
                    reference alignment.""")
    parser.add_argument("-t", dest="min_lhw", type=float, default=0.,
            help="""A value between 0 and 1, the minimal sum of likelihood weight of
                    an assignment to a specific rank. This value represents a confidence 
                    measure of the assignment, assignments below this value will be discarded. 
                    Default: 0 to output all possbile assignments.""")
    parser.add_argument("-o", dest="output_dir", default=None,
            help="""Directory for result files  (default: same as query file directory)""")
    parser.add_argument("-n", dest="output_name", default=None,
            help="""Run name, will be used to name result files.""")
    parser.add_argument("-p", dest="p_value", type=float, default=0.0,
            help="""P-value for branch length Erlang test. Default: 0 (filter off)\n""")
    parser.add_argument("-minalign", dest="minalign", type=float, default=0.9,
            help="""Minimal percent of the sites aligned to the reference alignment.  Default: 0.9""")
    parser.add_argument("-m", dest="method", default="1",
            help="""Assignment method 1 or 2
                    1: Max sum of likelihood weights (default)
                    2: Max likelihood weight placement""")
    parser.add_argument("-T", dest="num_threads", type=int, default=None,
            help="""Specify the number of CPUs.  Default: %d""" % multiprocessing.cpu_count())
    parser.add_argument("-v", dest="verbose", action="store_true",
            help="""Print classification results and additional info messages to the console.""")
    parser.add_argument("-debug", dest="debug", action="store_true",
            help="""Debug mode, intermediate files will not be cleaned up.""")
    parser.add_argument("-a", dest="align_only", action="store_true",
            help="""Alignment only: Do not perform classification, just build the combined alignment (RS+QS) (default: OFF)""")
    parser.add_argument("-j", dest="jplace_fname",
            help="""Do not call RAxML EPA, use existing .jplace file as input instead.""")
    parser.add_argument("-c", dest="config_fname", default=None,
            help="Config file name.")
    parser.add_argument("--ptp", 
            help = """Delimit species with ptp""",
            default = False,
            action="store_true")
    parser.add_argument("-x", dest="ignore_refalign", action="store_true",
            help="Query file contains complete alignment (query+reference), so use it and ignore reference alignment from .json file.")
    parser.add_argument("-tmpdir", dest="temp_dir", default=None,
            help="""Directory for temporary files.""")
    args = parser.parse_args()
    return args


def check_args(args):    
    if not args.ref_fname:
        print("Must specify the reference in json format!\n")
        print_options()
        sys.exit()
    
    if not os.path.exists(args.ref_fname):
        print("Input reference json file does not exists: %s" % args.ref_fname)
        print_options()
        sys.exit()
    
    if args.jplace_fname:
        if os.path.exists(args.jplace_fname):
            input_fname = args.jplace_fname 
        else: 
            print("Portable tree file does not exist: %s" % args.jplace_fname)
            sys.exit()
    elif args.query_fname:
        if os.path.exists(args.query_fname):
            input_fname = args.query_fname 
        else:
            print("Input query file does not exists: %s" % args.query_fname)
            sys.exit()
    else:
        print("Either query file or .jplace file must be specified!\n")
        print_options()
        sys.exit()
    
    query_dir, query_fname = os.path.split(os.path.abspath(input_fname))
    if not args.output_dir:
        args.output_dir = query_dir
        
    if not args.output_name:
        args.output_name = query_fname
        
    if args.min_lhw < 0 or args.min_lhw > 1.0:
         args.min_lhw = 0.0
    
    if args.p_value < 0:
        args.p_value = 0
    
    if not (args.method == "1" or args.method == "2"):
        args.method == "1"


def print_run_info(config, args):
    print ""
    config.print_version("SATIVA-classifier")

    config.log.info("\nEPA-classifier running with the following parameters:")
    config.log.info(" Reference:......................%s" % os.path.abspath(args.ref_fname))
    config.log.info(" Query:..........................%s" % os.path.abspath(args.query_fname))
    config.log.info(" Min percent of alignment sites..%s" % args.minalign)
    config.log.info(" Min likelihood weight:..........%f" % args.min_lhw)
    config.log.info(" Assignment method:..............%s" % args.method)
    config.log.info(" P-value for Erlang test:........%f" % args.p_value)
    config.log.info(" Number of threads:..............%d" % config.num_threads)
#    print("Result will be written to:")
#    print(args.output_fname)
    config.log.info("")

    if config.debug:
        config.log.debug("Running in DEBUG mode, temp files will be saved to: %s\n", os.path.abspath(config.temp_dir))


if __name__ == "__main__":
    if len(sys.argv) == 1: 
        sys.argv.append("-h")
    args = parse_args()
    
    if args.ptp:
        try:
            from epac.epa_util import epa_2_ptp
        except ImportError, e:
            print("PTP module (or its dependecies) not found, please see detail below:")
            print e
            sys.exit()
    
    check_args(args)
    config = EpacConfig(args)
    print_run_info(config, args)
    
    start_time = time.time()
   
    ec = EpaClassifier(config, args)
    ec.classify(query_fname = args.query_fname, method = args.method, minlw = args.min_lhw, pv = args.p_value, minp = args.minalign, ptp = args.ptp)
    
    config.clean_tempdir()
    
    config.log.info("Taxonomic assignment results were saved to: %s", os.path.abspath(ec.out_assign_fname))
    config.log.info("Execution log was saved to: %s\n", os.path.abspath(config.log_fname))

    elapsed_time = time.time() - start_time
    config.log.info("Classification completed successfully, elapsed time: %.0f seconds\n", elapsed_time)

