#!/usr/bin/env python

import os
import sys
import glob
import shutil
import datetime
import time
import logging
import random
import multiprocessing
import ConfigParser

from epac.version import SATIVA_BUILD,SATIVA_RELEASE_DATE,SATIVA_RAXML_VER

class DefaultedConfigParser(ConfigParser.SafeConfigParser):
    def get_param(self, section, option, ctype=str, default=None):
        if default is None:
            ret = self.get(section, option)
        else:
            confdict = self.__dict__.get('_sections')
            sec = confdict.get(section)
            if sec:
                if ctype == bool:
                    ret = self.getboolean(section, option)
                else:
                    ret = sec.get(option, default)
            else:
                ret = default
        return ctype(ret)

class EpacConfig:
    # this prefix will be added to every sequence name in reference to prevent 
    # name clashes with query sequences, which are coded with numbers
    REF_SEQ_PREFIX = "r_";
    QUERY_SEQ_PREFIX = "q_";

    CAT_LOWER_THRES   = 100
    CAT_GAMMA_THRES   = 500
    GAMMA_UPPER_THRES = 10000
    EPA_HEUR_THRES    = 1000
    
    SATIVA_INFO = \
    """%s %s, released on %s. Last version: https://github.com/amkozlov/sativa 
By A.Kozlov and J.Zhang, the Exelixis Lab. Based on RAxML %s by A.Stamatakis.\n"""\
    % ("%s", SATIVA_BUILD, SATIVA_RELEASE_DATE, SATIVA_RAXML_VER)
    
    
    @staticmethod
    def strip_prefix(seq_name, prefix):
        if seq_name.startswith(prefix):
            plen = len(prefix)
            return seq_name[plen:]
        else:
            return seq_name
        
    @staticmethod
    def strip_ref_prefix(seq_name):
        return EpacConfig.strip_prefix(seq_name, EpacConfig.REF_SEQ_PREFIX)
        
    @staticmethod
    def strip_query_prefix(seq_name):
        return EpacConfig.strip_prefix(seq_name, EpacConfig.QUERY_SEQ_PREFIX)

    def __init__(self, args=None): 
        self.basepath = os.path.dirname(os.path.abspath(__file__))
        self.epac_home = os.path.abspath(os.path.join(self.basepath, os.pardir)) + "/"

        self.set_defaults()

        if not args:
            return
        
        self.verbose = args.verbose
        self.debug = args.debug
        self.restart = args.restart
        self.refjson_fname = args.ref_fname 
               
        self.rand_seed = args.rand_seed 
        
        timestamp = "%d" % (time.time()*1000) #datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        if args.output_name:
            self.name = args.output_name
        else:
            self.name = str(random.randint(1, 99999))

        self.output_dir = args.output_dir

        if args.temp_dir:
            self.temp_dir = args.temp_dir
        else:
            self.temp_dir = os.path.join(self.epac_home, "tmp")

        if self.restart:
            tmpdirs = glob.glob(os.path.join(self.temp_dir, self.name + "_*"))
            if len(tmpdirs) > 0:
                tmpdirs.sort(key=os.path.getmtime, reverse=True)
                self.temp_dir = tmpdirs[0]
            else:
                self.exit_user_error("ERROR: Cannot resume execution: temp directory for the previous run not found in %s" % self.temp_dir)
        else:
            temp_name = self.name + "_" + timestamp
            self.temp_dir = os.path.join(self.temp_dir, temp_name)
            
        self.raxml_outdir = self.temp_dir
        self.raxml_outdir_abs = os.path.abspath(self.raxml_outdir)

        self.init_logger()
        
        if not args.config_fname:
            args.config_fname = os.path.join(self.epac_home, "sativa.cfg")
        
        self.config_path = os.path.dirname(os.path.abspath(args.config_fname))
        self.read_from_file(args.config_fname)

        # command line setting has preference over config file and default
        if args.num_threads:
            self.num_threads = args.num_threads        
        self.check_raxml()    
        
        if not self.restart:
            os.mkdir(self.temp_dir)

    def set_defaults(self):
        self.muscle_home = self.epac_home + "/epac/bin" + "/"
        self.hmmer_home = self.epac_home + "/epac/bin" + "/"
        self.raxml_home = self.epac_home + "/epac/bin" + "/"
        self.raxml_exec = "raxmlHPC-PTHREADS-SSE3"
        self.raxml_model = "AUTO"
        self.raxml_remote_host = ""
        self.raxml_remote_call = False        
        self.run_on_cluster = False
        self.cluster_epac_home = self.epac_home
        self.cluster_qsub_script = ""
        self.epa_load_optmod = True
        self.epa_use_heuristic = "AUTO"
        self.epa_heur_rate = 0.01
        self.min_confidence = 0.2
        self.num_threads = multiprocessing.cpu_count()
        self.compress_patterns = False
        self.use_bfgs = False
        self.save_memory = False
        self.taxa_ident_thres = 0.6
        self.debug = False
        self.restart = False
        self.verbose = False
        self.log = logging.getLogger('epac')
       
    def init_logger(self):
        self.log_fname = self.out_fname("%NAME%.log")
        if self.verbose or self.debug:
           log_lvl = logging.DEBUG
        else:
           log_lvl = logging.INFO

        # configure logger object
        self.log.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(message)s')

        # add console handler
        ch = logging.StreamHandler()
        ch.setLevel(log_lvl)
        ch.setFormatter(formatter)        
        self.log.addHandler(ch)

        # add file handler
        if self.restart:
            logf_mode = "a"
        else:
            logf_mode = "w"
        fh = logging.FileHandler(self.log_fname, mode=logf_mode)
        fh.setLevel(log_lvl)
        fh.setFormatter(formatter)        
        self.log.addHandler(fh)

    def resolve_auto_settings(self, tree_size):
        if self.raxml_model == "AUTO":
            if tree_size > EpacConfig.CAT_GAMMA_THRES:
                self.raxml_model = "GTRCAT"
            else:
                self.raxml_model = "GTRGAMMA"
        elif self.raxml_model == "GTRCAT" and tree_size < EpacConfig.CAT_LOWER_THRES:
            print "WARNING: You're using GTRCAT model on a very small dataset (%d taxa), which might lead to unreliable results!" % tree_size
            print "Please consider switching to GTRGAMMA model.\n"
        elif self.raxml_model == "GTRGAMMA" and tree_size > EpacConfig.GAMMA_UPPER_THRES:
            print "WARNING: You're using GTRGAMMA model on a very large dataset (%d taxa), which might lead to numerical issues!" % tree_size
            print "In case of problems, please consider switching to GTRCAT model.\n"

        if self.epa_use_heuristic == "AUTO": 
            if tree_size > EpacConfig.EPA_HEUR_THRES:
                self.epa_use_heuristic = "TRUE"
                self.epa_heur_rate = 0.5 * float(EpacConfig.EPA_HEUR_THRES) / tree_size
            else:
                self.epa_use_heuristic = "FALSE"

    def resolve_relative_path(self, rpath):
        if rpath.startswith("/"):
            return rpath
        else:
            return os.path.join(self.config_path, rpath)
        
    def check_raxml(self):
        self.raxml_exec_full = self.raxml_home + self.raxml_exec
        if self.raxml_remote_host in ["", "localhost"]:
            self.raxml_remote_call = False
            # if raxml_home is empty, raxml binary must be on PATH; otherwise check if file exists
            if self.raxml_home:  
                if not os.path.isdir(self.raxml_home):
                    self.exit_user_error("RAxML home directory not found: %s" % self.raxml_home)
                elif not os.path.isfile(self.raxml_exec_full):
                    self.exit_user_error("RAxML executable not found: %s" % self.raxml_exec_full)
        else:
            self.raxml_remote_call = True
        self.raxml_cmd = [self.raxml_exec_full, "-w", self.raxml_outdir_abs]
        if self.num_threads > 1:
            self.raxml_cmd += ["-T", str(self.num_threads)]
        
    def read_from_file(self, config_fname):
        if not os.path.exists(config_fname):
            self.exit_user_error("ERROR: Config file not found: %s" % config_fname)

        parser = DefaultedConfigParser() #ConfigParser.SafeConfigParser()
        parser.read(config_fname)
        
        self.raxml_home = parser.get_param("raxml", "raxml_home", str, self.raxml_home)
        if self.raxml_home:
            self.raxml_home = self.resolve_relative_path(self.raxml_home + "/")
        self.raxml_exec = parser.get_param("raxml", "raxml_exec", str, self.raxml_exec)
        self.raxml_remote_host = parser.get_param("raxml", "raxml_remote_host", str, self.raxml_remote_host)

        self.raxml_model = parser.get_param("raxml", "raxml_model", str, self.raxml_model).upper()
        self.num_threads = parser.get_param("raxml", "raxml_threads", int, self.num_threads)

        self.epa_use_heuristic = parser.get_param("raxml", "epa_use_heuristic", str, self.epa_use_heuristic).upper()
        self.epa_heur_rate = parser.get_param("raxml", "epa_heur_rate", float, self.epa_heur_rate)
        self.epa_load_optmod = parser.get_param("raxml", "epa_load_optmod", bool, self.epa_load_optmod)

        self.hmmer_home = self.resolve_relative_path(parser.get_param("hmmer", "hmmer_home", str, self.hmmer_home))
        self.muscle_home = self.resolve_relative_path(parser.get_param("muscle", "muscle_home", str, self.muscle_home))
        
        self.run_on_cluster = parser.get_param("cluster", "run_on_cluster", bool, self.run_on_cluster)
        self.cluster_epac_home = parser.get_param("cluster", "cluster_epac_home", str, self.cluster_epac_home) + "/"
        self.cluster_qsub_script = parser.get_param("cluster", "cluster_qsub_script", str, self.cluster_qsub_script)

        self.min_confidence = parser.get_param("assignment", "min_confidence", float, self.min_confidence)

        return parser

    def subst_name(self, in_str):
        """Replace %NAME% macros with an actual EPAC run name. Used to 
        generate unique run-specific identifiers (filenames, RAxML job names etc)"""
        return in_str.replace("%NAME%", self.name)
    
    def tmp_fname(self, fname):
        return os.path.join(self.temp_dir, self.subst_name(fname))

    def out_fname(self, fname):
        return os.path.join(self.output_dir, self.subst_name(fname))
        
    def clean_tempdir(self):
        if not self.debug and os.path.isdir(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def exit_fatal_error(self, msg=None):
        if msg:
            self.log.error(msg)
        sys.exit(13)    

    def exit_user_error(self, msg=None):
        if msg:
            self.log.error(msg)
        self.clean_tempdir()
        sys.exit(14)

    def print_version(self, progname):
        self.log.info(EpacConfig.SATIVA_INFO % progname)    
        if self.restart:
            self.log.info("Resuming %s execution using files from previous run found in: %s\n" % (progname, self.temp_dir))
    
class EpacTrainerConfig(EpacConfig):
    
    def __init__(self, args=None):
        EpacConfig.__init__(self, args)
        if args:
            self.taxonomy_fname = args.taxonomy_fname
            self.align_fname = args.align_fname
            self.no_hmmer = args.no_hmmer
            self.dup_rank_names  = args.dup_rank_names
            self.wrong_rank_count  = args.wrong_rank_count
            self.mfresolv_method = args.mfresolv_method
            self.taxcode_name = args.taxcode_name
            self.rep_num = args.rep_num
        
    def set_defaults(self):
        EpacConfig.set_defaults(self)
        self.no_hmmer = False
        # whether model parameters should be re-optimized from scratch on the best topology using "-f e"
        self.reopt_model = False
        self.dup_rank_names = "ignore"
        self.wrong_rank_count = "ignore"
        self.taxassign_method = "1"
        self.mfresolv_method = "thorough"
        self.taxcode_name = "bac"
        self.rep_num = 1
        # default settings below imply no taxonomy filtering, 
        # i.e. all sequences from taxonomy file will be included into reference tree
        self.reftree_min_rank = 0
        self.reftree_max_seqs_per_leaf = 1e6
        self.reftree_clades_to_include=[]
        self.reftree_clades_to_ignore=[]
        

    def read_from_file(self, config_fname):
        parser = EpacConfig.read_from_file(self, config_fname)
        
        self.reftree_min_rank = parser.get_param("reftree", "min_rank", int, self.reftree_min_rank)
        self.reftree_max_seqs_per_leaf = parser.get_param("reftree", "max_seqs_per_leaf", int, self.reftree_max_seqs_per_leaf)
        clades_str = parser.get_param("reftree", "clades_to_include", str, "")
        self.reftree_clades_to_include = self.parse_clades(clades_str)
        clades_str = parser.get_param("reftree", "clades_to_ignore", str, "")
        self.reftree_clades_to_ignore = self.parse_clades(clades_str)
        
    def parse_clades(self, clades_str):
        clade_list = []
        try:        
            if clades_str:
                clades = clades_str.split(",")
                for clade in clades:
                    toks = clade.split("|")
                    clade_list += [(int(toks[0]), toks[1])]
        except:
            print "Invalid format in config parameter: clades_to_include"
            sys.exit()

        return clade_list
        
class EpacClassifierConfig(EpacConfig):

    def __init__(self, args=None):
        EpacConfig.__init__(self, args)

        if args:
            self.taxassign_method = args.taxassign_method
            self.min_lhw = args.min_lhw
            self.brlen_pv = args.brlen_pv
        else:
            self.taxassign_method = "1"
            self.min_lhw = 0.
            self.brlen_pv = 0.

class SativaConfig(EpacTrainerConfig):
    
    def __init__(self, args):
        args.no_hmmer = True
        args.dup_rank_names = "ignore"
        args.wrong_rank_count = "ignore"
        args.taxassign_method = "1"

        EpacTrainerConfig.__init__(self, args)

        self.taxassign_method = args.taxassign_method
        self.min_lhw = args.min_lhw
        self.brlen_pv = args.brlen_pv
        self.ranktest = args.ranktest
        self.conf_cutoff = args.conf_cutoff
        self.jplace_fname = args.jplace_fname
        
        self.output_interim_files = True
        self.compress_patterns = True
        self.use_bfgs = True
        self.save_memory = False

        if self.refjson_fname:
            self.load_refjson = True
        else:
            self.load_refjson = False
            if self.output_interim_files:
                self.refjson_fname = self.out_fname("%NAME%.refjson")
            else:
                self.refjson_fname = self.tmp_fname("%NAME%.refjson")
                
        if self.restart and os.path.isfile(self.refjson_fname):
            self.load_refjson = True
        
    def set_defaults(self):
        EpacTrainerConfig.set_defaults(self)

    def read_from_file(self, config_fname):
        parser = EpacTrainerConfig.read_from_file(self, config_fname)
