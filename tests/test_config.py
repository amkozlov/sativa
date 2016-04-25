#!/usr/bin/env python
import os
import sys
import unittest
from argparse import Namespace

lib_path = os.path.abspath('..')
sys.path.append(lib_path)

from epac.taxonomy_util import Taxonomy, TaxCode
from epac.config import EpacConfig, EpacTrainerConfig, EpacClassifierConfig, SativaConfig
from epac.ete2 import Tree

class ConfigTests(unittest.TestCase):

    def setUp(self):
        self.test_dir = os.path.dirname(os.path.abspath(__file__))
        self.testfile_dir = os.path.join(self.test_dir, "testfiles")
        self.cfg_fname = os.path.join(self.testfile_dir, "test.cfg")
        self.tmp_dir = os.path.join(self.test_dir, "../tmp")
    
    def get_default_namespace(self):
        ns = Namespace()
        ns.verbose = False
        ns.debug = False
        ns.restart = False
        ns.ref_fname = None
        ns.rand_seed = None
        ns.output_name = None
        ns.output_dir = self.testfile_dir
        ns.temp_dir = self.tmp_dir
        ns.config_fname = None
        ns.num_threads = None
        return ns
        
    def get_trainer_namespace(self):
        ns = self.get_default_namespace()
        ns.taxonomy_fname = None
        ns.align_fname = None
        ns.no_hmmer = None
        ns.dup_rank_names  = None
        ns.wrong_rank_count  = None
        ns.mfresolv_method = None
        ns.taxcode_name = None
        ns.rep_num = None
        ns.synonym_fname = None
        return ns
        
    def get_classifier_namespace(self):
        ns = self.get_default_namespace()
        ns.taxassign_method = None
        ns.min_lhw = None
        ns.brlen_pv = None
        return ns

    def get_sativa_namespace(self):
        ns = self.get_trainer_namespace()
        ns.taxassign_method = None
        ns.min_lhw = None
        ns.brlen_pv = None
        ns.ranktest = None
        ns.jplace_fname = None
        ns.final_jplace_fname = None
        ns.conf_cutoff = None
        ns.save_memory = None
        return ns

    def check_common_config(self, cfg):
        self.assertTrue(len(cfg.name) > 0)
        self.assertTrue(os.path.isdir(cfg.temp_dir))
        cfg.clean_tempdir()
        self.assertFalse(os.path.isdir(cfg.temp_dir))
        self.assertTrue(cfg.raxml_model.lower() == "auto")
        self.assertTrue(cfg.epa_use_heuristic.lower() == "auto")
        cfg.resolve_auto_settings(10)
        self.assertFalse(cfg.raxml_model == "auto")
        self.assertFalse(cfg.epa_use_heuristic.lower() == "auto")
        cfg.set_defaults()
        self.assertTrue(cfg.raxml_model.lower() == "auto")
        self.assertTrue(cfg.epa_use_heuristic.lower() == "auto")
        cfg.resolve_auto_settings(10e6)
        self.assertFalse(cfg.raxml_model == "auto")
        self.assertFalse(cfg.epa_use_heuristic.lower() == "auto")
        os.remove(cfg.log_fname)
    
    def test_epac_config(self):
        args = self.get_default_namespace()
        cfg = EpacConfig(args)
        self.check_common_config(cfg)

    def test_trainer_config(self):
        args = self.get_trainer_namespace()
        cfg = EpacTrainerConfig(args)
        self.check_common_config(cfg)

    def test_classifier_config(self):
        args = self.get_classifier_namespace()
        cfg = EpacClassifierConfig(args)
        self.check_common_config(cfg)

    def test_sativa_config(self):
        args = self.get_sativa_namespace()
        cfg = SativaConfig(args)
        self.check_common_config(cfg)
        
if __name__ == '__main__':
    unittest.main()
