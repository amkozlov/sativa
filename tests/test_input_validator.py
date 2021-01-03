#!/usr/bin/env python3
import os
import sys
import unittest

lib_path = os.path.abspath('..')
sys.path.append(lib_path)

from epac.ete2 import SeqGroup
from epac.taxonomy_util import Taxonomy, TaxCode
from epac.config import EpacTrainerConfig,EpacConfig
from epa_trainer import InputValidator

class InputValidatorTests(unittest.TestCase):

    def setUp(self):
        cfg = EpacTrainerConfig()
        cfg.debug=True
        testfile_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "testfiles")
        tax_fname = os.path.join(testfile_dir, "test.tax")
        phy_fname = os.path.join(testfile_dir, "test.phy")
        tax = Taxonomy(EpacConfig.REF_SEQ_PREFIX, tax_fname)
        seqs = SeqGroup(sequences=phy_fname, format = "phylip")
        self.inval = InputValidator(cfg, tax, seqs, False)

        self.expected_mis_ids = ["Missing1", "Missing2"]
        self.expected_dups = ["DupSeq(01)", "DupSeq02"]
        self.expected_merges = [self.inval.taxonomy.seq_rank_id(sid) for sid in self.expected_dups]
    
    def test_seq_ids(self):
        mis_ids = self.inval.check_seq_ids()
        self.assertEqual(set(mis_ids), set(self.expected_mis_ids))
        
    def test_identical_seqs(self):
        count, dups = self.inval.check_identical_seqs()
        self.assertEqual(len(dups), 1)
        self.assertEqual(set(dups[0]), set(self.expected_dups))
        
    def test_identical_ranks(self):
        merged_ranks = self.inval.check_identical_ranks()
        self.assertEqual(len(merged_ranks), 1)
        self.assertEqual(set(merged_ranks.popitem()[1]), set(self.expected_merges))

    def test_validate(self):
        corr_ranks, corr_seqid, merged_ranks = self.inval.validate()
        self.assertEqual(len(merged_ranks), 1)
        merged_list = merged_ranks.popitem()[1]
        corr_merged_list = [self.inval.taxonomy.get_uncorr_rank_id(i) for i in merged_list]
        self.assertEqual(set(corr_merged_list), set(self.expected_merges))
        
if __name__ == '__main__':
    unittest.main()
