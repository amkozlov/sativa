#!/usr/bin/env python
import os
import sys
import unittest

lib_path = os.path.abspath('..')
sys.path.append(lib_path)

from epac.ete2 import SeqGroup
from epac.taxonomy_util import Taxonomy, TaxCode
from epac.config import EpacConfig
from epa_trainer import InputValidator

class InputValidatorTests(unittest.TestCase):

    def setUp(self):
        cfg = EpacConfig()
        testfile_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "testfiles")
        tax_fname = os.path.join(testfile_dir, "test.tax")
        phy_fname = os.path.join(testfile_dir, "test.phy")
        tax = Taxonomy(EpacConfig.REF_SEQ_PREFIX, tax_fname)
        seqs = SeqGroup(sequences=phy_fname, format = "phylip")
        self.inval = InputValidator(cfg, tax, seqs, False)
    
    def test_seq_ids(self):
        expected_mis_ids = ["Missing1", "Missing2"]
        mis_ids = self.inval.check_seq_ids()
        self.assertEqual(set(mis_ids), set(expected_mis_ids))
        
    def test_identical_seqs(self):
        expected_dups = ["DupSeq01", "DupSeq02"]
        count, dups = self.inval.check_identical_seqs()
        self.assertEqual(len(dups), 1)
        self.assertEqual(set(dups[0]), set(expected_dups))
        
    def test_identical_ranks(self):
        expected_dups = ["DupSeq01", "DupSeq02"]
        expected_merges = [self.inval.taxonomy.seq_rank_id(sid) for sid in expected_dups]
        merged_ranks = self.inval.check_identical_ranks()
        self.assertEqual(len(merged_ranks), 1)
        self.assertEqual(merged_ranks.popitem()[1], expected_merges)
        
if __name__ == '__main__':
    unittest.main()
