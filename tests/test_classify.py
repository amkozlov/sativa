#!/usr/bin/env python3
import os
import sys
import unittest

lib_path = os.path.abspath('..')
sys.path.append(lib_path)

from epac.taxonomy_util import Taxonomy, TaxCode, TaxTreeBuilder
from epac.config import EpacConfig, EpacClassifierConfig
from epac.ete2 import Tree
from epac.classify_util import TaxTreeHelper, TaxClassifyHelper
from epac.json_util import EpaJsonParser

class TaxTreeHelperTests(unittest.TestCase):

    def setUp(self):
        self.testfile_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "testfiles")
        self.tax_fname = os.path.join(self.testfile_dir, "test_clean.tax")
        self.taxonomy = Taxonomy(EpacConfig.REF_SEQ_PREFIX, self.tax_fname)   
        tax_map = self.taxonomy.get_map()    
        cfg = EpacConfig()
        self.taxtree_helper = TaxTreeHelper(cfg, tax_map)

        outgr_fname = os.path.join(self.testfile_dir, "outgroup.nw")
        self.expected_outgr = Tree(outgr_fname)
    
    def tearDown(self):
        self.taxonomy = None
        self.taxtree_helper = None
   
    def test_outgroup(self):
        mfu_tree_fname = os.path.join(self.testfile_dir, "taxtree.nw")
        mfu_tree = Tree(mfu_tree_fname)
        self.taxtree_helper.set_mf_rooted_tree(mfu_tree)
        outgr = self.taxtree_helper.get_outgroup()  
        self.assertEqual(outgr.get_leaf_names(), self.expected_outgr.get_leaf_names())

    def test_branch_labeling(self):
        bfu_tree_fname = os.path.join(self.testfile_dir, "resolved_tree.nw")
        bfu_tree = Tree(bfu_tree_fname)
        map_fname = os.path.join(self.testfile_dir, "bid_tax_map.txt")
        self.expected_map = {}
        with open(map_fname) as inf:
            for line in inf:
                bid, rank_id, rdiff, brlen = line.strip().split("\t")
                self.expected_map[bid] = (rank_id, int(rdiff), float(brlen))
                
        self.taxtree_helper.set_outgroup(self.expected_outgr)
        self.taxtree_helper.set_bf_unrooted_tree(bfu_tree)
        bid_tax_map = self.taxtree_helper.get_bid_taxonomy_map()
        self.assertEqual(len(bid_tax_map), 2 * len(bfu_tree) - 3)
        for bid in self.expected_map.keys():
            e_rec = self.expected_map[bid]
            rec = bid_tax_map[bid]
            self.assertEqual(e_rec[0], rec[0])
            self.assertEqual(e_rec[1], rec[1])
            self.assertAlmostEqual(e_rec[2], rec[2], 6)
                
class TaxClassifyHelperTests(unittest.TestCase):

    def setUp(self):
        self.testfile_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "testfiles")
        self.tax_fname = os.path.join(self.testfile_dir, "test_clean.tax")
        cfg = EpacClassifierConfig()

        map_fname = os.path.join(self.testfile_dir, "bid_tax_map2.txt")
        self.bid_tax_map = {}
        with open(map_fname) as inf:
            for line in inf:
                bid, rank_id, rdiff, brlen = line.strip().split("\t")
                self.bid_tax_map[bid] = (rank_id, int(rdiff), float(brlen))
        
        self.classify_helper = TaxClassifyHelper(cfg, self.bid_tax_map)

    def test_assign_taxonomy(self):
        assign_fname = os.path.join(self.testfile_dir, "true_assign.txt")
        expected_assign_map = {}
        with open(assign_fname) as inf:
            for line in inf:
                sid, ranks_str, lws = line.strip().split("\t")
                expected_assign_map[sid] = ranks_str.split(";")
    
        jplace_fname = os.path.join(self.testfile_dir, "test.jplace")
        parser = EpaJsonParser(jplace_fname)
        for p in parser.get_placement():
            sid = p["n"][0]
            edges = p["p"]
            ranks, conf = self.classify_helper.classify_seq(edges)
#            for e in edges: print self.bid_tax_map[str(e[0])], e[2]
#            print sid, "\t", ";".join(ranks) #, conf
            self.assertEqual(ranks, expected_assign_map[sid])
                
if __name__ == '__main__':
    unittest.main()
