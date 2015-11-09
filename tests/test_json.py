#!/usr/bin/env python
import os
import sys
import unittest

lib_path = os.path.abspath('..')
sys.path.append(lib_path)

from epac.ete2 import SeqGroup
from epac.taxonomy_util import Taxonomy, TaxCode
from epac.config import EpacConfig
from epac.json_util import EpaJsonParser, RefJsonParser, RefJsonBuilder
from epac.ete2 import Tree

class JsonTests(unittest.TestCase):

    def setUp(self):
        self.cfg = EpacConfig()
        self.testfile_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "testfiles")
        tax_fname = os.path.join(self.testfile_dir, "test.tax")
        phy_fname = os.path.join(self.testfile_dir, "test.phy")
        tax = Taxonomy(EpacConfig.REF_SEQ_PREFIX, tax_fname)
        seqs = SeqGroup(sequences=phy_fname, format = "phylip")
    
    def test_jplace_read(self):
        jplace_fname = os.path.join(self.testfile_dir, "test.jplace")
        parser = EpaJsonParser(jplace_fname)
        self.assertEquals(parser.get_raxml_version(), "8.2.3")
        t = Tree(parser.get_tree())
        t_len = len(t)
        self.assertEquals(t_len, 32)
        self.assertEquals(len(parser.get_placement()), 6)
        for p in parser.get_placement():
            self.assertFalse(p["n"][0] in t) 
            self.assertTrue(len(p["p"]) > 0) 
            for edge in p["p"]:
                branch = int(edge[0])
                lh = edge[1]
                lhw = edge[2]
                self.assertTrue(branch >= 0 and branch < (t_len * 2 - 3)) 
                self.assertTrue(lhw >= 0.0 and lhw <= 1.0) 
        
    def test_refjson_read(self):
        versions = ["1.4", "1.5"]
        for ver in versions:
            ref_fname = os.path.join(self.testfile_dir, "test.refjson.v%s" % ver)
            parser = RefJsonParser(ref_fname, ver)
            valid, errors = parser.validate()
            self.assertTrue(valid) 
        
        
if __name__ == '__main__':
    unittest.main()
