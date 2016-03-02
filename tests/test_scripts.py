#!/usr/bin/env python
import os
import sys
import unittest
import shutil
from subprocess import call,check_output,STDOUT,CalledProcessError

lib_path = os.path.abspath('..')
sys.path.append(lib_path)

from epac.taxonomy_util import Taxonomy, TaxCode, TaxTreeBuilder
from epac.config import EpacConfig, EpacClassifierConfig
from epac.ete2 import Tree
from epac.classify_util import TaxTreeHelper, TaxClassifyHelper
from epac.json_util import EpaJsonParser

class ScriptTests(unittest.TestCase):

    def setUp(self):
        self.testfile_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "testfiles")
        self.sativa_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
        self.example_dir = os.path.join(self.sativa_dir, "example")
        self.out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tmpout")
        os.mkdir(self.out_dir)
    
    def tearDown(self):
        shutil.rmtree(self.out_dir)
   
    def test_sativa(self):
        exec_script = os.path.join(self.sativa_dir, "sativa.py")
        phy_fname = os.path.join(self.testfile_dir, "ref.phy")
        tax_fname = os.path.join(self.testfile_dir, "full.tax")
        call_str = [exec_script, "-s", phy_fname, "-t", tax_fname, "-n", "testsat", "-x", "BAC", "-o", self.out_dir, "-T", "2"]
        try:
            out_str = check_output(call_str, stderr=STDOUT)
        except CalledProcessError as ex:
            print "\n\nCommand line: %s\n\nOutput:\n%s\n" % (ex.cmd, ex.output)
            self.assertTrue(False, msg="Error running SATIVA script")
            
        mis_fname = os.path.join(self.out_dir, "testsat.mis")
        self.assertTrue(os.path.isfile(mis_fname))
        
        
    def test_classifier(self):
        exec_script = os.path.join(self.sativa_dir, "epa_trainer.py")
        ali_fname = os.path.join(self.testfile_dir, "ref.phy")
        tax_fname = os.path.join(self.testfile_dir, "ref.tax")
        call_str = [exec_script, "-s", ali_fname, "-t", tax_fname, "-n", "testcl", "-x", "BAC", "-o", self.out_dir, "-no-hmmer", "-T", "2"]
        try:
            out_str = check_output(call_str, stderr=STDOUT)
#            call(call_str)
        except CalledProcessError as ex:
            print "\n\nCommand line: %s\n\nOutput:\n%s\n" % (ex.cmd, ex.output)
            self.assertTrue(False, msg="Error running epa_trainer.py script")

        refjson_fname = os.path.join(self.out_dir, "testcl.refjson")
        self.assertTrue(os.path.isfile(refjson_fname))
        
        exec_script = os.path.join(self.sativa_dir, "epa_classifier.py")
        query_fname = os.path.join(self.testfile_dir, "ref.phy")
        call_str = [exec_script, "-r", refjson_fname, "-q", query_fname, "-n", "testcl", "-x", "-o", self.out_dir, "-T", "2"]
        
        try:
            out_str = check_output(call_str, stderr=STDOUT)
#            call(call_str)
        except CalledProcessError as ex:
            print "\n\nCommand line: %s\n\nOutput:\n%s\n" % (ex.cmd, ex.output)
            self.assertTrue(False, msg="Error running epa_classifier.py script")        

        assign_fname = os.path.join(self.out_dir, "testcl.assignment.txt")
        self.assertTrue(os.path.isfile(assign_fname))
                
if __name__ == '__main__':
    unittest.main()
