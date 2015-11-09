#!/usr/bin/env python
import os
import sys
import unittest

lib_path = os.path.abspath('..')
sys.path.append(lib_path)

from epac.taxonomy_util import Taxonomy, TaxCode, TaxTreeBuilder
from epac.config import EpacConfig
from epac.ete2 import Tree

class TaxCodeTests(unittest.TestCase):

    def check_ranks(self, tests):
        for ranks_str, expected_levels_str in tests:
            ranks = ranks_str.split(";")
            expected_levels = expected_levels_str.split(";")
            for i in range(len(ranks)):
                guessed_level = self.taxcode.guess_rank_level_name(ranks, i)[0].lower()
#                print guessed_level
                self.assertEqual(guessed_level, expected_levels[i], \
                  msg="%s\n%s\n%s != %s" % (ranks_str, expected_levels_str, guessed_level, expected_levels[i]))
    
    def test_bac_taxcode(self):
        self.taxcode = TaxCode("BAC")
        tests = []
        ranks = "Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium;Clostridium rectum"
        expected_levels = "kingdom;phylum;class;order;family;genus;species"
        tests += [(ranks, expected_levels)]
        ranks = "Bacteria;Actinobacteria;Actinobacteria;Actinobacteridae;Actinomycetales;Micrococcineae;Micrococcaceae;Acaricomes;Acaricomes phytoseiuli"
        expected_levels = "kingdom;phylum;class;subclass;order;suborder;family;genus;species"
        tests += [(ranks, expected_levels)]
        ranks = "k__Bacteria; p__Acidobacteria; c__[Chloracidobacteria]; o__[Chloracidobacterales]; f__[Chloracidobacteraceae]; g__Candidatus Chloracidobacterium; s__"
        expected_levels = "kingdom;phylum;class;order;family;genus;species"
        tests += [(ranks, expected_levels)]
        self.check_ranks(tests)
        
    def test_bot_taxcode(self):
        self.taxcode = TaxCode("BOT")
        tests = []
        ranks = "Apocynoideae;Nerieae;Neriinae;Adenium;Adenium_swazicum"
        expected_levels = "subfamily;tribe;subtribe;genus;species"
        tests += [(ranks, expected_levels)]
        ranks = "Apocynoideae;Echiteae;Parsonsiinae;Parsonsia;Parsonsia_heterophylla"
        expected_levels = "subfamily;tribe;subtribe;genus;species"
        tests += [(ranks, expected_levels)]
#        ranks = "Eukaryota;Opisthokonta;Nucletmycea;Fungi;Dikarya;Ascomycota;Saccharomycotina;Saccharomycetes;Saccharomycetales;Metschnikowiaceae;Clavispora;Clavispora lusitaniae ATCC 42720"
        ranks = "Eukaryota;Fungi;Dikarya;Ascomycota;Saccharomycotina;Saccharomycetes;Saccharomycetales;Metschnikowiaceae;Clavispora;Clavispora lusitaniae ATCC 42720"
        expected_levels = "domain;kingdom;subkingdom;phylum;subphylum;class;order;family;genus;species"
        tests += [(ranks, expected_levels)]
        self.check_ranks(tests)
        
    def test_zoo_taxcode(self):
        self.taxcode = TaxCode("ZOO")
        tests = []
        ranks = "Animalia;Chordata;Mammalia;Rodentia;Muridae;Apodemus;agrarius"
        expected_levels = "kingdom;phylum;class;order;family;genus;species"
        tests += [(ranks, expected_levels)]
        ranks = "Animalia;Chordata;Mammalia;Primates;Hominidae;Homininae;Hominini;Hominina;Homo;sapiens"
        expected_levels = "kingdom;phylum;class;order;family;subfamily;tribe;subtribe;genus;species"
        tests += [(ranks, expected_levels)]
        self.check_ranks(tests)

class TaxonomyTests(unittest.TestCase):

    TAX_DICT = { "UnpC[Ceti]" : ["Bacteria", "Fusobacteria", "Fusobacteriia", "Fusobacteriales", "Fusobacteriaceae", "Cetobacterium", "Cetobacterium ceti"],
                 "UnpSomer," : ["Bacteria", "Fusobacteria", "Fusobacteriia", "Fusobacteriales", "Fusobacteriaceae", "Cetobacterium", "Cetobacterium somerae"],
                 "UpbRectu" : ["[Bacteria]", "'Firmicutes'", "Clostridia(1)", "Clostridiales", "Clostridiaceae", "Clostridium", "Clostridium rectum"]
    }
    
    def setUp(self):
        test_dir = os.path.dirname(os.path.abspath(__file__))
        self.tax_fname = os.path.join(test_dir, "test.tax")
        self.PREFIXED_TAX_DICT = {}
        with open(self.tax_fname, "w") as outf:
            for sid, ranks in self.TAX_DICT.iteritems():
                outf.write("%s\t%s\n" % (sid, ";".join(ranks)))
                self.PREFIXED_TAX_DICT[EpacConfig.REF_SEQ_PREFIX+sid] = ranks
        self.taxonomy = Taxonomy("", self.tax_fname)
    
    def tearDown(self):
        self.taxonomy = None
        os.remove(self.tax_fname)
    
    def test_load(self):
        self.assertEqual(self.TAX_DICT, self.taxonomy.seq_ranks_map)
        prefixed_tax = Taxonomy(EpacConfig.REF_SEQ_PREFIX, self.tax_fname)
        self.assertEqual(self.PREFIXED_TAX_DICT, prefixed_tax.seq_ranks_map)
        prefixed_tax = None
        
    def test_rank_uid(self):
        tax = self.taxonomy
        for sid in tax.get_map().iterkeys():
            self.assertEqual(tax.get_seq_ranks(sid), Taxonomy.split_rank_uid(tax.seq_rank_id(sid)))
            
    def test_normalize_seq_ids(self):
        tax = Taxonomy(tax_map=self.taxonomy.seq_ranks_map)
        self.assertTrue("UnpC[Ceti]" in tax.seq_ranks_map)
        self.assertTrue("UnpSomer," in tax.seq_ranks_map)
        tax.normalize_seq_ids()
        self.assertFalse("UnpC[Ceti]" in tax.seq_ranks_map)
        self.assertTrue("UnpC_Ceti_" in tax.seq_ranks_map)
        self.assertFalse("UnpSomer," in tax.seq_ranks_map)
        self.assertTrue("UnpSomer_" in tax.seq_ranks_map)
        
    def test_normalize_rank_names(self):
        tax = Taxonomy(tax_map=self.taxonomy.seq_ranks_map)
        ranks = tax.get_seq_ranks("UpbRectu")
        self.assertEqual(ranks[0], "[Bacteria]")
        self.assertEqual(ranks[1], "'Firmicutes'")
        self.assertEqual(ranks[2], "Clostridia(1)")
        corr_ranks = tax.normalize_rank_names()
        self.assertEqual(len(corr_ranks), 3)
        ranks = tax.get_seq_ranks("UpbRectu")
        self.assertEqual(ranks[0], "_Bacteria_")
        self.assertEqual(ranks[1], "_Firmicutes_")
        self.assertEqual(ranks[2], "Clostridia_1_")
        
    def test_merge_ranks(self):
        tax = Taxonomy(tax_map=self.taxonomy.seq_ranks_map)
        merge_sids = ["UnpC[Ceti]", "UnpSomer,"]
        rank_ids = [tax.seq_rank_id(sid) for sid in merge_sids]
        new_rank_id = tax.merge_ranks(rank_ids)
        self.assertEqual(merge_sids, tax.get_rank_seqs(new_rank_id))
        for sid in merge_sids:
            self.assertEqual(tax.seq_rank_id(sid), new_rank_id)
            
    def test_taxtree_builder(self):
        cfg = EpacConfig()
        testfile_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "testfiles")
        tax_fname = os.path.join(testfile_dir, "test.tax")
        tax = Taxonomy(EpacConfig.REF_SEQ_PREFIX, tax_fname)
        tree_fname = os.path.join(testfile_dir, "taxtree.nw")
        expected_tree = Tree(tree_fname, format=8)
        tb = TaxTreeBuilder(cfg, tax)
        tax_tree, seq_ids = tb.build()
        self.assertEqual(seq_ids, tax.get_map().keys())
        self.assertEqual(tax_tree.write(format=8), expected_tree.write(format=8))

if __name__ == '__main__':
    unittest.main()
