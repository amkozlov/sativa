#! /usr/bin/env python
import sys
import os
import json
import operator
import base64
from subprocess import call
from ete2 import Tree, SeqGroup
from taxonomy_util import TaxCode

class EpaJsonParser:
    """This class parses the RAxML-EPA json output file"""
    def __init__(self, jsonfin):
        self.jdata = json.load(open(jsonfin))
    
    def get_placement(self):
        return self.jdata["placements"]
        
    def get_tree(self):
        return self.jdata["tree"]
        
    def get_std_newick_tree(self):
        tree = self.jdata["tree"]
        tree = tree.replace("{", "[&&NHX:B=")
        tree = tree.replace("}", "]")
        return tree

    def get_raxml_version(self):
        return self.jdata["metadata"]["raxml_version"]

    def get_raxml_invocation(self):
        return self.jdata["metadata"]["invocation"]

class RefJsonChecker:
    def __init__(self, jsonfin= None, jdata = None):
        if jsonfin!=None:
            self.jdata = json.load(open(jsonfin))
        else:
            self.jdata = jdata
    
    def check_field(self, fname, ftype, fvals=None, fopt=False):
        if fname in self.jdata:
            field = self.jdata[fname]
            if isinstance(field, ftype):
                if not fvals or field in fvals:
                    return True
                else:
                    self.error = "Invalid value of field '%s': %s" % fname, repr(field)
                    return False
            else:
              self.error = "Field '%s' has a wrong type: %s (expected: %s)" % fname, type(field).__name__, ftype.__name__
              return False
        else:
            if fopt:
                return True
            else:
                self.error = "Field not found: %s" % fname
                return False
    
    def validate(self, ver = "1.6"):
        nver = float(ver)
        
        self.error = None

        valid = self.check_field("tree", unicode) \
                and self.check_field("raxmltree", unicode) \
                and self.check_field("rate", float) \
                and self.check_field("node_height", dict) \
                and self.check_field("origin_taxonomy", dict) \
                and self.check_field("sequences", list) \
                and self.check_field("binary_model", unicode) \
                and self.check_field("hmm_profile", list, fopt=True) 
                
        # check v1.1 fields, if needed
        if nver >= 1.1:
            valid = valid and \
                    self.check_field("ratehet_model", unicode)  # ["GTRGAMMA", "GTRCAT"]

        # check v1.2 fields, if needed
        if nver >= 1.2:
            valid = valid and \
                    self.check_field("tax_tree", unicode)

        # check v1.3 fields, if needed
        if nver >= 1.3:
            valid = valid and \
                    self.check_field("taxcode", unicode, TaxCode.TAX_CODE_MAP)
        
        # check v1.4 fields, if needed
        if nver >= 1.4:
            valid = valid \
                    and self.check_field("corr_seqid_map", dict) \
                    and self.check_field("corr_ranks_map", dict)

        # check v1.5 fields, if needed
        if nver >= 1.5:
            valid = valid \
                    and self.check_field("merged_ranks_map", dict)

        # "taxonomy" field has been renamed and its format was changed in v1.6
        if nver >= 1.6:
            valid = valid \
                    and self.check_field("branch_tax_map", dict)
        else:
            valid = valid \
                    and self.check_field("taxonomy", dict)

        return (valid, self.error)

class RefJsonParser:
    """This class parses the EPA Classifier reference json file"""
    def __init__(self, jsonfin):
        self.jdata = json.load(open(jsonfin))
        self.version = self.jdata["version"]
        self.nversion = float(self.version)
        self.corr_seqid = None
        self.corr_ranks = None
        self.corr_seqid_reverse = None
        
    def validate(self):
        jc = RefJsonChecker(jdata = self.jdata)
        return jc.validate(self.version)
    
    def get_version(self):
        return self.version

    def get_rate(self):
        return self.jdata["rate"]
    
    def get_node_height(self):
        return self.jdata["node_height"]
    
    def get_raxml_readable_tree(self, fout_name = None):
        tree_str = self.jdata["raxmltree"]
        #t.unroot()
        if fout_name != None:
            with open(fout_name, "w") as fout:
                fout.write(tree_str)
        else:
            return tree_str
    
    def get_reftree(self, fout_name = None):
        tree_str = self.jdata["tree"]
        if fout_name != None:
            with open(fout_name, "w") as fout:
                fout.write(tree_str)
        else:
            return Tree(tree_str, format=1)
    
    def get_tax_tree(self):
        t = Tree(self.jdata["tax_tree"], format=8)
        return t

    def get_outgroup(self):
        t = Tree(self.jdata["outgroup"], format=9)
        return t

    def get_branch_tax_map(self):
        if self.nversion >= 1.6:
            return self.jdata["branch_tax_map"]
        else:
            return None

    def get_taxonomy(self):
        if self.nversion < 1.6:
            return self.jdata["taxonomy"]
        else:
            return None

    def get_origin_taxonomy(self):
        return self.jdata["origin_taxonomy"]
    
    def get_alignment(self, fout):
        entries = self.jdata["sequences"]
        with open(fout, "w") as fo:
            for entr in entries:
                fo.write(">%s\n%s\n" % (entr[0], entr[1]))

        return fout
    
    def get_ref_alignment(self):
        entries = self.jdata["sequences"]
        alignment = SeqGroup()
        for entr in entries:
            alignment.set_seq(entr[0], entr[1])
        return alignment
    
    def get_alignment_list(self):
        return self.jdata["sequences"]
    
    def get_sequences_names(self):
        nameset = set()
        entries = self.jdata["sequences"]
        for entr in entries:
            nameset.add(entr[0])
        return nameset
    
    def get_alignment_length(self):
        entries = self.jdata["sequences"]
        return len(entries[0][1])
    
    def get_hmm_profile(self, fout):
        if "hmm_profile" in self.jdata:        
            lines = self.jdata["hmm_profile"]
            with open(fout, "w") as fo:
                for line in lines:
                    fo.write(line)
            return fout
        else:
            return None

    def get_binary_model(self, fout):
        model_str = self.jdata["binary_model"]
        with open(fout, "wb") as fo:
            fo.write(base64.b64decode(model_str))

    def get_ratehet_model(self):
        return self.jdata["ratehet_model"]

    def get_pattern_compression(self):
        if "pattern_compression" in self.jdata:
            return self.jdata["pattern_compression"]
        else:
            return False

    def get_taxcode(self):
        return self.jdata["taxcode"]
        
    def get_corr_seqid_map(self):
        if "corr_seqid_map" in self.jdata:
            self.corr_seqid = self.jdata["corr_seqid_map"]
        else:
            self.corr_seqid = {}
        return self.corr_seqid

    def get_corr_ranks_map(self):
        if "corr_ranks_map" in self.jdata:
            self.corr_ranks = self.jdata["corr_ranks_map"]
        else:
            self.corr_ranks = {}
        return self.corr_ranks

    def get_merged_ranks_map(self):
        if "merged_ranks_map" in self.jdata:
            self.merged_ranks = self.jdata["merged_ranks_map"]
        else:
            self.merged_ranks = {}
        return self.merged_ranks

    def get_metadata(self):
        return self.jdata["metadata"]
        
    def get_field_string(self, field_name):
        if field_name in self.jdata:
            return json.dumps(self.jdata[field_name], indent=4, separators=(',', ': ')).strip("\"")
        else:
            return None
            
    def get_uncorr_seqid(self, new_seqid):
        if not self.corr_seqid:
            self.get_corr_seqid_map()
        return self.corr_seqid.get(new_seqid, new_seqid)
        
    def get_corr_seqid(self, old_seqid):
        if not self.corr_seqid_reverse:
            if not self.corr_seqid:
                self.get_corr_seqid_map()
            self.corr_seqid_reverse = dict((reversed(item) for item in self.corr_seqid.items()))
        return self.corr_seqid_reverse.get(old_seqid, old_seqid)

    def get_uncorr_ranks(self, ranks):
        if not self.corr_ranks:
            self.get_corr_ranks_map()
        uncorr_ranks = list(ranks)
        for i in range(len(ranks)):
            uncorr_ranks[i] = self.corr_ranks.get(ranks[i], ranks[i])
        return uncorr_ranks        
                
class RefJsonBuilder:
    """This class builds the EPA Classifier reference json file"""
    def __init__(self, old_json=None):
        if old_json:
            self.jdata = old_json.jdata
        else:
            self.jdata = {}
            self.jdata["version"] = "1.6"
#            self.jdata["author"] = "Jiajie Zhang"
        
    def set_branch_tax_map(self, bid_ranks_map):
        self.jdata["branch_tax_map"] = bid_ranks_map

    def set_origin_taxonomy(self, orig_tax_map):
        self.jdata["origin_taxonomy"] = orig_tax_map

    def set_tax_tree(self, tr):
        self.jdata["tax_tree"] = tr.write(format=8)

    def set_tree(self, tr):
        self.jdata["tree"] = tr
        self.jdata["raxmltree"] = Tree(tr, format=1).write(format=5)
        
    def set_outgroup(self, outgr):
        self.jdata["outgroup"] = outgr.write(format=9)

    def set_sequences(self, seqs):    
        self.jdata["sequences"] = seqs
        
    def set_hmm_profile(self, fprofile):    
        with open(fprofile) as fp:
            lines = fp.readlines()
        self.jdata["hmm_profile"] = lines
       
    def set_rate(self, rate):    
        self.jdata["rate"] = rate
        
    def set_nodes_height(self, height):    
        self.jdata["node_height"] = height

    def set_binary_model(self, model_fname):  
        with open(model_fname, "rb") as fin:
            model_str = base64.b64encode(fin.read())
        self.jdata["binary_model"] = model_str

    def set_ratehet_model(self, model):
        self.jdata["ratehet_model"] = model

    def set_pattern_compression(self, value):
        self.jdata["pattern_compression"] = value

    def set_taxcode(self, value):
        self.jdata["taxcode"] = value

    def set_corr_seqid_map(self, seqid_map):
        self.jdata["corr_seqid_map"] = seqid_map

    def set_corr_ranks_map(self, ranks_map):
        self.jdata["corr_ranks_map"] = ranks_map

    def set_merged_ranks_map(self, merged_ranks_map):
        self.jdata["merged_ranks_map"] = merged_ranks_map

    def set_metadata(self, metadata):    
        self.jdata["metadata"] = metadata

    def dump(self, out_fname):
        self.jdata.pop("fields", 0)
        self.jdata["fields"] = self.jdata.keys()
        with open(out_fname, "w") as fo:
            json.dump(self.jdata, fo, indent=4, sort_keys=True)                


if __name__ == "__main__":
    if len(sys.argv) < 2: 
        print("usage: ./json_util.py jsonfile")
        sys.exit()
    jc = json_checker(jsonfin = sys.argv[1])
    if jc.valid():
        print("The json file is OK for EPA-classifer")
    else:
        print("!!!Invalid json file!!!")
    
