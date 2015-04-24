#!/usr/bin/env python

import sys
import re
from ete2 import Tree
#from epac.ete2 import Tree

class TaxCode:
    UNI_TAX_RANKS = {
         1: ("Kingdom", "k__"),
         2: ("Phylum", "p__"),
         3: ("Subphylum", "a__"),
         4: ("Class", "c__"),
         5: ("Subclass", "d__"),
         6: ("Superorder", "e__"),
         7: ("Order", "o__"),
         8: ("Suborder", "h__"),
         9: ("Infraorder", "i__"),
         10: ("Superfamily", "j__"),
         11: ("Epifamily", "l__"),
         12: ("Family", "f__"), 
         13: ("Subfamily", "m__"),
         14: ("Infrafamily", "n__"),
         15: ("Tribe", "t__"),
         16: ("Subtribe", "u__"),
         17: ("Infratribe", "v__"), 
         18: ("Genus", "g__"),
         19: ("Species", "s__"),
         20: ("Subspecies", "b__"),
         21: ("Strain", "r__"),
         22: ("Isolate", "q__")
        }        
    UNI_TAX_LEVELS = len(UNI_TAX_RANKS)
    
    BAC_TAX_CODE = {
         1: ((), ("bacteria", "archaea")),   # kingdom
         2: ((), ()),                        # phylum
         4: ((), ()),                        # class
         5: (("idae"), ()),                  # subclass
         7: (("ales"), ()),                  # order
         8: (("ineae"), ()),                 # suborder
         12: (("aceae"), ()),                # family
         13: (("oideae"), ()),               # subfamily
         15: (("eae"), ()),                  # tribe
         16: (("inae"), ()),                 # subtribe
         18: ((), ()),                       # genus
         19: ((), ()),                       # species
         20: ((), ()),                       # subspecies
         21: ((), ()),                       # strain
         22: ((), ())                        # isolate
    }

    BOT_TAX_CODE = {
         1: ((), ("plantae", "algae", "fungi")),                # kingdom
         2: (("phyta", "phycota", "mycota"), ()),               # phylum
         3: (("phytina", "phycotina", "mycotina"), ()),         # subphylum
         4: (("opsida", "phyceae", "mycetes"), ()),             # class
         5: (("idae", "phycidae", "mycetidae"), ()),            # subclass
         6: (("anae"), ()),                                     # superorder
         7: (("ales"), ()),                                     # order
         8: (("ineae"), ()),                                    # suborder
         9: (("aria"), ()),                                     # infraorder
         10: (("acea"), ()),                                    # superfamily
         12: (("aceae"), ()),                                   # family
         13: (("oideae"), ()),                                  # subfamily
         15: (("eae"), ()),                                     # tribe
         16: (("inae"), ()),                                    # subtribe
         18: ((), ()),                                          # genus
         19: ((), ()),                                          # species
         20: ((), ()),                                          # subspecies
         21: ((), ()),                                          # specimen
         22: ((), ())                                           # isolate
    }

    ZOO_TAX_CODE = {
         1: ((), ("animalia")),                                               # kingdom
         2: ((), ("chordata", "arthropoda", "mollusca", "nematoda")),         # phylum
         3: ((), ("vertebrata", "myriapoda", "crustacea", "hexapoda")),       # subphylum
         4: ((), ("mammalia", "aves", "reptilia", "amphibia", "insecta")),    # class
         5: ((), ()),                                                         # subclass
         6: ((), ()),                                                         # superorder
         7: ((), ()),                                                         # order
         8: ((), ()),                                                         # suborder
         9: ((), ()),                                                         # infraorder
         10: (("oidea"), ()),                                                 # superfamily
         11: (("oidae"), ()),                                                 # epifamily
         12: (("idae"), ()),                                                  # family
         13: (("inae"), ()),                                                  # subfamily
         13: (("odd"), ()),                                                   # infrafamily
         15: (("ini"), ()),                                                   # tribe
         16: (("ina"), ()),                                                   # subtribe
         17: (("ad", "iti"), ()),                                             # infratribe
         18: ((), ()),                                                        # genus
         19: ((), ()),                                                        # species
         20: ((), ()),                                                        # subspecies
         21: ((), ()),                                                        # specimen
         22: ((), ())                                                         # isolate
    }
    
    VIR_TAX_CODE = {
         1: ((), ("viruses")),      # kingdom
         4: ((), ()),               # class
         5: (("idae"), ()),         # subclass
         7: (("virales"), ()),      # order
         12: (("viridae"), ()),     # family
         13: (("virinae"), ()),     # subfamily
         18: (("virus"), ()),       # genus
         19: ((" virus"), ()),      # species
         21: ((), ()),              # strain
         22: ((), ())               # isolate
    }
    
    TAX_CODE_MAP = {
        "bac": BAC_TAX_CODE,
        "bot": BOT_TAX_CODE,
        "zoo": ZOO_TAX_CODE,
        "vir": VIR_TAX_CODE
    }

    def __init__(self, tax_code_name):
        self.tax_code = TaxCode.TAX_CODE_MAP.get(tax_code_name.lower(), None)
        if not self.tax_code:
            print "ERROR: Unknown taxonomic code: %s" % tax_code_name
            sys.exit()
        
    @staticmethod    
    def rank_level_name(uni_rank_level):
        return TaxCode.UNI_TAX_RANKS.get(uni_rank_level, ("Unknown", "?__"))
                            
    def guess_rank_level(self, ranks, rank_level):
        rank_name = re.sub("[\W_]+", "", ranks[rank_level].lower())
        
        sorted_tax_levels = sorted(self.tax_code.keys())
        
        real_level = 0
        
        # first, try to guess rank level based on its name or name suffix
        for lvl in sorted_tax_levels:
          (suffix, exact_match) = self.tax_code[lvl]
          if rank_name.endswith(suffix) or rank_name in exact_match:
            real_level = lvl
            break
        
        # if name-based identification failed, try to derive rank level from its parent
        if real_level == 0:
            if rank_level == 0:    # kingdom
                real_level = 1
            else:
                parent_level = self.guess_rank_level(ranks, rank_level-1)
                idx = sorted_tax_levels.index(parent_level)
                for i in range(idx+1, len(sorted_tax_levels)):
                  lvl = sorted_tax_levels[i]
                  suffix = self.tax_code[lvl][0]
                  # here we assume that if rank has a defined suffix, if would have been identified in the previous step
                  if len(suffix) == 0 or i == len(sorted_tax_levels):
                    real_level = lvl
                    break
                             
        return real_level
         
    def guess_rank_level_name(self, ranks, rank_level):
        real_level = self.guess_rank_level(ranks, rank_level)
        return TaxCode.rank_level_name(real_level)


class Taxonomy:
    EMPTY_RANK = "-"
    RANK_UID_DELIM = "@@"
    
    RANK_PLACEHOLDERS = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]    

    @staticmethod    
    def lineage_str(ranks):
        return ";".join(ranks).strip(';')

    @staticmethod    
    def lowest_assigned_rank_level(ranks):
        rank_level = len(ranks)-1
        while ranks[rank_level] == Taxonomy.EMPTY_RANK:
            rank_level -= 1
        return rank_level

    @staticmethod    
    def lowest_assigned_rank(ranks):
        rank_level = Taxonomy.lowest_assigned_rank_level(ranks)
        return ranks[rank_level]
        
    @staticmethod    
    def get_rank_uid(ranks, rank_level=-1):
        if rank_level == -1:
            rank_level = Taxonomy.lowest_assigned_rank_level(ranks)        
        return Taxonomy.RANK_UID_DELIM.join(ranks[:rank_level+1])

    @staticmethod    
    def split_rank_uid(rank_uid, min_lvls=0):
        ranks = rank_uid.split(Taxonomy.RANK_UID_DELIM)
        if len(ranks) < min_lvls:
            return ranks + [Taxonomy.EMPTY_RANK] * (min_lvls - len(ranks))
        else:
            return ranks
            
    def __init__(self, tax_fname="", prefix=""):
        self.prefix = prefix
        self.seq_ranks_map = {}
        tree_nodes = []
        self.common_ranks = set([])
        if tax_fname:
            self.load_taxonomy(tax_fname)

    def get_common_ranks(self):
        allranks = list(self.seq_ranks_map.items())
        numitems = len(allranks)
        if numitems > 0:
            self.common_ranks = set(allranks[0][1])
            for i in range(1, numitems):
                curr_set = set(allranks[i][1])
                inters = self.common_ranks & curr_set 
                self.common_ranks = inters
            return self.common_ranks
        else:
            return set([])

    def seq_count(self):
        return len(self.seq_ranks_map)

    def items(self):
        return self.seq_ranks_map.items()
        
    def iteritems(self):
        return self.seq_ranks_map.items()

    def get_map(self):
        return self.seq_ranks_map

    def remove_seq(self, seqid):
        del self.seq_ranks_map[seqid]

    def get_seq_ranks(self, seq_id):
        return self.seq_ranks_map[seq_id]
        
    def seq_lineage_str(self, seq_id):
        ranks = list(self.seq_ranks_map[seq_id])
        return Taxonomy.lineage_str(ranks)        

    def load_taxonomy(self, tax_fname):
        fin = open(tax_fname)
        for line in fin:
            line = line.strip()
            toks = line.split("\t")
            sid = self.prefix + toks[0]
            ranks_str = toks[1]
            ranks = ranks_str.split(";")
            for i in range(len(ranks)):
                rank_name = ranks[i].strip()
#                if rank_name in GGTaxonomyFile.rank_placeholders:
#                    rank_name = Taxonomy.EMPTY_RANK
                ranks[i] = rank_name
                
            self.seq_ranks_map[sid] = ranks     

        fin.close()

    def normalize_rank_name(self, rank, rank_name):
        invalid_chars = ['(', ')', ',', ';', ':']
        for ch in invalid_chars:
            rank_name = rank_name.replace(ch, "_")
        return rank_name

    def normalize_rank_names(self):
        for sid, ranks in self.seq_ranks_map.iteritems():
            for i in range(1, len(ranks)):
                ranks[i] = self.normalize_rank_name(i, ranks[i])
                
    def close_taxonomy_gaps(self):
        for sid, ranks in self.seq_ranks_map.iteritems():
            last_rank = None
            gap_count = 0
            for i in reversed(range(1, len(ranks))):
                if ranks[i] != Taxonomy.EMPTY_RANK:
                    last_rank = ranks[i]
                elif last_rank:
                    gap_count += 1
                    ranks[i] = "parent%d_" % gap_count + last_rank

    def check_for_duplicates(self, autofix=False):
        parent_map = {}
        dups = []
        old_fixed = {}
        old_ranks = {}
        for sid, ranks in self.seq_ranks_map.iteritems():
            for i in range(1, len(ranks)):
                if ranks[i] == Taxonomy.EMPTY_RANK:
                    break                
                parent = ranks[i-1]
                if not ranks[i] in parent_map:
                    parent_map[ranks[i]] = sid
                else:
                    old_sid = parent_map[ranks[i]]
                    if self.get_seq_ranks(old_sid)[i-1] != parent:
                        if autofix:
                            orig_name = self.lineage_str(sid)
                            orig_old_name = self.lineage_str(old_sid)

                            if old_sid not in old_fixed:
                                old_fixed[old_sid] = self.lineage_str(old_sid)
                                old_ranks[old_sid] = set([])
                                
                            if i not in old_ranks[old_sid]:
                                self.seq_ranks_map[old_sid][i] += "_" + self.seq_ranks_map[old_sid][i-1]
                                old_ranks[old_sid].add(i)

                            self.seq_ranks_map[sid][i] = ranks[i] + "_" + parent
                            dup_rec = (old_sid, orig_old_name, sid, orig_name, self.lineage_str(sid))
                        else:
                            dup_rec = (old_sid, self.lineage_str(old_sid), sid, self.lineage_str(sid))
                        dups.append(dup_rec)
                        
        if autofix:
            for sid, orig_name in old_fixed.iteritems():
                dup_rec = (sid, orig_name, sid, orig_name, self.lineage_str(sid))
                dups.append(dup_rec)

        return dups
        
    def check_for_disbalance(self, autofix=False):
        # the next block finds "orphan" ranks - could be used to decide which ranks to drop (not used now)
        if 0:
            child_map = {}
            for sid, ranks in self.seq_ranks_map.iteritems():
                for i in range(1, len(ranks)):
                    parent = "%d_%s" % (i-1, ranks[i-1])
                    if parent not in child_map:
                        child_map[parent] = set([ranks[i]])
                    else:
                        child_map[parent].add(ranks[i])
        
        errs = []
        for sid, ranks in self.seq_ranks_map.iteritems():
            if len(ranks) > 7:
                if autofix:
                    orig_name = self.lineage_str(sid)

                    dropq = []
                    keepq = []
                    restq = []
                    for i in range(1, len(ranks)):
                        # drop Subclass and Suborder, preserve Order and Family (based on common suffixes)
                        if ranks[i].endswith("dae") or ranks[i].endswith("neae"):
                            dropq += [i]
                        elif ranks[i].endswith("ceae") or ranks[i].endswith("ales"):
                            keepq += [i]
                        else:
                            restq += [i]

                    to_remove = dropq + restq + keepq
                    to_remove = to_remove[:len(ranks)-7]
                    
                    new_ranks = []
                    for i in range(len(ranks)):
                        if i not in to_remove:
                            new_ranks += [ranks[i]]
                            
                    self.seq_ranks_map[sid] = new_ranks
                    
                    err_rec = (sid, orig_name, self.lineage_str(sid))
                else:    
                    err_rec = (sid, self.lineage_str(sid))
                errs.append(err_rec)

        return errs
        
        
class TaxTreeBuilder:
    ROOT_LABEL = "<<root>>"

    def __init__(self, config, taxonomy):
        self.tree_nodes = {}
        self.leaf_count = {}
        self.config = config
        self.taxonomy = taxonomy

    def prune_unifu_nodes(self, tree):
        for node in tree.traverse("preorder"):
            if len(node.children) == 1:
                node.delete()

    def add_tree_node(self, tree, nodeId, ranks, rank_level):
        if rank_level >= 0:
            parent_level = rank_level            
            while ranks[parent_level] == Taxonomy.EMPTY_RANK:
                parent_level -= 1
            parentId = Taxonomy.get_rank_uid(ranks, parent_level)
        else:
            parentId = TaxTreeBuilder.ROOT_LABEL
            parent_level = -1

        if (parentId in self.tree_nodes):
            parentNode = self.tree_nodes[parentId]
        else:
            parentNode = self.add_tree_node(tree, parentId, ranks, parent_level-1)
            self.tree_nodes[parentId] = parentNode;

#        print "Adding node: %s, parent: %s, parent_level: %d" % (nodeId, parentId, parent_level)
        newNode = parentNode.add_child()
        newNode.add_feature("name", nodeId)
        return newNode        

    def build(self, min_rank=0, max_seqs_per_leaf=1e9, clades_to_include=[], clades_to_ignore=[]):

        if self.config.verbose:
            print "Number of nodes: %d" % self.taxonomy.seq_count()
        
        t0 = Tree()
        t0.add_feature("name", TaxTreeBuilder.ROOT_LABEL)
        self.tree_nodes[TaxTreeBuilder.ROOT_LABEL] = t0;
        self.leaf_count[TaxTreeBuilder.ROOT_LABEL] = 0;
        k = 0
        added = 0
        seq_ids = []
        # sequences are leafs of the tree, so they always have the lowest taxonomy level (e.g. "species"+1)        
        for sid, ranks in self.taxonomy.iteritems():
            k += 1
            if self.config.verbose and k % 1000 == 0:
                print "Processed nodes: ", k, ", added: ", added, ", skipped: ", k - added

            # filter by minimum rank level            
            if ranks[min_rank] == Taxonomy.EMPTY_RANK:
                continue       
    
            # filter by rank contraints (e.g. class Clostridia only)
            clade_is_ok = False

            # check against the inclusion list            
            if len(clades_to_include) > 0:
                for (rank_level, rank_name) in clades_to_include:            
                    if ranks[rank_level] == rank_name:
                        clade_is_ok = True
                        break
            else: # default: include all
                clade_is_ok = True

            # if sequence is about to be included, check it against the ignore list
            if clade_is_ok:
                for (rank_level, rank_name) in clades_to_ignore:
                    if ranks[rank_level] == rank_name:
                        clade_is_ok = False
                        break

            # final decision
            if not clade_is_ok:
                continue

            tax_seq_level = len(ranks)
            parent_level = tax_seq_level - 1
            while ranks[parent_level] == Taxonomy.EMPTY_RANK:
                parent_level -= 1
            parent_name = Taxonomy.get_rank_uid(ranks, parent_level)
            if parent_name in self.tree_nodes:
                parent_node = self.tree_nodes[parent_name]
#                max_seq_per_rank = max_seqs_per_leaf * (tax_seq_level - parent_level)                
                if parent_level == tax_seq_level - 1:
                    max_seq_per_rank = max_seqs_per_leaf # * (tax_seq_level - parent_level)                
                    if parent_name in self.leaf_count and self.leaf_count[parent_name] >= max_seq_per_rank:
                        continue

            self.leaf_count[parent_name] = self.leaf_count.get(parent_name, 0) + 1

            # all checks succeeded: add the sequence to the tree
            self.add_tree_node(t0, sid, ranks, parent_level)
            seq_ids += [sid]
            added += 1

        if self.config.verbose:
            print "Total nodes in resulting tree: ", added
        
        if self.config.debug:
            reftax_fname = self.config.tmp_fname("%NAME%_mf_unpruned.tre")
            t0.write(outfile=reftax_fname, format=8)

        self.prune_unifu_nodes(t0)
        return t0, seq_ids
