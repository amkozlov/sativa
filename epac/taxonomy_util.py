#!/usr/bin/env python

#from epac.ete2 import Tree
from ete2 import Tree

class Taxonomy:
    EMPTY_RANK = "-"
    RANK_UID_DELIM = "@@"
    
    def __init__(self):
        tree_nodes = []

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

    def get_seq_ranks(self, seq_id):
        return []

    def seq_count(self):
        return 0

    def max_rank_level(self):
        return 0

    def items(self):
        return 0

class GGTaxonomyFile(Taxonomy):
    rank_placeholders = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]    
    
    @staticmethod
    def strip_prefix(ranks):
        new_ranks = ['']*len(ranks);
        for i in range(len(ranks)):
            if ranks[i].startswith(GGTaxonomyFile.rank_placeholders[i]):
                plen = len(GGTaxonomyFile.rank_placeholders[i])
                new_ranks[i] = ranks[i][plen:]
            else:
                new_ranks[i] = ranks[i]
        return new_ranks

    @staticmethod
    def add_prefix(ranks):
        new_ranks = ['']*len(ranks);
        for i in range(len(ranks)):
            new_ranks[i] = GGTaxonomyFile.add_rank_prefix(ranks[i], i)
        return new_ranks

    @staticmethod
    def add_rank_prefix(rank_name, rank_level):
        prefix = GGTaxonomyFile.rank_placeholders[rank_level]
        if rank_name.startswith(prefix):
            return rank_name
        else:
            return prefix + rank_name


    @staticmethod
    def lineage_str(ranks, strip_prefix=False):
        if strip_prefix:
            ranks = GGTaxonomyFile.strip_prefix(ranks)
        return Taxonomy.lineage_str(ranks);


    def __init__(self, tax_fname, prefix=""):
        self.tax_fname = tax_fname        
        self.prefix = prefix
        self.seq_ranks_map = {}
        self.common_ranks = set([])
        self.load_taxonomy()

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

    # zero-based! 
    def max_rank_level(self):
        return 6   # level of species in standard 7-level taxonomy
    
    def rank_level_name(self, rank_level):
        return { 0: "Kingdom",
                 1: "Phylum",
                 2: "Class",
                 3: "Order",
                 4: "Family",
                 5: "Genus",
                 6: "Species"
                }[rank_level]

    def get_seq_ranks(self, seq_id):
        return self.seq_ranks_map[seq_id]

    def lineage_str(self, seq_id, use_placeholders=False):
        ranks = list(self.seq_ranks_map[seq_id])
        if use_placeholders:        
            for i in range(len(ranks)):
                if ranks[i] == Taxonomy.EMPTY_RANK:
                    ranks[i] = GGTaxonomyFile.rank_placeholders[i]
        return Taxonomy.lineage_str(ranks)

    def lowest_assigned_rank_level(self, seq_id):
        ranks = self.seq_ranks_map[seq_id]
        return Taxonomy.lowest_assigned_rank_level(ranks)

    def items(self):
        return self.seq_ranks_map.items()
        
    def iteritems(self):
        return self.seq_ranks_map.items()

    def get_map(self):
        return self.seq_ranks_map

    def remove_seq(self, seqid):
        del self.seq_ranks_map[seqid]

    def make_binomial_name(self, ranks):
        last_rank = len(ranks)-1
        if ranks[last_rank] != Taxonomy.EMPTY_RANK:
            genus_name = ranks[last_rank-1][3:]
            sp_name = ranks[last_rank][3:]
            prefix = ranks[last_rank][:3]
            if not sp_name.startswith(genus_name):
                ranks[last_rank] = prefix + genus_name + "_" + sp_name

    def normalize_rank_name(self, rank, rank_name):
        invalid_chars = ['(', ')', ',', ';', ':']
        for ch in invalid_chars:
            rank_name = rank_name.replace(ch, "_")
#        rank_prefix = self.rank_placeholders[rank] if rank < 7 else ""
#        if not rank_name.startswith(rank_prefix):
#            rank_name = rank_prefix + rank_name
        return rank_name

    def normalize_rank_names(self):
        for sid, ranks in self.seq_ranks_map.iteritems():
            for i in range(1, len(ranks)):
                ranks[i] = self.normalize_rank_name(i, ranks[i])
#            self.make_binomial_name(ranks);

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
            
    def load_taxonomy(self):
        fin = open(self.tax_fname)
        for line in fin:
            line = line.strip()
            toks = line.split("\t")
            sid = self.prefix + toks[0]
            ranks_str = toks[1]
            ranks = ranks_str.split(";")
            for i in range(len(ranks)):
                rank_name = ranks[i].strip()
                if rank_name in GGTaxonomyFile.rank_placeholders:
                    rank_name = Taxonomy.EMPTY_RANK
                ranks[i] = rank_name
                
            if len(ranks) < 7:
                ranks += [Taxonomy.EMPTY_RANK] * (7 - len(ranks))
#                print "WARNING: sequence " + sid + " has incomplete taxonomic annotation. Missing ranks were considered empty (%s)" % Taxonomy.EMPTY_RANK 

            self.seq_ranks_map[sid] = ranks     

        fin.close()

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
        tax_seq_level = self.taxonomy.max_rank_level() + 1
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

            parent_level = tax_seq_level - 1            
            while ranks[parent_level] == Taxonomy.EMPTY_RANK:
                parent_level -= 1
            parent_name = Taxonomy.get_rank_uid(ranks, parent_level)
            if parent_name in self.tree_nodes:
                parent_node = self.tree_nodes[parent_name]
                # filter by max number of seqs (threshold depends from rank level, 
                # i.e. for genus there can be more seqs than for species)
                max_seq_per_rank = max_seqs_per_leaf * (tax_seq_level - parent_level)                
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
  
if __name__ == "__main__":
    gt = GGTaxonomyFile("/home/zhangje/GIT/EPA-classifier/example/training_tax.txt")
    print gt.get_common_ranks()
