#! /usr/bin/env python
import math
from ete2 import Tree

class erlang:
    def __init__(self):
        pass
    
    def one_tail_test(self, rate, k, x):
        """rate: estimated branching rate from reference tree
           k: node height
           x: placement branch length"""
        p = 0.0
        for n in range(k):
            p = p + (1.0/float(math.factorial(n))) * math.exp((-rate)*x) * math.pow(rate*x, n)
        return p

class tree_param:
    def __init__(self, tree, origin_taxonomy):
        """tree: rooted and branch labled tree in newick format
            origin_taxonomy: a dictionary of leaf name and taxonomy ranks"""
        self.tree = tree 
        self.taxonomy = origin_taxonomy 
    
    def get_speciation_rate(self):
        #pruning the input tree such that each species only appear once
        species = set()
        keepseqs = []
        for name in self.taxonomy.keys():
            ranks = self.taxonomy[name]
            sp = ranks[-1]
            if sp == "-":
                keepseqs.append(name)
            else:
                if not sp in species:
                    keepseqs.append(name)
                    species.add(sp)
        root = Tree(self.tree)
        root.prune(keepseqs, preserve_branch_length=True)
        sumbr = 0.0
        cnt = 0.0 
        for node in root.traverse(strategy = "preorder"):
            sumbr = sumbr + node.dist
            cnt = cnt + 1.0
        return float(cnt) / float(sumbr)
       
    def get_speciation_rate_fast(self):
        """ETE2 prune() function is extremely slow on large trees, so
        this function don't use it and instead just removes "redundant"
        species-level nodes one-by-one"""

        species = set()
        root = Tree(self.tree)

        name2node = {}
        for node in root.traverse(strategy = "postorder"):
          if node.is_leaf():
              name2node[node.name] = node

        #pruning the input tree such that each species only appear once
        for name in self.taxonomy.keys():
            ranks = self.taxonomy[name]
            sp = ranks[-1]
            if sp != "-":
                if sp in species:
                    node = name2node.get(name, None)
                    if node:
                        node.delete(preserve_branch_length=True)
                    else:
                        raise ValueError("Node names not found in the tree: " + name)
                else:
                    species.add(sp)

        # traverse the pruned tree, counting the number of speciation events and 
        # summing up the branch lengths
        sumbr = 0.0
        cnt = 0
        for node in root.traverse(strategy = "preorder"):
            sumbr += node.dist
            cnt += 1

        # sp_rate = number_of_sp_events / sum_of_branch_lengts
        return float(cnt) / float(sumbr)

    def get_nodesheight(self):
        root = Tree(self.tree)
        nh_map = {}
        for node in root.traverse(strategy = "preorder"):
            if hasattr(node, "B"):
                height = node.get_closest_leaf(topology_only=True)
                #height = node.get_farthest_leaf(topology_only=True)
                nh_map[node.B] = height[1] + 1
        
        return nh_map

if __name__ == "__main__":
    print("This is erlang.py main")
    el = erlang()
    print el.one_tail_test(rate = 17, k = 1, x = 0.221977)
