#! /usr/bin/env python
from epac.taxonomy_util import Taxonomy
from epac.erlang import erlang
import math

class TaxTreeHelper:
    def __init__(self, tax_map, cfg):
        self.origin_taxonomy = tax_map
        self.cfg = cfg
        self.outgroup = None
        self.mf_rooted_tree = None
        self.bf_rooted_tree = None
        self.tax_tree = None
        self.bid_taxonomy_map = None
    
    def set_mf_rooted_tree(self, rt):
        self.mf_rooted_tree = rt
        self.save_outgroup()
        self.bf_rooted_tree = None
        self.tax_tree = None
        self.bid_taxonomy_map = None
    
    def set_bf_unrooted_tree(self, ut):
        self.restore_rooting(ut)
        
    def get_outgroup(self):
        return self.outgroup

    def get_tax_tree(self):
        if not self.tax_tree:
            self.label_bf_tree_with_ranks()
        return self.tax_tree
    
    def get_bid_taxonomy_map(self):
        self.get_tax_tree()
        if not self.bid_taxonomy_map:
            self.build_bid_taxonomy_map()
        return self.bid_taxonomy_map
    
    def save_outgroup(self):
        rt = self.mf_rooted_tree
        
        # remove unifurcation at the root
        if len(rt.children) == 1:
            rt = rt.children[0]

        if len(rt.children) > 1:
            outgr = rt.children[0]    
            outgr_size = len(outgr.get_leaves())
            for child in rt.children:
                if child != outgr:
                    child_size = len(child.get_leaves())
                    if child_size < outgr_size:
                        outgr = child
                        outgr_size = child_size
        else:
            raise AssertionError("Invalid tree: unifurcation at the root node!")
        
        self.outgroup = outgr
    
    def restore_rooting(self, utree):
        outgr_leaves = self.outgroup.get_leaf_names()
        # check if outgroup consists of a single node - ETE considers it to be root, not leaf
        if not outgr_leaves:
            outgr_root = utree&outgr.name
        elif len(outgr_leaves) == 1:
            outgr_root = utree&outgr_leaves[0]
        else:
            # Even unrooted tree is "implicitely" rooted in ETE representation.
            # If this pseudo-rooting happens to be within the outgroup, it cause problems
            # in the get_common_ancestor() step (since common_ancestor = "root")
            # Workaround: explicitely root the tree outside from outgroup subtree
            for node in utree.iter_leaves():
                if not node.name in outgr_leaves:
                    tmp_root = node.up
                    if not utree == tmp_root:
                        utree.set_outgroup(tmp_root)
                        break
            
            outgr_root = utree.get_common_ancestor(outgr_leaves)

        # we could be so lucky that the RAxML tree is already correctly rooted :)
        if outgr_root != utree:
            utree.set_outgroup(outgr_root)

        self.bf_rooted_tree = utree

    def label_bf_tree_with_ranks(self):
        """labeling inner tree nodes with taxonomic ranks"""
        if not self.bf_rooted_tree:
            raise AssertionError("self.bf_rooted_tree is not set: TaxTreeHelper.set_bf_unrooted_tree() must be called before!")
            
        for node in self.bf_rooted_tree.traverse("postorder"):
            if node.is_leaf():
                seq_ranks = self.origin_taxonomy[node.name]
                rank_level = Taxonomy.lowest_assigned_rank_level(seq_ranks)
                node.add_feature("rank_level", rank_level)
                node.add_feature("ranks", seq_ranks)
                node.name += "__" + seq_ranks[rank_level]
            else:
                if len(node.children) != 2:
                    raise AssertionError("FATAL ERROR: tree is not bifurcating!")
                lchild = node.children[0]
                rchild = node.children[1]
                rank_level = min(lchild.rank_level, rchild.rank_level)
                while rank_level >= 0 and lchild.ranks[rank_level] != rchild.ranks[rank_level]:
                    rank_level -= 1
                node.add_feature("rank_level", rank_level)
                node_ranks = [Taxonomy.EMPTY_RANK] * 7
                if rank_level >= 0:
                    node_ranks[0:rank_level+1] = lchild.ranks[0:rank_level+1]
                    node.name = lchild.ranks[rank_level]
                else:
                    node.name = "Undefined"
                    if hasattr(node, "B") and self.cfg.verbose:
                        print "INFO: no taxonomic annotation for branch %s (reason: children belong to different kingdoms)" % node.B

                node.add_feature("ranks", node_ranks)
        
        self.tax_tree = self.bf_rooted_tree

    def build_bid_taxonomy_map(self):
        self.bid_taxonomy_map = {}
        for node in self.tax_tree.traverse("postorder"):
            if not node.is_root() and hasattr(node, "B"):                
                parent = node.up                
                self.bid_taxonomy_map[node.B] = parent.ranks
#                bid_taxonomy_map[node.B] = node.ranks

    
class TaxClassifyHelper:
    def __init__(self, cfg, bid_taxonomy_map, brlen_pv = 0., sp_rate = 0., node_height = []):
        self.cfg = cfg
        self.bid_taxonomy_map = bid_taxonomy_map
        self.brlen_pv = brlen_pv
        self.sp_rate = sp_rate
        self.node_height = node_height
        self.erlang = erlang()

    def classify_seq(self, edges, method = "1", minlw = 0.):
        if len(edges) > 0:
            edges = self.erlang_filter(edges)
            if method == "1":
                ranks, lws = self.assign_taxonomy_maxsum(edges, minlw)
            else:
                ranks, lws = self.assign_taxonomy_maxlh(edges)
            return ranks, lws
        else:
            return None, None      
            
    def erlang_filter(self, edges):
        if self.brlen_pv == 0.:
            return edges
            
        newedges = []
        for edge in edges:
            edge_nr = str(edge[0])
            pendant_length = edge[4]
            pv = self.erlang.one_tail_test(rate = self.sp_rate, k = int(self.node_height[edge_nr]), x = pendant_length)
            if pv >= self.brlen_pv:
                newedges.append(edge)
        
        if len(newedges) == 0:
            return newedges
        
        # adjust likelihood weights -> is there a better way ???        
        sum_lh = 0
        max_lh = float(newedges[0][1])
        for edge in newedges:
            lh = float(edge[1])
            sum_lh += math.exp(lh - max_lh)
                        
        for edge in newedges:
            lh = float(edge[1])
            edge[2] = math.exp(lh - max_lh) / sum_lh

        return newedges
     
    def assign_taxonomy_maxlh(self, edges):
        #Calculate the sum of likelihood weight for each rank
        taxonmy_sumlw_map = {}
        for edge in edges:
            edge_nr = str(edge[0])
            lw = edge[2]
            taxonomy = self.bid_taxonomy_map[edge_nr]
            for rank in taxonomy:
                if rank == "-":
                    taxonmy_sumlw_map[rank] = -1
                elif rank in taxonmy_sumlw_map:
                    oldlw = taxonmy_sumlw_map[rank]
                    taxonmy_sumlw_map[rank] = oldlw + lw
                else:
                    taxonmy_sumlw_map[rank] = lw
        
        #Assignment using the max likelihood placement
        ml_edge = edges[0]
        edge_nr = str(ml_edge[0])
        maxlw = ml_edge[2]
        ml_ranks = self.bid_taxonomy_map[edge_nr]
        ml_ranks_copy = []
        for rk in ml_ranks:
            ml_ranks_copy.append(rk)
        lws = []
        cnt = 0
        for rank in ml_ranks:
            lw = taxonmy_sumlw_map[rank]
            if lw > 1.0:
                lw = 1.0
            lws.append(lw)
            if rank == "-" and cnt > 0 :                
                for edge in edges[1:]:
                    edge_nr = str(edge[0])
                    taxonomy = self.bid_taxonomy_map[edge_nr]
                    newrank = taxonomy[cnt]
                    newlw = taxonmy_sumlw_map[newrank]
                    higherrank_old = ml_ranks[cnt -1]
                    higherrank_new = taxonomy[cnt -1]
                    if higherrank_old == higherrank_new and newrank!="-":
                        ml_ranks_copy[cnt] = newrank
                        lws[cnt] = newlw
            cnt = cnt + 1
            
        return ml_ranks_copy, lws

    def assign_taxonomy_maxsum(self, edges, minlw):
        """this function sums up all LH-weights for each rank and takes the rank with the max. sum """
        # in EPA result, each placement(=branch) has a "weight"
        # since we are interested in taxonomic placement, we do not care about branch vs. branch comparisons,
        # but only consider rank vs. rank (e. g. G1 S1 vs. G1 S2 vs. G1)
        # Thus we accumulate weights for each rank, there are to measures:
        # "own" weight  = sum of weight of all placements EXACTLY to this rank (e.g. for G1: G1 only)
        # "total" rank  = own rank + own rank of all children (for G1: G1 or G1 S1 or G1 S2)
        rw_own = {}
        rw_total = {}
        rb = {}
        
        ranks = [Taxonomy.EMPTY_RANK]
        
        for edge in edges:
            br_id = str(edge[0])
            lweight = edge[2]
            lowest_rank = None

            if lweight == 0.:
                continue
            
            # accumulate weight for the current sequence                
            ranks = self.bid_taxonomy_map[br_id]
            for i in range(len(ranks)):
                rank = ranks[i]
                rank_id = Taxonomy.get_rank_uid(ranks, i)
                if rank != Taxonomy.EMPTY_RANK:
                    rw_total[rank_id] = rw_total.get(rank_id, 0) + lweight
                    lowest_rank = rank_id
                    if not rank_id in rb:
                        rb[rank_id] = br_id
                else:
                    break

            if lowest_rank:
                rw_own[lowest_rank] = rw_own.get(lowest_rank, 0) + lweight
                rb[lowest_rank] = br_id
            elif self.cfg.verbose:
                print "WARNING: no annotation for branch ", br_id
            
        # if all branches have empty ranks only, just return this placement
        if len(rw_total) == 0:
            return ranks, [1.] * len(ranks)
        
        # we assign the sequence to a rank, which has the max "own" weight AND 
        # whose "total" weight is greater than a confidence threshold
        max_rw = 0.
        s_r = None
        for r in rw_own.iterkeys():
            if rw_own[r] > max_rw and rw_total[r] >= minlw:
                s_r = r
                max_rw = rw_own[r] 
        if not s_r:
            s_r = max(rw_total.iterkeys(), key=(lambda key: rw_total[key]))

        a_br_id = rb[s_r]
        a_ranks = self.bid_taxonomy_map[a_br_id]

        # "total" weight is considered as confidence value for now
        a_conf = [0.] * len(a_ranks)
        for i in range(len(a_conf)):
            rank = a_ranks[i]
            if rank != Taxonomy.EMPTY_RANK:
                rank_id = Taxonomy.get_rank_uid(a_ranks, i)
                a_conf[i] = rw_total[rank_id]

        return a_ranks, a_conf
    
