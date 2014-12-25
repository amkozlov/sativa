#! /usr/bin/env python
import sys
import os
import json
import operator
import time
from ete2 import Tree, SeqGroup
from config import EpacConfig
from subprocess import call

class hmmer:
    def __init__(self, config, refalign = None, query = None, refprofile = None, discard = None, seqs = None, minp = 0.9):
        self.cfg = config
        self.refalign = refalign
        self.query = query
        self.refprofile = refprofile
        self.hmmbuildpath = config.hmmer_home + "/hmmbuild"
        self.hmmalignpath = config.hmmer_home + "/hmmalign"
        if self.refprofile == None:
            self.refprofile = config.tmp_fname("%NAME%.hmm")
        self.stockname = config.tmp_fname("%NAME%.stock")
        self.trimed = config.tmp_fname("%NAME%.trimed.afa")
        self.output = config.tmp_fname("%NAME%.aligned.afa")
        self.merged = config.tmp_fname("%NAME%.merged.afa")
        self.discardpath = discard
        self.seqs = seqs 
        self.minp = minp
        self.minl = 100
    
    def _remove(self, filename):
        if os.path.exists(filename):
            os.remove(filename)
    
    def __del__(self):
        self._remove(self.stockname)
        self._remove(self.trimed)
        self._remove(self.output)
    
    def __processHMMseq(self, seqin):
        newseq = ""
        for s in seqin:
            if s == ".":
                pass
            elif s == "-":
                newseq = newseq + s
            elif s.isupper():
                newseq = newseq + s
        return newseq
    
    def build_hmm_profile(self, informat="afa"):
        #hmmbuild --informat afa refotu.hmm ref_outs_547.fas
        call_str = [self.hmmbuildpath, "--symfrac", "0.0", "--informat", informat, self.refprofile, self.refalign]
        if self.cfg.verbose:
            print "\n" + ' '.join(call_str) + "\n"
        call(call_str) #, stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
        return self.refprofile

    def hmm_align(self):
        #hmmalign -o 454.stock refotu.hmm 454input.fna.min100.fasta
        call_str = [self.hmmalignpath,"-o", self.stockname, self.refprofile, self.query]
        if self.cfg.verbose:
            print "\n" + ' '.join(call_str) + "\n"
        call(call_str) #, stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
        return self.stockname
    
    def get_hmm_refalignment(self):
        sites = []
        hmp = open(self.refprofile)
        l = hmp.readline()
        start = False
        while l!="":
            if l.startswith("//"):
                break
            if start:
                l = l.strip()
                ll = l.split()
                usedsite = int(ll[5])
                sites.append(usedsite)
                l = hmp.readline()
                l = hmp.readline()
            else:
                if l.startswith("HMM "):
                    start = True
                    l = hmp.readline()
                    l = hmp.readline()
                    l = hmp.readline()
                    l = hmp.readline()
            l = hmp.readline()
        hmp.close()
        align = SeqGroup(self.refalign)
        fout = open(self.trimed, "w")
        for entr in align.get_entries():
            fout.write(">" + entr[0] + "\n")
            for pos in sites:
                fout.write(entr[1][pos-1])
            fout.write("\n")
        fout.close()
        return self.trimed, len(sites)

    def parse_HMM(self, l_ref):
        """stock format"""
        cnt = 0
        fin = open(self.stockname)
        line = fin.readline()
        seqs = {}
        while line!="":
            if line.startswith("//"):
                break
            elif line.startswith("#=GS"):
                cnt = 0
                pass
            elif line.startswith("#"):
                pass
            elif line.startswith("\n"):
                cnt = cnt + 1
            else:
                line = line.strip()
                if cnt == 1:
                    l2 = line.split()
                    ss = self.__processHMMseq(l2[1])
                    seqs[l2[0]] = ss
                else:
                    l2 = line.split()
                    seq = seqs[l2[0]]
                    ss = self.__processHMMseq(l2[1])
                    seqs[l2[0]] = seq + ss 
            line = fin.readline()
        fin.close()
        fout = open(self.output, "w")
        foutdiscard = open(self.discardpath, "w")
        for key in seqs.keys():
            
            #key is the sequence name which is the id
            numleft = count_non_gap(seqs[key])
            numall = len(self.seqs.get_seq(key))
            
            if numall > 0:
                pleft = float(numleft) / float(numall)
            else:
                pleft = 0
            
            if numleft >= self.minl and pleft >= self.minp:
                fout.write(">" + key + "\n")
                seq = seqs[key]
                lappd = l_ref - len(seq)
                if lappd > 0:
                    appd = "-" * lappd
                    seq = seq + appd
                elif lappd < 0:
                    print("Warning: query sequence > ref sequence")
            
                fout.write(seq + "\n")
            else:
                foutdiscard.write(">" + key + "\n")
        fout.close()
        return self.output

    def hmm_alignment(self, ref_align, query, outfolder, lmin = 100):
        if not os.path.exists(self.refprofile):
            self.build_hmm_profile()
        self.hmm_align()
        final_ref, ref_len = self.trim_refalign_hmm()
        final_query = self.parse_HMM(l_ref = ref_len, minl = lmin)

    def align(self):
        #aquire reference alignment that hmm would use
        refaln, numsite = self.get_hmm_refalignment()
        #alignment
        self.hmm_align()
        #process alignment
        queryaln = self.parse_HMM(l_ref = numsite)
        #merge refrence and query alignment
        merge_alignment(aln1 = refaln, aln2 = queryaln, fout = self.merged, numsites = numsite)
        return self.merged


class muscle:
    def __init__(self, config):
        self.cfg = config 
        self.musclepath = config.muscle_home + "/muscle"
        self.outname = config.tmp_fname("%NAME%.afa")
    
    def merge(self, aln1, aln2):
        #muscle -profile -in1 existing_msa.afa -in2 new_seq.fa -out combined.afa
        call_str = [self.musclepath,"-profile", "-in1", aln1, "-in2", aln2, "-out", self.outname]
        if self.cfg.debug:
            print "\n" + ' '.join(call_str) + "\n"
        call(call_str)
        return self.outname


def merge_alignment(aln1, aln2, fout, numsites):
    seqs1 = SeqGroup(aln1)
    seqs2 = SeqGroup(aln2)
    if len(seqs1) == 0 or len(seqs2) == 0:
        print("No sequences aligned! ")
        sys.exit()
    with open(fout, "w") as fo:
        for seq in seqs1.iter_entries():
            if len(seq[1].strip()) == numsites:
                fo.write(">" + seq[0] + "\n" + seq[1] + "\n")
            else:
                print("Error in alignment ....")
                sys.exit()
        for seq in seqs2.iter_entries():
            if len(seq[1].strip()) == numsites:
                fo.write(">" + seq[0] + "\n" + seq[1] + "\n")
            else:
                print("Error in alignment ....")
                sys.exit()


def count_non_gap(seqin):
    cnt = 0
    for s in seqin:
        if s!="-":
            cnt = cnt + 1
    return cnt


if __name__ == "__main__":
    print("This is main")
    cfg = EpacConfig()
    hm = hmmer(config = cfg, refalign = "example/t1_trimed.fa")
    #trimed = hm.process_ref_alignment()
    pf = hm.build_hmm_profile()
    print(pf)
