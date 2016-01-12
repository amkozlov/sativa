#!/usr/bin/env python

import os
import sys
import glob
import shutil
import datetime
import random
import re
from subprocess import call,STDOUT
from json_util import EpaJsonParser

class FileUtils:

    @staticmethod    
    def normalize_dir(dir_str):
        if dir_str and not dir_str.endswith("/"):
           dir_str += "/"
        return dir_str

    @staticmethod    
    def remove_if_exists(fname):
        try:
            os.remove(fname)
        except:
            pass

    @staticmethod    
    def rebase(fname, old_basedir, new_basedir):
        return fname.replace(old_basedir, new_basedir)    

class RaxmlWrapper:

    def __init__(self, config): 
        self.cfg = config
        if config.rand_seed:
            random.seed(config.rand_seed)
    
    def make_raxml_fname(self, stem, job_name, absolute=True):
        fname = "RAxML_" + stem + "." + job_name
        if absolute:
            return os.path.join(self.cfg.raxml_outdir, fname)
        else:
            return fname            

    def make_raxml_wildcard(self, job_name):
        return self.make_raxml_fname("*", job_name)

    def cleanup(self, job_name, remove_jplace=False):
        raxml_out_mask = self.make_raxml_wildcard(job_name)
        for fl in glob.glob(raxml_out_mask):
            os.remove(fl)
        if remove_jplace:
            jplace_fname = self.make_raxml_fname("portableTree", job_name) + ".jplace"
            FileUtils.remove_if_exists(jplace_fname)

    def reduce_alignment(self, align_fname, job_name="reduce"):
        reduced_fname = align_fname + ".reduced"
        # we don't have to do anything in restart mode
        if self.cfg.restart and os.path.isfile(reduced_fname):
            return reduced_fname
        else:
            FileUtils.remove_if_exists(reduced_fname)
            raxml_params = ["-f", "c", "-s", align_fname]
            raxml_params += ["--no-dup-check"]
            self.run(job_name, raxml_params)
            self.cleanup(job_name)
            if os.path.isfile(reduced_fname):
                return reduced_fname
            else:
                return align_fname

    def run_epa(self, job_name, align_fname, reftree_fname, optmod_fname="", silent=True, mode="epa", subtree_fname=None,\
    lhw_acc_threshold=0.999):
        raxml_params = ["-s", align_fname, "-t", reftree_fname]
        # assume that by the time we call EPA reference has been cleaned already (e.g. with previous reduce_alignment call)
        raxml_params += ["--no-seq-check"]
        if mode == "l1o_seq":
            raxml_params += ["-f", "O"]
            result_file_stem = "leaveOneOutResults"
        elif mode == "l1o_subtree":
            raxml_params += ["-f", "P", "-z", subtree_fname]
            result_file_stem = "subtreePlacement"
        elif mode == "epa":
            raxml_params += ["-f", "v"]
            result_file_stem = "portableTree"
        elif mode == "epa_mp":
            raxml_params += ["-f", "y"]
            result_file_stem = "portableTree"
        else:
            print "ERROR: Invalid RAxML-EPA running mode: %s" % mode
            sys.exit()
            
        raxml_params += ["--epa-accumulated-threshold", str(lhw_acc_threshold)]

        if self.cfg.epa_use_heuristic in ["TRUE", "YES", "1"]:
            raxml_params += ["-G", str(self.cfg.epa_heur_rate)]            

        if self.cfg.epa_load_optmod and optmod_fname:
            if os.path.isfile(optmod_fname):
                raxml_params += ["-R", optmod_fname]
                if self.cfg.raxml_model.startswith("GTRCAT") and not self.cfg.compress_patterns:
                    raxml_params +=  ["-H"]
            else:
                print "WARNING: Binary model file not found: %s" % optmod_fname
                print "WARNING: Model parameters will be estimated by RAxML"
                
        self.run(job_name, raxml_params, silent)
        
        jp = None
        failed = False
        if mode == "l1o_subtree":
            jp = []
            i = 0
            while True:
                jp_fname = self.make_raxml_fname(result_file_stem, job_name) + ".%d.jplace" % (i+1)
                if not os.path.isfile(jp_fname):
                    break
                jp.append(EpaJsonParser(jp_fname))
                i += 1
            failed = i == 0    
        else:
            jp_fname = self.make_raxml_fname(result_file_stem, job_name) + ".jplace"
            if os.path.isfile(jp_fname):
                jp = EpaJsonParser(jp_fname)
            else:
                failed = True
            
        if failed:
            print "RAxML EPA run failed, please examine the log for details:\n %s" \
                    % self.make_raxml_fname("output", job_name)
            sys.exit()
        else:        
            return jp

    def run(self, job_name, params, silent=True, chkpoint_fname=None):
        if self.cfg.raxml_model == "AUTO":
            print "ERROR: you should have called EpacConfig.resolve_auto_settings() in your script!\n"
            sys.exit()

        self.cleanup(job_name)
        
        lparams  = []
        lparams += params
        lparams += ["-m", self.cfg.raxml_model, "-n", job_name]

        if not self.cfg.use_bfgs:
            lparams += ["--no-bfgs"]
            
        if self.cfg.save_memory:
            lparams += ["-U"]
        
        if self.cfg.verbose:
            lparams += ["--verbose"]
        
        if not "-p" in lparams:
            seed = random.randint(1, 32000)
            lparams += ["-p", str(seed)]
            
        if chkpoint_fname:
            lparams += ["-Z", chkpoint_fname]

        if self.cfg.run_on_cluster:
            self.run_cluster(lparams)
            return;        

        if self.cfg.raxml_remote_call:
            call_str = ["ssh", self.cfg.raxml_remote_host]
        else:
            call_str = []
        call_str += self.cfg.raxml_cmd + lparams
        if silent:        
            self.cfg.log.debug(' '.join(call_str) + "\n")
            out_fname = self.make_raxml_fname("output", job_name)
            with open(out_fname, "w") as fout:
                call(call_str, stdout=fout, stderr=STDOUT)
        else:        
            call(call_str)

        return ' '.join(call_str)
        
    def run_multiple(self, job_name, params, repnum, silent=True):    
        best_lh = float("-inf")
        best_jobname = None
        check_old_jobs = self.cfg.restart
        
        for i in range(repnum):
            call_raxml = True
            chkpoint_fname = None

            rep_jobname = "%s.%d" % (job_name, i)
            
            if check_old_jobs:
                # in resume mode, we have to check where we have stopped before 
                next_jobname = "%s.%d" % (job_name, i+1)
                next_info = self.info_fname(next_jobname)
                # if RAxML_info file for the next job exists, current job has had finished -> skip it
                if os.path.isfile(next_info):
                    call_raxml = False
                else:
                    # use RAxML checkpoints if there are any
                    old_chkpoint_fname = self.checkpoint_fname(rep_jobname)
                    if not os.path.isfile(old_chkpoint_fname):
                        old_chkpoint_fname = self.bkup_checkpoint_fname(rep_jobname)
                        if not os.path.isfile(old_chkpoint_fname):
                            old_chkpoint_fname = None

                    FileUtils.remove_if_exists(next_info)
                    if old_chkpoint_fname:
                        chkpoint_fname = self.checkpoint_fname("last_chkpoint")
                        shutil.move(old_chkpoint_fname, chkpoint_fname)
                    
                    # this was the last RAxML run from previous SATIVA invocation -> proceed without checks from now on
                    check_old_jobs = False
                  
            if call_raxml:
                invoc_str = self.run(rep_jobname, params, silent, chkpoint_fname)
            lh = self.get_tree_lh(rep_jobname, "GAMMA")
            if lh > best_lh:
                best_lh = lh
                best_jobname = rep_jobname
            self.cfg.log.debug("Tree %d GAMMA-based logLH: %s\n" % (i, str(lh)))
        
        best_fname = self.info_fname(best_jobname)
        dst_fname = self.info_fname(job_name)
        shutil.copy(best_fname, dst_fname)

        best_fname = self.result_fname(best_jobname)
        dst_fname = self.result_fname(job_name)
        shutil.copy(best_fname, dst_fname)
        
        return invoc_str
        
    def run_cluster(self, params):
        if self.cfg.raxml_remote_call:
            qsub_call_str = ["ssh", self.cfg.raxml_remote_host]
        else:
            qsub_call_str = []
        
        raxml_call_cmd = self.cfg.raxml_cmd + params        
        for i in range(len(raxml_call_cmd)):
            if isinstance(raxml_call_cmd[i], basestring):
                raxml_call_cmd[i] = FileUtils.rebase(raxml_call_cmd[i], self.cfg.epac_home, self.cfg.cluster_epac_home)
        raxml_call_str = ' '.join(raxml_call_cmd)
                
        script_fname = self.cfg.tmp_fname("%NAME%_sub.sh")
        FileUtils.remove_if_exists(script_fname)
        shutil.copy(self.cfg.cluster_qsub_script, script_fname)
        qsub_job_name = "epa"        
        with open(script_fname, "a") as fout:
            fout.write("#$ -N %s\n" % qsub_job_name)
            fout.write("\n")            
            fout.write(raxml_call_str + "\n")

        cluster_script_fname = FileUtils.rebase(script_fname, self.cfg.epac_home, self.cfg.cluster_epac_home)
        qsub_call_str += ["qsub", "-sync", "y", cluster_script_fname]

        print raxml_call_str + "\n"
        print ' '.join(qsub_call_str) + "\n"
#        sys.exit()

        call(qsub_call_str)
        if not self.cfg.debug:
            FileUtils.remove_if_exists(script_fname)
            
    def get_tree_lh(self, job_name, ratehet="GAMMA"):
        info_fname = self.info_fname(job_name)
        with open(info_fname, "r") as info_file:
            info_str = info_file.read()
        
        lh_patterns = [ "Final %s-based Score of best tree " % ratehet,
                       "Final %s  likelihood: " % ratehet,
                       "%s-based likelihood " % ratehet]
        
        m = None
        for pat in lh_patterns:
            m = re.search('(?<=%s)[0-9.\-]+' % pat, info_str)
            if m:
                lh = float(m.group(0))
                return lh
                
        return None
        
    def get_invocation_str(self, job_name):
        info_fname = self.info_fname(job_name)
        with open(info_fname, "r") as info_file:
            info_str = info_file.read()
        
        pattern = "RAxML was called as follows:\n\n"
        
        m = re.search('(?<=%s).*' % pattern, info_str)
        if m:
            return m.group(0)
        else:
            return ""
    
    def result_fname(self, job_name):
        return self.make_raxml_fname("result", job_name)
    
    def besttree_fname(self, job_name):
        return self.make_raxml_fname("bestTree", job_name)

    def info_fname(self, job_name):
        return self.make_raxml_fname("info", job_name)

    def checkpoint_fname(self, job_name):
        return self.make_raxml_fname("binaryCheckpoint", job_name)

    def bkup_checkpoint_fname(self, job_name):
        return self.make_raxml_fname("binaryCheckpointBackup", job_name)

    def result_exists(self, job_name):
        if os.path.isfile(self.result_fname(job_name)):
            return True
        else:
            return False

    def besttree_exists(self, job_name):
        if os.path.isfile(self.besttree_fname(job_name)):
            return True
        else:
            return False

    def epa_result_exists(self, job_name):
        if os.path.isfile(self.make_raxml_fname("labelledTree", job_name)):
            return True
        else:
            return False

    def copy_result_tree(self, job_name, dst_fname):
        src_fname = self.result_fname(job_name)
        shutil.copy(src_fname, dst_fname)

    def copy_best_tree(self, job_name, dst_fname):
        src_fname = self.besttree_fname(job_name)
        shutil.copy(src_fname, dst_fname)

    def copy_optmod_params(self, job_name, dst_fname):
        src_fname = self.make_raxml_fname("binaryModelParameters", job_name)
        shutil.copy(src_fname, dst_fname)

    def copy_epa_orig_tree(self, job_name, dst_fname):
        src_fname = self.make_raxml_fname("originalLabelledTree", job_name)
        shutil.copy(src_fname, dst_fname)
        
    def copy_epa_result_tree(self, job_name, dst_fname):
        src_fname = self.make_raxml_fname("labelledTree", job_name)
        shutil.copy(src_fname, dst_fname)

    def copy_epa_jplace(self, job_name, dst_fname, mode="epa", move=False):
        if mode == "l1o_seq":
            result_file_stem = "leaveOneOutResults"
        elif mode == "l1o_subtree":
            result_file_stem = "subtreePlacement"
        elif mode == "epa":
            result_file_stem = "portableTree"
        else:
            return
        src_fname = self.make_raxml_fname(result_file_stem, job_name) + ".jplace"
        if move:
            shutil.move(src_fname, dst_fname)
        else:
            shutil.copy(src_fname, dst_fname)
