import os
import sys
import json
import logging
import subprocess

from .utils import which


class BCL(object):
    def __init__(self, fcdir, outdir, bcfile, cpu=60, bcl2fastq="", mis=1):
        self.fcdir = fcdir
        self.outdir = outdir
        self.bcfile = bcfile
        self.mis = str(mis)
        self.index = 1
        self.nproc = str(cpu)
        self.bcl2fastq = bcl2fastq or which("bcl2fastq")
        self.samplesheet = os.path.join(self.outdir, "sample-sheet.csv")

    def create_samplesheet(self):
        idx = []
        with open(self.bcfile) as fi:
            for line in fi:
                if not line.strip() or line.startswith("#"):
                    continue
                line = line.split()
                if len(line) == 2:
                    self.index = 1
                    idx.append((line))
                elif len(line) == 3:
                    self.index = 2
                    idx.append((line))
                else:
                    raise IOError("illegal barcode file %s" % self.bcfile)
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir)
        with open(self.samplesheet, "w") as fo:
            if self.index == 1:
                fo.write("[Data]\nSample_ID,Sample_Name,index\n")
            else:
                fo.write("[Data]\nSample_ID,Sample_Name,index,index2\n")
            for line in idx:
                fo.write(",".join([line[0]] + line) + "\n")

    def make_bcl_cmd(self):
        cmd = [self.bcl2fastq, "-o", self.outdir, "-R", self.fcdir,
               "-i", os.path.join(self.fcdir, "Data/Intensities/BaseCalls"),
               "--no-lane-splitting",
               "--barcode-mismatches", self.mis,
               "-p", self.nproc,
               "--sample-sheet", self.samplesheet]
        return " ".join(cmd)

    @staticmethod
    def call(cmd, run=True, verbose=False):
        if not run:
            if verbose:
                print(cmd)
            return
        if verbose:
            print(cmd)
            subprocess.check_call(cmd, shell=True, stdout=sys.stdout,
                                  stderr=sys.stderr)
        else:
            with open(os.devnull, "w") as fo:
                subprocess.check_call(cmd, shell=True, stdout=fo, stderr=fo)

    def stats(self, j):
        with open(j) as fi:
            stat = json.load(fi)
        info = stat["ConversionResults"][0]
        out = []
        for i in info["DemuxResults"]:
            sn = i["SampleName"]
            rdn = i["NumberReads"]
            out.append((sn, rdn))
        out.append(("Unknow", info["Undetermined"]["NumberReads"]))
        return out

    def run(self):
        self.create_samplesheet()
        cmd = self.make_bcl_cmd()
        self.cmd = cmd
        # self.logs.info(cmd)
        self.logs.info("start bcl2fastq.")
        self.call(cmd, run=True)
        #sms = self.stats(os.path.join(self.outdir, "Stats/Stats.json"))
        # self.logs.info(sms)
        # self.logs.info("Success")

    @property
    def logs(self):
        return logging.getLogger()
