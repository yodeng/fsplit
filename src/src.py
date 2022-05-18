#!/usr/bin/env python
# coding:utf-8

import os
import sys
import math
import time
import json
import gzip
import logging
import argparse
import subprocess
import editdistance

from .version import __version__
from . import sh


class Zopen(object):
    def __init__(self, name,  mode="rb", gzip=False):
        self.name = name
        self.mode = mode
        self.gzip = gzip
        self.handler = None

    def __enter__(self):
        if self.name.endswith(".gz"):
            if self.gzip:
                p = subprocess.Popen(
                    ["gzip", "-c", "-d", self.name], stdout=subprocess.PIPE)
                self.handler = p.stdout
            else:
                self.handler = gzip.open(self.name, self.mode)
        else:
            self.handler = open(self.name, self.mode)
        return self.handler

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.handler:
            self.handler.close()


class MultiZipHandle(object):
    def __init__(self, **infiles):
        self.info = infiles
        self.handler = {}

    def __enter__(self):
        for sn, f in self.info.items():
            if f.endswith(".gz"):
                self.handler[sn] = gzip.open(f, "wb")
            else:
                self.handler[sn] = open(f, "wb")
        return self.handler

    def __exit__(self, exc_type, exc_val, exc_tb):
        for sn, h in self.handler.items():
            h.close()


def chunk(lst, size=1, group=None):
    if group:
        size = int(math.ceil(len(lst)/float(group)))
    return [lst[i:i + size] for i in range(0, len(lst), size)]


class FastqIndex(object):

    def __init__(self, fastqfile=""):
        self.fq = os.path.abspath(fastqfile)
        self.fai = self.fq + ".fai"

    def index_pos(self):
        idx = []
        with open(self.fai) as fi:
            h = fi.readline().strip()
            if h.isdigit():
                idx.append(0)
                for line in fi:
                    idx.append(line.strip())
                self.logs.info(
                    "read fastq index (%s) done, %d sequences for input.", self.fai, len(idx))
                return idx
            h = h.split()
            if len(h) == 6:
                idx.append(0)
                e = int(h[5])
                for line in fi:
                    line = line.split()
                    idx.append(e+int(line[4]))
                    e = int(line[5])
        self.logs.info(
            "read fastq index (%s) done, %d sequences for input.", self.fai, len(idx))
        return idx

    def index_crt(self):
        with Zopen(self.fq, gzip=True) as fi, open(self.fai, "w") as fo:
            cur = 0
            for n, line in enumerate(fi):
                if n % 4 == 0:
                    fo.write("%d\n" % cur)
                cur += len(line)
            self.logs.info("create fastq index (%s) done.", self.fai)

    @classmethod
    def createindex(cls, fq=""):
        fqi = cls(fq)
        fqi.index_crt()

    @classmethod
    def getIndexPos(cls, fq=""):
        fqi = cls(fq)
        if os.path.isfile(fqi.fai):
            return fqi.index_pos()
        idx = []
        with Zopen(fqi.fq, gzip=True) as fi, open(fqi.fai, "w") as fo:
            cur = 0
            for n, line in enumerate(fi):
                if n % 4 == 0:
                    fo.write("%d\n" % cur)
                    idx.append(cur)
                cur += len(line)
            fqi.logs.info("create fastq index done.")
        return idx

    @property
    def logs(self):
        return logging.getLogger()


def parseArg():
    subargs = ["index", "split"]
    main_parser = argparse.ArgumentParser(
        description="split a mix fastq by barcode index.",)
    mode_parser = main_parser.add_argument_group("Commands options")
    mode_parser.add_argument("mode",  metavar="{%s}" % ",".join(
        subargs), help="command to run.", choices=subargs)
    args = main_parser.parse_args(sys.argv[1:2])
    mode = args.mode
    if mode == "index":
        des = "index fastq file for reading in multi processing, can be instead by `samtools fqidx <fqfile>`."
    else:
        des = "split sequence data by barcode."
    parser = argparse.ArgumentParser(
        description=des, prog=" ".join(sys.argv[0:2]))
    general_parser = parser.add_argument_group("General options")
    general_parser.add_argument("mode", metavar=mode, choices=subargs)
    general_parser.add_argument("-i", "--input", type=str, help="input fastq file or BCL flowcell directory, required",
                                required=True, metavar="<str>")
    general_parser.add_argument('-v', '--version',
                                action='version', version="v" + __version__)
    if mode == "index":
        parser_index = parser.add_argument_group("Options")
    elif mode == "split":
        parser_split = parser.add_argument_group("Options")
        parser_split.add_argument("-b", "--barcode", type=str,
                                  help='barcode and sample file, required', required=True, metavar="<file>")
        parser_split.add_argument('-m', "--mismatch", help="mismatch allowed for barcode search, 0 by default",
                                  type=int, default=0, metavar="<int>")
        parser_split.add_argument('-t', "--threads", help="threads core, 2 by default",
                                  type=int, default=2, metavar="<int>")
        parser_split.add_argument('-o', "--output", help="output directory, required",
                                  type=str, required=True, metavar="<str>")
        parser_split.add_argument("-d", '--drup',   action='store_true',
                                  help="drup barcode sequence in output if set",  default=False)
        parser_split.add_argument('--bcl2fq', metavar="<str>",
                                  help="bcl2fastq path is necessary, if not set, auto detected")
        parser_split.add_argument("--output-gzip",   action='store_true',
                                  help="gzip output fastq file, this will make your process slower", default=False)
    return parser.parse_args()


def splitFastq(fq, s, e, outQ, barcode, mis=0, drup=False):
    drup_pos = dict.fromkeys(barcode.keys(), 0)
    if drup:
        for bc in barcode:
            drup_pos[bc] = len(bc)
    bc_len = {b: len(b) for b in barcode.keys()}
    with Zopen(fq) as fi:
        fi.seek(s)
        out = {}
        while True:
            if fi.tell() > e:
                break
            name = fi.readline()
            seq = fi.readline()
            flag = fi.readline()
            qual = fi.readline()
            for b in barcode:
                bl = bc_len[b]
                bar = seq[:bl]
                dis = editdistance.eval(b, bar)
                if dis <= mis:
                    dp = drup_pos[b]
                    sn = barcode[b]
                    seq = seq[dp:]
                    qual = qual[dp:]
                    out.setdefault(sn, []).extend((name, seq, flag, qual))
                    break
            # else:
            #    out.setdefault("Unknow", []).extend((name, seq, flag, qual))
        outQ.put(out)


def timeRecord(func):
    def wrapper(*args, **kwargs):
        s = time.time()
        value = func(*args, **kwargs)
        sys.stdout.write("\nTime elapse: %d sec.\n" % int(time.time() - s))
        return value
    return wrapper


def load_bcl_stats(j):
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


def log(logfile=None, level="info"):
    logger = logging.getLogger()
    if level.lower() == "info":
        logger.setLevel(logging.INFO)
        f = logging.Formatter(
            '[%(levelname)s %(asctime)s] %(message)s')
    elif level.lower() == "debug":
        logger.setLevel(logging.DEBUG)
        f = logging.Formatter(
            '[%(levelname)s %(threadName)s %(asctime)s %(funcName)s(%(lineno)d)] %(message)s')
    if logfile is None:
        h = logging.StreamHandler(sys.stdout)
    else:
        h = logging.FileHandler(logfile, mode='w')
    h.setFormatter(f)
    logger.addHandler(h)
    return logger


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
