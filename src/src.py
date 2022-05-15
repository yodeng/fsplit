#!/usr/bin/env python
# coding:utf-8

import os
import sys
import math
import time
import gzip
import argparse
import editdistance

from .version import __version__


class Zopen(object):
    def __init__(self, name,  mode="rb"):
        self.name = name
        self.mode = mode
        self.handler = None

    def __enter__(self):
        if self.name.endswith(".gz"):
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
        self.fq = fastqfile
        self.fai = self.fq + ".fai"

    def index_pos(self):
        idx = []
        with open(self.fai) as fi:
            h = fi.readline().strip()
            if h.isdigit():
                idx.append(0)
                for line in fi:
                    idx.append(line.strip())
                return idx
            h = h.split()
            if len(h) == 6:
                idx.append(0)
                e = int(h[5])
                for line in fi:
                    line = line.split()
                    idx.append(e+int(line[4]))
                    e = int(line[5])
        return idx

    def index_crt(self):
        with Zopen(self.fq) as fi, open(self.fai, "w") as fo:
            while True:
                cur = fi.tell()
                _ = fi.readline()
                if cur == fi.tell():
                    break
                fo.write("%d\n" % cur)
                fi.readline()
                fi.readline()
                fi.readline()

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
        with Zopen(fqi.fq) as fi, open(fqi.fai, "w") as fo:
            while True:
                cur = fi.tell()
                _ = fi.readline()
                if cur == fi.tell():
                    break
                fo.write("%d\n" % cur)
                idx.append(cur)
                fi.readline()
                fi.readline()
                fi.readline()
        return idx


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
        des = "split fastq by barcode."
    parser = argparse.ArgumentParser(
        description=des, prog=" ".join(sys.argv[0:2]))
    general_parser = parser.add_argument_group("General options")
    general_parser.add_argument("mode", metavar=mode, choices=subargs)
    general_parser.add_argument("-i", "--input", type=str, help="input fastq file, required",
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
