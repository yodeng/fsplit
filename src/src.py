#!/usr/bin/env python
# coding:utf-8

import os
import sys
import math
import gzip
import argparse
import editdistance


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
            self.handler[sn] = gzip.open(f, "wb")

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
            h = fi.readline().split()
            if len(h) == 6:
                e = int(h[5])
            idx.append(0)
            for line in fi:
                line = line.split()
                if len(line) == 6:
                    idx.append(e+int(line[4]))
                    e = int(line[5])
                else:
                    idx.append(int(line[0]))
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
    return parser.parse_args()


def splitFastq(fq, s, e, outQ, barcode, mis, drup=False):
    drup_pos = dict.fromkeys(barcode.keys(), 0)
    if drup:
        for bc in barcode:
            drup_pos[bc] = len(bc)
    with Zopen(fq) as fi:
        fi.seek(s)
        while True:
            if fi.tell() > e:
                break
            name = fi.readline()
            seq = fi.readline()
            flag = fi.readline()
            qual = fi.readline()
            for b in barcode:
                bc_len = len(b)
                bar = seq[:bc_len]
                dis = editdistance.eval(b, bar)
                if dis <= mis:
                    dp = drup_pos[b]
                    sn = barcode[b]
                    seq = seq[dp:]
                    qual = qual[dp:]
                    outQ.put((sn, name, seq, flag, qual))
                    break
            else:
                outQ.put(("Unknow", name, seq, flag, qual))
        outQ.put(None)
