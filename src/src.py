#!/usr/bin/env python
# coding:utf-8

import os
import sys
import gzip
import logging
import argparse
import subprocess

from random import uniform

from .version import __version__
from .bcl import *
from .utils import *


class Zopen(object):
    def __init__(self, name,  mode="rb", gzip=False):
        self.name = name
        self.mode = mode
        self.gzip = gzip
        self.handler = None

    def __enter__(self):
        if self.name.endswith(".gz"):
            if self.gzip and "r" in self.mode:
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
                    "read fastq index (%s) done", self.fai)
                return idx
            h = h.split()
            if len(h) == 6:
                idx.append(0)
                e = int(h[5])
                for line in fi:
                    line = line.split()
                    idx.append(e+int(line[4]))
                    e = int(line[5])
                else:
                    idx.append(e+int(line[4]))
        self.logs.info(
            "read fastq index (%s) done.", self.fai)
        return idx

    def index_crt(self):
        with Zopen(self.fq, gzip=True) as fi, open(self.fai, "w") as fo:
            cur = 0
            fo.write("%d\n" % cur)
            for n, line in enumerate(fi):
                if n % 4 == 0:
                    if uniform(0, 1) < 0.001:
                        fo.write("%d\n" % cur)
                cur += len(line)
            fo.write("%d\n" % cur)
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
            fo.write("%d\n" % cur)
            idx.append(cur)
            for n, line in enumerate(fi):
                if n % 4 == 0:
                    if uniform(0, 1) < 0.001:
                        fo.write("%d\n" % cur)
                        idx.append(cur)
                cur += len(line)
            fo.write("%d\n" % cur)
            idx.append(cur)
            fqi.logs.info("create fastq index done.")
        return idx

    @property
    def logs(self):
        return logging.getLogger()


def parseArg():
    parser = argparse.ArgumentParser(
        description="split a mix fastq or BCL by barcode index.",)
    parser.add_argument("-v", '--version',
                        action='version', version="v" + __version__)
    parent_parser = argparse.ArgumentParser(add_help=False)
    general_parser = parent_parser.add_argument_group("common options")
    general_parser.add_argument("-i", "--input", type=str, help="input fastq file or BCL flowcell directory, required",
                                required=True, metavar="<str>")
    subparsers = parser.add_subparsers(
        metavar="command", dest="mode")
    parser_index = subparsers.add_parser(
        'index', parents=[parent_parser],  help="index fastq file for reading in multi processing, can be instead by `samtools fqidx <fqfile>`.")
    parser_split = subparsers.add_parser(
        'split', help="split sequence data by barcode.")
    parser_split.add_argument("-i", "--input", type=str, help="input fastq file, required",
                              required=True, metavar="<file>")
    parser_split.add_argument("-I", "--Input", type=str, help="input paired fastq file",
                              required=False, metavar="<file>")
    parser_split.add_argument("-b", "--barcode", type=str,
                              help='sample and barcode sequence info, two or three columns like "sampleName barcodeSeq1 barcodeSeq2", required', required=True, metavar="<file>")
    parser_split.add_argument('-m', "--mismatch", help="mismatch allowed for barcode search, 0 by default",
                              type=int, default=0, metavar="<int>")
    parser_split.add_argument('-o', "--output", help="output directory, required",
                              type=str, required=True, metavar="<str>")
    parser_split.add_argument("-d", '--drup',   action='store_true',
                              help="drup barcode sequence in output if set",  default=False)
    parser_split.add_argument("-rc1", "--rc-bc1", action="store_true", default=False,
                              help='reverse complement barcode1')
    parser_split.add_argument("-rc2", "--rc-bc2", action="store_true", default=False,
                              help='reverse complement barcode2')
    parser_split.add_argument("--output-gzip",   action='store_true',
                              help="gzip output fastq file, this will make your process slower", default=False)
    parser_bcl2fq = subparsers.add_parser(
        'bcl2fq', parents=[parent_parser], help="split flowcell bcl data to fastq.")
    parser_bcl2fq.add_argument('-t', "--threads", help="threads core, 10 by default",
                               type=int, default=10, metavar="<int>")
    parser_bcl2fq.add_argument("-s", "--sample", type=str,
                               help='sample index file, two or three columns like "sample index1(i7) index2(i5)", required', required=True, metavar="<file>")
    parser_bcl2fq.add_argument('-m', "--mismatch", help="mismatch allowed for barcode search, 1 by default",
                               type=int, default=1, metavar="<int>")
    parser_bcl2fq.add_argument('-o', "--output", help="output directory, required",
                               type=str, required=True, metavar="<str>")
    parser_bcl2fq.add_argument("-rc1", "--rc-index1", action="store_true", default=False,
                               help='reverse complement index1(i7)')
    parser_bcl2fq.add_argument("-rc2", "--rc-index2", action="store_true", default=False,
                               help='reverse complement index2(i5)')
    parser_bcl2fq.add_argument('--bcl2fq', metavar="<str>",
                               help="bcl2fastq path if necessary, if not set, auto detected")
    return parser.parse_args()


def splitFastq(fq, s, e, outQ, barcode, mis=0, drup=False, outfile=None, filelock=None):
    drup_pos = dict.fromkeys(barcode.keys(), 0)
    total = 0
    if drup:
        for bc in barcode:
            drup_pos[bc] = len(bc)
    out = {}
    with Zopen(fq) as fi:
        fi.seek(s)
        while True:
            if fi.tell() == e:
                break
            name = fi.readline()
            seq = fi.readline()
            flag = fi.readline()
            qual = fi.readline()
            for b in barcode:
                d = 0
                for n, i in enumerate(b):
                    if i != seq[n]:
                        d += 1
                    if d > mis:
                        break
                else:
                    dp = drup_pos[b]
                    sn = barcode[b]
                    seq = seq[dp:]
                    qual = qual[dp:]
                    out.setdefault(sn, []).extend((name, seq, flag, qual))
                    break
            # else:
            #    out.setdefault("Unknow", []).extend((name, seq, flag, qual))
            total += 1
    snms = {}
    for sn, seqs in out.items():
        with filelock[sn]:
            with Zopen(outfile[sn], "ab") as fo:
                snms[sn] = len(seqs)/4
                for line in seqs:
                    fo.write(line)
    outQ.put((total, snms))
