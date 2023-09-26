#!/usr/bin/env python
# coding:utf-8

from .utils import *


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
