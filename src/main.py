#!/usr/bin/env python
# coding:utf-8

import os
import multiprocessing as mp


from collections import Counter

from .src import *


def main():
    args = parseArg()
    if args.mode == "index":
        FastqIndex.createindex(args.input)
        return
    infq = args.input
    outdir = args.output
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    group = args.threads * 2

    barcode = {}
    outfile = {"Unknow": os.path.join(outdir, "Unknow.fq.gz")}
    with open(args.barcode) as fi:
        for line in fi:
            if not line.strip() or line.strip().startswith("#"):
                continue
            line = line.split()
            sn, bc = line[0], line[1]
            barcode[bc.encode()] = sn
            outfile[sn] = os.path.join(outdir, sn+".fq.gz")

    idx = FastqIndex.getIndexPos(infq)
    size = int(math.ceil(len(idx)/float(group)))
    l = len(idx)
    pos = []
    for i in range(0, len(idx), size):
        e = min(l, i+size)
        pos.append((idx[i], idx[e-1]))
    outQ = mp.Manager().Queue()
    p = mp.Pool(args.threads)
    for s, e in pos:
        p.apply_async(splitFastq,  args=(infq, int(s), int(e),
                      outQ, barcode, args.mismatch, args.drup))

    sms = Counter()
    ivs = len(pos)
    count = 0
    with MultiZipHandle(**outfile) as fh:
        while count < ivs:
            res = outQ.get()
            if res is None:
                count += 1
            else:
                sn, name, seq, flag, qual = res
                fh[sn].writelines([name, seq, flag, qual])
                sms[sn] += 1
    total_seq = sum(list(sms.values()))
    for sn in sorted(barcode.values()):
        num = sms[sn]
        sys.stdout.write("%s: %d(%.2f%%)\n" %
                         (sn, num, round(num/float(total_seq)*100, 2)))
    sys.stdout.write("Unknow: %d(%.2f%%)\n" % (
        sms["Unknow"], round(sms["Unknow"]/float(total_seq)*100, 2)))


if __name__ == "__main__":
    main()
