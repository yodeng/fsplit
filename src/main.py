#!/usr/bin/env python
# coding:utf-8

import os
import time
import multiprocessing as mp


from collections import Counter

from .src import *


@timeRecord
def main():
    args = parseArg()
    if not os.path.isfile(args.input):
        sys.exit("No such file %s" % args.input)
    if args.mode == "index":
        FastqIndex.createindex(args.input)
        return
    infq = args.input
    outdir = args.output
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    group = args.threads * 5

    barcode = {}
    # outfile = {"Unknow": os.path.join(outdir, "Unknow.fq")}
    outfile = {}
    with open(args.barcode) as fi:
        for line in fi:
            if not line.strip() or line.strip().startswith("#"):
                continue
            line = line.split()
            sn, bc = line[0], line[1]
            barcode[bc.encode()] = sn
            outfile[sn] = os.path.join(outdir, sn+".fq")
    if args.output_gzip:
        for sn in outfile:
            outfile[sn] += ".gz"
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
            res_p = outQ.get()
            for sn, seqs in res_p.items():
                h = fh[sn]
                sms[sn] += len(seqs) / 4
                for seq in seqs:
                    h.write(seq)
            count += 1
    total_seq = l
    sms["Unknow"] = total_seq - sum(sms.values())
    for sn in sorted(barcode.values()):
        num = sms[sn]
        sys.stdout.write("%s: %d(%.2f%%)\n" %
                         (sn, num, round(num/float(total_seq)*100, 2)))
    sys.stdout.write("Unknow: %d(%.2f%%)\n" % (
        sms["Unknow"], round(sms["Unknow"]/float(total_seq)*100, 2)))


if __name__ == "__main__":
    main()
