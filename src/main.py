#!/usr/bin/env python
# coding:utf-8

import os
import time
import multiprocessing as mp


from collections import Counter

from .src import *


@timeRecord
def main():
    logs = log()
    args = parseArg()
    infq = args.input
    outdir = os.path.abspath(args.output)
    if not os.path.exists(infq):
        IOError("No such file or directory: %s" % infq)
    logs.info("start fsplit")
    if args.mode == "index":
        FastqIndex.createindex(args.input)
        return
    if os.path.isdir(infq):
        if os.path.isfile(os.path.join(infq, "RTAComplete.txt")):
            bcl2fastq = args.bcl2fq or sh.which("bcl2fastq") or os.path.join(
                sys.prefix, "bin", "bcl2fastq")
            if not os.path.isfile(bcl2fastq):
                sys.exit("bcl2fastq not found, exists")
            bcl = BCL(infq, outdir, args.barcode, args.threads,
                      bcl2fastq=bcl2fastq, mis=args.mismatch)
            bcl.run()
            logs.info("Success")
            js = os.path.join(outdir, "Stats/Stats.json")
            if os.path.isfile(js):
                stats = load_bcl_stats(js)
                logs.info(stats)
            rmpdir = ["Stats", "Reports", "sample-sheet.csv"]
            rmfile = [os.path.join(outdir, i) for i in rmpdir]
            rmcmd = "rm -fr %s" % " ".join(rmfile)
            call(rmcmd)
            return
        else:
            sys.exit("flowcell directory '%s' incorrection" % infq)
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
    idx = FastqIndex.getIndexPos(infq)
    size = int(math.ceil(len(idx)/float(group)))
    l = len(idx)
    pos = []
    for i in range(0, len(idx), size):
        e = min(l-1, i+size)
        pos.append((idx[i], idx[e]))
    outQ = mp.Manager().Queue()
    p = mp.Pool(args.threads)
    for s, e in pos:
        p.apply_async(splitFastq,  args=(infq, int(s), int(e),
                      outQ, barcode, args.mismatch, args.drup, outdir))

    sms = Counter()
    ivs = len(pos)
    count = 0
    total_seq = 0
    for _, f in outfile.items():
        if os.path.isfile(f):
            call("rm -fr " + f)
    for s, e in pos:
        res = outQ.get()
        total, snms, snoutf = res
        sms.update(snms)
        total_seq += total
        for sn, _of in snoutf.items():
            call("cat %s >> %s" % (_of, outfile[sn]))
            call("rm -fr " + _of)
    logs.info("Success")
    sys.stdout.write("\n")
    sms["Unknow"] = total_seq - sum(sms.values())
    for sn in sorted(barcode.values()):
        num = sms[sn]
        sys.stdout.write("%s: %d(%.2f%%)\n" %
                         (sn, num, round(num/float(total_seq)*100, 2)))
    sys.stdout.write("Unknow: %d(%.2f%%)\n" % (
        sms["Unknow"], round(sms["Unknow"]/float(total_seq)*100, 2)))


if __name__ == "__main__":
    main()
