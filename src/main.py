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
    if not os.path.exists(args.input):
        IOError("No such file or directory: %s" % args.input)
    if args.mode == "index":
        logs.info("start fsplit")
        FastqIndex.createindex(args.input)
        return
    infq = args.input
    outdir = os.path.abspath(args.output)
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
    logs.info("start fsplit")
    if os.path.isdir(infq):
        if os.path.isfile(os.path.join(infq, "RTAComplete.txt")):
            with open(os.path.join(outdir, "sample-sheet.csv"), "w") as fo:
                fo.write("[Data]\nSample_ID,Sample_Name,index\n")
                for bcs, bcn in barcode.items():
                    bcs = bcs.decode()
                    fo.write(",".join([bcn, bcn, bcs]) + "\n")
            bcl2fastq = args.bcl2fq or sh.which("bcl2fastq") or os.path.join(
                sys.prefix, "bin", "bcl2fastq")
            if not os.path.isfile(bcl2fastq):
                sys.exit("bcl2fastq not found, exists")
            cmd = [bcl2fastq, "-o", outdir, "-R", infq, "--no-lane-splitting"]
            cmd.extend(["--barcode-mismatches", str(args.mismatch)])
            cmd.extend(["-p", str(args.threads)])
            cmd.extend(
                ["--sample-sheet", os.path.join(outdir, "sample-sheet.csv")])
            cmd.extend(
                ["-i", os.path.join(infq, "Data/Intensities/BaseCalls")])
            logs.info(" ".join(cmd))
            call(" ".join(cmd))
            logs.info("Success")
            sys.stdout.write("\n")
            js = os.path.join(outdir, "Stats/Stats.json")
            if os.path.isfile(js):
                stats = load_bcl_stats(js)
                for sn, rdn in stats:
                    sys.stdout.write("%s: %d\n" % (sn, rdn))
            rmpdir = ["Stats", "Reports", ]
            rmfile = [os.path.join(outdir, i) for i in rmpdir]
            rmcmd = "rm -fr %s" % " ".join(rmfile)
            call(rmcmd)
            return
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
    logs.info("Success")
    sys.stdout.write("\n")
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
