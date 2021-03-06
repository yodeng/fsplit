#!/usr/bin/env python
# coding:utf-8

import os
import multiprocessing as mp


from collections import Counter

from .src import *


@timeRecord
def main_idx_multi():
    logs = log()
    args = parseArg()
    infq = args.input
    if not os.path.exists(infq):
        raise IOError("No such file or directory: %s" % infq)
    logs.info("start fsplit")
    if args.mode == "index":
        FastqIndex.createindex(args.input)
        return
    outdir = os.path.abspath(args.output)
    if os.path.isdir(infq):
        if os.path.isfile(os.path.join(infq, "RTAComplete.txt")):
            bcl2fastq = args.bcl2fq or which("bcl2fastq") or os.path.join(
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
            clean(*rmfile)
            return
        else:
            sys.exit("flowcell directory '%s' incorrection" % infq)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    group = args.threads * 5

    mg = mp.Manager()
    barcode = {}
    # outfile = {"Unknow": os.path.join(outdir, "Unknow.fq")}
    filelock = {}
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
                outfile[sn] += ".gz"
            if os.path.isfile(outfile[sn]):
                os.remove(outfile[sn])
            filelock[sn] = mg.Lock()
    idx = FastqIndex.getIndexPos(infq)
    size = int(math.ceil(len(idx)/float(group)))
    l = len(idx)
    pos = []
    for i in range(0, len(idx), size):
        e = min(l-1, i+size)
        pos.append((idx[i], idx[e]))
    outQ = mg.Queue()
    p = mp.Pool(args.threads)
    for s, e in pos:
        p_args = tuple([infq, int(s), int(e), outQ, barcode,
                       args.mismatch, args.drup, outfile, filelock])
        p.apply_async(splitFastq,  args=p_args)

    sms = Counter()
    ivs = len(pos)
    count = 0
    total_seq = 0
    for s, e in pos:
        res = outQ.get()
        total, snms = res
        sms.update(snms)
        total_seq += total
    logs.info("Success")
    sys.stdout.write("\n")
    sms["Unknow"] = total_seq - sum(sms.values())
    for sn in sorted(barcode.values()):
        num = sms[sn]
        sys.stdout.write("%s: %d(%.2f%%)\n" %
                         (sn, num, round(num/float(total_seq)*100, 2)))
    sys.stdout.write("Unknow: %d(%.2f%%)\n" % (
        sms["Unknow"], round(sms["Unknow"]/float(total_seq)*100, 2)))


@timeRecord
def main():
    logs = log()
    args = parseArg()
    infq = args.input
    if not os.path.exists(infq):
        raise IOError("No such file or directory: %s" % infq)
    logs.info("start fsplit")
    if args.mode == "index":
        FastqIndex.createindex(args.input)
        return
    outdir = os.path.abspath(args.output)
    if os.path.isdir(infq):
        if os.path.isfile(os.path.join(infq, "RTAComplete.txt")):
            bcl2fastq = args.bcl2fq or which("bcl2fastq") or os.path.join(
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
            clean(*rmfile)
            return
        else:
            sys.exit("flowcell directory '%s' incorrection" % infq)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

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
                outfile[sn] += ".gz"
            if os.path.isfile(outfile[sn]):
                os.remove(outfile[sn])

    drup_pos = dict.fromkeys(barcode.keys(), 0)
    if args.drup:
        for bc in barcode:
            drup_pos[bc] = len(bc)

    mis = args.mismatch
    sms = Counter()
    total_seq = 0
    with MultiZipHandle(mode="wb", **outfile) as fh:
        with Zopen(infq, gzip=True) as fi:
            seq = [None] * 4
            for n, line in enumerate(fi):
                ni = n % 4
                seq[ni] = line
                if ni == 3:
                    for b, sn in barcode.items():
                        d = 0
                        for n, i in enumerate(b):
                            if i != seq[1][n]:
                                d += 1
                            if d > mis:
                                break
                        else:
                            dp = drup_pos[b]
                            seq[1] = seq[1][dp:]
                            seq[3] = seq[3][dp:]
                            fh[sn].writelines(seq)
                            sms[sn] += 1
                            break
                    total_seq += 1
    logs.info("Success")
    sys.stdout.write("\n")
    sms["Unknow"] = total_seq - sum(sms.values())
    for sn in sorted(set(barcode.values())):
        num = sms[sn]
        sys.stdout.write("%s: %d(%.2f%%)\n" %
                         (sn, num, round(num/float(total_seq)*100, 2)))
    sys.stdout.write("Unknow: %d(%.2f%%)\n" % (
        sms["Unknow"], round(sms["Unknow"]/float(total_seq)*100, 2)))


if __name__ == "__main__":
    main()
