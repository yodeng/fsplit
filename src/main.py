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
        if os.path.isdir(os.path.join(infq, "Data/Intensities/BaseCalls")):
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
    if args.mode == "bcl2fq":
        if os.path.isdir(os.path.join(infq, "Data/Intensities/BaseCalls")):
            bcl2fastq = args.bcl2fq or which("bcl2fastq")
            if not (bcl2fastq and os.path.isfile(bcl2fastq)):
                sys.exit("bcl2fastq not found, exit")
            kw = {
                "cpu": args.threads,
                "bcl2fastq": bcl2fastq,
                "mis": args.mismatch,
                "bcl2fastq": bcl2fastq,
                "rc_i7": args.rc_index1,
                "rc_i5": args.rc_index2,
                "print_cmd": False,
            }
            bcl = BCL(infq, outdir, args.sample, **kw)
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
    barcode_paired = {}
    # outfile = {"Unknow": os.path.join(outdir, "Unknow.fq")}
    outfile = {}
    outfile_paired = {}
    with open(args.barcode) as fi:
        for line in fi:
            if not line.strip() or line.strip().startswith("#"):
                continue
            line = line.split()
            sn, bc = line[0], line[1]
            bc1 = args.rc_bc1 and rc_seq(bc).encode() or bc.encode()
            barcode[bc1] = sn
            bc2 = len(line) == 3 and line[2] or bc
            barcode_paired[bc1] = args.rc_bc2 and rc_seq(
                bc2).encode() or bc2.encode()
            if args.Input:
                outfile[sn] = os.path.join(outdir, sn+".R1.fq")
                outfile_paired[sn] = os.path.join(outdir, sn+".R2.fq")
            else:
                outfile[sn] = os.path.join(outdir, sn+".fq")
            if args.output_gzip:
                outfile[sn] += ".gz"
                if sn in outfile_paired:
                    outfile_paired[sn] += ".gz"

    drup_pos = dict.fromkeys(list(barcode.keys()) +
                             list(barcode_paired.values()), 0)
    if args.drup:
        for bc in barcode.keys():
            drup_pos[bc] = len(bc)
        for bcp in barcode_paired.values():
            drup_pos[bcp] = len(bcp)

    mis = args.mismatch
    sms = Counter()
    total_seq = 0
    if not args.Input:
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
    else:
        with MultiZipHandle(mode="wb", **outfile) as f1:
            with MultiZipHandle(mode="wb", **outfile_paired) as f2:
                seq1 = [None] * 4
                seq2 = [None] * 4
                with Zopen(infq, gzip=True) as fq1, Zopen(args.Input, gzip=True) as fq2:
                    for n, line in enumerate(fq1):
                        ni = n % 4
                        seq1[ni] = line
                        seq2[ni] = fq2.readline()
                        if ni == 3:
                            for b, sn in barcode.items():
                                d = 0
                                for n, i in enumerate(b):
                                    if i != seq1[1][n]:
                                        d += 1
                                    if d > mis:
                                        break
                                else:
                                    d2 = 0
                                    b2 = barcode_paired[b]
                                    for n, i in enumerate(b2):
                                        if i != seq2[1][n]:
                                            d2 += 1
                                        if d2 > mis:
                                            break
                                    else:
                                        seq1[1] = seq1[1][drup_pos[b]:]
                                        seq1[3] = seq1[3][drup_pos[b]:]
                                        seq2[1] = seq2[1][drup_pos[b2]:]
                                        seq2[3] = seq2[3][drup_pos[b2]:]
                                        f1[sn].writelines(seq1)
                                        f2[sn].writelines(seq2)
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


def gsplit():
    gsplit_exe = os.path.join(os.path.dirname(__file__), "gsplit")
    if not os.path.isfile(gsplit_exe) or not os.access(gsplit_exe, os.R_OK | os.X_OK):
        raise IOError("gsplit not implemented yet")
    subprocess.check_call([gsplit_exe] + sys.argv[1:])


if __name__ == "__main__":
    main()
