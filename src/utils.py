import os
import sys
import json
import gzip
import time
import math
import shutil
import fnmatch
import logging
import argparse
import subprocess

import multiprocessing as mp

from random import uniform
from collections import defaultdict, Counter

from .version import __version__

PY3 = sys.version[0] == "3"

ambiguous_dna_letters = "GATCRYWSMKHBVDN"


class MultiZipHandle(object):

    def __init__(self, mode="rb", **infiles):
        self.info = infiles
        self.handler = {}
        self.mode = mode

    def __enter__(self):
        for sn, f in self.info.items():
            if f.endswith(".gz"):
                self.handler[sn] = gzip.open(f, self.mode)
            else:
                self.handler[sn] = open(f, self.mode)
        return self.handler

    def __exit__(self, exc_type, exc_val, exc_tb):
        for sn, h in self.handler.items():
            h.close()


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


def canonicalize(path):
    return os.path.abspath(os.path.expanduser(path))


def chunk(lst, size=1, group=None):
    if group:
        size = int(math.ceil(len(lst)/float(group)))
    return [lst[i:i + size] for i in range(0, len(lst), size)]


def log(logfile=None, level="info"):
    logger = logging.getLogger()
    if level.lower() == "info":
        logger.setLevel(logging.INFO)
        f = logging.Formatter(
            '[%(levelname)s %(asctime)s] %(message)s')
    elif level.lower() == "debug":
        logger.setLevel(logging.DEBUG)
        f = logging.Formatter(
            '[%(levelname)s %(threadName)s %(asctime)s %(funcName)s(%(lineno)d)] %(message)s')
    if logfile is None:
        h = logging.StreamHandler(sys.stdout)
    else:
        h = logging.FileHandler(logfile, mode='a')
    h.setFormatter(f)
    logger.addHandler(h)
    return logger


def call(cmd, run=True, verbose=False):
    if not run:
        if verbose:
            print(cmd)
        return
    if verbose:
        print(cmd)
        subprocess.check_call(cmd, shell=True, stdout=sys.stdout,
                              stderr=sys.stderr)
    else:
        with open(os.devnull, "w") as fo:
            subprocess.check_call(cmd, shell=True, stdout=fo, stderr=fo)


def run_exe(name):
    exe = os.path.join(os.path.dirname(__file__), name)
    if not is_exe(exe):
        raise IOError("%s not implemented yet" % name)
    subprocess.check_call([exe] + sys.argv[1:])


def clean(*path):
    for p in path:
        if os.path.isfile(p):
            os.remove(p)
        elif os.path.isdir(p):
            shutil.rmtree(p)


def clean_fnmatch(directory, *patterns):
    for root, dirnames, filenames in os.walk(directory, followlinks=False):
        for pat in patterns:
            for filename in fnmatch.filter(filenames, pat):
                filepath = os.path.join(root, filename)
                try:
                    os.remove(filepath)
                except:
                    pass
            for dirname in fnmatch.filter(dirnames, pat):
                dirpath = os.path.join(root, dirname)
                try:
                    shutil.rmtree(dirpath)
                except:
                    pass


if PY3:
    TRANS = str.maketrans("ATGCRYMKSWHBVDN", "TACGYRKMSWDVBHN")
else:
    import string
    TRANS = string.maketrans("ATGCRYMKSWHBVDN", "TACGYRKMSWDVBHN")


def rc_seq(seq):
    return seq.encode().decode().upper().translate(TRANS)[::-1]


def which(program, paths=None):
    ex = os.path.dirname(sys.executable)
    found_path = None
    fpath, fname = os.path.split(program)
    if fpath:
        program = canonicalize(program)
        if is_exe(program):
            found_path = program
    else:
        if is_exe(os.path.join(ex, program)):
            return os.path.join(ex, program)
        paths_to_search = []
        if isinstance(paths, (tuple, list)):
            paths_to_search.extend(paths)
        else:
            env_paths = os.environ.get("PATH", "").split(os.pathsep)
            paths_to_search.extend(env_paths)
        for path in paths_to_search:
            exe_file = os.path.join(canonicalize(path), program)
            if is_exe(exe_file):
                found_path = exe_file
                break
    return found_path


def is_exe(file_path):
    return (
        os.path.exists(file_path)
        and os.access(file_path, os.X_OK)
        and os.path.isfile(os.path.realpath(file_path))
    )


def load_bcl_stats(j):
    with open(j) as fi:
        stat = json.load(fi)
    info = stat["ConversionResults"][0]
    out = []
    for i in info["DemuxResults"]:
        sn = i["SampleName"]
        rdn = i["NumberReads"]
        out.append((sn, rdn))
    out.append(("Unknow", info["Undetermined"]["NumberReads"]))
    return out


def nestdict():
    return defaultdict(nestdict)


def timeRecord(func):
    def wrapper(*args, **kwargs):
        s = time.time()
        value = func(*args, **kwargs)
        sys.stdout.write("Time elapse: %d sec.\n" % int(time.time() - s))
        return value
    return wrapper


def is_ambiguous_dna(seq):
    for i in seq.upper():
        if i not in ambiguous_dna_letters:
            return False
    return True


def parseArg():
    parser = argparse.ArgumentParser(
        description="split a mix fastq or BCL by barcode index.",)
    parser.add_argument("-v", '--version',
                        action='version', version="v" + __version__)
    parent1_parser = argparse.ArgumentParser(add_help=False)
    parent1_parser.add_argument("-i", "--input", type=str, help="input fastq file or BCL flowcell directory, required",
                                required=True, metavar="<str>")
    parent2_parser = argparse.ArgumentParser(add_help=False)
    parent2_grp = parent2_parser.add_argument_group("logging arguments")
    parent2_grp.add_argument('--debug', action='store_true',
                             help='logging debug', default=False)
    parent2_grp.add_argument('--log', type=str,
                             help='append log to file, stdout by default', metavar="<file>")
    subparsers = parser.add_subparsers(
        metavar="command", dest="command")
    parser_index = subparsers.add_parser(
        'index', parents=[parent1_parser, parent2_parser],  help="index fastq file for reading in multi processing, can be instead by `samtools fqidx <fqfile>`.")
    parser_split = subparsers.add_parser(
        'split', parents=[parent2_parser], help="split sequence data by barcode.")
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
        'bcl2fq', parents=[parent1_parser, parent2_parser], help="split flowcell bcl data to fastq.")
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
    parser_bcl2fq.add_argument("--mode", choices=["single", "paired", "auto"], default="auto",
                               help="barcode mode, single-end or paired-end sequence, 'auto' by default")
    parser_bcl2fq.add_argument('--bcl2fq', metavar="<str>",
                               help="bcl2fastq path if necessary, if not set, auto detected")
    return parser.parse_args()
