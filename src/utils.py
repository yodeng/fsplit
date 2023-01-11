import os
import sys
import json
import gzip
import time
import math
import shutil
import logging
import subprocess

from collections import defaultdict

PY3 = sys.version[0] == "3"


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
        h = logging.FileHandler(logfile, mode='w')
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


def clean(*path):
    for p in path:
        if os.path.isfile(p):
            os.remove(p)
        elif os.path.isdir(p):
            shutil.rmtree(p)


if PY3:
    TRANS = str.maketrans("ATGCRYMKSWHBVDN", "TACGYRKMSWDVBHN")
else:
    import string
    TRANS = string.maketrans("ATGCRYMKSWHBVDN", "TACGYRKMSWDVBHN")


def rc_seq(seq):
    return seq.upper().translate(TRANS)[::-1]


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
