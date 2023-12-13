from .utils import *


class BCL(object):

    def __init__(self, fcdir, outdir, bcfile, cpu=60, bcl2fastq="", mis=1, rc_i7=False, rc_i5=False, mode="auto", logfile=None, **kw):
        self.fcdir = fcdir
        self.outdir = outdir
        self.bcfile = bcfile
        self.mis = str(mis)
        self.index = 1
        self.nproc = str(cpu)
        self.rc_i7 = rc_i7
        self.rc_i5 = rc_i5
        self.mode = mode
        self.logfile = logfile
        self.bcl2fastq = bcl2fastq or which("bcl2fastq")
        self.samplesheet = os.path.join(self.outdir, "sample-sheet.csv")

    def create_samplesheet(self):
        idx = []
        ignore_index2 = None
        with open(self.bcfile) as fi:
            for line in fi:
                if not line.strip() or line.startswith("#"):
                    continue
                line = line.split()[:3]
                if len(line) < 2:
                    raise IOError("illegal barcode input file")
                elif len(line) == 2:
                    if not is_ambiguous_dna(line[1]):
                        raise IOError("illegal barcode base %s" % line[1])
                    self.index = 1
                    idx.append((line))
                elif len(line) >= 3:
                    if not is_ambiguous_dna(line[1]) or self.mode == "paired" and not is_ambiguous_dna(line[2]):
                        raise IOError("illegal barcode base %s" % line)
                    elif not is_ambiguous_dna(line[2]) or len(line[2]) < 6:
                        ignore_index2 = 1
                    self.index = 2
                    idx.append((line))
                # else:
                #    raise IOError("illegal barcode file %s" % self.bcfile)
        if self.mode == "single":
            self.index = 1
        elif self.mode == "paired":
            self.index = 2
        else:
            if ignore_index2:
                self.logs.warning(
                    "ignore column 3 in barcode file, using single index")
                self.index = 1
        idx = [i[:self.index+1] for i in idx[:]]
        if len(sum(idx, [])) != len(idx) * (self.index+1):
            raise IOError("illegal barcode input file")
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir)
        with open(self.samplesheet, "w") as fo:
            if self.index == 1:
                fo.write("[Data]\nSample_ID,index\n")
            else:
                fo.write("[Data]\nSample_ID,index,index2\n")
            for line in idx:
                if self.rc_i7:
                    line[1] = rc_seq(line[1])
                if self.index == 2 and self.rc_i5:
                    line[2] = rc_seq(line[2])
                fo.write(",".join(line) + "\n")

    def make_bcl_cmd(self):
        cmd = [self.bcl2fastq, "-o", self.outdir, "-R", self.fcdir,
               "-i", os.path.join(self.fcdir, "Data/Intensities/BaseCalls"),
               "--no-lane-splitting",
               "--barcode-mismatches", self.mis,
               "-p", self.nproc,
               "--sample-sheet", self.samplesheet]
        return cmd

    def call(self, cmd, run=True):
        self.logs.debug(cmd)
        if not run:
            return
        try:
            out = sys.stdout
            if self.logfile:
                out = open(self.logfile, "a")
            subprocess.check_call(cmd, shell=isinstance(
                cmd, str) and True or False, stdout=out if self.logs.level == 10 else -3, stderr=-2)
        finally:
            out.close()

    def stats(self, j):
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

    def run(self):
        self.create_samplesheet()
        cmd = self.make_bcl_cmd()
        self.cmd = cmd
        # self.logs.info(cmd)
        self.logs.info("start bcl2fastq.")
        self.call(cmd, run=True)
        # sms = self.stats(os.path.join(self.outdir, "Stats/Stats.json"))
        # self.logs.info(sms)
        # self.logs.info("Success")

    @property
    def logs(self):
        return logging.getLogger()
