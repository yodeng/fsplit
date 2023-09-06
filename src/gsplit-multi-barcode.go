package main

import (
	"fmt"
	"io"
	"os"
	"regexp"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/yodeng/go-kit/hflag"
	"github.com/yodeng/xopen"
)

const VERSION = "v2023.09.05 15:00"

type SplitFlags struct {
	Fqfile      []string `hflag:"--input, -i; required; usage: input fastq file, *.gz/xz/zst or uncompress allowed, multi-input can be separated by ',' or whitespace, required"`
	Barcodefile string   `hflag:"--barcode, -b; required; usage: barcode and sample file, 'samplename barcode1 barcode1_pos barcode2 barcode2_pos ...', required"`
	Output      string   `hflag:"--output, -o; required; usage: output directory, will create if not exists, required"`
	Threads     int      `hflag:"--threads, -t; default: 10; usage: threads core, 10 by default"`
	Mismatch    int      `hflag:"--mismatch, -m; default: 0; usage: mismatch allowed for barcode search, 0 by default"`
	Split       int      `hflag:"--splitn, -s; default: 1; usage: split barcode output into N files, 1 by default"`
	Drup        bool     `hflag:"--drup, -d; default: false; usage: drup barcode sequence in output if set"`
	Nogz        bool     `hflag:"--no-gzip, -n; default: false; usage: do not gzip output fastq file"`
	Version     bool     `hflag:"--version, -v; usage: show version and exit"`
}

func checkError(err error) {
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func splitSample(seq string, bc map[string]string, mis int) (name, barcode string) {
BCLOOP:
	for b, sn := range bc {
		d := 0
		for n := range b {
			if seq[n] != b[n] {
				d++
			}
			if d > mis {
				continue BCLOOP
			}
		}
		return sn, b
	}
	return
}

func splitBarcode(seq string, bc []string, mis int) (barcode string) {
BCLOOP:
	for _, b := range bc {
		d := 0
		for n := range b {
			if seq[n] != b[n] {
				d++
			}
			if d > mis {
				continue BCLOOP
			}
		}
		barcode = b
		break BCLOOP
	}
	return
}

func min(x, y int) int {
	if x < y {
		return x
	}
	return y
}

func write_seq(seq *[4]string, s int, e int, out *xopen.Writer) {
	if e <= s && e > 0 {
		return
	}
	if e >= len(seq[1]) || e == 0 {
		out.WriteString(seq[0])
		out.WriteString(seq[1][s:])
		out.WriteString(seq[2])
		out.WriteString(seq[3][s:])
	} else {
		out.WriteString(seq[0])
		out.WriteString(seq[1][s:e] + "\n")
		out.WriteString(seq[2])
		out.WriteString(seq[3][s:e] + "\n")
	}
}

func clean_empty_file(files []string) {
	for _, path := range files {
		if isExist(path) {
			_, err := xopen.Ropen(path)
			if err == xopen.ErrNoContent {
				os.Remove(path)
			}
		}
	}
}

func ReadSeq(ch1 chan []string, infile string) {
	r, err := xopen.Ropen(infile)
	checkError(err)
	defer r.Close()
	seq := make([]string, 5, 5)
	for n := 0; ; n++ {
		line, err := r.ReadString('\n')
		if err == io.EOF {
			break
		}
		i := n % 4
		seq[i] = line
		if i == 3 {
			ch1 <- seq
			seq = make([]string, 5, 5)
		}
	}
	close(ch1)
}

func SplitSeq(ch1, ch2 chan []string, barcode map[string]string, mis int, drup bool, wg *sync.WaitGroup) {
	defer wg.Done()
	drup_pos := make(map[string]int, 10)
	for k, _ := range barcode {
		if drup {
			drup_pos[k] = len(k)
		} else {
			drup_pos[k] = 0
		}
	}
	for v := range ch1 {
	bcloop:
		for b, sn := range barcode {
			d := 0
			for n, _ := range b {
				if v[1][n] != b[n] {
					d++
				}
				if d > mis {
					continue bcloop
				}
			}
			dp := drup_pos[b]
			v[1] = v[1][dp:]
			v[3] = v[3][dp:]
			v[4] = sn
			ch2 <- v
			break bcloop
		}
	}
}

func isExist(path string) bool {
	_, err := os.Stat(path)
	if err != nil {
		if os.IsExist(err) {
			return true
		}
		if os.IsNotExist(err) {
			return false
		}
		return false
	}
	return true
}

func strslice2map(sl []string) map[string]struct{} {
	set := make(map[string]struct{}, len(sl))
	for _, v := range sl {
		set[v] = struct{}{}
	}
	return set
}

func inSlice(sl []string, s string) bool {
	m := strslice2map(sl)
	_, ok := m[s]
	return ok
}

func split_fq(args *SplitFlags) (samples []string, sms map[string]map[string]int, total int, outfiles []string) {
	barcode := make(map[string]map[string][2]int)
	drup_pos := make(map[string]map[string][2]int)

	bcf, err := xopen.Ropen(args.Barcodefile)
	checkError(err)
	defer bcf.Close()

	fout := make(map[string]map[string][]*xopen.Writer)

	if !isExist(args.Output) {
		err := os.MkdirAll(args.Output, os.ModePerm)
		checkError(err)
	}
	samples = []string{}
	sms = make(map[string]map[string]int)
	outfiles = []string{}
	for {
		line, err := bcf.ReadString('\n')
		if err == io.EOF {
			break
		}
		line = strings.TrimSpace(line)
		if len(line) == 0 || strings.HasPrefix(line, "#") {
			continue
		}
		reg := regexp.MustCompile(`\s+`)
		line_s := reg.Split(line, -1)
		sn := line_s[0]
		bc_info := make(map[string][2]int)
		b_ := make([]string, 0)
		p_ := make([]int, 0)
		for n, v := range line_s[1:] {
			if n%2 > 0 {
				p, _ := strconv.Atoi(v)
				p_ = append(p_, p-1)
			} else {
				b_ = append(b_, v)
			}
		}
		p_ = append(p_, 0)
		if len(p_) == 0 {
			fmt.Printf("read %s error", args.Barcodefile)
			os.Exit(1)
		}
		if p_[0] > 0 {
			bc_info[""] = [2]int{0, p_[0]}
		}
		for n, b := range b_ {
			bc_info[b] = [2]int{p_[n], p_[n+1]}
		}
		dp := make(map[string][2]int)
		for b, pos := range bc_info {
			if args.Drup {
				dp[b] = [2]int{pos[0] + len(b), pos[1]}
			} else {
				dp[b] = pos
			}
		}
		barcode[sn] = bc_info
		drup_pos[sn] = dp
		if !inSlice(samples, sn) {
			samples = append(samples, sn)
			sms[sn] = make(map[string]int)
			fo := make(map[string][]*xopen.Writer, 10)
			bn := 0
			for bc, _ := range bc_info {
				for n := 0; n < args.Split; n++ {
					outf := args.Output + "/" + sn + fmt.Sprintf("_b%d_%05d", bn+1, n+1) + ".fq"
					if !args.Nogz {
						outf += ".gz"
					}
					fh, _ := xopen.Wopen(outf)
					defer fh.Close()
					outfiles = append(outfiles, outf)
					fo[bc] = append(fo[bc], fh)
				}
				bn += 1
			}
			fout[sn] = fo
		}
	}

	mis := args.Mismatch
	seq := [4]string{}

	for _, fqfile := range args.Fqfile {
		r, err := xopen.Ropen(fqfile)
		defer r.Close()
		checkError(err)
		linefo := 0
		for {
			line, err := r.ReadString('\n')
			if err == io.EOF {
				break
			}
			i := linefo % 4
			seq[i] = line
			if i == 3 {
			BCLOOP:
				for sn, bc := range barcode {
					for b, pos := range bc {
						if b == "" {
							continue
						}
						d := 0
						for n := range b {
							if seq[1][pos[0]+n] != b[n] {
								d++
							}
							if d > mis {
								continue BCLOOP
							}
						}
					}
					for b, _ := range bc {
						s := drup_pos[sn][b][0]
						e := drup_pos[sn][b][1]
						bc_chunk_num := sms[sn][b] % args.Split
						write_seq(&seq, s, e, fout[sn][b][bc_chunk_num])
						sms[sn][b] += 1
					}
					break BCLOOP
				}
				total += 1
			}
			linefo++
		}
	}
	return
}

func main() {
	args := &SplitFlags{}
	if err := hflag.Bind(args); err != nil {
		panic(err)
	}
	err := hflag.AddDesc(fmt.Sprintf("for split fastq data from a mixed fastq by barcode/index of each sample. (version: %v)\n", VERSION))
	err = hflag.Parse()
	if len(os.Args) == 1 {
		fmt.Println(hflag.Usage())
		return
	}
	if args.Version {
		fmt.Println(VERSION)
		return
	}
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
	runtime.GOMAXPROCS(args.Threads)
	t := time.Now()
	samples, sms, total, outfiles := split_fq(args)
	clean_empty_file(outfiles)
	fmt.Println()
	snm := 0
	for _, sn := range samples {
		count := 0
		for _, c := range sms[sn] {
			count = c
		}
		snm += count
		fmt.Printf("%v: %v(%.2f%%)\n", sn, count, float64(count)/float64(total)*100)
	}
	fmt.Printf("Unknow: %v(%.2f%%)\n", total-snm, float64(total-snm)/float64(total)*100)
	d := time.Since(t)
	fmt.Printf("\nTime elapse: %v sec.\n", d)
}
