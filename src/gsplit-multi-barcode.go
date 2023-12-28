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

const VERSION = "v2023.12.28 15:00"

type SplitFlags struct {
	Fqfile      []string `hflag:"--input, -i; required; usage: input fastq file, *.gz/xz/zst or uncompress allowed, multi-input can be separated by ',' or whitespace, required"`
	Barcodefile string   `hflag:"--barcode, -b; required; usage: barcode and sample file, 'samplename barcode1 barcode1_pos barcode2 barcode2_pos ...', required"`
	Pos         string   `hflag:"--pos, -p; usage: barcode parts pos info"`
	Output      string   `hflag:"--output, -o; required; usage: output directory, will create if not exists, required"`
	Threads     int      `hflag:"--threads, -t; default: 10; usage: threads core, 10 by default"`
	Mismatch    int      `hflag:"--mismatch, -m; default: 0; usage: mismatch allowed for barcode search, 0 by default"`
	Split       int      `hflag:"--splitn, -s; default: 1; usage: split barcode output into N files, 1 by default"`
	Trimleft    int      `hflag:"--left-trim, -l; default: 0; usage: trim N base in the first barcode output from left, 0 by default"`
	Trimright   int      `hflag:"--right-trim, -r; default: 0; usage: trim N base in the last barcode output from right, 0 by default"`
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

func write_seq(seq *[4]string, s int, e int, b string, p int, out *xopen.Writer) {
	if e <= s && e > 0 {
		return
	}
	if length := len(seq[1]); e == 0 || e >= length {
		name := fmt.Sprintf("%v barcode:%v,part:%d,pos:%d-%d\n", strings.TrimSpace(seq[0]), b, p+1, s+1, len(seq[1])+1)
		out.WriteString(name)
		out.WriteString(seq[1][s:])
		out.WriteString(seq[2])
		out.WriteString(seq[3][s:])
	} else {
		name := fmt.Sprintf("%v barcode:%v,part:%d,pos:%d-%d\n", strings.TrimSpace(seq[0]), b, p+1, s+1, e+1)
		out.WriteString(name)
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

func has_key(k string, dict map[string][][][2]int) bool {
	_, ok := dict[k]
	return ok
}

func split_fq(args *SplitFlags) (samples []string, sms map[string]map[string]map[int]int, total int, outfiles []string) {
	barcode := make(map[string]map[string][2]int)
	drup_pos := make(map[string]map[string][2]int)

	barcode_split := make(map[string]map[string][][2]int)

	bcf, err := xopen.Ropen(args.Barcodefile)
	checkError(err)
	defer bcf.Close()

	fout := make(map[string]map[string]map[int][]*xopen.Writer)

	if !isExist(args.Output) {
		err := os.MkdirAll(args.Output, os.ModePerm)
		checkError(err)
	}
	samples = []string{}
	sms = make(map[string]map[string]map[int]int)
	outfiles = []string{}
	sn2bc := make(map[string][]string)

	sn_bc_pos_info := make(map[string][][][2]int, 0)
	if len(args.Pos) > 0 {
		posfile, err := xopen.Ropen(args.Pos)
		checkError(err)
		defer posfile.Close()
		for {
			line, err := posfile.ReadString('\n')
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
			bc_pos_info := make([][][2]int, 0)
			for _, pos_info := range line_s[1:] {
				pos_pairs := make([][2]int, 0)
				reg_p := regexp.MustCompile(`,`)
				pos_pair_list := reg_p.Split(pos_info, -1)
				for _, pos_pair_info := range pos_pair_list {
					reg_p_p := regexp.MustCompile(`-`)
					pos_pair_s_e := reg_p_p.Split(pos_pair_info, -1)
					s, _ := strconv.Atoi(pos_pair_s_e[0])
					e, _ := strconv.Atoi(pos_pair_s_e[1])
					s = s - 1
					e = e - 1
					if e <= 0 {
						e = -1
					}
					if s <= 0 {
						s = 0
					}
					if s > e && e != -1 {
						fmt.Println("start pos great then end pos", s, e)
					}
					pos_pair := [2]int{s, e}
					pos_pairs = append(pos_pairs, pos_pair)
				}
				bc_pos_info = append(bc_pos_info, pos_pairs)
				sn_bc_pos_info[sn] = bc_pos_info
			}
		}
	}

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
		bc_list := make([]string, 0)
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
			bc_list = append(bc_list, b)
		}
		sn2bc[sn] = bc_list
		dp := make(map[string][2]int)
		for b, pos := range bc_info {
			if args.Drup {
				dp[b] = [2]int{pos[0] + len(b), pos[1]}
			} else {
				dp[b] = pos
			}
			if pos[0] == 0 {
				dp[b] = [2]int{dp[b][0] + args.Trimleft, dp[b][1]}
			}
			if pos[1] == 0 {
				dp[b] = [2]int{dp[b][0], -args.Trimright - 1}
			}
		}
		barcode[sn] = bc_info
		drup_pos[sn] = dp
		bc_parts_pos := make(map[string][][2]int, 0)
		for n, bc := range bc_list {
			if has_key(sn, sn_bc_pos_info) && n <= len(sn_bc_pos_info[sn])-1 {
				pos_info := sn_bc_pos_info[sn][n]
				for _, pos_p := range pos_info {
					if dp[bc][1] != -1 && (pos_p[0] < dp[bc][0] || pos_p[1] > dp[bc][1]) {
						err := fmt.Sprintf("sample:'%v', barcode:'%v', position %v out of barcode pos %v", sn, bc, pos_p, dp[bc])
						panic(err)
					}
				}
				bc_parts_pos[bc] = pos_info
			} else {
				pos_info := make([][2]int, 0)
				pos_info = append(pos_info, dp[bc])
				bc_parts_pos[bc] = pos_info
			}
		}
		barcode_split[sn] = bc_parts_pos

		if !inSlice(samples, sn) {
			samples = append(samples, sn)
			sms[sn] = make(map[string]map[int]int)
			fo := make(map[string]map[int][]*xopen.Writer, 10)
			bn := 0
			for bc, _ := range bc_info {
				sms[sn][bc] = make(map[int]int)
				fop := make(map[int][]*xopen.Writer, 10)
				for p, _ := range bc_parts_pos[bc] {
					for n := 0; n < args.Split; n++ {
						outf := args.Output + "/" + sn + fmt.Sprintf("_b%d_part%d_%05d", bn+1, p+1, n+1) + ".fq"
						if !args.Nogz {
							outf += ".gz"
						}
						fh, _ := xopen.Wopen(outf)
						defer fh.Close()
						outfiles = append(outfiles, outf)
						fop[p] = append(fop[p], fh)
					}
				}
				fo[bc] = fop
				bn += 1
			}
			fout[sn] = fo
		}
	}

	/*
	   fmt.Println("###", barcode_split)
	   fmt.Println("~~~", fout )
	*/

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
						for p, se := range barcode_split[sn][b] {
							s := se[0]
							e := se[1]
							if e < 0 {
								e += len(seq[1])
							}
							bc_chunk_num := sms[sn][b][p] % args.Split
							write_seq(&seq, s, e, b, p, fout[sn][b][p][bc_chunk_num])
							sms[sn][b][p] += 1
						}

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
		for _, b := range sms[sn] {
			for _, c := range b {
				count = c
			}
		}
		snm += count
		fmt.Printf("%v: %v(%.2f%%)\n", sn, count, float64(count)/float64(total)*100)
	}
	fmt.Printf("Unknow: %v(%.2f%%)\n", total-snm, float64(total-snm)/float64(total)*100)
	d := time.Since(t)
	fmt.Printf("\nTime elapse: %v sec.\n", d)
}
