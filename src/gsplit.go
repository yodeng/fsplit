package main

import (
	"fmt"
	"io"
	"os"
	"regexp"
	"runtime"
	"strings"
	"sync"
	"time"

	"github.com/yodeng/go-kit/hflag"
	"github.com/yodeng/xopen"
)

const VERSION = "v2022.07.08 15:00"

type SplitFlags struct {
	Fqfile      []string `hflag:"--input, -i; required; usage: input fastq file, *.gz/xz/zst or uncompress allowed, multi-input can be separated by ',' or whitespace, required"`
	Fqfile2     []string `hflag:"--Input, -I; usage: input read2 fastq file if there is, *.gz/xz/zst or uncompress allowed, multi-input can be separated by ',' or whitespace"`
	Barcodefile string   `hflag:"--barcode, -b; required; usage: barcode and sample file, 1st column for sample name and 2nd column for barcode sequence, required"`
	Output      string   `hflag:"--output, -o; required; usage: output directory, will create if not exists, required"`
	Threads     int      `hflag:"--threads, -t; default: 10; usage: threads core, 10 by default"`
	Mismatch    int      `hflag:"--mismatch, -m; default: 0; usage: mismatch allowed for barcode search, 0 by default"`
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

/*  run channel for results in pool
func main() {
	var wg sync.WaitGroup
	args := &SplitFlags{}
	if err := hflag.Bind(args); err != nil {
		panic(err)
	}
	if err := hflag.Parse(); err != nil {
		fmt.Println(hflag.Usage())
		panic(err)
	}

	runtime.GOMAXPROCS(args.Threads)
	t := time.Now()
	barcode := make(map[string]string, 10)

	bcf, err := xopen.Ropen(args.Barcodefile)
	checkError(err)
	defer bcf.Close()

	fout := make(map[string]*xopen.Writer, 10)

	if !isExist(args.Output) {
		err := os.MkdirAll(args.Output, os.ModePerm)
		checkError(err)
	}

	for {
		line, err := bcf.ReadString('\n')
		if err == io.EOF {
			break
		}
		reg := regexp.MustCompile(`\s+`)
		line_s := reg.Split(strings.TrimSpace(line), -1)
		sn, bc := line_s[0], line_s[1]
		barcode[bc] = sn
		outf := args.Output + "/" + sn + ".fq"
		if !args.Nogz {
			outf += ".gz"
		}
		fo, _ := xopen.Wopen(outf)
		defer fo.Close()
		fout[line_s[0]] = fo
	}

	ch1 := make(chan []string, 10000)
	res := make(chan []string, 10000)

	for _, f := range args.Fqfile {
		go ReadSeq(ch1, f)
	}

	for i := 1; i <= args.Threads; i++ {
		wg.Add(1)
		go SplitSeq(ch1, res, barcode, args.Mismatch, args.Drup, &wg)
	}
	go func() {
		wg.Wait()
		close(res)
	}()

	// sms := map[string]int{}

	var out sync.WaitGroup

	locks := make(map[string]*sync.Mutex)
	for _, sn := range barcode {
		locks[sn] = &sync.Mutex{} // lock for each file writing
	}

	for i := 0; i < len(barcode); i++ {
		out.Add(1)
		go func() {
			defer out.Done()
			for out := range res {
				sn := out[4]
				locks[sn].Lock()
				for _, line := range out[:4] {
					fout[sn].WriteString(line)
				}
				// sms[sn] += 1  // conrutine write map error
				locks[sn].Unlock()

			}
		}()
	}
	out.Wait()
	d := time.Since(t)

	// fmt.Println(sms)

	fmt.Printf("Time elapse: %v sec.\n", d)
}
*/

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
	barcode := make(map[string]string, 10)

	bcf, err := xopen.Ropen(args.Barcodefile)
	checkError(err)
	defer bcf.Close()

	fout := make(map[string][]*xopen.Writer, 10)

	if !isExist(args.Output) {
		err := os.MkdirAll(args.Output, os.ModePerm)
		checkError(err)
	}
	drup_pos := make(map[string]int, 10)
	samples := []string{}
	if args.Fqfile2 != nil {
		if len(args.Fqfile2) != len(args.Fqfile) {
			fmt.Println("Miss input fastq file of R1 or R2")
			os.Exit(1)
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
		sn, bc := line_s[0], line_s[1]
		barcode[bc] = sn
		if len(line_s) >= 3 {
			barcode[line_s[2]] = sn
		}
		if args.Drup {
			drup_pos[bc] = len(bc)
		}
		if !inSlice(samples, sn) {
			samples = append(samples, sn)
			if args.Fqfile2 != nil {
				outf1 := args.Output + "/" + sn + ".R1.fq"
				outf2 := args.Output + "/" + sn + ".R2.fq"
				if !args.Nogz {
					outf1 += ".gz"
					outf2 += ".gz"
				}
				fo1, _ := xopen.Wopen(outf1)
				defer fo1.Close()
				fo2, _ := xopen.Wopen(outf2)
				defer fo2.Close()
				fout[sn] = append(fout[sn], []*xopen.Writer{fo1, fo2}...)
			} else {
				outf := args.Output + "/" + sn + ".fq"
				if !args.Nogz {
					outf += ".gz"
				}
				fo, _ := xopen.Wopen(outf)
				defer fo.Close()
				fout[sn] = append(fout[sn], fo)
			}

		}
	}
	mis := args.Mismatch
	seq := [4]string{}
	sms := make(map[string]int, len(barcode))
	sample2barcode := make(map[string][]string, len(samples))
	for b, sn := range barcode {
		sample2barcode[sn] = append(sample2barcode[sn], b)
	}
	total := 0
	for fn, fqfile := range args.Fqfile {
		r, err := xopen.Ropen(fqfile)
		checkError(err)
		defer r.Close()
		linefo := 0
		if args.Fqfile2 != nil {
			r2, err := xopen.Ropen(args.Fqfile2[fn])
			checkError(err)
			defer r2.Close()
			seq2 := [4]string{}
			for {
				line, err := r.ReadString('\n')
				if err == io.EOF {
					break
				}
				line2, _ := r2.ReadString('\n')
				i := linefo % 4
				seq[i] = line
				seq2[i] = line2
				if i == 3 {
					sn, b1 := splitSample(seq[1], barcode, mis)
					if sn != "" {
						b2 := splitBarcode(seq2[1], sample2barcode[sn], mis)
						if b2 != "" {
							seq[1] = seq[1][drup_pos[b1]:]
							seq[3] = seq[3][drup_pos[b1]:]
							seq2[1] = seq2[1][drup_pos[b2]:]
							seq2[3] = seq2[3][drup_pos[b2]:]
							sms[sn] += 1
							for _, line := range seq {
								fout[sn][0].WriteString(line)
							}
							for _, line := range seq2 {
								fout[sn][1].WriteString(line)
							}
						}
					}
					total += 1
				}
				linefo++
			}
		} else {
			for {
				line, err := r.ReadString('\n')
				if err == io.EOF {
					break
				}
				i := linefo % 4
				seq[i] = line
				if i == 3 {
				BCLOOP:
					for b, sn := range barcode {
						d := 0
						for n := range b {
							if seq[1][n] != b[n] {
								d++
							}
							if d > mis {
								continue BCLOOP
							}
						}
						seq[1] = seq[1][drup_pos[b]:]
						seq[3] = seq[3][drup_pos[b]:]
						sms[sn] += 1
						for _, line := range seq {
							fout[sn][0].WriteString(line)
						}
						break BCLOOP
					}
					total += 1
				}
				linefo++
			}
		}
	}
	fmt.Println()
	snm := 0
	for _, sn := range samples {
		count := sms[sn]
		snm += count
		fmt.Printf("%v: %v(%.2f%%)\n", sn, count, float64(count)/float64(total)*100)
	}
	fmt.Printf("Unknow: %v(%.2f%%)\n", total-snm, float64(total-snm)/float64(total)*100)
	d := time.Since(t)
	fmt.Printf("\nTime elapse: %v sec.\n", d)
}
