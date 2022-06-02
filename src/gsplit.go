package main

import (
	"fmt"
	"io"
	"os"
	"regexp"
	"runtime"
	"strings"
	"time"

	"github.com/hpifu/go-kit/hflag"
	"github.com/shenwei356/xopen"
)

const VERSION = "v2022.06.02 15:00"

type SplitFlags struct {
	Fqfile      []string `hflag:"--input, -i; required; usage: input fastq file, *.gz/xz/zst or uncompress allowed, multi-input call be separated by ',', required"`
	Barcodefile string   `hflag:"--barcode, -b; required; usage: barcode and sample file, required"`
	Threads     int      `hflag:"--threads, -t; default: 10; usage: threads core, 10 by default"`
	Mismatch    int      `hflag:"--mismatch, -m; default: 0; usage: mismatch allowed for barcode search, 0 by default"`
	Output      string   `hflag:"--output, -o; required; usage: output directory, required"`
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

func ReadSeq(ch1 chan []string, infile string) {

	r, err := xopen.Ropen(infile)
	checkError(err)
	defer r.Close()
	seq := make([]string, 4, 5)
	n := 0
	for {
		line, err := r.ReadString('\n')
		if err == io.EOF {
			break
		}
		i := n % 4
		seq[i] = line
		if i == 3 {
			ch1 <- seq
			seq = make([]string, 4, 5)
		}
		n++
	}
	close(ch1)
}

func SplitSeq(ch1, ch2 chan []string, exch1 chan bool, barcode map[string]string, mis int, drup bool) {
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
			v = append(v, sn)
			ch2 <- v
			break bcloop
		}
	}
	exch1 <- true
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

/*  run channel for results in pool
func main() {

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

	bcf, _ := xopen.Ropen(args.Barcodefile)
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

	ch1 := make(chan []string, 10)
	res := make(chan []string, 5)
	exch := make(chan bool, 1)

	go ReadSeq(ch1, args.Fqfile)

	for i := 1; i <= args.Threads; i++ {
		go SplitSeq(ch1, res, exch, barcode, args.Mismatch, args.Drup)
	}

	go func() {
		for i := 1; i <= args.Threads; i++ {
			<-exch
		}
		close(res)
	}()

	for out := range res {
		sn := out[4]
		for _, line := range out[:4] {
			fout[sn].WriteString(line)
		}

	}
	d := time.Since(t)
	fmt.Printf("Time elapse: %v sec.\n", d)
}
*/

func main() {

	args := &SplitFlags{}
	if err := hflag.Bind(args); err != nil {
		panic(err)
	}
	err := hflag.Parse()
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

	bcf, _ := xopen.Ropen(args.Barcodefile)
	defer bcf.Close()

	fout := make(map[string]*xopen.Writer, 10)

	if !isExist(args.Output) {
		err := os.MkdirAll(args.Output, os.ModePerm)
		checkError(err)
	}
	drup_pos := make(map[string]int, 10)
	samples := []string{}
	for {
		line, err := bcf.ReadString('\n')
		if err == io.EOF {
			break
		}
		reg := regexp.MustCompile(`\s+`)
		line_s := reg.Split(strings.TrimSpace(line), -1)
		sn, bc := line_s[0], line_s[1]
		barcode[bc] = sn
		samples = append(samples, sn)
		outf := args.Output + "/" + sn + ".fq"
		if !args.Nogz {
			outf += ".gz"
		}
		if args.Drup {
			drup_pos[bc] = len(bc)
		} else {
			drup_pos[bc] = 0
		}
		fo, _ := xopen.Wopen(outf)
		defer fo.Close()
		fout[line_s[0]] = fo
	}
	/*
		unkonwfile := args.Output + "/" + "Unknow" + ".fq"
		if !args.Nogz {
			unkonwfile += ".gz"
		}
		unfo, _ := xopen.Wopen(unkonwfile)
		defer unfo.Close()
		fout["Unknow"] = unfo
	*/
	mis := args.Mismatch
	seq := [4]string{}
	sms := make(map[string]int, len(barcode))
	total := 0
	for _, fqfile := range args.Fqfile {
		r, err := xopen.Ropen(fqfile)
		checkError(err)
		defer r.Close()
		linefo := 0
		for {
			line, err := r.ReadString('\n')
			if err == io.EOF {
				break
			}
			i := linefo % 4
			seq[i] = line
			// unknow := true
			if i == 3 {
			BCLOOP:
				for b, sn := range barcode {
					d := 0
					for n, _ := range b {
						if seq[1][n] != b[n] {
							d++
						}
						if d > mis {
							continue BCLOOP
						}
					}
					// unknow = false
					seq[1] = seq[1][drup_pos[b]:]
					seq[3] = seq[3][drup_pos[b]:]
					sms[sn] += 1
					for _, line := range seq {
						fout[sn].WriteString(line)
					}
					break BCLOOP
				}
				/*
				   if unknow {
				       for _, line := range seq {
				           fout["Unknow"].WriteString(line)
				       }
				   }
				*/
				total += 1
			}
			linefo++
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
