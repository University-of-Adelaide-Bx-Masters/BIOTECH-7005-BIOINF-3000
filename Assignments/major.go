package main

import (
	"bufio"
	"fmt"
	"log"
	"math/rand"
	"os"
	"sort"
	"strconv"
	"strings"
)

var faSynthesisData = []string{
	"fae1_7-8",
	"fae1_9-10",
	"fae1_11-12",
	"CL37_7-8",
	"CL37_9-10",
	"CL37_11-12",
	"PDAT_7-8",
	"PDAT_9-10",
	"PDAT_11-12",
}

var arabidopsisData = []string{
	"hda5-1 vs wild-type (arabidopsis dataset)",
	"hda6-6 vs wild-type (arabidopsis dataset)",
	"hda9-1 vs wild-type (arabidopsis dataset)",
}

func pairwiseDatasets(data []string) []string {
	var (
		mutantSet = make(map[string]bool)
		timeSet   = make(map[string]bool)
	)
	for _, d := range data {
		f := strings.Split(d, "_")
		mutantSet[f[0]] = true
		timeSet[f[1]] = true
	}

	var mutants []string
	for m := range mutantSet {
		mutants = append(mutants, m)
	}
	sort.Strings(mutants)
	var times []string
	for t := range timeSet {
		times = append(times, t)
	}
	sort.Strings(times)

	var options []string
	for _, t := range times {
		for i, m1 := range mutants {
			for _, m2 := range mutants[i+1:] {
				options = append(options, fmt.Sprintf("%s vs %s at %s (fatty acid synthesis dataset)", m1, m2, t))
			}
		}
	}
	for _, m := range mutants {
		for i, t1 := range times {
			for _, t2 := range times[i+1:] {
				options = append(options, fmt.Sprintf("%s at %s vs %s (fatty acid synthesis dataset)", m, t1, t2))
			}
		}
	}

	return options
}

func main() {
	option := append(pairwiseDatasets(faSynthesisData), arabidopsisData...)
	sc := bufio.NewScanner(os.Stdin)
	for sc.Scan() {
		fields := strings.Fields(sc.Text())
		if len(fields) == 0 {
			continue
		}
		if fields[0][0] != 'a' {
			continue
		}
		a, err := strconv.Atoi(fields[0][1:])
		if err != nil {
			log.Fatalf("failed to parse id: %v", err)
		}
		rand.Seed(int64(a))
		fmt.Printf("|%s|%s|\n", fields[0], option[rand.Intn(len(option))])
	}
	err := sc.Err()
	if err != nil {
		log.Fatalf("error during scan: %v", err)
	}

}
