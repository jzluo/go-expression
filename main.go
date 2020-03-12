package main

import (
	//"bufio"
	//"flag"
	//"log"
	//"os"
	"fmt"
	"runtime"
	//"runtime/pprof"
	//"strings"
	"time"
)

//var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to `file`")
//var memprofile = flag.String("memprofile", "", "write memory profile to `file`")

func main() {
	//flag.Parse()
	//if *cpuprofile != "" {
	//	f, err := os.Create(*cpuprofile)
	//	if err != nil {
	//		log.Fatal("could not create CPU profile: ", err)
	//	}
	//	defer f.Close()
	//	if err := pprof.StartCPUProfile(f); err != nil {
	//		log.Fatal("could not start CPU profile: ", err)
	//	}
	//	defer pprof.StopCPUProfile()
	//}

	numProcs := runtime.NumCPU()
	fmt.Println(numProcs, "threads detected")
	fmt.Println()

	fmt.Println("----------------------------------------")
	fmt.Println("CO-EXPRESSION GRAPHS AND THEIR ANALYSIS")
	fmt.Println("----------------------------------------")

	fmt.Println("For the covariance tests, please enter P for Pearsons (faster) or B for BiWeightedCorrelation:")
	var correlation string
	fmt.Scanln(&correlation)
	for correlation != "B" && correlation != "P" && correlation != "b" && correlation != "p" {
		fmt.Println("Please input B or P")
		fmt.Scanln(&correlation)
	}

	fmt.Println("Enter S to build a Signed correlation network or U for Unsigned:")
	var sign string
	fmt.Scanln(&sign)
	for sign != "S" && sign != "U" && sign != "s" && sign != "u" {
		fmt.Println("Please input S or U")
		fmt.Scanln(&sign)
	}

	fmt.Println(time.Now().Format("2006-01-02 15:04:05"), "doing IO")
	filteredData := io()

	fmt.Println(time.Now().Format("2006-01-02 15:04:05"), "Building expression matrix")
	expMatrix, genes := ConvertMaptoGeneExpressionMatrix(filteredData)

	//fmt.Println(len(expMatrix))
	//write genes
	//writer, _ := os.Create("genes.csv")
	//defer writer.Close()
	//for i := range genes {
	//	fmt.Fprintf(writer, "%s,", genes[i])
	//}

	fmt.Println("----------------------------------------")
	fmt.Println(time.Now().Format("2006-01-02 15:04:05"), "Building covariance matrix")
	covarianceMatrix := make([][]float64, len(genes))
	for i := range genes {
		covarianceMatrix[i] = make([]float64, len(genes))
	}
	if correlation == "P" || correlation == "p" {
		covarianceMatrix = BuildCovarianceMatrixPearson(expMatrix, genes)
	} else if correlation == "B" || correlation == "b" {
		covarianceMatrix = BuildCovarianceMatrixBiWeight(expMatrix, genes)
	}

	fmt.Println("----------------------------------------")
	fmt.Println(time.Now().Format("2006-01-02 15:04:05"), "Building adjacency matrix")

	corrMatrix := make([][]float64, len(genes))
	for i := range genes {
		corrMatrix[i] = make([]float64, len(genes))
	}
	if sign == "S" || sign == "s" {
		corrMatrix = BuildSignedWeightedCorrelationNetwork(covarianceMatrix)
	} else if sign == "U" || sign == "u" {
		corrMatrix = BuildUnsignedWeightedCorrelationNetwork(covarianceMatrix)
	}
	fmt.Println("----------------------------------------")
	fmt.Println(time.Now().Format("2006-01-02 15:04:05"), "Building TOM")
	tom := TOM(corrMatrix)
	fmt.Println("----------------------------------------")
	fmt.Println(time.Now().Format("2006-01-02 15:04:05"), "Building DissTOM")
	//
	//
	dissTOM := DissTOM(tom)

	// read genes
	//file, err := os.Open("genes.csv")
	//if err != nil {
	//	log.Fatal(err)
	//	panic("Error: Issue Opening Raw Data Files.")
	//}
	//defer file.Close()
	//scanner := bufio.NewScanner(file)
	//scanner.Scan()
	//genes := make([]string, 0)
	//splitLine := strings.Split(scanner.Text(), ",")
	//for i := range splitLine {
	//	genes = append(genes, splitLine[i])
	//}

	//dissTOM := ReadDissTOM("dissTOM.csv")
	//fmt.Println(dissTOM)
	fmt.Println("----------------------------------------")
	fmt.Println(time.Now().Format("2006-01-02 15:04:05"), "Performing hierarchical clustering")
	tree := UPGMA(dissTOM, genes)

	fmt.Println(time.Now().Format("2006-01-02 15:04:05"), "Writing clusters to clusters.txt")
	nodes := make([]*Node, 100000) // can't get recursively appending to a slice to work
	for i := range nodes {
		newnode := Node{0.0, 0.0, "notnode",
			nil, nil}
		nodes[i] = &newnode
	}
	index := []int{0}
	//PrintTree(tree[len(tree)-1])
	TreeCut(tree[len(tree)-1], nodes, index)

	//fmt.Println(nodes[1].label)
	//fmt.Println(len(nodes))
	//fmt.Println(tree[len(tree)-1].label, tree[len(tree)-1].height)
	//fmt.Println(nodes[0].label, nodes[0].height)

	//clusterRoots := NonRecursiveCut(tree)

	//nodes=RemoveNullNodes(nodes)
	WriteClusters("clusters.txt", nodes)
	fmt.Println("Clusters written to clusters.txt")

	//if *memprofile != "" {
	//	f, err := os.Create(*memprofile)
	//	if err != nil {
	//		log.Fatal("could not create memory profile: ", err)
	//	}
	//	defer f.Close()
	//	runtime.GC() // get up-to-date statistics
	//	if err := pprof.WriteHeapProfile(f); err != nil {
	//		log.Fatal("could not write memory profile: ", err)
	//	}
	//}

}
