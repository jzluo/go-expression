package main

import (
	"fmt"
	"gonum.org/v1/gonum/mat"
	"os"
	"runtime"
	"strconv"
)

type Tree []*Node

type Cluster *Node

type Node struct {
	height         float64
	avg            float64
	label          string
	child1, child2 *Node
}

type MatrixElement struct {
	row   int
	col   int
	value float64
}

func TOM(adjacencyMatrix [][]float64) [][]float64 {
	/*
		Constructs a topological overlap matrix (TOM)from a weighted adjacency matrix.
		w_ij = (l_ij + a_ij) / (min{k_i, k_j} + 1 - a_ij)  eq 4 in WGCNA framework paper

		Splits the adjacency matrix into numProcs sub-matrices to construct the matrix in parallel
	*/
	tom := make([][]float64, len(adjacencyMatrix))
	for i := range tom {
		tom[i] = make([]float64, len(adjacencyMatrix))
		tom[i][i] = 1 // set diagonal to 1
	}

	numProcs := runtime.NumCPU()
	chunks := len(adjacencyMatrix) / numProcs
	c := make(chan int)

	for i := 0; i < numProcs; i++ {
		start := i * chunks
		end := (i + 1) * chunks
		if i == numProcs-1 {
			end = len(adjacencyMatrix) - 1
			go func() {
				//	fmt.Println("starting last goroutine", start, end)
				BuildTOM(adjacencyMatrix, tom, start, end, c)
			}()
		} else {
			go func() {
				//	fmt.Println("starting goroutine ", start, end)
				BuildTOM(adjacencyMatrix, tom, start, end, c)
			}()
		}
	}
	countGR := 0
	for i := 0; i < numProcs; i++ {
		countGR += <-c
		fmt.Printf("\r\tSub-matrices completed: %d/%d", countGR, numProcs)
	}
	fmt.Println()
	//
	//for i := range adjacencyMatrix {
	//	//a_iu := adjacencyMatrix[i]
	//	gene_a := mat.NewVecDense(len(adjacencyMatrix[0]), adjacencyMatrix[i])
	//	k_i := SumRow(adjacencyMatrix[i])
	//	for j := i + 1; j < len(adjacencyMatrix[0]); j++ { // since it's symmetric, calculate only for half of matrix
	//		//if tom[i][j] != 0 {
	//		//	continue
	//		//}
	//
	//		/*
	//			w_ij = (l_ij + a_ij) / (min{k_i, k_j} + 1 - a_ij)  eq 4 in WGCNA framework paper
	//			wnere w_ij = topological overlap of nodes i and j
	//			a_ij = direct connection between nodes i and j
	//			l_ij = \sum_u a_iu a_uj
	//			k_i = \sum_u a_iu "connectivity" of node i
	//		*/
	//		gene_b := mat.NewVecDense(len(adjacencyMatrix[0]), adjacencyMatrix[j])
	//		l_ij := mat.Dot(gene_a, gene_b)
	//		a_ij := adjacencyMatrix[i][j]
	//		a_ju := adjacencyMatrix[j]
	//		//l_ij := Dot(a_iu, a_ju)
	//		k_j := SumRow(a_ju)
	//		k_min := Min(k_i, k_j)
	//
	//		w_ij := (l_ij + a_ij) / (k_min + 1 - a_ij)
	//
	//		tom[i][j] = w_ij
	//		tom[j][i] = w_ij // s y m m e t r y
	//		fmt.Println(i, j)
	//	}
	//}
	return tom
}

func BuildTOM(AdjacencyMatrix [][]float64, tom [][]float64, start, end int, c chan int) {
	/*
		Builds the TOM for a provided submatrix. Uses gonum to perform the dot products instead of Dot() for performance
		for large sets of genes.

		w_ij = (l_ij + a_ij) / (min{k_i, k_j} + 1 - a_ij)  eq 4 in WGCNA framework paper
		wnere w_ij = topological overlap of nodes i and j
		a_ij = direct connection between nodes i and j
		l_ij = \sum_u a_iu a_uj
		k_i = \sum_u a_iu "connectivity" of node i
	*/

	//fmt.Println(len(AdjacencyMatrix), len(AdjacencyMatrix[0]), start)
	for i := start; i < end; i++ {
		gene_a := mat.NewVecDense(len(AdjacencyMatrix[0]), AdjacencyMatrix[i])
		k_i := SumRow(AdjacencyMatrix[i])
		for j := start + 1; j < len(AdjacencyMatrix[0]); j++ { // since it's symmetric, calculate only for half of matrix
			//if tom[i][j] != 0 {
			//	continue
			//}

			gene_b := mat.NewVecDense(len(AdjacencyMatrix[0]), AdjacencyMatrix[j])
			l_ij := mat.Dot(gene_a, gene_b)
			a_ij := AdjacencyMatrix[i][j]
			a_ju := AdjacencyMatrix[j]
			//l_ij := Dot(a_iu, a_ju)
			k_j := SumRow(a_ju)
			k_min := Min(k_i, k_j)

			w_ij := (l_ij + a_ij) / (k_min + 1 - a_ij)

			tom[i][j] = w_ij
			tom[j][i] = w_ij // s y m m e t r y
			//fmt.Println(i, j)
		}
	}
	c <- 1
}

func Dot(vi, vj []float64) float64 {
	// computes dot product of two slices of floats
	var product float64
	for idx, val := range vi {
		product += val * vj[idx]
	}
	return product
}

func SumRow(row []float64) float64 {
	// computes sum of a slice of floats
	var sum float64
	for _, val := range row {
		sum += val
	}
	return sum
}

//func CountSharedNeighbors(AdjacencyMatrix, i, j) int {
//	// Counts nodes that nodes i and j both have connections to
//	i_neighbors := make([]int, 0)
//	var numSharedNeighbors int
//	for u := range AdjacencyMatrix[i] {
//		if AdjacencyMatrix[i][u] == 1 {
//			i_neighbors = append(i_neighbors, u)
//		}
//	}
//
//	for _, u := range i_neighbors {
//		if AdjacencyMatrix[j][u] == 1 {
//			numSharedNeighbors++
//		}
//	}
//	return numSharedNeighbors
//}

//func CountNeighbors(AdjacencyMatrix, i) int {
//	// Count total neighbors for node i (overall "connectivity")
//	var numNeighbors int
//	for u := range AdjacencyMatrix[i] {
//		if AdjacencyMatrix[i][u] == 1 {
//			numNeighbors++
//		}
//	}
//	return numNeighbors
//}
//

func Min(a, b float64) float64 {
	if a < b {
		return a
	} else {
		return b
	}
}

func DissTOM(tom [][]float64) [][]float64 {
	// return a "dissimilarity" matrix from the TOM for clustering.
	//writer, _ := os.Create("dissTOM.csv")
	//defer writer.Close()

	for i := range tom {
		//fmt.Fprintln(writer)
		for j := range tom[i] {
			tom[i][j] = 1 - tom[i][j]
			//fmt.Fprintf(writer, "%f,", tom[i][j])
		}
	}
	return tom
}

//func ReadDissTOM(fileName string) [][]float64 {
//	file, err := os.Open(fileName)
//	if err != nil {
//		log.Fatal(err)
//		panic("Error: Issue Opening Raw Data Files.")
//	}
//	defer file.Close()
//	scanner := bufio.NewScanner(file)
//
//	dissTOM := make([][]float64, 12700)
//	row := make([]float64 ,0)
//	counter := 0
//	for scanner.Scan() {
//		//read current line
//		currentLine := scanner.Text()
//		// split the line on delimiter
//		splitLine := strings.Split(currentLine, ",")
//		for i := range splitLine {
//			x, err := strconv.ParseFloat(splitLine[i], 64)
//			if err == nil {
//				fmt.Println(splitLine[i])
//				panic("unknown float?")
//			}
//			row = append(row, x)
//		}
//		dissTOM[counter] = row
//		counter++
//	}
//	return dissTOM
//}

func UPGMA(dissTOM [][]float64, genes []string) Tree {
	/*
	 Performs hierarchical clustering on the TOM dissimilarity matrix, returning a tree.
	*/

	t := InitializeTree(genes)

	clusters := t.InitializeClusters() // all genes are leaves in beginning

	numLeaves := len(genes)
	for k := numLeaves; k < 2*numLeaves-1; k++ {
		//row, col, currentMin := NextMerge(dissTOM)  // find min element in parallel since it uses most of the runtime
		row, col, currentMin := GetMinimumElementSingleproc(dissTOM)
		t[k].height = currentMin // height of the node in the dendrogram, proportional to dissimilarity of children

		t[k].child1 = CompareNodes(clusters[row], clusters[col])[1] // place longer one on left
		t[k].child2 = CompareNodes(clusters[row], clusters[col])[0]

		dissTOM = AddRowCol(dissTOM, clusters, row, col)
		clusters = append(clusters, t[k])

		dissTOM = DelRowCol(dissTOM, row, col)
		clusters = DelClusters(clusters, row, col)
	}
	//fmt.Println(len(clusters))
	return t
}

func CompareNodes(a, b *Node) []*Node {
	if a.height < b.height {
		return []*Node{a, b}
	} else {
		return []*Node{b, a}
	}

}

func InitializeTree(genes []string) Tree {
	// Initialize a tree where in its initial state, each gene is a node (leaf)
	numLeaves := len(genes)

	// total number of nodes is 2*numLeaves - 1
	var tree Tree = make([]*Node, 2*numLeaves-1)

	for i := range tree {
		var vx Node

		if i < numLeaves {
			// if leaf, label with the gene name
			vx.label = genes[i]
		} else {
			vx.label = "Gene cluster " + strconv.Itoa(i)
		}
		// point tree[i] to the node vx
		tree[i] = &vx
	}
	return tree
}

func (t Tree) InitializeClusters() []Cluster {
	// returns a slice of pointers to the leaves of the Tree
	numNodes := len(t)
	numLeaves := (numNodes + 1) / 2

	clusters := make([]Cluster, numLeaves)

	for i := 0; i < numLeaves; i++ {
		clusters[i] = t[i]
	}
	return clusters
}

func NextMerge(dissTOM [][]float64) (int, int, float64) {
	// Coordinate parallelizing of finding minimum in a matrix by splitting into sub-matrices, then getting min of results
	numProcs := runtime.NumCPU()
	chunks := len(dissTOM) / numProcs
	c := make(chan MatrixElement)

	//if chunks == 0 {
	//	numProcs = 1
	//}

	for i := 0; i < numProcs; i++ {
		start := i * chunks
		end := (i + 1) * chunks
		if i == numProcs-1 {
			end = len(dissTOM) - 1
			go func() {
				//	fmt.Println("starting last goroutine", start, end)
				GetMinElement(dissTOM, start, end, c)
			}()
		} else {
			go func() {
				//	fmt.Println("starting goroutine ", start, end)
				GetMinElement(dissTOM, start, end, c)
			}()
		}
	}

	minimums := make([]MatrixElement, 0)
	for i := 0; i < numProcs; i++ {
		minimums = append(minimums, <-c)
	}

	return GetMESliceMin(minimums)
}

func GetMinElement(dissTOM [][]float64, start, end int, c chan MatrixElement) { //(int, int, float64) {
	// finds minimum element of sub-matrix
	fmt.Println(len(dissTOM), len(dissTOM[0]), start)
	minElement := MatrixElement{value: 999.0}

	for i := start; i < end; i++ {
		for j := start + 1; j < len(dissTOM[0]); j++ {
			if dissTOM[i][j] < minElement.value {
				minElement.row = i
				minElement.col = j
				minElement.value = dissTOM[i][j]
			}
		}

	}
	c <- minElement

}

func GetMinimumElementSingleproc(dissTOM [][]float64) (int, int, float64) {
	row := 0
	col := 1
	currentMin := dissTOM[row][col]

	for i := 0; i < len(dissTOM); i++ {
		for j := i + 1; j < len(dissTOM[i]); j++ {
			if dissTOM[i][j] < currentMin {
				row = i
				col = j
				currentMin = dissTOM[i][j]
			}
		}
	}

	return row, col, currentMin
}

func GetMESliceMin(minimums []MatrixElement) (int, int, float64) {
	// find minimum of a slice of MatrixElement

	min := MatrixElement{
		row:   minimums[0].row,
		col:   minimums[0].col,
		value: minimums[0].value,
	}
	for i := range minimums {
		if minimums[i].value < min.value {
			min.row = minimums[i].row
			min.col = minimums[i].col
			min.value = minimums[i].value
		}
	}
	return min.row, min.col, min.value
}

func AddRowCol(dissTOM [][]float64, clusters []Cluster, row, col int) [][]float64 {
	n := len(dissTOM)
	newRow := make([]float64, n+1)

	for j := 0; j < n; j++ {
		if j != row && j != col {

			// Use weighted average based on number of genes in cluster

			size1 := CountLeaves(clusters[row])
			size2 := CountLeaves(clusters[col])

			newRow[j] = (float64(size1)*dissTOM[row][j] + float64(size2)*dissTOM[col][j]) / float64(size1+size2)
		}
	}

	// append newRow to matrix
	dissTOM = append(dissTOM, newRow)

	// add last column to each row
	for i := 0; i < n; i++ {
		dissTOM[i] = append(dissTOM[i], newRow[i])
	}

	return dissTOM
}

func CountLeaves(vx *Node) int {
	// Counts the number of leaves (genes) in the tree rooted at the node.
	if vx.child1 == nil || vx.child2 == nil { // we are at a leaf
		return 1
	}
	// inductive step: call CountLeaves for children
	return CountLeaves(vx.child1) + CountLeaves(vx.child2)
}

func DelRowCol(dissTOM [][]float64, row, col int) [][]float64 {
	dissTOM = append(dissTOM[:col], dissTOM[col+1:]...)
	dissTOM = append(dissTOM[:row], dissTOM[row+1:]...)

	for i := range dissTOM {
		dissTOM[i] = append(dissTOM[i][:col], dissTOM[i][col+1:]...)
		dissTOM[i] = append(dissTOM[i][:row], dissTOM[i][row+1:]...)
	}

	return dissTOM
}

func DelClusters(clusters []Cluster, row, col int) []Cluster {
	clusters = append(clusters[:col], clusters[col+1:]...)
	clusters = append(clusters[:row], clusters[row+1:]...)
	return clusters
}

func WriteClusters(fileName string, clusterRoots []*Node) {
	// Writes the genes belonging to each cluster found to fileName, one cluster per line in tab-delimited format
	outfile, err := os.Create(fileName)
	if err != nil {
		fmt.Println("Could not create file", fileName)
		os.Exit(0)
	}
	defer outfile.Close()
	for i := range clusterRoots {
		if clusterRoots[i].label != "notnode" {
			labels := make([]string, 100000) // can't get it to recursively append to slice
			index := []int{0}
			TraverseCluster(clusterRoots[i], labels, index)
			//if node[i].child1 != nil || node[i].child2 != nil {
			//	fmt.Println("child", node[i].child1.label)
			//	fmt.Println("child", node[i].child2.label)
			//}
			labels = GetGeneLabels(labels)
			//fmt.Println("Number of genes in cluster ", i, ":", len(labels), clusterRoots[i].height)
			if len(labels) > 1 {
				for j := range labels {
					fmt.Fprintf(outfile, "%s\t", labels[j])
				}
				fmt.Fprintf(outfile, "\n")
			}
		}
	}
	//fmt.Println(len(clusterRoots))

}

func RemoveNullNodes(nodes []*Node) []*Node {
	for i, gene := range nodes {
		if gene != nil {
			return nodes[0 : i+1]
		}
	}
	return nodes
}

func GetGeneLabels(labels []string) []string {
	for i, gene := range labels {
		if gene == "" {
			return labels[0:i]
		}
	}
	return labels
}

func TraverseCluster(vx *Node, labels []string, index []int) {
	// recursively traverse the tree until reaching the genes (leaves)
	//fmt.Println(index[0])
	if vx.child1 == nil && vx.child2 == nil {
		labels[index[0]] = vx.label
		update := index[0] + 1
		index[0] = update
		return
	}

	if vx.child1 != nil {
		TraverseCluster(vx.child1, labels, index)
	}
	if vx.child2 != nil {
		TraverseCluster(vx.child2, labels, index)
	}
}

func TreeCut(vx *Node, nodes []*Node, index []int) { //(vx *Node) {
	// cuts the tree at a given height to give gene clusters
	// Recursively traverses the tree and returns the first nodes below the cut, which are the cluster roots
	if vx.height < 0.94 {
		//fmt.Println(vx.label, vx.age)
		////nodes = append(nodes, vx) //append([]*Node{},vx)...)
		//ages = append(ages, vx.age)
		//names = append(names, vx.label)
		nodes[index[0]] = vx
		update := index[0] + 1
		index[0] = update
		//fmt.Println(nodes)
		return // add the node to the list of clusters and stop traversing
	}

	TreeCut(vx.child1, nodes, index)
	TreeCut(vx.child2, nodes, index)
}

func NonRecursiveCut(tree Tree) []*Node {
	cut := 10.0
	children := make([]*Node, 0)
	roots := make([]*Node, 0)
	for i := len(tree) - 1; i > 0; i-- {
		if tree[i].height < cut {
			if len(roots) == 0 {
				roots = append(roots, tree[i])
				children = append(children, []*Node{tree[i].child1, tree[i].child2}...)
				continue
			}
			keep := true
			for i := range children {
				if tree[i] == children[i] {
					keep = false
					break
				}
			}
			if keep {
				roots = append(roots, tree[i])
				children = append(children, []*Node{tree[i].child1, tree[i].child2}...)
			}
		}
	}
	return roots
}
