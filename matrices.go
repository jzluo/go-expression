package main

import (
	//"log"
	"fmt"
	"math"
	"math/rand"
	"os"
	//"time"
	//"github.com/go-echarts/go-echarts/charts"
)

type edge struct {
	source, target string
	weight         float64
}

/*
Input: Struct of type edge that contains information involving source, destination gene and edge weight
Output: Writes out the structure onto a .tsv file to export to Cytoscape for Protein Network Visualization
*/

func WriteAdjacencyMatrixtoFile(goods []edge) {
	writer, _ := os.Create("adjacency_matrix.tsv") //Writes in file writer
	defer writer.Close()
	fmt.Fprintln(writer, "Source\tTarget\tWeight")
	for i := range goods {
		fmt.Fprintf(writer, "%s\t%s\t%f\n", goods[i].source, goods[i].target, goods[i].weight)
	}
}

/*
Input: Map[string] Map[string] float64
Mouse : Genes : Counts
Output: 2D Slice - Mice - Rows, Columns - Genes counts
*/

func ConvertMaptoGeneExpressionMatrix(MiceGeneCounts map[string]map[string]float64) ([][]float64, []string) {
	CoexpressionMatrix := make([][]float64, len(MiceGeneCounts))
	GeneLabels := make([]string, 0)
	index := 0
	first_time := true
	for _, mousegenemap := range MiceGeneCounts {
		row := make([]float64, 0)
		for gene, count := range mousegenemap {
			if first_time {
				GeneLabels = append(GeneLabels, gene)
			}

			row = append(row, count)
		}
		CoexpressionMatrix[index] = row
		first_time = false
		index++
	}

	return CoexpressionMatrix, GeneLabels
}

/*
Input: 2D CoexpressionMatrix - Mice - Rows, Gene- Columns
OUtput: 2D with correlation coefficents of each gene to each other
*/

func BuildCovarianceMatrixPearson(CoexpressionMatrix [][]float64, GeneLabels []string) [][]float64 {
	CovarianceMatrix := make([][]float64, len(GeneLabels))
	for i := range GeneLabels {
		CovarianceMatrix[i] = make([]float64, len(GeneLabels))
		for j := range CovarianceMatrix[i] {
			CovarianceMatrix[i][j] = -1000
		}
		CovarianceMatrix[i][i] = 1
	}
	for i := 0; i < len(CovarianceMatrix); i++ {
		for j := 0; j < len(CovarianceMatrix); j++ {
			if i != j && CovarianceMatrix[i][j] == -1000 {
				matrix1 := make([]float64, 0)
				matrix2 := make([]float64, 0)
				for k := 0; k < 12; k++ {
					matrix1 = append(matrix1, CoexpressionMatrix[k][i])
					matrix2 = append(matrix2, CoexpressionMatrix[k][j])
				}
				x := CalculatePearsonCorrelation(matrix1, matrix2)
				//fmt.Println(i," ",j," ",x)
				CovarianceMatrix[i][j] = x
				// fmt.Println(x)
				CovarianceMatrix[j][i] = x
			}

		}
	}
	return CovarianceMatrix
}

func BuildCovarianceMatrixBiWeight(CoexpressionMatrix [][]float64, GeneLabels []string) [][]float64 {
	CovarianceMatrix := make([][]float64, len(GeneLabels))
	for i := range GeneLabels {
		CovarianceMatrix[i] = make([]float64, len(GeneLabels))
		for j := range CovarianceMatrix[i] {
			CovarianceMatrix[i][j] = -1000
		}
		CovarianceMatrix[i][i] = 1
	}
	for i := 0; i < len(CovarianceMatrix); i++ {
		for j := 0; j < len(CovarianceMatrix); j++ {
			if i != j && CovarianceMatrix[i][j] == -1000 {
				matrix1 := make([]float64, 0)
				matrix2 := make([]float64, 0)
				for k := 0; k < 12; k++ {
					matrix1 = append(matrix1, CoexpressionMatrix[k][i])
					matrix2 = append(matrix2, CoexpressionMatrix[k][j])
				}
				x := BiWeightedMidCorrelation(matrix1, matrix2)
				//fmt.Println(i," ",j," ",x)
				CovarianceMatrix[i][j] = x
				// fmt.Println(x)
				CovarianceMatrix[j][i] = x
			}

		}
	}
	return CovarianceMatrix
}

/*
Input : Signed Weighted Correlation Network
      OR  Unsigned Weighted Correlation Network
Output : AdjacencyMatrix of 1 and 0s with a chosen threshold

*/

func BuildAdjacencyMatrix(a [][]float64) [][]float64 {
	threshold := 0.5
	AdjacencyMatrix := make([][]float64, len(a))
	for i := 0; i < len(a); i++ {
		AdjacencyMatrix[i] = make([]float64, len(a))
	}
	for i := 0; i < len(a); i++ {
		for j := 0; j < len(a); j++ {
			if a[i][j] >= threshold {
				fmt.Println(a[i][j])
				AdjacencyMatrix[i][j] = 1.0
			} else {
				AdjacencyMatrix[i][j] = 0.0
			}
		}
	}
	return AdjacencyMatrix
}

/*
Input: CovarianceMatrix
Output: Transformed continuously with a power adjacency matrix to a weighted network
Selected Optimum value of Beta
*/

func BuildUnsignedWeightedCorrelationNetwork(CovarianceMatrix [][]float64) [][]float64 {
	beta := 6.0
	AdjacencyMatrix := make([][]float64, len(CovarianceMatrix))
	for i := range AdjacencyMatrix {
		AdjacencyMatrix[i] = make([]float64, len(CovarianceMatrix))
		AdjacencyMatrix[i][i] = 0
	}
	for i := range CovarianceMatrix {
		for j := range CovarianceMatrix[0] {
			if i == j {
				AdjacencyMatrix[i][j] = math.Pow(CovarianceMatrix[i][j], beta)
			}

		}
	}
	return AdjacencyMatrix
}

/*
Input:CovarianceMatrix
Output: adjacency matrix dichotomized to arrive at unweighted network
*/

func BuildSignedWeightedCorrelationNetwork(CovarianceMatrix [][]float64) [][]float64 {
	beta := 12.0
	AdjacencyMatrix := make([][]float64, len(CovarianceMatrix))
	for i := range AdjacencyMatrix {
		AdjacencyMatrix[i] = make([]float64, len(CovarianceMatrix))
		AdjacencyMatrix[i][i] = 0
	}
	for i := range CovarianceMatrix {
		for j := range CovarianceMatrix[0] {
			if i != j {
				temp := 0.5 + 0.5*CovarianceMatrix[i][j]
				AdjacencyMatrix[i][j] = math.Pow(temp, beta)
			}

		}
	}
	return AdjacencyMatrix
}

//Calculates the Pearson correlation factor between the expression of two genes
//This is used to construct the similiarity matrix between genes across samples of mice
/*
Input:  two slices to type float64 corresponding to the two genes from samples of mice
Output: Pearson's correlation coefficient between those two genes
*/

func CalculatePearsonCorrelation(x, y []float64) float64 {
	var correlation_coeff float64
	correlation_coeff = 0.0
	var x_bar, y_bar, squaredeviation_x, squaredeviation_y float64
	x_bar = CalculateMean(x)
	y_bar = CalculateMean(y)
	deviation := CalculateDeviation(x, x_bar, y, y_bar)
	squaredeviation_x = CalculateSquareDeviation(x, x_bar)
	squaredeviation_y = CalculateSquareDeviation(y, y_bar)
	correlation_coeff = deviation
	correlation_coeff /= (math.Sqrt(squaredeviation_x) * math.Sqrt(squaredeviation_y))
	return correlation_coeff
}

//Calculates the mean value of a slice of normalized RNA counts
/*
Input: slice of type float64
Output: the average value of the slice
*/

func CalculateMean(x []float64) float64 {
	var sum, mean float64
	sum = 0.0
	mean = 0.0
	for i := range x {
		sum += x[i]
	}
	mean = sum / float64(len(x))
	return mean
}

//Calculate the standard deviation from the average of the gene count from the
/*
Input: slice of type float64
Output: the sum of standard deviation of elements of the slice from the mean
*/

func CalculateDeviation(a []float64, mean_a float64, b []float64, mean_b float64) float64 {
	var sumdeviation float64
	sumdeviation = 0.0
	for i := range a {
		sumdeviation += (a[i] - mean_a) * (b[i] - mean_b)
	}
	return sumdeviation
}

/*
Input: slice of type float64
Output: the sum of squares of standard deviation of elements of the slice from the mean
*/

func CalculateSquareDeviation(a []float64, mean float64) float64 {
	var sumdeviation float64 ///
	sumdeviation = 0.0
	for i := range a {
		sumdeviation += math.Pow(a[i]-mean, 2.0)
	}
	return sumdeviation
}

/*
Input : two slices to type float64 corresponding to the two genes from samples of mice
Output: Biweighted mid correlation coefficient between those two genes
*/
func BiWeightedMidCorrelation(x, y []float64) float64 {
	var med_x, mad_x, med_y, mad_y float64
	//fmt.Println(len(x),len(y))
	med_x = CalculateMedian(x)
	mad_x = CalculateMedianAbsoluteDeviation(x, med_x)
	med_y = CalculateMedian(y)
	mad_y = CalculateMedianAbsoluteDeviation(y, med_y)
	u := make([]float64, 0)
	v := make([]float64, 0)
	for i := range x {
		a := (x[i] - med_x) / (9 * mad_x)
		b := (y[i] - med_y) / (9 * mad_y)
		u = append(u, a)
		v = append(v, b)
	}
	w_x := CalculateWeightsCorrelation(u)
	w_y := CalculateWeightsCorrelation(v)

	x_i_dash := NormalizeWeights(x, w_x, med_x)
	y_i_dash := NormalizeWeights(y, w_y, med_y)
	bicor := 0.0
	for i := 0; i < len(x); i++ {
		bicor = x_i_dash[i] * y_i_dash[i]
	}
	return bicor
}

/*
Input: a slice of type float64 corresponding to one gene and the median value of that count
Output: the median of the standard deviation
*/

func CalculateMedianAbsoluteDeviation(a []float64, median float64) float64 {
	var median_absolute_deviation float64
	a_median := make([]float64, 0)
	for i := range a {
		x := math.Abs(a[i] - median)
		a_median = append(a_median, x)
	}
	median_absolute_deviation = CalculateMedian(a_median)
	return median_absolute_deviation
}

func CalculateMedian(a []float64) float64 {
	QuickSort(a)
	var middleIndex, n int
	n = len(a)
	var median float64
	if len(a)%2 != 0 {
		middleIndex = n / 2
		median = a[middleIndex]
	} else {
		middleIndex = n / 2
		median = (a[middleIndex] + a[middleIndex-1]) / 2
	}
	return median
}

func QuickSort(a []float64) {
	if len(a) < 20 {
		InsertionSort(a)
		return
	}
	pivot := FindPivot(a)
	lowerIndex, higherIndex := Partition(a, pivot)
	QuickSort(a[:lowerIndex])
	QuickSort(a[higherIndex:])
}

func FindPivot(a []float64) float64 {
	n := len(a)
	return Median(a[rand.Intn(n)],
		a[rand.Intn(n)],
		a[rand.Intn(n)])
}

func Median(a1, a2, a3 float64) float64 {
	if a1 < a2 {
		switch {
		case a2 < a3:
			return a2
		case a1 < a3:
			return a3
		default:
			return a1
		}
	}
	switch {
	case a1 < a3:
		return a1
	case a2 < a3:
		return a3
	default:
		return a2
	}
}

// Partition rearranges the elements of a such that:
// - all elements in a[:lowerIndex] are less than the pivot,
// - all elements in a[lowerIndex:higherIndex] are equal to the pivot,
// - all elements in a[higherIndex:] are greater than pivot.
func Partition(a []float64, pivot float64) (lowerIndex, higherIndex int) {
	lowerIndex, higherIndex = 0, len(a)
	for middleIndex := 0; middleIndex < higherIndex; {
		// Invariant:
		//  - a[:lowerIndex] < pivot
		//  - a[lowerIndex:middleIndex] = pivot
		//  - a[middleIndex:higherIndex] are unknown
		//  - a[higherIndex:] > pivot

		switch x := a[middleIndex]; {
		case x < pivot:
			a[middleIndex] = a[lowerIndex]
			a[lowerIndex] = x
			lowerIndex++
			middleIndex++
		case x == pivot:
			middleIndex++
		default: // x > pivot
			a[middleIndex] = a[higherIndex-1]
			a[higherIndex-1] = x
			higherIndex--
		}
	}
	return
}

func InsertionSort(a []float64) {
	for i := 1; i < len(a); i++ {
		// Invariant: a[:j] contains the same elements as
		// the original slice a[:j], but in sorted order.
		key := a[i]
		k := i - 1
		for k >= 0 && a[k] > key {
			a[k+1] = a[k]
			k--
		}
		a[k+1] = key
	}
}

func CalculateWeightsCorrelation(x []float64) []float64 {
	w := make([]float64, 0)
	var I_x float64
	for i := range x {
		a := 1 - (math.Pow(x[i], 2.0))
		if (1 - math.Abs(x[i])) > 0 {
			I_x = 1.0
		} else {
			I_x = 0.0
		}
		weight := math.Pow(a, 2.0) * I_x
		w = append(w, weight)
	}
	return w
}

func NormalizeWeights(a, w []float64, median float64) []float64 {
	normalized_w := make([]float64, 0)
	for i := 0; i < len(a); i++ {
		x := ((a[i] - median) * w[i]) / math.Pow(SumNormalizedWeightSquare(a, w, median), 2.0)
		normalized_w = append(normalized_w, x)
	}
	return normalized_w
}

func SumNormalizedWeightSquare(a, w []float64, median float64) float64 {
	var sum float64
	for i := 0; i < len(a); i++ {
		x := ((a[i] - median) * w[i])
		sum += math.Pow(x, 2.0)
	}
	return sum
}

func CalculateSum(x []float64) float64 {
	var sum float64
	sum = 0.0
	for i := range x {
		sum += x[i]
	}
	return sum
}

func CalculateSumTwo(x, y []float64) float64 {
	var sum float64
	sum = 0.0
	for i := range x {
		sum += x[i] * y[i]

	}
	return sum
}

func CalculateSumSquare(x []float64) float64 {
	var sum float64
	sum = 0.0
	for i := range x {
		sum += math.Pow(x[i], 2.0)
	}
	return sum
}
