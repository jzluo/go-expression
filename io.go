package main

import (
	"bufio"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
	//"gonum.org/v1/gonum/mat"
)

func io() map[string]map[string]float64 {
	//Input where all the files are
	dirName := "data/"
	files, err := ioutil.ReadDir(dirName)
	if err != nil {
		log.Fatal(err)
		panic("Error: Cannot Read Input Directory.")
	}

	allRawFiles := make([]string, 0)
	// putting all raw file names into a slice of strings
	for _, file := range files {
		if file.Name()[len(file.Name())-3:] == "txt" {
			allRawFiles = append(allRawFiles, dirName+file.Name())
			//fmt.Println("The user has imported:", file.Name())
		}
	}

	threshold := 1.5

	fmt.Print("Importing data... ")
	rawData, genes := ImportRawData(allRawFiles)

	fmt.Print(" done.\n")
	//fmt.Println(rawData)

	fmt.Print("Filtering lowly-expressed genes... ")
	highExpGenes := HighExpGenes(threshold, rawData, genes)
	//fmt.Println(highExpGenes)
	filteredData := FilterGene(rawData, highExpGenes)
	fmt.Print(" done.\n")
	//fmt.Println(filteredData)

	// expressionMatrix, geneLabels:= ConvertMaptoGeneExpressionMatrix(logData)

	// flipedMatrix := FlipMatrix(expressionMatrix)

	numGroups := 3
	numSamplePerGroup :=4
	alphaValue := 0.025
	geneLabels := GeneLabelsInFilteredData(filteredData)
	sortedSampleNames := SortSampleNames(filteredData)
	fmt.Print("Filtering genes with \"flat\" expression... ")
	significantData:= FilterSignificantGenes(alphaValue, numGroups, numSamplePerGroup, sortedSampleNames, geneLabels, filteredData )
	fmt.Print("done.\n")
	fmt.Print("Normalizing... ")
	logData := LogTransform(significantData)
	fmt.Print("done.\n")
	//fmt.Println(logData)
	return logData
}

/*
Input: A slice of string that contain all raw file names
Output: A Map of Map that contains the information about Sample Name(directory), Gene Name and RPKM number
*/
func ImportRawData(allRawFiles []string) (map[string]map[string]float64, []string) {
	// The Key for the map is the Sample Name, the Value for the map is another map of Gene name to RPKM Number
	rawData := make(map[string]map[string]float64)
	allGenes := make([]string, 0)
	//Open one file in the slice
	for i := 0; i < len(allRawFiles); i++ {
		file, err := os.Open(allRawFiles[i])
		if err != nil {
			log.Fatal(err)
			panic("Error: Issue Opening Raw Data Files.")
		}
		defer file.Close()

		scanner := bufio.NewScanner(file)
		// reads the first line and not saving anything
		scanner.Scan()
		// initiate the 2nd map evertime we read a file/a new sample
		rawData[allRawFiles[i]] = make(map[string]float64)
		//read one line at a time starting at the second line
		// A slice of string that holds the information after split them
		splitLine := make([]string, 0)
		for scanner.Scan() {
			//read current line
			currentLine := scanner.Text()
			// split the line into multiple string when seeing a tab
			splitLine = strings.Split(currentLine, "\t")
			geneName := splitLine[0]
			//append all gene names into a slice
			allGenes = append(allGenes, geneName)
			rawData[allRawFiles[i]][geneName], _ = strconv.ParseFloat(splitLine[len(splitLine)-1], 64)
		}
	}

	genes := make([]string, 0)
	genes = RemoveDuplicatesString(allGenes)
	//  for i,_ := range rawData{
	//    fmt.Println(i)
	//  }
	// fmt.Println(rawData)
	return rawData, genes
}

/*
Input: the slice of gene names contain each name numFile times
Output: the slice of gene names contain each name only once
*/
func RemoveDuplicatesString(allGenes []string) []string {
	geneMap := make(map[string]int)
	genes := make([]string, 0)
	for i := range allGenes {
		_, exist := geneMap[allGenes[i]]
		if !exist {
			geneMap[allGenes[i]] = 1
		} else {
			geneMap[allGenes[i]]++
		}
	}
	for key, _ := range geneMap {
		genes = append(genes, key)
	}
	return genes
}

/*
Input: A Map of Map that contains the information about Sample Name(directory), Gene Name and RPKM number and a given threshold for RPKM
Output: A Filtered map of map which contains only information after filtering.
*/
func FilterGene(rawData map[string]map[string]float64, highExpGenes []string) map[string]map[string]float64 {

	for sample, _ := range rawData {
		//fmt.Println(rawData[sample])
		for sample2, _ := range rawData[sample] {
			//fmt.Println(sample2)
			if NeedFilter(sample2, highExpGenes) == true {
				// fmt.Println(rawData[sample][sample2])
				delete(rawData[sample], sample2)
			}
		}
	}
	filteredData := rawData
	return filteredData
}

/*
Input: the secondary map and the list of genes with high expression
Output: a boolean value of whether the gene needs to be filtered or not
*/
func NeedFilter(geneName string, highExpGenes []string) bool {
	for i := range highExpGenes {
		if geneName == highExpGenes[i] {
			return false
		}
	}
	return true
}

/*
Input: Threshold, rawData and genes(a slice of strings which contan gene names)
Output: a slice of string that contain gene names who dosen't need filtering
*/
func HighExpGenes(threshold float64, rawData map[string]map[string]float64, genes []string) []string {
	preserve := make([]string, 0)
	for i := range genes {
		for sample, _ := range rawData {
			if rawData[sample][genes[i]] > threshold {
				//fmt.Println(rawData[sample][genes[i]])
				preserve = append(preserve, genes[i])
			}
		}
	}
	preservedGenes := RemoveDuplicatesString(preserve)
	//fmt.Println(removedGenes)
	//highExpGenes := DeleteLowExpGene(preservedGenes,genes)
	return preservedGenes
}

/*
Input: the map of map of filtered Data
Output: the map of map contains filtered Data with log tranformed RPKM
*/
func LogTransform(filteredData map[string]map[string]float64) map[string]map[string]float64 {
	for _, value := range filteredData {
		for key2, _ := range value {
			value[key2] = math.Log2(value[key2] + 1)
		}
	}
	return filteredData
}

/*
Input: coexpressionMatrix, which is a 2D slice of float64, that has gene as columns and mouse as names
Output: flip the coexpressionMatrix to the flipedMatrix with mouse as coulmns and gene as rows
*/
func FlipMatrix(CoexpressionMatrix [][]float64) [][]float64 {
	flipedMatrix := make([][]float64, len(CoexpressionMatrix[0]))
	for r := range flipedMatrix {
		flipedMatrix[r] = make([]float64, len(CoexpressionMatrix))
	}
	for i := range CoexpressionMatrix[0] {
		for j := range CoexpressionMatrix {
			flipedMatrix[i][j] = CoexpressionMatrix[j][i]
		}
	}
	return flipedMatrix
}

/*
Input: A float of how many groups we have(control,experimental....)(input in func main). We have 3 groups.
Output: A float6 of degree of freedom between the samples
*/
func CalculateDFBetween(numGroups int) float64 {
	return float64(numGroups - 1)
}

/*
Input: Input how many groups(input in func main) and how many sample per group(input in funcm ain)
Output: A integer of degree of freedom within the sample
*/
func CalculateDFWithin(numGroups, numSamplePerGroup int) float64 {
	g := float64(numSamplePerGroup)
	f := float64(numGroups)
	N := (g * f) - f
	return N
}

/*
Input: DF of Between and DF of within
Output: A integer of DF of total
*/
func CalculateDFTotal(numGroups, numSamplePerGroup int) float64 {
	return CalculateDFBetween(numGroups) + CalculateDFWithin(numGroups, numSamplePerGroup)
}

/*
Input: A float of critical value(input in func main)
Output: A string of name of which table is correct to use
*/
func DetermineCriticalTable(alphaValue float64) string {
	if alphaValue == 0.025 {
		return "0.05.txt"
	}
	if alphaValue == 0.05 {
		return "0.05.txt"
	}
	if alphaValue == 0.1 {
		return "0.1.txt"
	}
	if alphaValue == 0.01 {
		return "0.01.txt"
	}
	if alphaValue == 0.001 {
		return "0.001.txt"
	}
	return ("Error: Can Only take 0.001, 0.01, 0.05 and 0.1.")
}

/*
Input: Read in an appropriate file which contains a F Distribution Value
Output: A 2D slice of that table with the first row as DF1(between) and first column as DF2(within)
*/
func ReadCriticalTable(FTableFilesName string) [][]float64 {
	//declare a 2D slice
	FValues := make([][]float64, 36)
	for r := range FValues {
		FValues[r] = make([]float64, 20)
	}

	count := 0
	//Open the File
	file, err := os.Open(FTableFilesName)
	if err != nil {
		log.Fatal(err)
		panic("Error: Issue Opening F Distribution Table Files.")
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	splitLine := make([]string, 0)
	//Read one line at a time
	for scanner.Scan() {
		// which row it is reading
		count++
		//fmt.Println(count-1)
		currentLine := scanner.Text()
		splitLine = strings.Split(currentLine, "\t")
		for i := 0; i < len(splitLine); i++ {
			//readvar,_:=strconv.ParseFloat(splitLine[i], 64)
			FValues[count-1][i], _ = strconv.ParseFloat(splitLine[i], 64)
		}
	}
	return FValues
}

/*
Input: Input DF1 and DF2 from above functions
Output: The appropriate value stored in the F distribution table
*/
func FCritical(numGroups, numSamplePerGroup int, FValues [][]float64) float64 {
	FCriticalValue := 0.0
	dfBetween := CalculateDFBetween(numGroups)
	dfWithin := CalculateDFWithin(numGroups, numSamplePerGroup)

	if dfBetween == 0 {
		panic("Error: dfBetween is zero.")
	}
	if dfWithin == 0 {
		panic("Error: dfWithin is zero.")
	}

	df1InTableIndex := FindClosestDF1inTable(FValues, dfBetween)
	df2InTableIndex := FindClosestDF2inTable(FValues, dfWithin)

	FCriticalValue = FValues[df1InTableIndex][df2InTableIndex]
	return FCriticalValue
}

/*
Input: Input DF1(dfBetween) and the 2D slice contains the F values
Output: The closest DF index in the F values
*/
func FindClosestDF1inTable(FValues [][]float64, dfBetween float64) int {
	//if DF1 is greater than 120, then use tha last column of the table
	df1InTableIndex := 0
	min := math.Inf(+1)
	if dfBetween > 120 {
		//grap the largest df1 value in the table at row 0 and column 19
		df1InTableIndex = 19
		return df1InTableIndex
	}
	// make a slice of floats and contain every df provided by table other than the last one
	df1 := make([]float64, 0)
	//range over the first row of the matrix which contains DF1 values. The first position in the
	for c := 0; c < len(FValues[0])-1; c++ {
		df1 = append(df1, FValues[0][c])
	}
	// range over the slice
	// find the number that has smalles difference to dfBetween
	for i := range df1 {
		if math.Abs(df1[i]-dfBetween) < min {
			min = math.Abs(df1[i] - dfBetween)
			df1InTableIndex = i

		}
	}
	return df1InTableIndex
}

/*
Input: Input DF2(Within) and the 2D slice contains the F values
Output: The index for closes DF exist in the F values
*/
func FindClosestDF2inTable(FValues [][]float64, dfWithin float64) int {
	//if DF1 is greater than 120, then use tha last column of the table
	df2InTableIndex := 0
	min := math.Inf(+1)
	if dfWithin > 120 {
		//grap the largest df1 value in the table at row 0 and column 19
		df2InTableIndex = 34
		return df2InTableIndex
	}
	// make a slice of floats and contain every df provided by table other than the last one
	df2 := make([]float64, 0)
	//range over the first row of the matrix which contains DF1 values. The first position in the
	for c := 0; c < len(FValues)-1; c++ {
		df2 = append(df2, FValues[c][0])
	}
	// range over the slice
	// find the number that has smalles difference to dfBetween
	for i := range df2 {
		if math.Abs(df2[i]-dfWithin) < min {
			min = math.Abs(df2[i] - dfWithin)
			df2InTableIndex = i
		}
	}
	return df2InTableIndex
}

/*
Input: The map of Map contain filtered Data
Output: A slice of string contains sorted key for the map(sorted sample name)
Because in filteredData before this function, the sample(mouse) names are random.
We need to sort them in a correct order so every group's samples are adjacent.
This will be convinient to caculate the mean for one gene in each group for Anova Testing.
*/
func SortSampleNames(filteredData map[string]map[string]float64) []string {
	sortedSampleNames := make([]string, 0)
	for key := range filteredData {
		sortedSampleNames = append(sortedSampleNames, key)
	}
	return QuickSortStrings(sortedSampleNames)
}

/*
Input: A slice of string
Output: A slice of sorted string
*/
func QuickSortStrings(strings []string) []string {
	if len(strings) < 2 { // base case
		return strings
	}

	jon := strings[0] // set first element as "pivot"
	i := 1            // start from 2nd element
	j := len(strings) - 1
	for j-i > 0 { // stop when i, j meet
		for IsGreater(strings[i], jon) != true && strings[i] != jon { // if strings equal, swap anyway
			if i < j {
				i++
			} else {
				break
			}
		}

		for IsGreater(strings[j], jon) != false && strings[j] != jon {
			if j > i {
				j--
			} else {
				break
			}
		}
		strings[i], strings[j] = strings[j], strings[i] // swap the two and continue
	}
	strings[0], strings[i] = strings[i], strings[0] // move pivot to end of left section

	QuickSortStrings(strings[:i])
	QuickSortStrings(strings[i+1:])

	return strings
}

/*
Input: 2 strings
Output: return true if s1 has greater length than s2
*/
func IsGreater(s1, s2 string) bool {
	// Compare if s1 should come after s2 lexicographically
	s1 = strings.ToLower(s1)
	s2 = strings.ToLower(s2)

	for i := 0; i < len(s1) && i < len(s2); i++ {
		if s1[i] < s2[i] {
			return false
		} else if s2[i] < s1[i] {
			return true
		}
	}
	if len(s1) > len(s2) { // if identical up to end of the shorter string, longer > shorter
		return true
	} else {
		return false
	}
}

/*
Input: the filteredData
Output: a slice of string contains all gene names
*/
func GeneLabelsInFilteredData(filteredData map[string]map[string]float64) []string {
	geneLabels1 := make([]string, 0)
	for sample, _ := range filteredData {
		for sample2, _ := range filteredData[sample] {
			geneLabels1 = append(geneLabels1, sample2)
		}
	}
	geneLabels := RemoveDuplicatesString(geneLabels1)
	return geneLabels
}

/*
Input: filtereData, a map of map contains all data after low threshold filtering.
      sortedSample names and geneLabels(a list of genes exist in the map)
      also input numGroups and numSampple per group
Output: significantData, a map of map contains all genes that only shows
        a significant differences between control and experimental group.
*/

func FilterSignificantGenes(alphaValue float64, numGroups, numSamplePerGroup int, sortedSampleNames, geneLabels []string, filteredData map[string]map[string]float64) map[string]map[string]float64 {
	expValues := make([]float64, 0)
	meanValues := make([]float64, 0)
	count := 0

	for sample, _ := range filteredData {
		for sample2, _ := range filteredData[sample] {
			geneLabels = append(geneLabels, sample2)
		}
	}
	//loop through the slice of gene names
	for i := range geneLabels {
		// loop through the sorted sample names
		for j := range sortedSampleNames {
			//for that one gene and that one mouse, we append the expression into a slice
			expValues = append(expValues, filteredData[sortedSampleNames[j]][geneLabels[i]])
		}
		meanValues = MeanExpForOneGene(numGroups, numSamplePerGroup, expValues)
		if IsSignificant(alphaValue, numGroups, numSamplePerGroup, expValues, meanValues) == false {
			count++
			for sample, _ := range filteredData {
				delete(filteredData[sample], geneLabels[i])

			}
		}
	}
	return filteredData
}

func IsSignificant(alphaValue float64, numGroups, numSamplePerGroup int, expValues, meanValues []float64) bool {
	FTableName := DetermineCriticalTable(alphaValue)
	FValues := ReadCriticalTable(FTableName)
	FCriticalValue := FCritical(numGroups, numSamplePerGroup, FValues)
	dataFValues := CalculateDataFValue(numGroups, numSamplePerGroup, expValues, meanValues)
	if dataFValues < FCriticalValue {
		return false
	}
	return true
}

func CalculateSumOfSquaresWithin(numGroups, numSamplePerGroup int, expValues []float64) float64 {
	a := 0
	b := numSamplePerGroup
	sum := 0.0
	ssValues := make([]float64, 0)

	for k := 0; k < numGroups; k++ {
		l := numSamplePerGroup
		slice := expValues[a:b]
		mean := CalculateMean(expValues[a:b])
		ss := CalculateSumOfSquares(slice, mean)
		a = a + l
		b = a + l
		ssValues = append(ssValues, ss)
	}
	for i := range ssValues {
		sum = sum + ssValues[i]
	}
	return sum
}

func CalculateSumOfSquaresBetween(numGroups, numSamplePerGroup int, expValues, meanValues []float64) float64 {
	grandMean := CalculateGrandMean(meanValues)
	sumOfSquaresTotal := CalculateSumOfSquares(expValues, grandMean)
	sumOfSquaresWithin := CalculateSumOfSquaresWithin(numGroups, numSamplePerGroup, expValues)
	sumOfSquaresBetween := sumOfSquaresTotal - sumOfSquaresWithin

	return sumOfSquaresBetween
}

/*
Input: the mean of one gene expression of each group
Output: the sum of squares
*/
func CalculateSumOfSquares(expValues []float64, mean float64) float64 {
	sumOfSquares := 0.0
	for i := range expValues {
		sumOfSquares = sumOfSquares + (expValues[i]-mean)*(expValues[i]-mean)
	}
	return sumOfSquares
}

/*
Input: also input numGroups and numSampple per group and slice of gene expression values
Output: a slice contains only mean gene values for each group[meanValues]
*/
func MeanExpForOneGene(numGroups, numSamplePerGroup int, expValues []float64) []float64 {
	a := 0
	b := numSamplePerGroup
	meanValues := make([]float64, 0)

	for k := 0; k < numGroups; k++ {
		l := numSamplePerGroup
		mean := CalculateMean(expValues[a:b])
		a = a + l
		b = a + l
		meanValues = append(meanValues, mean)
	}
	if len(meanValues) == numGroups {
		return meanValues
	}
	panic("Error:meanExpForOneGene return more values than Number of Groups")
}

/*
Input: a slice of floats contain mean values
Output: one float represent the grand mean
*/
func CalculateGrandMean(meanValues []float64) float64 {
	grandMean := 0.0
	sum := 0.0
	for i := range meanValues {
		sum = sum + meanValues[i]
	}
	grandMean = sum / float64(len(meanValues))
	return grandMean
}

func CalculateVarianceBet(numGroups, numSamplePerGroup int, expValues, meanValues []float64) float64 {
	ssBetween := CalculateSumOfSquaresBetween(numGroups, numSamplePerGroup, expValues, meanValues)
	dfBetween := CalculateDFBetween(numGroups)
	varBetween := ssBetween / dfBetween
	return varBetween
}

func CalculateVarianceWit(numGroups, numSamplePerGroup int, expValues []float64) float64 {
	ssWithin := CalculateSumOfSquaresWithin(numGroups, numSamplePerGroup, expValues)
	dfWithin := CalculateDFWithin(numGroups, numSamplePerGroup)
	varWithin := ssWithin / dfWithin
	return varWithin
}

func CalculateDataFValue(numGroups, numSamplePerGroup int, expValues, meanValues []float64) float64 {
	varBetween := CalculateVarianceBet(numGroups, numSamplePerGroup, expValues, meanValues)
	varWithin := CalculateVarianceWit(numGroups, numSamplePerGroup, expValues)
	dataFValues := varBetween / varWithin
	return dataFValues
}
