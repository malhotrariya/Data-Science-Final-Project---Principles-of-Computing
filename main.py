#15-110 Protein Sequencing Project
#Name: Riya Malhotra
#AndrewID: riyamalh

### Import Libraries ###
import support as test
project = "Protein"

# Opens the file, reads its content, and removes newline characters
def readFile(filename):
    f1 = open(filename, "r")
    text = f1.read()
    text2 = text.replace('\n','')
    return text2

# Converts DNA sequences to RNA sequences and returns a list of codons
def dnaToRna(dna, startIndex):
    codonList = []
    for i in range(startIndex, len(dna), 3):
        codon = dna[i:i+3]
        if (('TAA' == codon) or ('TAG' == codon) or ('TGA' == codon)):
            TreplaceU = codon.replace('T','U')
            codonList.append(TreplaceU)
            break
        TreplaceU = codon.replace('T','U')
        codonList.append(TreplaceU)
    return codonList

# Reads a JSON file containing codon information and creates a dictionary
def makeCodonDictionary(filename):
    import json
    f = open(filename, "r")
    jsonFile = json.load(f)
    aminoAcid = {}
    for amino in jsonFile:
        for codon in jsonFile[amino]:
           TreplaceU = codon.replace('T','U')
           aminoAcid[TreplaceU] = amino
    return aminoAcid

# Generates a protein sequence based on codons and a codon dictionary
def generateProtein(codons, codonD):
    tempList = []
    startFlag = False
    for codon in codons:
        if codon == 'AUG':
            if startFlag == False:
                tempList.append("Start")
                startFlag = True
            else:
               tempList.append("Met")   
        elif (('TAA' == codon) or ('TAG' == codon) or ('TGA' == codon)):
            tempList.append("Stop")
        else:
            tempList.append(codonD[codon])
    return tempList

# Synthesizes proteins from DNA sequences using codon information
def synthesizeProteins(dnaFilename, codonFilename):
    import json
    dnaFile = readFile(dnaFilename)
    codonFile = makeCodonDictionary(codonFilename)
    dnaIndex = 0
    finalList = []
    countVal = 0
    proteinCount = 0
    while(dnaIndex < len(dnaFile)):
        if (len(dnaFile) - dnaIndex) > 3:
            dna = dnaFile[dnaIndex: dnaIndex + 3]
            if dna == 'ATG':
                proteinCount = proteinCount + 1
                codonList = dnaToRna(dnaFile, dnaIndex)
                protein = generateProtein(codonList, codonFile)
                finalList.append(protein)
                dnaIndex = dnaIndex + (3 * len(protein))
            else:
                dnaIndex = dnaIndex + 1
                countVal = countVal + 1
        else:
            break
    print("Total bases:",len(dnaFile))
    print("Unused-base count: ",countVal)
    print("Total synthesized protein:",proteinCount)
    return(finalList)

# Executes Week 1 tasks by synthesizing proteins from human and elephant DNA sequences.
def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


# Finds common proteins between two lists of proteins
def commonProteins(proteinList1, proteinList2):
    commonList = []
    for aminoList1 in proteinList1:
        if aminoList1 in proteinList2 and aminoList1 not in commonList:
            commonList.append(aminoList1)
    return commonList

# Combines a list of proteins into a single list
def combineProteins(proteinList):
    allAminoAcid = []
    for aminoList in proteinList:
        for amino in aminoList:
            allAminoAcid.append(amino)
    return allAminoAcid

# Creates a dictionary mapping amino acids to their occurrences
def aminoAcidDictionary(aaList):
    result = {}
    for amino in aaList:
        if amino not in result:
            result[amino] = 1
        else:
            result[amino] = result[amino] + 1
    return result

# Finds amino acid differences between two protein lists based on frequency
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    result = []
    list1 = combineProteins(proteinList1)
    list2 = combineProteins(proteinList2)
    dict1 = aminoAcidDictionary(list1)
    dict2 = aminoAcidDictionary(list2)
    totalLenList1 = len(list1)
    totalLenList2 = len(list2)
    freq1 = {}
    freq2 = {}
    freqDiff = {}
    for i in dict1:
        freq1[i] = dict1[i] / totalLenList1
    for i in dict2:
        freq2[i] = dict2[i] / totalLenList2
    for i in dict1:
        if i in ["Start","Stop"]:
            continue
        frq1 = freq1[i]
        if i in dict2:
            frq2 = freq2[i]
        else:
            frq2 = 0
        freqDiff[i] = abs(frq1-frq2)
        if abs(frq1-frq2) > cutoff:
            result.append([i,frq1, frq2])
    for i in dict2:
        if i in ["Start","Stop"]:
            continue
        frq2 = freq2[i]
        if i in dict1:
            frq1 = freq1[i]
        else:
            frq1 = 0
        freqDiff[i] = abs(frq1-frq2)
        if abs(frq1-frq2) > cutoff:
            if [i,frq1,frq2] not in result:
                result.append([i,frq1, frq2])
    return result

# Displays common proteins and amino acid differences
def displayTextResults(commonalities, differences):
    print("\nThe following proteins occurred in both DNA Sequences:")
    for protein in commonalities:
        if protein != ["Start", "Stop"]:
            print(protein[1])
    
    print("\nThe following amino acids occurred at very different rates in the two DNA sequences:")
    for amino in differences:
        print(amino[0], " : ", round(amino[1], 2), " in Seq1, ", round(amino[2], 2), " in Seq2" )

# Finding common proteins and differences between human and elephant DNA sequences
def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")
    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)



# Creates a list of amino acid labels based on two protein lists
def makeAminoAcidLabels(proteinList1, proteinList2):
    proteins1 = combineProteins(proteinList1)
    proteins2 = combineProteins(proteinList2)
    amino1 = aminoAcidDictionary(proteins1)
    amino2 = aminoAcidDictionary(proteins2)
    totalAmino = []
    for i in amino1:
        if i not in totalAmino:
            totalAmino.append(i)
    for i in amino2:
        if i not in totalAmino:
            totalAmino.append(i)
    totalAmino.sort()
    return totalAmino

# Sets up data for creating a chart based on labels and a protein list
def setupChartData(labels, proteinList):
    freqList = []
    combineProteinList = combineProteins(proteinList)
    aminoDict = aminoAcidDictionary(combineProteinList)
    totalAminoCount = len(list(combineProteinList))
    for i in labels:
        if i not in aminoDict:
            freqList.append(0)
        else:
            freqList.append(aminoDict[i] / totalAminoCount)
    return freqList

# Creates a chart using matplotlib with specified data
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    import numpy as np
    
    xLabelTuple = tuple(xLabels)
    geneLabel = {}
    geneLabel[label1] = freqList1
    geneLabel[label2] = freqList2
    
    totalFreq = freqList1 + freqList2

    x = np.arange(len(xLabelTuple))  # the label locations
    width = 0.4  # the width of the bars
    multiplier = 0

    fig, ax = plt.subplots(layout='constrained')

    for attribute, measurement in geneLabel.items():
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement, width, label=attribute, edgecolor=edgeList)
        multiplier += 1

    ax.set_xticks(x + width, xLabelTuple)
    ax.legend(loc='upper right', ncols=1)
    ax.set_ylim(0, max(totalFreq) + 0.1)
    fig.tight_layout(pad=3.0)
    
    plt.show()

# Creates a list of edge colors for highlighting differences in the chart
def makeEdgeList(labels, biggestDiffs):
    biggestDiffsAminoList = []
    for i in biggestDiffs:
        if i[0] not in biggestDiffsAminoList:
            biggestDiffsAminoList.append(i[0])
    commonList = []
    for i in labels:
        if i in biggestDiffsAminoList:
            commonList.append("black")
        else:
            commonList.append("white")
    return commonList

# Executes the full program by synthesizing proteins, finding commonalities and differences,
# and creating a chart based on human and elephant DNA sequences.
def runFullProgram():
    humanProteins = synthesizeProteins(r"D:\Work\CMU\Fall 2023 Courses\15-110 Principles of Computing\HW6\hw6_protein_starter\hw6_protein_starter\human_p53.txt",r"D:\Work\CMU\Fall 2023 Courses\15-110 Principles of Computing\HW6\hw6_protein_starter\hw6_protein_starter\codon_table.json")
    elephantProteins = synthesizeProteins(r"D:\Work\CMU\Fall 2023 Courses\15-110 Principles of Computing\HW6\hw6_protein_starter\hw6_protein_starter\elephant_p53.txt",r"D:\Work\CMU\Fall 2023 Courses\15-110 Principles of Computing\HW6\hw6_protein_starter\hw6_protein_starter\codon_table.json")
    print(" Running Full Program!")
    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)
    labels = makeAminoAcidLabels(humanProteins, elephantProteins)
    freqList1 = setupChartData(labels, humanProteins)
    freqList2 = setupChartData(labels, elephantProteins)
    edges = makeEdgeList(labels, differences)
    createChart(labels, freqList1, "Humans", freqList2, "Elephant", edgeList=edges)
    return


# Test Cases
if __name__ == "__main__":
    print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    test.week1Tests()
    print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    runWeek1()

    ## Uncomment these for Week 2 ##
    
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
    

    ## Uncomment these for Week 3 ##
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
