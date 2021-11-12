"""
Protein Sequencing Project
Name:
Roll Number:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    f=open(filename)
    t=f.read()
    text=''
    for i in t:
        if i!="\n":
            text+=i
    return text


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    stop=["UAA","UAG","UGA"]
    j=startIndex
    codons=[]
    while j in range(len(dna)):
        codon=""
        for i in range(j,j+3):
            if dna[i]=="T":
                codon+="U"
            else:
                codon+=dna[i]
        codons.append(codon)
        if codon in stop:
            return codons
        j+=3
    return codons


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    f=open(filename)
    j=json.load(f)
    d={}
    for amino in j:
        for k in j[amino]:
            codon=""
            for i in range(len(k)):
                if k[i]=="T":
                    codon+="U"
                else:
                    codon+=k[i]
            d[codon]=amino    
    return d


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    protein=[]
    for i in range(len(codons)):
        codon=codons[i]
        if i==0:
            protein.append("Start")
        else:
            protein.append(codonD[codon])
        if codonD[codon]=="Stop":
            return protein
    return protein


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    dna=readFile(dnaFilename)
    codonD=makeCodonDictionary(codonFilename)
    proteins=[]
    count=0
    i=0
    while i in range(len(dna)):
        if dna[i:i+3]=="ATG":
            codons=dnaToRna(dna,i)
            protein=generateProtein(codons,codonD)
            proteins.append(protein)
            i+=3*len(protein)
        else:
            i+=1
            count+=1
    print(len(proteins),count)
    return proteins


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    uniqueproteins=[]
    for i in proteinList1:
        if i not in uniqueproteins:
            uniqueproteins.append(i)
    commonproteins=[]
    for i in uniqueproteins:
        if i in proteinList2:
            if i not in commonproteins:
                commonproteins.append(i)
    return commonproteins


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    amino=[]
    for i in proteinList:
        for j in i:
            amino.append(j)
    return amino



'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    d={}
    for i in aaList:
        if i not in d:
            d[i]=0
        d[i]+=1
    return d

'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    list1=combineProteins(proteinList1)
    list2=combineProteins(proteinList2)
    d1=aminoAcidDictionary(list1)
    d2=aminoAcidDictionary(list2)
    d=[]
    aminoacids=[]
    for i in d1:
        d1[i]=d1[i]/len(list1)
        if i not in d:
            d.append(i)
    for i in d2:
        d2[i]=d2[i]/len(list2)
        if i not in d:
            d.append(i)
    for i in d:
        if i!="Start" and i!="Stop":
            if i not in d1:
                f1=0
            else:
                f1=d1[i]
            if i not in d2:
                f2=0
            else:
                f2=d2[i]
            if abs(f1-f2)>cutoff:
                lis=[i,f1,f2]
                aminoacids.append(lis)
    return aminoacids



'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("The following proteins occurred in both DNA Sequences:")
    for i in commonalities:
        for j in i:
            if j!="Start" and j!="Stop":
                print(j)
    print("The following amino acids occurred at very different rates in the two DNA sequences:")
    for i in differences:
        print(i[0],":",round(i[1]*100,2),"% in seq1",round(i[2]*100,2),"% in seq2")



def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    p1=combineProteins(proteinList1)
    p2=combineProteins(proteinList2)
    p=[]
    for i in p2:
        if i not in p:
            p.append(i)
    for i in p1:
        if i not in p:
            p.append(i)
    return sorted(p)


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    lis=combineProteins(proteinList)
    d=aminoAcidDictionary(lis)
    freq=[]
    for i in labels:
        if i not in d:
            freq.append(0)
        else:
            freq.append(d[i]/len(lis))
    return freq


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    w=0.35
    plt.bar(xLabels,freqList1,width=-w,align='edge',edgecolor=edgeList,label=label1)
    plt.bar(xLabels,freqList2,width=w,align='edge',edgecolor=edgeList,label=label2)
    plt.xticks(rotation="vertical")
    plt.legend()
    plt.show()


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    amino=[]
    edge=[]
    for i in biggestDiffs:
        amino.append(i[0])
    for i in labels:
        if i in amino:
            edge.append("black")
        else:
            edge.append("white")
    return edge


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    test.week1Tests()
    print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    runWeek1()
    # test.testReadFile()
    # test.testDnaToRna()
    # test.testMakeCodonDictionary()
    # test.testGenerateProtein()
    # test.testSynthesizeProteins()

    ## Uncomment these for Week 2 ##
    
    # print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    # test.week2Tests()
    # print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    # runWeek2()
    
    # test.testCommonProteins()
    # test.testCombineProteins()
    # test.testAminoAcidDictionary()
    # test.testFindAminoAcidDifferences()

    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """
    # test.testMakeAminoAcidLabels()
    # test.testSetupChartData()
    # test.testCreateChart()
    test.testMakeEdgeList()
