from operator import itemgetter
import matplotlib.pyplot as plt
import math
import numpy as np
import os
import csv
import bisect
import random
#import random

import ioLib
import approximateRepeatLib
import graphPlottingLib
import distanceComputeLib
### Checking
def testIndelStatistics():
    listOfData = ioLib.readFromMummerOutput("ecoli536")
    listOfData = ioLib.filterData(listOfData)
    
    listOfData =ioLib.transformMummerOutput(listOfData) 
    listOfData = sorted(listOfData)
    
    listOfData= approximateRepeatLib.findApproximateIndelRepeatStatistics(listOfData, "oneLine.fasta","oneLine2.fasta")
    print listOfData
    

def testreadFromMummerOutput():    
    listOfData = ioLib.readFromMummerOutput("ecoli536")
    
    listOfData = sorted(listOfData, key = itemgetter(2,0,1)) ## sort by order of 2, 0, 1 in the list 
    for eachitem in listOfData:
        print eachitem
    graphPlottingLib.plotExactRepeatStatistics(listOfData)

def checking(temp,genomeSource1,genomeSource2 ):
    print temp
    f1 = open(genomeSource1,'r')
    f2 = open(genomeSource2,'r')
    
    f1.seek(temp[0]-1)
    f2.seek(temp[1]-1)
    
    str1  = f1.read(temp[2])
    str2 = f2.read(temp[2])
    

    print "Hamming distance", distanceComputeLib.hammingDistance(str1,str2, len(str1))


    f1.seek(temp[0]-2)
    f2.seek(temp[1]-2)
    
    str1  = f1.read(temp[2]+2)
    str2 = f2.read(temp[2]+2)
    print "Hamming distance", distanceComputeLib.hammingDistance(str1,str2, len(str1))

    f1.close()
    f2.close()    



def tester(folderName):
    print "Alignment file"    
    fout =open('test.txt', 'w')
    fin = open( folderName+ "/oneline.fasta", 'r')
    
    # input format listOfApproxRepeat : (startIndex1, startIndex2, lapprox, mutationRate, lengthOfExactRepeat)
    listOfRepeat = [[2176433,2292744,5092,7,2862]]
    for eachitem in listOfRepeat :
        fout.write(">Repeat 1 \n")
        fin.seek(eachitem[0])
        temp = fin.read(eachitem[2])
        
        fout.write(temp+ '\n')
        
        fout.write(">Repeat 2 \n")
        fin.seek(eachitem[1])
        temp = fin.read(eachitem[2])
        
        fout.write(temp+ '\n')
        
    
    fin.close()
    fout.close()
    
    
def testCase(folderName):
    
    
    print "Debugging: preprocessingRepeatInfo ---------------"
    dataList = ioLib.preprocessingRepeatInfo(folderName)
    #dataList= ioLib.filterData(dataList)   
    testRepeat = dataList[2]
    # E.g. ['EcoliK12/long_repeats.txt', 3465034, 2724834, 3027]
    
    ### Plot distance between consecutive rise
    print "Debugging: findapproxrepeatLength ---------------"
    startIndex1, startIndex2, lapprox, mutationRate = approximateRepeatLib.findapproxrepeatLength(folderName+ "/oneLine.fasta", folderName+ "/oneLine2.fasta" , testRepeat[1], testRepeat[2], testRepeat[3])
    # E.g. [2724708 3464908 3154 0.000317057704502]
    testApproxRepeat = [startIndex1, startIndex2, lapprox, mutationRate,testRepeat[3], testRepeat[1], testRepeat[2]]
    print "Debugging: plotdetailedPatternForTopRepeats2 ---------------"
    #approximateRepeatLib.plotdetailedPatternForTopRepeats(folderName, 1, [testApproxRepeat],inverted = False)
    approximateRepeatLib.plotdetailedPatternForTopRepeats(folderName, 1, [testApproxRepeat],inverted = False)

    ### Plot sliding hamming window
    print "Debugging: extractApproxRepeatStatistics ---------------"    
    plt.figure(3)    
    graphPlottingLib.plotPatternBeyondRepeat( folderName+ "/oneLine.fasta",folderName+ "/oneLine2.fasta",  startIndex1 - int( 1.5*lapprox), startIndex1 - int( 1.5*lapprox),  startIndex2 - int( 1.5*lapprox) ,  startIndex2 - int( 1.5*lapprox) ,int(lapprox*5/10))
    
    ### Outputing values
    print "Debugging: reportPatternBeyondRepeat ---------------"  
    approximateRepeatLib.reportPatternBeyondRepeat(folderName+ "/oneLine.fasta",folderName+ "/oneLine2.fasta",  startIndex1 - int( 1.5*lapprox), startIndex1 - int( 1.5*lapprox),  startIndex2 - int( 1.5*lapprox) ,  startIndex2 - int( 1.5*lapprox) ,int(lapprox*5/10),folderName+"\outputResult.txt")
      
    ### Edit distance Plot    
    print "Debugging: findapproxrepeatLengthEditDistance ---------------"  
    plt.figure(4)
    approximateRepeatLib.findapproxrepeatLengthEditDistance(folderName+ "/oneLine.fasta",folderName+ "/oneLine2.fasta", testRepeat[1],testRepeat[2], testRepeat[3])
    ### 

    plt.show()	
testCase("EcoliK12")	































