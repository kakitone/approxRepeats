 
from operator import itemgetter
import matplotlib.pyplot as plt
import math
import numpy as np
import os
import csv
import bisect
import random

import distanceComputeLib
import ioLib
import approximateRepeatLib


#import random
### Sliding window View Plot
def readAndPlotStat(file1,file2,folderName, index):
    defaultsize = 28
    
    f1 = open(file1,'r')
    f2 = open(file2,'r')
    
    windowSize = 10
    
    list1 = []
    list2 = []
    plt.figure(figsize=(6*3.13,4*3.13))
    
    temp1 = f1.readline()
    while (len(temp1) > 0):
        list1.append(int(temp1[0:len(temp1)-1]))
        temp1 = f1.readline()
        
    temp2 = f2.readline()
    
    while (len(temp2) > 0):
        list2.append(int(temp2[0:len(temp2)-1]))
        temp2 = f2.readline()    
        
    print list1
    plt.subplot(111)
    plt.plot(range(0,len(list1)*windowSize,windowSize), list1)
    
    

    plt.title(folderName +" right-sided repeat :" +str(index),fontsize=defaultsize)
    #plt.xlabel('Position',fontsize=defaultsize)
    plt.ylabel('Hamming distance \n in a window of 10',fontsize=defaultsize)    
    plt.tick_params(axis='both', which='major', labelsize=defaultsize)
    plt.tick_params(axis='both', which='minor', labelsize=defaultsize)
    
#    plt.subplot(212)
#    plt.plot(range(0,len(list2)*windowSize,windowSize), list2)
        
    f1.close()
    f2.close()

#    plt.title(folderName +" left-sided repeat :" +str(index),fontsize=defaultsize)
#    plt.xlabel('Position',fontsize=defaultsize)
#    plt.ylabel('Hamming distance \n in a window of 10',fontsize=defaultsize)
#    plt.tick_params(axis='both', which='major', labelsize=defaultsize)
#    plt.tick_params(axis='both', which='minor', labelsize=defaultsize)
    

    #plt.show()                                    
    plt.savefig(folderName +"/repeat " +str(index)+ ".png") 
    #plt.show()
    
def plotPatternBeyondRepeat(genomeSource1,genomeSource2, startpt1, endpt1, startpt2, endpt2,plotRange):
    f1 = open(genomeSource1,'r')
    f2 = open(genomeSource2,'r')
    
    pointerLocation1 = endpt1
    pointerLocation2 = endpt2
    
    windowSize = 10
    distanceList = []
    
    for index in range(plotRange):
        f1.seek(pointerLocation1)
        f2.seek(pointerLocation2)
    
        str1 = f1.read(windowSize)
        str2 = f2.read(windowSize)
        
        pointerLocation1 = pointerLocation1 + windowSize
        pointerLocation2 = pointerLocation2 + windowSize
        
        distance = distanceComputeLib.hammingDistance(str1, str2, min(windowSize,len(str1),len(str2)))
        distanceList.append(distance)        
        
    plt.subplot(211)
    plt.plot(range(0,len(distanceList)*windowSize,windowSize), distanceList)

    
    windowSize = 10
    distanceList = []
    
    pointerLocation1 = startpt1 - windowSize
    pointerLocation2 = startpt2 - windowSize


#    for index in range(plotRange):
#        f1.seek(pointerLocation1)
#        f2.seek(pointerLocation2)
    
#        str1 = f1.read(windowSize)
#        str2 = f2.read(windowSize)
        
#        pointerLocation1 = pointerLocation1 - windowSize
#        pointerLocation2 = pointerLocation2 - windowSize
        
#        distance = hammingDistance(str1, str2, windowSize)
#        distanceList.append(distance)        
        
#    plt.subplot(212)
#    plt.plot(range(0,len(distanceList)*windowSize,windowSize), distanceList)
        
    #plt.show()
    
    f1.close()
    f2.close()

def batchPatternBeyondRepeat(sourceFile,genomeSource1, genomeSource2,outputfile):
    listOfData = ioLib.readFromMummerOutput(sourceFile)

    listOfData= sorted(listOfData,key = itemgetter(2), reverse = True)

    if listOfData[0][0] == 1 and listOfData[0][1] == 1:
        listOfData.pop(0)

    start1 = listOfData[0][0]    
    start2 = listOfData[0][1]    
    length = listOfData[0][2]
    
    print listOfData
    print start1
    print start2
    print length
    graphPlottingLib.plotPatternBeyondRepeat(genomeSource1,genomeSource2, start1 + int( length/2) ,start1 + int( length/2) +1 , start2 + int( length/2) ,start2 + int( length/2) +1 , int(length*5/10))
#    reportPatternBeyondRepeat(genomeSource1,genomeSource2, start1 + int( length/2) ,start1 + int( length/2) +1 , start2 + int( length/2) ,start2 + int( length/2) +1 , int(length*5/10),outputfile)

def batchComputationOfStatistics(folders):
    for eachfolder in folders:
        ioLib.generateOneLineFile(eachfolder, eachfolder+ "/genome.fasta")
        batchPatternBeyondRepeat(eachfolder + "/long_repeats.txt", eachfolder + "/oneLine.fasta", eachfolder + "/oneLine2.fasta" , eachfolder+ "/result")
  

### Whole Repeat Spectrum  Plot
def plotExactRepeatStatistics(listOfData):
    countList = []
    for eachitem in listOfData:
        countList.append(eachitem[2])
    countList.pop()
    
    plt.hist(countList,bins = 6000,log=True )
    plt.title("Exact Repeat Statistics")
    plt.xlabel("Length of repeat")
    plt.ylabel("Frequency")
    plt.show()

    
def batchGenerateApproximateRepeatStat(genomeSource1,genomeSource2,mummerOutput,HDrange):
    arrayOfListOfData = []    
    
    listOfData = ioLib.readFromMummerOutput(mummerOutput)
    listOfData = sorted(listOfData,key = itemgetter(2))
    listOfData = ioLib.filterData(listOfData)
    
    arrayOfListOfData.append(listOfData)
    
    
    ### Checking    
    temp = listOfData[len(listOfData)-1]
    #debuggingLib.checking(temp,genomeSource1,genomeSource2 )
    ### End Checking

    for index in range(1,HDrange):
        listOfData = approximateRepeatLib.findApproxRepeatStatistics(listOfData, genomeSource1, genomeSource2)
        listOfData = ioLib.filterData(listOfData)
        arrayOfListOfData.append(listOfData)
        
        temp = listOfData[len(listOfData)-1]
        print listOfData[len(listOfData)-5:len(listOfData)-1]
        #debuggingLib.checking(temp,genomeSource1,genomeSource2 )
        #print listOfData[len(listOfData)-2:len(listOfData)]
        
        
    #plotgraph(arrayOfListOfData, HDrange) 
    
    
def batchGenerateApproximateRepeatStatIndel(genomeSource1,genomeSource2,mummerOutput,HDrange):
    arrayOfListOfData = []    
   
    listOfData = ioLib.readFromMummerOutput(mummerOutput)
    listOfData = sorted(listOfData,key = itemgetter(2))
    listOfData = ioLib.filterData(listOfData)    
    listOfData = ioLib.transformMummerOutput(listOfData)  

    arrayOfListOfData.append(listOfData)
    
    for index in range(1,HDrange):
        listOfData= approximateRepeatLib.findApproximateIndelRepeatStatistics(listOfData, "oneLine.fasta","oneLine2.fasta")
        listOfData = ioLib.filterData(listOfData)
        arrayOfListOfData.append(listOfData)
        temp = listOfData[len(listOfData)-1]
        print "Approx repeat indel"
        print listOfData[len(listOfData)-5:len(listOfData)-1]
        #checking(temp,genomeSource1,genomeSource2 )
        print listOfData[len(listOfData)-2:len(listOfData)]

    plotgraph(arrayOfListOfData, HDrange) 
    
def plotgraph(arrayOfListOfData,HDrange):

    for index in range(0, len(arrayOfListOfData)): 
        dataList = arrayOfListOfData[index]

        if index == 0:        
            ax1 = plt.subplot(HDrange * 100+11)   
            countList = []
            
            for eachitem in dataList:
                if len(eachitem) ==4:
                    countList.append(max(eachitem[2],eachitem[3]))
                else:
                    countList.append(eachitem[2])
                    
            ax1.hist(countList,bins = 6000,log=True)                

        else:
            typeOfPlot = HDrange * 100 +10 + (index+1)
            ax2 = plt.subplot(typeOfPlot,sharex=ax1)      
            countList = []
            
            for eachitem in dataList:
                if len(eachitem) ==4:
                    countList.append(max(eachitem[2],eachitem[3]))
                else:
                    countList.append(eachitem[2])
                    
            ax2.hist(countList,bins = 6000,log=True)            
            
    plt.xlabel("Length of approximate repeat")
    plt.ylabel("Frequency")
    plt.show()

### Classification of approx repeats plot  
def plotClusteringOfRepeat(listOfApproxRepeat,folderName):
    print "clustering"
    xList = []
    yList = []
    lapproxList = []    
    
    ### startIndex1, startIndex2, lapprox, mutationRate, lengthOfExactRepeat
    for eachitem in listOfApproxRepeat:
        xList.append(eachitem[3]*float(eachitem[2]) / (float(max(eachitem[2] - eachitem[4],1))))
        yList.append(float(eachitem[2])/eachitem[4])
        lapproxList.append(eachitem[2])
    
    fig = plt.figure(figsize=(16.0, 10))
    ax = fig.add_subplot(211)
    
    defaultsize = 14
    plt.ylim((0.9,10))
    plt.yscale('log')
    
    plt.xlim((pow(10, -5),1))
    plt.xscale('log')
    
    plt.scatter(xList,yList, c=lapproxList,alpha=0.5, s = 100)
    plt.colorbar()

    plt.axhline(y= 1.25, linewidth=2, color='k', ls= 'dotted')    
        
    
    plt.title("Classification of approximate repeat for " +  folderName,fontsize=defaultsize)
    plt.xlabel("Mutation Rate = number of mutation / ( l_approx - l_exact) ",fontsize=defaultsize)
    plt.ylabel(" l_approx / l_exact ",fontsize=defaultsize)
    
    
def plotApproxRepeatSpectrum(listOfApproxRepeat,longestExactRepeat,folderName):
    
    nonhomologous = []
    homologous = []
    
    longestlapprox = 0 
    for eachitem in listOfApproxRepeat:
        #startIndex1, startIndex2, lapprox, mutationRate, lengthOfExactRepeat = eachitem
        if float(eachitem[2])/eachitem[4]  > 1.25:
            homologous.append(eachitem[2])
        else:
            nonhomologous.append(eachitem[2])
            
        if eachitem[2] > longestlapprox :
            longestlapprox = eachitem[2]
            


    print homologous
    print nonhomologous
    
    plt.subplot(212)
    
    if len(homologous) > 0:    
        homologous = sorted(homologous)
        plt.hist(homologous,bins = range(0,longestlapprox+25,20),color='b')
    if len(nonhomologous) > 0:
        nonhomologous = sorted(nonhomologous)
        plt.hist(nonhomologous, bins  = range(0,longestlapprox+25,20), color = 'r', alpha = 0.5)
    

    plt.axvline(x= longestExactRepeat, linewidth=4, color='g', ls= 'dotted')

    defaultsize = 14
    
    plt.xlim((0, 8500))
    plt.ylim((0,10))
    
    plt.title("Approximate Repeat Spectrum for "+ folderName,fontsize=defaultsize)
    plt.xlabel('Length of approximate repeat',fontsize=defaultsize)
    plt.ylabel('Number of approximate repeat',fontsize=defaultsize)    


### Effect Of Interleaving
def effectOfInterleaving(genomeSource1,Lrepeat, Liid):
        
    f1 = open(genomeSource1, 'r')

    G = len(f1.read())
    print G
    
    oldh = Liid
    
    totalNumberOfRounds = 100000
    
    for numberOfRounds in range(totalNumberOfRounds):
        i = random.randint(0,G- Lrepeat - Liid -1)
        j = random.randint(i+1, G- Lrepeat - Liid)

        f1.seek(i)
        substring1 = f1.read(Liid)
        f1.seek(j)
        substring2= f1.read(Liid)
        
        h1= distanceComputeLib.hammingDistance(substring1,substring2, len(substring1))
        
        f1.seek(i + Lrepeat)
        substring1 = f1.read(Liid)
        f1.seek(j+ Lrepeat)
        substring2= f1.read(Liid)
        
        h2= distanceComputeLib.hammingDistance(substring1,substring2, len(substring1))            

        h = max(h1,h2)
        print i,j, h
        
        if h < oldh:
            oldh = h
                
    print "Minimum hamming distance over " ,totalNumberOfRounds, " is ", oldh, "\n"
    f1.close()
    