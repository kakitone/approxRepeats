from operator import itemgetter
import matplotlib.pyplot as plt
import math
import numpy as np
import os
import csv
import bisect
import random
#import random

### Subroutine
def hammingDistance(str1, str2, k):
    count  =0
    for index in range(k):
        if str1[index] != str2[index]:
            count = count +1
    return count

def filterData(listOfData):
    listOfData = sorted(listOfData)
    oldData = []
    index =0
    while (index < len(listOfData)):
        if listOfData[index][0] >= listOfData[index][1]: 
            listOfData.pop(index)
        elif listOfData[index] == oldData:
            listOfData.pop(index)
        else:
            oldData = listOfData[index]
            index = index +1
    #print listOfData
    listOfData = sorted(listOfData,key = itemgetter(2))
    return listOfData
        
def readFromMummerOutput(sourceFile):
    f = open(sourceFile, 'r')
    genomeName = f.readline()
    print genomeName
    testline = 'NOT NULL YET'
    listOfData = []
    while len(testline) > 0 :
        testline = f.readline()
        testList = testline.split()
        testList = map(int, testList)
        if len(testList) > 0 :
            listOfData.append(testList)
    f.close()
    return listOfData           
            
### Approximate repeat statistics (Hamming type)
def findApproxRepeatStatistics(listOfData, genomeSource1, genomeSource2):    
    listOfData = sorted(listOfData)
    approximateRepeatList = []
    sourceFile = genomeSource1
    sourceFile2 = genomeSource2
    
    f = open(sourceFile, 'r') 
    f2 = open(sourceFile2, 'r')
        
    for eachitem in listOfData:

        start1  = eachitem[0] -1
        start2  = eachitem[1] -1
        end1 = eachitem[0] + eachitem[2] -2        
        end2 = eachitem[1] + eachitem[2] -2
        
        check  = 1 
        counter = 2
        loc1 =  end1 +counter
        loc2 =  end2 +counter
     
        while (check and loc1 < len(sourceFile) and loc2 < len(sourceFile2) ) : 
            f.seek(loc1)
            f2.seek(loc2)
            t1 = f.read(1)
            t2 = f2.read(1)
            
            if  t1!=t2 :
                check  =0 
            else: 
                counter = counter + 1
                loc1 =end1 +counter
                loc2 = end2 +counter

                
                
        newentry = [start1 + 1, start2 + 1, eachitem[2] + counter -1]     
        approximateRepeatList.append(newentry)
        
        check  = 1 
        counter = 2
        loc1 =  start1 -counter
        loc2 =  start2 -counter
        

        while (check and loc1>=0  and loc2 >=0 ) : 
            f.seek(loc1)
            f2.seek(loc2)
            t1 = f.read(1)
            t2 = f2.read(1)

            if  t1!=t2 :
                check  =0 
            else: 
                counter = counter + 1
                loc1 = start1 -counter
                loc2 =  start2 -counter    
        
        newentry = [loc1 + 2, loc2+2, eachitem[2] + counter -1]      
        approximateRepeatList.append(newentry)

    print "Approx repeat "    
    approximateRepeatList = sorted(approximateRepeatList,key = itemgetter(2))

    f.close()
    f2.close()
    
    return approximateRepeatList

### Approximate repeat statistics (edit distance type)
def transformMummerOutput(listOfData):
    newListOfData = []
    
    for eachitem in listOfData :
        temp = eachitem
        temp.append(eachitem[2])
        newListOfData.append(temp)
    return newListOfData

def findApproximateIndelRepeatStatistics(listOfData, genomeSource1, genomeSource2):    
    
    listOfData = sorted(listOfData)
    approximateRepeatList = []
    sourceFile = genomeSource1
    sourceFile2 = genomeSource2
    
    f = open(sourceFile, 'r') 
    f2 = open(sourceFile2, 'r')
    
    for eachitem in listOfData:
        start1  = eachitem[0] -1
        start2  = eachitem[1] -1
        end1 = eachitem[0] + eachitem[2] -2        
        end2 = eachitem[1] + eachitem[2] -2
        
        
        ### Modify part 1
        check  = 1 
        counter = 2
        loc1 =  end1 +counter -1
        loc2 =  end2 +counter
     
        while (check and loc1 < len(sourceFile) and loc2 < len(sourceFile2) ) : 
            f.seek(loc1)
            f2.seek(loc2)
            t1 = f.read(1)
            t2 = f2.read(1)
            
            if  t1!=t2 :
                check  =0 
            else: 
                counter = counter + 1
                loc1 =end1 +counter -1
                loc2 = end2 +counter
    
                
                
        newentry = [start1 + 1, start2 + 1, eachitem[2] + counter -2,  eachitem[2] + counter -1]     
        approximateRepeatList.append(newentry)
        
        
        ### Modify part 2
        check  = 1 
        counter = 2
        loc1 =  end1 +counter 
        loc2 =  end2 +counter -1
     
        while (check and loc1 < len(sourceFile) and loc2 < len(sourceFile2) ) : 
            f.seek(loc1)
            f2.seek(loc2)
            t1 = f.read(1)
            t2 = f2.read(1)
            
            if  t1!=t2 :
                check  =0 
            else: 
                counter = counter + 1
                loc1 =end1 +counter 
                loc2 = end2 +counter -1
    
                
                
        newentry = [start1 + 1, start2 + 1, eachitem[2] + counter -1,  eachitem[2] + counter -2]     
        approximateRepeatList.append(newentry)
        
        ### Modify part 3 
        check  = 1 
        counter = 2
        loc1 =  start1 -counter  +1
        loc2 =  start2 -counter
        
    
        while (check and loc1>=0  and loc2 >=0 ) : 
            f.seek(loc1)
            f2.seek(loc2)
            t1 = f.read(1)
            t2 = f2.read(1)
    
            if  t1!=t2 :
                check  =0 
            else: 
                counter = counter + 1
                loc1 = start1 -counter +1
                loc2 =  start2 -counter    
        
        newentry = [loc1 + 2, loc2+2, eachitem[2] + counter ,eachitem[2] + counter -1]      
        approximateRepeatList.append(newentry)
        
        ### Modify part 4
        check  = 1 
        counter = 2
        loc1 =  start1 -counter  
        loc2 =  start2 -counter +1
        
    
        while (check and loc1>=0  and loc2 >=0 ) : 
            f.seek(loc1)
            f2.seek(loc2)
            t1 = f.read(1)
            t2 = f2.read(1)
    
            if  t1!=t2 :
                check  =0 
            else: 
                counter = counter + 1
                loc1 = start1 -counter 
                loc2 =  start2 -counter +1    
        
        newentry = [loc1 + 2, loc2+2, eachitem[2] + counter-1 ,eachitem[2] + counter ]      
        approximateRepeatList.append(newentry)
    
        
        
    f.close()
    f2.close()
    return approximateRepeatList
    
### Checking
def testIndelStatistics():
    listOfData = readFromMummerOutput("ecoli536")
    listOfData = filterData(listOfData)
    
    listOfData =transformMummerOutput(listOfData) 
    listOfData = sorted(listOfData)
    
    listOfData= findApproximateIndelRepeatStatistics(listOfData, "oneLine.fasta","oneLine2.fasta")
    print listOfData
    

def testreadFromMummerOutput():    
    listOfData = readFromMummerOutput("ecoli536")
    
    listOfData = sorted(listOfData, key = itemgetter(2,0,1)) ## sort by order of 2, 0, 1 in the list 
    for eachitem in listOfData:
        print eachitem
    plotExactRepeatStatistics(listOfData)

def checking(temp,genomeSource1,genomeSource2 ):
    print temp
    f1 = open(genomeSource1,'r')
    f2 = open(genomeSource2,'r')
    
    f1.seek(temp[0]-1)
    f2.seek(temp[1]-1)
    
    str1  = f1.read(temp[2])
    str2 = f2.read(temp[2])
    

    print "Hamming distance", hammingDistance(str1,str2, len(str1))


    f1.seek(temp[0]-2)
    f2.seek(temp[1]-2)
    
    str1  = f1.read(temp[2]+2)
    str2 = f2.read(temp[2]+2)
    print "Hamming distance", hammingDistance(str1,str2, len(str1))

    f1.close()
    f2.close()    

### Graph plotting   

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
    
    listOfData = readFromMummerOutput(mummerOutput)
    listOfData = sorted(listOfData,key = itemgetter(2))
    listOfData = filterData(listOfData)
    
    arrayOfListOfData.append(listOfData)
    
    
    ### Checking    
    temp = listOfData[len(listOfData)-1]
    checking(temp,genomeSource1,genomeSource2 )
    ### End Checking

    for index in range(1,HDrange):
        listOfData = findApproxRepeatStatistics(listOfData, genomeSource1, genomeSource2)
        listOfData = filterData(listOfData)
        arrayOfListOfData.append(listOfData)
        
        temp = listOfData[len(listOfData)-1]
        print listOfData[len(listOfData)-5:len(listOfData)-1]
        checking(temp,genomeSource1,genomeSource2 )
        #print listOfData[len(listOfData)-2:len(listOfData)]
        
        
    #plotgraph(arrayOfListOfData, HDrange) 
    
    
def batchGenerateApproximateRepeatStatIndel(genomeSource1,genomeSource2,mummerOutput,HDrange):
    arrayOfListOfData = []    
   
    listOfData = readFromMummerOutput(mummerOutput)
    listOfData = sorted(listOfData,key = itemgetter(2))
    listOfData = filterData(listOfData)    
    listOfData =transformMummerOutput(listOfData)  

    arrayOfListOfData.append(listOfData)
    
    for index in range(1,HDrange):
        listOfData= findApproximateIndelRepeatStatistics(listOfData, "oneLine.fasta","oneLine2.fasta")
        listOfData = filterData(listOfData)
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
        
        distance = hammingDistance(str1, str2, min(windowSize,len(str1),len(str2)))
        distanceList.append(distance)        
        
    plt.subplot(111)
    plt.plot(range(0,len(distanceList)*windowSize,windowSize), distanceList)

    
    windowSize = 10
    defaultsize = 24
    
    lims = plt.xlim()    
    
    plt.tick_params(axis='both', which='major', labelsize=defaultsize)
    plt.tick_params(axis='both', which='minor', labelsize=defaultsize)
    
       
    #plt.title("Analysis of region beyond repeat in ")
    plt.ylabel("Distance within a window of length 10 ",fontsize=defaultsize)
    plt.xlabel("Genomics location ( unit : base )",fontsize=defaultsize)

    
    plt.show()

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
    listOfData = readFromMummerOutput(sourceFile)

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
    plotPatternBeyondRepeat(genomeSource1,genomeSource2, start1 + int( length/2) ,start1 + int( length/2) +1 , start2 + int( length/2) ,start2 + int( length/2) +1 , int(length*5/10))
#    reportPatternBeyondRepeat(genomeSource1,genomeSource2, start1 + int( length/2) ,start1 + int( length/2) +1 , start2 + int( length/2) ,start2 + int( length/2) +1 , int(length*5/10),outputfile)
   
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
        
        h1= hammingDistance(substring1,substring2, len(substring1))
        
        f1.seek(i + Lrepeat)
        substring1 = f1.read(Liid)
        f1.seek(j+ Lrepeat)
        substring2= f1.read(Liid)
        
        h2= hammingDistance(substring1,substring2, len(substring1))            

        h = max(h1,h2)
        print i,j, h
        
        if h < oldh:
            oldh = h
                
    print "Minimum hamming distance over " ,totalNumberOfRounds, " is ", oldh, "\n"
    f1.close()
    
### Batch generation of statistics
def detailRepeatPattern(genomeSource1,genomeSource2, startpt1, startpt2, lrepeat,mutationTolerant,outputResult):
    
    print startpt1, startpt2, lrepeat
    f1 = open(genomeSource1,'r')
    f2 = open(genomeSource2,'r')
    
    lastTime = 'N'
    l1 =     lrepeat
    l2 =     lrepeat 
    deltal1 = 0
    deltal2  =0 
    
    outputList = []	
    runningIndex  =0 
    outputList.append([runningIndex,l1+l2 - lrepeat])
    
    while runningIndex < mutationTolerant:
        #print deltal1, deltal2
        if lastTime == 'L' or lastTime == 'N':
            # extend L
            runningForwarder1 = startpt1+lrepeat-l1 - 2
            runningForwarder2 = startpt2+lrepeat-l1  -2
            f1.seek(runningForwarder1)
            f2.seek(runningForwarder2)
            # Skip the mutated base
            temp = f1.read(1)
            temp2 = f2.read(1)    
            # Check 
            #print "Should be false", temp == temp2, temp,temp2
            
            temp = f1.read(1)
            temp2 = f2.read(1)    
            # Check 2            
            #print "Should be True", temp == temp2, temp,temp2
            # continue searching
            runningForwarder1 = runningForwarder1 -1
            runningForwarder2 = runningForwarder2 -1
                        
            f1.seek(runningForwarder1)
            f2.seek(runningForwarder2)            
            
            temp = f1.read(1)
            temp2 = f2.read(1)
            
            deltal1 = 0
            while (temp == temp2):
                deltal1 = deltal1 + 1
                
                runningForwarder1 = runningForwarder1 - 1
                runningForwarder2 = runningForwarder2 -1
                #print temp == temp2, runningForwarder1,runningForwarder2
                f1.seek(runningForwarder1)
                f2.seek(runningForwarder2) 
            
                temp = f1.read(1)
                temp2 = f2.read(1)
            
            
            deltal1 = deltal1 + 1
            #print deltal1
			
        if  lastTime == 'R'or lastTime == 'N':
            # extend R
           # print  f1.read(1)
            
            f1.seek(startpt1 -1 + l2)
            f2.seek(startpt2 -1  + l2)
            # Skip the mutated base
            temp = f1.read(1)
            temp2 = f2.read(1)
            # Check 
            #print "Should be false", temp == temp2, temp,temp2
            
           # print temp, temp2
            # continue searching
            temp = f1.read(1)
            temp2 = f2.read(1)
            # Check 2            
            #print "Should be True", temp == temp2, temp,temp2
                        
            deltal2 = 0
            while (temp == temp2):
                deltal2 = deltal2 + 1
                temp = f1.read(1)
                temp2 = f2.read(1)
            
            deltal2 = deltal2 + 1
           # print deltal2
        runningIndex = runningIndex +1           
        
        if deltal1 > deltal2 : 
            lastTime = 'L'
            l1 = l1 + deltal1
            outputList.append([runningIndex,l1+l2 - lrepeat])
            
        elif deltal1 < deltal2:
            lastTime = 'R'
            l2 = l2 + deltal2
            outputList.append([runningIndex,l1+l2 - lrepeat])
            
        elif deltal1 == deltal2: 
            breakTie = random.random()
            if breakTie < 0.5:
                lastTime = 'L'
                l1 = l1 + deltal1
                outputList.append([runningIndex,l1+l2 - lrepeat])
            else:
                lastTime = 'R'
                l2 = l2 + deltal2
                outputList.append([runningIndex,l1+l2 - lrepeat])
   
    f1.close()
    f2.close()
    return outputList
    
def plotdetailedPatternForTopRepeats(folderName, topX, dataList):
    #generateOneLineFile(folderName, folderName+ "/genome.fasta")  
    #dataList = preprocessingRepeatInfo(folderName)
    #print dataList
    ### Input format : [startIndex1, startIndex2, lapprox, mutationRate,lengthOfExactRepeat,start1,start2])
    print dataList    
    lengthOfRepeatList = []
    for plotindex  in range(2):
        plt.figure(figsize=(16.0, 10))
        plt.grid(which = 'both')
        if plotindex  == 0:
            plt.xlim(0,0.45)
            plt.ylim(0.9,8)
        
            
        for eachitem in dataList:
            outputList =     detailRepeatPattern(folderName+"/oneLine.fasta",folderName+"/oneLine2.fasta", eachitem[5], eachitem[6], eachitem[4],50,folderName+"/outputResult.txt")
            lrepeat = eachitem[4]
            
            #print "outputList"
            plotY = []
            plotX = []
            for inneritem in outputList:
                #print eachitem
                plotY.append(inneritem[1]/float(lrepeat))
                plotX.append(inneritem[0]/float(lrepeat))
            
            plotY2 = []
            for elementindex in range(len(plotY)-1):
                plotY2.append( (plotY[elementindex+1] - plotY[elementindex] )/ float(plotX[elementindex+1] - plotX[elementindex]))
            
            #print "plotY2", plotY2
            #print "plotY", plotY
            #print "plotX", plotX
            
            #print "mean(plotY2)",np.mean(plotY2)
            #plt.yscale('log')
            #plt.plot(plotX[0:-1], plotY2,linewidth=1, marker = 'x')
            plt.plot(plotX, plotY,linewidth=1, marker = 'x')
            lengthOfRepeatList.append(str(eachitem[2])+" , "+ str(eachitem[4]))
        defaultsize = 24
        
        lims = plt.xlim()    
    
        plt.tick_params(axis='both', which='major', labelsize=defaultsize)
        plt.tick_params(axis='both', which='minor', labelsize=defaultsize)
    
       
        plt.title("Analysis of region beyond repeat in " +  folderName +" for top"+str(topX) +" Exact Repeats",fontsize=defaultsize)
        plt.xlabel("Mutation Rate = number of mutation / l_exactRepeat ",fontsize=defaultsize)
        plt.ylabel("Stretch =  l_approxRepeat / l_exactRepeat ",fontsize=defaultsize)
        plt.legend(lengthOfRepeatList)
 
        if plotindex  == 0:       
            plt.savefig("regionBeyondRepeat/FixScale/"+folderName+"_approxRepeatAnalysisPlot.png") 
            
             
            
        
        else:
            plt.savefig("regionBeyondRepeat/bestFit/"+folderName+"_approxRepeatAnalysisPlot.png") 
           # plt.show()
        

    

def reportPatternBeyondRepeat(genomeSource1,genomeSource2, startpt1, endpt1, startpt2, endpt2,plotRange,outputResult):
    Gchecker = open(genomeSource1,'r')
    G = len(Gchecker.read())
    print G
    Gchecker.close()    
    
    f1 = open(genomeSource1,'r')
    f2 = open(genomeSource2,'r')
    

    pointerLocation1 = endpt1
    pointerLocation2 = endpt2
    
    windowSize = 10
    distanceList = []
    
    for index in range(plotRange):
        if pointerLocation1 < G and pointerLocation2 < G:
            f1.seek(pointerLocation1)
            f2.seek(pointerLocation2)
        
            str1 = f1.read(windowSize)
            str2 = f2.read(windowSize)
            
            pointerLocation1 = pointerLocation1 + windowSize
            pointerLocation2 = pointerLocation2 + windowSize
            
            distance = hammingDistance(str1, str2, windowSize)
            distanceList.append(distance)        
        
    f = open(outputResult,'w')
    for eachitem in distanceList:  
        f.write(str(eachitem) + "\n" )
    f.close()
    
    #print "Forward List",  distanceList
    
#    windowSize = 10
#    distanceList = []
    
#    pointerLocation1 = startpt1 - windowSize
#    pointerLocation2 = startpt2 - windowSize


#    for index in range(plotRange):
#        if pointerLocation1 > 0 and pointerLocation2 >0:
#            f1.seek(pointerLocation1)
#            f2.seek(pointerLocation2)
        
#            str1 = f1.read(windowSize)
#            str2 = f2.read(windowSize)
            
#            pointerLocation1 = pointerLocation1 - windowSize
#            pointerLocation2 = pointerLocation2 - windowSize
            
#            distance = hammingDistance(str1, str2, windowSize)
#            distanceList.append(distance)        
        
#    f = open(outputResult+'2','w')
#    for eachitem in distanceList:    
#        f.write(str(eachitem) + "\n")
#    f.close()

    #print "Backward List", distanceList
    
    f1.close()
    f2.close()    

def generateOneLineFile(foldername, filepath):
    f = open(filepath, 'r') 
    
    f2 = open(foldername +"/oneLine.fasta", 'w')  
    f3= open(foldername +"/oneLine2.fasta", 'w' )
    
    temp = f.readline()
    
    while (len(temp) > 0):
        temp = f.readline()
        f2.write(temp[0:len(temp)-1])
        f3.write(temp[0:len(temp)-1])
    
    
    f.close()
    f2.close()
    f3.close()        

def batchComputationOfStatistics(folders):
    for eachfolder in folders:
        generateOneLineFile(eachfolder, eachfolder+ "/genome.fasta")
        batchPatternBeyondRepeat(eachfolder + "/long_repeats.txt", eachfolder + "/oneLine.fasta", eachfolder + "/oneLine2.fasta" , eachfolder+ "/result")
  


def extractApproxRepeatStatistics(folderName, topX):
    generateOneLineFile(folderName, folderName+ "/genome.fasta")
    listOfData = readFromMummerOutput(folderName+ "/long_repeats.txt")
    listOfData= sorted(listOfData,key = itemgetter(2), reverse = True)

    if listOfData[0][0] == 1 and listOfData[0][1] == 1:
        listOfData.pop(0)

    length = listOfData[0][2]

    for index in range(min(topX,len(listOfData))):
        start1 = listOfData[index][0]    
        start2 = listOfData[index][1]    

        
        print "Starting Location 1:" , start1
        print "Starting Location 2:" , start2
        print ""

        #reportPatternBeyondRepeat( folderName+ "/oneLine.fasta", folderName+ "/oneLine2.fasta", start1 + int( length/2) ,start1 + int( length/2) +1 , start2 + int( length/2) ,start2 + int( length/2) +1 , int(length*5/10),folderName+"/outputfile" + str(index)+".txt")
        reportPatternBeyondRepeat( folderName+ "/oneLine.fasta", folderName+ "/oneLine2.fasta", start1 - int( 1.5*length) ,start1 - int(  1.5*length) -1 , start2 - int( 1.5* length) ,start2 - int( 1.5* length) -1 , int(length*5/10),folderName+"/outputfile" + str(index)+".txt")
        
        readAndPlotStat(folderName+"/outputfile" + str(index)+".txt",folderName+"/outputfile" + str(index)+".txt2",folderName, index)

    
### Clustering of approx repeat 
def editDistance(str1, str2, windowSize, threshold):
    #print "Compute Edit Distance"
    rowLimit = len(str1)
    colLimit = len(str2)
    
    K = min(rowLimit,colLimit )
    listOfDiagonalPoints = []
    stoppingposition = 0
    
    lastRow = [0]
    lastCol = [0]
    
    currentRow = []
    currentCol = []
    
    
    for index in range(1,K):
        currentRow = [index]
        currentCol = [index]     
        
        for index2 in range(1,index+1):
            currentRow.append(-1)
            currentCol.append(-1)
            
            if str1[index] == str2[index2] :
                checkCurrentBase = 0 
            else:
                checkCurrentBase = 1
            
            if index2 != index:
                currentRow[index2] = min( currentRow[index2-1] +1 ,lastRow[index2] +1 , lastRow[index2-1] + checkCurrentBase )            
                currentCol[index2] = min( currentCol[index2-1] +1, lastCol[index2] +1, lastCol[index2-1]+ checkCurrentBase   )
            else:
                currentRow[index2] = min(currentRow[index2-1] +1 ,currentCol[index2-1] +1 , lastRow[index2-1] + checkCurrentBase )            
                currentCol[index2] = min(currentCol[index2-1] +1, currentRow[index2-1] +1, lastCol[index2-1]+ checkCurrentBase )
                
        
        listOfDiagonalPoints.append(currentRow[index])
         
        lastRow = []
        lastCol = []
            
        for eachitem in currentRow:
            lastRow.append(eachitem)
        for eachitem in currentCol:
            lastCol.append(eachitem)
        
        stoppingpositionFromStart = index
        numberOfMutation = listOfDiagonalPoints[-1]
        
        if len(listOfDiagonalPoints) > windowSize:
            if listOfDiagonalPoints[-1] -  listOfDiagonalPoints[-1-windowSize] > threshold:
                stoppingpositionFromStart = index
                numberOfMutation = listOfDiagonalPoints[-1]
                break
    
    return stoppingpositionFromStart, numberOfMutation

def reverseString(str1):
    str2 = ""    
    for index in range(len(str1)-1, -1, -1):
        str2 = str2 + str1[index]
    return str2
            
            
def findapproxrepeatLengthEditDistance(filename1,filename2, start1, start2, lengthOfExactRepeat):
    ### Finding editDistance RHS
    f1 = open(filename1, 'r')
    f2 = open(filename2, 'r')
    
    #G = f1.read()
    
    totalNumberOfError = 0
    
    windowSize = 100
    threshold = 30 

    lookaheadRage = 10000    
    ### Decision rule : if > 30 error in the latest window of 100, then stop counting 
    
    numberOfError= 0     
    
    ### Compute the RHS
    f1.seek(start1+ lengthOfExactRepeat- windowSize-1)
    f2.seek(start2+ lengthOfExactRepeat- windowSize-1)
    
    temp1  = f1.read(lookaheadRage)
    temp2 = f2.read(lookaheadRage)
        
    stoppingpositionFromStart, numberOfError = editDistance(temp1, temp2, windowSize, threshold)
    totalNumberOfError = totalNumberOfError + numberOfError
    
    ### Finding editDistance LHS
    endIndex1 = start1+ lengthOfExactRepeat- windowSize-1   +  stoppingpositionFromStart- int(threshold* 1.3333)
    endIndex2 = start1+ lengthOfExactRepeat- windowSize-1   +  stoppingpositionFromStart- int(threshold* 1.3333)


    numberOfError= 0     
    
    ### Compute the LHS
    f1.seek(start1-lookaheadRage )
    f2.seek(start2-lookaheadRage)
    
    f1.read(windowSize)
    f2.read(windowSize)
    
    temp1  = f1.read(lookaheadRage)
    temp2 = f2.read(lookaheadRage)

    temp1 = reverseString(temp1)
    temp2 = reverseString(temp2)
    
    stoppingpositionFromStart, numberOfError = editDistance(temp1, temp2, windowSize, threshold)
    #print     "stoppingpositionFromStart, numberOfError",stoppingpositionFromStart, numberOfError  
    totalNumberOfError = totalNumberOfError + numberOfError
    
    ### Combining results
    startIndex1 = start1- (- windowSize-1   +  stoppingpositionFromStart- int(threshold* 1.3333))
    startIndex2 = start2- (- windowSize-1   +  stoppingpositionFromStart- int(threshold* 1.3333))
    
    
    lapprox = endIndex1-startIndex1
    
    print lapprox, totalNumberOfError ,threshold,lengthOfExactRepeat
    mutationRate = (totalNumberOfError - 2* threshold )/float(lapprox)
    
    if mutationRate == 0:
        mutationRate = 1/float(lapprox)
    
    print "mutationRate", mutationRate, (totalNumberOfError - 2* threshold )
    
    
    if startIndex1 > startIndex2 :
        dummy = startIndex1
        startIndex1 = startIndex2 
        startIndex2 = dummy
    
    if lapprox <= lengthOfExactRepeat:
        return start1+1, start2+1, lengthOfExactRepeat+1, 1/ float(lengthOfExactRepeat)
    else:
        return  startIndex1, startIndex2, lapprox, mutationRate 


def findapproxrepeatLength(filename1, filename2, start1, start2, lengthOfExactRepeat):
    f1 = open(filename1, 'r')
    f2 = open(filename2, 'r')
    
    totalNumberOfError = 0
    
    windowSize = 100
    threshold = 30 
    ### Decision rule : if > 50 error in the latest window of 100, then stop counting 
    
    numberOfError= 0     
    
    ### Compute the RHS
    f1.seek(start1+ lengthOfExactRepeat- windowSize-1)
    f2.seek(start2+ lengthOfExactRepeat- windowSize-1)
    
    temp1  = f1.read(windowSize)
    temp2 = f2.read(windowSize)

    
    lastPosition1 = start1+ lengthOfExactRepeat- windowSize-1
    lastPosition2 = start2+ lengthOfExactRepeat- windowSize-1
    
    print "CheckPoint 1 : ", hammingDistance(temp1, temp2, len(temp1))
    numberOfError = hammingDistance(temp1, temp2, len(temp1))
    totalNumberOfError = totalNumberOfError + numberOfError
    
    while (numberOfError < threshold):
        
        f1.seek(lastPosition1)
        char1 = f1.read(1)
        f2.seek(lastPosition2)
        char2 = f2.read(1)
        #print numberOfError, char1, char2
        if char1 != char2 : 
          # print "Mutation"
           numberOfError = numberOfError - 1 
        
        f1.seek(lastPosition1  + windowSize)
        f2.seek(lastPosition2 + windowSize)
        
        char1 = f1.read(1)
        char2 = f2.read(1)
        
        if char1 != char2:
            numberOfError = numberOfError + 1
            totalNumberOfError= totalNumberOfError + 1
           # print ":-)"
        
        lastPosition1 = lastPosition1 + 1
        lastPosition2 = lastPosition2 + 1
    
    
    endIndex1 = lastPosition1 + windowSize - int(threshold* 1.3333)
    endIndex2 = lastPosition2 + windowSize - int(threshold* 1.3333)


    numberOfError= 0     
    
    ### Compute the LHS
    f1.seek(start1)
    f2.seek(start2)
    
    temp1  = f1.read(windowSize)
    temp2 = f2.read(windowSize)

    
    lastPosition1 = start1 + windowSize -1
    lastPosition2 = start2  + windowSize -1
    
    numberOfError = hammingDistance(temp1, temp2, len(temp1))
    print "checkPoint2 :",  hammingDistance(temp1, temp2, len(temp1))
    totalNumberOfError = totalNumberOfError + numberOfError
    
    while (numberOfError < threshold):
        
        f1.seek(lastPosition1)
        char1 = f1.read(1)
        f2.seek(lastPosition2)
        char2 = f2.read(1)
        

        if char1 != char2 : 
           numberOfError = numberOfError - 1 
        
        f1.seek(lastPosition1  - windowSize)
        f2.seek(lastPosition2 - windowSize)
        
        char1 = f1.read(1)
        char2 = f2.read(1)
        
        if char1 != char2:
            numberOfError = numberOfError + 1
            totalNumberOfError = totalNumberOfError +1
        
        
        lastPosition1 = lastPosition1 - 1
        lastPosition2 = lastPosition2 - 1
    
    
    startIndex1 = lastPosition1 - windowSize + int(threshold* 1.3333)
    startIndex2 = lastPosition2 - windowSize + int(threshold* 1.3333)
    
    lapprox = endIndex1-startIndex1 
    print lapprox, totalNumberOfError ,threshold,lengthOfExactRepeat
    mutationRate = (totalNumberOfError - 2* threshold )/float(lapprox)
    if mutationRate == 0:
        mutationRate = 1/float(lapprox)
    
    print "mutationRate", mutationRate, (totalNumberOfError - 2* threshold )
    
    
    if startIndex1 > startIndex2 :
        dummy = startIndex1
        startIndex1 = startIndex2 
        startIndex2 = dummy
    
    if lapprox <= lengthOfExactRepeat:
        return start1, start2, lengthOfExactRepeat+1, 1/ float(lengthOfExactRepeat)
    else:
        return  startIndex1, startIndex2, lapprox, mutationRate
    
    
    
def  filterSameApproxRepeat(listOfApproxRepeat):
    runningIndex = 0
    
    while ( runningIndex < len(listOfApproxRepeat) -1):
        if listOfApproxRepeat[runningIndex][0:2] == listOfApproxRepeat[runningIndex+1][0:2] or ( abs(listOfApproxRepeat[runningIndex][0] - listOfApproxRepeat[runningIndex+1][0] ) <5 and abs(listOfApproxRepeat[runningIndex][1] - listOfApproxRepeat[runningIndex+1][1] )<5 ):
            #print len(listOfApproxRepeat), runningIndex
            listOfApproxRepeat.pop(runningIndex)
        else:
            runningIndex = runningIndex + 1
            
    return listOfApproxRepeat 
    
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
    
    #runningindex  =0
    #for i,j in zip(xList,yList):
    #    ax.annotate(str(lapproxList[runningindex]),xy=(i,j ))
    #    runningindex = runningindex +1
    
    plt.axhline(y= 1.25, linewidth=2, color='k', ls= 'dotted')    
        
    
    plt.title("Classification of approximate repeat for " +  folderName,fontsize=defaultsize)
    plt.xlabel("Mutation Rate = number of mutation / ( l_approx - l_exact) ",fontsize=defaultsize)
    plt.ylabel("Stretch = \n   l_approx / l_exact ",fontsize=defaultsize)
    
    
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
    
    #fig = plt.figure()
    plt.subplot(212)
    
    if len(homologous) > 0:    
        homologous = sorted(homologous)
        #plt.bar(homologous,np.ones(len(homologous)), color='b')
        plt.hist(homologous,bins = range(0,longestlapprox+25,20),color='b')
    if len(nonhomologous) > 0:
        nonhomologous = sorted(nonhomologous)
        plt.hist(nonhomologous, bins  = range(0,longestlapprox+25,20), color = 'r', alpha = 0.5)
        #plt.bar(nonhomologous,np.ones(len(nonhomologous)),  color='r')
    

    plt.axvline(x= longestExactRepeat, linewidth=4, color='g', ls= 'dotted')

    defaultsize = 14
    
    plt.xlim((0, 8500))
    plt.ylim((0,10))
    
    plt.title("Approximate Repeat Spectrum for "+ folderName,fontsize=defaultsize)
    plt.xlabel('Length of approximate repeat',fontsize=defaultsize)
    plt.ylabel('Number of approximate repeat',fontsize=defaultsize)    


def clusteringRepeatMain(folderName, listOfData, topX, typeOfRepeat, errorType, longestExactRepeat, premapping):
    listOfApproxRepeat = []
    
    if len(listOfData) > 0:
        if listOfData[0][0] == 1 and listOfData[0][1] == 1:
            listOfData.pop(0)
    
    
    if longestExactRepeat == -1:
        longestExactRepeat = listOfData[0][2]
    print folderName
    print "longestExactRepeat",longestExactRepeat
    
    for index in range(min(topX,len(listOfData))):
        start1 = listOfData[index][0]    
        start2 = listOfData[index][1]    
        lengthOfExactRepeat = listOfData[index][2]
        
        print "Starting Location 1:" , start1
        print "Starting Location 2:" , start2
        print "length Of Exact Repeat", lengthOfExactRepeat

        
        if typeOfRepeat == 'c':
            inversePrefix = "/reverseoneLine.fasta"    
        else:
            inversePrefix = "/oneLine2.fasta"
        
        if errorType == 'h':
            startIndex1, startIndex2, lapprox, mutationRate = findapproxrepeatLength( folderName+ "/oneLine.fasta", folderName+inversePrefix , start1, start2, lengthOfExactRepeat)
        elif errorType == 'e':
            startIndex1, startIndex2, lapprox, mutationRate = findapproxrepeatLengthEditDistance( folderName+ "/oneLine.fasta", folderName+ inversePrefix, start1, start2, lengthOfExactRepeat)
            
            
        if len(premapping) == 0:
            listOfApproxRepeat.append([startIndex1, startIndex2, lapprox, mutationRate,lengthOfExactRepeat, start1, start2])
        else:
            searchtemp = ['eee',start1, start2,lengthOfExactRepeat ]
            tempindex = bisect.bisect(premapping,searchtemp)
            temptarget = premapping[tempindex]
            listOfApproxRepeat.append(temptarget[1:6])
            
            
        
    
    listOfApproxRepeat = sorted(listOfApproxRepeat)
    listOfApproxRepeat = filterSameApproxRepeat(listOfApproxRepeat)
    
    print "len(listOfApproxRepeat)", len(listOfApproxRepeat)
    ### Only plot the spectrum if we are doing simple repeat, ignore triple/intereleave repeats     
    if len(premapping) == 0 and typeOfRepeat == 'r':
        listOfApproxRepeat = sorted(listOfApproxRepeat, key = itemgetter(2),reverse = True)
        plotdetailedPatternForTopRepeats(folderName,topX, listOfApproxRepeat)  
        
    listOfApproxRepeat = sorted(listOfApproxRepeat)  

    plotClusteringOfRepeat(listOfApproxRepeat,folderName)
    
    plotApproxRepeatSpectrum(listOfApproxRepeat,longestExactRepeat,folderName) 
    
    folderHeader = ""
    if errorType == 'h':
        folderHeader = "dataHammingDistance/"+folderName
    elif errorType == 'e':
        folderHeader = "dataEditDistance/"+folderName
    
    if typeOfRepeat == 'r':
        folderHeader = folderHeader + "_1approxrepeatstat.png"
    elif typeOfRepeat == 'i':
        folderHeader = folderHeader + "_2approxinterleaverepeatstat.png"
    elif typeOfRepeat == 't': 
        folderHeader = folderHeader + "_3approxtriplerepeatstat.png"        
    elif typeOfRepeat == 'c':
        folderHeader= folderHeader + "_4approxinvertedrepeatstat.png"
    
    plt.savefig(folderHeader) 
    #plt.show()    
    
    ### Output to csv file
    # input format listOfApproxRepeat : (startIndex1, startIndex2, lapprox, mutationRate, lengthOfExactRepeat)
    # output format : (repeat index, exact length, approximate length, number mutations)


    ofile  = open(folderHeader+'.csv', "wb")
    #writer = csv.writer(ofile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_ALL)
    writer = csv.writer(ofile)   
    writer.writerow(["repeat index 1 ","repeat index 2", "exact length", "approximate length", "number of mutations"])
    for eachitem in listOfApproxRepeat:
        #print "eachitem", eachitem
        
        templist = [eachitem[0], eachitem[1], eachitem[4], eachitem[2], int(round(eachitem[3]*eachitem[2]))]
        #print "templist", templist
        writer.writerow(templist)

    ofile.close() 
    return listOfApproxRepeat

### Extract Triple/interleave Repeat
def preprocessingRepeatInfo(folderName):
    #dirList=os.listdir(folderName)
    dirList = [folderName+'/long_repeats.txt']
    dataList = []
    
    for fname in dirList:
        #f= open(folderName+"/"+fname, 'r')
        f= open(fname,'r')        
        f.readline()
        f.readline()
        temp = f.readline()
        while (len(temp) > 0 ):
            index1, index2, score = temp[0:len(temp)-1].split()
            
            index1ToAdd = int(index1)
            index2ToAdd = int(index2)
            scoreToAdd = int(score)

            dataList.append([fname,index1ToAdd,index2ToAdd,scoreToAdd])
            temp = f.readline()
        f.close()
    
    dataList = sorted(dataList, key= itemgetter(3), reverse = True)
     
    return dataList

def computeInterleaveRepeat(currentCluster):
    
   if len(currentCluster) == 0:
         return []
   else:
        filename = currentCluster[0][0] 
        lengthOfCluster = len(currentCluster)
        
        ### Output Format  : [file_name, x1, y1, repeat_length_1, x2, y2, repeat_length_2]
        interleaveRepeatList = []
        for i in range(lengthOfCluster-1):
            (x1, y1, repeat1) = currentCluster[i][1], currentCluster[i][2], currentCluster[i][3]
            
            for j in range(i+1, lengthOfCluster): 
                (x2, y2, repeat2) = currentCluster[j][1], currentCluster[j][2], currentCluster[j][3]
                #print x1, x2, y1, y2
                if ( x1 < x2 and x2 < y1 and y1  <y2)  or (x2 < x1 and x1  < y2 and y2 < y1 ):                    
                    interleaveRepeatList.append([filename,x1, y1, repeat1, x2, y2, repeat2, min(repeat1, repeat2) ])
        
        #print "Interleave Repeat"
        
        #for eachitem in interleaveRepeatList:
        #    print eachitem
            
        return interleaveRepeatList

def computeTripleRepeat(currentCluster):
   if len(currentCluster) == 0:
         return []
   else:
        filename = currentCluster[0][0] 
        lengthOfCluster = len(currentCluster)
        
        ### Output Format  : [file_name, x1, y1, repeat_length_1, x2, y2, repeat_length_2]
        tripleRepeatList = []
        for i in range(lengthOfCluster-1):
            (x1, y1, repeat1) = currentCluster[i][1], currentCluster[i][2], currentCluster[i][3]
            
            for j in range(i+1, lengthOfCluster): 
                (x2, y2, repeat2) = currentCluster[j][1], currentCluster[j][2], currentCluster[j][3]
                #print x1, x2, y1, y2
                if (x2 >= y1 and x2 + repeat2 <=y1 + repeat1 and  (x2 != y1 or y2!= x1))  or (x2 >= x1 and x2 + repeat2 <= x1+ repeat1 and  (x2 != x1 or y2 != y1)):
                    tripleRepeatList.append([filename,x1, y1, repeat1, x2, y2, repeat2, min(repeat1, repeat2) ])
        
        #print "Triple Repeat: "
        #for eachitem in tripleRepeatList:
            #print eachitem
        return tripleRepeatList


def findDifferentTypesOfRepeat(dataList,topX):

    dataList = sorted(dataList)

    ### Input Data Format : [ file_name , startLoc_1 , startLoc_2, repeat_length]    
    ### Find Interleave repeat and triple repeat
    
    metaListForInterleaveRepeat = []
    metaListForTripleRepeat = []
    prevName = ""
    currentCluster = [] 
    for eachitem in dataList:
        currentFileName = eachitem[0]
        if currentFileName != prevName:            
            templistOfInterRepeat =      computeInterleaveRepeat(currentCluster)
            if len(templistOfInterRepeat) > 0: 
                metaListForInterleaveRepeat= metaListForInterleaveRepeat + templistOfInterRepeat
                
            templistOfTripleRepeat = computeTripleRepeat(currentCluster)
            if len(templistOfTripleRepeat) > 0: 
                metaListForTripleRepeat = metaListForTripleRepeat + templistOfTripleRepeat
            
            currentCluster = []
            currentCluster.append(eachitem)
            prevName = currentFileName
            
        else:
            currentCluster.append(eachitem)
        
            
    templistOfInterRepeat =      computeInterleaveRepeat(currentCluster)
    if len(templistOfInterRepeat) > 0: 
        metaListForInterleaveRepeat= metaListForInterleaveRepeat + templistOfInterRepeat
        
    templistOfTripleRepeat = computeTripleRepeat(currentCluster)
    if len(templistOfTripleRepeat) > 0: 
        metaListForTripleRepeat = metaListForTripleRepeat + templistOfTripleRepeat
    
    metaListForInterleaveRepeat = sorted(metaListForInterleaveRepeat, key= itemgetter(-1), reverse = True )
    metaListForTripleRepeat = sorted(metaListForTripleRepeat, key= itemgetter(-1), reverse = True) 
    
    print "metaListForInterleaveRepeat"
    for index in range(len(metaListForInterleaveRepeat)):
        tempdata = metaListForInterleaveRepeat[index] 
        print tempdata
    print "metaListForTripleRepeat"
    for index in  range(len(metaListForTripleRepeat)):
        print metaListForTripleRepeat[index]
        
    return metaListForInterleaveRepeat[0:min(len(metaListForInterleaveRepeat),topX)], metaListForTripleRepeat[0:min(len(metaListForTripleRepeat),topX)]



def extractInterleaveRepeat(maxInterleaveRepeatLength):
    listOfInterleaveRepeat = []
    for eachitem in maxInterleaveRepeatLength:

        if eachitem[3] == eachitem[7]:
            if eachitem[1] < eachitem[2]:
                listOfInterleaveRepeat.append(eachitem[1:4])
            else: 
                listOfInterleaveRepeat.append([eachitem[2], eachitem[1], eachitem[3]])
        if eachitem[6] == eachitem[7]:
            if eachitem[4] < eachitem[5]:
                listOfInterleaveRepeat.append(eachitem[4:7])
            else:
                listOfInterleaveRepeat.append([eachitem[5], eachitem[4], eachitem[6]])
                
    
    listOfInterleaveRepeat = sorted(listOfInterleaveRepeat)
    listOfInterleaveRepeat = filterSameApproxRepeat(listOfInterleaveRepeat)
    
    listOfInterleaveRepeat = sorted(listOfInterleaveRepeat,key = itemgetter(2) ,reverse = True)
    
    return listOfInterleaveRepeat
            
def extractExactRepeats(folderName):
    generateOneLineFile(folderName, folderName+ "/genome.fasta")    
    genomeInversion(folderName)
    
    listOfData = readFromMummerOutput(folderName+ "/long_repeats.txt")
    listOfData= sorted(listOfData,key = itemgetter(2), reverse = True)
    
    return listOfData
    
        

def combinedStatisticsPlot(folderName,topX,typeOferror):
    
    dataList = preprocessingRepeatInfo(folderName)
    #print "dataList" , dataList
    metaListForInterleaveRepeat, metaListForTripleRepeat = findDifferentTypesOfRepeat(dataList,topX)
       
    maxInterleaveRepeatLength = metaListForInterleaveRepeat[0][-1]
    maxTripleRepeatLength = metaListForTripleRepeat[0][-1]
    
    print "maxInterleaveRepeatLength,maxTripleRepeatLength", maxInterleaveRepeatLength,maxTripleRepeatLength  

    listOfRepeat = extractExactRepeats(folderName)
    print listOfRepeat
    listOfApproxRepeat = clusteringRepeatMain(folderName, listOfRepeat, topX, 'r', typeOferror , -1, [])



    for eachitem in listOfApproxRepeat:
        #startIndex1, startIndex2, lapprox, mutationRate,lengthOfExactRepeat, start1, start2
        if eachitem[2] < eachitem[4]:
            print "Errorrrrrrrrrrrrrrrr", eachitem[2] , eachitem[4]
        eachitem.insert(0,'eee')
     
    print "listOfApproxRepeat", listOfApproxRepeat
    metaListForInterleaveRepeat, metaListForTripleRepeat = findDifferentTypesOfRepeat(listOfApproxRepeat,topX)    
    print "metaListForInterleaveRepeat, metaListForTripleRepeat",metaListForInterleaveRepeat, metaListForTripleRepeat
    
    listOfInterleaveRepeat = extractInterleaveRepeat(metaListForInterleaveRepeat)
    clusteringRepeatMain(folderName, listOfInterleaveRepeat, topX, 'i',typeOferror, maxInterleaveRepeatLength,listOfApproxRepeat)

    listOfTripleRepeat = extractInterleaveRepeat(metaListForTripleRepeat)
    clusteringRepeatMain(folderName, listOfTripleRepeat, topX, 't',typeOferror,maxTripleRepeatLength,listOfApproxRepeat)

    #listOfRepeat = readFromMummerOutput(folderName+ "/long_repeats_inverted.txt")
    #listOfRepeat= sorted(listOfRepeat,key = itemgetter(2), reverse = True)
    #listOfInvertedRepeat = clusteringRepeatMain(folderName, listOfRepeat, topX, 'c', typeOferror , -1, [])   




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

### Inverted Repeat Treatment
def genomeInversion(folderName):
    fin = open(folderName  + "/oneLine.fasta", 'r')
    fout = open(folderName +"/reverseoneLine.fasta", 'w')

    temp = fin.read()   
    G = len(temp)
    
    for index in range(G-1, -1, -1):    
        fin.seek(index)
        temp = fin.read(1)
        if temp == 'A':
            fout.write('T')
            
        elif temp == 'C':
            fout.write('G')
            
        elif temp == 'T':
            fout.write('A')
            
        elif temp == 'G':
            fout.write('C')
        else:
            print "Error", temp
            
    
    fout.close()
    fin.close()
#detailRepeatPattern("Ecoli536/oneLine.fasta","Ecoli536/oneLine2.fasta", 228619, 4419727, 3353,50,"outputResult.txt")
#detailRepeatPattern("testGenome.txt","testGenome2.txt", 20, 48, 9,10,"outputResult.txt")
def captureSegment():
    f = open("segmentEK12\\oneLine.fasta", 'r')
    fout = open("segmentEK12\\genome.fasta", 'w')    
    start , end =2724717 ,3468075 + 2167   
    bufferR = 100    
    f.seek( start- bufferR)
    temp = f.read(end-start+bufferR)
    fout.write(">Seg1 \n")    
    fout.write(temp)    
    f.close()

def captureSegment2():
    f = open("LactoAcid\\oneLine.fasta",'r')
    f2 = open("LactoAcid\\reverseOneLine.fasta",'r')
    fout = open("segmentedLactoAcid.fasta", 'w')    


    offset = 50000
    start2, end2 = 363225 - offset ,363225 + offset
    start1, end1 = 469122 - offset ,469122 + offset 
    
    fout.write(">Seg1\n")
        
    f.seek(start1 )
    temp = f.read(2*offset)
    fout.write(temp)
    
    f2.seek(start2)
    temp = f2.read(2*offset)
    fout.write(temp)
    
    f.close()
    f2.close()
    
    fout.close()
    
def generateCombineComplement():
    f = open("EcoliK12\\oneLine.fasta", 'r')
    f2 = open("EcoliK12\\reverseoneLine.fasta", 'r')
    fout = open("ecolik12Combined.fasta",'w')
    
    temp = f2.read()
    fout.write(temp)    
    temp = f.read() 
    fout.write(temp)
    
    f.close()
    f2.close()
    fout.close()


#folders = ["LactoAcid","Ecoli536","Salmonella"]    
#folders = ["segmentEK12"]




#for item in folders:
#    print item
#    genomeInversion(item)
#    plotdetailedPatternForTopRepeats(item, 20)       
#    combinedStatisticsPlot(item,30,'h')
    #combinedStatisticsPlot(item,20,'e')




def pipeLine():
    folderList = range(1)
    
    for eachitem in folderList:
        foldername = str(eachitem) 
        os.system("cp "+foldername +"\genome.fasta oneLine.fasta")
        os.system("tr -d '\n' < "+foldername + "\oneLine.txt")
        os.system("cp "+foldername + "\oneLine.txt"+ " oneLine2.txt")    
        
        path = ""
        os.system(path+"mummer -maxrepeat genome.fasta genome.fasta > " + foldername+"\long_repeats.txt")
        
        combinedStatisticsPlot(foldername,30,'h')
    
#pipeLine().
def viewdata():
    dataset = readFromMummerOutput("special//output.txt")
    dataset = sorted(dataset)
    for eachitem in dataset:
        print eachitem

#viewdata()
#generateCombineComplement()   
#extractApproxRepeatStatistics("LactoAcid", 10)    
#plotPatternBeyondRepeat("pacBioEcoli/oneLine.fasta","pacBioEcoli/oneLine.fasta", 281-200, 281 - 200 + 3990, 9281 -200, 9281 -200 + 3990,100)
#folders = ["Yesnina"]  
#folders = ['HumanChr14','HumanChr22','Buchnera','HumanChr19','SAureus','Salmonella','Ecoli536','Pmarinus','Sislandicus','Heli51','Yesnina','RSPHAEROIDES', 'LactoAcid'] 
#folders = ["LactoAcid","Pmarinus", "SAureus","Ecoli536","Salmonella","Yesnina","Chromosome22"]

#folders = ["LactoAcid"]
#batchComputationOfStatistics(folders)

#    combinedStatisticsPlot(item, 20)

#tester("SAureus")
#batchPatternBeyondRepeat("pacBioEcoli/mummer", "pacBioEcoli/oneLine.fasta", "pacBioEcoli/oneLine2.fasta", "pacBioEcoli/result")
#batchGenerateApproximateRepeatStat("oneLine.fasta","oneLine2.fasta","ecoli536",2)
start1, start2, lengthRep , offset =1193350,2656552,1217, 2000
plotPatternBeyondRepeat("Meiothermus\\oneLine.fasta","Meiothermus\\oneLine2.fasta", start1-offset, start1-offset + 1, start2-offset, start2-offset + 1,1000)

#plotPatternBeyondRepeat("pacBioEcoli/oneLine.fasta","pacBioEcoli/oneLine2.fasta", 959 - 900, 959 - 900 + 1, 4643451- 900 , 4643451 - 900 +1,1000)
#batchGenerateApproximateRepeatStatIndel("3\oneLine.fasta","3\oneLine2.fasta","ecoli536",8)
#testIndelStatistics()