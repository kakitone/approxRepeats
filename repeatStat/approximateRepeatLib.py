from operator import itemgetter
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
import numpy as np
import os
import csv
import bisect
import random

import graphPlottingLib
import distanceComputeLib
import ioLib


### Complete Spectrum of Approximate repeat(Hamming distance type and edit distance type)

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
    
    
### Finding Approx Repeat 
def editDistance(str1, str2, windowSize, threshold):
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
    
    return stoppingpositionFromStart, numberOfMutation,listOfDiagonalPoints
            
def findapproxrepeatLengthEditDistance(filename1,filename2, start1, start2, lengthOfExactRepeat):
    ### Finding editDistance RHS
    f1 = open(filename1, 'r')
    f2 = open(filename2, 'r')
    
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
        
    stoppingpositionFromStart, numberOfError,listOfDiagonalPointsRight = editDistance(temp1, temp2, windowSize, threshold)
    

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

    temp1 = ioLib.reverseString(temp1)
    temp2 = ioLib.reverseString(temp2)
    
    stoppingpositionFromStart, numberOfError,listOfDiagonalPointsLeft = editDistance(temp1, temp2, windowSize, threshold)
 
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
    
    
    outputList = [[0,0]]
    
    ### Prefilter

    runningIndex = 0
    while runningIndex <= len(outputList) -1:
        if listOfDiagonalPointsLeft[runningIndex] == 0 :
            listOfDiagonalPointsLeft.pop(runningIndex)
        else:
            runningIndex = runningIndex +1
     
    runningIndex = 0
    while runningIndex <= len(outputList) -1:
        if listOfDiagonalPointsRight[runningIndex] == 0 :
            listOfDiagonalPointsRight.pop(runningIndex)
        else:
            runningIndex = runningIndex +1
            
            
    
    for index in range(len(listOfDiagonalPointsLeft)):
        outputList.insert(0,[-listOfDiagonalPointsLeft[index] ,index])
        
    for index in range(len(listOfDiagonalPointsRight)):
        outputList.append([listOfDiagonalPointsRight[index],index ])
    
     ### Filtering outputList
    runningIndex = 1
    while runningIndex <= len(outputList) -1:
        if outputList[runningIndex][0] == outputList[runningIndex-1][0] and outputList[runningIndex] != [0,0]:
            outputList.pop(runningIndex)
        else:
            runningIndex = runningIndex +1
     
        
    if startIndex1 > startIndex2 :
        dummy = startIndex1
        startIndex1 = startIndex2 
        startIndex2 = dummy
    
    if lapprox <= lengthOfExactRepeat:
        return start1+1, start2+1, lengthOfExactRepeat+1, 1/ float(lengthOfExactRepeat),outputList
    else:
        return  startIndex1, startIndex2, lapprox, mutationRate ,outputList


def findapproxrepeatLength(filename1, filename2, start1, start2, lengthOfExactRepeat):
    f1 = open(filename1, 'r')
    f2 = open(filename2, 'r')
    
    totalNumberOfError = 0
    
    windowSize = 100
    threshold = 25 
    ### Decision rule : if > 50 error in the latest window of 100, then stop counting 
    
    numberOfError= 0     
    
    ### Compute the RHS
    f1.seek(start1+ lengthOfExactRepeat- windowSize-1)
    f2.seek(start2+ lengthOfExactRepeat- windowSize-1)
    
    temp1  = f1.read(windowSize)
    temp2 = f2.read(windowSize)

    
    lastPosition1 = start1+ lengthOfExactRepeat- windowSize-1
    lastPosition2 = start2+ lengthOfExactRepeat- windowSize-1
    
    print "CheckPoint 1 : ", distanceComputeLib.hammingDistance(temp1, temp2, len(temp1))
    numberOfError = distanceComputeLib.hammingDistance(temp1, temp2, len(temp1))
    totalNumberOfError = totalNumberOfError + numberOfError
    
    while (numberOfError < threshold):
        
        f1.seek(lastPosition1)
        char1 = f1.read(1)
        f2.seek(lastPosition2)
        char2 = f2.read(1)

        if char1 != char2 : 
           numberOfError = numberOfError - 1 
        
        f1.seek(lastPosition1  + windowSize)
        f2.seek(lastPosition2 + windowSize)
        
        char1 = f1.read(1)
        char2 = f2.read(1)
        
        if char1 != char2:
            numberOfError = numberOfError + 1
            totalNumberOfError= totalNumberOfError + 1

        
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
    
    numberOfError = distanceComputeLib.hammingDistance(temp1, temp2, len(temp1))
    print "checkPoint2 :",  distanceComputeLib.hammingDistance(temp1, temp2, len(temp1))
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
    #runningIndex = 0
    
    #while ( runningIndex < len(listOfApproxRepeat) -1):
    #    if listOfApproxRepeat[runningIndex][0:2] == listOfApproxRepeat[runningIndex+1][0:2] or ( abs(listOfApproxRepeat[runningIndex][0] - listOfApproxRepeat[runningIndex+1][0] ) <20 and abs(listOfApproxRepeat[runningIndex][1] - listOfApproxRepeat[runningIndex+1][1] )<20 ):
    #        #print len(listOfApproxRepeat), runningIndex
    #        listOfApproxRepeat.pop(runningIndex)
    #    else:
    #        runningIndex = runningIndex + 1
    
    listOfApproxRepeatNew = []   
    print "IMPORTANT", listOfApproxRepeat
    
    for index1 in range(len(listOfApproxRepeat))       :
        count = 0
        for index2 in range(index1 +1,len(listOfApproxRepeat) ) : 
             if listOfApproxRepeat[index1][0:2] == listOfApproxRepeat[index2][0:2] or ( abs(listOfApproxRepeat[index1][0] - listOfApproxRepeat[index2][0] ) <20 and abs(listOfApproxRepeat[index1][1] - listOfApproxRepeat[index2][1] )<20 ):
                count = count +1
        if count == 0:
            listOfApproxRepeatNew.append(listOfApproxRepeat[index1])
            
    print "IMPORTANT2 ", listOfApproxRepeatNew
    return listOfApproxRepeatNew 

def clusteringRepeatMain(folderName, listOfData, topX, typeOfRepeat, errorType, longestExactRepeat, premapping):
    
    
    print "typeOfRepeat", typeOfRepeat
    listOfApproxRepeat = []

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
            startIndex1, startIndex2, lapprox, mutationRate, outputList = findapproxrepeatLengthEditDistance( folderName+ "/oneLine.fasta", folderName+ inversePrefix, start1, start2, lengthOfExactRepeat)
            
            
        if len(premapping) == 0:
            listOfApproxRepeat.append([startIndex1, startIndex2, lapprox, mutationRate,lengthOfExactRepeat, start1, start2])
        else:
            searchtemp = ['eee',start1, start2,lengthOfExactRepeat ]
            tempindex = bisect.bisect(premapping,searchtemp)
            temptarget = premapping[tempindex]
            listOfApproxRepeat.append(temptarget[1:6])
            
    listOfApproxRepeat = sorted(listOfApproxRepeat, key = itemgetter(4))
    listOfApproxRepeat = filterSameApproxRepeat(listOfApproxRepeat)
    
    print "len(listOfApproxRepeat)", len(listOfApproxRepeat)
    
    if len(listOfApproxRepeat) == 0:
        return listOfApproxRepeat
    ### Only plot the spectrum if we are doing simple repeat, ignore triple/intereleave repeats     
    if len(premapping) == 0 and typeOfRepeat == 'r' and errorType == 'h' :
        listOfApproxRepeat = sorted(listOfApproxRepeat, key = itemgetter(2),reverse = True)
        plotdetailedPatternForTopRepeats(folderName,topX, listOfApproxRepeat,False)  
        
        plotdetailedPatternForTopRepeats2(folderName,topX, listOfApproxRepeat,False)  
        
    if len(premapping) == 0 and typeOfRepeat == 'c' and errorType == 'h':
        listOfApproxRepeat = sorted(listOfApproxRepeat, key = itemgetter(2),reverse = True)
        plotdetailedPatternForTopRepeats(folderName,topX, listOfApproxRepeat, True)  
        plotdetailedPatternForTopRepeats2(folderName,topX, listOfApproxRepeat, True)  

    if len(premapping) == 0 and typeOfRepeat == 'r' and errorType == 'e' :
        listOfApproxRepeat = sorted(listOfApproxRepeat, key = itemgetter(2),reverse = True)
        plotdetailedPatternForTopRepeats2Edit(folderName,topX, listOfApproxRepeat,False)  
        
    if len(premapping) == 0 and typeOfRepeat == 'c' and errorType == 'e':
        listOfApproxRepeat = sorted(listOfApproxRepeat, key = itemgetter(2),reverse = True)
        plotdetailedPatternForTopRepeats2Edit(folderName,topX, listOfApproxRepeat, True)  

        
    listOfApproxRepeat = sorted(listOfApproxRepeat)  

    graphPlottingLib.plotClusteringOfRepeat(listOfApproxRepeat,folderName)
    
    graphPlottingLib.plotApproxRepeatSpectrum(listOfApproxRepeat,longestExactRepeat,folderName) 
    
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
    writer = csv.writer(ofile)   
    writer.writerow(["repeat index 1 ","repeat index 2", "exact length", "approximate length", "number of mutations"])
    for eachitem in listOfApproxRepeat:        
        templist = [eachitem[0], eachitem[1], eachitem[4], eachitem[2], int(round(eachitem[3]*eachitem[2]))]
        writer.writerow(templist)

    ofile.close() 
    return listOfApproxRepeat

### Viewpoint : Distance between consecutive rise
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
def detailRepeatPattern2(genomeSource1,genomeSource2, startpt1, startpt2, lrepeat,mutationTolerant,outputResult):
    
    print startpt1, startpt2, lrepeat
    f1 = open(genomeSource1,'r')
    f2 = open(genomeSource2,'r')
    
    n1  =0 

    rightSideList = [0]
    rightStartPoint1 = startpt1 + lrepeat 
    rightStartPoint2 = startpt2 + lrepeat
    runningextension = 1
    
    while n1 < mutationTolerant:
        f1.seek(rightStartPoint1)
        f2.seek(rightStartPoint2)
        temp1 = f1.read(1)
        temp2 = f2.read(1)
        
        if temp1 != temp2 :
            n1 = n1 +1
            rightSideList.append(runningextension)
            
            
        rightStartPoint1 = rightStartPoint1 +1
        rightStartPoint2 = rightStartPoint2 +1
        runningextension = runningextension +1


    n2 = 0    
    leftSideList = [0]
    
    leftStartPoint1 = startpt1  -1
    leftStartPoint2 = startpt2 -1
    runningextension = 1
    while n2 < mutationTolerant:
        f1.seek(leftStartPoint1)
        f2.seek(leftStartPoint2)
        temp1 = f1.read(1)
        temp2 = f2.read(1)
        
        if temp1 != temp2 :
            n2 = n2 +1
            leftSideList.append(runningextension)
            
            
        leftStartPoint1 = leftStartPoint1 -1
        leftStartPoint2 = leftStartPoint2 -1
        runningextension = runningextension +1
    
    
    outputList = []
    
    for index in range(1,len(leftSideList)):
        outputList.insert(0,[-index,leftSideList[index]])
    
    for index in range(len(rightSideList)):
        outputList.append([index,rightSideList[index]])

    #plt.plot(range( 0,-len(leftSideList),-1), leftSideList)
    #plt.plot(range(len(rightSideList)), rightSideList)
    
    f1.close()
    f2.close()
    return outputList        
    
    
def plotdetailedPatternForTopRepeats2(folderName, topX, dataList,inverted = False):
    
    print dataList    
    ### Input format : [startIndex1, startIndex2, lapprox, mutationRate,lengthOfExactRepeat,start1,start2])
    dataList = sorted(dataList, key = itemgetter(2), reverse = True)
    dataList[min(topX, len(dataList)):]    = []
    dataList = sorted(dataList, key = itemgetter(2))

    
    lengthOfRepeatList = []
    for plotindex  in range(2):
        plt.figure(figsize=(16.0, 10))
        plt.grid(which = 'both')
        if plotindex  == 0:
            plt.xlim(-0.25,0.25)
            plt.ylim(0,6)
        
        #assert(1==2)
        colors= cm.rainbow(np.linspace(0, 1, min(topX,len(dataList))))
        #colors  = []
        #for index in range(len(reversecolors)-1, -1,-1):
        #    colors.append(reversecolors[index])
        ax = plt.subplot(1,1,1)
        counterindex  =0 
        
        if inverted:
            folderName2 = folderName + "_inverted_"
        else:
            folderName2 = folderName + "_noninverted_"
                        
        
        ofile  = open("regionBeyondRepeat/data_HD_twoSides_"+folderName2+"_"+str(counterindex)+".csv", "wb")
        writer = csv.writer(ofile)   
        writer.writerow(["Repeat ID","Extension ID", "Length Of Exact Repeat", "Length Of Approx Repeat"])
        
        for eachitem ,c in zip( dataList ,colors):
            if inverted == False :
                outputList =     detailRepeatPattern2(folderName+"/oneLine.fasta",folderName+"/oneLine2.fasta", eachitem[5], eachitem[6], eachitem[4],100,folderName+"/outputResult.txt")
            else:
                outputList =     detailRepeatPattern2(folderName+"/oneLine.fasta",folderName+"/reverseoneLine.fasta", eachitem[5], eachitem[6], eachitem[4],100,folderName+"/outputResult.txt")
            lrepeat = eachitem[4]
            
            plotY = []
            plotX = []
            for inneritem in outputList:
                plotY.append(inneritem[1]/float(lrepeat))
                plotX.append(inneritem[0]/float(lrepeat))
            
            plotY2 = []
            for elementindex in range(len(plotY)-1):
                plotY2.append( (plotY[elementindex+1] - plotY[elementindex] )/ float(plotX[elementindex+1] - plotX[elementindex]))
            
            lengthOfRepeatList.append(str(eachitem[2])+" , "+ str(eachitem[4]))
            plt.plot(plotX, plotY,linewidth=1, marker = 'x', color = c,label= str(eachitem[2])+" , "+ str(eachitem[4]))
         
            for index in range(len(plotX)):        
                writer.writerow([counterindex, int(plotX[index]*lrepeat), lrepeat,int(plotY[index]*lrepeat)])
            
            
            counterindex = counterindex +1 
        
        ofile.close()  
        defaultsize = 24
        
        lims = plt.xlim()    
    
        plt.tick_params(axis='both', which='major', labelsize=defaultsize)
        plt.tick_params(axis='both', which='minor', labelsize=defaultsize)
    
       
        plt.title("Analysis of region beyond repeat in " +  folderName +" for top"+str(topX) +" Exact Repeats",fontsize=defaultsize)
        plt.xlabel("Cumulative Hamming distance / l_repeat ",fontsize=defaultsize)
        plt.ylabel("Stretch =  l_extended-repeat / l_repeat ",fontsize=defaultsize)
        handles, labels = ax.get_legend_handles_labels()

    # reverse the order
        ax.legend(handles[::-1], labels[::-1],prop={'size':8})
        #plt.legend(lengthOfRepeatList)

 
        if inverted == False:
            if plotindex  == 0:       
                plt.savefig("regionBeyondRepeat/FixScale_twoSides/"+folderName+"_approxRepeatAnalysisPlot.png") 
            else:
                plt.savefig("regionBeyondRepeat/bestFit_twoSides/"+folderName+"_approxRepeatAnalysisPlot.png") 

        else:
            if plotindex  == 0:       
                plt.savefig("regionBeyondRepeat/FixScale_twoSides/"+folderName+"_inverted_approxRepeatAnalysisPlot.png") 
            else:
                plt.savefig("regionBeyondRepeat/bestFit_twoSides/"+folderName+"_inverted_approxRepeatAnalysisPlot.png") 
  
  
  


       
def plotdetailedPatternForTopRepeats(folderName, topX, dataList,inverted = False):
    ### Input format : [startIndex1, startIndex2, lapprox, mutationRate,lengthOfExactRepeat,start1,start2])
    print dataList    
    ### Input format : [startIndex1, startIndex2, lapprox, mutationRate,lengthOfExactRepeat,start1,start2])
    dataList = sorted(dataList, key = itemgetter(2), reverse = True)
    dataList[min(topX, len(dataList)):]    = []
    dataList = sorted(dataList, key = itemgetter(2))    
    
    lengthOfRepeatList = []
    for plotindex  in range(2):
        plt.figure(figsize=(16.0, 10))
        plt.grid(which = 'both')
        
        if plotindex  == 0:
            plt.xlim(0,0.45)
            plt.ylim(0,8)
        
        counterindex = 0
        
        
        if inverted:
            folderName2 = folderName + "_inverted_"
        else:
            folderName2 = folderName + "_noninverted_"
        ofile  = open("regionBeyondRepeat/data_HD_greedy_"+folderName2+"_"+str(counterindex)+".csv", "wb")
        writer = csv.writer(ofile)   
        writer.writerow(["Repeat ID","Mutation ID", "Length Of Exact Repeat", "Length Of Approx Repeat"])
        
        colors= cm.rainbow(np.linspace(0, 1, min(topX,len(dataList))))
        #colors  = []
        #for index in range(len(reversecolors)-1, -1,-1):
        #    colors.append(reversecolors[index])
        ax = plt.subplot(1,1,1)
        counterindex  =0 
        for eachitem ,c in zip( dataList ,colors):        
            if inverted == False :
                outputList =     detailRepeatPattern(folderName+"/oneLine.fasta",folderName+"/oneLine2.fasta", eachitem[5], eachitem[6], eachitem[4],50,folderName+"/outputResult.txt")
            else:
                outputList =     detailRepeatPattern(folderName+"/oneLine.fasta",folderName+"/reverseoneLine.fasta", eachitem[5], eachitem[6], eachitem[4],50,folderName+"/outputResult.txt")
            lrepeat = eachitem[4]
            
            plotY = []
            plotX = []
            for inneritem in outputList:
                plotY.append(inneritem[1]/float(lrepeat)-1)
                plotX.append(inneritem[0]/float(lrepeat))
            
            plotY2 = []
            for elementindex in range(len(plotY)-1):
                plotY2.append( (plotY[elementindex+1] - plotY[elementindex] )/ float(plotX[elementindex+1] - plotX[elementindex]))

            plt.plot(plotX, plotY,linewidth=1, marker = 'x',color = c,label= str(eachitem[2])+" , "+ str(eachitem[4]))
            lengthOfRepeatList.append(str(eachitem[2])+" , "+ str(eachitem[4]))
            
            for index in range(len(plotX)):        
                writer.writerow([counterindex, index, lrepeat,int(plotY[index]*lrepeat)])
            

            counterindex = counterindex +1 
        
        ofile.close() 
        defaultsize = 24
        
        lims = plt.xlim()    
    
        plt.tick_params(axis='both', which='major', labelsize=defaultsize)
        plt.tick_params(axis='both', which='minor', labelsize=defaultsize)
    
       
        plt.title("Analysis of region beyond repeat in " +  folderName +" for top"+str(topX) +" Exact Repeats",fontsize=defaultsize)
        plt.xlabel("Cumulative Hamming distance / l_repeat ",fontsize=defaultsize)
        plt.ylabel("Stretch =  l_extended-repeat / l_repeat ",fontsize=defaultsize)
        handles, labels = ax.get_legend_handles_labels()

    # reverse the order
        ax.legend(handles[::-1], labels[::-1],prop={'size':8})
        #plt.legend(lengthOfRepeatList)

        if inverted == False:
            if plotindex  == 0:       
                plt.savefig("regionBeyondRepeat/FixScale_greedy/"+folderName+"_approxRepeatAnalysisPlot.png") 
            else:
                plt.savefig("regionBeyondRepeat/bestFit_greedy/"+folderName+"_approxRepeatAnalysisPlot.png") 

        else:
            if plotindex  == 0:       
                plt.savefig("regionBeyondRepeat/FixScale_greedy/"+folderName+"_inverted_approxRepeatAnalysisPlot.png") 
            else:
                plt.savefig("regionBeyondRepeat/bestFit_greedy/"+folderName+"_inverted_approxRepeatAnalysisPlot.png") 



def  plotdetailedPatternForTopRepeats2Edit(folderName, topX, dataList,inverted = False):
    print dataList    
    ### Input format : [startIndex1, startIndex2, lapprox, mutationRate,lengthOfExactRepeat,start1,start2])
    dataList = sorted(dataList, key = itemgetter(2), reverse = True)
    dataList[min(topX, len(dataList)):]    = []
    dataList = sorted(dataList, key = itemgetter(2))
    
    lengthOfRepeatList = []
    for plotindex  in range(2):
        plt.figure(figsize=(16.0, 10))
        plt.grid(which = 'both')
        if plotindex  == 0:
            plt.xlim(-0.25,0.25)
            plt.ylim(0,6)
        
        colors= cm.rainbow(np.linspace(0, 1, min(topX,len(dataList))))

        ax = plt.subplot(1,1,1)
        
        counterindex =0 
        
        
        if inverted:
            folderName2 = folderName + "_inverted_"
        else:
            folderName2 = folderName + "_noninverted_"            
        
        ofile  = open("regionBeyondRepeat/data_ED_twoSides_"+folderName2+"_"+str(counterindex)+".csv", "wb")
        writer = csv.writer(ofile)   
        writer.writerow(["Repeat ID","Mutation ID", "Length Of Exact Repeat", "Length Of Approx Repeat"])
        
        
        for eachitem ,c in zip( dataList ,colors):
            if inverted == False :
                startIndex1, startIndex2, lapprox, mutationRate ,outputList =     findapproxrepeatLengthEditDistance(folderName+"/oneLine.fasta",folderName+"/oneLine2.fasta", eachitem[5], eachitem[6], eachitem[4])
            else:
                startIndex1, startIndex2, lapprox, mutationRate ,outputList =     findapproxrepeatLengthEditDistance(folderName+"/oneLine.fasta",folderName+"/reverseoneLine.fasta", eachitem[5], eachitem[6], eachitem[4])
            lrepeat = eachitem[4]
            
            plotY = []
            plotX = []
            
            for inneritem in outputList:
                plotY.append(inneritem[1]/float(lrepeat))
                plotX.append(inneritem[0]/float(lrepeat))
            
#            plotY2 = []
#            for elementindex in range(len(plotY)-1):
#                plotY2.append( (plotY[elementindex+1] - plotY[elementindex] )/ float(plotX[elementindex+1] - plotX[elementindex]))
            
            lengthOfRepeatList.append(str(eachitem[2])+" , "+ str(eachitem[4]))
            plt.plot(plotX, plotY,linewidth=1, marker = 'x', color = c,label= str(eachitem[2])+" , "+ str(eachitem[4]))
            

            
            for index in range(len(plotX)):        
              writer.writerow([counterindex, int(plotX[index]*lrepeat), lrepeat,int(plotY[index]*lrepeat)])
            

            counterindex = counterindex +1 
        
        ofile.close() 
        defaultsize = 24
        
        lims = plt.xlim()    
    
        plt.tick_params(axis='both', which='major', labelsize=defaultsize)
        plt.tick_params(axis='both', which='minor', labelsize=defaultsize)
    
       
        plt.title("Analysis of region beyond repeat in " +  folderName +" for top"+str(topX) +" Exact Repeats",fontsize=defaultsize)
        plt.xlabel("Cumulative Hamming distance / l_repeat ",fontsize=defaultsize)
        plt.ylabel("Stretch =  l_extended-repeat / l_repeat ",fontsize=defaultsize)
        handles, labels = ax.get_legend_handles_labels()

    # reverse the order
        ax.legend(handles[::-1], labels[::-1],prop={'size':8})

        #plt.legend(lengthOfRepeatList)

 
        if inverted == False:
            if plotindex  == 0:       
                plt.savefig("regionBeyondRepeat/FixScale_twoSides/"+folderName+"_approxRepeatAnalysisPlotEditDistance.png") 
            else:
                plt.savefig("regionBeyondRepeat/bestFit_twoSides/"+folderName+"_approxRepeatAnalysisPlotEditDistance.png") 

        else:
            if plotindex  == 0:       
                plt.savefig("regionBeyondRepeat/FixScale_twoSides/"+folderName+"_inverted_approxRepeatAnalysisPlotEditDistance.png") 
            else:
                plt.savefig("regionBeyondRepeat/bestFit_twoSides/"+folderName+"_inverted_approxRepeatAnalysisPlotEditDistance.png") 
  
  

  
### Viewpoint : Sliding window
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
            
            distance = distanceComputeLib.hammingDistance(str1, str2, windowSize)
            distanceList.append(distance)        
        
    f = open(outputResult,'w')
    for eachitem in distanceList:  
        f.write(str(eachitem) + "\n" )
    f.close()
    
    
    f1.close()
    f2.close()    


def extractApproxRepeatStatistics(folderName, topX):
    ioLib.generateOneLineFile(folderName, folderName+ "/genome.fasta")
    listOfData = ioLib.readFromMummerOutput(folderName+ "/long_repeats.txt")
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

        #reportPatternBeyondRepeat( folderName+ "/oneLine.fasta", folderName+ "/oneLine2.fasta", start1 - int( 1.5*length) ,start1 - int(  1.5*length) -1 , start2 - int( 1.5* length) ,start2 - int( 1.5* length) -1 , int(length*5/10),folderName+"/outputfile" + str(index)+".txt")
        #graphPlottingLib.readAndPlotStat(folderName+"/outputfile" + str(index)+".txt",folderName+"/outputfile" + str(index)+".txt2",folderName, index)
        graphPlottingLib.plotPatternBeyondRepeat( folderName+ "/oneLine.fasta",folderName+ "/oneLine2.fasta",  start1 - int( 1.5*length), start1 - int(  1.5*length) -1, start2 - int( 1.5* length) , start2 - int( 1.5* length) ,int(length*5/10))