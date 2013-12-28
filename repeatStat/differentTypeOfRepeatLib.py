from operator import itemgetter
import matplotlib.pyplot as plt
import math
import numpy as np
import os
import csv
import bisect
import random
#import random

import approximateRepeatLib
### Extract Triple/interleave Repeat

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
    
    '''
    print "metaListForInterleaveRepeat"
    for index in range(len(metaListForInterleaveRepeat)):
        tempdata = metaListForInterleaveRepeat[index] 
        print tempdata
    print "metaListForTripleRepeat"
    for index in  range(len(metaListForTripleRepeat)):
        print metaListForTripleRepeat[index]
    '''
      
    return metaListForInterleaveRepeat[0:min(len(metaListForInterleaveRepeat),topX)], metaListForTripleRepeat[0:min(len(metaListForTripleRepeat),topX)]


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

                if ( x1 < x2 and x2 < y1 and y1  <y2)  or (x2 < x1 and x1  < y2 and y2 < y1 ):                    
                    interleaveRepeatList.append([filename,x1, y1, repeat1, x2, y2, repeat2, min(repeat1, repeat2) ])
        
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

                if (x2 >= y1 and x2 + repeat2 <=y1 + repeat1 and  (x2 != y1 or y2!= x1))  or (x2 >= x1 and x2 + repeat2 <= x1+ repeat1 and  (x2 != x1 or y2 != y1)):
                    tripleRepeatList.append([filename,x1, y1, repeat1, x2, y2, repeat2, min(repeat1, repeat2) ])
        
        return tripleRepeatList


def extractShorterCopy(maxInterleaveRepeatLength):
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
    listOfInterleaveRepeat = approximateRepeatLib.filterSameApproxRepeat(listOfInterleaveRepeat)
    
    listOfInterleaveRepeat = sorted(listOfInterleaveRepeat,key = itemgetter(2) ,reverse = True)
    
    return listOfInterleaveRepeat
            

  
