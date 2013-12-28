from operator import itemgetter

import approximateRepeatLib
import differentTypeOfRepeatLib
import ioLib



def combinedStatisticsPlot(folderName,topX,typeOferror):
    dataList = ioLib.preprocessingRepeatInfo(folderName)

    metaListForInterleaveRepeat, metaListForTripleRepeat = differentTypeOfRepeatLib.findDifferentTypesOfRepeat(dataList,topX)
       
    maxInterleaveRepeatLength = metaListForInterleaveRepeat[0][-1]
    maxTripleRepeatLength = metaListForTripleRepeat[0][-1]
    
    print "maxInterleaveRepeatLength,maxTripleRepeatLength", maxInterleaveRepeatLength,maxTripleRepeatLength  

    listOfRepeat = ioLib.extractExactRepeats(folderName)
    print listOfRepeat
    listOfApproxRepeat = approximateRepeatLib.clusteringRepeatMain(folderName, listOfRepeat, topX, 'r', typeOferror , -1, [])


    for eachitem in listOfApproxRepeat:
        #startIndex1, startIndex2, lapprox, mutationRate,lengthOfExactRepeat, start1, start2
        if eachitem[2] < eachitem[4]:
            print "Errorrrrrrrrrrrrrrrr", eachitem[2] , eachitem[4]
        eachitem.insert(0,'eee')
     
    print "listOfApproxRepeat", listOfApproxRepeat
    metaListForInterleaveRepeat, metaListForTripleRepeat = differentTypeOfRepeatLib.findDifferentTypesOfRepeat(listOfApproxRepeat,topX)    
    print "metaListForInterleaveRepeat, metaListForTripleRepeat",metaListForInterleaveRepeat, metaListForTripleRepeat
    
    listOfInterleaveRepeat = differentTypeOfRepeatLib.extractShorterCopy(metaListForInterleaveRepeat)
    
    if len(listOfInterleaveRepeat) > 0:
        approximateRepeatLib.clusteringRepeatMain(folderName, listOfInterleaveRepeat, topX, 'i',typeOferror, maxInterleaveRepeatLength,listOfApproxRepeat)

    listOfTripleRepeat = differentTypeOfRepeatLib.extractShorterCopy(metaListForTripleRepeat)
    
    if len(listOfTripleRepeat) > 0:
        approximateRepeatLib.clusteringRepeatMain(folderName, listOfTripleRepeat, topX, 't',typeOferror,maxTripleRepeatLength,listOfApproxRepeat)

    listOfRepeat = ioLib.readFromMummerOutput(folderName+ "/long_repeats_inverted.txt")
    print listOfRepeat
    
    fForG = open(folderName+ "/oneLine.fasta", 'r')
    
    G = len(fForG.read())
    
    fForG.close()
    listOfRepeat = ioLib.filterDataInverted(listOfRepeat, G)
    print listOfRepeat
    listOfRepeat= sorted(listOfRepeat,key = itemgetter(2), reverse = True)
    
    
    
    print "Analyzing inverted repeats", listOfRepeat
    if len(listOfRepeat) > 0:
        listOfInvertedRepeat = approximateRepeatLib.clusteringRepeatMain(folderName, listOfRepeat, topX, 'c', typeOferror , -1, [])   

folders = ["LactoAcid","EcoliK12","Ecoli536","Salmonella","Pmarinus" ]    

#folders = ["EcoliK12"]
for item in folders:
    combinedStatisticsPlot(item,100,'h')
    combinedStatisticsPlot(item,100,'e')



