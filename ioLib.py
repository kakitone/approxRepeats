from operator import itemgetter

### MUMMER 
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

    listOfData = sorted(listOfData,key = itemgetter(2))
    return listOfData


def filterDataInverted(listOfData, G):
    listOfData = sorted(listOfData)
    oldData = []
    index =0
    while (index < len(listOfData)):
        if listOfData[index][0] >= G -  listOfData[index][1] - listOfData[index][2]: 
            listOfData.pop(index)
        elif listOfData[index] == oldData:
            listOfData.pop(index)
        else:
            oldData = listOfData[index]
            index = index +1

    listOfData = sorted(listOfData,key = itemgetter(2))
    return listOfData


        
def readFromMummerOutput(sourceFile, writeNew= False):
    f = open(sourceFile, 'r')
    
    if writeNew :
        f2 = open(sourceFile+"temp.txt", 'w')
        print "here"
    
    genomeName = f.readline()
    print genomeName
    testline = 'NOT NULL YET'
    listOfData = []
    while len(testline) > 0 :
        testline = f.readline()
        testList = testline.split()
        testList = map(int, testList)
        if len(testList) > 0 :
            
            if int(testList[2]) > 45:
                listOfData.append(testList)

            if writeNew and int(testList[2]) > 100:
                print testline
                f2.write(testline)
                
    f.close()
    
    if writeNew:
        f2.close()
    return listOfData           
            
def transformMummerOutput(listOfData):
    newListOfData = []
    
    for eachitem in listOfData :
        temp = eachitem
        temp.append(eachitem[2])
        newListOfData.append(temp)
    return newListOfData


def preprocessingRepeatInfo(folderName):

    dirList = [folderName+'/long_repeats.txt']
    dataList = []
    
    for fname in dirList:
        f= open(fname,'r')        
        f.readline()
        f.readline()
        temp = f.readline()
        while (len(temp) > 0 ):
            index1, index2, score = temp[0:len(temp)-1].split()
            
            index1ToAdd = int(index1)
            index2ToAdd = int(index2)
            scoreToAdd = int(score)
            if index1ToAdd <= index2ToAdd:
                dataList.append([fname,index1ToAdd,index2ToAdd,scoreToAdd])
                
            temp = f.readline()
        f.close()
    
    dataList = sorted(dataList, key= itemgetter(3), reverse = True)

    return dataList
    
def extractExactRepeats(folderName):    
    listOfData = readFromMummerOutput(folderName+ "/long_repeats.txt", True)
    listOfData= sorted(listOfData,key = itemgetter(2), reverse = True)
    
    return listOfData

### Inverted Repeat Treatment
def genomeInversion(folderName):
    fin = open(folderName  + "/oneLine.fasta", 'r')
    fout = open(folderName +"/reverseoneLine.fasta", 'w')

    temp = fin.read()   
    G = len(temp)
    
    fout.write(">Seq2 \n")
    
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

def reverseString(str1):
    str2 = ""    
    for index in range(len(str1)-1, -1, -1):
        str2 = str2 + str1[index]
    return str2
            
### File Formatting
def generateFilesInDifferentFormat(foldername):
    generateOneLineFile(foldername, foldername+"/genome.fasta" )
    genomeInversion(foldername)


def generateOneLineFile(foldername, filepath):
    f = open(filepath, 'r') 
    
    f2 = open(foldername +"/oneLine.fasta", 'w')  
    f3= open(foldername +"/oneLine2.fasta", 'w' )
    
    f4 = open(foldername + "/forward.fasta", 'w')
    
    temp = f.readline()
    
    while (len(temp) > 0):
        temp = f.readline()
        f2.write(temp[0:len(temp)-1])
        f3.write(temp[0:len(temp)-1])
        f4.write(temp[0:len(temp)-1])
    
    
    f.close()
    f2.close()
    f3.close()        
    f4.close()

    
    
### Analysis of lambda
def analysisOfLambda():
    f = open('lambda.fasta', 'r')
    
    
    for window in range(25, 10, -1):        
        listOfItems = []
        f.seek(0)
    
        for index in range(48502 - window  +1 ):
            f.seek(index)
            listOfItems.append( [f.read(window), index] )
     
        listOfItems = sorted(listOfItems)
                
        counter = 0
        for index in range(1, len(listOfItems)):
            if listOfItems[index][0] == listOfItems[index -1][0]:
                counter = counter +1 
                
        print "Window :"+ str(window) + "---Counter :" + str( counter)
                
    f.close()
    
#analysisOfLambda()
















    
    