import matplotlib.pyplot as plt
import os
import sys


def hammingDistance(str1, str2, k):
    count  =0
    for index in range(k):
        if str1[index] != str2[index]:
            count = count +1
    return count
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

def slidingWindowPlot(filepath, option,start1, start2 , plotRange, offset):

    myGenomeFile = open(filepath,'r')
    header = myGenomeFile.readline()
    myInfo = header.split('|')
    
    myname = myInfo[4].split(',')
    if myname[0][0] == ' ':
        
        folderName = myname[0][1:len(myname[0])]
    else:
        folderName = myname[0]
    
    newFolder = ""
    for eachChar in folderName: 
        if eachChar == ' ':
            newFolder= newFolder + "_"
        else:
            newFolder = newFolder + eachChar
            
        
    myGenomeFile.close()
    
    if option == "setup":
        os.system("mkdir "+newFolder)
        os.system("cp " + filepath + " " + newFolder+ "/genome.fasta")
        generateOneLineFile(newFolder, filepath)
 
    
    #start1, start2 , plotRange, offset =227625 , 4418733 ,10000,  2000
    plotPatternBeyondRepeat(newFolder+"/oneLine.fasta",newFolder+"/oneLine2.fasta",  start1-offset , start2-offset ,plotRange)


def plotPatternBeyondRepeat(genomeSource1,genomeSource2,  start1,  start2, plotRange):
    f1 = open(genomeSource1,'r')
    f2 = open(genomeSource2,'r')
    
    pointerLocation1 = start1
    pointerLocation2 = start2
    
    windowSize = 10
    
    defaultsize = 24
    distanceList = []
    plotRange = plotRange/windowSize
    
    for index in range(plotRange):
        f1.seek(pointerLocation1)
        f2.seek(pointerLocation2)
    
        str1 = f1.read(windowSize)
        str2 = f2.read(windowSize)
        
        pointerLocation1 = pointerLocation1 + windowSize
        pointerLocation2 = pointerLocation2 + windowSize
        
        distance = hammingDistance(str1, str2, min(windowSize,len(str1),len(str2)))
        distanceList.append(distance)        
        
    
    #plt.subplot(111)
    plt.plot(range(0,len(distanceList)*windowSize, windowSize), distanceList)
    
    plt.tick_params(axis='both', which='major', labelsize=defaultsize)
    plt.tick_params(axis='both', which='minor', labelsize=defaultsize)
    
    plt.xlabel("Genomics location of Copy 1 ( unit : base )", fontsize=defaultsize)    
    plt.ylabel("Hamming distance within a window of length 10 ",fontsize=defaultsize)
    
    #print plt.xticks()
    locs =range(0,len(distanceList)*windowSize+1, len(distanceList)*windowSize/4)
    tick_lbls =  range(start1, start1+plotRange*windowSize+1, plotRange*windowSize/4)
    plt.xticks(locs, tick_lbls, fontsize=defaultsize)
    
    ax2 = plt.twiny()
    plt.tick_params(axis='both', which='major', labelsize=defaultsize)
    plt.tick_params(axis='both', which='minor', labelsize=defaultsize)
    
    plt.xlabel("Genomics location of Copy 2 ( unit : base )", fontsize=defaultsize)    
    plt.ylabel("Hamming distance within a window of length 10 ",fontsize=defaultsize)
    
    tick_locs = range(0, plotRange*windowSize+1, plotRange*windowSize/4)
    tick_lbls = range(start2, start2+plotRange*windowSize+1, plotRange*windowSize/4)
    plt.xticks(tick_locs, tick_lbls,fontsize=defaultsize)
  
    plt.show()

    f1.close()
    f2.close()



def main(argv):
    slidingWindowPlot(argv[0],argv[1],int(argv[2]) ,int(argv[3]) ,int(argv[4]),  int(argv[5])) 
if __name__ == "__main__":
   main(sys.argv[1:])