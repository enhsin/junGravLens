import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
import os
import numpy as np

def readTextFile(fileList):
    concList = []
    for file in fileList:
        f = open(file, 'r')
        for line in f.readlines():
            temp = line.split()
            if(len(temp) > 1):
                concList.append(temp)
    return concList



def getFileList(path, prefix) :
    fileList = []
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    for file in onlyfiles:
        if file.startswith(prefix):
            fileList.append(file)
    return fileList


def findMin(finalList):
    maxObjFunc = [0, 0, 0, 1000000000]
    maxIndex = [0, 0, 0, 0]

    curIndex = 0
    for line in finalList:
        if(float(line[-3]) > maxObjFunc[0]):
            maxObjFunc[0] = float(line[-3])
            maxIndex[0] = curIndex
        if(float(line[-2]) > maxObjFunc[1]):
            maxObjFunc[1] = float(line[-2])
            maxIndex[1] = curIndex
        #if(float(line[-1]) > maxObjFunc[2]):
        #    maxObjFunc[2] = float(line[-1])
        #    maxIndex[2] = curIndex
        if(float(line[-1]) < maxObjFunc[3]):
            maxObjFunc[3] = float(line[-1])
            maxIndex[3] = curIndex


        curIndex += 1
    return maxIndex


def printModels(finalList, maxIndex):
    bestScatter = finalList[maxIndex[3]]
    outStr = "scatter: \t" + str(bestScatter[-1]) + "\n"
    for i in range(len(bestScatter)):
        outStr += (bestScatter[i] + "\t")
        if (i+1)%8 == 1:
            outStr  += "\n"
    # outStr = "Zeroth order: \t" + str(bestZeroth[-3]) + "\n"
    # for i in range(len(bestZeroth)):
    #     outStr += (bestZeroth[i] + "\t")
    #     if (i+1)%8 == 1:
    #         outStr  += "\n"
    # print outStr
    #
    # bestGrad = finalList[maxIndex[1]]
    # outStr = "Gradient order: \t" + str(bestGrad[-2]) + "\n"
    # for i in range(len(bestGrad)):
    #     outStr += (bestGrad[i] + "\t")
    #     if (i+1)%8 == 1:
    #         outStr  += "\n"
    # print outStr
    #
    # bestCubic = finalList[maxIndex[1]]
    # outStr = "Cubic order: \t" + str(bestCubic[-2]) + "\n"
    # for i in range(len(bestCubic)):
    #     outStr += (bestCubic[i] + "\t")
    #     if (i+1)%8 == 1:
    #         outStr  += "\n"
    print outStr

def main():
    prefix = "output."

    path = os.path.dirname(os.path.realpath(__file__))
    fileList = getFileList(path, prefix)
    finalList = readTextFile(fileList)
    maxIndex = findMin(finalList)
    printModels(finalList, maxIndex)


if __name__=='__main__':
    main()
