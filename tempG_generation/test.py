#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

from math import *
from collections import *
from functools import *
from itertools import *
import datetime
import os

datasetPath = "/Users/siqingzhang/Desktop/dataset"
chosenFile = "wiki-talk-temporal.txt"

datafile = os.path.join(datasetPath, chosenFile)


f = open(datafile)
line = f.readlines()
allData = list()

# mark=0
for edge in line:
    # mark+=1
    allData.append(edge.split())
    # if mark>20:
    #     break

def my_compare(x,y):
    if x[2] > y[2]:
        return 1
    elif x[2] < y[2]:
        return -1
    elif x[2] == y[2]:
        if x[0] > y[0]:
            return 1
        elif x[0] < y[0]:
            return -1
        elif x[1] == y[1]:
            if x[1] > y[1]:
                return 1
            elif x[1] < y[1]:
                return -1
    return 0

# allData.sort(key=lambda x:x[2])

allData.sort(key=cmp_to_key(my_compare))


"""To show sorted info"""
# outputfile = chosenFile[:-4]+"-sorted"+".txt"
# fo = open(outputfile, "w+")
# for _ in allData:
#     line = " ".join(_)
#     line += "\n"
#     # print(line)
#     fo.write(line)
# fo.close()
#
"""To show sorted info with real-word time"""
# timefile = chosenFile[:-4]+"-time"+".txt"
# fd = open(timefile, "w+")
# for _ in allData:
#     realTime = datetime.datetime.fromtimestamp(int(_[2]))
#     _.append(str(realTime))
#
#     line = " ".join(_)
#     line += "\n"
#     # print(line)
#     fd.write(line)
# fd.close()

"""To show example"""
nVertices = 24818
exmaplefile = chosenFile[:-4]+"-"+str(nVertices)+"-example"+".txt"
fe = open(exmaplefile, "w+")
for _ in allData:
    if int(_[0]) == int(_[1]) or int(_[0]) > nVertices or int(_[1]) > nVertices:
        continue
    realTime = datetime.datetime.fromtimestamp(int(_[2]))
    _.append(str(realTime))
    line = " ".join(_)
    line += "\n"
    # print(line)
    fe.write(line)
fe.close()


