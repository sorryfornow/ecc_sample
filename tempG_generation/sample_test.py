#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

from math import *
from collections import *
from functools import *
from itertools import *
import datetime
import os

chosenFile = "wiki-talk-temporal.txt"
nVertices = 24818   # edges outcome: 1325803
datafile = chosenFile[:-4]+"-"+str(nVertices)+"-example"+".txt"


f = open(datafile)
line = f.readlines()
allData = list()
for edge in line:
    allData.append(edge.split())


Dge = defaultdict(list)
timeStamp = dict()
timeCount = 1
arrayEdges = set()

""" edge[3] represent the date"""
for edge in allData:
    Dge[edge[0]] += (edge[1], edge[3])
    Dge[edge[1]] += (edge[0], edge[3])
    if edge[3] not in timeStamp:
        timeStamp[edge[3]] = timeCount
        timeCount += 1

    if int(edge[0]) < int(edge[1]):
        arrayEdges.add((edge[0],edge[1],str(timeStamp[edge[3]])))
    else:
        arrayEdges.add((edge[1], edge[0], str(timeStamp[edge[3]])))

print("time interval:", 1, "-", timeCount-1)
c = []
for i in Dge.keys():
    c.append(int(i))
c.sort()
print("vertices:", c)

arrayEdges = list(arrayEdges)
arrayEdges.sort(key=lambda x: int(x[2]))
print("number of edges:", len(arrayEdges))
# print(arrayEdges)


outputfile = chosenFile[:-4]+"-fine-"+str(nVertices)+"-example"+".txt"
fo = open(outputfile, "w+")
for _ in arrayEdges:
    line = " ".join(_)
    line += "\n"
    # print(line)
    fo.write(line)
fo.close()



