import csv
import itertools
import random

nodeNum = 12
edgeNum = 24

with open("/Users/anakin/Downloads/data_test/test11.nodes.csv", "w") as csvFile:
    fieldNames = ["id", "x", "y"]
    writer = csv.DictWriter(csvFile, fieldnames=fieldNames)

    writer.writeheader()

    for i in range(0, nodeNum):
        writer.writerow({'id':str(i), 'x':0, 'y':0})


with open("/Users/anakin/Downloads/data_test/test11.edges.csv", "w") as csvFile:
    edges = list(itertools.permutations(range(0,nodeNum),2))
    edges_slice = random.sample(edges, edgeNum)
    fieldNames = ["source", "target"]
    writer = csv.DictWriter(csvFile, fieldnames=fieldNames)

    writer.writeheader()

    for i in range(0, edgeNum):
        writer.writerow({'source':str(edges_slice[i][0]), 'target':str(edges_slice[i][1])})
