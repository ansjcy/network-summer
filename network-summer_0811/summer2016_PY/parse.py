import csv

readerEdge = csv.reader(file("/Users/anakin/Downloads/data/karate.edges.csv", "rb"))
readerNode = csv.reader(file("/Users/anakin/Downloads/data/karate.nodes.csv", "rb"))



edges = []

for line in readerNode:
    if readerNode.line_num == 1:
        continue
    maxNode = line[0]

for line in readerEdge:
    if readerEdge.line_num == 1:
        continue
    edges.append([line[0], line[1]])

outputFile = open("/Users/anakin/Desktop/haha.json", "w")
outputFile.write("{\n\"nodes\":[\n")
for i in range(0, int(maxNode)):
    outputFile.write("{\"id\":\""+str(i)+"\", \"group\": 1},\n")
outputFile.write("{\"id\":\""+str(int(maxNode))+"\", \"group\": 1}\n],\n\"links\":[\n")
for i, edge in enumerate(edges):
    if i == len(edges) - 1:
        outputFile.write("{\"source\": \""+str(edge[0])+"\", \"target\":\""+str(edge[1])+"\", \"value\":1}\n]\n}")
    else:
        outputFile.write("{\"source\": \""+str(edge[0])+"\", \"target\":\""+str(edge[1])+"\", \"value\":1},\n")






