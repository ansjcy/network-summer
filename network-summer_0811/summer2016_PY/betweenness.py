import networkx as nx
import csv

readerEdge = csv.reader(file("/Users/anakin/Downloads/data_test/adjnoun.edges.csv", "rb"))
readerNode = csv.reader(file("/Users/anakin/Downloads/data_test/adjnoun.nodes.csv", "rb"))

# readerEdge = csv.reader(file("/Users/anakin/Downloads/data_test/test12.edges.csv", "rb"))
# readerNode = csv.reader(file("/Users/anakin/Downloads/data_test/test12.nodes.csv", "rb"))

# readerEdge = csv.reader(file("/Users/anakin/Downloads/data_test/karate.edges.csv", "rb"))
# readerNode = csv.reader(file("/Users/anakin/Downloads/data_test/karate.nodes.csv", "rb"))


G = nx.DiGraph()

# rowNum = sum(1 for row in readerNode)
#
# for line in readerNode:
#     if readerNode.line_num == 1:
#         continue
#     maxNode = line[0]

# line = readerNode.get

# G.add_nodes_from([0, 7])
# G.add_edges_from([(2,3),(7,3),(7,1),(2,1),(7,2),(5,2),(3,4),(5,6),(1,5)])

# G.add_nodes_from([0, 3])
# G.add_edges_from([(0,2),(1,2),(1,3),(0,1)])

for line in readerEdge:
    if readerEdge.line_num == 1:
         continue
    G.add_edge(int(line[0]), int(line[1]), weight=2.0)
    # G.add_edge(int(line[1]), int(line[0]), weight=2.0)

# G.add_edge(1,0,weight=2)
# G.add_edge(10,11,weight=2)
# G.add_edge(3,0,weight=1)
G[1][0]['weight'] = 1
G[11][10]['weight'] = 1
# G[0][3]['weight'] = 1.66667
# G[0][5]['weight'] = 1.66667
# G[1][0]['weight'] = 1.66667
# G[3][0]['weight'] = 1.0
# G[5][0]['weight'] = 1.5

betweenness = nx.betweenness_centrality(G, normalized=False, weight='weight')
for i in betweenness:
    print str(i) + ': '+str(betweenness[i])





