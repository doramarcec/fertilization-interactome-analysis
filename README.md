# Network analysis of the fertilization pathway

## Data accessibility

The data can be accessed on the [STRING](https://string-db.org/) database by searching for the term 'fertilization', selecting 'MAP-1187000 (Reactome)' from the list of terms, and choosing Homo sapiens as an organism. The network can be exported as tabular text output, which should download as the 'string_interactions.tsv' file. 

## Network analysis

First, we will install the NetworkX package in the terminal:
```
pip install networkx
```

Then, we will then import the required packages:
```
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
```

The quickest way to import the protein interaction data is to export a .tsv file from the STRING database as described and read it using pandas, as demonstrated below. Once imported, we will convert the data into the data frame and select columns of interest.
```
data = pd.read_csv("data/string_interactions.tsv", delimiter = "\t")

df = pd.DataFrame(data, columns = ['#node1', 'node2'])
```

We will read the created data frame, which acts as an edge list, into a new variable. 
```
fertilization = nx.from_pandas_edgelist(df, source = "#node1", target = "node2")
```

From this variable, we can draw a simple graph. 
```
nx.draw_networkx(fertilization)
```

The above line of code plots the following graph:<br/>

![fertilization-network](https://github.com/user-attachments/assets/77ffe8e6-78e4-43bb-b1c5-7c346314ecdc)

From this, we can see that our network is divided into two groups of proteins that are loosely connected. In other words, our network has a modular structure which is very common for biological networks, and the two modules encode different biological functions. To understand this better and take this analysis further, we will introduce a few concepts from the Graph Theory.

The Graph Theory states that:<br/>
>A *graph* consists of a set *V* of vertices and a set *E* of edges, each edge being associated with two vertices.

Vertices are also known as nodes, and they are the blue interconnected points on the graph, representing different proteins. A line connecting the two nodes is called an edge, representing an interaction between the two proteins.

The first and simplest step in analysing any network is to count the number of nodes and edges in a graph.
```
nx.number_of_nodes(fertilization) 
nx.number_of_edges(fertilization)
```

From this, we learn that our graph has 26 nodes and 124 edges, which tells us that our network consists of 26 proteins with 124 interactions. 

Next we will identify the connectivity, also known as the degree of all nodes, which represents a number of links (i.e. edges) a node has to other nodes. 
```
nx.degree(fertilization)
```

From this, we get the following output, and learn that ADAM2 node has the highest degree of 17 links with other nodes.<br/>
>DegreeView({'ACR': 11, 'SPAM1': 12, 'ADAM2': 17, 'ZP1': 14, 'IZUMO2': 6, 'IZUMO3': 9, 'CD9': 7, 'ZP4': 13, 'ZP3': 13, 'IZUMO4': 6, 'IZUMO1': 10, 'ZP2': 12, 'CATSPER1': 11, 'OVGP1': 10, 'B4GALT1': 10, 'ADAM30': 10, 'ADAM21': 10, 'ADAM20': 10, 'CATSPERB': 7, 'CATSPER3': 8, 'HVCN1': 6, 'KCNU1': 6, 'CATSPER2': 8, 'CATSPER4': 8, 'CATSPERD': 8, 'CATSPERG': 6})

We can visualise this information in a degree distribution:
```
plt.plot(nx.degree_histogram(fertilization))
```

This plots the histogram below, where X-axis represents a degree, and y-axis represents a frequency of a degree:<br/>

![degree-histogram](https://github.com/user-attachments/assets/469278b1-77df-4041-a25b-f3f4788f56f3)

From this histogram, we can understand that the minimum and maximum number of links a node in our network has is 6 and 17, respectively, including 5 nodes containing 6 links (in other words, 5 nodes with a degree of 6), and 1 node with a degree of 17. Moreover, 10 is the most frequent degree, as 6 nodes have a degree of 10.

Nodes with high degrees are known as hubs, which are essential nodes in a network that often have a crucial role in information-transfer. As we are working with quite a small network, we easily identified the hub from the above output. However, imagine we had a much larger network that consists of hundreds of nodes... it would take too much time to manually search for the hub. In that case, this is how we could go about identifying our hub:
```
degrees = dict(nx.degree(fertilization))
hub = max(degrees, key = degrees.get)
print("This network contains a hub", hub, "with a degree of", degrees[hub])
```

This prints the following output:<br/>
>This network contains a hub ADAM2 with a degree of 17

From the analysis we have done so far, we understand that our network of 26 proteins has two modules representing two biological functions with a total of 124 interactions. We identified the most connected protein in the network, which is ADAM2 with links to 17 other proteins, which forms an essential node with a likely role in information-transfer between the two modules. 

Now that we have identified the degree of our nodes and the hub, we will assess the closeness and betweenness centrality. Closeness centrality can simply be defined as a measure of the distance between a node and all other nodes of the network. We can calculate it using the code below:
```
nx.closeness_centrality(fertilization)
```
And get the output below:
>{'ACR': 0.5434782608695652,<br/>
>'SPAM1': 0.5434782608695652,<br/>
>'ADAM2': 0.7575757575757576,<br/>
>'ZP1': 0.6944444444444444,<br/>
>'IZUMO2': 0.49019607843137253,<br/>
>'IZUMO3': 0.5952380952380952,<br/>
>'CD9': 0.5,<br/>
>'ZP4': 0.6756756756756757,<br/>
>'ZP3': 0.5681818181818182,<br/>
>'IZUMO4': 0.49019607843137253,<br/>
>'IZUMO1': 0.5319148936170213,<br/>
>'ZP2': 0.5434782608695652,<br/>
>'CATSPER1': 0.6410256410256411,<br/>
>'OVGP1': 0.5208333333333334,<br/>
>'B4GALT1': 0.5208333333333334,<br/>
>'ADAM30': 0.5208333333333334,<br/>
>'ADAM21': 0.5208333333333334,<br/>
>'ADAM20': 0.5208333333333334,<br/>
>'CATSPERB': 0.44642857142857145,<br/>
>'CATSPER3': 0.5102040816326531,<br/>
>'HVCN1': 0.43859649122807015,<br/>
>'KCNU1': 0.43103448275862066,<br/>
>'CATSPER2': 0.45454545454545453,<br/>
>'CATSPER4': 0.45454545454545453,<br/>
>'CATSPERD': 0.45454545454545453,<br/>
>'CATSPERG': 0.43859649122807015}<br/>

This output tells us that ADAM2 node has the highest closeness centrality of 0.7576, followed by ZP1 (0.6944) and ZP4 (0.6757). Nodes with a high closeness centrality tend to be highly influential within a network. 

Betweenness centrality, on the other hand, measures the number of shortest paths between all pairs of nodes that pass through that node, and can be calculated as follows:
```
nx.betweenness_centrality(fertilization)
```
Which gives us the following betweenness centrality measures for each node:<br/>
>{'ACR': 0.01630952380952381,<br/>
>'SPAM1': 0.014111111111111116,<br/>
>'ADAM2': 0.21144841269841266,<br/>
>'ZP1': 0.10678174603174607,<br/>
>'IZUMO2': 0.0,<br/>
>'IZUMO3': 0.0855,<br/>
>'CD9': 0.003781746031746032,<br/>
>'ZP4': 0.09408730158730162,<br/>
>'ZP3': 0.017753968253968255,<br/>
>'IZUMO4': 0.0,<br/>
>'IZUMO1': 0.012527777777777778,<br/>
>'ZP2': 0.006587301587301587,<br/>
>'CATSPER1': 0.3495277777777778,<br/>
>'OVGP1': 0.0,<br/>
>'B4GALT1': 0.0,<br/>
>'ADAM30': 0.0,<br/>
>'ADAM21': 0.0,<br/>
>'ADAM20': 0.0,<br/>
>'CATSPERB': 0.0016388888888888892,<br/>
>'CATSPER3': 0.06983333333333334,<br/>
>'HVCN1': 0.0009722222222222222,<br/>
>'KCNU1': 0.0005555555555555556,<br/>
>'CATSPER2': 0.0028611111111111116,<br/>
>'CATSPER4': 0.0028611111111111116,<br/>
>'CATSPERD': 0.0028611111111111116,<br/>
>'CATSPERG': 0.0}<br/>

From this output, we can see that CATSPER1 node has the highest betweenness centrality with a value of 0.34953, followed by ADAM2 with a value of 0.21145. Nodes with the highest betweenness centrality tend to represent bridges or bottlenecks for information-transfer through the network.

Based on that, we can assume that the interaction between the ADAM2 protein and CATSPER1 protein connect the two biological functions in a network. As our network has a modular structure, it may be worth looking deeper into individual modules. To do so, we will explore the concept of cliques. 

A clique is a complete subgraph of a graph in which every node (i.e. protein) has edges connecting it to every other node that is a part of the defined subgraph, in other words, it represents a subset of proteins in which every single protein interacts with each other. 

To count the number of cliques every node is a part of, we can take the following approach:
```
print(nx.number_of_cliques(fertilization))
```

Which gives us the following output:<br/>
>{'ACR': 5, 'SPAM1': 3, 'ADAM2': 7, 'ZP1': 5, 'IZUMO2': 1, 'IZUMO3': 3, 'CD9': 2, 'ZP4': 4, 'ZP3': 4, 'IZUMO4': 1, 'IZUMO1': 3, 'ZP2': 3, 'CATSPER1': 5, 'OVGP1': 1, 'B4GALT1': 1, 'ADAM30': 1, 'ADAM21': 1, 'ADAM20': 1, 'CATSPERB': 2, 'CATSPER3': 3, 'HVCN1': 2, 'KCNU1': 2, 'CATSPER2': 4, 'CATSPER4': 4, 'CATSPERD': 4, 'CATSPERG': 1}

From this, we deduct that ADAM2 is a part of 7 cliques, and CATSPER1 is a part of 5 cliques, just like ZP1, making these three proteins highly central and explaining their high scores on closeness and betweenness centrality. 

However, we can get more information about the cliques, by finding the largest clique and identifying all of its nodes, as below:
```
print(max(nx.find_cliques(fertilization), key = len))
```

From that, we get the list of nodes:<br/>
>['ADAM2', 'ZP1', 'ZP3', 'ZP2', 'ZP4', 'SPAM1', 'ADAM30', 'ADAM21', 'ADAM20', 'B4GALT1', 'OVGP1']

This tells us that aside from being a hub, ADAM2 forms a clique with 10 other proteins, with two clearly defined protein families, ZP1, ZP2, ZP3, and ZP4, and ADAM2, ADAM20, ADAM21, and ADAM30. 

Finally, to visualise all of our findings more effectively, we can plot our network by colouring the nodes by their degree and altering their size by their betweenness centrality. 
```
pos = nx.spring_layout(fertilization)
betCent = nx.betweenness_centrality(fertilization, normalized=True, endpoints=True)
node_color = [20000.0 * fertilization.degree(v) for v in fertilization]
node_size =  [v * 10000 for v in betCent.values()]
plt.figure(figsize=(20,20))
nx.draw_networkx(fertilization, pos=pos, with_labels=True,
                 node_color=node_color,
                 node_size=node_size)
```

From the above code chunk, we get the following output:

![coloured-network-graph](https://github.com/user-attachments/assets/6abe66fe-2efc-4890-b92d-ac0c95ed6ff9)

From this, we can see the module separation more clearly, and we can see that actually there may be not two, but three modules. 

## Analysis in the context of literature

ZP1, ZP2, ZP3, and ZP4 are glycoproteins composing zona pellucida, which is an oocyte-surrounding extracellular matrix that mediates sperm binding. ZP1 is specifically responsible for maintaining the structural integrity of zona pellucida, ZP3 is involved in the formation of the zona matrix and binding of the sperm, whereas ZP2 and ZP4 do not have clearly defined roles but they are speculated to act as sperm receptors. 

On the other hand, ADAM2, ADAM20, ADAM21, and ADAM30 are domains of a sperm membrane glycoprotein known as fertilin, which has an important role in migration of the spermatozoa through the oviduct, and binding both the zona pellucida and eventually the oocyte. 

CATSPER1, with the greatest betweenness centrality, is a subunit of a Catsper channel in spermatozoa, which is a Ca<sup>2+</sup> channel that mediates the motility of the spermatozoa. 



