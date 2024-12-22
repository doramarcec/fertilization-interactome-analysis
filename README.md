# Network analysis of the fertilization pathway

## Data accessibility

The data can be accessed on the [STRING](https://string-db.org/) database by searching for the fertilization term, selecting 'MAP-1187000 (Reactome)' from the list of terms, and choosing Homo sapiens as an organism. The network can be exported as tabular text output, which should download as the 'string_interactions.tsv' file. 

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
df
```

We will read the created data frame, which acts as an edge list, into a new variable. 
```
fertilization = nx.from_pandas_edgelist(df, source = "#node1", target = "node2")
```

From this variable, we can draw a simple graph. However, our understanding of this graph will be limited without understanding several basic concepts from Graph Theory. 
```
nx.draw_networkx(fertilization)
```

The graph theory states that:<br/>
>A *graph* consists of a set *V* of vertices and a set *E* of edges, each edge being associated with two vertices.

Vertices are also known as nodes, and they are the blue interconnected points on the graph. A line connecting the two nodes is called an edge. In the context of our fertilization pathway that we are about to analyse, nodes represent proteins involved in the pathway and edges represent the interactions between those proteins. 

The first and simplest step in analysing any network is to count the number of nodes and edges in a graph.
```
nx.number_of_nodes(fertilization) 
nx.number_of_edges(fertilization)
```

From this, we learn that our graph has 26 nodes and 124 edges, which tells us that our network consists of 26 proteins with 124 interactions. 

Let's identify the diameter of our network, which is defined as the "longest shortest path between any two nodes in the network." In the context of biological networks, this refers to the maximum number of interactions required for two proteins to indirectly affect one another.
```
nx.diameter(fertilization)
```

Networks with small diameters are usually highly interconnected, with efficient interaction between proteins, whereas networks with larger diameters often contain more isolated nodes representative of specialised proteins that require more intermediate nodes to connect with the rest of the network. 

The diameter of our network is 3





