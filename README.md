# Protein interaction network analysis in fertilization

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

The quickest way to import the protein interaction data is by exporting a .tsv file from the STRING database as described and reading it using pandas, as demonstrated below. Once imported, we will convert the data into the data frame and select columns of interest.
```
data = pd.read_csv("data/string_interactions.tsv", delimiter = "\t")

df = pd.DataFrame(data, columns = ['#node1', 'node2'])
df
```

We will read the created data frame, which acts as an edge list, into a new variable. 
```
fertilization = nx.from_pandas_edgelist(df, source = "#node1", target = "node2")
```

From this variable, we can draw a simple network. 
```
nx.draw_networkx(fertilization)
```






