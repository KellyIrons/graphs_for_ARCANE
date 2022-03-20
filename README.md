# Overview

This set of scripts interfaces with ARCANE to add the ability to nicely visualize chemical mechanisms. The main script is A_Graph_Creation.py, which calls graph-specific functions in ARCANE's drgep module to write a .graphml file. A_Find_EP_paths.py is called within the main script to determine the "path" between each species and the target species that is taken when calculating the error propagation coefficient. 


# Workflow 

1. Download the three files mentioned above. The repository also includes an example .graphml file but that is not needed elsewhere.
2. Change lines 23-26 of A_Graph_Creation.py to the parameters you want to govern your graph.
3. Change lines 29-40 of A_Graph_Creation.py to align with the mechanism you want to analyze.
4. Run (the script)!
5. Visualize the mechanism using the method outlined in the next section.


# Visualization

The code creates a .graphml file using the NetworkX python package. The graph object used to create the file contains a lot of information on its own (see https://networkx.org/ if you are interested in doing more graph theory style analysis) but visualization is an easy way to get a good sense of the mechanism with little additional effort.

Multiple tools can be used to visualize .graphml files, but the preferred one by the author is Cytoscape because it is very user-friendly and makes sophisticated looking graphs. Cytoscape can be downloaded for free https://cytoscape.org/ and the following process can be used to visualize your file:

1. Select the "Import Network From File System" button (the arrow with the little graph/molecule) and select your file.
