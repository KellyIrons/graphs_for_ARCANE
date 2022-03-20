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
2. Select the "Style" tab and under the top drop-down choose a pre-set style if desired. Node and edge display options can be changed in their respective tabs shown near the bottom left of the window.
      a. Use "continuous mapping" to vary different display options by parametere values, such as having edge with increse with increasing DIC or making    species with higher error propagation values darker in color.
3. Select "Layout" in the top toolbar to change the distribution of nodes. Some useful layouts are Attribute Circle Layout (by "EP") and Prefuse Force Directed Layout (by "weight"). You can also drag individual nodes to the desired location.
4. Select "Layout Tools" in the bottom left and drag the scale slider to make the nodes closer together or farther apart while preserving the same distribution.
5. To visualize the paths between individual species and the target species, select "Filter" --> "Column Filter" --> "Edge: path_name" --> "is". To use the filter, type the exact name of the species you want to see (as it is written on the node) and the edges in the path should be highlighted in red.

That is it! Cytoscape has lots of options to customize your graph, but following the steps above is a very good start.

For questions, please contact Kelly Irons at kei5@cornell.edu
