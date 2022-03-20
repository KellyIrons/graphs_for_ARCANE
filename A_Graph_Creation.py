"""
ARCANE-friendly Graph Creation
"""


''' # Things for Kelly to run code
import sys
sys.path.append('/Users/kirons/Desktop/RDX_Project/Python_Files/arcane-master/')
import ARCANE
'''

# Import required packages
import numpy as np
import networkx as nx
import ARCANE.drgep as drgep
import ARCANE.cases as cases
import ARCANE.mechanisms as mechanisms
import ARCANE.sampling as sampling
from A_Find_EP_Paths import A_Find_EP_Paths

# Graph-specific user inputs
graph_time = 0.5  # time at which to calculate DICs for the graph
DIC_cutoff = 0.25 # minimum DIC to be included on graph
DIC_tol = 1E-15 # maximum difference between actual EP for a given species and EP
                # calculated from species path
name = 'gri211_w_paths' # name of the graph file to be output

# Set up the mechanism to be drawn (pulled from other demos to test code)
cti_init = "cti/gri211.cti"
graph_mechanism = mechanisms.Mechanism(cti_init, name='GRI-Mech211')
air = "X/O2/0.21/N2/0.79"
fuel = "X/CH4/1"
case_targets = ["CH4"]
temperature_range_0D = "1000"
temperature_range_1D = "300"
pressure_range = "1e5"
phi_range = "1"
targets = case_targets
caselist = []
caselist.extend(cases.create_case(reactor="0DIsochoric",
                                  mechanism=graph_mechanism,
                                  fuel=fuel,
                                  oxidizer=air,
                                  pressure=pressure_range,
                                  temperature=temperature_range_0D,
                                  phi=phi_range,
                                  targets=case_targets,
                                  error_dict={'tig': 0.05}))


# Run Cases
cases.run_cases(caselist, graph_mechanism, overwrite=True)
sdb = sampling.samples_database(caselist, graph_mechanism)

############## Create the graph!  ##########################################
reduction_type = 'species'

grids = []
for samples in sdb:
    grids.append(samples.grid)

# Find the sample corresponding to the desired time
if graph_time > sdb[-1].grid or graph_time < sdb[0].grid:
    print('Time out of range, middle time chosen')
    graph_sample = sdb[int(np.floor(len(sdb)/2))]
else:
    array = np.asarray(grids)
    idx = (np.abs(array - graph_time)).argmin()
    graph_sample = sdb[idx]

# Compute DIC and EP
DIC_spec, DIC_reac = drgep.compute_DIC(graph_mechanism, graph_sample, reduction_type)
EP = drgep.error_propagation(graph_mechanism, sdb, targets, reduction_type)

# Build the graph
G = nx.MultiDiGraph()

Slist = graph_mechanism.species_names
J = len(Slist)
added = [False for i in range(len(Slist))]

# For each pait of species, add either species as a node if it hasn't already 
# been added, and add an edge bewteen them if the DIC is above DIC_cutoff 
for i in range(J):
    for j in range(J):
            if i!=j and DIC_spec[i,j] > DIC_cutoff:
                if not added[i]:
                    G.add_nodes_from([(i, {"name": Slist[i], 'EP':EP[i]})])
                    added[i] = True
                if not added[j]:
                    G.add_nodes_from([(j, {"name": Slist[j],'EP':EP[j]})])
                    added[j] = True
                if added[i] and added[j]:
                    G.add_edge(i, j, weight=DIC_spec[i,j]) # arrow from j to i

# Find the path from each each species to the target
EP_paths = A_Find_EP_Paths(graph_mechanism, targets, sdb, DIC_tol)

# Include paths of all species that are in the graph
for i in range(len(graph_mechanism.species_names)):
    if added[i]:
        nx.add_path(G, EP_paths[i], path_name = graph_mechanism.species_names[i], path_num = i)

# Write the graph to a .graphml file
nx.write_graphml(G, name + '.graphml')









