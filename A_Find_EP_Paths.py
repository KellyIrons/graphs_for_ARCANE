"""
ARCANE-friendly EP Path Finder
"""

def A_Find_EP_Paths(mechanism, targets, sdb, DIC_tol):
    
    
    """
    Finds the 'paths' taken by ARCANE when calculating the error propagation between
    the targe and each species

    :param mechanism: Mechanism object
    :param targets: list of reduction targets
    :param sdb: sample database on which DRGEP is applied
    :param DIC_tol: maximum difference between actual EP for a given species 
        and EP calculated from species path

    :return EP_paths: list of lists of species paths, by species index

    Created: 22/3/3 [KI]
    """
    
    # Import required packages
    import numpy as np
    import ARCANE.drgep as drgep 
    
    S = len(mechanism.species_names)
    net = mechanism.network
    
    # Find indices of targets 
    target_indices = []
    species_names = mechanism.species_names
    for target in targets:
        index = 0;
        while index <= (S-1) and target != species_names[index]:
            index = index + 1
        target_indices.append(index)
    
    # Initialize paths
    EP_paths = []
    
    # Calculate EP
    [EP, ind_list, DIC_spec, DIC_ind, coeff_list, alpha_norm] = drgep.graph_error_propagation(mechanism, sdb, targets)
    
    # Fins path for each species
    for i in range(S):
        
        #For multiple targets, would have to figure out which one is the one 
        # that the EP coefficient comes from. For now, just take the one target
        
        # Start out by adding the initial point (target) and the final point (species i)
        this_path = [target_indices[0], i]    
        
        
        # If the list isn't empty, then the target and the species have a path
        #  of more than one edge     
        if len(ind_list[i]) > 0:
            
            # Take the necessary data subsets
            DIC = DIC_spec[:, :, int(DIC_ind[i])]
            coeff_list_i = coeff_list[int(DIC_ind[i])]
            alpha_norm_loc = alpha_norm[0,int(DIC_ind[i])]
            
            # Initialize a list of possible paths
            paths = [this_path]
            
            
            # Take the last entry of ind_list as the second to last node, add it to the path
            node = ind_list[i][-1]
            paths[0].insert(1,node)
            
            final_EP = []
            
            flag = True
            # If the numer of iterations for this species in local_error_propagation was 1..
            if coeff_list_i[node,0] == -1:
                flag = False
                this_EP = DIC[target_indices[0], node]*DIC[node,i]
                flag2 = abs(EP[i] - this_EP*alpha_norm_loc) <= DIC_tol
                if flag2:
                    final_EP = paths[0]
                else:
                    print('Error! Try increasing DIC_tol.')
            else:
                count = 0 # need to count the number of levels we go through
                flag3 = False
                while flag and not flag3:
                    # Add a level
                    count = count+1
                    if not flag3:
                        for j in range(len(paths)):    # for each possible path          
                        
                            # For each level, look at the adjacent species as possible next paths
                            next_node = net.indsj[net.indsi == node].copy() 
                            
                            # shouldn't this be looking for species that are adjascent to the middle node?
                            for ind in paths[j]: 
                                if ind in next_node:
                                    del_ind = np.where(next_node == ind)[0][0]
                                    next_node = np.delete(next_node, del_ind)
                            
                            # Add potential next nodes to the path
                            for k in range(len(next_node)):
                                temp_path = paths[j].copy()
                                temp_path.insert(1, next_node[k])
                                
                                # Calculate EP
                                if coeff_list_i[next_node[k],count-1] == -1 and not flag3 :
                                    this_EP = 1
                                    for ind in range(len(temp_path)-1):
                                        this_EP = this_EP*DIC[temp_path[ind], temp_path[ind+1]]
                                    
                                    # Check if calculated EP is correct
                                    if abs(EP[i] - this_EP*alpha_norm_loc) <= DIC_tol:
                                        final_EP = temp_path.copy()
                                        flag3 = True
                                        continue
                                else:
                                    paths.append(temp_path)


        else:
            final_EP = this_path

        EP_paths.append(final_EP) 
     
        
    # Checking step - recalculate EP from paths
    for i in range(S):
        DIC_i = DIC_ind[i]
        DIC = DIC_spec[:, :, int(DIC_i)]
        path = EP_paths[i]
        
        coeff = 1
        for ind in range(len(path)-1):
            coeff = coeff*DIC[path[ind], path[ind+1]]
            
        alpha_norm_loc = alpha_norm[0,int(DIC_ind[i])]


        flag = abs(EP[i] - coeff*alpha_norm_loc) <= 1e-10
    
        if not flag:
            print('Error: Bad path for species', i)   
    
    
    return EP_paths

