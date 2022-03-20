"""Module implementing DRGEP for species and reactions reduction"""

import ARCANE.sampling as sampling
import ARCANE.display as display
import ARCANE.kwdict as kwdict
import ARCANE.analysis as analysis

import numpy as np

import sys


logger = display.Logger()
logger.set_log('logReduction')
kwdict = kwdict.Kwdict()


def drgep_ranked(caselist, mechanism, reduction_type, integrity=True, sensitivity=0):
    """Sort species using DRGEP method
    Sorts species using DRGEP method
    DRGEP is applied to all samples in sdb, uses mechanism in mechanism

    :param caselist: list of Case objects
    :param mechanism: Mechanism object
    :param reduction_type: type of reduction (species/S or reactions/R)
    :param integrity: check if there is no mass sink or source
    :param sensitivity: Activate aiding DRGEP with SA (Sensitivity analysis)
    (0 to be not activated, a value to be the sensor on the DRGEP coefficients. See sensitivity_analysis in analysis.)
    :return: array containing ranked species

    Created: 17/11/14 [PP]
    Last modified: 19/08/30 [JW]
    """

    if reduction_type not in ['species', 'S', 'reactions', 'R']:
        logger.error(" The type of reduction must be either 'species' (or 'S') or  'reactions' (or 'R')")
        sys.exit()

    ctmech = mechanism.ctmech
    species_names = ctmech.species_names
    reactions = ctmech.reactions()

    myns = mechanism.ns
    mynr = mechanism.nr

    if reduction_type in ['species', 'S']:
        EP = np.zeros(myns, 'd')
    else:
        EP = np.zeros(mynr, 'd')

    # Putting together cases with the same targets
    cases_dict = {}
    for case in caselist:
        targets = repr(case.targets)
        if targets not in list(cases_dict.keys()):
            cases_dict[targets] = [case]
        else:
            cases_dict[targets].append(case)

    cases_tuples = list(zip(cases_dict.keys(), cases_dict.values()))

    for targets, cases in cases_tuples:

        targets = eval(targets)

        # Create sample database on cases
        sdb = sampling.samples_database(cases, mechanism)

        EPtmp = error_propagation(mechanism, sdb, targets, reduction_type)
        EP = np.fmax(EP, EPtmp)

    if integrity and reduction_type in ['species', 'S']:
        no_more_consumed, no_more_produced, remove_if_removed = integrity_check(mechanism)
        # Integrity check
        if integrity:
            for spec_removed in no_more_produced:
                for spec_impacted in no_more_produced[spec_removed]:
                    EP[ctmech.species_index(spec_impacted)] = EP[ctmech.species_index(spec_removed)]

            for spec_removed in no_more_consumed:
                for spec_impacted in no_more_consumed[spec_removed]:
                    max_EP_spec_impacted = [EP[ctmech.species_index(spec_impacted)] \
                                            for spec_impacted in no_more_consumed[spec_removed]]
                    max_value = np.max([EP[ctmech.species_index(spec_removed)],
                                        np.max(max_EP_spec_impacted)])
                    EP[ctmech.species_index(spec_removed)] = max_value
                    EP[ctmech.species_index(spec_impacted)] = max_value

    if reduction_type in ['species', 'S']:
        # Keep important species
        important_species = important_species_identification(caselist, mechanism)

        # Add sensitivity to the important species
        # Species are switched according to an arbitrary coefficient defined in sensitivity
        # => If the species is a target, the DRGEP coefficient will be 1.
        # => If the species is under the sensitivity coefficient, treatment is done to keep the DRGEP in this order.
        # => If the species is above the sensitivity coefficient, treatment is done to put the species coefficient
        # bigger than the DRGEP.
        if sensitivity:
            logger.info('Sensitivity applied for species coefficient above ' + str(sensitivity))
            for id_case, case in enumerate(caselist):
                result_analysis = analysis.sensitivity_analysis(case, reduction_type)[case.mechanism.name][case.myid]
                for id_err, err in enumerate(list(case.error_dict.keys())):
                    if id_case + id_err == 0:
                        sa_species = result_analysis[err]
                    else:
                        sa_species = {key_dic: max(value_dic, result_analysis[err][key_dic])\
                                      for key_dic, value_dic in sa_species.items()}
            for id_spec, spec in enumerate(species_names):
                if spec in important_species:
                    EP[ctmech.species_index(spec)] = 1
                elif sa_species[id_spec] < sensitivity:
                    EP[ctmech.species_index(spec)] = EP[ctmech.species_index(spec)]*max(min(sa_species.values()), 1e-30)\
                                                     /max(sa_species.values())
                else:
                    EP[ctmech.species_index(spec)] = sa_species[id_spec]/max(sa_species.values())

                #logger.info('Sensible species :' + ', '.join([sa_species[spec] for spec in species_names if sa_species[spec] < sensitivity]))
        # Classical method
        else:
            EP = [1.0 if spec in important_species else EP[ctmech.species_index(spec)] for spec in species_names]

        drgep_dict = dict(zip(species_names, EP))

    elif reduction_type in ['reactions', 'R']:

        # Add sensitivity to the important reactions
        # Species are switched according to an arbitrary coefficient defined in sensitivity
        # => If the reaction is a target, the DRGEP coefficient will be 1.
        # => If the reaction is under the sensitivity coefficient, treatment is done to keep the DRGEP in this order.
        # => If the reaction is above the sensitivity coefficient, treatment is done to put the species coefficient
        # bigger than the DRGEP.

        if sensitivity:
            logger.info('Sensitivity applied for reactions coefficient above ' + str(sensitivity))
            for id_case, case in enumerate(caselist):
                result_analysis = analysis.sensitivity_analysis(case, reduction_type)[case.mechanism.name][case.myid]
                for id_err, err in enumerate(list(case.error_dict.keys())):
                    if id_case + id_err == 0:
                        sa_reactions = result_analysis[err]
                    else:
                        sa_reactions = {key_dic: max(value_dic, result_analysis[err][key_dic])\
                                      for key_dic, value_dic in sa_reactions.items()}
            for id_reac, reac in enumerate(reactions):
                if sa_reactions[reac.equation] < sensitivity:
                    EP[id_reac] = EP[id_reac]*max(min(sa_reactions.values()), 1e-30)\
                                                     /max(sa_reactions.values())
                else:
                    EP[id_reac] = sa_reactions[reac.equation]/max(sa_reactions.values())

        drgep_dict = dict(zip(reactions, EP))

    else:
        drgep_dict = {}
        logger.error('DRGEP is not applicable on this case !')

    return drgep_dict


def error_propagation(mechanism, sdb, targets, reduction_type):
    """
    Computes error propagation coefficients

    :param mechanism: Mechanism object
    :param sdb: sample database on which DRGEP is applied
    :param targets: list of reduction targets
    :param reduction_type: type of reduction performed (species or reactions)

    :return: dictionary of species or reactions with their corresponding coefficient

    Created: 17/11/14 [PP]
    Last modified: 18/04/03 [QC]
    """

    myns = mechanism.ns
    mynr = mechanism.nr

    if reduction_type in ['species', 'S']:
        EP = np.zeros(myns, 'd')
    elif reduction_type in ['reactions', 'R']:
        EP = np.zeros(mynr, 'd')
    else:
        EP = []
        logger.error('Wrong reduction type !')

    alpha_norm = scaling_coefficient(mechanism, sdb, targets)

    logger.terminator('\r')
    for i_sample, mysample in enumerate(sdb):
        logger.info('####### Samples treated : ' + str(i_sample + 1) + ' / ' + str(len(sdb)) + ' #######')

        EPtmp = local_error_propagation(mechanism, mysample, targets, reduction_type, alpha_norm[:, i_sample])

        EP = np.fmax(EP, EPtmp)
    logger.info('\n')
    logger.terminator('\n')

    return EP


def local_error_propagation(mechanism, mysample, targets, reduction_type, alpha_norm_loc):
    """
    Computes error propagation coefficients on one sample

    :param mechanism: Mechanism object
    :param mysample: sample on which error propagation is applied
    :param targets: list of reduction targets
    :param reduction_type: type of reduction performed (species or reactions)
    :param alpha_norm_loc: scaling coefficient for the given sample

    :return: list of error coefficients

    Created: 17/11/14 [PP]
    Last modified: 19/01/28 [QC]
    """

    myns = mechanism.ns
    mynr = mechanism.nr
    net = mechanism.network
    ctmech = mechanism.ctmech
    species_names = ctmech.species_names

    # Parameters
    EPmin = 1e-7
    EPcomp = np.zeros(myns, 'd')
    EP_spec = np.zeros(myns, 'd')
    EP_reac = np.zeros(mynr, 'd')
    EP_spec_dict = {}

    # Indices of targets
    target_indices = [species_names.index(target) for target in targets 
                            if target not in ['HeatRelease', 'HeatRelease']]
    if 'HeatRelease' in targets or 'HeatRelease' in targets:
        target_indices.append(myns)

    DIC_spec, DIC_reac = compute_DIC(mechanism, mysample, reduction_type)

    # ------------
    # Error Propagation (EP)
    # ------------

    # Go through each target
    for index_target_local, index_target_global in enumerate(target_indices):

        # Initialize working array - Works for HR as well
        EPtmp = DIC_spec[index_target_global, :] * alpha_norm_loc[index_target_local]

        # Initialize EP arrays
        array_up = np.zeros(myns, 'd')
        array_down = np.zeros(myns, 'd')

        # Initialize array_up
        array_up[:] = -1.0
        for i in range(myns):
            if i == index_target_global or EPtmp[i] < EPmin:
                continue

            array_up[i] = EPtmp[i]

        # Iterate until all relevant coefficients have been included

        flag = True
        while flag:

            # Init inner loop
            flag = False
            array_down[:] = -1.0

            # Loop over array_up
            for i in range(myns):
                indj = net.indsj[net.indsi == i]
                # If coeff is positive, store coeff and continue

                if array_up[i] > 0.0:
                    coeff_up = array_up[i]

                    # Loop over all species
                    for j in indj:
                        coeff_down = DIC_spec[i, j] * coeff_up

                        # New coeff counts if i!=j and > EPmin
                        if i != j and coeff_down > EPmin:
                            flag = True
                            # Update EPtmp and array_down for next iteration
                            if coeff_down > EPtmp[j]:
                                EPtmp[j] = coeff_down
                                array_down[j] = coeff_down
                

            if list(EPcomp) == list(EPtmp):
                flag = False
            else:
                EPcomp = EPtmp

            if targets[index_target_local] != 'HeatRelease':
                EP_spec_dict[species_names[index_target_global]] = EPtmp
            else:
                EP_spec_dict['HeatRelease'] = EPtmp

            array_up[:] = array_down[:]
        # Adjust maximum coefficient for target/species pair
        EP_spec = np.fmax(EP_spec, EPtmp)

    if reduction_type in ['reactions', 'R']:
        EPtmp = np.zeros(mynr, 'd')

        for index_target_local, index_target_global in enumerate(target_indices):
            # Initialize working array - Works for HR as well
            R_target = np.zeros(mynr, 'd')
            for index_spec in range(myns):
                if targets[index_target_local] == 'HeatRelease':
                    Rtmp = EP_spec_dict['HeatRelease'][index_spec] * DIC_reac[index_spec, :]
                else:
                    Rtmp = EP_spec_dict[species_names[index_target_global]][index_spec] \
                           * DIC_reac[index_spec, :]

                R_target = np.fmax(R_target, Rtmp)

            EPtmp[:] = R_target[:]

            # Adjust maximum coefficient for target/species pair
            EP_reac = np.fmax(EP_reac, EPtmp)

    if reduction_type in ['species', 'S']:
        EP = EP_spec
    else:
        EP = EP_reac

    return EP


def compute_DIC(mechanism, mysample, reduction_type):
    """
    Compute Direct interaction coefficients (DIC) for a given sample

    :param mechanism: Mechanism object
    :param mysample: sample for which DIC is computed
    :param reduction_type: type of reduction performed (species or reactions)

    :return: number_of_species * number_of_species matrix of coefficients for type species
             number_of_species * number_of_reactions matrix of coefficients for type reactions

    Created: 17/11/14 [PP]
    Last modified: 18/04/03 [QC]
    """

    # Parameters
    myns = mechanism.ns
    mynr = mechanism.nr
    net = mechanism.network
    ctmech = mechanism.ctmech
    nup = net.nup
    nur = net.nur
    nu = nup - nur

    # Get reaction rates for each reaction
    ctmech.TPY = float(mysample.T), float(mysample.P), [float(myY) for myY in mysample.Y]

    omega = ctmech.net_rates_of_progress

    # Get production and consumption rates for species
    PA = ctmech.creation_rates
    CA = ctmech.destruction_rates

    # Get enthalpy from production and consumption of species in reactions
    HR_species = species_heat_release(mechanism, mysample)
    HR_reactions = reactions_heat_release(mechanism, mysample)
    # HR_prod_spec = sum([HR_spec for HR_spec in HR_species if HR_spec > 0])
    # HR_cons_spec = sum([HR_spec for HR_spec in HR_species if HR_spec < 0])
    HR_prod_reac = sum([HR_reac for HR_reac in HR_reactions if HR_reac > 0])
    HR_cons_reac = sum([HR_reac for HR_reac in HR_reactions if HR_reac < 0])
    i_HR = myns

    # Workspace
    DIC_spec = np.zeros([myns + 1, myns], 'd')

    # Evaluate DIC(i,j)
    for i in range(myns):
        # reactions containing species i
        booli = net.indr[net.inds == i]
        indj = net.indsj[net.indsi == i]

        for j in indj:
            boolj = net.indr[net.inds == j]  # reactions containing species j
            indk = np.intersect1d(booli, boolj)  # reactions containing species i and j

            # Compute the DIC
            for k in indk:
                DIC_spec[i, j] += (nup[i, k] - nur[i, k]) * omega[k]

            # Normalize
            DIC_spec[i, j] = abs(DIC_spec[i, j]) / max(max(PA[i], CA[i]), 1e-60)

        # Heat release term
        DIC_spec[i_HR, i] += abs(HR_species[i])

        # Normalize
        DIC_spec[i_HR, i] = abs(DIC_spec[i_HR, i]) / max(max(abs(HR_prod_reac), abs(HR_cons_reac)), 1e-60)

    if reduction_type in ['reactions', 'R']:
        # Workspace
        DIC_reac = np.zeros([myns + 1, mynr], 'd')

        # Evaluate DIC(i,j)
        for i in range(myns):
            # reactions containing species i
            indk = net.indr[net.inds == i]  # reactions containing species i

            # Compute the DIC
            for k in indk:
                DIC_reac[i, k] += (nup[i, k] - nur[i, k]) * omega[k]

                # Normalize
                DIC_reac[i, k] = abs(DIC_reac[i, k]) / max(max(PA[i], CA[i]), 1e-60)

        for i in range(mynr):

            # Heat release term
            DIC_reac[i_HR, i] += abs(HR_reactions[i])

            # Normalize
            DIC_reac[i_HR, i] = abs(DIC_reac[i_HR, i]) / max(max(abs(HR_prod_reac), abs(HR_cons_reac)), 1e-60)

    else:
        DIC_reac = None

    return DIC_spec, DIC_reac


def integrity_check(mechanism):
    """
    Dictionaries with intel about integrity

    Gives dictionaries telling:
    - which species are no more consumed if a given one (the key of the dict) is removed
    - which species are no more produced if a given one (the key of the dict) is removed
    - which species it would be wise to remove if a given one (the key of the dict) is removed

    :param mechanism: Mechanism object

    :return the 3 dictionaries previously described in the same order

    Created: 18/04/03 [QC]
    """
    myns = mechanism.ns
    net = mechanism.network
    ctmech = mechanism.ctmech
    species_names = ctmech.species_names
    nup = net.nup
    nur = net.nur
    nu = nup - nur

    # Integrity dict: grouping species that induce truncated path
    no_more_produced = {}
    no_more_consumed = {}
    # Global dictionary giving the species that needs to be removed if key is removed
    remove_if_removed = {}

    for i in range(myns):
        nu_integrity = nu.copy()
        # reactions containing species i
        indj = net.indsj[net.indsi == i]

        for j in indj:

            # Integrity check
            if [x for x in nu_integrity[j, :] if x != 0]:
                if np.max([x for x in nu_integrity[j, :] if
                           x != 0]) < 0:  # if species i is removed species j is no more produced
                    if species_names[i] not in no_more_produced:
                        no_more_produced[species_names[i]] = []
                    else:
                        if species_names[j] not in no_more_produced[species_names[i]]:
                            no_more_produced[species_names[i]].append(species_names[j])
                    # Global dict
                    if species_names[i] not in remove_if_removed:
                        remove_if_removed[species_names[i]] = []
                    else:
                        if species_names[j] not in remove_if_removed[species_names[i]]:
                            remove_if_removed[species_names[i]].append(species_names[j])

                elif np.min([x for x in nu_integrity[j, :] if
                             x != 0]) > 0:  # if species i is removed species j is no more consumed
                    if species_names[i] not in no_more_consumed:
                        no_more_consumed[species_names[i]] = []
                    else:
                        if species_names[j] not in no_more_consumed[species_names[i]]:
                            no_more_consumed[species_names[i]].append(species_names[j])
                    # Global dict
                    if species_names[i] not in remove_if_removed:
                        remove_if_removed[species_names[i]] = []
                    else:
                        if species_names[j] not in remove_if_removed[species_names[i]]:
                            remove_if_removed[species_names[i]].append(species_names[j])

    return no_more_consumed, no_more_produced, remove_if_removed


def scaling_coefficient(mechanism, sdb, targets):
    """
    Computes scaling coefficient for DRGEP

    :param mechanism: Mechanism object
    :param sdb: sample database on which DRGEP is applied
    :param targets: list of reduction targets

    :return: number_of_targets * number_of_samples matrix of coefficients

    Created: 18/04/03 [QC]
    """

    ctmech = mechanism.ctmech
    species_names = ctmech.species_names

    elements = ctmech.element_names

    # Indices of targets
    target_indices = [species_names.index(target) for target in targets if target not in ['HeatRelease', 'HeatRelease']]
    target_indices.append(len(species_names))

    # Heat release
    try:
        i_HR = targets.index('HeatRelease')
    except ValueError:
        i_HR = None

    # Store composition of species in a n_elements*n_species matrix
    element_in_species = np.zeros([len(elements), len(species_names)])
    for index_spec, spec in enumerate(species_names):
        comp = ctmech.species(index_spec).composition
        for index_element, element in enumerate(elements):
            if element in comp:
                element_in_species[index_element, index_spec] = comp[element]

    # Identifying the different cases in the samples
    cases = list(set([sample.case.myid for sample in sdb]))
    n_cases = len(cases)

    # Initialize workspace
    n_samples = len(sdb)
    n_elements = len(elements)
    n_targets = len(targets)

    ind_e, ind_s = np.nonzero(element_in_species)
    alpha = np.zeros([n_elements, n_targets, n_samples])
    alpha_temp = np.zeros([n_elements, n_targets, n_samples])
    max_alpha = np.ones([n_elements, n_targets, n_cases])  # Ones instead of zeros to avoid zero division
    alpha_norm = np.zeros([n_targets, n_samples])

    for i_sample, mysample in enumerate(sdb):

        # Sample case index
        case_index = cases.index(mysample.case.myid)

        # Get reaction rates for each reaction
        ctmech.TPY = float(mysample.T), float(mysample.P), [float(myY) for myY in mysample.Y]

        # Get production and consumption rates for species
        PA = ctmech.creation_rates
        CA = ctmech.destruction_rates

        # Heat Release
        HR = float(mysample.HR)

        # ------------
        # Scaling coefficient
        # ------------

        # Element production
        for index_element, element in enumerate(elements):
            P_element = 0.0
            # Species containing the element
            spec_with_element = ind_s[ind_e == index_element]

            for index_spec in spec_with_element:
                # P_e = Sum(N_e,s * max(P_s - C_s, 0))
                P_element += max(element_in_species[index_element, index_spec] * (PA[index_spec] - CA[index_spec]),
                                 1e-60)

            for index_target, target in enumerate(targets):
                if target != 'HeatRelease':
                    # alpha_e,t = (N_e,s * |P_t - C_t|) / P_e
                    alpha[index_element, index_target, i_sample] = element_in_species[index_element, int(
                        target_indices[index_target])] * abs(
                            PA[target_indices[index_target]] - CA[target_indices[index_target]]) / P_element

                    if alpha[index_element, index_target, i_sample] > max_alpha[index_element, index_target,
                                                                                case_index]:

                        max_alpha[index_element, index_target, case_index] = alpha[index_element, index_target, i_sample]

        if i_HR is not None:
            alpha[0, i_HR, i_sample] = abs(HR)
            if alpha[0, i_HR, i_sample] > max_alpha[0, i_HR, case_index]:
                max_alpha[0, i_HR, case_index] = alpha[0, i_HR, i_sample]

    # Reconstructing
    max_alpha_temp = np.zeros([n_elements, n_targets, n_samples])
    for index, sample in enumerate(sdb):
        # Sample case index
        case_index = cases.index(sample.case.myid)
        max_alpha_temp[:, :, index] = max_alpha[:, :, case_index]

    max_alpha = max_alpha_temp

    for index_element, element in enumerate(elements):
        for index_target, target in enumerate(targets):
            alpha_temp[index_element, index_target, :] = alpha[index_element, index_target, :] \
                                                             / max_alpha[index_element, index_target, :]

    # Normalization: alpha_t = max_e(alpha_e,t/max_sample(alpha_e,t))
    for index_target, target in enumerate(targets):
        alpha_norm[index_target, :] = np.max(alpha_temp[:, index_target, :], 0)

    return alpha_norm


def important_species_identification(caselist, mechanism):
    """
    Identifies which the important species (fuel, oxidizer, products and diluents)
    Those species are the one with a mass faction representing more than 1% of the fresh and burnt gases

    :param caselist: list of Case objects
    :param mechanism: Mechanism object

    :return: list of important species
    """

    species_names = mechanism.ctmech.species_names

    important_species = []

    for case in caselist:
        data = case.extract_profile(mechanism)
        data_names_dict = case.names_dictionary(mechanism)
        for spec in species_names:
            if spec not in important_species:
                y_start = data[0, data_names_dict[spec]]
                y_end = data[-1, data_names_dict[spec]]

                if y_start > 1e-3 or y_end > 1e-2:
                    important_species.append(spec)

                if spec in case.targets:
                    important_species.append(spec)

    return important_species


def species_heat_release(mechanism, sample):
    """
    Computes the heat release rate from each species

    :param mechanism: Mechanism object
    :param sample: Sample object

    :return hr_species: list of heat release rate of each species
    """

    ctmech = mechanism.ctmech
    net = mechanism.network

    species_names = ctmech.species_names

    ctmech.TPY = float(sample.T), float(sample.P), [float(myY) for myY in sample.Y]

    hr_species = np.zeros((len(species_names)))
    hr_reactions = reactions_heat_release(mechanism, sample)
    for i_spec, spec in enumerate(species_names):
        # reactions containing species i
        booli = net.indr[net.inds == i_spec]
        for i_reac in booli:
            hr_species[i_spec] = hr_species[i_spec] + hr_reactions[i_reac]

    return hr_species


def reactions_heat_release(mechanism, sample):
    """
    Computes the heat release rate from each reaction

    :param mechanism: Mechanism object
    :param sample: Sample object

    :return hr_reactions: list of heat release rate of each reaction
    """

    ctmech = mechanism.ctmech

    n_reac = len(ctmech.reactions())

    ctmech.TPY = float(sample.T), float(sample.P), [float(myY) for myY in sample.Y]

    hr_reactions = [- ctmech.net_rates_of_progress[i_reac] * ctmech.delta_enthalpy[i_reac] for i_reac in range(n_reac)]

    return hr_reactions



############### BELOW IS KELLY!!! #######################################################

def graph_error_propagation(mechanism, sdb, targets):
    """
    Computes error propagation coefficients and outputs parameters used to 
    create mechanism graph with "paths" between target and each species.
    Modified version of error_propagation.

    :param mechanism: Mechanism object
    :param sdb: sample database on which DRGEP is applied
    :param targets: list of reduction targets

    :return EP: dictionary of species or reactions with their corresponding coefficient
    :return ind_list_f: list of species that might have contributed to given species' EP
    :return DIC: 3d matrix composed of DIC_spec matrix at each sample time
    :return DIC_ind: list of the sample index that each EP was calculated at
    :return coeffs: list of intermediate coefficient matrix at each sample time
    :return alpha_norm: list of scaling coefficient at each sample time


    Created: 22/3/3 [KI]
    """

    myns = mechanism.ns

    EP = np.zeros(myns, 'd')

    alpha_norm = scaling_coefficient(mechanism, sdb, targets)

    logger.terminator('\r')
    
    # Count number of samples
    N = 0
    for i_sample, mysample in enumerate(sdb):
        N = N+1
    
    # Initialize working arrays/lists
    EP_all = np.zeros((myns,N))
    DIC_ind = np.zeros((myns))
    DIC = np.zeros((myns,myns, N))
    ind_list = []
    ind_list_f = []
    coeffs = []
    
    for i_sample, mysample in enumerate(sdb):
        logger.info('####### Samples treated : ' + str(i_sample + 1) + ' / ' + str(len(sdb)) + ' #######')

        # Calculate EPs for this sample
        [EP_i, ind_list_i,DIC_i, coeffs_i] = graph_local_error_propagation(mechanism, mysample, targets, alpha_norm[:, i_sample])
    
        # Fill matrix/lists with info from this samplee
        EP_all[:,i_sample] = EP_i
        ind_list.append(ind_list_i)
        DIC[:,:,i_sample] = DIC_i[:-1,:]
        coeffs.append(coeffs_i[:,1:])
        
    # Take the max EP for each species over all samples    
    EP = np.max(EP_all, axis = 1)
        
    # Make final ind_list 
    for i in range(myns):
        ind_list_f.append([])
    
    # Fill out remaininfg outputs for each species
    for i in range(myns):
        if sum(EP[i] == EP_all[i,:]) == 1:
            # Find the sample that this species' EP came from
            ind = np.where(EP[i] == EP_all[i,:])
            
            # Store the ind_list from this sample
            ind_list_f[i] = ind_list[int(ind[0])][i]
            
            # Store the index of this sample
            DIC_ind[i] = int(ind[0])          
            
        elif sum(EP[i] == EP_all[i,:]) > 1:
            # If the EP could have been calculated from multiple samples (value
            #  is likely 0 or 1) take the last sample it could have been
            #  calculated from
            print('Multiple Matches: ', mechanism.species_names[i])
            ind_list_f[i] = ind_list[int(np.where(EP[i] == EP_all[i,:])[0][-1])][i]
            DIC_ind[i] = np.where(EP[i] == EP_all[i,:])[0][-1]
        else:
            print('Error: ',  mechanism.species_names[i])
        
        
        
        
    logger.info('\n')
    logger.terminator('\n')

    return [EP, ind_list_f, DIC, DIC_ind, coeffs, alpha_norm]


def graph_local_error_propagation(mechanism, mysample, targets, alpha_norm_loc):
    """
    Computes error propagation coefficients on one sample and outputs parameters
    used to create mechanism graph with "paths" between target and each species.
    Modified version of local_error_propagation

    :param mechanism: Mechanism object
    :param mysample: sample on which error propagation is applied
    :param targets: list of reduction targets
    :param alpha_norm_loc: scaling coefficient for the given sample

    :return EP: list of error coefficients
    :return ind_list: list of species that might have contributed to given 
        species' EP (for this sample)
    :return DIC_spec: DIC used to calculate EP
    :return coeffs: intermediate coefficients in EP calculation

    Created: 22/3/3 [KI]

    """

    myns = mechanism.ns
    net = mechanism.network
    ctmech = mechanism.ctmech
    species_names = ctmech.species_names
    reduction_type = 'species'

    # Parameters
    EPmin = 1e-7
    EPcomp = np.zeros(myns, 'd')
    EP_spec = np.zeros(myns, 'd')
    EP_spec_dict = {}

    # Indices of targets
    target_indices = [species_names.index(target) for target in targets 
                            if target not in ['HeatRelease', 'HeatRelease']]
    if 'HeatRelease' in targets or 'HeatRelease' in targets:
        target_indices.append(myns)

    DIC_spec, DIC_reac = compute_DIC(mechanism, mysample, reduction_type)
    coeffs = np.zeros((myns, 1))
    
    # Initialize list to track species that change a given species' EP at some point
    #  --> the last species in this list is adjascent to the given species in its path
    ind_list = []
    for i in range(myns):
        ind_list.append([])

    # ------------
    # Error Propagation (EP)
    # ------------

    # Go through each target
    for index_target_local, index_target_global in enumerate(target_indices):

        # Initialize working array - Works for HR as well
        EPtmp = DIC_spec[index_target_global, :] * alpha_norm_loc[index_target_local]

        # Initialize EP arrays
        array_up = np.zeros(myns, 'd')
        array_down = np.zeros(myns, 'd')

        # Initialize array_up
        array_up[:] = -1.0
        for i in range(myns):
            if i == index_target_global or EPtmp[i] < EPmin:
                continue

            array_up[i] = EPtmp[i]


        # Iterate until all relevant coefficients have been included
        flag = True
        while flag:

            # Init inner loop
            flag = False
            array_down[:] = -1.0

            # Loop over array_up
            for i in range(myns):
                indj = net.indsj[net.indsi == i]
                # If coeff is positive, store coeff and continue

                if array_up[i] > 0.0:
                    coeff_up = array_up[i]

                    # Loop over all species
                    for j in indj:
                        coeff_down = DIC_spec[i, j] * coeff_up

                        # New coeff counts if i!=j and > EPmin
                        if i != j and coeff_down > EPmin:
                            flag = True
                            # Update EPtmp and array_down for next iteration
                            if coeff_down > EPtmp[j]:
                                EPtmp[j] = coeff_down
                                array_down[j] = coeff_down
                                ind_list[j].append(i) 
                                # record that species i affected species j's EP


            if list(EPcomp) == list(EPtmp):
                flag = False
            else:
                EPcomp = EPtmp

            if targets[index_target_local] != 'HeatRelease':
                EP_spec_dict[species_names[index_target_global]] = EPtmp
            else:
                EP_spec_dict['HeatRelease'] = EPtmp

            array_up[:] = array_down[:]
            
            # Store intermediate coefficients - this helps determine how many 
            # "levels" are in the path of a given species
            old_coeffs = array_down[:].reshape(array_down.shape[0],1)
            coeffs = np.append(coeffs, old_coeffs,1)
            
        # Adjust maximum coefficient for target/species pair
        EP_spec = np.fmax(EP_spec, EPtmp)

    EP = EP_spec

    return [EP, ind_list, DIC_spec, coeffs]
                                
                                
                