import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis import *
from MDAnalysis.analysis.align import *
from MDAnalysis.lib.distances import distance_array
from tqdm import tqdm
import seaborn as sns
import pandas as pd
from scipy.signal import savgol_filter
from itertools import groupby
import networkx as nx
import re

def longest_contiguous_stretch_of_ones(input_string: str):
    """
    Function inpired by and modified from: https://stackoverflow.com/a/75069897
    Input: A string of the form 000011111100, where 1's represent residues with beta-contacts and 0's represent those without.

    Returns the longest contiguous stretch of 1's in the given input.
    """

    # groups the input into substrings: contiguous stretches of 0's only or 1's only
    substring_generator = ("".join(g) for _, g in groupby(input_string))

    # Sort by length
    sorted_substrings = sorted(
        # filter for substrings whose first (or in fact, any) element is 1
        [sub for sub in list(substring_generator) if sub[0] == '1'],
        key=lambda x: -len(x)
    )

    result = [list(g) for k, g in groupby(
        sorted_substrings, key=lambda x: len(x))][0][0]

    result_span = re.search(pattern=result, string=input_string).span()
    return result, result_span


def contiguous_beta_sheet(a, min_contiguous_contacts=4):

    if np.sum(a) < 1:
        return (False, (0, 0), (0, 0))

    last_relative_diagonal = a.shape[0]-min_contiguous_contacts
    for k in range(-last_relative_diagonal, last_relative_diagonal+1):
        # print(a)
        # break
        dia = np.diag(a, k=k)
        # print(dia)
        if np.sum(dia) >= min_contiguous_contacts:
            # Convert the diagonal into a string (s1),
            s1 = "".join(
                map(
                    str,
                    dia,
                )
            )

            contig_stretch, contig_span = longest_contiguous_stretch_of_ones(
                s1)
            if len(contig_stretch) >= min_contiguous_contacts:
                # print(k, contig_span)
                if k >= 0:
                    pep_i_span = (contig_span[0], contig_span[1])
                    pep_j_span = (k+contig_span[0], k+contig_span[1])
                elif k < 0:
                    pep_i_span = (abs(k)+contig_span[0], abs(k)+contig_span[1])
                    pep_j_span = (contig_span[0], contig_span[1])
                # print(pep_i_span, pep_j_span)
                # if this condition is met, there is no need to check any other diagonals
                return (True, pep_i_span, pep_j_span)

        anti_dia = np.diag(np.flipud(a), k=k)
        if np.sum(anti_dia) >= min_contiguous_contacts:
            # Convert the anti-diagonal into a string (s2),
            s2 = "".join(
                map(
                    str,
                    anti_dia,
                )
            )
            # print(s2)
            contig_stretch, contig_span = longest_contiguous_stretch_of_ones(
                s2)

            if len(contig_stretch) >= min_contiguous_contacts:
                # print(k, contig_span)
                if k >= 0:
                    pep_i_span = (contig_span[0], contig_span[1])
                    pep_j_span = (k+contig_span[0], k+contig_span[1])
                elif k < 0:
                    pep_i_span = (abs(k)+contig_span[0], abs(k)+contig_span[1])
                    pep_j_span = (contig_span[0], contig_span[1])
                # if this condition is met, there is no need to check any other diagonals
                return (True, pep_i_span, pep_j_span)

    return (False, (0, 0), (0, 0))


def e2e_distance_calculation(u, ts, selection_string="name BB"):

    all_peps = u.select_atoms(selection_string)
    BB_BB_distances = []
    for pep in all_peps.fragments:
        pep = pep.select_atoms('name BB')
        dists = distance_array(pep.positions, pep.positions, box=ts.dimensions)
        BB_BB_distances += [dists[0, -1]]

    return BB_BB_distances


def agg_num(u, ts, selection_string="name BB S1 S2 S3", cutoff=7):
    """
    Returns a graph object of aggregates, identified by contacts between BB groups
    Nodes of returned graph are numbered by peptide id
    """

    n_beads_per_peptide = None
    prot_bb = u.select_atoms(selection_string)
    pep_num = len(prot_bb.fragments)
    for temp in prot_bb.fragments[:1]:
        n_beads_per_peptide = len(temp.select_atoms(selection_string))
        continue

    dist_array = distance_array(
        prot_bb.positions, prot_bb.positions, np.copy(ts.dimensions))
    ct_ids = np.where(dist_array < cutoff)  # id's of beads/groups in contact

    # converts BB-bead index to peptide index
    E = ct_ids[0]//n_beads_per_peptide
    # converts BB-bead index to peptide index
    E1 = ct_ids[1]//n_beads_per_peptide
    G = nx.Graph()
    for i in range(len(E)):
        G.add_edge(int(E[i]), int(E1[i]))  # created edge between peptides

    G.remove_edges_from(nx.selfloop_edges(G))

    agg_num_dict = {pep_id: 1 for pep_id in range(pep_num)}
    connected_comps = [list(c) for c in sorted(
        nx.connected_components(G), key=len)]

    for comp in connected_comps:
        agg_size = len(comp)
        for pep_id in comp:
            agg_num_dict[pep_id] = agg_size

    agg_num_list = [agg_num for (pep_id, agg_num)
                                 in sorted(agg_num_dict.items())]

    return agg_num_list


def is_helical_dipole_contacts(u, ts, helical_contact_cutoff=3):

    result = []
    pep_num = len(u.select_atoms("name BBp").fragments)
    n_res = len(u.select_atoms(
        "name BBp").fragments[0].select_atoms("name BBp").resids)
    bbp = u.select_atoms("name BBp").positions.reshape((pep_num, n_res, 3))
    bbm = u.select_atoms("name BBm").positions.reshape((pep_num, n_res, 3))

    for aa_i in range(n_res-4):
        first_bbm_positions = bbm[:, aa_i, :]
        first_bbp_positions = bbp[:, aa_i, :]
        last_bbm_positions = bbm[:, aa_i+4, :]
        last_bbp_positions = bbp[:, aa_i+4, :]

        # To find euclidean distance between a1 and a2 (both are lists of position vectors (n x 3))
        # distance = np.linalg.norm(a1 - a2, axis=1)

        alpha_contact = np.any([np.linalg.norm(first_bbm_positions-last_bbp_positions, axis=1) < helical_contact_cutoff,
                                np.linalg.norm(first_bbp_positions-last_bbm_positions, axis=1) < helical_contact_cutoff],
                                axis=0)  # operate over the length of the resulting Boolean arrays

        result.append(alpha_contact.astype(int))

    return result


def rgyr(u, ts, selection_string="name BB BBm BBp S1 S1m S1p S2 S2m S2p S3 S4", cutoff=7, pep_num=50):

    all_peps = u.select_atoms(selection_string)
    rgyrs = []
    for pep in all_peps.fragments:
        rgyrs += [pep.radius_of_gyration()]

    # agg_num_list  = [ agg_num for (pep_id, agg_num) in sorted(agg_num_dict.items())]

    return rgyrs


def beta_components_time_series_dictionary(
    u,
    start = 0,
    stop = None,
    step = 100,
    bb_bb_cutoff = 7.,
    min_contig = 4,
    n_res = 7
) -> dict:

    index_thing=[]
    bb=u.select_atoms('name BB')
    beads_per_pep=None
    for i, pep in enumerate(bb.fragments):
        pep_group=pep.select_atoms('name BB')
        index_thing += [f'index {pep_group[0].index} to {pep_group[-1].index+4}']
        beads_per_pep=len(pep_group)
        # print(id, f"(index {pep.select_atoms('name BB').indices[0]} to {pep.select_atoms('name BB').indices[-1]+4})",)
        # break
    pep_num=len(bb.fragments)

    times=[]
    pep_ids=[]
    beta_spans=[]
    beta_spans_size=[]
    cluster_strings=[]
    cluster_lengths=[]
    vmd_group=[]
    e2e_ds=[]
    pep_rgyrs=[]
    agg_nums=[]
    # helical_contacts=[[]]*(n_res-4)
    # dist_to_hep=[]
    # polar_order=[]
    # nematic_order=[]
    with tqdm(total=len(u.trajectory[start:stop:step])) as pbar:
        for ts in u.trajectory[start:stop:step]:

            ts_times=[None]*pep_num
            ts_pep_ids=[None]*pep_num
            ts_beta_spans=[None]*pep_num
            ts_beta_spans_size=[None]*pep_num
            ts_cluster_strings=['']*pep_num
            ts_cluster_lengths=[0]*pep_num
            ts_e2e_ds=e2e_distance_calculation(u, ts)
            ts_pep_rgyrs=rgyr(u, ts)
            ts_agg_nums=agg_num(u, ts)

            # pol, nem=nematic_and_polar_order(u, ts, pep_num = pep_num)
            # ts_polar_order=[pol]*pep_num
            # ts_nematic_order=[nem]*pep_num

            # tmp_spans = [[]]*pep_num
            tmp_spans=[[] for pep in range(pep_num)]
            # Distance array of all BB's to all other BB's
            bb_bb=distance_array(bb.positions, bb.positions, ts.dimensions)

            # Two modifications to this array: removal of intrapeptide contacts,
            # and removal of all lower triangle elements (array is symmetric, any one triangle is sufficient)
            for pep_id in range(pep_num):
                # to remove these elements from the analysis, just add any number > bb_bb-cutoff
                bb_bb[
                    pep_id*beads_per_pep: pep_id*beads_per_pep+beads_per_pep,
                    : pep_id*beads_per_pep+beads_per_pep] += bb_bb_cutoff+10

            # Convert to a boolean array -> a_ij = 1 if contact exists, else 0.
            bb_bb_bool= (bb_bb <= bb_bb_cutoff).astype(int)

            # beta-sheet array
            beta_sheet_array= np.zeros((pep_num, pep_num), dtype=int)

            pep_i= 0
            while pep_i < pep_num:
                # To parse upper triangle, skip all pep_ids <= pep_i
                pep_j= pep_i+1
                while pep_j < pep_num:
                    ij_contacts = bb_bb_bool[pep_i*beads_per_pep: (pep_i+1)*beads_per_pep,
                                            pep_j*beads_per_pep: (pep_j+1)*beads_per_pep
                                            ]
                    # print(pep_i, pep_j)
                    contig_sheet_exists, ispan, jspan= contiguous_beta_sheet(ij_contacts, min_contiguous_contacts=min_contig)
                    # print(ispan, jspan)
                    beta_sheet_array[pep_i, pep_j]= int(contig_sheet_exists)
                    tmp_spans[pep_i] += [ispan]
                    tmp_spans[pep_j] += [jspan]
                    pep_j += 1
                pep_i += 1

            for pep_n in range(pep_num):
                ts_times[pep_n]= ts.time//1000
                # sorry... this is complex
                # sort the spans from smallest to largest for all the possible beta-sheet regions
                # add the longest span to the dataframe. This may need to be changed if the smaller fragments are distinct and significant
                # print(ts_beta_spans)
                # print(ts_beta_spans[pep_n])
                longest_span = sorted([(len(range(*span)), span) for span in tmp_spans[pep_n] if span != None])[: : -1][0]

                ts_beta_spans[pep_n]= f'({longest_span[1][0]}_{longest_span[1][1]})'
                ts_beta_spans_size[pep_n]= len(range(*longest_span[1]))

                ts_pep_ids[pep_n]= pep_n


            # Make a graph object of the beta sheet array
            G= nx.Graph(beta_sheet_array)
            # Detect connected components, i.e., beta sheets connected via BB contacts
            beta_components= [list(c) for c in sorted(
                nx.connected_components(G), key=len, reverse=True)]

            for comp in beta_components:
                # print(comp, type(comp))
                if len(comp)>1: comp.sort()
                cluster_string= '(' + '_'.join(map(str, comp)) + ')'
                cluster_length= len(comp)
                for pep_n in comp:
                    ts_cluster_lengths[pep_n]= cluster_length
                    ts_cluster_strings[pep_n]= cluster_string

            # helical_result = is_helical_dipole_contacts(u, ts)
            # for i, t in enumerate(helical_result):
            #     helical_contacts[i]= np.append(helical_contacts[i], t)

            times += ts_times
            pep_ids += ts_pep_ids
            beta_spans += ts_beta_spans
            beta_spans_size += ts_beta_spans_size
            cluster_strings += ts_cluster_strings
            cluster_lengths += ts_cluster_lengths
            vmd_group += index_thing
            e2e_ds += ts_e2e_ds
            pep_rgyrs += ts_pep_rgyrs
            agg_nums += ts_agg_nums
            # dist_to_hep += ts_dist_to_hep
            # polar_order += ts_polar_order
            # nematic_order += ts_nematic_order

            pbar.update(1)


    data = {
        'time': pd.Series(times), 
        'pep_id': pd.Series(pep_ids, dtype=str), 
        'e2e_ds': pd.Series(e2e_ds, dtype=float).round(decimals=3), 
        'pep_rgyrs': pd.Series(pep_rgyrs, dtype=float).round(decimals=3), 
        # 'dist_to_hep': pd.Series(dist_to_hep, dtype=float).round(decimals=3), 
        'agg_nums': pd.Series(agg_nums, dtype=str), 
        'beta_span': pd.Series(beta_spans, dtype=str), 
        'beta_size': pd.Series(beta_spans_size, dtype=int), 
        'cluster_lengths': pd.Series(cluster_lengths, dtype=int), 
        'cluster_strings': pd.Series(cluster_strings, dtype=str), 
        'vmd_group': pd.Series(vmd_group, dtype=str), 
        # 'P1': pd.Series(polar_order), 
        # 'P2': pd.Series(nematic_order)
    }
    # data.update(
    #     {
    #         f'helical{i+1}': pd.Series(helical) for i, helical in enumerate(helical_contacts)
    #     }
    # )
    return data

def in_contact_with_heparin(u, ts,
                          cutoff=7.,
                          hep_selection="name BG",
                          prot_selection="name BB S1 S2 S3",
                          ):

    hep = u.select_atoms(hep_selection)
    prot = u.select_atoms(prot_selection)

    min_dist = []

    if len(hep) > 0:
        for pep in prot.fragments:
            pep = pep.select_atoms(prot_selection)

            pep_hep = distance_array(
                pep.center_of_mass(), hep.positions, np.copy(ts.dimensions))

            min_dist += [np.min(pep_hep)]

    in_contact = np.array(min_dist) <= cutoff

    return in_contact

def hep_contact_per_peptide(
    u,
    ts,
    pep_particles = "name BB S1 S2 S3",
    hep_particles = "name BG B2 B3 B6",
    cutoff=7,
) -> list:

    """
    Returns a list, of number of contacts between each peptide and a heparin strand.
    """
    pep = u.select_atoms(pep_particles)
    hep = u.select_atoms(hep_particles)

    n_pep = u.select_atoms("resname LYS and name BB").n_atoms
    
    pep_hep_distances = distance_array(
        pep.positions,
        hep.positions,
        box=ts.dimensions
    )

    # Convert the array to boolen format, 
    # where an A_ij = 1 if the the distance between i and j are <= cutoff
    pep_hep_distances = (pep_hep_distances <= cutoff).astype(int)
    
    # Split the array into sub-arrays - one per peptide
    # Get the sum of each sub-array to get the number of heparin-peptide contacts
    n_contacts =  [np.sum(x) for x in np.array_split(pep_hep_distances, n_pep)]
    
    return n_contacts