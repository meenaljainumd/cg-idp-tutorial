#!/usr/local/bin/python2.7
import numpy as np
import os
import sys
# sys.argv.append('-test')
import math
import time
from MDAnalysis import *
import pandas as pd
# from Bio.PDB import PDBParser
# from Bio.PDB.DSSP import DSSP
# p = PDBParser()

import argparse

# Create the parser and add arguments
parser = argparse.ArgumentParser()
parser.add_argument('-o', '--output_dir', type=str,
					default='/home/sgotla/heparin/polarizable_heparin_cg/make_ff_itp/protein_itp_outputs/', 
                    help='Output directory')
parser.add_argument('-f', '--input_seq', type=str,
					default='seq.txt', help='Input Sequence')
# parser.add_argument('-c', '--charged', type=int, default=1,
# 					help='Where to start numbering charged beads')
# parser.add_argument('-a', '--aromatic', type=int, default=1,
# 					help='Where to start numbering aromatic beads')
# parser.add_argument('-n', '--protno', type=int,
# 					default=50, help='Peptide number')
parser.add_argument('-count', '--peptide_count', type=int,
					default=50, help='Number of peptides in the system')
parser.add_argument('-bbdih', '--backbone_dih_type', type=int,
					default=6, help='Type of backbone dihedral. Default: 6 (double well potential)')
parser.add_argument('-bbdih_fc', '--backbone_dih_force_constant', type=int,
					default=2, help='Force constant for backbone dihedrals. Default: 2 (kJ/mol)')
parser.add_argument('-cys', '--make_cys_bridge', type=int,
					default=0, help='Assuming there are 2 cysteines in the peptide, turns cysteine bridges on or off.')
# Parse and print the results
args = parser.parse_args()
# print(args.argument1)

input_seq = args.input_seq
with open(input_seq, 'r') as f:
    seq = f.readlines()[0].strip()
    print(seq)
output_itp = args.output_dir
# cb = args.charged
# arb = args.aromatic
n_peptides = args.peptide_count
out_dir = args.output_dir
bb_dih_type = args.backbone_dih_type
bb_dih_fc = args.backbone_dih_force_constant
make_cys_bridge = bool(args.make_cys_bridge)
print(make_cys_bridge)
# if make_cys_bridge == 0:
#     make_cys_bridge = False
# elif make_cys_bridge == 'yes':
#     make_cys_bridge = True
# else:
#     exit()
# cys_resids = [resid+1 for resid, aa in enumerate(seq) if aa=='C']

# "/home/sgotla/heparin/polarizable_heparin_cg/make_ff_itp/protein_itp_outputs"


######################
# Main
######################
# If PDB is the input -- use that; else use seq.txt if exception arise.


seq_to_bead = {
	'A': ['BB', 'BBp', 'BBm'],
	'G': ['BB', 'BBp', 'BBm'],
	'S': ['BB', 'BBp', 'BBm', 'S1', 'S1p', 'S1m'],
	'I': ['BB', 'BBp', 'BBm', 'S1'],
	'L': ['BB', 'BBp', 'BBm', 'S1'],
	'V': ['BB', 'BBp', 'BBm', 'S1'],
	'C': ['BB', 'BBp', 'BBm', 'S1'],
	'M': ['BB', 'BBp', 'BBm', 'S1'],
	'D': ['BB', 'BBp', 'BBm', 'S1'],
	'Q': ['BB', 'BBp', 'BBm', 'S1', 'S1p', 'S1m'],
	'N': ['BB', 'BBp', 'BBm', 'S1', 'S1p', 'S1m'],
	'P': ['BB', 'BBp', 'BBm', 'S1'],
	'T': ['BB', 'BBp', 'BBm', 'S1', 'S1p', 'S1m'],
	'K': ['BB', 'BBp', 'BBm', 'S1', 'S2'],
	'R': ['BB', 'BBp', 'BBm', 'S1', 'S2'],
	'E': ['BB', 'BBp', 'BBm', 'S1', 'S2'],
	'W': ['BB', 'BBp', 'BBm', 'S1', 'S2', 'S2p', 'S2m', 'S3', 'S4'],
	'F': ['BB', 'BBp', 'BBm', 'S1', 'S2', 'S3'],
	'Y': ['BB', 'BBp', 'BBm', 'S1', 'S2', 'S3', 'S3p', 'S3m'],
	'H': ['BB', 'BBp', 'BBm', 'S1', 'S2', 'S2p', 'S2m', 'S3', 'S3p', 'S3m'],
	'Rsd': ['BB', 'BBp', 'BBm', 'S1', 'S1p', 'S1m', 'S2', 'S2p', 'S2m']
}
seq_to_type = {
	'A': ['PP5', 'D2', 'D2'],
	'G': ['PP5', 'D2', 'D2'],
	'S': ['PP5', 'D2', 'D2', 'PN0', 'D2', 'D2'],
	'I': ['PP5', 'D2', 'D2', 'C1'],
	'L': ['PP5', 'D2', 'D2', 'C1'],
	'V': ['PP5', 'D2', 'D2', 'C3'],
	'C': ['PP5', 'D2', 'D2', 'C3'],
	'M': ['PP5', 'D2', 'D2', 'C3'],
	'D': ['PP5', 'D2', 'D2', 'Qa'],
	'Q': ['PP5', 'D2', 'D2', 'PNda', 'D2', 'D2'],
	'N': ['PP5', 'D2', 'D2', 'PNda', 'D2', 'D2'],
	'P': ['APP5', 'D2', 'D2', 'PAC1'],
	'T': ['PP5', 'D2', 'D2', 'PNda', 'D2', 'D2'],
	'K': ['PP5', 'D2', 'D2', 'C3', 'Qd'],
	'R': ['PP5', 'D2', 'D2', 'N0', 'Qd'],
	'E': ['PP5', 'D2', 'D2', 'C3', 'Qa'],
	'W': ['PP5', 'D2', 'D2', 'AR', 'ARP', 'D2A', 'D2A', 'AR', 'AR'],
	'F': ['PP5', 'D2', 'D2', 'AR', 'AR', 'AR'],
	'Y': ['PP5', 'D2', 'D2', 'AR', 'AR', 'ARP', 'D2A', 'D2A'],
	'H': ['PP5', 'D2', 'D2', 'AR', 'ARP', 'D2', 'D2', 'ARP', 'D2', 'D2'],
	'Rsd': ['PP5', 'D2', 'D2', 'PP5', 'D2', 'D2', 'PP5', 'D2', 'D2']
}
seq_to_cgnr = {
	'A': [1, 0, 1],
	'G': [1, 0, 1],
	'S': [1, 0, 1, 1, 0, 1],
	'I': [1, 0, 1, 1],
	'L': [1, 0, 1, 1],
	'V': [1, 0, 1, 1],
	'C': [1, 0, 1, 1],
	'M': [1, 0, 1, 1],
	'D': [1, 0, 1, 1],
	'Q': [1, 0, 1, 1, 0, 1],
	'N': [1, 0, 1, 1, 0, 1],
	'P': [1, 0, 1, 1],
	'T': [1, 0, 1, 1, 0, 1],
	'K': [1, 0, 1, 1, 1],
	'R': [1, 0, 1, 1, 1],
	'E': [1, 0, 1, 1, 1],
	'W': [1, 0, 1, 1, 1, 1, 0, 1, 1],
	'F': [1, 0, 1, 1, 1, 1],
	'Y': [1, 0, 1, 1, 1, 1, 0, 1],
	'H': [1, 0, 1, 1, 1, 0, 1, 1, 0, 1],
	'Rsd': [1, 0, 1, 1, 0, 1, 1, 0, 1]
}
seq_to_charge = {
	'A': [0, 0.34, -0.34],
	'G': [0, 0.34, -0.34],
	'S': [0, 0.34, -0.34,  0.00, 0.144, -0.144],
	'I': [0, 0.34, -0.34,  0.00],
	'L': [0, 0.34, -0.34,  0.00],
	'V': [0, 0.34, -0.34,  0.00],
	'C': [0, 0.34, -0.34,  0.00],
	'M': [0, 0.34, -0.34,  0.00],
	'D': [0, 0.34, -0.34, -1.00],
	'Q': [0, 0.34, -0.34,  0.00,  0.256, -0.256],
	'N': [0, 0.34, -0.34,  0.00,  0.256, -0.256],
	'P': [0, 0.34, -0.34,  0.00],
	'T': [0, 0.34, -0.34,  0.00,  0.153, -0.153],
	'K': [0, 0.34, -0.34,  0.00,  1.00],
	'R': [0, 0.34, -0.34,  0.00,  1.00],
	'E': [0, 0.34, -0.34,  0.00, -1.00],
	'W': [0, 0.34, -0.34,  0.00,  0.00, 0.138, -0.138, 0.00, 0.00],
	'F': [0, 0.34, -0.34,  0.00,  0.00, 0.00],
	'Y': [0, 0.34, -0.34,  0.00,  0.00,  0.00, 0.138, -0.138],
	'H': [0, 0.34, -0.34,  0.00,  0.00,  0.2, -0.2, 0.00,  0.136, -0.136],
	'Rsd': [0, 0.34, -0.34, 0, 0.34, -0.34, 0, 0.34, -0.34]
}
seq_to_resname = {
	'A': 'ALA',
	'G': 'GLY',
	'S': 'SER',
	'I': 'ILE',
	'L': 'LEU',
	'V': 'VAL',
	'C': 'CYS',
	'M': 'MET',
	'D': 'ASP',
	'Q': 'GLN',
	'N': 'ASN',
	'P': 'PRO',
	'T': 'THR',
	'K': 'LYS',
	'R': 'ARG',
	'E': 'GLU',
	'W': 'TRP',
	'F': 'PHE',
	'Y': 'TYR',
	'H': 'HIS',
	'Rsd': 'RSD'
}


# Dictionary reversal To allow compatinility to PDB files.
# resname_to_seq = {v: k for k, v in seq_to_resname.iteritems()}
resname_to_seq = {v: k for k, v in seq_to_resname.items()}

# try:
# 	u = Universe(input_txt)
# 	residue_names = u.select_atoms("protein").residues.resnames
# 	seq = ''.join([resname_to_seq[i] for i in residue_names])
# 	secstr = 0


# except:
# 	with open(input_txt, 'r') as f:
# 		seq = f.read()
# 		seq = seq.strip()
# 		secstr = 0

tableno = {'G': 1, 'A': 2, 'L': 3, 'I': 4, 'F': 5, 'R': 6, 'K': 7, 'D': 8, 'E': 9, 'N': 10,
		   'T': 11, 'S': 12, 'W': 13, 'Y': 14, 'V': 15, 'P': 16, 'M': 17, 'C': 18, 'H': 19, 'Q': 20}
tablenos1 = {'G': 21, 'A': 22, 'L': 23, 'I': 24, 'F': 25, 'R': 26, 'K': 27, 'D': 28, 'E': 29,
			 'N': 30, 'T': 31, 'S': 32, 'W': 33, 'Y': 34, 'V': 35, 'P': 36, 'M': 37, 'C': 38, 'H': 39, 'Q': 40}


for pepnum in range(1,n_peptides+1):

    out = []
    comment = [
        f';This is the itp file for the sequence {seq}',
        f'; Dihedral type: {bb_dih_type}, force constant: {bb_dih_fc}\n; Cysteine bridges active: {make_cys_bridge}',
    ]
    
    # if make_cys_bridge:
    #     comment.append('; Cysteine bridges active')

    
    # out.append(comment)
    out+=comment
    out.append('')
    ###############################################################################
    out.append('[moleculetype]')
    out.append(';molname		nrexcl')
    out.append('Protein' + str(pepnum) + '		 1')
    out.append('')

    ###############################################################################
    out.append('[atoms]')
    out.append(';id 		 type 		 resnr 		 residue 	atom 		 cgnr 		 charge ')
    aid = 1
    rid = 1
    charge_num = 1

    for aa in seq:
        # rid += 1
        resname = seq_to_resname[aa]
        beads = seq_to_bead[aa]
        types = seq_to_type[aa]
        cgnrs = seq_to_cgnr[aa]
        charge = seq_to_charge[aa]

        for (bead, tp, cgnr, cg) in zip(beads, types, cgnrs, charge):
            
            if (resname in ['LYS', 'ARG']) and (bead == 'S2'):
                # cb += 1
                out.append('%-17d%-16s%-16d%-14s%-17s%-17d%-6.4f' %
                        (aid, tp+f'_r{rid}p{pepnum}', rid, resname, bead, charge_num, cg))
                charge_num += cgnr
            elif (resname in ['PHE', 'TRP', 'TYR', 'HIS']) and (bead in ['S1', 'S2', 'S3', 'S4']):
                    # if (rid != arcount+1):
                    #     arb += 1
                    #     arcount = rid-1
                out.append('%-17d%-16s%-16d%-14s%-17s%-17d%-6.4f' %
                    (aid, tp+f'_r{rid}p{pepnum}', rid, resname, bead, charge_num, cg))

                charge_num += cgnr
            else:
                # arb += 0
                out.append('%-17d%-16s%-16d%-14s%-17s%-17d%-6.4f' %
                            (aid, tp, rid, resname, bead, charge_num, cg))
                charge_num += cgnr
            
            aid += 1

        rid+=1

    # Sort out PRO
    out_split = np.array([i.split() for i in out[8:]])
    df = pd.DataFrame(out_split)
    df2 = df.loc[(df[4] == 'BB')]
    df3 = np.array(df2.iloc[np.flatnonzero(df2[3] == 'PRO') - 1][0])
    for i in df3:
        df.at[int(i)-1, 1] = 'APP5'

    out_update = np.array(df)
    New = []
    for i in out_update:
        New.append('%-17s%-16s%-16s%-14s%-17s%-17s%-6s' %
                (i[0], i[1], i[2], i[3], i[4], i[5], i[6]))
    out[8:] = New

    ###############################################################################
    out.append('')
    out.append('[bonds]')
    out.append(';bonds between BB and dummies')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        bond_1 = backbone + 0  # BB
        bond_2 = backbone + 1  # BBp
        out.append('%-16d%-13d1				0.14		   5000' % (bond_1, bond_2))
        bond_2 = backbone + 2  # BBm
        out.append('%-16d%-13d1				0.14		   5000' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')


    out.append(';bonds between S1 and dummies')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['Q', 'T', 'N', 'S']:
            bond_1 = backbone + 3  # S1
            bond_2 = backbone + 4  # S1p
            out.append('%-16d%-13d1				0.14		   5000' % (bond_1, bond_2))
            bond_1 = backbone + 3  # S1
            bond_2 = backbone + 5  # S1m
            out.append('%-16d%-13d1				0.14		   5000' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';bonds between S2 and dummies')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['W', 'H']:
            bond_1 = backbone + 4  # S2
            bond_2 = backbone + 5  # S2p
            out.append('%-16d%-13d1				0.14		   5000' % (bond_1, bond_2))
            bond_2 = backbone + 6  # S2m
            out.append('%-16d%-13d1				0.14		   5000' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';bonds between S3 and dummies')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['Y']:
            bond_1 = backbone + 5  # S2
            bond_2 = backbone + 6  # S2p
            out.append('%-16d%-13d1				0.14		   5000' % (bond_1, bond_2))
            bond_2 = backbone + 7  # S2m
            out.append('%-16d%-13d1				0.14		   5000' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';bonds between S3 and dummies')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['H']:
            bond_1 = backbone + 7  # S2
            bond_2 = backbone + 8  # S2p
            out.append('%-16d%-13d1				0.14		   5000' % (bond_1, bond_2))
            bond_2 = backbone + 9  # S2m
            out.append('%-16d%-13d1				0.14		   5000' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')

    ################################################################################

    out.append(';bonds between S1 bead and BB bead')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa not in ['A', 'G']:
            bond_1 = backbone + 0  # BB
            bond_2 = backbone + 3  # S1
            out.append('%-16d%-13d1				0.25		   5000' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';bonds between S1 bead and S2 bead')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['K', 'R', 'E']:
            bond_1 = backbone + 3  # S1
            bond_2 = backbone + 4  # S2
            out.append('%-16d%-13d1				0.28		   5000' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';bonds between BB and BB beads')
    backbone = 1
    for i in range(0, len(seq)-1):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        bond_1 = backbone
        bond_2 = backbone + bead_num
        out.append('%-16d%-13d1				0.385		   7500' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')
    ###############################################################################
    out.append(';Bond between PRO-S1 and BB-Prev')
    backbone = 1
    for aa in seq:
        if aa in ['P']:
            bond_1 = backbone + 3  # S1
            bond_2 = backbone - prev  # S2
            out.append('%-16d%-13d1				0.385		   5000' % (bond_1, bond_2))
        backbone += len(seq_to_bead[aa])
        prev = len(seq_to_bead[aa])
    out.append('')
    ###############################################################################
    out.append(';Cys-Cys bridges: Cys_a-S1 to Cys_b-S1')
    if make_cys_bridge:
        backbone = 1
        bond_atoms = [] # atomids for the first and second cysteine S1 beads
        for aa in seq:
            if aa in ['C']:
                bond_atoms.append(backbone + 3)  # S1
                # bond_2 = backbone - prev  # S2
                # out.append('%-16d%-13d1				0.385		   5000' % (bond_1, bond_2))
            backbone += len(seq_to_bead[aa])
        out.append('%-16d%-13d1				0.385		   5000' % (bond_atoms[0], bond_atoms[1]))
        out.append('')
    else:
        out.append('')
    ###############################################################################
    out.append('[constraints]')

    '''
    out.append(';atom1/BB  	atom2/dummy 	func 		 distance');
    ###############################################################################
    out.append(';bonds between BB and dummies')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        print(backbone)
        bond_1 = backbone + 0 # BB
        bond_2 = backbone + 1 # BBp
        out.append('%-16d%-13d1				0.14' % (bond_1, bond_2))
        bond_2 = backbone + 2 # BBm
        out.append('%-16d%-13d1				0.14' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')
    '''


    out.append(';constraints between S1 and S2 beads for W Y F H')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['W', 'Y', 'F', 'H']:
            bond_1 = backbone + 3  # S1
            bond_2 = backbone + 4  # S2
            out.append('%-16d%-13d1				0.27' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';constraints between S2 and S3 beads for Y F')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['Y', 'F']:
            bond_1 = backbone + 4  # S2
            bond_2 = backbone + 5  # S3
            out.append('%-16d%-13d1				0.27' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';constraints between S1 and S3 beads for Y F')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['Y', 'F']:
            bond_1 = backbone + 3  # S1
            bond_2 = backbone + 5  # S3
            out.append('%-16d%-13d1				0.27' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';constraints between S1 and S3 beads for W')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['W']:
            bond_1 = backbone + 3  # S1
            bond_2 = backbone + 7  # S3
            out.append('%-16d%-13d1				0.27' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';constraints between S2 and S3 beads for W')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['W']:
            bond_1 = backbone + 4  # S2
            bond_2 = backbone + 7  # S3
            out.append('%-16d%-13d1				0.27' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';constraints between S2 and S4 beads for W')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['W']:
            bond_1 = backbone + 4  # S2
            bond_2 = backbone + 8  # S4
            out.append('%-16d%-13d1				0.27' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';constraints between S3 and S4 beads for W')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['W']:
            bond_1 = backbone + 7  # S3
            bond_2 = backbone + 8  # S4
            out.append('%-16d%-13d1				0.27' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')
    ######## HIS ADDITIONS ##########################
    out.append(';constraints between S1 and S3 beads for H')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['H']:
            bond_1 = backbone + 3  # S1
            bond_2 = backbone + 7  # S3
            out.append('%-16d%-13d1				0.27' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')
    out.append(';constraints between S2 and S3 beads for H')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['H']:
            bond_1 = backbone + 4  # S2
            bond_2 = backbone + 7  # S3
            out.append('%-16d%-13d1				0.27' % (bond_1, bond_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append('[angles]')
    out.append(';atom1/dummy  	atom2/BB 	atom3/dummy 	func 		 angle 	force.cons')
    ###############################################################################
    out.append(';angles between BB and dummies')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        angle_1 = backbone + 1  # BBp
        angle_2 = backbone + 0  # BB
        angle_3 = backbone + 2  # BBm
        out.append('%-16d%-16d%-16d2 		 180	  7.2' % (angle_1, angle_2, angle_3))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';angles between S1 and dummies')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['Q', 'T', 'N', 'S']:
            angle_1 = backbone + 4  # S1p
            angle_2 = backbone + 3  # S1
            angle_3 = backbone + 5  # S1m
            out.append('%-16d%-16d%-16d2 		 180	  7.2' %
                    (angle_1, angle_2, angle_3))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';angles between S2 and dummies - W H')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['W', 'H']:
            angle_1 = backbone + 5  # S2p
            angle_2 = backbone + 4  # S2
            angle_3 = backbone + 6  # S2m
            out.append('%-16d%-16d%-16d2 		 180	  7.2' %
                    (angle_1, angle_2, angle_3))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';angles between S3 and dummies - Y')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['Y']:
            angle_1 = backbone + 6  # S3p
            angle_2 = backbone + 5  # S3
            angle_3 = backbone + 7  # S3m
            out.append('%-16d%-16d%-16d2 		 180	  7.2' %
                    (angle_1, angle_2, angle_3))
        backbone += bead_num
    out.append('')

    out.append(';angles between S3 and dummies - H')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['H']:
            angle_1 = backbone + 8  # S3p
            angle_2 = backbone + 7  # S3
            angle_3 = backbone + 9  # S3m
            out.append('%-16d%-16d%-16d2 		 180	  7.2' %
                    (angle_1, angle_2, angle_3))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';angles between BB and BB beads')
    backbone = 1
    for i in range(0, len(seq)-2):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        aa_next = seq[i+1]
        bead_num_next = len(seq_to_bead[aa_next])
        angle_1 = backbone
        angle_2 = backbone + bead_num
        angle_3 = backbone + bead_num + bead_num_next
        if aa_next == 'P':
            out.append('%-16d%-16d%-16d2 		 98	  100' %
                    (angle_1, angle_2, angle_3))
        else:
            out.append('%-16d%-16d%-16d2 		 109	  75' %
                    (angle_1, angle_2, angle_3))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';angles between BB-BB-S1')
    backbone = 1
    for i in range(0, len(seq)-1):
        aa = seq[i]
        aa_next = seq[i+1]
        bead_num = len(seq_to_bead[aa])
        if aa_next not in ['A', 'G', 'P']:
            bead_num_next = len(seq_to_bead[aa_next])
            angle_1 = backbone
            angle_2 = backbone + bead_num
            angle_3 = backbone + bead_num + 3
            out.append('%-16d%-16d%-16d8   %-16d 	1' %
                    (angle_1, angle_2, angle_3, tablenos1[aa_next]))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';angles between BB\'-BB-S1')
    out.append('')

    ###############################################################################
    out.append(';angles between BB-S1-S2')
    backbone = 1
    for i in range(0, len(seq)):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        if (aa in ['E']):
            angle_1 = backbone + 0  # BB
            angle_2 = backbone + 3  # S1
            angle_3 = backbone + 4  # S2
            out.append('%-16d%-16d%-16d2			   180	   25.0' %
                    (angle_1, angle_2, angle_3))
        elif (aa in ['K', 'R']):  # 09/21/16 Silvina requested it from one of Sai's RVVGE itp files
            angle_1 = backbone + 0  # BB
            angle_2 = backbone + 3  # S1
            angle_3 = backbone + 4  # S2
            out.append('%-16d%-16d%-16d2			   151	   25.0' %
                    (angle_1, angle_2, angle_3))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';angles between BB-S1-S2 for F Y H')
    backbone = 1
    for i in range(0, len(seq)):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        if (aa in ['F', 'Y', 'H']):
            angle_1 = backbone + 0  # BB
            angle_2 = backbone + 3  # S1
            angle_3 = backbone + 4  # S2
            out.append('%-16d%-16d%-16d2			   150	   50.0' %
                    (angle_1, angle_2, angle_3))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';angles between BB-S1-S3 for F Y')
    backbone = 1
    for i in range(0, len(seq)):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        if (aa in ['F', 'Y']):
            angle_1 = backbone + 0  # BB
            angle_2 = backbone + 3  # S1
            angle_3 = backbone + 5  # S3
            out.append('%-16d%-16d%-16d2			   150	   50.0' %
                    (angle_1, angle_2, angle_3))
        backbone += bead_num
    out.append('')

    out.append(';angles between BB-S1-S3 for H')
    backbone = 1
    for i in range(0, len(seq)):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        if (aa in ['H']):
            angle_1 = backbone + 0  # BB
            angle_2 = backbone + 3  # S1
            angle_3 = backbone + 7  # S3
            out.append('%-16d%-16d%-16d2			   150	   50.0' %
                    (angle_1, angle_2, angle_3))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';angles between BB-S1-S2 for W')
    backbone = 1
    for i in range(0, len(seq)):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        if (aa in ['W']):
            angle_1 = backbone + 0  # BB
            angle_2 = backbone + 3  # S1
            angle_3 = backbone + 4  # S2
            out.append('%-16d%-16d%-16d2			   210	   50.0' %
                    (angle_1, angle_2, angle_3))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';angles between BB-S1-S3 for W')
    backbone = 1
    for i in range(0, len(seq)):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        if (aa in ['W']):
            angle_1 = backbone + 0  # BB
            angle_2 = backbone + 3  # S1
            angle_3 = backbone + 7  # S3
            out.append('%-16d%-16d%-16d2			   90	   50.0' %
                    (angle_1, angle_2, angle_3))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append('[dihedrals]')
    out.append('; i   j  k  l  funct  tableno  k')
    ###############################################################################

    '''
    Helix: Table_D0
    BetaSheet: Table_D1
    3-10: Table_D2  
    '''

    out.append(';tabulated dihedrals between BB and BB beads')
    backbone = 1
    for i in range(0, len(seq)-3):
        aa_0 = seq[i]
        bead_num_0 = len(seq_to_bead[aa_0])
        aa_1 = seq[i+1]
        bead_num_1 = len(seq_to_bead[aa_1])
        aa_2 = seq[i+2]
        bead_num_2 = len(seq_to_bead[aa_2])
        dih_1 = backbone + 0
        dih_2 = backbone + bead_num_0
        dih_3 = backbone + bead_num_0 + bead_num_1
        dih_4 = backbone + bead_num_0 + bead_num_1 + bead_num_2

        # out.append('%-16d%-16d%-16d%-16d8 		  6	 2' % (dih_1, dih_2, dih_3, dih_4))
        out.append(f'{dih_1} \t{dih_2} \t{dih_3} \t{dih_4} \t8 \t{bb_dih_type} \t{bb_dih_fc}')

        backbone += bead_num_0
    out.append('')

    ###############################################################################
    out.append(';Dihedral between BB-S2-S3-S1 for F Y')
    backbone = 1
    for i in range(0, len(seq)):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        if (aa in ['F', 'Y']):
            dih_1 = backbone + 0  # BB
            dih_2 = backbone + 4  # S2
            dih_3 = backbone + 5  # S3
            dih_4 = backbone + 3  # S1
            out.append('%-16d%-16d%-16d%-16d2 		  0	 50' %
                    (dih_1, dih_2, dih_3, dih_4))
        backbone += bead_num
    out.append('')


    out.append(';Dihedral between BB-S2-S3-S1 for H')
    backbone = 1
    for i in range(0, len(seq)):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        if (aa in ['H']):
            dih_1 = backbone + 0  # BB
            dih_2 = backbone + 4  # S2
            dih_3 = backbone + 7  # S3
            dih_4 = backbone + 3  # S1
            out.append('%-16d%-16d%-16d%-16d2 		  0	 50' %
                    (dih_1, dih_2, dih_3, dih_4))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';Dihedral between BB-S2-S3-S1 for W')
    backbone = 1
    for i in range(0, len(seq)):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        if (aa in ['W']):
            dih_1 = backbone + 0  # BB
            dih_2 = backbone + 4  # S2
            dih_3 = backbone + 7  # S3
            dih_4 = backbone + 3  # S1
            out.append('%-16d%-16d%-16d%-16d2 		  0	 50' %
                    (dih_1, dih_2, dih_3, dih_4))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';Dihedral between S1-S2-S4-S3 for W ... preserve planirity')
    backbone = 1
    for i in range(0, len(seq)):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        if (aa in ['W']):
            dih_1 = backbone + 0  # S1
            dih_2 = backbone + 4  # S2
            dih_3 = backbone + 8  # S4
            dih_4 = backbone + 7  # S3
            out.append('%-16d%-16d%-16d%-16d2 		  0	 200' %
                    (dih_1, dih_2, dih_3, dih_4))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append('[exclusions]')
    ###############################################################################
    out.append(';backbone dummy within bead')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        ex_1 = backbone + 1
        ex_2 = backbone + 2
        out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';backbone main site to main site 1-3')
    backbone = 1
    for i in range(0, len(seq)-2):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        aa_next = seq[i+1]
        bead_num_next = len(seq_to_bead[aa_next])
        ex_1 = backbone
        ex_2 = backbone + bead_num + bead_num_next
        out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')


    ###############################################################################
    out.append(';backbone main site to dummy 1-2')
    # Forward
    backbone = 1
    for i in range(0, len(seq)-1):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        ex_1 = backbone
        ex_2 = backbone + bead_num + 1
        out.append('%-10d%d' % (ex_1, ex_2))
        ex_1 = backbone
        ex_2 = backbone + bead_num + 2
        out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    # Backward
    backbone = 1
    for i in range(0, len(seq)-1):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        ex_1 = backbone + 1
        ex_2 = backbone + bead_num
        out.append('%-10d%d' % (ex_1, ex_2))
        ex_1 = backbone + 2
        ex_2 = backbone + bead_num
        out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';backbone dummy to dummy 1-2 ;replaced by tuning fudgeQQ')
    backbone = 1
    for i in range(0, len(seq)-1):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        ex_1 = backbone + 1
        ex_2 = backbone + bead_num + 1
        out.append('%-10d%d' % (ex_1, ex_2))
        ex_1 = backbone + 1
        ex_2 = backbone + bead_num + 2
        out.append('%-10d%d' % (ex_1, ex_2))
        ex_1 = backbone + 2
        ex_2 = backbone + bead_num + 1
        out.append('%-10d%d' % (ex_1, ex_2))
        ex_1 = backbone + 2
        ex_2 = backbone + bead_num + 2
        out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')

    ###############################################################################

    out.append(';backbone dummy to dummy 1-3')
    '''
    backbone = 1
    for i in xrange(0, len(seq)-2):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        aa_next = seq[i+1]
        bead_num_next = len(seq_to_bead[aa_next])
        ex_1 = backbone + 1
        ex_2 = backbone + bead_num + bead_num_next + 1
        out.append('%-10d%d' % (ex_1, ex_2))
        ex_1 = backbone + 1
        ex_2 = backbone + bead_num + bead_num_next + 2
        out.append('%-10d%d' % (ex_1, ex_2))
        ex_1 = backbone + 2
        ex_2 = backbone + bead_num + bead_num_next + 1
        out.append('%-10d%d' % (ex_1, ex_2))
        ex_1 = backbone + 2
        ex_2 = backbone + bead_num + bead_num_next + 2
        out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    '''
    out.append('')

    ###############################################################################
    out.append(';sidechain S1 to backbone dummy 1-2')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa not in ['A', 'G']:
            ex_1 = backbone + 3  # S1
            ex_2 = backbone + 1  # BBp
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 3  # S1
            ex_2 = backbone + 2  # BBm
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';sidechain S1 main bead to backbone dummy 1-3 (Optional)')
    # Forward
    backbone = 1
    for i in range(0, len(seq)-1):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        if aa in ['D']:
            ex_1 = backbone + 3  # S1
            ex_2 = backbone + bead_num + 1  # next BBp
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 3  # S1
            ex_2 = backbone + bead_num + 2  # next BBm
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    # Backward
    backbone = 1
    for i in range(0, len(seq)-1):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        aa_next = seq[i+1]
        if aa_next in ['D']:
            ex_1 = backbone + 1  # BBp
            ex_2 = backbone + bead_num + 3  # next S1
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 2  # BBm
            ex_2 = backbone + bead_num + 3  # next S1
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';side chain S1 dummy within bead')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['Q', 'T', 'N', 'S']:
            ex_1 = backbone + 4  # S1p
            ex_2 = backbone + 5  # S1m
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';side chain S1 dummy to backbone bead')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['Q', 'T', 'N', 'S']:
            ex_1 = backbone + 4  # S1p
            ex_2 = backbone + 0  # BB
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 5  # S1m
            ex_2 = backbone + 0  # BB
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';side chain S1 dummy to backbone dummy 1-2 ;replaced by tuning fudgeQQ')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['Q', 'T', 'N', 'S']:
            ex_1 = backbone + 4  # S1p
            ex_2 = backbone + 1  # BBp
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 4  # S1p
            ex_2 = backbone + 2  # BBm
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 5  # S1m
            ex_2 = backbone + 1  # BBp
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 5  # S1m
            ex_2 = backbone + 2  # BBm
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';sidechain S1 dummy to backbone dummy 1-3')
    '''
    # Forward
    backbone = 1
    for i in xrange(0, len(seq)-1):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        if aa in ['Q', 'P', 'T', 'N']:
            ex_1 = backbone + 4 # S1p
            ex_2 = backbone + bead_num + 1 # next BBp
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 4 # S1p
            ex_2 = backbone + bead_num + 2 # next BBm
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 5 # S1m
            ex_2 = backbone + bead_num + 1 # next BBp
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 5 # S1m
            ex_2 = backbone + bead_num + 2 # next BBm
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    # Backward
    backbone = 1
    for i in xrange(0, len(seq)-1):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        aa_next = seq[i+1]
        if aa_next in ['Q', 'P', 'T', 'N']:
            ex_1 = backbone + 1 # BBp
            ex_2 = backbone + bead_num + 4 # next S1p
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 2 # BBm
            ex_2 = backbone + bead_num + 4 # next S1p
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 1 # BBp
            ex_2 = backbone + bead_num + 5 # next S1m
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 2 # BBm
            ex_2 = backbone + bead_num + 5 # next S1m
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    '''
    out.append('')

    ###############################################################################
    out.append(';sidechain backbone S2 to sidechain S1 dummy')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if (aa in ['Rsd']):
            ex_1 = backbone + 6  # S2
            ex_2 = backbone + 4  # S1p
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 6  # S2
            ex_2 = backbone + 5  # S1m
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')


    ###############################################################################
    out.append(';sidechain S2 to backbone BB 1-3')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if (aa in ['K', 'R', 'E']):
            ex_1 = backbone + 0  # BB
            ex_2 = backbone + 4  # S2
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';sidechain S2 main bead to backbone dummy 1-3')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if (aa in ['K', 'R', 'E']):
            ex_1 = backbone + 4
            ex_2 = backbone + 1
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 4
            ex_2 = backbone + 2
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')

    ###############################################################################
    out.append(';Specific Exclusions for Disulfide/dicysteine bridges')
    out.append('; Scheme: BB--S1--S1\'--BB\'')
    out.append('; CYS bridge BB<--x-->S1\' and BB\'<--x-->S1 exclusions 1-3')
    backbone = 1
    if make_cys_bridge:
        BB_exclusion = []
        S1_exclusion = []
        for aa in seq:
            """
            WARNING: If multiple cys-bridges exist, this will NOT WORK.
            """
            bead_num = len(seq_to_bead[aa])
            
            if aa == 'C':
                BB_exclusion.append(backbone + 0)  # BB
                S1_exclusion.append(backbone + 3)  # S1
                # out.append('%-10d%d' % (ex_1, ex_2))
            backbone += bead_num
        
        # AGAIN, only assumes 1 cys bridge
        out.append(f"{BB_exclusion[0]}\t{S1_exclusion[1]}   ; BB<--x-->S1\'")
        out.append(f"{BB_exclusion[1]}\t{S1_exclusion[0]}   ; BB\'<--x-->S1")
    out.append('')
    out.append('; CYS bridge BB dummies<--x-->S1\' and BB\' dummies<--x-->S1 exclusions 1-3')
    backbone = 1
    if make_cys_bridge:
        BBm_exclusion = []
        BBp_exclusion = []
        S1_exclusion = []
        # print(S1_exclusion)
        for aa in seq:
            """
            WARNING: If multiple cys-bridges exist, this will NOT WORK.
            """
            bead_num = len(seq_to_bead[aa])
            
            if aa == 'C':
                BBp_exclusion.append(backbone + 1)  # BBp
                BBm_exclusion.append(backbone + 2)  # BBm
                S1_exclusion.append(backbone + 3)  # S1
                # out.append('%-10d%d' % (ex_1, ex_2))
            backbone += bead_num
        # print(S1_exclusion)
        # AGAIN, only assumes 1 cys bridge
        out.append(f"{BBp_exclusion[0]}\t{S1_exclusion[1]}   ; BBp<--x-->S1\'")
        out.append(f"{BBm_exclusion[0]}\t{S1_exclusion[1]}   ; BBm<--x-->S1\'")
        out.append(f"{BBp_exclusion[1]}\t{S1_exclusion[0]}   ; BBp\'<--x-->S1")
        out.append(f"{BBm_exclusion[1]}\t{S1_exclusion[0]}   ; BBm\'<--x-->S1")
    out.append('')


    ###############################################################################
    out.append(';Specific Exclusions for aromatic residues')

    out.append('')
    out.append(';For F')
    out.append(';BB to sidechain main site S2')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['F']:
            ex_1 = backbone + 0  # BB
            ex_2 = backbone + 4  # S2
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')
    out.append(';BB to sidechain main site S3')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['F']:
            ex_1 = backbone + 0  # BB
            ex_2 = backbone + 5  # S3
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')


    ###############################################################################
    out.append(';For Y')
    out.append(';side chain S3 dummy within bead')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['Y']:
            ex_1 = backbone + 6  # S3p
            ex_2 = backbone + 7  # S3m
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append(';side chain S3 dummy to sidechain main site S1')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['Y']:
            ex_1 = backbone + 6  # S3p
            ex_2 = backbone + 3  # S1
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 7  # S3m
            ex_2 = backbone + 3  # S1
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')
    out.append(';side chain S3 dummy to sidechain main site S2')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['Y']:
            ex_1 = backbone + 6  # S3p
            ex_2 = backbone + 4  # S2
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 7  # S3m
            ex_2 = backbone + 4  # S2
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')
    out.append(';side chain S3 dummy to sidechain main site S1')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['Y']:
            ex_1 = backbone + 6  # S3p
            ex_2 = backbone + 3  # S1
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 7  # S3m
            ex_2 = backbone + 3  # S1
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')
    out.append(';BB to sidechain main site S2')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['Y']:
            ex_1 = backbone + 0  # BB
            ex_2 = backbone + 4  # S2
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')
    out.append(';BB to sidechain main site S3')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['Y']:
            ex_1 = backbone + 0  # BB
            ex_2 = backbone + 5  # S3
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')
    ###############################################################################
    out.append(';For W')
    out.append(';side chain S2 dummy within bead')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['W']:
            ex_1 = backbone + 5  # S2p
            ex_2 = backbone + 6  # S2m
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append(';side chain S2 dummy to sidechain main site S1')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['W']:
            ex_1 = backbone + 5  # S2p
            ex_2 = backbone + 3  # S1
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 6  # S2m
            ex_2 = backbone + 3  # S1
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')
    out.append(';side chain S2 dummy to sidechain main site S3')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['W']:
            ex_1 = backbone + 5  # S2p
            ex_2 = backbone + 7  # S3
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 6  # S2m
            ex_2 = backbone + 7  # S3
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')
    out.append(';side chain S2 dummy to sidechain main site S4')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['W']:
            ex_1 = backbone + 5  # S2p
            ex_2 = backbone + 8  # S4
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 6  # S2m
            ex_2 = backbone + 8  # S4
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')
    out.append(';BB to sidechain main site S2')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['W']:
            ex_1 = backbone + 0  # BB
            ex_2 = backbone + 4  # S2
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')
    out.append(';BB to sidechain main site S3')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['W']:
            ex_1 = backbone + 0  # BB
            ex_2 = backbone + 7  # S3
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')
    out.append(';S4 to sidechain main site S1')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['W']:
            ex_1 = backbone + 3  # S1
            ex_2 = backbone + 8  # S4
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')


    ###############################################################################
    out.append(';For H')
    out.append(';side chain S2 dummy within bead')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['H']:
            ex_1 = backbone + 5  # S2p
            ex_2 = backbone + 6  # S2m
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num

    out.append(';For H')
    out.append(';side chain S3 dummy within bead')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['H']:
            ex_1 = backbone + 8  # S3p
            ex_2 = backbone + 9  # S3m
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num

    out.append(';side chain S2 dummy to sidechain main site S1')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['H']:
            ex_1 = backbone + 5  # S2p
            ex_2 = backbone + 3  # S1
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 6  # S2m
            ex_2 = backbone + 3  # S1
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')
    out.append(';side chain S3 dummy to sidechain main site S1')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['H']:
            ex_1 = backbone + 8  # S3p
            ex_2 = backbone + 3  # S1
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 9  # S3m
            ex_2 = backbone + 3  # S1
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')
    out.append(';side chain S2 dummy to sidechain main site S3')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['H']:
            ex_1 = backbone + 5  # S2p
            ex_2 = backbone + 7  # S3
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 6  # S2m
            ex_2 = backbone + 7  # S3
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')
    out.append(';side chain S2 dummy to sidechain main site S3-dummies')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['H']:
            ex_1 = backbone + 5  # S2p
            ex_2 = backbone + 8  # S3p
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 5  # S2p
            ex_2 = backbone + 9  # S3m
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num

    out.append('')
    out.append(';side chain S2 dummy to sidechain main site S3-dummies')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['H']:
            ex_1 = backbone + 6  # S2m
            ex_2 = backbone + 8  # S3p
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 6  # S2m
            ex_2 = backbone + 9  # S3m
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')

    out.append(';side chain S2 dummy to sidechain main site BB-dummies')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['H']:
            ex_1 = backbone + 5  # S2p
            ex_2 = backbone + 1  # BBp
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 5  # S2p
            ex_2 = backbone + 2  # BBm
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num

    out.append('')
    out.append(';side chain S2 dummy to sidechain main site BB-dummies')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['H']:
            ex_1 = backbone + 6  # S2m
            ex_2 = backbone + 1  # BBp
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 6  # S2m
            ex_2 = backbone + 2  # BBm
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')
    out.append(';side chain S3 dummy to sidechain main site BB-dummies')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['H']:
            ex_1 = backbone + 8  # S3p
            ex_2 = backbone + 1  # BBp
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 8  # S3p
            ex_2 = backbone + 2  # BBm
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num

    out.append('')
    out.append(';side chain S3 dummy to sidechain main site BB-dummies')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['H']:
            ex_1 = backbone + 9  # S3m
            ex_2 = backbone + 1  # BBp
            out.append('%-10d%d' % (ex_1, ex_2))
            ex_1 = backbone + 9  # S3m
            ex_2 = backbone + 2  # BBm
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')

    out.append(';BB to sidechain main site S2')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['H']:
            ex_1 = backbone + 0  # BB
            ex_2 = backbone + 4  # S2
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')
    out.append(';BB to sidechain main site S3')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['H']:
            ex_1 = backbone + 0  # BB
            ex_2 = backbone + 7  # S3
            out.append('%-10d%d' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')


    ###############################################################################
    out.append('[pairs]')
    ###############################################################################
    out.append(';1-2 dummy-dummy interactions')
    ###############################################################################
    out.append(';backbone dummy to dummy 1-2')
    backbone = 1
    for i in range(0, len(seq)-1):
        aa = seq[i]
        bead_num = len(seq_to_bead[aa])
        ex_1 = backbone + 1
        ex_2 = backbone + bead_num + 1
        out.append('%-10d%-10d1         0.0          0.0' % (ex_1, ex_2))
        ex_1 = backbone + 1
        ex_2 = backbone + bead_num + 2
        out.append('%-10d%-10d1         0.0          0.0' % (ex_1, ex_2))
        ex_1 = backbone + 2
        ex_2 = backbone + bead_num + 1
        out.append('%-10d%-10d1         0.0          0.0' % (ex_1, ex_2))
        ex_1 = backbone + 2
        ex_2 = backbone + bead_num + 2
        out.append('%-10d%-10d1         0.0          0.0' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')
    ###############################################################################
    out.append(';side chain S1 dummy to backbone dummy 1-2 ;replaced by tuning fudgeQQ')
    backbone = 1
    for aa in seq:
        bead_num = len(seq_to_bead[aa])
        if aa in ['Q', 'T', 'N', 'S']:
            ex_1 = backbone + 4  # S1p
            ex_2 = backbone + 1  # BBp
            out.append('%-10d%-10d1         0.0          0.0' % (ex_1, ex_2))
            ex_1 = backbone + 4  # S1p
            ex_2 = backbone + 2  # BBm
            out.append('%-10d%-10d1         0.0          0.0' % (ex_1, ex_2))
            ex_1 = backbone + 5  # S1m
            ex_2 = backbone + 1  # BBp
            out.append('%-10d%-10d1         0.0          0.0' % (ex_1, ex_2))
            ex_1 = backbone + 5  # S1m
            ex_2 = backbone + 2  # BBm
            out.append('%-10d%-10d1         0.0          0.0' % (ex_1, ex_2))
        backbone += bead_num
    out.append('')


    outdata = '\n'.join(out)
    with open(f"{out_dir}/pep{pepnum}.itp", 'w') as f:
        f.write(outdata)
