# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 15:45:43 2021

@author: MDI0316
"""
import os, sys
from pathlib import Path
from common_functions import src_dir

proj_name = 'MOFS_QML'
tml_runs_dir = Path( 'C:', '/Users', 'MDI0316', 'Desktop', 'Turbomol_runs',  proj_name)
ref_tml_runs_dir = os.path.join( tml_runs_dir, 'References', 'completed' )
qml_tml_runs_dir = os.path.join( tml_runs_dir, 'candidates', 'completed' )

proj_src_dir = os.path.join( src_dir, proj_name )

if proj_src_dir not in sys.path:
    sys.path.insert(0, proj_src_dir)

import pandas as pd
import ase
from ase import Atoms
from ase.build import molecule
from ase.io import read
from ase.visualize import view
import ase.io.png
from IPython.display import Image
from pysmiles import read_smiles
import networkx as nx
import matplotlib.pyplot as plt
import shutil
import math
#!pip install imolecule 
# or
#!conda install -c patrickfuller imolecule
import imolecule

Ha2kjmol=2600
Ang2Bohr=1.8897161646320724

import PARSE_FILES as PF
import GAMESS as GMS
from known_smiles import aromatic_ring_smiles, rejected_smiles

def read_tmol_coords( tmol_file ):
    if os.path.exists(tmol_file):
        xyz_file = Tmol_dft_crd_file + '.xyz'
        if not os.path.exists(xyz_file):
            print('reading: ', tmol_file)
            with open(tmol_file, 'r') as f:
                coords = f.read().split('$')[1].split('coord')[1]
            coords_list = coords.strip().split('\n')
            print('writing: ', xyz_file)
            with open(xyz_file, 'w+') as f:
                f.write('{}\n\n'.format(len(coords_list)))
                for crd in coords_list:
                    x,y,z,e = crd.split()
                    f.write('{}\t{}\t{}\t{}\t\n'.format(e,float(x)/Ang2Bohr,float(y)/Ang2Bohr,float(z)/Ang2Bohr))
    else:
        print( 'missing file: ', tmol_file )

def read_tmol_dft_ene( tmol_file ):
    if os.path.exists(tmol_file):
        print('reading: ', tmol_file)
        with open(tmol_file, 'r') as f:
            lines = f.readlines()
        lines = [f for f in lines if f != '' and not f.strip().startswith('$') ]
        return float(lines[-1].split()[1])
    else:
        print( 'missing file: ', tmol_file )
        return math.nan
    
def read_tmol_mp2_ene( tmol_file ):
    if os.path.exists(tmol_file):
        print('reading: ', tmol_file)
        with open(tmol_file, 'r') as f:
            lines = f.readlines()
        lines = [f for f in lines if f != '' and not f.strip().startswith('$') ]
        return float(lines[-1].split()[1]) + float(lines[-1].split()[4]) 
    else:
        print( 'missing file: ', tmol_file )
        return math.nan

def xyz_2_coord( xyz_file, coord_file ):
    if os.path.exists(xyz_file):
        print('reading: ', xyz_file)
        with open(xyz_file, 'r') as f:
            lines = f.readlines()
        with open(coord_file, 'w+') as f:
            f.write('$coords\n')
            for line in lines[2:]:
                e,x,y,z = line.split()
                f.write('{}\t{}\t{}\t{}\n'.format(x,y,z,e))
            f.write('$ends\n')
        print( lines )
        print( sfdpfk ) 
    else:
        print( 'missing file: ', xyz_file )

def dict_to_xyz(xyz_dict, xyz_file, comment=''):
    with open( xyz_file, 'w+' ) as f:
        f.write('{}\n'.format(len(xyz_dict)))
        f.write('{}\n'.format(comment))
        for k,v in xyz_dict.items():
            f.write('{} {} {} {}\n'.format(v['elem'], *v['coords'] ))
    return

def coord_2_xyz( coord_file, xyz_file ):
    if os.path.exists(coord_file):
        if not os.path.exists(xyz_file):
            print('reading: ', coord_file)
            with open(coord_file, 'r') as f:
                coords = f.read().split('$')[1].split('coord')[1]
            coords_list = coords.strip().split('\n')
            print('writing: ', xyz_file)
            with open(xyz_file, 'w+') as f:
                f.write('{}\n\n'.format(len(coords_list)))
                for crd in coords_list:
                    x,y,z,e = crd.split()
                    f.write('{}\t{}\t{}\t{}\t\n'.format(e,float(x)/Ang2Bohr,float(y)/Ang2Bohr,float(z)/Ang2Bohr))
    else:
        print( 'missing file: ', coord_file )

#PLOT WITH NETWORKX
def plot_networx(qml_label, df):
    gamess = df.loc[df[('qml_label','0','0')] == qml_label ].head(1)
    mol = read_smiles(fin_smi)
    elements = nx.get_node_attributes(mol, name = "element")
    nx.draw(mol, with_labels=True, labels = elements, pos=nx.spring_layout(mol))
    plt.gca().set_aspect('equal')
    plt.show()

#plot with ASE
def plot_ase(qml_label, df):
    gamess = df.loc[df[('qml_label','0','0')] == qml_label ]
    conf_label = gamess[('H2_conf_label', '0', '0')].values[0]
    atoms = read('XYZ/closed_ring/{}_{}.xyz'.format(qml_label,conf_label))
    Image(filename='PNG/closed_ring/{}_{}.png'.format(qml_label,conf_label)) 
    
def plot_gamess_structure(qml_label, df):
    gamess = df.loc[ final_df[('qml_label','0','0') ] == qml_label ]
    conf_label, H2_conf_label = gamess[[('conf_label','0','0'),('H2_conf_label','0','0')]].values[0]
    
    GMS_conf_xyz = os.path.join( 'XYZ', 'closed_ring', '{}_{}.xyz'.format(qml_label, conf_label))
    GMS_H2_conf_xyz = os.path.join( 'XYZ', 'closed_ring', '{}_{}.xyz'.format(qml_label, H2_conf_label))
    
    imolecule.draw(GMS_H2_conf_xyz)
    
def plot_turbomole_structure(qml_label, df):
    gamess = df.loc[ df[('qml_label','0','0') ] == qml_label ]
    conf_label, H2_conf_label = gamess[[('conf_label','0','0'),('H2_conf_label','0','0')]].values[0]

    TML_conf = os.path.join( qml_tml_runs_dir, '{}_{}'.format(qml_label, conf_label), 'coord' )
    TML_H2_conf = os.path.join( qml_tml_runs_dir, '{}_{}'.format(qml_label, H2_conf_label), 'coord' )

    TML_conf_xyz = os.path.join( qml_tml_runs_dir, '{}_{}'.format(qml_label, conf_label), 'coord.xyz' )
    TML_H2_conf_xyz = os.path.join( qml_tml_runs_dir, '{}_{}'.format(qml_label, H2_conf_label), 'coord.xyz' )
    
    for tmol, xyz in zip( [TML_conf, TML_H2_conf], [TML_conf_xyz, TML_H2_conf_xyz] ):
        if not os.path.exists(xyz):
            tmol_2_xyz( tmol, xyz )
        
    imolecule.draw(TML_H2_conf_xyz)