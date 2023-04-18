__author__ = "Marco Di Gennaro"
__email__ = "marco.di.gennaro@external.toyota-europe.com"

import os, re, sys
import numpy as np
import pickle as pck
import pandas as pd
import ast

import matplotlib.pyplot as plt
import seaborn as sns

from mpl_toolkits.mplot3d import Axes3D
from itertools import product, cycle

from mendeleev import element
import math

import scipy.linalg as LA
from scipy.optimize import least_squares

from ase import Atoms

#from converter import Converter
import subprocess as sp
import scipy.constants as const
import warnings
import getpass
import copy


Ha2eV = const.value('hartree-electron volt relationship') #27.211
Ang2Bohr = 1.8897259886
Ang2m = 1e-10
F = const.value('Faraday constant')
elem_charge = const.value(u'elementary charge') #1.6021766208e-19
eps0 = const.value('electric constant') #8.854187817620389e-12

joule_to_ev = const.value('joule-electron volt relationship') #6.241509126e+18
ha_to_joule = const.value('hartree-joule relationship') #4.35974434e-18 J

def last_slurm_id(folder):
    return max([ int(f.replace('slurm-','').replace('.out','')) for f in os.listdir(folder) if f.startswith('slurm') ])

def find_last_log(folder):
    if os.path.exists( folder ):
       all_logs = [ f for f in os.listdir(folder) if f.startswith('log.') ]
       if len(all_logs) == 0:
          print( 'NO out file yet ', folder )
          return False
       else:
          if len(all_logs) == 1:
             return all_logs[0]
          else:
             print_tab( 1,  'WARNING: {} log files found in {}'.format(len(all_logs), folder) )
             last_slurm = last_slurm_id(folder)
             return [ log for log in all_logs if log.endswith(str(last_slurm)) ][0]
             #warnings.warn('{} log files found in {}'.format(len(all_logs), folder.replace('/data/mdi0316/WORK/','')))
             #print_tab( 2, [ 'all  logs = ' ] + all_logs +
             #              [  'last logs = {}'.format( last_log ) ] )
    else:
      return False

def running_jobs():
    user = 'mdi0316' #getpass.getuser()
    running_jobs = sp.getoutput('squeue -u {} -o "%i %.150j"'.format(user))
    if len(running_jobs) > 1:
      all_ids    = [ item.split()[0] for item in running_jobs.split('\n')[1:] ]
      all_labels = [ item.split()[1] for item in running_jobs.split('\n')[1:] ]
      #all_status = [ item.split()[2] for item in running_jobs.split('\n')[1:] ]
      #all_nodes  = [ item.split()[2] for item in running_jobs.split('\n')[1:] ]
      return all_ids, all_labels 
    else:
      return [],[]

def running_label( label ):
    run_idx, run_labels = running_jobs()
    running = False
    for tmp_lab in run_labels:
        if label in tmp_lab:
           running = True
           break
    return( running )

def print_tab( ntab, msg_to_print ):
    if   isinstance(msg_to_print,str):
         print( '{}{}'.format(ntab*'\t', msg_to_print) )
    elif isinstance(msg_to_print,int):
         print( '{}{}'.format(ntab*'\t', msg_to_print) )
    elif isinstance(msg_to_print,float):
         print( '{}{}'.format(ntab*'\t', msg_to_print) )
    elif isinstance(msg_to_print,list):
         for line in msg_to_print:
             print( '{}{}'.format(ntab*'\t', line) )
    elif isinstance(msg_to_print,dict):
         for kk, vv in msg_to_print.items():
             print( '{}{}: {}'.format(ntab*'\t', kk, vv) )
    elif isinstance(msg_to_print,pd.core.series.Series):
         print( msg_to_print )
    elif isinstance(msg_to_print,pd.DataFrame):
         for row, line in msg_to_print.iterrows():
             print( '{}{}'.format(ntab*'\t', line.values) )
    else:
         print( 'Unrecognized object: {} '.format(type(msg_to_print) ))
         print( wrg )

def cart2polar( x,y,z ):
    r     = np.sqrt(x**2 + y**2 + z**2)
    theta = math.acos(z/r)
    phi   = math.atan2(y,x)
    return( r, phi, theta )

def polar2cart( r, phi, theta ):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return( x,y,z )

def perpendicular_vector(v):
    if v[1] == 0 and v[2] == 0:
        if v[0] == 0:
            raise ValueError('zero vector')
        else:
            return np.cross(v, [0, 1, 0])
    return np.cross(v, [1, 0, 0])

def perpendicular( a ):
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

def df_to_repres( read_dir, file_name, energy_df, show_distr = False ):
    X_npy, Y_npy = '{}/{}_X.npy'.format(read_dir,file_name), '{}/{}_Y.npy'.format(read_dir, file_name)

    if os.path.exists( X_npy ) and os.path.exists( Y_npy ):
       print_tab( 0, ['Representation files exists:', X_npy, Y_npy] )
       X = np.load( X_npy )
       Y = np.load( Y_npy )
    else:
       print_tab( 0, ['Writing representation files:', X_npy, Y_npy] )
       X = np.empty((len(energy_df),6))
       Y = np.empty((len(energy_df),1))
       for row, index in energy_df.iterrows():
           ns1 = re.sub( '\[\s+',  '[', index['ZGP'] )
           ns2 = re.sub(    '\n',  ' ', ns1 )
           ns3 = re.sub(   '\s+', ', ', ns2 )
           chem_array = np.array( ast.literal_eval( ns3 ) )
           ns1 = re.sub( '\[\s+',  '[', index['cart'] )
           ns2 = re.sub(    '\n',  ' ', ns1 )
           ns3 = re.sub(   '\s+', ', ', ns2 )
           coords_array = np.array( ast.literal_eval( ns3 ) )
           for chem, cart in zip( chem_array, coords_array):
               Z,G,P = chem
               x,y,z = cart
               elem = element(int(Z))
               tmp_array = np.array( [ Z, elem.period, elem.group_id, x, y, z ] )
           X[row] = tmp_array
           Y[row] = index['D.ener.']
       np.save( X_npy, X )
       np.save( Y_npy, Y )

    Y = np.array( [ y-min(Y) for y in Y ] )

    if show_distr:
       fig = plt.figure()
       ax1 = fig.add_subplot()
       sns.distplot(Y, ax=ax1)
       plt.ylabel('counts')
       plt.xlabel('energy')
       plt.show()
       plt.close()

    return( X, Y )


def rotation_matrix(axis, angle):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by angle radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(angle / 2.0)
    b, c, d = -axis * math.sin(angle / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def average(lst): 
    return sum(lst) / len(lst) 

def rotate_molecule( cc_array, idx0=0, idx1=1, idx2=2 ):
    """
    Input: np cartesian coordinates
    Returns: rotated matrix:
      cc_array[idx0] at the origin,
      cc_array[idx1] on the y axis
      cc_array[idx2] on the xy plane
    """
    ax_x, ax_y, ax_z = [1,0,0], [0,1,0], [0,0,1]

    # translation: at.0 => origin
    cc_array -= cc_array[idx0]

    # rotation: at.1 => xz plane
    x1, y1, z1 = cc_array[idx1]
    ang1 = np.arctan( x1/y1 )
    rot_mat1 = rotation_matrix(ax_z, ang1)
    cc_array = np.dot( rot_mat1, cc_array.T ).T
    
    # rotation: at.1 => x axis
    x1, y1, z1 = cc_array[idx1]
    ang2 = -np.arctan( z1/y1 )
    rot_mat2 = rotation_matrix(ax_x,ang2)
    cc_array = np.dot( rot_mat2, cc_array.T ).T

    # rotation: at.2 => xy plane
    x2, y2, z2 = cc_array[idx2]
    ang3 = np.arctan( z2/x2 )
    rot_mat3 = rotation_matrix(ax_y, ang3)
    cc_array = np.dot( rot_mat3, cc_array.T ).T
  
    # rotation: at.1 => positive y
    x1, y1, z1 = cc_array[idx1]
    if y1 < 0:
       rot_mat4 = rotation_matrix(ax_z, math.pi)
       cc_array = np.dot( rot_mat4, cc_array.T ).T

    # rotation: at.2 => positive x  
    x2, y2, z2 = cc_array[idx2]
    if x2 < 0:
       rot_mat5 = rotation_matrix(ax_y, math.pi)
       cc_array = np.dot( rot_mat5, cc_array.T ).T

    return( cc_array )

def rotate_ase_object( ase_object, print_out = False ):
    """
    Input: ase object
    Returns: rotated ase object with:
      atom1 at the origin,
      atom2 on the y axis
      atom3 on the xy plane
    Options:
      print cartesian and Zmat coordinates
      plot 3D figure
    """
    # translate
    coord = ase_object.get_positions()
    ase_object.translate(-coord[0])
    transl_coord = ase_object.get_positions()
    # rotate
    atom1_coord = transl_coord[1]
    radius, phi, theta = cart2polar( *atom1_coord )
    Delta_theta = math.pi/2 - theta
    atom1_xy_proj = atom1_coord[:2]/LA.norm( atom1_coord[:2] )
    atom1_xy_perp = np.array( [ *perpendicular( atom1_xy_proj ), 0 ] )
    ase_object.rotate( math.degrees(Delta_theta), atom1_xy_perp )
    # rotate
    rot1_coord   = ase_object.get_positions()
    atom1_coord  = rot1_coord[1]
    atom1_y_proj = atom1_coord[:2]
    rot_angle    = math.atan2( atom1_coord[0], atom1_coord[1] )
    ase_object.rotate( math.degrees(rot_angle), np.array([0,0,1]) )
    # rotate
    rot2_coord  = ase_object.get_positions()
    atom2_coord = rot2_coord[2]
    rot_angle   = math.atan2( atom2_coord[2], atom2_coord[0] )
    ase_object.rotate( math.degrees(rot_angle), np.array([0,1,0]) )
    rot3_coord  = ase_object.get_positions()

    if print_out:
       cart_coord_file = os.path.join( 'cart_ase_object.dat' )
       zmat_coord_file = os.path.join( 'zmat_ase_object.dat' )
       with open( cart_coord_file, 'w+' ) as f:
            f.write('24\n\n')
            for line in X_array:
                Z,g,p,x,y,z = line
                elem = element( int(Z) )
                new_line = '{:<3} {:<1.16e} {:<1.16e} {:1.16e}\n'.format(elem.symbol,x,y,z)
                if len(new_line) > 79:
                   print( Warning_GAMESS_line_too_long )
                f.write(new_line)

       a = Converter()
       b = a.read_cartesian( input_file = cart_coord_file )
       a.run_cartesian(      input_file = cart_coord_file, output_file = zmat_coord_file )

    return( ase_object )

def ase_to_ZGPxyz( ase_object ):
    """
    Input: ase object
    Output: np array each line being = [Z,G,P,x,y,z]
    """
    atomic_nums = ase_object.get_atomic_numbers()
    coordinates = ase_object.get_positions()
    ZPG_xyz = np.empty((len(coordinates),6))
    for pos, Z, coords in zip( range(len(coordinates)), atomic_nums, coordinates ):
        elem = element( int(Z) )
        x, y, z = coordinates[pos]
        ZPG_xyz[pos] = np.array( [ Z, elem.period, elem.group_id, x, y, z ] )

    return( ZPG_xyz )

def center_of_mass( cart_coords, length_factor = 1. ):
    com = [0., 0., 0.]
    tot_weight = 0.
    for kk, coord in cart_coords.items():
        elem = element( coord['elem'].upper() )
        com[0] += elem.atomic_weight * length_factor * float(coord['x'])
        com[1] += elem.atomic_weight * length_factor * float(coord['y'])
        com[2] += elem.atomic_weight * length_factor * float(coord['z'])
        tot_weight += elem.atomic_weight
    com = [ xx/tot_weight for xx in com ]
    return( list(com) )



