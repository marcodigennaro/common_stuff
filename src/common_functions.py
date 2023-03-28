
### common input start
import os, sys, re
import numpy as np
import pandas as pd
import shutil
import subprocess as sp
import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt
import time
import math
from pathlib import Path
from matplotlib.lines import Line2D

from numpy import linalg as LA
import scipy.constants as const
from mpl_toolkits.mplot3d import axes3d

desktop_dir = Path( 'C:', '/Users', 'MDI0316', 'Desktop' )
data_dir = Path( 'C:', '/Users', 'MDI0316', 'Desktop', 'DATA' )
codes_dir = Path( 'C:', '/Users', 'MDI0316', 'Desktop', 'CODES' )
work_dir = Path( 'C:', '/Users', 'MDI0316', 'OneDrive - Toyota Motor Europe', 'WORK' )
src_dir = Path( work_dir, 'src' )

import datetime
now = datetime.datetime.now()

month_num = now.month
datetime_object = datetime.datetime.strptime(str(month_num), "%m")
month_str = datetime_object.strftime("%b").upper()

month_label = '{:02d}_{}'.format(month_num, month_str)

def project_dirs( proj_name ):
    proj_work_dir = Path( work_dir, month_label, proj_name )
    proj_data_dir = Path( data_dir, proj_name )
    proj_src_dir = Path( src_dir, proj_name )
    jupy_dir = Path( proj_work_dir, 'jupyter' )
    figs_dir = Path( proj_work_dir, 'FIGS' )
    os.makedirs( figs_dir, exist_ok=True )
    return proj_work_dir, proj_data_dir, proj_src_dir, jupy_dir, figs_dir

linestyle_dict = {'solid': 'solid',
 'dotted': (0, (1, 1)),
 'dashed': (0, (5, 5)),
 'dashdot': 'dashdot',
 'loosely dotted': (0, (1, 10)),
 'densely dotted': (0, (1, 1)),
 'loosely dashed': (0, (5, 10)),
 'densely dashed': (0, (5, 1)),
 'loosely dashdotted': (0, (3, 10, 1, 10)),
 'dashdotted': (0, (3, 5, 1, 5)),
 'densely dashdotted': (0, (3, 1, 1, 1)),
 'dashdotdotted': (0, (3, 5, 1, 5, 1, 5)),
 'loosely dashdotdotted': (0, (3, 10, 1, 10, 1, 10)),
 'densely dashdotdotted': (0, (3, 1, 1, 1, 1, 1))}

Ha2eV = const.value('hartree-electron volt relationship') #27.211
Ang2Bohr = 1.8897259886
Ang2m = 1e-10
F = const.value('Faraday constant')
elem_charge = const.value(u'elementary charge') #1.6021766208e-19
epd0 = const.value('electric constant') #8.854187817620389e-12
joule_to_ev = const.value('joule-electron volt relationship') #6.241509126e+18
ha_to_joule = const.value('hartree-joule relationship') #4.35974434e-18 J

def Coulomb_Energy( x, q1 = 1, q2 = -1, unit='Hartree'):
    coulomb_energy = q1*q2/x
    coulomb_energy *= 1/(4*math.pi*epd0) * 1/Ang2m * elem_charge**2  # in Joules
    if unit == 'Hartree':
        coulomb_energy /= ha_to_joule
    elif unit == 'eV':
        coulomb_energy *= joule_to_ev
    return coulomb_energy

def Coulomb_Force( x, q1 = 1, q2 = -1, unit='Hartree'):
    coulomb_force = q1*q2/x**2
    coulomb_force *= 1/(4*math.pi*epd0) * 1/Ang2m * elem_charge**2  # in Joules
    if unit == 'Hartree':
        coulomb_force /= ha_to_joule
    elif unit == 'eV':
        coulomb_force *= joule_to_ev
    return coulomb_force

def calculate_R2(xx,yy):
    yy_aver = np.mean(yy)
    return 1 - np.sum( (xx-yy)**2 )/np.sum( (xx-yy_aver)**2 )

#this is for arrays
def Energy_minimum_coordinates(dd,ee):
    return [ (d,e) for (d,e) in zip(dd,ee) if e==ee.min() ][0]
##this is for lists
#def Energy_minimum_coordinates(dd,ee):
#    return [ (d,e) for (d,e) in zip(dd,ee) if e==min(ee) ][0]

def chi_square( E1, E2 ):
    return np.sum( (E1-E2)**2 )

def calculate_C_vdW( r, E_nC ):
    return( np.sum(E_nC/r**6)/np.sum(1/r**12) )
    
def Dispersive_energy( x, C_vdW ):
    return C_vdW/x**6

def Dispersive_force( x, C_vdW ):
    return 6*C_vdW/x**7

def Dispersive_chi_square( x, E_nC ):
    C_vdW = calculate_C_vdW( x, E_nC )
    return np.sum( (E_nC - C_vdW/x**6)**2 )

def LJ_nm_energy( x, d0, e0, n, m ):
    return e0/(m-n)*( m*(d0/x)**n - n*(d0/x)**m )

def LJ_nm_force( x, d0, e0, n, m ):
    return m*n/(m-n)*e0/d0*((d0/x)**(n+1)-(d0/x)**(m+1)) 

def LJ_nm_chi_square( d, E, n, m, d0=False, e0=False ):
    if any([d0,e0]) == False:
        (d0,e0) = Energy_minimum_coordinates( d, E )
    E_LJ_nm = LJ_nm_energy( d, d0, e0, n, m )
    R2 = calculate_R2( E, E_LJ_nm )
    chi2 = chi_square( E, E_LJ_nm )
    #return np.sum( (E - E_LJ_nm)**2 ) 
    return R2, chi2

def LJ_AB_energy( x, n, m, A, B ):
    return A/x**n - B/x**m

def LJ_AB_force( x, n, m, A, B ):
    return n*A/x**(n+1)-m*B/x**(m+1)

def LJ_AB_chi_square( d, E, n, m, A, B ):
    lj_AB_energy = LJ_AB_energy( d, n, m, A, B )
    R2 = calculate_R2( E, lj_AB_energy )
    chi2 = chi_square( E, lj_AB_energy )
    return R2, chi2
    #return np.sum( (E - A/d**n + B/d**m )**2 )

def LJ_ab_energy( x, n, a, b ):
    return (a - b*np.log(x))/x**n

def LJ_ab_force( x, n, a, b ):
    #return n*a/x**(n+1) + b/x**(n+1)*(1-n*np.log(x))
    return (n*a + b - b*n*np.log(x) )/x**(n+1)

def LJ_ab_chi_square( d, E, n, a, b ):
    lj_ab_energy = LJ_ab_energy( d, n, a, b )
    R2 = calculate_R2( E, lj_ab_energy )
    chi2 = chi_square( E, lj_ab_energy )
    return R2, chi2
    #return np.sum( (E - a/d**n + b*np.log(d)/d**n )**2 )

def ab_from_Ed( e0, d0 ):
    a_LJ =  e0*(1-n*np.log(d0))*d0**n
    b_LJ = -e0*n*d0**n
    return a_LJ, b_LJ

def d_E_from_AB( n, m, A, B ):
    d = ( (n/m)*(A/B) )**(1/(n-m))
    E_A = (m-n)/m*A*d**(-n)
    E_B = (m-n)/n*B*d**(-m)
    if E_A == E_B:
        E = E_A
    else:
        print('Warning: EA={}, EB={}'.format(E_A,E_B) )
        E = E_A
    return d, E

def d_E_from_ab( n, a, b ):
    d = np.exp( 1/n + a/b)
    E = -b/n*np.exp(-(1+n*a/b))
    return d, E

def LJ_force_Ed(x, e0, d0, n):
    return -n**2*e0*np.log(d0/x)*d0**n/x**(n+1)

def calculate_AB(x,E,n,m):
    
    S1 = np.sum(E/x**n)
    S2 = np.sum(E/x**m)
    S12 = - np.sum(1/x**(n+m))
    S21 = - S12
    S11 =   np.sum(1/x**(2*n))
    S22 = - np.sum(1/x**(2*m))
    
    DELTA = S11*S22 - S12*S21
    
    A = ( S1*S22 - S2*S12 )/DELTA
    B = ( S2*S11 - S1*S21 )/DELTA
    
    return A,B

def calculate_ab(x,E,n):
    
    S1 = np.sum(E/x**n)
    S2 = np.sum(E*np.log(x)/x**n)
    S12 = - np.sum(np.log(x)/x**(2*n))
    S21 = - S12
    S11 = np.sum(1/x**(2*n))
    S22 = - np.sum(np.log(x)**2/x**(2*n))
    
    DELTA = S11*S22 - S12*S21
    
    a = ( S1*S22 - S2*S12 )/DELTA
    b = ( S2*S11 - S1*S21 )/DELTA
    
    return a,b

def angle_df( tmp_df, theta, phi, column_min = False ):
    tmp_df = tmp_df.loc[ tmp_df['Theta']==theta ].loc[tmp_df['Phi']==phi]
    tmp_df.sort_values(by=['Radius'], inplace=True)
    if column_min:
        line_min = tmp_df.loc[ tmp_df[column_min] == tmp_df[column_min].min() ]
    else:
        line_min = False
    return tmp_df, line_min

def create_subplot( tmp_df, tmp_theta, tmp_phi, tmp_column, tmp_marker, tmp_color,
                    label=None, ax=None, tmp_ls = '-', tmp_fs = 'full', tmp_size = 5):

    if ax is None:
        ax = plt.gca()

    sorted_df = angle_df( tmp_df, tmp_theta, tmp_phi )[0]
    tmp_rad = sorted_df['Radius'].values
    tmp_val = sorted_df[tmp_column].values

    ax.axhline(y=0.0, color = 'k', ls='--')
    sub_plt = ax.plot( tmp_rad, tmp_val, label=label, color=tmp_color, marker=tmp_marker,
                       linestyle=tmp_ls, fillstyle = tmp_fs, ms=tmp_size )

    return sub_plt
