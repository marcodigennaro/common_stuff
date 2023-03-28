#!/home/mdi0316/anaconda3/bin/python

from common_stuff import *

work_label = 'MOFS_RGI2972'
month_label = '03_MAR'
mofs_data_dir  = Path( data_dir, 'MOFS_RGI2972' )
mofs_work_dir  = Path( work_dir, month_label, 'MOFS_RGI2972' )

ML_dir = os.path.join( mofs_work_dir, 'ML' )
Regression_dir = os.path.join( mofs_work_dir, 'ML', 'Regression' )

#import json
#import ast
#import joblib

from sklearn import preprocessing
import ast, math, random
from numpy import std
from sklearn.metrics import r2_score

#from pathlib import Path

#import matplotlib.pyplot as plt
import mendeleev
#import shutil
from itertools import cycle
#import numpy as np
#import pandas as pd
from scipy import linalg as LA

import ase
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones
from numpy import std

from sklearn.preprocessing import StandardScaler

from sklearn.gaussian_process.kernels import WhiteKernel, ExpSineSquared
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_absolute_error as mae
from sklearn.metrics import pairwise
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score as r2s
from sklearn.model_selection import RandomizedSearchCV
from sklearn.svm import SVR
from sklearn.metrics import  pairwise

#import qml
#from qml.fchl import get_local_symmetric_kernels
#from qml.math import cho_solve
#from qml.kernels import gaussian_kernel
#from qml.kernels import get_local_kernels_gaussian
#
#from dscribe.descriptors import ACSF
#from dscribe.descriptors import SOAP

from sklearn.ensemble import RandomForestRegressor

#def xyz_to_ase(xyz_file):
#
#    with open( xyz_file, 'r' ) as f:
#        xyz_lines = f.readlines()
#
#    for count, line in enumerate(xyz_lines):
#        if re.search('\*\^', line):
#            xyz_lines[count] = line.replace('*^','E')
#            print(count,xyz_lines[count])
#
#    nat=int(xyz_lines[0])
#
#    species = [ line.split()[0] for line in xyz_lines[2:nat+2] ]
#    atomic_numbers = [ mendeleev.element(line.split()[0]).atomic_number for line in xyz_lines[2:nat+2] ]
#    xx = [ float(line.split()[1]) for line in xyz_lines[2:nat+2] ]
#    yy = [ float(line.split()[2]) for line in xyz_lines[2:nat+2] ]
#    zz = [ float(line.split()[3]) for line in xyz_lines[2:nat+2] ]
#    positions = [ [x,y,z] for (x,y,z) in zip(xx, yy, zz) ]
#    molecule = Atoms( species, positions = positions )
#
#    return molecule



def print_distribution( data, destination ):
    print('Plotting distribution:', destination )
    if not os.path.exists( destination ):
        fig, ax = plt.subplots()
        ax.hist(data, bins=100, density=True)
        ax.set_ylabel('distribution')
        ax.set_xlabel('Interaction Energy (Ha)')
        plt.savefig( destination )

def standardize_df( ene_col = 'MP2.INT.EN.' ):
    dime_results_csv = os.path.join( mofs_dir, 'MOFS_CSV/dime_results.csv' )
    dime_results_df = pd.read_csv( dime_results_csv , index_col = 0 )

    os.makedirs('Regression/FIGS', exist_ok=True )

    start_distribution_pdf = os.path.join( ML_dir, 'Regression/FIGS/start_distribution.pdf' )
    print_distribution( dime_results_df[ene_col], start_distribution_pdf )
    dime_results_df.drop( dime_results_df[ dime_results_df['FG'] == '01_H' ].index, inplace=True )
    dime_results_df = dime_results_df.loc[ ( dime_results_df[ene_col] > -0.003)  &
                                                        ( dime_results_df[ene_col] < 0.) ]

    shifted_distribution_pdf = os.path.join( ML_dir, 'ML/FIGS/shifted_distribution.pdf' )
    shifted_df = pd.DataFrame(dime_results_df)
    min_ene = shifted_df[ene_col].min()
    shifted_df[ene_col] = shifted_df[ene_col] - min_ene
    print_distribution( shifted_df[ene_col], shifted_distribution_pdf )

    scaler = StandardScaler()
    scaled_df = pd.DataFrame(dime_results_df)
    scaled_df[ene_col] = scaler.fit_transform(dime_results_df[[ene_col]])
    scaled_distribution_pdf = os.path.join( ML_dir, 'ML/FIGS/scaled_distribution.pdf' )
    print_distribution( scaled_df[ene_col], scaled_distribution_pdf )

    return shifted_df

def center_of_mass( masses, coordinates ):

    comx, comy, comz = 0, 0, 0
    totm = sum(masses)
    for m,[x,y,z] in zip(masses, coordinates):
        comx += m*x
        comy += m*y
        comz += m*z
    comx /= totm
    comy /= totm
    comz /= totm
    com = np.array([comx, comy, comz])
    return com

def coord_to_xyz( labels, coords, out_file ):
    with open( out_file, 'w+' ) as f:
        f.write(str(len(coords)))
        f.write('\n\n')
        for label, coord in zip(labels, coords):
            f.write('{} {} {} {}\n'.format(label, *coord))
    return

def xyz_to_df(xyz_file):
    with open( xyz_file, 'r' ) as f:
        xyz_lines = f.readlines()
    for count, line in enumerate(xyz_lines):
        if re.search('\*\^', line):
            xyz_lines[count] = line.replace('*^','E')
            print(count,xyz_lines[count])
    nat=int(xyz_lines[0])
    species = [ line.split()[0] for line in xyz_lines[2:nat+2] ]
    atomic_numbers = [ mendeleev.element(line.split()[0]).atomic_number for line in xyz_lines[2:nat+2] ]
    xx = [ float(line.split()[1]) for line in xyz_lines[2:nat+2] ]
    yy = [ float(line.split()[2]) for line in xyz_lines[2:nat+2] ]
    zz = [ float(line.split()[3]) for line in xyz_lines[2:nat+2] ]
    molecule_df = pd.DataFrame()
    for specie, at_num, x, y, z in zip( species, atomic_numbers, xx, yy, zz ):
        molecule_df = molecule_df.append( [ { 'specie' : specie, 'at.number' : at_num,
                                              'x' : x, 'y' : y, 'z' : z } ], ignore_index = True )
    return molecule_df

def read_xyz_file(xyz_file):
    with open( xyz_file, 'r' ) as f:
        xyz_lines = f.readlines()
    for count, line in enumerate(xyz_lines):
        if re.search('\*\^', line):
            xyz_lines[count] = line.replace('*^','E')
            print(count,xyz_lines[count])
    nat=int(xyz_lines[0])
    species = [ line.split()[0] for line in xyz_lines[2:nat+2] ]
    atomic_numbers = [ mendeleev.element(line.split()[0]).atomic_number for line in xyz_lines[2:nat+2] ]
    xx = [ float(line.split()[1]) for line in xyz_lines[2:nat+2] ]
    yy = [ float(line.split()[2]) for line in xyz_lines[2:nat+2] ]
    zz = [ float(line.split()[3]) for line in xyz_lines[2:nat+2] ]
    return np.array(species), np.array(atomic_numbers), np.array([xx,yy,zz]).T

def fchl_cv_kkr( representation, energy ):

    alphas = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12]
    sigmas = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]

    print('representation.shape =', representation.shape)
    X_train, X_test, y_train, y_test = train_test_split(representation, energy, test_size=0.33, random_state=42)

    print('X_train.shape =', X_train.shape)
    print('X_test.shape  =', X_test.shape)
    print('y_train.shape =', y_train.shape)
    print('y_test.shape  =', y_test.shape)

    K_FCHL = get_local_symmetric_kernels(X_train, sigmas, cut_distance=10.0)

    fchl_results = pd.DataFrame()
    best_mae = 1e4
    for alpha in alphas:
      for sigma, K in zip( sigmas, K_FCHL ):
          K[np.diag_indices_from(K)] += alpha
          coeffs = cho_solve(K, y_train)

          #K_train = gaussian_kernel(X_train, X_train, sigma)
          pred_train = np.dot(K, coeffs)
          mae_train = mae(y_train, pred_train)
          std_train = std(y_train - pred_train)
          r2s_train = r2s(y_train, pred_train)

          K_test = gaussian_kernel(X_test, X_train, sigma)
          pred_test = np.dot(K_test, coeffs)
          mae_test = mae(y_test, pred_test)
          std_test = std(y_test - pred_test)
          r2s_test = r2s(y_test, pred_test)

          fchl_results = fchl_results.append( [ { 'alpha' : alpha, 'sigma' : sigma,
                                                  'mae.train' : mae_train, 'std.train' : std_train, 'r2s.train' : r2s_train,
                                                  'mae.test'  : mae_test,  'std.test'  : std_test,  'r2s.test'  : r2s_test  } ],
                                                  ignore_index = True )

          if mae_test < best_mae:
             best_mae = mae_test
             best_sigma = sigma
             best_alpha = alpha
             best_coeffs = coeffs

             best_pred_train = pred_train
             best_pred_test = pred_test

    train = [ X_train, y_train, best_pred_train ]
    test  = [ X_test,  y_test,  best_pred_test ]
    return fchl_results, best_coeffs, train, test

def cv_krr( representation, energy ):

    print('representation.shape =', representation.shape)
    X_train, X_test, y_train, y_test = train_test_split(representation, energy, test_size=0.2, random_state=42)

    print('X_train.shape =', X_train.shape)
    print('X_test.shape  =', X_test.shape)
    print('y_train.shape =', y_train.shape)
    print('y_test.shape  =', y_test.shape)

    CV_kernel = GridSearchCV(   KernelRidge(), scoring='neg_mean_absolute_error',
                                            param_grid=KRR_param_grid, cv=3 )
    CV_kernel.fit(X_train, y_train)

    best_alpha = CV_kernel.best_params_['alpha']
    best_gamma = CV_kernel.best_params_['gamma']
    best_kernel = CV_kernel.best_params_['kernel']

    print(KRR_param_grid)

    print( 'best_alpha  =', best_alpha )
    print( 'best_gamma  =', best_gamma )
    print( 'best_kernel =', best_kernel )

    regr = KernelRidge(alpha=best_alpha, gamma=best_gamma, kernel=best_kernel)
    regr.fit(X_train, y_train)

    pred_train = regr.predict(X_train)
    pred_test = regr.predict(X_test)

    plt.scatter(pred_train, y_train)
    plt.scatter(pred_test, y_test)
    plt.show()

    return { 'alpha' : best_alpha, 'gamma' : best_gamma, 'kernel' : best_kernel }

def sklearn_svr( X_train, X_test, y_train, y_test, parameters_file = None ):

    #default parameters
    kernel  = 'rbf'
    degree  = 3
    gamma   = 'auto'
    C       = 1
    epsilon = 0.1
         
    if os.path.exists( parameters_file ):
       print( 'Read parameters from: ', parameters_file )
       params_df = pd.read_csv( parameters_file, index_col = 0 )
       if 'kernel' in params_df.columns:
           kernel = params_df['kernel'].values[0]
       if 'degree' in params_df.columns:
           degree = params_df['degree'].values[0]
       if 'gamma' in params_df.columns:
           gamma = params_df['gamma'].values[0]
       if 'C' in params_df.columns:
           C = params_df['C'].values[0]
       if 'epsilon' in params_df.columns:
           epsilon = params_df['epsilon'].values[0]

    model = SVR( kernel = kernel, degree = degree, gamma = gamma, C = C, epsilon = epsilon )
    print( model )
    model.fit(X_train, y_train)
    pred_train = model.predict(X_train)
    pred_test  = model.predict(X_test)
    train = [ X_train, y_train, pred_train ]
    test  = [ X_test,  y_test,  pred_test  ]

    return train, test, model 

def sklearn_rfr( X_train, X_test, y_train, y_test, parameters_file = None ):

    #default parameters
    n_estimators = 100
    max_depth = None
    min_samples_split = 2
    min_samples_leaf = 1
    max_features = 'auto'
    bootstrap = True 

    if os.path.exists( parameters_file ):
       print( 'Read parameters from: ', parameters_file )
       params_df = pd.read_csv( parameters_file, index_col = 0 )
       if 'n_estimators' in params_df.columns:
           n_estimators = params_df['n_estimators'].values[0]
       if 'max_depth' in params_df.columns:
           max_depth  = params_df['max_depth'].values[0]
       if 'max_features' in params_df.columns:
           max_features = params_df['max_features'].values[0]
       if 'min_samples_leaf' in params_df.columns:
           min_samples_leaf = params_df['min_samples_leaf'].values[0]
       if 'min_samples_split' in params_df.columns:
           min_samples_split = params_df['min_samples_split'].values[0]
       if 'bootstrap' in params_df.columns:
           bootstrap  = params_df['bootstrap'].values[0]

    if math.isnan( max_depth ):
       max_depth = None

    model = RandomForestRegressor( criterion = 'mae',
                                   bootstrap = bootstrap,
                                   max_depth = max_depth,
                                   max_features = max_features,
                                   min_samples_leaf = min_samples_leaf,
                                   min_samples_split = min_samples_split,
                                   n_estimators = n_estimators )

    print( model )
    model.fit(X_train, y_train)

    pred_train = model.predict(X_train)
    pred_test  = model.predict(X_test)

    train = [ X_train, y_train, pred_train ]
    test  = [ X_test,  y_test,  pred_test  ]

    return train, test, model 

def sklearn_krr( X_train, X_test, y_train, y_test, parameters_file = None ):

    #default parameters
    alpha=1e-9 
    gamma=0.2
    kernel='rbf'

    if os.path.exists( parameters_file ):
       print( 'Read parameters from: ', parameters_file )
       params_df = pd.read_csv( parameters_file, index_col = 0 )
       if 'alpha' in params_df.columns:
           alpha = params_df['alpha'].values[0]
       if 'gamma' in params_df.columns:
           gamma = params_df['gamma'].values[0]
       if 'kernel' in params_df.columns:
           kernel = params_df['kernel'].values[0]

    model = KernelRidge(alpha=alpha, gamma=gamma, kernel=kernel)
    print( model )

    model.fit(X_train, y_train)
    pred_train = model.predict(X_train)
    pred_test = model.predict(X_test)

    train = [ X_train, y_train, pred_train ]
    test  = [ X_test,  y_test,  pred_test  ]

    return train, test, model 

def my_kkr( representation, energy, alpha=1e-9, gamma=0.2, kernel='rbf' ):

    X_train, X_test, y_train, y_test = train_test_split(representation, energy, test_size=0.1, random_state=42)

    print(X_train.shape)
    K_train = pairwise.pairwise_kernels(X_train, X_train)
    K_test = pairwise.pairwise_kernels(X_test, X_train)

    K_train[np.diag_indices_from(K_train)] += alpha
    coeffs = cho_solve(K_train, y_train)

    pred_train = np.dot(K_train, coeffs)
    pred_test = np.dot(K_test, coeffs)

    plt.scatter( y_train, pred_train )
    plt.scatter( y_test, pred_test )
    plt.show()
    exit()
    train = [ X_train, y_train, pred_train ]
    test  = [ X_test,  y_test,  pred_test  ]

    return train, test, clf.dual_coef_

SVR_param_grid = { 'kernel'  : ['linear', 'rbf'],
                   'C'       : np.logspace(-1, 2, 25),
                   'epsilon' : np.logspace(-5,-4, 25),
                   'gamma'   : np.logspace(-5,-2, 25),
                   'degree' : [2,3,4,5] }

KRR_param_grid = { 'alpha'  : np.logspace(-12, -6, 7),
                   "gamma"  : np.logspace( -5, -2, 50),
                   "kernel" : ['rbf', 'laplacian']}

RFR_param_grid = { 'n_estimators': [ int(x) for x in np.linspace( 800, 1200, 18 ) ],
                   'min_samples_split': list(range(1,11)),
                   'min_samples_leaf':  list(range(1,11)),
                   'max_depth': [ math.floor(x) for x in np.logspace( 0, 3 , 25) ],
                   'max_features': ['auto', 'sqrt'],
                   'bootstrap': [True, False] }

RFR_param_grid['max_depth'].append(None)
