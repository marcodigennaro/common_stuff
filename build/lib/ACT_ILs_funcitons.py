#!/home/mdi0316/anaconda3/bin/python

from common_stuff import *

act_data_dir  = Path( data_dir, 'ACT_ILs' )
mono_data_dir = Path( act_data_dir, 'MONOMERS' )
dime_data_dir = Path( act_data_dir, 'DIMERS' )

funct_dict = { 'PBE0'  : ['c','s'],
               'B3LYP' : ['r','o'],
               'TPSS'  : ['y','p'],
               'M06'   : ['m','h'],
               'M11'   : ['g','^'],
               'wB97x-D' : ['b','D'] }

basis_dict = { 'STO'  :      ['c','s'],
               'N311' :      ['r','o'],
               'APCseg-1'  : ['m','p'],
               'CCQ'   :     ['g','h'] }

full_funct_list = list(funct_dict.keys())
full_gbasis_list = [ 'STO', 'N21', 'N31', 'N311', 'G3L',
                     'DZV', 'TZV',
                     'CCD', 'CCT', 'CCQ',
                     'PCseg-0', 'PCseg-1', 'PCseg-2',
                     #'PCseg-3', 'PCseg-4',
                     'APCseg-0', 'APCseg-1',
                     #'APCseg-3', 'APCseg-4'
                     ]

reduced_gbasis_list = [ 'STO', 'N311', 'APCseg-1', 'CCQ']
reduced_funct_list = [ 'PBE0', 'B3LYP', 'M11' , 'wB97x-D' ]

def get_dimer_df( dime_label, basis, funct ):

    dft_ene_csv = Path( dime_data_dir, dime_label, '{}_{}_scan_dft_ene.csv'.format(basis, funct) )
    #dft_eda_csv = Path( dime_data_dir, dime_label, '{}_{}_scan_dft_eda.csv'.format(basis, funct) )
    mp2_ene_csv = Path( dime_data_dir, dime_label, '{}_scan_mp2_ene.csv'.format(basis) )
    #mp2_eda_csv = Path( dime_data_dir, dime_label, '{}_scan_mp2_eda.csv'.format(basis) )

    #opt_df = pd.read_csv( opt_csv, index_col = 0 )
    #err_df = pd.read_csv( err_csv, index_col = 0 )
    dft_ene_df = pd.read_csv( dft_ene_csv, index_col = 0 )
    #dft_eda_df = pd.read_csv( dft_eda_csv, index_col = 0 )
    mp2_ene_df = pd.read_csv( mp2_ene_csv, index_col = 0 )
    #mp2_eda_df = pd.read_csv( mp2_eda_csv, index_col = 0 )
    return( dft_ene_df, None, mp2_ene_df, None )

def get_mono_line( mono_label, basis=False, funct=False ):
    mono_csv = os.path.join( mono_data_dir, '{}_energy.csv'.format(mono_label) )
    mono_df = pd.read_csv( mono_csv, index_col = 0 )
    if basis==False and funct==False:
        return mono_df
    elif basis==False and funct!=False:
        return mono_df.loc[ mono_df['FUNCT'] == funct ]
    elif basis!=False and funct==False:
        return mono_df.loc[ mono_df['BASIS'] == basis ]
    else:
        return mono_df.loc[ (mono_df['BASIS'] == basis) & (mono_df['FUNCT'] == funct) ]

def get_dime_line( basis=False, funct=False ):
    dft_csv = os.path.join( dime_data_dir, 'CONVERGENCE', 'dft_results.csv' )
    mp2_csv = os.path.join( dime_data_dir, 'CONVERGENCE', 'mp2_results.csv' )
    dft_df = pd.read_csv( dft_csv, index_col = 0 )
    mp2_df = pd.read_csv( mp2_csv, index_col = 0 )

    dft_df.loc[dft_df['BASIS'].isin(full_gbasis_list)].loc[dft_df['FUNCT'].isin(full_funct_list)]
    mp2_df.loc[mp2_df['BASIS'].isin(full_gbasis_list)].loc[mp2_df['FUNCT'].isin(full_funct_list)]

    #null_dft_df = dft_df[dft_df.isnull().any(axis=1)]
    #null_mp2_df = dft_df[mp2_df.isnull().any(axis=1)]

    if basis==False and funct==False:
        return dft_df, mp2_df
    elif basis==False and funct!=False:
        dft_line = dft_df.loc[ dft_df['FUNCT'] == funct ]
        mp2_line = mp2_df.loc[ mp2_df['FUNCT'] == funct ]
        return dft_line, mp2_line
    elif basis!=False and funct==False:
        dft_line = dft_df.loc[ dft_df['BASIS'] == basis ]
        mp2_line = mp2_df.loc[ mp2_df['BASIS'] == basis ]
        return dft_line, mp2_line
    else:
        dft_line = dft_df.loc[ (dft_df['BASIS'] == basis) & (dft_df['FUNCT'] == funct) ]
        mp2_line = mp2_df.loc[ (mp2_df['BASIS'] == basis) & (mp2_df['FUNCT'] == funct) ]
        return dft_line, mp2_line

#mpl.rcParams["figure.figsize"] = 10, 20
mpl.rcParams["font.size"] = 20
mpl.rcParams["lines.markersize"] = 10 #12
mpl.rcParams["lines.linewidth"] = 3
mpl.rcParams["axes.linewidth"] = 3
plt.rcParams["xtick.major.size"] = 8
plt.rcParams["xtick.minor.size"] = 4
plt.rcParams["ytick.major.size"] = 8
plt.rcParams["ytick.minor.size"] = 4


linestyle_dict = {
     'solid':  'solid',
     'dotted': 'dotted',    # Same as (0, (1, 1)) or '.'
     'dashed': 'dashed',    # Same as '--'
     'dashdot': 'dashdot',  # Same as '-.'
     
     
     'loosely dotted' :       (0, (1, 10)),
     'dotted':                (0, (1, 1)),
     'densely dotted':        (0, (1, 1)),

     'loosely dashed':        (0, (5, 10)),
     'dashed':                (0, (5, 5)),
     'densely dashed':        (0, (5, 1)),

     'loosely dashdotted':   (0, (3, 10, 1, 10)),
     'dashdotted':            (0, (3, 5, 1, 5)),
     'densely dashdotted':    (0, (3, 1, 1, 1)),

     'dashdotdotted':         (0, (3, 5, 1, 5, 1, 5)),
     'loosely dashdotdotted': (0, (3, 10, 1, 10, 1, 10)),
     'densely dashdotdotted': (0, (3, 1, 1, 1, 1, 1))}