# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 11:47:06 2020

@author: MDI0316
"""

from common_stuff import *

post_processing_dir = Path( codes_dir, 'abinit-9.2.2', 'scripts', 'post_processing' )


def read_edos( dos_file ):
    dos_df = pd.DataFrame()
    with open(dos_file,'r') as f:
        dos_lines = f.readlines()
    for d_line in dos_lines:
        if d_line.startswith('# Fermi energy :'):
            Fermi_En = float(d_line.split()[-1])
        elif not d_line.startswith('#'):
            energy, DOS, Int_DOS, DOS_tsmear_1, DOS_tsmear_2 = d_line.split()
            tmp_dict = {
                    'energy' : float(energy), 
                    'DOS' : float(DOS), 
                    'Integr. DOS' : float(Int_DOS), 
                    'DOS (tsmear/2)' : float(DOS_tsmear_1), 
                    'DOS (tsmear*2)' : float(DOS_tsmear_2)
                    }
            dos_df = dos_df.append( [tmp_dict], ignore_index=True )
    return dos_df, Fermi_En

def plot_edos( dos_file ):
    dos_df, Fermi_En = read_edos(dos_file)
    dos_df = dos_df.apply(pd.to_numeric)
    energy = dos_df['energy'].values
    DOS = dos_df['DOS'].values
    Int_DOS = dos_df['Integr. DOS'].values
    print(DOS)
    print(energy)
    fig, ax = plt.subplots()
    ax.plot(energy, DOS)
    ax.axvline(x=Fermi_En, ls='--')
    ax.set_ylabel('DOS')
    ax.set_xlabel('energy')
    plt.show()
    
proj_work_dir, proj_data_dir, jupy_dir, figs_dir = project_dirs( '12_DEC', 'CH4_eDOS' )
dos_file = os.path.join(proj_data_dir, 'xo_DOS')

plot_edos(dos_file)
