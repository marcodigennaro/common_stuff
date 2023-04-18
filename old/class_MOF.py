# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 09:53:56 2021

@author: MDI0316
"""
from common_functions import project_dirs
from MOFS_QML import tml_runs_dir, ref_tml_runs_dir, qml_tml_runs_dir, proj_name, read_tmol_dft_ene, read_tmol_mp2_ene, coord_2_xyz
import os

proj_work_dir, proj_data_dir, proj_src_dir, jupy_dir, figs_dir = project_dirs( proj_name )

class MOF_conf:
    def __init__(self, qml_label, conf_label, reference=False ):
        self.qml_label = qml_label
        
        if reference == True:
            self.full_label = qml_label
            self.tbm_run_dir = os.path.join(ref_tml_runs_dir, self.full_label)
        else:
            self.conf_label = conf_label #self.H2_label = H2_label
            self.full_label = '{}_{}'.format(self.qml_label, self.conf_label )
            self.tbm_run_dir = os.path.join(qml_tml_runs_dir, self.full_label)
        
        self.gms_out_file = os.path.join(proj_data_dir, 'B3LYP', self.full_label) + '.out'
        self.gms_xyz_file = os.path.join(proj_data_dir, 'B3LYP', self.full_label) + '.xyz'
            
        self.tbm_coord_file = os.path.join( self.tbm_run_dir, 'coord')
        self.tbm_xyz_file = os.path.join( self.tbm_run_dir, '{}.xyz'.format(self.full_label))
        
        if os.path.exists( self.tbm_coord_file ) and not os.path.exists( self.tbm_xyz_file ):
            coord_2_xyz( self.tbm_coord_file, self.tbm_xyz_file )
        
        self.tbm_dft_ene_file = os.path.join( self.tbm_run_dir, 'job_0000', 'energy' )
        self.tbm_mp2_ene_file = os.path.join( self.tbm_run_dir, 'job_0001', 'energy' )
    
        
    def calc_dft_ene( self ):
        self.dft_ene = read_tmol_dft_ene( self.tbm_dft_ene_file )
        return self.dft_ene
    
    def calc_mp2_ene( self ):
        self.mp2_ene = read_tmol_mp2_ene( self.tbm_mp2_ene_file )
        return self.mp2_ene