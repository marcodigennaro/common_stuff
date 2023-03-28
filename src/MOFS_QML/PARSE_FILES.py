__author__ = "Marco Di Gennaro"
__email__ = "marco.di.gennaro@external.toyota-europe.com"

"""
general function to parse files
"""

import math
import numpy as np
import re
import os
import shutil
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist
import ase
import json
import mendeleev
import subprocess as sp
import AVOGADRO as AVO
import GAMESS as GMS
import Functions as FUN

obabel='/home/mdi0316/CODES/openbabel-install/bin/obabel'

def read_xyz_file(xyz_file):
    with open( xyz_file, 'r' ) as f:
        xyz_lines = f.readlines()
    for count, line in enumerate(xyz_lines):
        if re.search('\*\^', line):
            xyz_lines[count] = line.replace('*^','E')
    nat=int(xyz_lines[0])
    species = [ line.split()[0] for line in xyz_lines[2:nat+2] ]
    xx = [ float(line.split()[1]) for line in xyz_lines[2:nat+2] ]
    yy = [ float(line.split()[2]) for line in xyz_lines[2:nat+2] ]
    zz = [ float(line.split()[3]) for line in xyz_lines[2:nat+2] ]
    return np.array(species), np.array([xx,yy,zz]).T

def dict_to_xyz(xyz_dict, xyz_file, comment=''):
    with open( xyz_file, 'w+' ) as f:
      f.write('{}\n'.format(len(xyz_dict)))
      f.write('{}\n'.format(comment))
      for k,v in xyz_dict.items():
        f.write('{}   {}   {}   {}\n'.format(v['elem'], v['coords'][0], v['coords'][1], v['coords'][2]))
    return

def xyz_to_dict(xyz_file):
    species, coords = read_xyz_file(xyz_file)
    xyz_dict = {}
    for idx,(s,c) in enumerate(zip(species, coords)):
      charge = mendeleev.element(s).atomic_number
      xyz_dict[idx] = { 'elem' : s, 
                        'char' : charge, 
#                        'x' : float(c[0]) , 'y' : float(c[1]), 'z' : float(c[2]) }
                        'coords' : [ float(c[0]) , float(c[1]), float(c[2]) ] }
    return xyz_dict

def xyz_to_smi(xyz_file):
    smi_file = xyz_file.replace('.xyz', '.smi')
    run_dir  = '/'.join(xyz_file.split('/')[:-1])
    os.chdir(run_dir)
    sp.call( '{} -i xyz {} -o smi > {} '.format(obabel, xyz_file, smi_file), shell=True )
    cmd = "cat {} ".format( smi_file )
    smi_out = sp.getoutput( cmd ).split()[0]
    return smi_out

def read_smi_file(smi_file):
    with open( smi_file, 'r' ) as f:
        smi_lines = f.readlines()
    return smi_lines[0].split()[0]

#def append_file_2_list(path_to_file, avo_list_file):
#    with open(avo_list_file, 'a') as f:
#      f.write('{}\n'.format(path_to_file))

def adjust_xyz( xyz_file, qml9_original = False ):
    xyz_lines = open(xyz_file,'r').readlines()
    if xyz_lines[0].startswith('Full Formula'):
       with open( xyz_file, 'w+' ) as f:
           f.write( '%d\n#pymatgen generated\n' %(len(xyz_lines)-4) )
           for line in xyz_lines[4:]:
               f.write( ' '.join(line.strip().split()[1:]) )
               f.write( '\n' )
    elif qml9_original:
       natoms = int(xyz_lines[0])
       with open( xyz_file, 'w+' ) as f:
           f.write( '{}\n'.format(natoms) )
           f.write( '#smi: {}\n'.format(xyz_lines[-2].split()[0]) )
           for line in xyz_lines[2:2+natoms]:
               f.write( line )
    return

def find_all_hs(xyz_file):
    species, coords = read_xyz_file(xyz_file)
    h_idxs = [ idx for (idx,s) in enumerate(species) if s == 'H']
    h_coordinates = coords[h_idxs]
    return h_idxs, h_coordinates

def decorated_dict( **kwargs ):
    label = ''
#    label_dict = {}
    for k,v in kwargs.items():
      label += k.replace('_idx', '').upper()
      label += '{:02d}_'.format(int(v))
#      label_dict[k] = int(v)
    return label[:-1] #, label_dict )

def sort_coordinates(index, species, coordinates):
    ''' index = index of an atom
        rearrange species and coordinates according to
        the distance from atom corresponding to index  '''
    dist_matrix = squareform(pdist(coordinates))  # matrix of distances
    dist_array = dist_matrix[index]
    new_species = [ s for (d,s) in sorted(zip(dist_array, species)) ]
    new_coords = [ c for (d,c) in sorted(zip(dist_array, coordinates)) ]
    return( new_species, new_coords )

def rotate_functional_group( ase_obj, o_idx = 1 ):
    if o_idx == 1:
        pass
    elif o_idx == 2:
        ase_obj.rotate( 90, 'x')
    elif o_idx == 3:
        ase_obj.rotate( 180, 'x')
    elif o_idx == 4:
        ase_obj.rotate( 270, 'x')
    else:
        print('unknown o_idx = {}'.format(o_idx))
        print(error)
    return( ase_obj )

def shift_functional_group( ase_obj, p_idx = 1, d_idx = 1 ):
    p0 = ase_obj.positions
    if p_idx == 1:
        p0 += [d_idx, 0, 0]
    elif p_idx == 2:
        p0 += [0, d_idx, 0]
    elif p_idx == 3:
        p0 += [0, -d_idx, 0]
    elif p_idx == 4:
        p0 += [0, 0, d_idx]
    elif p_idx == 5:
        p0 += [0, 0, -d_idx]
    else:
        print('unknown p_idx = {}'.format(p_idx))
        print(error)
    shifted_ase = ase.Atoms( ase_obj.symbols, p0 )
    return( shifted_ase )

def add_H2_molecule( species, coordinates, h2_position = 0 ):
    benzene_com = calculate_com( species[:11], coordinates[:11] )
    if h2_position == 0:
       h2_1 = benzene_com + [ 0., 0.,  2.0 ]
       h2_2 = benzene_com + [ 0., 0.,  2.7 ]
    elif h2_position == 1:
       h2_1 = benzene_com + [ 0., 0., -2.0 ]
       h2_2 = benzene_com + [ 0., 0., -2.7 ]
    elif h2_position == 2:
       h2_1 = benzene_com + [ 0.,  2.0, 0. ]
       h2_2 = benzene_com + [ 0.,  2.7, 0. ]
    elif h2_position == 3:
       h2_1 = benzene_com + [ 0., -2.0, 0. ]
       h2_2 = benzene_com + [ 0., -2.7, 0. ]
    elif h2_position == 4:
       h2_1 = benzene_com + [  2.0, 0., 0. ]
       h2_2 = benzene_com + [  2.7, 0., 0. ]
    elif h2_position == 5:
       h2_1 = benzene_com + [ -2.0, 0., 0. ]
       h2_2 = benzene_com + [ -2.7, 0., 0. ]
    species = np.hstack( [ species, [ 'H', 'H' ] ] )
    coordinates = np.vstack( [ coordinates, h2_1, h2_2 ] )
    return species, coordinates

def calculate_com(species, coordinates):
    com_x = 0
    com_y = 0
    com_z = 0
    tot_m = 0
    for s, [x,y,z] in zip(species, coordinates):
        mass = mendeleev.element(s).mass
        com_x += mass*x
        com_y += mass*y
        com_z += mass*z
        tot_m += mass
    com_x /= tot_m
    com_y /= tot_m
    com_z /= tot_m
    return np.array( [ com_x, com_y, com_z ] )

def xyz_to_file(species, coordinates, xyz_file, comment_line=''):
    with open( xyz_file, 'w+' ) as f:
      f.write( '{}\n'.format( len(species) ) )
      f.write( '{}\n'.format( comment_line ) )
      for s, [x, y, z] in zip(species, coordinates):
          f.write( '{} {} {} {}\n'.format(s, x, y, z) )
    return

complete_dict = { 'PM3' : {'status' : 'TERMINATED.NORMALLY', 'geometry' : 'Converged' }}
def read_key_from_dict( key, read_dict, theory_level ):
    if key in read_dict.keys():
       conf_dict = read_dict[key][theory_level]
       if 'status' in conf_dict.keys():
          key_status = conf_dict['status']
          if key_status == 'TERMINATED.NORMALLY':
             FUN.print_tab(3, '{} complete'.format(key))
             return 0
          elif key_status in [ 'TERMINATED.ABNORMALLY', 'GMS.ERROR.FOUND', 'Known.Error', 'Unknown.Error']:
             FUN.print_tab(3, '{} {}'.format(key, key_status))
             return 0 
          elif key_status in [ 'missing_out' ]:
             return 1 
          else:
             print( key_status )
             print( unknown, key_status )
       else:
          FUN.print_tab(3, '{} unknown status {}'.format(key, conf_dict.keys()))
    else:
        FUN.print_tab(  3, '{} missing'.format(key) )
    return 1

def read_qml_files( qml_obj, qml_level ):
    FUN.print_tab( 3, 'read_qml_files in: {}'.format(qml_obj.run_dir) )
    out_dict = {} 
 
    if qml_obj.theory_level == 'PM3' and not os.path.exists( qml_obj.init_xyz_file ):
       out_dict['status'] = 'missing_inp'
       if qml_level == 'isolated_fungr':
            qml_obj.write_isolated_fungr()
       elif qml_level == 'default_decorated_benzene':
            qml_obj.default_decorate()
       elif qml_level == 'modified_decorated_benzene':
            qml_obj.modify_relative_position()
       elif qml_level == 'H2_modified_decorated_benzene':
            qml_obj.insert_H2()
       elif qml_level in ['B3LYP']:
            pass
       else:
            print(qml_obj) 
            print(error, unknownfdjadjsm)
 
    if os.path.exists( qml_obj.out_file ):
       ## UFF/AVOGADRO
       if qml_obj.theory_level == 'UFF':
          out_dict['status'],  out_dict['runtime'], \
          out_dict['nstep'],   out_dict['tot_ene'] = AVO.read_out_file( qml_obj.out_file )
          if out_dict['status'] == 'cannot.read.inp.file':
             FUN.print_tab( 3, 'Reomoving folder: {}'.format(out_dict['status']) )
             shutil.rmtree( qml_obj.run_dir )
             return out_dict
          else:
             print( check_possibilites_UFF )
          out_dict['fin_xyz'] = xyz_to_dict( qml_obj.fin_xyz_file )

       ## GAMESS
       elif qml_obj.theory_level in ['PM3', 'B3LYP']:
          out_dict['run_dir'] = qml_obj.run_dir
          out_dict['status'], out_dict['runtime'] = GMS.read_execution( qml_obj.out_file )
          if out_dict['status'] in [ 'TERMINATED.ABNORMALLY', 'GMS.ERROR.FOUND' ]:
            out_dict['error'] = GMS.read_error( qml_obj.out_file )
            print( 'Warning: {} error found'.format(out_dict['error']) )
            if not os.path.exists(qml_obj.out_file):
               out_dict['status'] = 'missing_out'
            ## LIST of known errors:
            elif out_dict['error'] in GMS.known_error_dict.values():
               out_dict['status'] = 'Known.Error'
            elif out_dict['error'] in GMS.unknown_error_dict.values():
               print( 'Warning: removing folder due to ', out_dict['error'])
               out_dict['status'] = 'Unknown.Error'
            elif out_dict['error'] == None:
               out_dict['status'] = 'missing_out'
            else:
               print( qml_obj ) 
               print( out_dict['run_dir'] )
               print( out_dict['error'] )
               print( unknonw )
            return out_dict
          elif out_dict['status'] == 'MISSING.OUTPUT.FILE': 
            print( 'Warning: \'MISSING.OUTPUT.FILE\'' )
            out_dict['status'] = 'missing_out'
          else:
            out_dict['geometry'], out_dict['nsteps'], out_dict['tot_ene'], fin_xyz_dict = \
            GMS.read_geometry( qml_obj.out_file )
            out_dict['fin_xyz'] = fin_xyz_dict
            if out_dict['geometry'] == 'Not.Converged':
               return out_dict
       if (os.path.exists(qml_obj.smi_file) and os.stat(qml_obj.smi_file).st_size == 0) or not os.path.exists(qml_obj.smi_file):
          xyz_to_smi(qml_obj.xyz_file)
       out_dict['fin_smi'] = read_smi_file( qml_obj.smi_file )
    else:
       out_dict['status'] = 'missing_out'
    return out_dict
