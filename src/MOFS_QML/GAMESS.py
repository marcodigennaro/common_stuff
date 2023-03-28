__author__ = "Marco Di Gennaro"
__email__ = "marco.di.gennaro@external.toyota-europe.com"

"""
parser for GAMESS US files
"""

import os
import re
import mendeleev
import math
import Functions as FUN

## input information
run_type_msg        = 'INPUT CARD>  RUNTYP='
## general execution
tot_ener_msg        = 'TOTAL ENERGY ='
## Mulliken charges
mulliken_msg        = 'TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS'
distance_msg        = 'INTERNUCLEAR DISTANCES'
## scf cycle
mp2_no_scf_conv_msg = 'SCF DID NOT CONVERGE...NO MPLEVL=2 CALCULATION'
no_scf_conv_msg     = 'SCF IS UNCONVERGED'
density_conv_msg    = 'DENSITY CONVERGED'
energy_conv_msg     = 'ENERGY CONVERGED'
diis_conv_msg       = 'DIIS CONVERGED'
## geometry optimization
ok_geom_conv_msg    = 'EQUILIBRIUM GEOMETRY LOCATED'
no_geom_conv_msg    = 'THE GEOMETRY SEARCH IS NOT CONVERGED!'
fail_geom_msg       = 'FAILURE TO LOCATE STATIONARY POINT'
fail_geom_steps     = 'FAILURE TO LOCATE STATIONARY POINT, TOO MANY STEPS TAKEN'
fail_geom_scf       = 'FAILURE TO LOCATE STATIONARY POINT, SCF HAS NOT CONVERGED'
## Abortion messages:
opt_abort_msg       = '0OPTIMIZATION ABORTED.'
grad_abort_msg      = '-- GRADIENT OUT OF RANGE'
forceabort_msg      = '-- MAXIMUM ALLOWED FORCE (FMAXT) =   10.000'
internuclear_msg    = 'INTERNUCLEAR DISTANCES (ANGS.)'
atoms_apart_msg     = 'APART, QUITTING'
## geometry related msgs
beg_geom_msg        = 'BEGINNING GEOMETRY SEARCH POINT NSERCH='
cart_coords_msg     = "COORDINATES OF ALL ATOMS ARE \(ANGS\)"
inte_coords_msg     = "INTERNAL COORDINATES"
atomic_coords_msg   = "ATOM      ATOMIC                      COORDINATES \(BOHR\)"
end_geom_msg        = 'NSERCH:'
zmat_msg            = 'NO.   TYPE    I  J  K  L  M  N'
zmat_msg_bohr       = "NO.   TYPE    I  J  K  L  M  N        \(BOHR,RAD\)       \(ANG,DEG\)"   ## very begin and very end of the out file
                                                                                               ## plus begin of each NSERCH
zmat_msg_ang        = "NO.   TYPE    I  J  K  L  M  N         \(ANG,DEG\)     \(H\/B,H\/RAD\)" ## end of each NSERCH
cart_msg            = "#          XYZ \[AU\] "                                                 ## end of each NSERCH
zmat_label_dict     = { 'STRETCH' : 0, 'BEND' : 1, 'TORSION' : 2 }
inv_label_dict      = { 0 : 'STRETCH', 1 : 'BEND', 2 : 'TORSION' }
## EBS
occ_orbs_alpha_msg  =  'NUMBER OF OCCUPIED ORBITALS \(ALPHA\)'
occ_orbs_beta_msg   =  'NUMBER OF OCCUPIED ORBITALS \(BETA \)'
## force
cart_force_matrix_msg  =  'CARTESIAN FORCE CONSTANT MATRIX'
## execution messages
wall_time_msg       = 'TOTAL WALL CLOCK TIME'
cpus_time_msg       = 'CPU timing information for all processes'

## dict of known errors:
exec_dict = { 'EXECUTION OF GAMESS TERMINATED NORMALLY'        : 'TERMINATED.NORMALLY',
              'EXECUTION OF GAMESS TERMINATED -ABNORMALLY-'    : 'TERMINATED.ABNORMALLY',
              'Execution terminated due to error\(s\).'        : 'GMS.ERROR.FOUND' ,
              'Input file .inp does not exist.'                : 'MISSING.INPUT.FILE',
              'Fortran runtime error: End of file'             : 'FORT.EOF',
              'Fortran runtime error: No space left on device' : 'FORT.NO.SPACE.LEFT',
              'Fortran runtime error: Unformatted file structure has been corrupted' : 'FORT.CORRUPTED.FILE'
            }


unknown_error_dict = { 'check system limit for sysv semaphores.'       : 'semaphores.error.',
                       'THIS IS A SERIOUS FAILURE'                     : 'Serious.Failure',
                       'Error EINVAL'                                  : 'Error.Einval',
                       'ABORT FORCED.'                                 : 'Forced.Abort',
                       'semop errno=EINVAL'                            : 'Error.Einval',
                       'shmget returned an error'                      : 'Error.Shmget' ,
                       'dawrit'                                        : 'Error.Dawrit' }
#                       'terminated upon request'                       : 'DDI.Process.terminated' }

known_error_dict = { #'LOOSEN -QMTTOL' : 'reduce.qmttol',
                     'GRADIENT OUT OF RANGE'                         : 'Gradient.out.of.range',
                     'THERE ARE ATOMS LESS THAN'                     : 'atoms.too.close',
                     'MAXIMUM ANGULAR MOMENTUM EXCEEDED'             : 'max. ang. momentum exceeded',
                     'SCF DID NOT CONVERGE\.\.\.NO MPLEVL=2 CALCULATION' : 'SCF not converge (MP2)',
                     'SCF DID NOT CONVERGE\.\.\.NO CCTYP=CCSD\(T\)  CALCULATION' : 'SCF not converge (CCSDT)',
                     'ERROR! THERE ARE NOT 5 OR 6 TRANS/ROT MODES NUM T/R=' : 'TRANS/ROT MODES error',
                     'NO FORCE FIELD, SCF DOES NOT CONVERGE AT VIB0 POINT' : 'no force field',
                     'THE ORBITAL HESSIAN IS NOT POSITIVE DEFINITE.' : 'orbital hessian error',
                     'INSUFFICIENT DISTRIBUTED MEMORY REQUESTED'     : 'insufficient.distributed.memory',
                     'INSUFFICIENT REPLICATED MEMORY REQUESTED'      : 'insufficient.replicated.memory',
                     'MEMORY REQUEST EXCEEDS AVAILABLE MEMORY'       : 'memory.request.exceeds.available.memory',
                     'NOT ENOUGH REPLICATED MEMORY FOR PARALLEL '    : 'not.enough.replicated.memory' ,
                     'AMPLITUDE ITERATIONS DID NOT CONVERGE'         : 'amplitude did not converge',
                     'execution of EXETYP=CHECK is recommended'      : 'check.recommended',
                     'SPELLING OR LOGIC MISTAKE'                     : 'spelling.or.logig.mistake',
                     'FAILURE TO LOCATE STATIONARY POINT'            : 'Stationary.Point.Location.failed',
                     'NUMERICAL GRADIENT CANNOT PROCEED'             : 'Error.Numerical.Gradient',
                     'GLDIAG FAILURE IN -BKRNR-'                     : 'Error.BKRNR',
                     'THE REMAPPING OF H/I SHELLS FROM GAMESS ORDER' : 'Error.shell.remapping' }

error_dict = { **known_error_dict , **unknown_error_dict }

def write_inp_file(new_file, xyz_file, theory_level):

  with open(xyz_file, 'r') as f:
    xyz_lines = f.readlines()

  with open(new_file, 'w+') as f:
    f.write( ' $CONTRL SCFTYP=RHF UNITS=ANGSTROM MAXIT=200 \n' )
    f.write( '         EXETYP=RUN RUNTYP=OPTIMIZE' )
    if theory_level == 'B3LYP':
       f.write( ' DFTTYP=B3LYP')
    f.write( ' $END\n' )
    #f.write( ' $SYSTEM TIMLIM=600000 MWORDS=20 $END\n' )
    f.write( ' $STATPT NSTEP=1000   $END\n' )
    f.write( ' $SCF    NCONV=5      $END\n' )
    if theory_level == 'PM3':
       f.write( ' $BASIS  GBASIS=PM3   $END\n' )
    elif theory_level == 'B3LYP':
       f.write( ' $BASIS  GBASIS=N311  NDFUNC=1  NPFUNC=1  DIFFSP=.T.  DIFFS=.T.  NGAUSS=6 $END\n' )
       f.write( ' $DFT    IDCVER=4     $END\n' )
       f.write( ' $SCF    DIRSCF=.T.   DIIS=.T.  DAMP=.T.  $END\n' )
    f.write( ' $GUESS  GUESS=HUCKEL $END\n' )
    f.write( ' $DATA\n' )
    f.write( '  {}'.format(xyz_lines[1]) )
    f.write( '  C1\n' )
    for ll in xyz_lines[2:]:
      elem, xx, yy, zz = ll.split()
      cc = mendeleev.element(elem).atomic_number
      f.write( '  {}  {}  {}  {}  {}\n'.format(elem, cc, xx, yy, zz) )
    f.write( ' $END\n' )

def read_execution(out_file):
    execution = False
    wall_time = False
    if not os.path.exists( out_file ):
       execution = 'MISSING.OUTPUT.FILE'
       wall_time = math.nan
    else:
       out_lines = open( out_file, 'r', encoding = "ISO-8859-1" ).readlines()
       if len(out_lines) > 1e6:
          FUN.print_tab( 4, 'out len = {}'.format(len(out_lines)) )
          remove = input('Remove? (Y)')
          if remove.lower() in ['y', 'yes']:
             os.remove(out_file)
             json_file = input('Remove json file?')
             if json_file.lower() not in ['n', 'no']:
                os.remove(json_file)
             return 'MISSING.OUTPUT.FILE', math.nan
       out_lines.reverse()
       for o_count, o_line in enumerate(out_lines):
         if o_line.strip().startswith(wall_time_msg):
           read_idx = o_count #CPU time count
           wall_time = o_line.split()[4]
           break
       for o_line in out_lines: #[read_idx-10: read_idx+10]:
         for msg, status in exec_dict.items():
            exec_match = re.search( msg, o_line )
            if exec_match:
               execution = status
               break
    if not execution:
       FUN.print_tab( 1, ['could not find execution in', out_file ] )
       os.remove(out_file)
       return 'MISSING.OUTPUT.FILE', math.nan
    return execution, wall_time

def read_error(out_file):
    out_lines = open( out_file, 'r', encoding = "ISO-8859-1" ).readlines()
    out_lines.reverse()
    for o_line in out_lines:
      for k_err, v_err in known_error_dict.items():
        err_match = re.search( k_err, o_line )
        if err_match:
           FUN.print_tab( 4, 'Error found in GMS.read_error(): {}'.format(v_err) )
           return v_err
    for o_line in out_lines:
      for k_err, v_err in unknown_error_dict.items():
        err_match = re.search( k_err, o_line )
        if err_match:
           FUN.print_tab( 4, 'Error found in GMS.read_error(): {}'.format(v_err) )
           return v_err

def read_geometry(out_file):
    out_lines = open( out_file, 'r', encoding = "ISO-8859-1" ).readlines()
    xyz_file = out_file.replace('out', 'xyz')
    smi_file = out_file.replace('out', 'smi')
    for o_line in out_lines:
        if o_line.strip().startswith( 'TOTAL NUMBER OF ATOMS' ):
           natoms = int(o_line.split()[-1])
           break
    out_lines.reverse()
    for o_count, o_line in enumerate(out_lines):
      if ok_geom_conv_msg in o_line:
        status  = 'Converged'
        read_line = out_lines[o_count+3].replace('NSERCH:','')
        nstep   = read_line.split()[0]
        tot_ene = read_line.split()[2]
        fin_xyz = {}
        with open( xyz_file, 'w+' ) as f:
          f.write( '{}\n#gamess optimized structure\n'.format(natoms) )
          for at_idx, rev_idx in enumerate(range(4,natoms+4)):
            ee, cc, xx, yy, zz = out_lines[o_count-rev_idx].split()
            fin_xyz[at_idx] = { 'elem' : ee, 
                                'char' : float(cc), 
                                #'x' : float(xx), 'y' : float(yy), 'z' : float(zz) }
                                'coords' : [ float(xx) , float(yy), float(zz) ] }
            f.write( '{} {} {} {}\n'.format(ee, xx, yy, zz) )
        break
      elif any( [ msg in o_line for msg in [no_geom_conv_msg, fail_geom_msg] ] ):
        status  = 'Not.Converged'
        nstep   = math.nan
        tot_ene = math.nan
        fin_xyz = math.nan
        break
    return status, nstep, tot_ene, fin_xyz 

