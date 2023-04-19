__author__ = "Marco Di Gennaro"
__email__ = "marco.di.gennaro@external.toyota-europe.com"

"""
parser for AVOGADRO files
"""

import math

def read_opt_file( opt_file ):
    try:
      with open( opt_file, 'r' ) as f:
         opt_log_lines = f.readlines()
      opt_log_lines.reverse()
      # find Time line
      for tmp_idx, tmp_line in enumerate(opt_log_lines[:10]):
        if 'cannot read input file!' in tmp_line:
           return 'cannot.read.inp.file', math.nan, math.nan, math.nan
        elif tmp_line.strip() == '==============================':
           pass
        elif tmp_line.strip() == '*** Open Babel Warning  in ReadMolecule':
           pass
        elif tmp_line.strip() == 'Problems reading an XYZ file: Cannot read the first line.':
           pass
        elif tmp_line.startswith('Time'):
           time_idx = tmp_idx
           break
      dtime = opt_log_lines[time_idx].split()[1].replace('seconds.','')
      if opt_log_lines[time_idx+1].strip() == 'STEEPEST DESCENT HAS CONVERGED':
         status = 'Converged'
         nstep, tot_ene, prev_ene = opt_log_lines[time_idx+2].split()
      else:
         status = 'Not.Converged'
         tot_ene, nstep = 2*[math.nan]
      return status, float(dtime), nstep, float(tot_ene)
    except(UnboundLocalError):
      print( 'UnboundLocalError in ', opt_file ) 
      return 'UnboundLocalError', math.nan, math.nan, math.nan



"""
d_ene_out_file( self, calc_type = 'fin' ):
 49          t0 = time.time()
 50          if calc_type == 'init':
 51             ene_xyz_file = self.init_xyz_file
 52             ene_out_file = self.init_ene_out_file
 53             ene_log_file = self.init_ene_log_file
 54          elif calc_type == 'fin':
 55             ene_xyz_file = self.fin_xyz_file
 56             ene_out_file = self.tot_ene_out_file
 57             ene_log_file = self.tot_ene_log_file
 58          else:
 59             print( unknown_calc_type )
 60          print_tab( 1, 'reading ene_file: {}'.format(ene_out_file) )
 61          cmd = 'grep "TOTAL ENERGY" {}'.format(ene_out_file)
 62          total_energy = float(sp.getoutput(cmd).split()[-2])
 63          t1 = time.time()
 64          print_tab( 1, '-------- ends ({} secs.)'.format(t1-t0) )
 65          return total_energy

48      def slurm_submit( self, running_labels = [], exclude_nodes = False ):
 49          if self.init_xyz_name in running_labels:
 50             print( 'Job running in the queue' )
 51          else:
 52             os.chdir( self.run_dir )
 53             shutil.copy( '/home/mdi0316/SCRIPTS/submit_avogadro.sh', './' )
 54             if exclude_nodes:
 55                cmd = 'sbatch --exclude={} -J {} submit_avogadro.sh {}'.format( exclude_nodes, self.init_xyz_name, self.init_xyz_name )
 56             else:
 57                cmd = 'sbatch -J {} submit_avogadro.sh {}'.format( self.init_xyz_name, self.init_xyz_name )
 58             sp.call( cmd , shell=True )
 59          return
"""
