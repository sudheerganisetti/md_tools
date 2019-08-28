
"""
@Author : Sudheer Ganisetti
@Date   : So 18. Aug 22:06:23 CEST 2019
This code is to compute MSD of given list of atoms
"""

import sys
import ganisetti_tools
import numpy as np

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if __name__=="__main__":

  tot_argv = len(sys.argv)
  if tot_argv != 4:
      print("************************************** S. Ganisetti **************************************")
      print("Error: usage is wrong")
      print("correct usage is: python3.6 this_program  Ref_ChkptFile   Cur_ChkptFile   AtomsList.data")
      print("******************************************************************************************")
      exit()

  REF_LAMMPS_FILE = sys.argv[1]
  CUR_LAMMPS_FILE = sys.argv[2]
  ATOMS_DATA_FILE = sys.argv[3]
  ganisetti_tools.get_ganisetti_tools_version()
  ref_atoms = ganisetti_tools.get_atoms_info_from_lammps(REF_LAMMPS_FILE)
  cur_atoms = ganisetti_tools.get_atoms_info_from_lammps(CUR_LAMMPS_FILE)
  atoms_list=np.loadtxt(ATOMS_DATA_FILE)

  #boxx2=pow(ref_atoms.box_xx,2)
  #boxy2=pow(ref_atoms.box_yy,2)
  #boxz2=pow(ref_atoms.box_zz,2)
  box_x=ref_atoms.box_xx
  box_y=ref_atoms.box_yy
  box_z=ref_atoms.box_zz
  MSD=0.0
  MSD_x=0.0
  MSD_y=0.0
  MSD_z=0.0

  for i in atoms_list:
    dx = ref_atoms.posx[i] - cur_atoms.posx[i]
    dy = ref_atoms.posy[i] - cur_atoms.posy[i]
    dz = ref_atoms.posz[i] - cur_atoms.posz[i]
    if abs(dx) > box_x/2:
      dx = abs(dx) - box_x
    if abs(dy) > box_y/2:
      dy = abs(dy) - box_y
    if abs(dz) > box_z/2:
      dz = abs(dz) - box_z
    dx2 = dx*dx
    dy2 = dy*dy
    dz2 = dz*dz
    dr2 = dx2+dy2+dz2
    MSD_x = MSD_x + dx2
    MSD_y = MSD_y + dy2
    MSD_z = MSD_z + dz2
    MSD   = MSD   + dr2

  print("%lf\t%lf\t%lf\t%lf\n" %(MSD_x,MSD_y,MSD_z,MSD/len(atoms_list)))



