
"""
@Author : Sudheer Ganisetti
@Date   : Mi 28. Aug 22:54:43 CEST 2019
This code is to compute MSD of each atom and write this information into imd format file
"""

import sys
import ganisetti_tools
import numpy as np

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if __name__=="__main__":

  tot_argv = len(sys.argv)
  if tot_argv != 3:
      print("************************************** S. Ganisetti **************************************")
      print("Error: usage is wrong")
      print("correct usage is: python3.6 this_program  Ref_ChkptFile   Cur_ChkptFile")
      print("******************************************************************************************")
      exit()

  REF_LAMMPS_FILE = sys.argv[1]
  CUR_LAMMPS_FILE = sys.argv[2]
  ganisetti_tools.get_ganisetti_tools_version()
  ref_atoms = ganisetti_tools.get_atoms_info_from_lammps(REF_LAMMPS_FILE)
  cur_atoms = ganisetti_tools.get_atoms_info_from_lammps(CUR_LAMMPS_FILE)

  box_x=ref_atoms.box_xx
  box_y=ref_atoms.box_yy
  box_z=ref_atoms.box_zz
  MSD=0.0
  MSD_x=0.0
  MSD_y=0.0
  MSD_z=0.0
  output1=open(CUR_LAMMPS_FILE+str("_msd.atoms"),'w')
  ganisetti_tools.write_imd_header_custom_property(output1,cur_atoms.box,"msd")

  for i in ref_atoms.id:
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
    ganisetti_tools.write_imd_atom_custom_property(output1,i,cur_atoms,dr2)
  output1.close()


