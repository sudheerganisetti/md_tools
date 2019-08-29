
"""
@Author : Sudheer Ganisetti
@Date   : Do 29. Aug 22:14:17 CEST 2019
This code is to extract the selected atoms from a lammps dump chkpt
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
      print("correct usage is: python3.6 this_program  ChkptFile   SelectedAtomsList")
      print("******************************************************************************************")
      exit()

  LAMMPS_CHKPT_FILE 	= sys.argv[1]
  SELECTED_ATOMS_FILE	= sys.argv[2]
  ganisetti_tools.get_ganisetti_tools_version()
  ref_atoms 		= ganisetti_tools.get_atoms_info_from_lammps(LAMMPS_CHKPT_FILE)
  selected_atoms1 	= open(SELECTED_ATOMS_FILE)
  selected_atoms	= []
  for line in selected_atoms1:
    elements=line.strip().split()
    for each_element in elements:
      selected_atoms.append(each_element)

  output1=open(LAMMPS_CHKPT_FILE+str("_selectedatoms.atoms"),'w')
  ganisetti_tools.write_imd_header_custom_property(output1,ref_atoms.box,"EMPTY")

  for i in selected_atoms:
    #print(ref_atoms.id[int(i)])
    ganisetti_tools.write_imd_atom_custom_property(output1,int(i),ref_atoms,0)
  output1.close()

