"""
@Author : Sudheer Ganisetti
@Date   : Mo 27. Jul 16:36:19 CEST 2020
This code is to extract the given list of atoms from a chkpt
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
      print("correct usage is: python3.6 this_program  ChkptFile_Prefix ChkptFile_PostFix AtomsList.data")
      print("******************************************************************************************")
      exit()
  CUR_LAMMPS_FILE_PREFIX    = sys.argv[1]
  CUR_LAMMPS_FILE_POSTFIX   = sys.argv[2]
  ATOMS_DATA_FILE           = sys.argv[3]
  ganisetti_tools.get_ganisetti_tools_version()

  atoms_list1=np.loadtxt(ATOMS_DATA_FILE,dtype="int")
  atoms_list=[]
  collected_data={}
  try:
        for i in atoms_list1:
          atoms_list.append(i)
  except TypeError:
        atoms_list1=open(ATOMS_DATA_FILE,"r")
        atoms_list2=atoms_list1.readline().strip().split(" ")[0]
        atoms_list.append(int(atoms_list2))

  CUR_LAMMPS_FILE             = str(CUR_LAMMPS_FILE_PREFIX)+str(".")+str(CUR_LAMMPS_FILE_POSTFIX)
  config1                   = ganisetti_tools.get_atoms_info_from_lammps(CUR_LAMMPS_FILE)
  output1=open(CUR_LAMMPS_FILE_PREFIX+str("_SelectedAtoms.atoms"),'w')
  ganisetti_tools.write_imd_header_basic(output1,config1.box)
  for atom_id in atoms_list:
  	output1.write("%d  %d  %lf  %lf  %lf\n" %(atom_id,config1.type[atom_id],config1.posx[atom_id],config1.posy[atom_id],config1.posz[atom_id]))
  output1.close()





