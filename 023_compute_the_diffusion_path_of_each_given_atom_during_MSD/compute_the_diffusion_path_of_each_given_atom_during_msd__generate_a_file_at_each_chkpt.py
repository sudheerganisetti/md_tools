
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
  if tot_argv != 6:
      print("************************************** S. Ganisetti **************************************")
      print("Error: usage is wrong")
      print("correct usage is: python3.6 this_program  ChkptFile_Prefix  Start_Chkpt_Number  End_Chkpt_Number Interval_Chkpt_Number  AtomsList.data")
      print("******************************************************************************************")
      exit()

  #REF_LAMMPS_FILE           = sys.argv[1]
  CUR_LAMMPS_FILE_PREFIX    = sys.argv[1]
  START_CHKPT_NUMBER        = int(sys.argv[2])
  END_CHKPT_NUMBER          = int(sys.argv[3])
  INTERVAL_CHKPT_NUMBER     = int(sys.argv[4])
  ATOMS_DATA_FILE           = sys.argv[5]
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

  atom_count1=1
  for j in range(START_CHKPT_NUMBER,END_CHKPT_NUMBER,INTERVAL_CHKPT_NUMBER):
    k=("%08d" %(j))
    CUR_LAMMPS_FILE 	= str(CUR_LAMMPS_FILE_PREFIX)+str(".")+str(k)+str(".dump")
    cur_atoms 			= ganisetti_tools.get_atoms_info_from_lammps(CUR_LAMMPS_FILE)
    output1=open(CUR_LAMMPS_FILE_PREFIX+str("_selected_atoms_from_all_files.")+str(k)+str(".atoms"),'w')
    ganisetti_tools.write_imd_header_custom_property(output1,cur_atoms.box,"molecule_id")
    if j != START_CHKPT_NUMBER:
      for i in collected_data.keys():
        output1.write("%d  %d  %lf  %lf  %lf  %d \n" %(collected_data[i]["id"],collected_data[i]["type"],collected_data[i]["posx"],collected_data[i]["posy"],collected_data[i]["posz"],collected_data[i]["molecule_id"]))

    molecule_id=1
    for i in atoms_list:
      output1.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_count1,cur_atoms.type[i],cur_atoms.posx[i],cur_atoms.posy[i],cur_atoms.posz[i],molecule_id))
      temp1={"id":atom_count1,"type":cur_atoms.type[i],"posx":cur_atoms.posx[i],"posy":cur_atoms.posy[i],"posz":cur_atoms.posz[i],"molecule_id":molecule_id}
      temp2={atom_count1:temp1}
      collected_data.update(temp2)
      atom_count1 = atom_count1 + 1
      molecule_id = molecule_id + 1
    output1.close()
