
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
  total_atoms=len(atoms_list)
  for j in range(START_CHKPT_NUMBER,END_CHKPT_NUMBER,INTERVAL_CHKPT_NUMBER):
    k=("%08d" %(j))
    CUR_LAMMPS_FILE 	= str(CUR_LAMMPS_FILE_PREFIX)+str(".")+str(k)+str(".dump")
    cur_atoms 			= ganisetti_tools.get_atoms_info_from_lammps(CUR_LAMMPS_FILE)
    output1=open(CUR_LAMMPS_FILE_PREFIX+str("_selected_atoms_from_all_files.")+str(k)+str(".atoms"),'w')
    output1.write("# MSD processed file; ITEM: ATOMS id type x y z\n\n")
    output1.write("%d atoms\n" %(total_atoms))
    output1.write("%d bonds\n" %(total_atoms-len(atoms_list)))
    output1.write("0 angles\n0 dihedrals\n0 impropers\n\n")
    output1.write("10 atom types\n1 bond types\n\n")
    output1.write("%lf  %lf  xlo  xhi\n" %(cur_atoms.box[0][0],cur_atoms.box[0][1]))
    output1.write("%lf  %lf  ylo  yhi\n" %(cur_atoms.box[1][0],cur_atoms.box[1][1]))
    output1.write("%lf  %lf  zlo  zhi\n\n" %(cur_atoms.box[2][0],cur_atoms.box[2][1]))
    output1.write("Atoms \n\n")

    if j != START_CHKPT_NUMBER:
      for i in collected_data.keys():
        output1.write("%d  %d  %d  1.0  %lf  %lf  %lf  \n" %(collected_data[i]["id"],collected_data[i]["molecule_id"],collected_data[i]["type"],collected_data[i]["posx"],collected_data[i]["posy"],collected_data[i]["posz"]))

    molecule_id=1
    for i in atoms_list:
      output1.write("%d  %d  %d  1.0  %lf  %lf  %lf\n" %(atom_count1,molecule_id,cur_atoms.type[i],cur_atoms.posx[i],cur_atoms.posy[i],cur_atoms.posz[i]))
      temp1={"id":atom_count1,"type":cur_atoms.type[i],"posx":cur_atoms.posx[i],"posy":cur_atoms.posy[i],"posz":cur_atoms.posz[i],"molecule_id":molecule_id}
      temp2={atom_count1:temp1}
      collected_data.update(temp2)
      atom_count1 = atom_count1 + 1
      molecule_id = molecule_id + 1

    if j != START_CHKPT_NUMBER:
      bond_count1=1
      output1.write("\nBonds\n\n")
      for i in range(len(atoms_list)+1,atom_count1,1):
        output1.write("%d  1  %d  %d\n" %(bond_count1,i-len(atoms_list),i))
        bond_count1 = bond_count1 + 1
    total_atoms = total_atoms+len(atoms_list)
    output1.close()
