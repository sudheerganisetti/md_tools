#!/usr/bin/python
"""
@Author : Sudheer Ganisetti
@Date   : Wed Nov  9 15:56:04 CET 2016
This code is to compute nnlchkpt
"""
import numpy as np
import sys
""" **************** Main Function **************** """  
if __name__=="__main__":

  tot_argv=len(sys.argv)
  if tot_argv != 4:
     print "************* S. Ganisetti *************"
     print "Error: usage is wrong"
     print "./this_program   cur_chkpt   ref_nnlFile   cur_nnlFile"
     print "this program is to compute nnlchkpt"
     print "adds TotalNeighbours, NewBonds, SwitchedBonds, BrokenBonds & BondStatus information to cur_chkpt"
     print "bond_status values: 0=no_change; 1=new_bond; 2=switched_bond; 3=broken_bond"
     print "****************************************"
     sys.exit(0)
     
  MAX_ATOMS=72000
  MAX_NEIGHBOURS=11
  
  #REF_CHKPT_FILENAME=str(sys.argv[1])
  CUR_CHKPT_FILENAME=str(sys.argv[1])
  REF_NNL_FILENAME=str(sys.argv[2])
  CUR_NNL_FILENAME=str(sys.argv[3])

  OutputFileName1=CUR_CHKPT_FILENAME+str(".nnlchkpt")
  file1=open(OutputFileName1,'w')

  ref_nnl_data=np.loadtxt(REF_NNL_FILENAME,dtype='int')
  cur_nnl_data=np.loadtxt(CUR_NNL_FILENAME,dtype='int')
  
  #ref_atoms=np.loadtxt(REF_CHKPT_FILENAME)
  #ref_atoms_data=open(REF_CHKPT_FILENAME,'r')

  cur_atoms=np.loadtxt(CUR_CHKPT_FILENAME)
  cur_atoms_data=open(CUR_CHKPT_FILENAME,'r')

  mass=[]
  mass.append(0.000030)
  mass.append(0.000017)

  #ref_atom_id=[-1 for i in xrange(MAX_ATOMS)]
  #ref_atom_type=[-1 for i in xrange(MAX_ATOMS)]
  #ref_atom_posx=[-1.0 for i in xrange(MAX_ATOMS)]
  #ref_atom_posy=[-1.0 for i in xrange(MAX_ATOMS)]
  #ref_atom_posz=[-1.0 for i in xrange(MAX_ATOMS)]

  cur_atom_id=[-1 for i in xrange(MAX_ATOMS)]
  cur_atom_type=[-1 for i in xrange(MAX_ATOMS)]
  cur_atom_posx=[-1.0 for i in xrange(MAX_ATOMS)]
  cur_atom_posy=[-1.0 for i in xrange(MAX_ATOMS)]
  cur_atom_posz=[-1.0 for i in xrange(MAX_ATOMS)]

  ref_nnl=[[ -1 for i in xrange(MAX_NEIGHBOURS)] for j in xrange(MAX_ATOMS)]
  cur_nnl=[[ -1 for i in xrange(MAX_NEIGHBOURS)] for j in xrange(MAX_ATOMS)]

  tot_nn=[]
  switched_bonds=[]
  new_bonds=[]
  broken_bonds=[]
  bond_status=[]
  for i in range(MAX_ATOMS):
    tot_nn.append(-1)
    switched_bonds.append(-1)
    new_bonds.append(-1)
    broken_bonds.append(-1)
    bond_status.append(-1)

  #####################################
  # Reading Current Chkpt Box Dimension
  for line in cur_atoms_data:
    data=line.strip().split()
    if data[0] == "#X":
      cur_box_xx=float(data[1])
      cur_box_xy=float(data[2])
      cur_box_xz=float(data[3])
    if data[0] == "#Y":
      cur_box_yx=float(data[1])
      cur_box_yy=float(data[2])
      cur_box_yz=float(data[3])
    if data[0] == "#Z":
      cur_box_zx=float(data[1])
      cur_box_zy=float(data[2])
      cur_box_zz=float(data[3])

  ######################################
  # Raeding Current Chkpt Atoms
  for i in range(MAX_ATOMS):
    global_atom_id=int(cur_atoms[i][0])
    cur_atom_id[global_atom_id]=global_atom_id
    cur_atom_type[global_atom_id]=int(cur_atoms[i][1])
    cur_atom_posx[global_atom_id]=float(cur_atoms[i][3])
    cur_atom_posy[global_atom_id]=float(cur_atoms[i][4])
    cur_atom_posz[global_atom_id]=float(cur_atoms[i][5])

  ######################################
  # Reading Reference nnl Information
  for i in range(len(ref_nnl_data)):
    global_atom_id=int(ref_nnl_data[i][0])
    a1=list(ref_nnl_data[i])
    a2=set(a1)
    # Remove -1 and its own number from the list
    if -1 in a2:
      a2.remove(-1)
    a2.remove(global_atom_id)
    ref_nnl[global_atom_id]=list(a2)
  
  ######################################
  # Reading Current nnl Information
  for i in range(len(cur_nnl_data)):
    global_atom_id=int(cur_nnl_data[i][0])
    a1=list(cur_nnl_data[i])
    a2=set(a1)
    # Remove -1 and its own number from the list
    if -1 in a2:
      a2.remove(-1)
    a2.remove(global_atom_id)
    cur_nnl[global_atom_id]=list(a2)


  #######################################
  # Main loop to compute broken bond information
  for i in range(MAX_ATOMS):
    tot_nn[i]=len(cur_nnl[i])
    a1=set(ref_nnl[i]).difference(set(cur_nnl[i]))
    a2=set(cur_nnl[i]).difference(set(ref_nnl[i]))
    b1=len(a1)
    b2=len(a2)
    if len(ref_nnl[i]) == len(cur_nnl[i]):
      new_bonds[i]=b2
      switched_bonds[i]=b1
      broken_bonds[i]=b1
      if a1:
         bond_status[i]=2
      else:
         bond_status[i]=0
    elif len(ref_nnl[i]) < len(cur_nnl[i]):
      new_bonds[i]=b2
      switched_bonds[i]=b1
      broken_bonds[i]=b1
      bond_status[i]=1
    elif len(cur_nnl[i]) < len(ref_nnl[i]):
      new_bonds[i]=b2
      switched_bonds[i]=b2
      broken_bonds[i]=b1
      bond_status[i]=3

  #######################################
  # writing out to current chkpt
      # writing header
  file1.write("#F A 1 1 1 3 5 0\n#C number type mass x y z tot_nn new_bonds switched_bonds broken_bonds bond_status\n")
  file1.write("#X \t %lf \t %lf \t %lf\n"%(cur_box_xx,cur_box_xy,cur_box_xz))
  file1.write("#Y \t %lf \t %lf \t %lf\n"%(cur_box_yx,cur_box_yy,cur_box_yz))
  file1.write("#Z \t %lf \t %lf \t %lf\n"%(cur_box_zx,cur_box_zy,cur_box_zz))
  file1.write("## computed chkpt using cur_chkpt=%s, ref_nnl=%s & cur_nnl=%s \n"%(CUR_CHKPT_FILENAME,REF_NNL_FILENAME,CUR_NNL_FILENAME))
  file1.write("## bond_status values: 0=no_change; 1=new_bond; 2=switched_bond; 3=broken_bond \n")
  file1.write("#E\n")

  for i in range(MAX_ATOMS):
    file1.write("%d \t %d  %lf \t %lf \t %lf \t %lf \t %d  %d  %d  %d  %d \n" %(cur_atom_id[i],cur_atom_type[i],mass[cur_atom_type[i]],cur_atom_posx[i],cur_atom_posy[i],cur_atom_posz[i],tot_nn[i],new_bonds[i],switched_bonds[i],broken_bonds[i],bond_status[i]))


