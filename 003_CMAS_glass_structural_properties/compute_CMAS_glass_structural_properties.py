#!/usr/bin/python
"""
@Author : Sudheer Ganisetti
@Date   : Mo 25. Dez 18:34:13 CET 2017
This code is to compute nnl, nnlchk,  connectivity(Qn) and A-O-A triplets for the sample (SiO2 + Al2O3 + CaO + MgO)
"""

import numpy as np
import sys
import math as mt
import itertools as it
import subprocess
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
def write_imd_header1(output,box,atype1,r1,atype2,r2):
  #output
  output.write("#F A 1 1 3 1 0\n")
  output.write("#C number type x y z tot_nn \n")
  output.write("#X %lf 0.0 0.0 \n" %(box[0]))
  output.write("#Y 0.0 %lf 0.0 \n" %(box[1]))
  output.write("#Z 0.0 0.0 %lf \n" %(box[2]))
  output.write("## rcut (%s-O) = %lf and rcut (%s-O) = %lf \n" %(str(atype1),r1,str(atype2),r2))
  output.write("#E \n")

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
def write_imd_header2(output,box,atype1,r1,atype2,r2):
  #output
  output.write("#F A 1 1 3 3 0\n")
  output.write("#C number type x y z tot_nn n_in_QnAlm m_in_QnAlm\n")
  output.write("#X %lf 0.0 0.0 \n" %(box[0]))
  output.write("#Y 0.0 %lf 0.0 \n" %(box[1]))
  output.write("#Z 0.0 0.0 %lf \n" %(box[2]))
  output.write("## rcut (%s-O) = %lf and rcut (%s-O) = %lf \n" %(str(atype1),r1,str(atype2),r2))
  output.write("#E \n")

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
def write_imd_header3(output,box,atype1,r1):
  #output
  output.write("#F A 1 1 3 1 0\n")
  output.write("#C number type x y z tot_nn \n")
  output.write("#X %lf 0.0 0.0 \n" %(box[0]))
  output.write("#Y 0.0 %lf 0.0 \n" %(box[1]))
  output.write("#Z 0.0 0.0 %lf \n" %(box[2]))
  output.write("## rcut of %s-O = %lf, %s-O = %lf, %s-O = %lf, %s-O = %lf \n" %(str(atype1[2]),r1[2][1],str(atype1[3]),r1[3][1],str(atype1[4]),r1[4][1],str(atype1[5]),r1[5][1]))
  output.write("#E \n")

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
def write_imd_header4(output,box,atype1,r1):
  #output
  output.write("#F A 1 1 3 1 1 0\n")
  output.write("#C number type x y z tot_nn status\n")
  output.write("#X %lf 0.0 0.0 \n" %(box[0]))
  output.write("#Y 0.0 %lf 0.0 \n" %(box[1]))
  output.write("#Z 0.0 0.0 %lf \n" %(box[2]))
  output.write("## rcut of %s-O = %lf, %s-O = %lf, %s-O = %lf, %s-O = %lf \n" %(str(atype1[2]),r1[2][1],str(atype1[3]),r1[3][1],str(atype1[4]),r1[4][1],str(atype1[5]),r1[5][1]))
  output.write("#E \n")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if __name__=="__main__":
  
  tot_argv=len(sys.argv)
  if tot_argv != 2:
     #subprocess.call("sudheer_banner")
     print "************* S. Ganisetti *************"
     print "Error: usage is wrong"
     print "./this_program  ChkptFile"
     print "The program is prepared for a maximum of 5 atom types"
     print "By default we are assuming following atoms types and cutoffs"
     print "O(1), Si(2), Al(3), Ca(4), Mg(5) => Si-O= 1.85 ; Al-O= 2.0 ; CaO= 2.9 ; MgO= 2.7 "
     print "Please specify the cutoff radii of all required pair of atoms"
     print "If you want more atom types please modify MAX_ATOM_TYPES and also give corresponding atom pair cutoff radii"
     print "Since the code is not optimized there are alot of unwanted files are getting produced for analysis but these are not needed everytime"
     print "Therefore, enable the code as phases using the following controls:   Enable_phase1_calculations; Enable_phase2_calculations"
     print "****************************************"
     sys.exit(0)

  # --------------------------------------------------------------------
  # The cutoff radii are hard coaded, please edit according to your need
  MAX_NEIGHBOURS=15
  rc_SiO=2.34
  rc_AlO=2.34
  rc_CaO=3.0
  rc_MgO=2.75
  MAX_ATOM_TYPES=5
  Enable_O_env_status_loop="yes"
  Enable_phase1_calculations="no"
  Enable_phase2_calculations="yes"
  atype=["","O","Si","Al","Ca","Mg"]
  rc=[[0.0 for i in xrange(MAX_ATOM_TYPES+1)] for j in xrange(MAX_ATOM_TYPES+1)]
  rc[2][1] = rc[1][2] = rc_SiO
  rc[3][1] = rc[1][3] = rc_AlO
  rc[4][1] = rc[1][4] = rc_CaO
  rc[5][1] = rc[1][5] = rc_MgO
  rcut=0.0
  for i in range(MAX_ATOM_TYPES+1):
    for j in range(MAX_ATOM_TYPES+1):
       rcut=max(rcut,rc[i][j])
  # for debugging
  temp87=0
  if Enable_phase1_calculations == "yes" and Enable_O_env_status_loop == "no":
    print "\n\n*****************************************************************************************************************"
    print "*                                                                                                               *"
    print "*  Please note that computing the type of Oxygen based on its environment is time consuming and is not enabled  *"
    print "*                                                                                                               *"
    print "*****************************************************************************************************************\n\n"
  if Enable_phase1_calculations == "yes" and Enable_O_env_status_loop == "yes":
    print "\n\n*****************************************************************************************************************"
    print "*                                                                                                               *"
    print "* Beware that you have enabled 'Enable_O_env_status_loop'                                                       *"
    print "* computing the type of Oxygen based on its environment is time consuming and requires a lot of memory          *"
    print "* if the memory exceeds the system memory resulting the computer to hang                                        *"
    print "* if you want to avoid your system to hang please watch the memory with 'htop' command in another terminal and  *"
    print "* terminate the program if the memory is exceeding 90% of the systems memory                                    *"
    print "*                                                                                                               *"
    print "*****************************************************************************************************************\n\n"

  print "\n\n*****************************************************************************************************************"
  print "*                                                                                                               *"
  print "* Enable_phase1_calculations = %s										*" %(str(Enable_phase1_calculations))
  print "* Enable_phase2_calculations = %s										*" %(str(Enable_phase2_calculations))
  print "*                                                                                                               *"
  print "*****************************************************************************************************************\n\n"
  


  # --------------------------------------------------------------------

  BASE_FILE=str(sys.argv[1])
  LAMMPS_DUMP_FILE=BASE_FILE+str('.dump')
  foutput1=BASE_FILE+str('.nnl')
  output1=open(foutput1,'w')
  output1.write("# nnl file: used cutoffs are r_SiO = %lf ; r_AlO = %lf ; r_CaO = %lf ; r_MgO = %lf\n" %(rc[1][2],rc[1][3],rc[1][4],rc[1][5]))  

  atoms_data1=np.loadtxt(LAMMPS_DUMP_FILE,dtype='int,int,float,float,float,float,float,float,float',skiprows=9)
  atoms_data2=open(LAMMPS_DUMP_FILE,'r')
  TotalAtoms=len(atoms_data1)
  TotalProperties=len(atoms_data1[0])

  # Find the max id of the atom number
  atoms_data3=[]
  for line in atoms_data1:
    for i in line:
      atoms_data3.append(i)
  atoms_data3=np.reshape(atoms_data3,(TotalAtoms,TotalProperties))

  # define parameters to store atom properties
  MAX_ATOM_NUMBER=int(max(atoms_data3.T[0]))
  atom_status=[0 for i in xrange(TotalAtoms)]
  atom_local_id=[-1 for i in xrange(MAX_ATOM_NUMBER+1)]		# atom_local_id[global_id]=local_id
  atom_global_id=[-1 for i in xrange(TotalAtoms)]		# atom_global_id[local_id]=global_id
  atom_type=[-1 for i in xrange(TotalAtoms)]
  atom_charge=[0.0 for i in xrange(TotalAtoms)]
  atom_posx=[-1.0 for i in xrange(TotalAtoms)]
  atom_posy=[-1.0 for i in xrange(TotalAtoms)]
  atom_posz=[-1.0 for i in xrange(TotalAtoms)]
  
  # define parameters to store nnl related information
  nnl=[[-1 for j in xrange(MAX_NEIGHBOURS)] for i in xrange(TotalAtoms)]
  nnl_count=[0 for i in xrange(TotalAtoms)]
  #storing atom types of each neighbour in the nnl_types for getting things easy to work with
  nnl_type=[[-1 for j in xrange(MAX_NEIGHBOURS)] for i in xrange(TotalAtoms)]
  #store each pair distance to use them in later
  nnl_each_pair_distance=[[-1 for j in xrange(MAX_NEIGHBOURS)] for i in xrange(TotalAtoms)]


  # Reading box dimensions
  count=1
  for line in atoms_data2:
    data=line.strip().split()
    if count < 10:
      if count == 4:
        if TotalAtoms != int(data[0]):
          print "WARNING: The total number of atoms mentioned in the file is not equal to the total number of atoms in the file\n"
      elif count == 6:
        box_xx_lo=float(data[0])
        box_xx_hi=float(data[1])
        box_xx=box_xx_hi-box_xx_lo
      elif count == 7:
        box_yy_lo=float(data[0])
        box_yy_hi=float(data[1])
        box_yy=box_yy_hi-box_yy_lo
      elif count == 8:
        box_zz_lo=float(data[0])
        box_zz_hi=float(data[1])
        box_zz=box_zz_hi-box_zz_lo
    count=count+1
  imd_box=np.array([box_xx,box_yy,box_zz])

  # Reading atoms position
  count=0
  for i in range(TotalAtoms):
    atom_local_id[atoms_data1[i][0]]=i				# atom_local_id[global_id]=local_id
    atom_global_id[i]=atoms_data1[i][0]				# atom_global_id[local_id]=global_id
    atom_type[i]=atoms_data1[i][1]
    atom_charge[i]=atoms_data1[i][2]
    # atom positions are not shifted if the box is not shifted
    #atom_posx[i]=atoms_data1[i][3]
    #atom_posy[i]=atoms_data1[i][4]
    #atom_posz[i]=atoms_data1[i][5]
    # atom positions are also shifted if the box is shifted
    atom_posx[i]=atoms_data1[i][3] - box_xx_lo
    atom_posy[i]=atoms_data1[i][4] - box_yy_lo
    atom_posz[i]=atoms_data1[i][5] - box_zz_lo

  # Devide the box into voxels
  MaxCells_X=int(box_xx/rcut)
  MaxCells_Y=int(box_yy/rcut)
  MaxCells_Z=int(box_zz/rcut)

  # LEC= Length of Each Cell
  LEC_X=float(box_xx/MaxCells_X)
  LEC_Y=float(box_yy/MaxCells_Y)
  LEC_Z=float(box_zz/MaxCells_Z)
  
  cell_atoms_count=[[[0 for i in xrange(MaxCells_Z)] for j in xrange(MaxCells_Y)] for k in xrange(MaxCells_X)]				# cell_atoms_count[X][Y][Z]
  cell_atoms_list=[[[[0 for i in xrange(200)] for j in xrange(MaxCells_Z)] for k in xrange(MaxCells_Y)] for l in xrange(MaxCells_X)]	# cell_atoms_list[X][Y][Z][200]

  #assign the atoms into respective cell number
  for i in range(TotalAtoms):
    # If the atoms are outside of the box then bring them into other side according to PBC
    if atom_posx[i] >= box_xx:
      atom_posx[i]=atom_posx[i]-box_xx
    if atom_posy[i] >= box_yy:
      atom_posy[i]=atom_posy[i]-box_yy
    if atom_posz[i] >= box_zz:
      atom_posz[i]=atom_posz[i]-box_zz
    
    cellX=int(atom_posx[i]/LEC_X)
    cellY=int(atom_posy[i]/LEC_Y)
    cellZ=int(atom_posz[i]/LEC_Z)
    cell_atoms_list[cellX][cellY][cellZ][cell_atoms_count[cellX][cellY][cellZ]]=i
    cell_atoms_count[cellX][cellY][cellZ]=cell_atoms_count[cellX][cellY][cellZ]+1

  total_atoms_of_type=[0 for i in xrange(MAX_ATOM_TYPES+1)]
  for i in xrange(TotalAtoms):
    if atom_type[i] == 1:
      total_atoms_of_type[1] = total_atoms_of_type[1] + 1
    elif atom_type[i] == 2:
      total_atoms_of_type[2] = total_atoms_of_type[2] + 1
    elif atom_type[i] == 3:
      total_atoms_of_type[3] = total_atoms_of_type[3] + 1
    elif atom_type[i] == 4:
      total_atoms_of_type[4] = total_atoms_of_type[4] + 1
    elif atom_type[i] == 5:
      total_atoms_of_type[5] = total_atoms_of_type[5] + 1

  # Main loop to compute nnl and other information
  # (i,j,k) => selecting atom which is belongs to one cell
  for i in range(MaxCells_X):
    for j in range(MaxCells_Y):
      for k in range(MaxCells_Z):
        for atom1 in range(cell_atoms_count[i][j][k]):
          atom1_id=cell_atoms_list[i][j][k][atom1]
          pos_x1=atom_posx[atom1_id]
          pos_y1=atom_posy[atom1_id]
          pos_z1=atom_posz[atom1_id]

          # (l,m,n) => selecting neighbour cell 
          for l in np.arange(0,2,1):
            new_i=i+l
            an=new_i
            pbc_x=0.0
            # -----------------------------------------------------
            # Dealing with periodic boundary condition along x-axis
            if new_i < 0:
              # But the conditions will not lead to excute this case
              new_i=MaxCells_X-1
              print "WARNING!!! There is somethong wrong, please check it well !!!\n"
            if new_i > MaxCells_X-1:
              new_i=0
              pbc_x=box_xx
            # -----------------------------------------------------
            for m in np.arange(-1*l,2,1):
              new_j=j+m
              bn=new_j
              pbc_y=0.0
              # -----------------------------------------------------
              # Dealing with periodic boundary condition along y-axis
              if new_j < 0:
                new_j=MaxCells_Y-1
                pbc_y=-1.0*box_yy
              if new_j > MaxCells_Y-1:
                new_j=0
                pbc_y=box_yy
              # -----------------------------------------------------
              n1=-1
              if l==m:
                n1=-1*l
              for n in np.arange(n1,2,1):
                new_k=k+n
                cn=new_k
                pbc_z=0.0
                # -----------------------------------------------------
                # Dealing with periodic boundary conditions along z-axis
                if new_k < 0:
                  new_k=MaxCells_Z-1
                  pbc_z=-1.0*box_zz
                if new_k > MaxCells_Z-1:
                  new_k = 0
                  pbc_z=box_zz
                # debugging
                #if i==0 and j==0 and k==0:
                #  print ("%d %d %d %d %d %d   atom1_id %d") %(l,m,n,an,bn,cn,atom1_id)

                # -----------------------------------------------------
                # This is a trick to avoid double counting the interactions if the two atoms are in the same cell
                temp=0
                if i == new_i and j == new_j and k == new_k:
                  temp=atom1+1
                # selecting the neighbour atom
                for atom2 in range(temp,cell_atoms_count[new_i][new_j][new_k]):
                  atom2_id=cell_atoms_list[new_i][new_j][new_k][atom2]
                  pos_x2=atom_posx[atom2_id] + pbc_x
                  pos_y2=atom_posy[atom2_id] + pbc_y
                  pos_z2=atom_posz[atom2_id] + pbc_z

                  r=np.linalg.norm(np.array([pos_x1,pos_y1,pos_z1])-np.array([pos_x2,pos_y2,pos_z2]))
                  if r <= rc[atom_type[atom1_id]][atom_type[atom2_id]]:
                    nnl[atom1_id][nnl_count[atom1_id]]=atom2_id
                    nnl[atom2_id][nnl_count[atom2_id]]=atom1_id
                    nnl_each_pair_distance[atom1_id][nnl_count[atom1_id]]=r
                    nnl_each_pair_distance[atom2_id][nnl_count[atom2_id]]=r
                    nnl_count[atom1_id]=nnl_count[atom1_id]+1
                    nnl_count[atom2_id]=nnl_count[atom2_id]+1
                    # debugging 
                    #if atom1_id==713:
                    #  print "-------------"
                    #  print "%d %d %d %d %d %d %d %d\n" %(atom1_id,i,j,k,atom2_id,new_i,new_j,new_k)
                    #  print "-------------"
                    #if i==0 and j==0 and k==0:
                    #   print ("%d %d %d %d %d %d") %(i,j,k,new_i,new_j,new_k)
                    #temp80=nnl[atom1_id]
                    #temp81=nnl[atom2_id]
                    #temp82=[]
                    #for temp83 in temp80:
                    #  if temp83 != -1:
                    #    temp82.append(temp83)
                    #temp84=[]
                    #for temp85 in temp81:
                    #  if temp85 != -1:
                    #    temp84.append(temp85)
                    #if (len(temp82) != len(set(temp82))) and temp87 ==-1 :
                    #  print "%d %d %d %d %d %d %d %d\n" %(atom1_id,i,j,k,atom2_id,new_i,new_j,new_k)
                    #  print "max cells %d %d %d" %(MaxCells_X,MaxCells_Y,MaxCells_Z)
                    #  print "atoms in 1st cell "
                    #  for temp86 in range(cell_atoms_count[i][j][k]):
                    #      print "%d  " %(cell_atoms_list[i][j][k][temp86])
                    #  print "nnl of %d" %(atom1_id)
                    #  for temp86 in nnl[atom1_id]:
                    #    if temp86 != -1:
                    #      print "%d " %(temp86)
                    #  print "atoms in 2nd cell "
                    #  for temp86 in range(cell_atoms_count[new_i][new_j][new_k]):
                    #      print "%d  " %(cell_atoms_list[new_i][new_j][new_k][temp86])
                    #  print "nnl of %d" %(atom2_id)
                    #  for temp86 in nnl[atom2_id]:
                    #    if temp86 != -1:
                    #      print "%d " %(temp86)
                    #  temp87=temp87+1
  # **************************************************************************************
  # removing -1 from nnl[i]
  atom_nnl_type=[[-1 for j in xrange(MAX_NEIGHBOURS)] for i in xrange(TotalAtoms)]
  for i in range(TotalAtoms):
    temp1=[]
    temp2=[]
    temp3=[]
    for j in range(len(nnl[i])):
      if nnl[i][j] !=-1:
        temp1.append(nnl[i][j])
        temp2.append(atom_type[nnl[i][j]])
        temp3.append(nnl_each_pair_distance[i][j])
    nnl[i]=temp1
    nnl_type[i]=temp2
    nnl_each_pair_distance[i]=temp3
    for j in xrange(nnl_count[i]):
      atom_nnl_type[i][j]=atom_type[nnl[i][j]]
  # **************************************************************************************
  # writing out all nnl list
  for i in range(0,TotalAtoms):
    output1.write("%d " %(atom_global_id[i]))
    for j in range(MAX_NEIGHBOURS):
      if j < nnl_count[i]:
        output1.write("%d  " %(atom_global_id[nnl[i][j]]))
      else:
        output1.write("-1  ")
    output1.write("\n")
  output1.close()

  # **************************************************************************************
  # writing out chkpt with nnl status
  output=open(BASE_FILE+str("_withnnlstatus.atoms"),'w')
  write_imd_header3(output,imd_box,atype,rc)
  for i in range(TotalAtoms):
    output.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
  output.close()

  # **************************************************************************************
  # writing out individual cations and their anions
  for i in range(2,MAX_ATOM_TYPES+1):
    output=open(BASE_FILE+str("_")+str(atype[i])+str("-O.atoms"),'w')
    write_imd_header1(output,imd_box,str(atype[i]),rc[i][1],str(atype[i]),rc[i][1])
    SelectAtoms1=[0 for ii in xrange(TotalAtoms)]
    for j in range(TotalAtoms):
      if i == atom_type[j]:
        SelectAtoms1[j]=1
        for k in range(nnl_count[j]):
          SelectAtoms1[nnl[j][k]]=1
    for j in range(TotalAtoms):
      if SelectAtoms1[j] == 1:
        output.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
    output.close()

  if Enable_phase1_calculations == "yes":
    # **************************************************************************************
    # writing out bridging and non-bridging oxygen
    NBOS=[0 for i in xrange(TotalAtoms)]									# Non-Bridging Oxygen Status
    OCa_count=[0 for i in xrange(TotalAtoms)]								# Number of Ca connected to each oxygen
    OMg_count=[0 for i in xrange(TotalAtoms)]								# Number of Mg connected to each oxygen
    BOCa_2Si_count=0											# Number of Ca connected to bridging oxygen and 2 Si
    BOCa_2Al_count=0											# Number of Ca connected to bridging oxygen and 2 Al
    BOCa_SiAl_count=0											# Number of Ca connected to bridging oxygen and 1 Si+  1Al
    BOMg_2Si_count=0                                                                                      # Number of Mg connected to bridging oxygen and 2 Si
    BOMg_2Al_count=0                                                                                      # Number of Mg connected to bridging oxygen and 2 Al
    BOMg_SiAl_count=0                                                                                     # Number of Mg connected to bridging oxygen and 1 Si+  1Al
    NBOCa_Si_count=0                                                                                      # Number of Ca connected to non-bridging oxygen and 2 Si
    NBOCa_Al_count=0                                                                                      # Number of Ca connected to non-bridging oxygen and 2 Al
    NBOMg_Si_count=0
    NBOMg_Al_count=0
    TRI_3Si_0Al=TRI_2Si_1Al=TRI_1Si_2Al=TRI_0Si_3Al=0
    TRI_1Si_2Al_list=[]
    TRI_0Si_3Al_list=[]
    TRI_Al_atoms=[-1 for i in xrange(TotalAtoms)]                                                       	# Al atoms involved in tri-clusters
    TRI_1Si_2Al_tet_atoms=[-1 for i in xrange(TotalAtoms)]						# tri-cluster tetra atoms
    TRI_0Si_3Al_tet_atoms=[-1 for i in xrange(TotalAtoms)]						# tri-cluster tetra atoms
    NBOSi_count	=NBOAl_count	=NBO3Si_count	=0
    BOSiOSi_count	=BOSiOAl_count	=BOAlOAl_count	=0
    NBOI=[[[[0 for i in xrange(8)] for j in xrange(8)] for k in xrange(8)] for l in xrange(8)]		# Non-Bridging Oxygen Information
    NBOI1=[[[0 for i in xrange(8)] for j in xrange(8)] for k in xrange(8)]				# Non-bridging Oxygen Information
    BOI=[[[[0 for i in xrange(8)] for j in xrange(8)] for k in xrange(8)] for l in xrange(8)]             # Bridging Oxygen Information
    BOI1=[[[0 for i in xrange(8)] for j in xrange(8)] for k in xrange(8)]                                 # Bridging Oxygen Information
    AlOAl_NBO_Count1=[0 for i in xrange(10)]
    AlOAl_NBO_Count2=[0 for i in xrange(10)]
    AlOAl_NBO_Count3=[0 for i in xrange(10)]
    AlOAl_NBO_Count4=[0 for i in xrange(10)]
    output1=open(BASE_FILE+str("_NBO_info.data"),'w')
    output2=open(BASE_FILE+str("_BO_info.data"),'w')
    output3=open(BASE_FILE+str("_BO_additional_info.data"),'w')
    output4=open(BASE_FILE+str("_tri-cluster_info.data"),'w')
    output5=open(BASE_FILE+str("_tri-cluster_1Si2Al_OnlyOxygen.atoms"),'w')
    output6=open(BASE_FILE+str("_tri-cluster_1Si2Al_totaltets.atoms"),'w')
    output7=open(BASE_FILE+str("_tri-cluster_0Si3Al_OnlyOxygen.atoms"),'w')
    output8=open(BASE_FILE+str("_tri-cluster_0Si3Al_totaltets.atoms"),'w')
    output11=open(BASE_FILE+str("_NBO.atoms"),'w')
    output12=open(BASE_FILE+str("_NBO_and_its_tetra.atoms"),'w')
    output21=open(BASE_FILE+str("_BO.atoms"),'w')
    output22=open(BASE_FILE+str("_BO_and_its_tetra.atoms"),'w')
    output31=open(BASE_FILE+str("_BO_with2Si.atoms"),'w')
    output32=open(BASE_FILE+str("_BO_with2Si_and_its_tetra.atoms"),'w')
    output33=open(BASE_FILE+str("_BO_with2Si0Al0R.atoms"),'w')
    output41=open(BASE_FILE+str("_BO_with2Al.atoms"),'w')
    output42=open(BASE_FILE+str("_BO_with2Al_and_it_tetra.atoms"),'w')
    output51=open(BASE_FILE+str("_BO_with1Si1Al.atoms"),'w')
    output52=open(BASE_FILE+str("_BO_with1Si1Al_and_it_tetra.atoms"),'w')
    output61=open(BASE_FILE+str("_NBO_with0Si1AlnR.atoms"),'w')
    output62=open(BASE_FILE+str("_NBO_with0Si1AlnR_and_its_tetra.atoms"),'w')
    write_imd_header3(output11,imd_box,atype,rc)
    write_imd_header3(output12,imd_box,atype,rc)
    write_imd_header3(output21,imd_box,atype,rc)
    write_imd_header3(output22,imd_box,atype,rc)
    write_imd_header3(output31,imd_box,atype,rc)
    write_imd_header3(output32,imd_box,atype,rc)
    write_imd_header3(output33,imd_box,atype,rc)
    write_imd_header3(output41,imd_box,atype,rc)
    write_imd_header3(output42,imd_box,atype,rc)
    write_imd_header3(output51,imd_box,atype,rc)
    write_imd_header3(output52,imd_box,atype,rc)
    write_imd_header3(output61,imd_box,atype,rc)
    write_imd_header3(output62,imd_box,atype,rc)
    write_imd_header3(output5,imd_box,atype,rc)
    write_imd_header3(output6,imd_box,atype,rc)
    write_imd_header3(output7,imd_box,atype,rc)
    write_imd_header3(output8,imd_box,atype,rc)
    output1.write("# non-bridging_oxygen information \n")
    output1.write("# Si	Al	Ca	Mg	Number_Of_Units \n")
    output2.write("# bridging_oxygen information \n")
    output2.write("# Si   Al      Ca      Mg      Number_Of_Units \n")
    for i in range(TotalAtoms):
      #if atom_type[i] == 1 and nnl_count[i] != 2 :    
      if atom_type[i] == 1 :
        neighbour_atom_types=[]
        for j in nnl[i]:
          neighbour_atom_types.append(atom_type[j])
        #if (4 in set(neighbour_atom_types)) or (5 in set(neighbour_atom_types)):
        #  output1.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
        #else:
        #  output2.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
        number_of_Si=neighbour_atom_types.count(2)
        number_of_Al=neighbour_atom_types.count(3)
        number_of_Ca=neighbour_atom_types.count(4)
        number_of_Mg=neighbour_atom_types.count(5)
        if number_of_Si+number_of_Al < 2:
          output11.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          output12.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          NBOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]=NBOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]+1
          NBOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]=NBOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]+1
          for j in nnl[i]:
            output11.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
            output12.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
            if atom_type[j] == 2 or atom_type[j] == 3:
              for k in nnl[j]:
                if k != i :
                  output12.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
          NBOS[i]=1
          OCa_count[i]=number_of_Ca
          OMg_count[i]=number_of_Mg
          if number_of_Si == 1:
            NBOCa_Si_count=NBOCa_Si_count+number_of_Ca
            NBOMg_Si_count=NBOMg_Si_count+number_of_Mg
            NBOSi_count=NBOSi_count+1
          elif number_of_Al == 1:
            NBOCa_Al_count=NBOCa_Al_count+number_of_Ca
            NBOMg_Al_count=NBOMg_Al_count+number_of_Mg
            NBOAl_count=NBOAl_count+1
        elif number_of_Si == 2:
          output21.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          output31.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          output32.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          BOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]=BOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]+1
          BOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]=BOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]+1
          for j in nnl[i]:
            output21.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
            output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
            output31.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
            output32.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
            if atom_type[j] == 2:
              for k in nnl[j]:
                if k != i :
                  output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
                  output32.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
                  if number_of_Al == 0 and number_of_Ca == 0 and number_of_Mg == 0:
                    output33.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
            if number_of_Al == 0 and number_of_Ca == 0 and number_of_Mg == 0:
              output33.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
          if number_of_Al == 0 and number_of_Ca == 0 and number_of_Mg == 0:
            output33.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          OCa_count[i]=number_of_Ca
          OMg_count[i]=number_of_Mg
          BOCa_2Si_count=BOCa_2Si_count+number_of_Ca
          BOMg_2Si_count=BOMg_2Si_count+number_of_Mg
          BOSiOSi_count=BOSiOSi_count+1
        elif number_of_Al == 2:
          output21.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          output41.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          output42.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          BOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]=BOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]+1
          BOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]=BOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]+1
          for j in nnl[i]:
            output21.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
            output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
            output41.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
            output42.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
            if atom_type[j] == 3:
              for k in nnl[j]:
                if k != i :
                  output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
                  output42.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
          OCa_count[i]=number_of_Ca
          OMg_count[i]=number_of_Mg
          BOCa_2Al_count=BOCa_2Al_count+number_of_Ca
          BOMg_2Al_count=BOMg_2Al_count+number_of_Mg
          BOAlOAl_count=BOAlOAl_count+1
        elif number_of_Si == 1 and number_of_Al == 1:
          output21.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          output51.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          output52.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          BOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]=BOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]+1
          BOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]=BOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]+1
          for j in nnl[i]:
            output21.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
            output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
            output51.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
            output52.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
            if atom_type[j] == 2 or atom_type[j] == 3:
              for k in nnl[j]:
                if k != i :
                  output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
                  output52.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
          OCa_count[i]=number_of_Ca
          OMg_count[i]=number_of_Mg
          BOCa_SiAl_count = BOCa_SiAl_count + number_of_Ca
          BOMg_SiAl_count = BOMg_SiAl_count + number_of_Mg
          BOSiOAl_count=BOSiOAl_count+1
        elif number_of_Si+number_of_Al > 2:
          NBO3Si_count=NBO3Si_count+1
          #if number_of_Al == 3 or number_of_Al=2 and :  
          NBOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]=NBOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]+1
          NBOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]=NBOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]+1
          NBOS[i]=1
        # this loop is to compute tri-cluster information
        if number_of_Si == 3 and number_of_Al == 0:
          TRI_3Si_0Al = TRI_3Si_0Al+1
        elif number_of_Si == 2 and number_of_Al == 1:
          TRI_2Si_1Al = TRI_2Si_1Al+1
        elif number_of_Si == 1 and number_of_Al == 2:
          TRI_1Si_2Al = TRI_1Si_2Al+1
          TRI_1Si_2Al_list.append(i)
          output5.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          output6.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          for j in nnl[i]:
            if atom_type[j] == 2 or atom_type[j] == 3:
              TRI_1Si_2Al_tet_atoms[j]=1
              for k in nnl[j]:
                TRI_1Si_2Al_tet_atoms[k]=1
            if atom_type[j] == 3:
              TRI_Al_atoms[j]=1
        elif number_of_Si == 0 and number_of_Al == 3:
          TRI_0Si_3Al = TRI_0Si_3Al+1
          TRI_0Si_3Al_list.append(i)
          output7.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          output8.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          for j in nnl[i]:
            if atom_type[j] == 2 or atom_type[j] == 3:
              TRI_0Si_3Al_tet_atoms[j]=1
              for k in nnl[j]:
                TRI_0Si_3Al_tet_atoms[k]=1
            if atom_type[j] == 3:
              TRI_Al_atoms[j]=1
        
    for i in range(8):
      for j in range(8):
        for k in range(8):
          for l in range(8):
            if NBOI[i][j][k][l] != 0:
              output1.write("%d Si + %d Al + %d Ca + %d Mg = %d \n" %(i,j,k,l,NBOI[i][j][k][l]))
            if BOI[i][j][k][l] != 0:
              output2.write("%d Si + %d Al + %d Ca + %d Mg = %d \n" %(i,j,k,l,BOI[i][j][k][l]))
    # counting total Al invovled in tri-clusters
    count1=0
    for i in range(TotalAtoms):
      if TRI_Al_atoms[i] == 1:
        count1=count1+1
    # write the tri-cluster information
    output4.write("# tri-cluster information\n")
    output4.write("0 Si \t 3 Al = %d \n" %(TRI_0Si_3Al))
    output4.write("1 Si \t 2 Al = %d \n" %(TRI_1Si_2Al))
    output4.write("2 Si \t 1 Al = %d \n" %(TRI_2Si_1Al))
    output4.write("3 Si \t 0 Al = %d \n" %(TRI_3Si_0Al))
    output4.write("\n\nTotal Al invovled in tri-clusters (O-(Al+Al+Al), O-(Al+Al+Si)) = %d \n" %(count1))
    output4.close()

    # write out the tri-cluster atoms data
    for i in range(TotalAtoms):
      if TRI_1Si_2Al_tet_atoms[i] == 1:
        output6.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
      if TRI_0Si_3Al_tet_atoms[i] == 1:
        output8.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
    output5.close()
    output6.close()
    output7.close()
    output8.close()

    # computing total Ca and Mg connected to Non-Bridiging oxygen and Bridging Oxygen
    sum1=0
    sum2=0
    sum3=0
    sum4=0
    sum11=0
    sum31=0
    for i in range(TotalAtoms):
      # for Non bridging oxygen
      if atom_type[i] == 1 and NBOS[i]==1 :
        sum1=sum1+OCa_count[i]
        sum2=sum2+OMg_count[i]
        sum11=sum11+1
      # For bridging oxygen
      elif atom_type[i] == 1 and NBOS[i]==0 :
        sum3=sum3+OCa_count[i]
        sum4=sum4+OMg_count[i]
        sum31=sum31+1

    output1.write("\n\n\n\n***************************************\n")
    output1.write("# Si	Al 	(R=either Ca or Mg)	Number_of_units\n")
    output2.write("\n\n\n\n***************************************\n")
    output2.write("# Si   Al      (R=either Ca or Mg)     Number_of_units\n")
    for i in range(8):
      for j in range(8):
        for k in range(8):
          if NBOI1[i][j][k] != 0:
            output1.write("%d Si + %d Al + %d (R) = %d \n" %(i,j,k,NBOI1[i][j][k]))
          if BOI1[i][j][k] != 0:
            output2.write("%d Si + %d Al + %d (R) = %d \n" %(i,j,k,BOI1[i][j][k]))

    output1.write("\n\nTotal number of Non-Bridging oxygen atoms 		= %d \n" %(sum11))
    output1.write("Number of Ca Connected to Non-Bridging oxygen  	= %d \n" %(sum1))
    output1.write("Number of Mg Connected to Non-Bridging oxygen		= %d \n" %(sum2))
    output1.write("Number of Ca Connected to Non-Bridging oxygen (Si) 	= %d \n" %(NBOCa_Si_count))
    output1.write("Number of Ca Connected to Non-Bridging oxygen (Al) 	= %d \n" %(NBOCa_Al_count))
    output1.write("Number of Mg Connected to Non-Bridging oxygen (Si) 	= %d \n" %(NBOMg_Si_count))
    output1.write("Number of Mg Connected to Non-Bridging oxygen (Al) 	= %d \n" %(NBOMg_Al_count))
    output1.write("\n\nNon-Bridging Oxygen SiO 				= %d \n" %(NBOSi_count))
    output1.write("Non-Bridging Oxygen AlO                            	= %d \n" %(NBOAl_count))
    output1.write("Non-Bridging Oxygen 3(Si+Al)                           = %d \n" %(NBO3Si_count))

    output2.write("\n\nTotal number of Bridging oxygen atoms		= %d \n" %(sum31))
    output2.write("Number of Ca Connected to Bridging oxygen		= %d \n" %(sum3))
    output2.write("Number of Mg Connected to Bridging oxygen 		= %d \n" %(sum4))
    output2.write("Number of Ca Connected to Bridging oxygen (2Si)	= %d \n" %(BOCa_2Si_count))
    output2.write("Number of Ca Connected to Bridging oxygen (1Si+1Al)	= %d \n" %(BOCa_SiAl_count))
    output2.write("Number of Ca Connected to Bridging oxygen (2Al)	= %d \n" %(BOCa_2Al_count))
    output2.write("Number of Mg Connected to Bridging oxygen (2Si)	= %d \n" %(BOMg_2Si_count))
    output2.write("Number of Mg Connected to Bridging oxygen (1Si+1Al)	= %d \n" %(BOMg_SiAl_count))
    output2.write("Number of Mg Connected to Bridging oxygen (2Al) 	= %d \n" %(BOMg_2Al_count))
    output2.write("\n\nBridging Oxygen SiOSi                            	= %d \n" %(BOSiOSi_count))
    output2.write("Bridging Oxygen SiOAl					= %d \n" %(BOSiOAl_count))
    output2.write("Bridging Oxygen AlOAl					= %d \n" %(BOAlOAl_count))

    # computing number of non-bridging oxygens of the cations which are connected to O in Al-O-Al
    for i in range(TotalAtoms):
      if atom_type[i] == 1 :
        neighbour_atom_types=[]
        Si_atoms=[]
        Al_atoms=[]
        Ca_atoms=[]
        Mg_atoms=[]
        for j in nnl[i]:
          neighbour_atom_types.append(atom_type[j])
          if atom_type[j] == 2:
            Si_atoms.append(j)
          elif atom_type[j] == 3:
            Al_atoms.append(j)
          elif atom_type[j] == 4:
            Ca_atoms.append(j)
          elif atom_type[j] == 5:
            Mg_atoms.append(j)
        number_of_Si=neighbour_atom_types.count(2)
        number_of_Al=neighbour_atom_types.count(3)
        number_of_Ca=neighbour_atom_types.count(4)
        number_of_Mg=neighbour_atom_types.count(5)
        if number_of_Al == 2 and number_of_Si == 0:
          if (number_of_Ca == 1 or number_of_Ca == 2) and number_of_Mg == 0:
            for j in Ca_atoms:
              NBO_Count1=0
              for k in nnl[j]:
                NBO_Count1=NBO_Count1+NBOS[k]
                #if NBOS[k]==1:
                #  NBO_Count1=NBO_Count1+1
            AlOAl_NBO_Count1[NBO_Count1]=AlOAl_NBO_Count1[NBO_Count1]+1
          elif number_of_Ca == 1  and number_of_Mg == 1:
            for j in Ca_atoms:
              NBO_Count2=0
              for k in nnl[j]:
                NBO_Count2=NBO_Count2+NBOS[k]
            AlOAl_NBO_Count2[NBO_Count2]=AlOAl_NBO_Count2[NBO_Count2]+1
            for j in Mg_atoms:
              NBO_Count3=0
              for k in nnl[j]:
                NBO_Count3=NBO_Count3+NBOS[k]
            AlOAl_NBO_Count3[NBO_Count3]=AlOAl_NBO_Count3[NBO_Count3]+1
          elif number_of_Ca == 0  and number_of_Mg == 1:
            for j in Mg_atoms:
              NBO_Count4=0
              for k in nnl[j]:
                NBO_Count4=NBO_Count4+NBOS[k]
            AlOAl_NBO_Count4[NBO_Count4]=AlOAl_NBO_Count4[NBO_Count4]+1

        # Print Al non-bridging oxigen and its neighbourings
        if number_of_Al == 1 and number_of_Si == 0 and (number_of_Ca+number_of_Mg) > 0 :
          output61.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          output62.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          for j in nnl[i]:
            output61.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
            output62.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
            if atom_type[j] == 3:
               for k in nnl[j]:
                 if k != i:
                   output62.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
  
    output3.write("# additional information of Al-O-Al links\n")
    output3.write("# omputing number of non-bridging oxygens of the cations which are connected to O in Al-O-Al \n")
    output3.write("# total number of oxygen atoms (in the links Al-O-Al) such that x non-bridging oxygens are connected to the cations linked to oxygen in Al-O-Al links\n")
    output3.write("# 0 Si + 2 Al + 1 or Ca + 0 Mg ==> \n")
    for i in range(10):
      if AlOAl_NBO_Count1[i] != 0:
        output3.write("\t \t %d NBO = %d \n" %(i,AlOAl_NBO_Count1[i]))
    output3.write("# 0 Si + 2 Al + 1 Ca + 1 Mg ==> \n     Ca nonbridging oxygens \n")
    for i in range(10):
      if AlOAl_NBO_Count2[i] != 0:
        output3.write("\t \t %d NBO = %d \n" %(i,AlOAl_NBO_Count2[i]))
    output3.write("# 0 Si + 2 Al + 1 Ca + 1 Mg ==> \n     Mg nonbridging oxygens \n")
    for i in range(10):
      if AlOAl_NBO_Count3[i] != 0:
        output3.write("\t \t %d NBO = %d \n" %(i,AlOAl_NBO_Count3[i]))
    output3.write("# 0 Si + 2 Al + 0 or Ca + 1 Mg ==> \n")
    for i in range(10):
      if AlOAl_NBO_Count4[i] != 0:
        output3.write("\t \t %d NBO = %d \n" %(i,AlOAl_NBO_Count4[i]))
 
    output1.close()
    output2.close()
    output3.close()
    output11.close()
    output12.close()
    output21.close()
    output22.close()
    output31.close()
    output32.close()
    output33.close()
    output41.close()
    output42.close()
    output51.close()
    output52.close()
    output61.close()
    output62.close()
    """
    # writing out triplets
    output2=open(BASE_FILE+str("_SiOSi.data"),'w')
    output3=open(BASE_FILE+str("_AlOAl.data"),'w')
    output4=open(BASE_FILE+str("_SiOAl.data"),'w')
    output5=open(BASE_FILE+str("_SiOCa.data"),'w')
    output6=open(BASE_FILE+str("_SiOMg.data"),'w')
  
    for i in range(0,TotalAtoms):
      #neighbour_atom_types=[]
      if atom_type[i] == 1:
         #for j in nnl[i]:
         #  if j != -1:
         #    neighbour_atom_types.append(atom_type[j])
         #neighbour_atom_types_set=set(neighbour_atom_types)
         if nnl_count[i]==2:
           # writing out Si-O-Si
           if atom_type[nnl[i][0]] == 2 and atom_type[nnl[i][1]] == 2:
              output2.write("%d %d %d\n" %(atom_global_id[i],atom_global_id[nnl[i][0]],atom_global_id[nnl[i][1]]))
           # writing out Al-O-Al
           elif atom_type[nnl[i][0]] == 3 and atom_type[nnl[i][1]] == 3:
              output3.write("%d %d %d\n" %(atom_global_id[i],atom_global_id[nnl[i][0]],atom_global_id[nnl[i][1]]))
           # writing out Si-O-Al
           elif (atom_type[nnl[i][0]] == 2 and atom_type[nnl[i][1]] == 3 ) or (atom_type[nnl[i][0]] == 3 and atom_type[nnl[i][1]] == 2 ):
              output4.write("%d %d %d\n" %(atom_global_id[i],atom_global_id[nnl[i][0]],atom_global_id[nnl[i][1]]))
           # writing out Si-O-Ca
           elif (atom_type[nnl[i][0]] == 2 and atom_type[nnl[i][1]] == 4 ) or (atom_type[nnl[i][0]] == 4 and atom_type[nnl[i][1]] == 2 ):
              output5.write("%d %d %d\n" %(atom_global_id[i],atom_global_id[nnl[i][0]],atom_global_id[nnl[i][1]]))
           # writing out Si-O-Mg
           elif (atom_type[nnl[i][0]] == 2 and atom_type[nnl[i][1]] == 5 ) or (atom_type[nnl[i][0]] == 5 and atom_type[nnl[i][1]] == 2 ):
              output6.write("%d %d %d\n" %(atom_global_id[i],atom_global_id[nnl[i][0]],atom_global_id[nnl[i][1]]))
         else :
           print "Oxygen atom %d have nnls %d" %(atom_global_id[i],nnl_count[i])
    """
    # **************************************************************************************
    # writing out triplets
    triplet=[[ 0 for i in xrange(0,MAX_ATOM_TYPES+1)] for j in xrange(0,MAX_ATOM_TYPES+1)]
    triplet_atoms1=[[[ 0 for i in range(50000)] for j in xrange(0,MAX_ATOM_TYPES+1)] for k in xrange(0,MAX_ATOM_TYPES+1)]
    triplet_atoms2=[[[ 0 for i in range(50000)] for j in xrange(0,MAX_ATOM_TYPES+1)] for k in xrange(0,MAX_ATOM_TYPES+1)]
    triplet_atoms3=[[[ 0 for i in range(50000)] for j in xrange(0,MAX_ATOM_TYPES+1)] for k in xrange(0,MAX_ATOM_TYPES+1)]
    output_triplets=open(BASE_FILE+str("_triplets.data"), 'w')
    for i in range(0,TotalAtoms):
      neighbour_atoms=[]
      neighbour_atoms_type=[]
      if atom_type[i] == 1:
         for j in nnl[i]:
           if j != -1:
             neighbour_atoms.append(j)
             neighbour_atoms_type.append(atom_type[j])
         combs=list(it.combinations(neighbour_atoms,2))
         for pair in combs:
           a1=pair[0]
           a2=pair[1]
           b1=atom_type[a1]
           b2=atom_type[a2]
           if b1 < b2:
             b3=b1
             b4=b2
             a3=a1
             a4=a2
           elif b1 > b2:
             b3=b2
             b4=b1
             a3=a2
             a4=a1
           else :
             b3=b1
             b4=b2
             a3=a1
             a4=a2
           triplet_atoms1[b3][b4][triplet[b3][b4]]=a3
           triplet_atoms2[b3][b4][triplet[b3][b4]]=a4
           triplet_atoms3[b3][b4][triplet[b3][b4]]=i
           triplet[b3][b4]=triplet[b3][b4]+1

    for i in range(1,MAX_ATOM_TYPES+1):
      for j in range(i,MAX_ATOM_TYPES+1):
        #if i != j:
        output_triplets.write("%s-O-%s  %d \n" %(atype[i],atype[j],triplet[i][j]))
        if triplet[i][j] != 0:
          output_triplet_atoms_id=open(BASE_FILE+str("_triplet_atoms_")+str(atype[i])+str("O")+str(atype[j])+str(".atomsid"),'w')
          output_triplet_atoms_id.write("# triplet atom ids \n")
          output_triplet_atoms=open(BASE_FILE+str("_triplet_atoms_")+str(atype[i])+str("O")+str(atype[j])+str(".atoms"),'w')
          output_triplet_oxygen_atoms=open(BASE_FILE+str("_triplet_atoms_")+str(atype[i])+str("O")+str(atype[j])+str("_oxygen.atoms"),'w')
          write_imd_header1(output_triplet_atoms,imd_box,str(atype[i]),rc[i][1],str(atype[j]),rc[j][1])
          write_imd_header1(output_triplet_oxygen_atoms,imd_box,str(atype[i]),rc[i][1],str(atype[j]),rc[j][1])
          for k in range(triplet[i][j]):
            a1=triplet_atoms1[i][j][k]
            a2=triplet_atoms2[i][j][k]
            a3=triplet_atoms3[i][j][k]
            output_triplet_atoms_id.write("%d %d %d\n" %(a3,a1,a2))
            output_triplet_atoms.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[a1],atom_type[a1],atom_posx[a1],atom_posy[a1],atom_posz[a1],nnl_count[a1]))
            output_triplet_atoms.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[a2],atom_type[a2],atom_posx[a2],atom_posy[a2],atom_posz[a2],nnl_count[a2]))
            output_triplet_oxygen_atoms.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[a3],atom_type[a3],atom_posx[a3],atom_posy[a3],atom_posz[a3],nnl_count[a3]))
            a4=set(nnl[a1]) | set(nnl[a2])
            #a4.remove(-1)
            for l in a4:
              output_triplet_atoms.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[l],atom_type[l],atom_posx[l],atom_posy[l],atom_posz[l],nnl_count[l]))
          output_triplet_atoms.close()
          output_triplet_oxygen_atoms.close()
          output_triplet_atoms_id.close()

    output_triplet_oxygen_in_SiOSi=open(BASE_FILE+str("_triplet_atoms_oxygen_in_Si4OSi4.atoms"),'w')
    write_imd_header3(output_triplet_oxygen_in_SiOSi,imd_box,atype,rc)
    # **************************************
    # triplet of Si-O-Si which are perfectly tetrahedra
    Si_O_Si=0
    for i in range(0,TotalAtoms):
      neighbour_atoms=[]
      neighbour_atoms_type=[]
      if atom_type[i] ==1:
         for j in nnl[i]:
           if j != -1:
             neighbour_atoms.append(j)
             neighbour_atoms_type.append(atom_type[j])
         if list(neighbour_atoms_type).count(2) == 2:
           output_triplet_oxygen_in_SiOSi.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
           for j in nnl[i]:
             if j != -1:
               if atom_type[j] == 2 and nnl_count[j] == 4:
                 Si_O_Si=Si_O_Si+1
    output_triplets.write("\n\n# Si-O-Si count of perfect tetrahedras\n")
    output_triplets.write("Si-O-Si = %d \n" %(Si_O_Si/2))

    output_triplet_oxygen_in_SiOAl=open(BASE_FILE+str("_triplet_atoms_oxygen_in_Si4OAl4.atoms"),'w')
    write_imd_header3(output_triplet_oxygen_in_SiOAl,imd_box,atype,rc)
    # **************************************
    # triplet of Si-O-Al which are perfectly tetrahedra
    Si_O_Al=0
    Select_Al_Atoms_In_Triplets=[0 for i in xrange(TotalAtoms)]
    for i in range(0,TotalAtoms):
      neighbour_atoms=[]
      neighbour_atoms_type=[]
      if atom_type[i] ==1:
         for j in nnl[i]:
           if j != -1:
             neighbour_atoms.append(j)
             neighbour_atoms_type.append(atom_type[j])
         if list(neighbour_atoms_type).count(2) == 1 and list(neighbour_atoms_type).count(3) == 1:
           output_triplet_oxygen_in_SiOAl.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
           for j in nnl[i]:
             if j != -1:
               if (atom_type[j] == 2 or atom_type[j] == 3) and nnl_count[j] == 4:
                 Si_O_Al=Si_O_Al+1
                 if atom_type[j] == 3:
                   Select_Al_Atoms_In_Triplets[j]=1
    Al_Atoms_In_Triplets1=0
    for i in range(0,TotalAtoms):
      if atom_type[i] == 3:
        if Select_Al_Atoms_In_Triplets[i] == 1:
          Al_Atoms_In_Triplets1=Al_Atoms_In_Triplets1+1
 
    output_triplets.write("\n\n# Si-O-Al count of perfect tetrahedras\n")
    output_triplets.write("Si-O-Al = %d \n" %(Si_O_Al/2))

    output_triplet_oxygen_in_AlOAl=open(BASE_FILE+str("_triplet_atoms_oxygen_in_Al4OAl4.atoms"),'w')
    write_imd_header3(output_triplet_oxygen_in_AlOAl,imd_box,atype,rc)
    # **************************************
    # triplet of Al-O-Al which are perfectly tetrahedra
    Al_O_Al=0
    for i in range(0,TotalAtoms):
      neighbour_atoms=[]
      neighbour_atoms_type=[]
      if atom_type[i] ==1:
         for j in nnl[i]:
           if j != -1:
             neighbour_atoms.append(j)
             neighbour_atoms_type.append(atom_type[j])
         if list(neighbour_atoms_type).count(3) == 2:
           output_triplet_oxygen_in_AlOAl.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
           for j in nnl[i]:
             if j != -1:
               if atom_type[j] == 3 and nnl_count[j] == 4:
                 Al_O_Al=Al_O_Al+1
                 if atom_type[j] == 3:
                   Select_Al_Atoms_In_Triplets[j]=1
    Al_Atoms_In_Triplets2=0
    for i in range(0,TotalAtoms):
      if atom_type[i] == 3:
        if Select_Al_Atoms_In_Triplets[i] == 1:
          Al_Atoms_In_Triplets2=Al_Atoms_In_Triplets2+1

    output_triplets.write("\n\n# Al-O-Al count of perfect tetrahedras\n")
    output_triplets.write("Al-O-Al = %d \n" %(Al_O_Al/2))
    output_triplets.write("\n\n# total Al atoms which are in Si-O-Al = %d" %(Al_Atoms_In_Triplets1))
    output_triplets.write("\n\n# total Al atoms which are in Al-O-Al = %d" %(Al_Atoms_In_Triplets2))
    output_triplet_oxygen_in_SiOSi.close()
    output_triplet_oxygen_in_SiOAl.close()
    output_triplet_oxygen_in_AlOAl.close()
    # **************************************************************************************
    # finding Q0, Q1, Q2, Q3, Q4
    Q=[[0 for i in xrange(5)] for j in xrange(5)]
    Q_Si_atoms=[[[0 for k in xrange(50000)] for j in xrange(5)] for i in xrange(5)]
    Q_Si_status=[-1 for i in xrange(TotalAtoms)]
    Q_Al_status=[-1 for i in xrange(TotalAtoms)]
    SelectAtoms1=[-1 for i in xrange(TotalAtoms)]
    SelectAtoms2=[-1 for i in xrange(TotalAtoms)]
    Q_Al_list=[0 for i in xrange(TotalAtoms)]
    output_q=open(BASE_FILE+str("_QSi_info.data"), 'w')
    output_q_temp1=open(BASE_FILE+str("_QSi_additional_info.data"), 'w')
    output_q.write("# QSi information\n")
    output_q_temp1.write("#Special Si which has 4 neighboring oxygen atoms but the neighbouring oxygen atoms has bonded with more than 4 cations\n")
    # looping all atoms
    for Si1 in xrange(TotalAtoms):
      # go only Si atoms with 4 coordination
      if atom_type[Si1] == 2 and nnl_count[Si1] == 4:
        count_Si=0
        count_Al=0
        tri_clusters_status=0
        # first neighbours of Si (i.e O1)
        for O1 in nnl[Si1]:
          if tri_clusters_status == 0:
            if (nnl_type[O1].count(2)+nnl_type[O1].count(3)) > 2:
              tri_clusters_status = 1
              #print "%d %d %d" %(O1,tri_clusters_status,(nnl_type[O1].count(2)+nnl_type[O1].count(3)))
        if tri_clusters_status == 0:
          for O1 in nnl[Si1]:
            # second neighbours of Si (i.e cations)
            for Si2 in nnl[O1]:
              # count if Si2 is not equal to Si1 also Si2 is 4 coordinated
              if Si2 != Si1:
                # increase the count if the other cation is Si
                if atom_type[Si2] == 2 and nnl_count[Si2] == 4:
                  count_Si=count_Si+1
                # increase the count if the other cation is Al
                elif atom_type[Si2] == 3 and nnl_count[Si2] == 4:
                  count_Al=count_Al+1
                  Q_Al_list[Si2]=1
          if count_Si+count_Al < 5:
            Q_Si_atoms[count_Si+count_Al][count_Al][Q[count_Si+count_Al][count_Al]]=Si1
            Q_Si_status[Si1]=count_Si+count_Al
            Q_Al_status[Si1]=count_Al
            Q[count_Si+count_Al][count_Al] = Q[count_Si+count_Al][count_Al]+1
        else:
          output_q_temp1.write("%d \n" %(Si1))
 
    for i in range(5):
      for j in range(5):
        if j<=i:
          output_q.write("QSi%d (%dAl) = %d \n" %(i,j,Q[i][j]))
          if j == 0:
            for k in Q_Si_atoms[i][j]:
              SelectAtoms1[k]=1
              for l in range(nnl_count[k]):
                SelectAtoms1[nnl[k][l]]=1
          else:
            for k in Q_Si_atoms[i][j]:
              SelectAtoms2[k]=1
              for l in range(nnl_count[k]):
                SelectAtoms2[nnl[k][l]]=1
  
    output_q_atoms1=open(BASE_FILE+str("_QnAl0_Si_tets.atoms"),'w')
    output_q_atoms2=open(BASE_FILE+str("_QnAlm_Si_tets.atoms"),'w')
    write_imd_header2(output_q_atoms1,imd_box,str(atype[2]),rc[2][1],str(atype[3]),rc[3][1])
    write_imd_header2(output_q_atoms2,imd_box,str(atype[2]),rc[2][1],str(atype[3]),rc[3][1])
    for k in range(TotalAtoms):
      if SelectAtoms1[k] == 1:
        output_q_atoms1.write("%d  %d  %lf  %lf  %lf %d %d %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k],Q_Si_status[k],Q_Al_status[k]))
    for k in range(TotalAtoms):
      if SelectAtoms2[k] == 1:
        output_q_atoms2.write("%d  %d  %lf  %lf  %lf %d %d %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k],Q_Si_status[k],Q_Al_status[k]))

    output_q_Al_list=open(BASE_FILE+str("_QnAlm_Al_tets.atoms"),'w')
    write_imd_header2(output_q_Al_list,imd_box,str(atype[2]),rc[2][1],str(atype[3]),rc[3][1])
    for i in range(TotalAtoms):
      if Q_Al_list[i]==1:
        output_q_Al_list.write("%d  %d  %lf  %lf  %lf %d %d %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i],Q_Si_status[i],Q_Al_status[i]))
        for j in range(nnl_count[i]):
          output_q_Al_list.write("%d  %d  %lf  %lf  %lf %d %d %d \n" %(atom_global_id[nnl[i][j]],atom_type[nnl[i][j]],atom_posx[nnl[i][j]],atom_posy[nnl[i][j]],atom_posz[nnl[i][j]],nnl_count[nnl[i][j]],Q_Si_status[nnl[i][j]],Q_Al_status[nnl[i][j]]))

    output_q.close()
    output_q_temp1.close()
    output_q_atoms1.close()
    output_q_atoms2.close()
    output_q_Al_list.close()

    # **************************************************************************************
    # finding Q0, Q1, Q2, Q3, Q4 based on Al
    # you should definately know what you are doing otherwise it is very difficult to understand this code
    # Since I dont have a lot of time, I am not cleaning the code, cleaning of code is needed for better understanding
    output1=open(BASE_FILE+str("_QnSim_SomeSelected_AlOAl_tets.atoms"),'w')
    output2=open(BASE_FILE+str("_QnSim_SomeSelected_AlOAl_Neighbortets.atoms"),'w')
    write_imd_header3(output1,imd_box,atype,rc)
    write_imd_header3(output2,imd_box,atype,rc)
    Q=[[0 for i in xrange(5)] for j in xrange(5)]
    Q_Al_atoms=[[[0 for k in xrange(50000)] for j in xrange(5)] for i in xrange(5)]
    Q_SiAl_status=[-1 for i in xrange(TotalAtoms)]
    Q_Si_status=[-1 for i in xrange(TotalAtoms)]
    #SelectAtoms1=[-1 for i in xrange(TotalAtoms)]
    #SelectAtoms2=[-1 for i in xrange(TotalAtoms)]
    #Q_Si_list=[0 for i in xrange(TotalAtoms)]
    output_q=open(BASE_FILE+str("_QAl_info.data"), 'w')
    output_q_temp1=open(BASE_FILE+str("_QAl_additional_info.data"), 'w')
    output_q.write("# QAl information\n")
    output_q_temp1.write("#Special Al which has 4 neighboring oxygen atoms but the neighbouring oxygen atoms has bonded with more than 4 cations\n")
    # looping all atoms
    for Al1 in range(TotalAtoms):
      # go only Al atoms with 4 coordination
      if atom_type[Al1] == 3 and nnl_count[Al1] == 4:
        count_Si=0
        count_Al=0
        tri_clusters_status=0
        # first neighbours of Si (i.e O)
        for O1 in nnl[Al1]:
          if tri_clusters_status == 0:
            if (nnl_type[O1].count(2)+nnl_type[O1].count(3)) > 2:
              tri_clusters_status = 1
              #print "%d %d %d" %(O1,tri_clusters_status,(nnl_type[O1].count(2)+nnl_type[O1].count(3)))
        if tri_clusters_status == 0:
          for O1 in nnl[Al1]:
            # second neighbours of Si (i.e cations)
            for Al2 in nnl[O1]:
              # go only to other cations and also which is 4 coordinated
              if Al2 != Al1:
                # increase the count if the other cation is Si
                if atom_type[Al2] == 2 and nnl_count[Al2] == 4:
                  count_Si=count_Si+1
                  #Q_Si_list[nnl[O1][k]]=1
                # increase the count if the other cation is Al
                # elif atom_type[nnl[O1][k]] == 3 and nnl_count[nnl[O1][k]] == 4:
                elif atom_type[Al2] == 3 and nnl_count[Al2] == 4:
                  count_Al=count_Al+1
                  #Q_Al_list[nnl[O1][k]]=1
          if count_Si+count_Al < 5:
            #Q_Al_atoms[count_Al+count_Si][count_Si][Q[count_Al+count_Si][count_Si]]=i
            Q_SiAl_status[Al1]=count_Si+count_Al
            Q_Si_status[Al1]=count_Si
            Q[count_Si+count_Al][count_Si] = Q[count_Si+count_Al][count_Si]+1
            if (count_Si+count_Al == 3 and count_Si == 1) or (count_Si+count_Al == 3 and count_Si == 2) or \
               (count_Si+count_Al == 4 and count_Si == 1) or (count_Si+count_Al == 4 and count_Si == 2) or \
               (count_Si+count_Al == 4 and count_Si == 3) :
               output1.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[Al1],atom_type[Al1],atom_posx[Al1],atom_posy[Al1],atom_posz[Al1],nnl_count[Al1]))
               output2.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[Al1],atom_type[Al1],atom_posx[Al1],atom_posy[Al1],atom_posz[Al1],nnl_count[Al1]))
               for j in nnl[Al1]:
                 output1.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
                 #output2.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
                 for k in nnl[j]:
                   if i != k:
                     output2.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
                     for l in nnl[k]:
                       output2.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[l],atom_type[l],atom_posx[l],atom_posy[l],atom_posz[l],nnl_count[l]))
        else:
          output_q_temp1.write("%d \n" %(Al1))
  
    for i in range(5):
      for j in range(5):
        if j<=i:
          output_q.write("QAl%d (%dSi) = %d \n" %(i,j,Q[i][j]))
          #if j == 0:
          #  for k in Q_Si_atoms[i][j]:
          #    SelectAtoms1[k]=1
          #    for l in range(nnl_count[k]):
          #      SelectAtoms1[nnl[k][l]]=1
          #else:
          #  for k in Q_Si_atoms[i][j]:
          #    SelectAtoms2[k]=1
          #    for l in range(nnl_count[k]):
          #      SelectAtoms2[nnl[k][l]]=1

    #output_q_atoms1=open(BASE_FILE+str("_QnAl0_Si_tets.atoms"),'w')
    #output_q_atoms2=open(BASE_FILE+str("_QnAlm_Si_tets.atoms"),'w')
    #write_imd_header2(output_q_atoms1,imd_box,str(atype[2]),rc[2][1],str(atype[3]),rc[3][1])
    #write_imd_header2(output_q_atoms2,imd_box,str(atype[2]),rc[2][1],str(atype[3]),rc[3][1])
    #for k in range(TotalAtoms):
    #  if SelectAtoms1[k] == 1:
    #    output_q_atoms1.write("%d  %d  %lf  %lf  %lf %d %d %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k],Q_Si_status[k],Q_Al_status[k]))
    #for k in range(TotalAtoms):
    #  if SelectAtoms2[k] == 1:
    #    output_q_atoms2.write("%d  %d  %lf  %lf  %lf %d %d %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k],Q_Si_status[k],Q_Al_status[k]))

    #output_q_Al_list=open(BASE_FILE+str("_QnAlm_Al_tets.atoms"),'w')
    #write_imd_header2(output_q_Al_list,imd_box,str(atype[2]),rc[2][1],str(atype[3]),rc[3][1])
    #for i in range(TotalAtoms):
    #  if Q_Al_list[i]==1:
    #    output_q_Al_list.write("%d  %d  %lf  %lf  %lf %d %d %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i],Q_Si_status[i],Q_Al_status[i]))
    #    for j in range(nnl_count[i]):
    #      output_q_Al_list.write("%d  %d  %lf  %lf  %lf %d %d %d \n" %(atom_global_id[nnl[i][j]],atom_type[nnl[i][j]],atom_posx[nnl[i][j]],atom_posy[nnl[i][j]],atom_posz[nnl[i][j]],nnl_count[nnl[i][j]],Q_Si_status[nnl[i][j]],Q_Al_status[nnl[i][j]]))

    output_q.close()
    output_q_temp1.close()
    #output_q_atoms1.close()
    #output_q_atoms2.close()
    #output_q_Al_list.close()
    output1.close()
    output2.close()

    # **************************************************************************************
    # compute bond lengths based on Si-O-Si, Si-O-Al, Si-O-Ca, Si-O-Mg environments
    output11=open(BASE_FILE+str("_triplets_average_bond_lengths.data"),'w')
    sum_bond_length=[[[[[ 0 for e in xrange(MAX_ATOM_TYPES+1)] for a in xrange(2)] for b in xrange(2)] for c in xrange(2)] for d in xrange(2)]
    cnt_bond_length=[[[[[ 1 for e in xrange(MAX_ATOM_TYPES+1)] for a in xrange(2)] for b in xrange(2)] for c in xrange(2)] for d in xrange(2)]

    for i in range(TotalAtoms):
      if atom_type[i] == 1:
        temp_Si = 0
        temp_Al = 0
        temp_Ca = 0
        temp_Mg = 0
        if 2 in nnl_type[i]:
          temp_Si = 1
        if 3 in nnl_type[i]:
          temp_Al = 1
        if 4 in nnl_type[i]:
          temp_Ca = 1
        if 5 in nnl_type[i]:
          temp_Mg = 1
        for j in range(nnl_count[i]):
          sum_bond_length[temp_Si][temp_Al][temp_Ca][temp_Mg][atom_type[nnl[i][j]]] = sum_bond_length[temp_Si][temp_Al][temp_Ca][temp_Mg][atom_type[nnl[i][j]]] + nnl_each_pair_distance[i][j]
          cnt_bond_length[temp_Si][temp_Al][temp_Ca][temp_Mg][atom_type[nnl[i][j]]] = cnt_bond_length[temp_Si][temp_Al][temp_Ca][temp_Mg][atom_type[nnl[i][j]]] + 1

    for i in range(2):
      for j in range(2):
        for k in range(2):
          for l in range(2):
            temp_Si=""
            temp_Al=""
            temp_Ca=""
            temp_Mg=""
            if i == 1:
              temp_Si ="Si"
            if j == 1:
              temp_Al ="Al"
            if k == 1:
              temp_Ca ="Ca"
            if l == 1:
              temp_Mg ="Mg"
            #output11=open(BASE_FILE+str(_)+str()+str()+str()+str()+str("_bond_length.data"))
            output11.write("\n\n# %s%s%s%s \n" %(temp_Si,temp_Al,temp_Ca,temp_Mg))
            for m in range(2,MAX_ATOM_TYPES+1):
              if cnt_bond_length[i][j][k][l][m] > 1:
                output11.write("average bond length of %s-O = %lf , for total bonds of %d \n" %(atype[m],sum_bond_length[i][j][k][l][m]/cnt_bond_length[i][j][k][l][m],cnt_bond_length[i][j][k][l][m]))
    output11.close()

    # **************************************************************************************
    # coordination of each atom
    Coord=[[0 for j in range(MAX_NEIGHBOURS+1)] for i in range(MAX_ATOM_TYPES+1)]		# Coord[atom_type][coordination]
    for i in range(TotalAtoms):
      Coord[atom_type[i]][nnl_count[i]]=Coord[atom_type[i]][nnl_count[i]]+1
    output_coord1=open(BASE_FILE+str("_average_coord"),'w')
    output_coord1.write("# Chemical_Element\tNumber_of_Atoms\tAverageCoordination\n")
 
    MAX_LENGTH_OF_COORD=0
    for i in range(1,MAX_ATOM_TYPES+1):
      avg_coord1=0
      avg_coord2=0
      output_coord=open(BASE_FILE+str("_")+str(atype[i])+str("_coord"),'w')
      for j in range(0,MAX_NEIGHBOURS+1):
        output_coord.write("%d %d \n" %(j,Coord[i][j]))
        MAX_LENGTH_OF_COORD=max(MAX_LENGTH_OF_COORD,Coord[i][j])
        avg_coord1=avg_coord1+(j*Coord[i][j])
        avg_coord2=avg_coord2+Coord[i][j]
      output_coord.close()
      if avg_coord2 != 0.0:
        output_coord1.write("%s  %d  %lf \n" %(atype[i],total_atoms_of_type[i],1.0*avg_coord1/avg_coord2))
    output_coord1.close()

    Coordination=[[[0 for k in range(MAX_LENGTH_OF_COORD)] for j in range(MAX_NEIGHBOURS+1)] for i in range(MAX_ATOM_TYPES+1)] 	# Coordination[atom_type][coordination][count] = atom_id
    Coordination_count=[[0 for j in range(MAX_NEIGHBOURS+1)] for i in range(MAX_ATOM_TYPES+1)]					# Coordination_count[atom_type][coordination]
    for i in range(TotalAtoms):
      Coordination[atom_type[i]][nnl_count[i]][Coordination_count[atom_type[i]][nnl_count[i]]]=i
      Coordination_count[atom_type[i]][nnl_count[i]]=Coordination_count[atom_type[i]][nnl_count[i]]+1

    for i in range(1,MAX_ATOM_TYPES+1):
      for j in range(MAX_NEIGHBOURS+1):
        if Coord[i][j] !=0:
          output=open(BASE_FILE+str("_")+str(atype[i])+str("_")+str(j)+str("_coord.atomsid"),'w')
          for k in range(Coordination_count[i][j]):
            l=Coordination[i][j][k]
            output.write("%d  " %(atom_global_id[l]))
            for m in nnl[l]:
              if m != -1:
                output.write("%d  " %(atom_global_id[m]))
            output.write("\n")
          output.close


    # **************************************************************************************
    # tri-cluster (1 oxygen with 3 Al)
    SelectAtoms=[0 for i in xrange(TotalAtoms)]
    output2=open(BASE_FILE+str("_tri-cluster_3Al1O.atoms"),'w')
    output1=open(BASE_FILE+str("_tri-cluster_OnlyOxygen.atoms"),'w')
    write_imd_header3(output1,imd_box,atype,rc)
    write_imd_header3(output2,imd_box,atype,rc)
    for i in range(TotalAtoms):
      if atom_type[i] == 1:
        temp_atom_type=[]
        for j in nnl[i]:
          temp_atom_type.append(atom_type[j])
        if list(temp_atom_type).count(3) == 3 :
          output1.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
          for j in nnl[i]:
            if atom_type[j] == 3:
              SelectAtoms[j]=1
              for k in nnl[j]:
                SelectAtoms[k]=1
    for i in range(TotalAtoms):
      if SelectAtoms[i]==1:
        output2.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
    output1.close()
    output2.close()
  
    # **************************************************************************************
  #  # compute %Al in Al(IV)-O-Al, Al(IV)-O-Al
  #  for i in range(TotalAtoms):						# iterating all atoms for Al(1st)
  #    if ((atom_type[i] == 3) and (nnl_count[i] == 4)):			# now i = Al(1st)
  #      for j in nnl[i]:							# iterating Al(1st) neighbours for O(1st)
  #        if atom_type[j] == 1:						# now j = O(1st)
  #          for k in nnl[j]:						# iterating O(1st) neighbours for Al(2nd)
  #            if i != k and atom_type[k] == 3 and nnl_count[k] == 4:	# 
  
    # **************************************************************************************
    # Study based on Al-O-Al 
  #  output1=open(BASE_FILE+str("_AlOAl_additional_info.data"),'w')
  #  write_imd_header3(output1,imd_box,atype,rc)
  #  for i in range(TotalAtoms):
  #    if atom_type[i] == 1:
  #      neighbour_atom_types=[]
  #      for j in nnl[i]:
  #        neighbour_atom_types.append(atom_type[j])
  #      number_of_Si=neighbour_atom_types.count(2)
  #      number_of_Al=neighbour_atom_types.count(3)
  #      number_of_Ca=neighbour_atom_types.count(4)
  #      number_of_Mg=neighbour_atom_types.count(5)
  #      if number_of_Al == 2 and number_of_Si == 0 :
        
    # ***************************************************************************************
    # Getting out Ca and Mg neighbours information based on NBO or BO also seprate if it is NBO whether it is Si or Al and if it is BO then whether it is 2Si or SiAl or 2Al
    output1=open(BASE_FILE+str("_Ca_neighbours_info.data"),'w')
    output2=open(BASE_FILE+str("_Mg_neighbours_info.data"),'w')
    output1.write("# Ca neighbours information\n")
    output1.write("# NBO_atoms \t BO_atoms \t count\n")
    output2.write("# Mg neighbours information\n")
    output2.write("# NBO_atoms \t BO_atoms \t count\n")
    Ca_NBRS=[[0 for i in xrange(10)] for j in xrange(10)]
    Mg_NBRS=[[0 for i in xrange(10)] for j in xrange(10)]
    Ca_NBRS1=[[[[[0 for i1 in xrange(10)] for i2 in xrange(10)] for i3 in xrange(10)] for i4 in xrange(10)] for i5 in xrange(10)]
    Mg_NBRS1=[[[[[0 for i1 in xrange(10)] for i2 in xrange(10)] for i3 in xrange(10)] for i4 in xrange(10)] for i5 in xrange(10)]
    Ca_TOTAL_NBO=[0 for i in xrange(TotalAtoms)]
    Ca_TOTAL_NBO_Si=[0 for i in xrange(TotalAtoms)]
    Ca_TOTAL_NBO_Al=[0 for i in xrange(TotalAtoms)]
    Ca_TOTAL_BO=[0 for i in xrange(TotalAtoms)]
    Ca_TOTAL_BO_2Si=[0 for i in xrange(TotalAtoms)]
    Ca_TOTAL_BO_SiAl=[0 for i in xrange(TotalAtoms)]
    Ca_TOTAL_BO_2Al=[0 for i in xrange(TotalAtoms)]
    Mg_TOTAL_NBO=[0 for i in xrange(TotalAtoms)]
    Mg_TOTAL_NBO_Si=[0 for i in xrange(TotalAtoms)]
    Mg_TOTAL_NBO_Al=[0 for i in xrange(TotalAtoms)]
    Mg_TOTAL_BO=[0 for i in xrange(TotalAtoms)]
    Mg_TOTAL_BO_2Si=[0 for i in xrange(TotalAtoms)]
    Mg_TOTAL_BO_SiAl=[0 for i in xrange(TotalAtoms)]
    Mg_TOTAL_BO_2Al=[0 for i in xrange(TotalAtoms)]

    for i in range(TotalAtoms):
      if atom_type[i] == 4 :
        Ca_NBO_count=0
        Ca_NBO_Si_count=0
        Ca_NBO_Al_count=0
        Ca_BO_count=0
        Ca_BO_2Si_count=0
        Ca_BO_SiAl_count=0
        Ca_BO_2Al_count=0
        for j in nnl[i]:
          if NBOS[j] == 1:
            Ca_NBO_count=Ca_NBO_count+1
            if nnl_type[j].count(2) == 1 and nnl_type[j].count(3) == 0:
              Ca_NBO_Si_count=Ca_NBO_Si_count+1
            elif nnl_type[j].count(2) == 0 and nnl_type[j].count(3) == 1:
              Ca_NBO_Al_count=Ca_NBO_Al_count+1
          elif NBOS[j] == 0:
            Ca_BO_count=Ca_BO_count+1
            if nnl_type[j].count(2) == 2:
              Ca_BO_2Si_count=Ca_BO_2Si_count+1
            elif nnl_type[j].count(2) == 1 and nnl_type[j].count(3) == 1:
              Ca_BO_SiAl_count=Ca_BO_SiAl_count+1
            elif nnl_type[j].count(3) == 2:
              Ca_BO_2Al_count=Ca_BO_2Al_count+1
        Ca_NBRS[Ca_NBO_count][Ca_BO_count]=Ca_NBRS[Ca_NBO_count][Ca_BO_count]+1
        Ca_TOTAL_NBO[i]     = Ca_NBO_count
        Ca_TOTAL_NBO_Si[i]  = Ca_NBO_Si_count
        Ca_TOTAL_NBO_Al[i]  = Ca_NBO_Al_count
        Ca_TOTAL_BO[i]      = Ca_BO_count
        Ca_TOTAL_BO_2Si[i]  = Ca_BO_2Si_count
        Ca_TOTAL_BO_SiAl[i] = Ca_BO_SiAl_count
        Ca_TOTAL_BO_2Al[i]  = Ca_BO_2Al_count
        Ca_NBRS1[Ca_NBO_Si_count][Ca_NBO_Al_count][Ca_BO_2Si_count][Ca_BO_SiAl_count][Ca_BO_2Al_count]=Ca_NBRS1[Ca_NBO_Si_count][Ca_NBO_Al_count][Ca_BO_2Si_count][Ca_BO_SiAl_count][Ca_BO_2Al_count]+1
  
    # writing out Ca neighbours
    for i in range(8):
      for j in range(8):
        if Ca_NBRS[i][j] != 0:
          output1.write("%d \t %d \t %d\n" %(i,j,Ca_NBRS[i][j]))
    output1.write("\n\n# NBO_atoms \t BO_atoms \t count\n")
    output1.write("# NBO_Si \t NBO_Al \t BO_2Si \t BO_SiAl \t BO_2Al \t count\n")
    for k in range(10):
      output1.write("### total_neighbours = %d \n" %(k))
      for i in range(8):
        for j in range(8):
          if i+j == k:
            if Ca_NBRS[i][j] != 0:
              output1.write("## NBO's and BO's %d \t %d \t %d\n" %(i,j,Ca_NBRS[i][j]))
              if Ca_NBRS[i][j] > 100:
                for i1 in range(8):
                  for i2 in range(8):
                    for i3 in range(8):
                      for i4 in range(8):
                        for i5 in range(8):
                          if Ca_NBRS1[i1][i2][i3][i4][i5] != 0 and i1+i2 == i and i3+i4+i5 == j:
                            output1.write("%d  %d  %d  %d  %d %d \n" %(i1,i2,i3,i4,i5,Ca_NBRS1[i1][i2][i3][i4][i5]))

    # list Mg neighbours
    for i in range(TotalAtoms):
      if atom_type[i] == 5 :
        Mg_NBO_count=0
        Mg_NBO_Si_count=0
        Mg_NBO_Al_count=0
        Mg_BO_count=0
        Mg_BO_2Si_count=0
        Mg_BO_SiAl_count=0
        Mg_BO_2Al_count=0
        for j in nnl[i]:
          if NBOS[j] == 1:
            Mg_NBO_count=Mg_NBO_count+1
            if nnl_type[j].count(2) == 1 and nnl_type[j].count(3) == 0:
              Mg_NBO_Si_count=Mg_NBO_Si_count+1
            elif nnl_type[j].count(2) == 0 and nnl_type[j].count(3) == 1:
              Mg_NBO_Al_count=Mg_NBO_Al_count+1
          elif NBOS[j] == 0:
            Mg_BO_count=Mg_BO_count+1
            if nnl_type[j].count(2) == 2:
              Mg_BO_2Si_count=Mg_BO_2Si_count+1
            elif nnl_type[j].count(2) == 1 and nnl_type[j].count(3) == 1:
              Mg_BO_SiAl_count=Mg_BO_SiAl_count+1
            elif nnl_type[j].count(3) == 2:
              Mg_BO_2Al_count=Mg_BO_2Al_count+1
        Mg_NBRS[Mg_NBO_count][Mg_BO_count]=Mg_NBRS[Mg_NBO_count][Mg_BO_count]+1
        Mg_TOTAL_NBO[i]     = Mg_NBO_count
        Mg_TOTAL_NBO_Si[i]  = Mg_NBO_Si_count
        Mg_TOTAL_NBO_Al[i]  = Mg_NBO_Al_count
        Mg_TOTAL_BO[i]      = Mg_BO_count
        Mg_TOTAL_BO_2Si[i]  = Mg_BO_2Si_count
        Mg_TOTAL_BO_SiAl[i] = Mg_BO_SiAl_count
        Mg_TOTAL_BO_2Al[i]  = Mg_BO_2Al_count
        Mg_NBRS1[Mg_NBO_Si_count][Mg_NBO_Al_count][Mg_BO_2Si_count][Mg_BO_SiAl_count][Mg_BO_2Al_count]=Mg_NBRS1[Mg_NBO_Si_count][Mg_NBO_Al_count][Mg_BO_2Si_count][Mg_BO_SiAl_count][Mg_BO_2Al_count]+1
    # writing out Mg neighbours
    for i in range(8):
      for j in range(8):
        if Mg_NBRS[i][j] != 0:
          output2.write("%d \t %d \t %d\n" %(i,j,Mg_NBRS[i][j]))
    output2.write("\n\n# NBO_atoms \t BO_atoms \t count\n")
    output2.write("# NBO_Si \t NBO_Al \t BO_2Si \t BO_SiAl \t BO_2Al \t count\n")
    for k in range(10):
      output2.write("### total_neighbours = %d \n" %(k))
      for i in range(8):
        for j in range(8):
          if i+j == k :
            if Mg_NBRS[i][j] != 0:
              output2.write("## NBO's and BO's  %d \t %d \t %d\n" %(i,j,Mg_NBRS[i][j]))
              if Mg_NBRS[i][j] > 40:
                for i1 in range(8):
                  for i2 in range(8):
                    for i3 in range(8):
                      for i4 in range(8):
                        for i5 in range(8):
                          if Mg_NBRS1[i1][i2][i3][i4][i5] != 0 and i1+i2 == i and i3+i4+i5 == j:
                            output2.write("%d  %d  %d  %d  %d %d \n" %(i1,i2,i3,i4,i5,Mg_NBRS1[i1][i2][i3][i4][i5]))
 
    output1.close()
    output2.close()
  
    # output all atoms in imd format to avoid the bounding box problem
    output1=open(BASE_FILE+str(".imd"),'w')
    write_imd_header3(output1,imd_box,atype,rc)
    for i in range(TotalAtoms):
      output1.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
    output1.close()

    #######################################################################################################
    # working on bridging and non-bridging oxygens again
    # May be I am rewriting the part of existing code again
    #atom_nnl_type=[[-1 for j in xrange(MAX_NEIGHBOURS)] for i in xrange(TotalAtoms)]
    atom_nnl_type_count=[[-1 for i in xrange(6)] for j in xrange(TotalAtoms)]	
    ca1=8
    count_BO=[[[[ 0 for i in xrange(ca1)] for j in xrange(ca1)] for k in xrange(ca1)] for l in range(ca1)]	# count_BO[Si][Al][Ca][Mg]
    count_BO1=[[[ 0 for i in xrange(2*ca1)] for j in xrange(ca1)] for k in xrange(ca1)]        			# count_BO[Si][Al][R]
    O_env_status=[0 for i in xrange(TotalAtoms)]									# Giving all the oxygen a status based on the environment
    # 1=SiAlR; 2=2SiR; 3=2AlR; 4=SiAl2R; 5=2Al2R; 6=Si2R; 7=Si3R; 8=Al3R; 9=Al2R; 10=2Al2R; 11=Al2R; 12=SiR
    SelectAtoms1_O=[0 for i in xrange(TotalAtoms)]
    count_O_connected_to_Si=0
    ca2=4
    ca3=10			# 1 to 10 or 0 to 9
    cb=[0 for i in xrange(ca3)]
    cc=[0 for i in xrange(ca3)]
    cd=[0 for i in xrange(ca3+3)]
    atom_nnl_O_env=[[-1 for i in xrange(MAX_NEIGHBOURS)] for j in xrange(TotalAtoms)]
    SiOSi_count1=SiOAl_count1=AlOAl_count1=0

    for i in xrange(TotalAtoms):
      count_O=count_Si=count_Al=count_Ca=count_Mg=0
      #for j in xrange(nnl_count[i]):
      #  atom_nnl_type[i][j]=atom_type[nnl[i][j]]
      for j in atom_nnl_type[i]:
        if j == 1:
          count_O  = count_O  + 1
        elif j == 2:
          count_Si = count_Si + 1
        elif j == 3:
          count_Al = count_Al + 1
        elif j == 4:
          count_Ca = count_Ca + 1
        elif j == 5:
          count_Mg = count_Mg + 1
      atom_nnl_type_count[i][1]=count_O
      atom_nnl_type_count[i][2]=count_Si
      atom_nnl_type_count[i][3]=count_Al
      atom_nnl_type_count[i][4]=count_Ca
      atom_nnl_type_count[i][5]=count_Mg

    for i in xrange(TotalAtoms):
      a2=a3=a4=a5=0
      if atom_type[i] == 1:
        a2=atom_nnl_type_count[i][2]
        a3=atom_nnl_type_count[i][3]
        a4=atom_nnl_type_count[i][4]
        a5=atom_nnl_type_count[i][5]
        count_BO[a2][a3][a4][a5] = count_BO[a2][a3][a4][a5]+1
        count_BO1[a2][a3][a4+a5] = count_BO1[a2][a3][a4+a5]+1
        # 1=SiAlR; 2=2SiR; 3=2AlR; 4=SiAl2R; 5=2Al2R; 6=Si2R; 7=Si3R; 8=Al3R; 9=Al2R;
        if a2 == 1 and a3 == 1 and (a4+a5) == 1:
          O_env_status[i]=1
        elif a2 == 2 and a3 == 0 and (a4+a5) == 1:
          O_env_status[i]=2
        elif a2 == 0 and a3 == 2 and (a4+a5) == 1:
          O_env_status[i]=3
        elif a2 == 1 and a3 == 1 and (a4+a5) == 2:
          O_env_status[i]=4
        elif a2 == 0 and a3 == 2 and (a4+a5) == 2:
          O_env_status[i]=5
        elif a2 == 1 and a3 == 0 and (a4+a5) == 2:
          O_env_status[i]=6
        elif a2 == 1 and a3 == 0 and (a4+a5) == 3:
          O_env_status[i]=7
        elif a2 == 0 and a3 == 1 and (a4+a5) == 3:
          O_env_status[i]=8
        elif a2 == 0 and a3 == 1 and (a4+a5) == 2:
          O_env_status[i]=9
        #elif a2 == 1 and a3 == 0 and (a4+a5) == 1:
        #  O_env_status[i]=12
        elif a2 == 2 and a3 == 0 and (a4+a5) == 0:
          O_env_status[i]=11
        elif a2 == 1 and a3 == 1 and (a4+a5) == 0:
          O_env_status[i]=12
        elif a2 == 0 and a3 == 2 and (a4+a5) == 0:
          O_env_status[i]=13
        else:
          O_env_status[i]=10
        cd[O_env_status[i]-1] = cd[O_env_status[i]-1] + 1
 
    output1=open(BASE_FILE+str("_Oxygen_Environment_v01.data"),'w')
    output1.write("# count of each oxygen type based on the environment it has\n")
    output1.write("# 1=SiAlR; 2=2SiR; 3=2AlR; 4=SiAl2R; 5=2Al2R; 6=Si2R; 7=Si3R; 8=Al3R; 9=Al2R; 10=other; 11=SiOSi; 12=SiOAl; 13=AlOAl \n")
    output1.write("# types 1 to 5 are bridging oxygens and 6 to 9 are non bridging oxygens \n")
    for i in xrange(ca3+3):
      output1.write("%d \t %d \t %2.2lf\n" %(i+1,cd[i],cd[i]*100.0/total_atoms_of_type[1]))
    output1.close()
    for i in xrange(ca3):
      output1=open(BASE_FILE+str("_Oxygen_Environment_type")+str(i+1)+str(".atoms"),'w')
      output2=open(BASE_FILE+str("_Oxygen_Environment_type")+str(i+1)+str("_neighbours.atoms"),'w')
      write_imd_header3(output1,imd_box,atype,rc)
      write_imd_header3(output2,imd_box,atype,rc)
      SelectAtoms1=[0 for k in xrange(TotalAtoms)]
      for j in xrange(TotalAtoms):
        if O_env_status[j] == i+1:
          output1.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
          SelectAtoms1[j]=1
          for k in nnl[j]:
            SelectAtoms1[k]=1
      for j in xrange(TotalAtoms):
        if SelectAtoms1[j] == 1:
          output2.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
    output1.close()
    output2.close()
 
    output1=open(BASE_FILE+str("_BO_additional_info_v01.data"),'w')
    output1.write("# bridging oxygen additional information\n")
    output1.write("# O_atom_number\t NumberOf_O_neighbours\tNumberOf_Si_neighbours\tNumberOf_Al_neighbours\tNumberOf_Ca_neighbours\tNumberOf_Mg_neighbours\n")
    for i in xrange(TotalAtoms):
      if atom_type[i] == 1:
        if atom_nnl_type_count[i][2] == 2:
          output1.write("%d %d %d %d %d %d \n" %(atom_global_id[i],atom_nnl_type_count[i][0],atom_nnl_type_count[i][1],atom_nnl_type_count[i][2],atom_nnl_type_count[i][3],atom_nnl_type_count[i][4]))
    output1.close()

    output1=open(BASE_FILE+str("_BO_NBO_info_v01.data"),'w')
    output1.write("# bridging and non-bridging oxygen additional information\n")
    output2=open(BASE_FILE+str("_O_environment_info_of_Ca_Mg_v01.data"),'w')
    output2.write("# Ca or Mg ( nnl ) ( atom_number ) => each neighboring oxygen environment type\n")
    for a2 in xrange(ca1):
      for a3 in xrange(ca1):
        for a4 in xrange(2*ca1):
          if count_BO1[a2][a3][a4] != 0:
            output1.write("%d Si + %d Al + %d R = %d \n" %(a2,a3,a4,count_BO1[a2][a3][a4]))
    output1.write("\n\n##############################################\n\n")
    for a2 in xrange(ca1):
      for a3 in xrange(ca1):
        for a4 in xrange(ca1):
          for a5 in xrange(ca1):
            if count_BO[a2][a3][a4][a5] != 0:
               output1.write("%d Si + %d Al + %d Ca + %d Mg = %d \n" %(a2,a3,a4,a5,count_BO[a2][a3][a4][a5]))
    output1.write("\n\n##############################################\n\n")
    if Enable_O_env_status_loop == "yes":
      # compute the maximum size of each type of oxygen to define the matrix
      for i in xrange(TotalAtoms):
        if atom_type[i] == 4 or atom_type[i] == 5:
          a1=[0 for i1 in xrange(ca3)]
          for j in nnl[i]:		# j = neighbouring atom number
            a1[O_env_status[j]-1] = a1[O_env_status[j]-1] + 1	# counting how many O atoms of atom i have the same environment
          for j in xrange(ca3):
            if a1[j] > cb[j]-1:
              cb[j]=a1[j]+1
      print "\n\n*****************************************************************************************************************"
      print "* The max elements of matrix for computing type of Oxygen based on its environment is                           *"
      print "* count_O_env_status_of_R",
      for i in xrange(ca3):
        print "[%d]" %(cb[i]),
      print "                                              *"
      print "*****************************************************************************************************************\n\n"
      for i in cb:
        if i >= 10:
          print "\n\n*****************************************************************************************************************\n\n"
          print " The dimension of a big matrix to compute the type of Oxygen based on its environment is too big and thus the system \n"
          print " can hang if there is not enough ram therefore please enable this loop if you are sure that you have enough memory   \n\n"
          print "*****************************************************************************************************************\n\n"
          Enable_O_env_status_loop == "no"
  
    if Enable_O_env_status_loop == "yes":
      count_O_env_status_of_R=[[[[[[[[[[0 for i1 in xrange(cb[ca3-1])] for i2 in xrange(cb[ca3-2])] for i3 in xrange(cb[ca3-3])] for i4 in xrange(cb[ca3-4])] for i5 in xrange(cb[ca3-5])] for i6 in xrange(cb[ca3-6])] for i7 in xrange(cb[ca3-7])] for i8 in xrange(cb[ca3-8])] for i9 in xrange(cb[ca3-9])] for i10 in xrange(cb[ca3-10])] 
      for i in xrange(TotalAtoms):
        if atom_type[i] == 4 or atom_type[i] == 5:
          a1=[0 for i1 in xrange(ca3)]
          for j in nnl[i]:          # j = neighbouring atom number
            a1[O_env_status[j]-1] = a1[O_env_status[j]-1] + 1       # counting how many O atoms of atom i have the same environment 
          count_O_env_status_of_R[a1[0]][a1[1]][a1[2]][a1[3]][a1[4]][a1[5]][a1[6]][a1[7]][a1[8]][a1[9]]=count_O_env_status_of_R[a1[0]][a1[1]][a1[2]][a1[3]][a1[4]][a1[5]][a1[6]][a1[7]][a1[8]][a1[9]]+1

      for i1 in xrange(cb[0]):
       for i2 in xrange(cb[1]):
        for i3 in xrange(cb[2]):
         for i4 in xrange(cb[3]):
          for i5 in xrange(cb[4]):
           for i6 in xrange(cb[5]):
            for i7 in xrange(cb[6]):
             for i8 in xrange(cb[7]):
              for i9 in xrange(cb[8]):
               for i10 in xrange(cb[9]):
                   #for i11 in xrange(ca4):
                   # for i12 in xrange(ca4):
                   ca5=0
                   if Enable_O_env_status_loop == "yes":
                     ca5=count_O_env_status_of_R[i1][i2][i3][i4][i5][i6][i7][i8][i9][i10]
                   #ca5=count_O_env_status_of_R[i1][i2][i3][i4][i5][i6][i7][i8][i9][i10][i11]
                   #ca5=count_O_env_status_of_R[i1][i2][i3][i4][i5][i6][i7][i8][i9][i10][i11][i12]
                   if ca5 != 0:
                     output1.write("%d %d %d %d %d %d %d %d %d %d  = %d\n" %(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,ca5))
                     #output1.write("%d %d %d %d %d %d %d %d %d %d %d  = %d\n" %(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,ca5))
                     #output1.write("%d %d %d %d %d %d %d %d %d %d %d %d  = %d\n" %(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,ca5))                 
                     if i1 != 0:
                       cc[0]=cc[0]+ca5
                     if i2 != 0:
                       cc[1]=cc[1]+ca5
                     if i3 != 0:
                       cc[2]=cc[2]+ca5
                     if i4 != 0:
                       cc[3]=cc[3]+ca5
                     if i5 != 0:
                       cc[4]=cc[4]+ca5
                     if i6 != 0:
                       cc[5]=cc[5]+ca5
                     if i7 != 0:
                       cc[6]=cc[6]+ca5
                     if i8 != 0:
                       cc[7]=cc[7]+ca5
                     if i9 != 0:
                       cc[8]=cc[8]+ca5
                     if i10 != 0:
                       cc[9]=cc[9]+ca5
      # writing out percentage of R atoms which does not have a particular type of environment
      for i in xrange(ca3):
        output1.write("Number of R which has type %d oxygen environment = %d and their percentage is %2.2lf \n" %(i+1,cc[i],cc[i]*100.0/(total_atoms_of_type[4]+total_atoms_of_type[5])))
      output1.close()
      # writing out type of oxygen environment of neighbouring Ca and Mg atoms
      for i in xrange(TotalAtoms):
        if atom_type[i] != 1:
          for j in xrange(nnl_count[i]):
            atom_nnl_O_env[i][j]=O_env_status[nnl[i][j]]
          temp_list=list(atom_nnl_O_env[i])
          for j in xrange(nnl_count[i],MAX_NEIGHBOURS):
            temp_list.remove(-1)
          atom_nnl_O_env[i]=np.array(temp_list) 
          atom_nnl_O_env[i].sort()
      for i in xrange(TotalAtoms):
        if atom_type[i] == 4 or atom_type[i] == 5:
          output2.write("%s ( %d ) ( %d ) => " %(atype[atom_type[i]],nnl_count[i],atom_global_id[i]))
          for j in atom_nnl_O_env[i]:
            output2.write("%d  " %(j))
          output2.write("\n")
      output2.close()

    count_2Si=count_SiAl=count_2Al=count_Si=count_Al=count_3Al=count_Si2Al=count_missing_oxygen_atoms=0
    missing_oxygen_atoms=[]
    output1=open(BASE_FILE+str("_BO_NBO_info_v02.data"),'w')
    for i in xrange(TotalAtoms):
      if atom_type[i] == 1:
        if atom_nnl_type[i].count(2) == 2 and atom_nnl_type[i].count(3) == 0:
          count_2Si=count_2Si+1
        elif atom_nnl_type[i].count(2) == 1 and atom_nnl_type[i].count(3) == 1:
          count_SiAl=count_SiAl+1
        elif atom_nnl_type[i].count(2) == 0 and atom_nnl_type[i].count(3) == 2:
          count_2Al=count_2Al+1
        elif atom_nnl_type[i].count(2) == 1 and atom_nnl_type[i].count(3) == 0:
          count_Si=count_Si+1
        elif atom_nnl_type[i].count(2) == 0 and atom_nnl_type[i].count(3) == 1:
          count_Al=count_Al+1
        elif atom_nnl_type[i].count(2) == 0 and atom_nnl_type[i].count(3) == 3:
          count_3Al=count_3Al+1
        elif atom_nnl_type[i].count(2) == 1 and atom_nnl_type[i].count(3) == 2:
          count_Si2Al=count_Si2Al+1
        else:
          count_missing_oxygen_atoms=count_missing_oxygen_atoms+1
          missing_oxygen_atoms.append(i)
    temp_BO=count_2Si+count_SiAl+count_2Al
    # total oxygen atoms which are connected to any of Si-O-Si,Si-O-Al,Si-NBO
    for i in xrange(TotalAtoms):
      if atom_type[i] == 1:
        if atom_nnl_type[i].count(2) != 0:
          SelectAtoms1_O[i]=1
    for i in xrange(TotalAtoms):
      if SelectAtoms1_O[i] == 1:
        count_O_connected_to_Si=count_O_connected_to_Si+1
    output1.write("# BO and NBO information\n")
    output1.write("# Used cutoff radii for Si-O = %lf , Al-O = %lf , Ca-O = %lf , and Mg-O = %lf \n" %(rc_SiO,rc_AlO,rc_CaO,rc_MgO))
    output1.write("# Bridging Oxygen     (BO)   = %d \t %lf \n" %((count_2Si+count_SiAl+count_2Al),(count_2Si+count_SiAl+count_2Al)*100.0/total_atoms_of_type[1]))
    output1.write("# Non-Bridging Oxygen (NBO)  = %d \t %lf \n\n" %((count_Si+count_Al),(count_Si+count_Al)*100.0/total_atoms_of_type[1]))
    output1.write("# Si-BO-Si                   = %d \t %lf \n" %(count_2Si,count_2Si*100.0/total_atoms_of_type[1]))
    output1.write("# Si-BO-Al                   = %d \t %lf \n" %(count_SiAl,count_SiAl*100.0/total_atoms_of_type[1]))
    output1.write("# Al-BO-Al                   = %d \t %lf \n" %(count_2Al,count_2Al*100.0/total_atoms_of_type[1]))
    output1.write("# Si-NBO                     = %d \t %lf \n" %(count_Si,count_Si*100.0/total_atoms_of_type[1]))
    output1.write("# Al-NBO                     = %d \t %lf \n" %(count_Al,count_Al*100.0/total_atoms_of_type[1]))
    output1.write("# 3Al                        = %d \t %lf \n" %(count_3Al,count_3Al*100.0/total_atoms_of_type[1]))
    output1.write("# Si2Al                      = %d \t %lf \n" %(count_Si2Al,count_Si2Al*100.0/total_atoms_of_type[1]))
    output1.write("# Any other type oxygen      = %d \t %lf \n" %(count_missing_oxygen_atoms,count_missing_oxygen_atoms*100.0/total_atoms_of_type[1]))
    output1.write("\n\n# Total number of oxygen which are connected to more than 0 number of Si = %d \n\n" %(count_O_connected_to_Si))
    output1.write("# Following are the nnl of any other type of atoms which are listed above\n")
    for i in missing_oxygen_atoms:
      output1.write("\n%d => \t" %(atom_global_id[i]))
      for j in nnl[i]:
        output1.write("%d ( %d ) \t " %(atom_global_id[j],atom_type[j]))
    output1.close()

    # following loop is to write out Al-NBO atoms
    output1=open(BASE_FILE+str("_AlNBO.atoms"),'w')
    output2=open(BASE_FILE+str("_AlNBO_neighbours.atoms"),'w')
    write_imd_header3(output1,imd_box,atype,rc)
    write_imd_header3(output2,imd_box,atype,rc)
    SelectAtoms1=[0 for i in xrange(TotalAtoms)]
    SelectAtoms2=[0 for i in xrange(TotalAtoms)]
    for i in xrange(TotalAtoms):
      if atom_type[i] == 1:
        if atom_nnl_type_count[i][2] == 0 and atom_nnl_type_count[i][3] == 1:
          SelectAtoms1[i]=1
          for j in nnl[i]:
            #if atom_type[j] == 3 or atom_type[j] == 4 atom_type[j] == 5:
            SelectAtoms1[j]=1
    for i in xrange(TotalAtoms):
      if atom_type[i] != 1 and SelectAtoms1[i] == 1:
        SelectAtoms2[i] = 1
        for j in nnl[i]:
          SelectAtoms2[j] = 1
    for i in xrange(TotalAtoms):
      if SelectAtoms1[i] == 1:
        output1.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
      if SelectAtoms2[i] == 1:
        output2.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
    output1.close()
    output2.close()

    ########################################################################################################
    # computing Si-O-Al, Si-O-Si and Al-O-Al links again (omiting tri-clusters)
    SelectAtoms1=[0 for i in xrange(TotalAtoms)]
    SelectAtoms2=[0 for i in xrange(TotalAtoms)]
    SelectAtoms3=[0 for i in xrange(TotalAtoms)]
    SelectAtoms4=[0 for i in xrange(TotalAtoms)]
    SelectAtoms5=[0 for i in xrange(TotalAtoms)]
    SelectAtoms6=[0 for i in xrange(TotalAtoms)]
    output1=open(BASE_FILE+str("_SiOAl.atoms"),'w')
    output2=open(BASE_FILE+str("_SiOAl_SiAndAlTets.atoms"),'w')
    output3=open(BASE_FILE+str("_SiOSi.atoms"),'w')
    output4=open(BASE_FILE+str("_SiOSi_SiTets.atoms"),'w')
    output5=open(BASE_FILE+str("_AlOAl.atoms"),'w')
    output6=open(BASE_FILE+str("_AlOAl_AlTets.atoms"),'w')
    write_imd_header3(output1,imd_box,atype,rc)
    write_imd_header3(output2,imd_box,atype,rc)
    write_imd_header3(output3,imd_box,atype,rc)
    write_imd_header3(output4,imd_box,atype,rc)
    write_imd_header3(output5,imd_box,atype,rc)
    write_imd_header3(output6,imd_box,atype,rc)
    for i in xrange(TotalAtoms):
      if atom_type[i] == 1:
         # Si-O-Al
         if atom_nnl_type[i].count(2) == 1 and atom_nnl_type[i].count(3) == 1:
           SelectAtoms1[i]=1
           SelectAtoms2[i]=1
           for j in nnl[i]:
             if atom_type[j] == 2 or atom_type[j] == 3:
               SelectAtoms1[j]=1
               SelectAtoms2[j]=1
               for k in nnl[j]:
                 SelectAtoms2[k]=1
         elif atom_nnl_type[i].count(2) == 2 and atom_nnl_type[i].count(3) ==0:
           SelectAtoms3[i]=1
           SelectAtoms4[i]=1
           for j in nnl[i]:
             if atom_type[j] == 2:
               SelectAtoms3[j]=1
               SelectAtoms4[j]=1
               for k in nnl[j]:
                 SelectAtoms4[k]=1
         elif atom_nnl_type[i].count(2) == 0 and atom_nnl_type[i].count(3) == 2:
           SelectAtoms5[i]=1
           SelectAtoms6[i]=1
           for j in nnl[i]:
             if atom_type[j] == 3:
               SelectAtoms5[j]=1
               SelectAtoms6[j]=1
               for k in nnl[j]:
                 SelectAtoms6[k]=1
    for i in xrange(TotalAtoms):
      if SelectAtoms1[i] == 1:
        output1.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
      if SelectAtoms2[i] == 1:
        output2.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
      if SelectAtoms3[i] == 1:
        output3.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
      if SelectAtoms4[i] == 1:
        output4.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
      if SelectAtoms5[i] == 1:
        output5.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
      if SelectAtoms6[i] == 1:
        output6.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
    output1.close()
    output2.close()
    output3.close()
    output4.close()
    output5.close()
    output6.close()

    ################################################################################
    # Writing a part of the code again for Ca and Mg environment based on BO and NBO
    Ca_SiO=[0 for i1 in xrange(TotalAtoms)]
    Ca_AlO=[0 for i1 in xrange(TotalAtoms)]
    Ca_SiOSi=[0 for i1 in xrange(TotalAtoms)]
    Ca_SiOAl=[0 for i1 in xrange(TotalAtoms)]
    Ca_AlOAl=[0 for i1 in xrange(TotalAtoms)]
    Ca_Other=[0 for i1 in xrange(TotalAtoms)]
    Mg_SiO=[0 for i1 in xrange(TotalAtoms)]
    Mg_AlO=[0 for i1 in xrange(TotalAtoms)]
    Mg_SiOSi=[0 for i1 in xrange(TotalAtoms)]
    Mg_SiOAl=[0 for i1 in xrange(TotalAtoms)]
    Mg_AlOAl=[0 for i1 in xrange(TotalAtoms)]
    Mg_Other=[0 for i1 in xrange(TotalAtoms)]
    R_SiO=[0 for i1 in xrange(TotalAtoms)]
    R_AlO=[0 for i1 in xrange(TotalAtoms)]
    R_SiOSi=[0 for i1 in xrange(TotalAtoms)]
    R_SiOAl=[0 for i1 in xrange(TotalAtoms)]
    R_AlOAl=[0 for i1 in xrange(TotalAtoms)]
    R_Other=[0 for i1 in xrange(TotalAtoms)]
    SelectAtoms_R_1SiO=[0 for i in xrange(TotalAtoms)]
    SelectAtoms_R_2SiO=[0 for i in xrange(TotalAtoms)]
    SelectAtoms_R_3SiO=[0 for i in xrange(TotalAtoms)]
    SelectAtoms_R_4SiO=[0 for i in xrange(TotalAtoms)]
    SelectAtoms_R_5SiO=[0 for i in xrange(TotalAtoms)]
    SelectAtoms_R_6SiO=[0 for i in xrange(TotalAtoms)]
    count_SiNBO_per_Ca=count_AlNBO_per_Ca=count_SiOSi_per_Ca=count_SiOAl_per_Ca=count_AlOAl_per_Ca=count_Other_per_Ca=0
    count_SiNBO_per_Mg=count_AlNBO_per_Mg=count_SiOSi_per_Mg=count_SiOAl_per_Mg=count_AlOAl_per_Mg=count_Other_per_Mg=0

    for i in xrange(TotalAtoms):
      # For Ca atoms
      Ca_SiO_count = Ca_AlO_count = Ca_SiOSi_count = Ca_SiOAl_count = Ca_AlOAl_count = Ca_Other_count = 0
      if atom_type[i] == 4:
        for j in nnl[i]:
          if atom_nnl_type[j].count(2) == 1 and atom_nnl_type[j].count(3) == 0:
            Ca_SiO_count = Ca_SiO_count + 1
          elif atom_nnl_type[j].count(2) == 0 and atom_nnl_type[j].count(3) == 1:
            Ca_AlO_count = Ca_AlO_count + 1
          elif atom_nnl_type[j].count(2) == 2 and atom_nnl_type[j].count(3) == 0:
            Ca_SiOSi_count = Ca_SiOSi_count + 1
          elif atom_nnl_type[j].count(2) == 1 and atom_nnl_type[j].count(3) == 1:
            Ca_SiOAl_count = Ca_SiOAl_count + 1
          elif atom_nnl_type[j].count(2) == 0 and atom_nnl_type[j].count(3) == 2:
            Ca_AlOAl_count = Ca_AlOAl_count + 1
          else :
            Ca_Other_count = Ca_Other_count + 1
        Ca_SiO[i]=Ca_SiO_count
        Ca_AlO[i]=Ca_AlO_count
        Ca_SiOSi[i]=Ca_SiOSi_count
        Ca_SiOAl[i]=Ca_SiOAl_count
        Ca_AlOAl[i]=Ca_AlOAl_count
        Ca_Other[i]=Ca_Other_count
        count_SiNBO_per_Ca = count_SiNBO_per_Ca + Ca_SiO_count
        count_AlNBO_per_Ca = count_AlNBO_per_Ca + Ca_AlO_count
        count_SiOSi_per_Ca = count_SiOSi_per_Ca + Ca_SiOSi_count
        count_SiOAl_per_Ca = count_SiOAl_per_Ca + Ca_SiOAl_count
        count_AlOAl_per_Ca = count_AlOAl_per_Ca + Ca_AlOAl_count
        count_Other_per_Ca = count_Other_per_Ca + Ca_Other_count

      # For Mg atoms
      Mg_SiO_count = Mg_AlO_count = Mg_SiOSi_count = Mg_SiOAl_count = Mg_AlOAl_count = Mg_Other_count = 0
      if atom_type[i] == 5:
        for j in nnl[i]:
          if atom_nnl_type[j].count(2) == 1 and atom_nnl_type[j].count(3) == 0:
            Mg_SiO_count = Mg_SiO_count + 1
          elif atom_nnl_type[j].count(2) == 0 and atom_nnl_type[j].count(3) == 1:
            Mg_AlO_count = Mg_AlO_count + 1
          elif atom_nnl_type[j].count(2) == 2 and atom_nnl_type[j].count(3) == 0:
            Mg_SiOSi_count = Mg_SiOSi_count + 1
          elif atom_nnl_type[j].count(2) == 1 and atom_nnl_type[j].count(3) == 1:
            Mg_SiOAl_count = Mg_SiOAl_count + 1
          elif atom_nnl_type[j].count(2) == 0 and atom_nnl_type[j].count(3) == 2:
            Mg_AlOAl_count = Mg_AlOAl_count + 1
          else :
            Mg_Other_count=Mg_Other_count+1
        Mg_SiO[i]=Mg_SiO_count
        Mg_AlO[i]=Mg_AlO_count
        Mg_SiOSi[i]=Mg_SiOSi_count
        Mg_SiOAl[i]=Mg_SiOAl_count
        Mg_AlOAl[i]=Mg_AlOAl_count
        Mg_Other[i]=Mg_Other_count
        count_SiNBO_per_Mg = count_SiNBO_per_Mg + Mg_SiO_count
        count_AlNBO_per_Mg = count_AlNBO_per_Mg + Mg_AlO_count
        count_SiOSi_per_Mg = count_SiOSi_per_Mg + Mg_SiOSi_count
        count_SiOAl_per_Mg = count_SiOAl_per_Mg + Mg_SiOAl_count
        count_AlOAl_per_Mg = count_AlOAl_per_Mg + Mg_AlOAl_count
        count_Other_per_Mg = count_Other_per_Mg + Mg_Other_count

      # For Ca or Mg
      R_SiO_count = R_AlO_count = R_SiOSi_count = R_SiOAl_count = R_AlOAl_count = R_Other_count = 0
      if atom_type[i] == 4 or atom_type[i] == 5:
        for j in nnl[i]:
          if atom_nnl_type[j].count(2) == 1 and atom_nnl_type[j].count(3) == 0:
            R_SiO_count = R_SiO_count + 1
          elif atom_nnl_type[j].count(2) == 0 and atom_nnl_type[j].count(3) == 1:
            R_AlO_count = R_AlO_count + 1
          elif atom_nnl_type[j].count(2) == 2 and atom_nnl_type[j].count(3) == 0:
            R_SiOSi_count = R_SiOSi_count + 1
          elif atom_nnl_type[j].count(2) == 1 and atom_nnl_type[j].count(3) == 1:
            R_SiOAl_count = R_SiOAl_count + 1
          elif atom_nnl_type[j].count(2) == 0 and atom_nnl_type[j].count(3) == 2:
            R_AlOAl_count = R_AlOAl_count + 1
          else :
            R_Other_count = R_Other_count + 1
        R_SiO[i]=R_SiO_count
        R_AlO[i]=R_AlO_count
        R_SiOSi[i]=R_SiOSi_count
        R_SiOAl[i]=R_SiOAl_count
        R_AlOAl[i]=R_AlOAl_count
        R_Other[i]=R_Other_count

    Ca_SiO_count=[0 for i in xrange(10)]
    Ca_AlO_count=[0 for i in xrange(10)]
    Ca_SiOSi_count=[0 for i in xrange(10)]
    Ca_SiOAl_count=[0 for i in xrange(10)]
    Ca_AlOAl_count=[0 for i in xrange(10)]
    Ca_Other_count=[0 for i in xrange(10)]
    Mg_SiO_count=[0 for i in xrange(10)]
    Mg_AlO_count=[0 for i in xrange(10)]
    Mg_SiOSi_count=[0 for i in xrange(10)]
    Mg_SiOAl_count=[0 for i in xrange(10)]
    Mg_AlOAl_count=[0 for i in xrange(10)]
    Mg_Other_count=[0 for i in xrange(10)]
    R_SiO_count=[0 for i in xrange(10)]
    R_AlO_count=[0 for i in xrange(10)]
    R_SiOSi_count=[0 for i in xrange(10)]
    R_SiOAl_count=[0 for i in xrange(10)]
    R_AlOAl_count=[0 for i in xrange(10)]
    R_Other_count=[0 for i in xrange(10)]
 
    for i in xrange(TotalAtoms):
      if atom_type[i] == 4:
        Ca_SiO_count[Ca_SiO[i]]     = Ca_SiO_count[Ca_SiO[i]] + 1
        Ca_AlO_count[Ca_AlO[i]]     = Ca_AlO_count[Ca_AlO[i]] + 1
        Ca_SiOSi_count[Ca_SiOSi[i]] = Ca_SiOSi_count[Ca_SiOSi[i]] + 1
        Ca_SiOAl_count[Ca_SiOAl[i]] = Ca_SiOAl_count[Ca_SiOAl[i]] + 1
        Ca_AlOAl_count[Ca_AlOAl[i]] = Ca_AlOAl_count[Ca_AlOAl[i]] + 1
        Ca_Other_count[Ca_Other[i]] = Ca_Other_count[Ca_Other[i]] + 1
      if atom_type[i] == 5:
        Mg_SiO_count[Mg_SiO[i]]     = Mg_SiO_count[Mg_SiO[i]]+1
        Mg_AlO_count[Mg_AlO[i]]     = Mg_AlO_count[Mg_AlO[i]]+1
        Mg_SiOSi_count[Mg_SiOSi[i]] = Mg_SiOSi_count[Mg_SiOSi[i]]+1
        Mg_SiOAl_count[Mg_SiOAl[i]] = Mg_SiOAl_count[Mg_SiOAl[i]]+1
        Mg_AlOAl_count[Mg_AlOAl[i]] = Mg_AlOAl_count[Mg_AlOAl[i]]+1
        Mg_Other_count[Mg_Other[i]] = Mg_Other_count[Mg_Other[i]]+1
      if atom_type[i] == 4 or atom_type[i] == 5:
        R_SiO_count[R_SiO[i]]     = R_SiO_count[R_SiO[i]]+1
        R_AlO_count[R_AlO[i]]     = R_AlO_count[R_AlO[i]]+1
        R_SiOSi_count[R_SiOSi[i]] = R_SiOSi_count[R_SiOSi[i]]+1
        R_SiOAl_count[R_SiOAl[i]] = R_SiOAl_count[R_SiOAl[i]]+1
        R_AlOAl_count[R_AlOAl[i]] = R_AlOAl_count[R_AlOAl[i]]+1
        R_Other_count[R_Other[i]] = R_Other_count[R_Other[i]]+1

    for i in xrange(TotalAtoms):
      if R_SiO[i] == 1:
        SelectAtoms_R_1SiO[i]=1
      elif R_SiO[i] == 2:
        SelectAtoms_R_2SiO[i]=1
      elif R_SiO[i] == 3:
        SelectAtoms_R_3SiO[i]=1
      elif R_SiO[i] == 4:
        SelectAtoms_R_4SiO[i]=1
      elif R_SiO[i] == 5:
        SelectAtoms_R_5SiO[i]=1
      elif R_SiO[i] == 6:
        SelectAtoms_R_6SiO[i]=1

    output1=open(BASE_FILE+str("_R_BO_and_NBO_info_v01.data"),'w')
    output1.write("# BO and NBO of Ca info: Number of Ca atoms\t\t %of Ca atoms\n\n")
    for i in xrange(10):
      if Ca_SiO_count[i] != 0:
        output1.write("Number of Ca atoms with %d SiO (NBO) = %d    %2.2lf\n" %(i,Ca_SiO_count[i],100.0*Ca_SiO_count[i]/total_atoms_of_type[4]))
    output1.write("\n")
    for i in xrange(10):
      if Ca_AlO_count[i] != 0:
        output1.write("Number of Ca atoms with %d AlO (NBO) = %d    %2.2lf\n" %(i,Ca_AlO_count[i],100.0*Ca_AlO_count[i]/total_atoms_of_type[4]))
    output1.write("\n")
    for i in xrange(10):
      if Ca_SiOSi_count[i] != 0:
        output1.write("Number of Ca atoms with %d SiOSi (BO) = %d    %2.2lf\n" %(i,Ca_SiOSi_count[i],100.0*Ca_SiOSi_count[i]/total_atoms_of_type[4]))
    output1.write("\n")
    for i in xrange(10):
      if Ca_SiOAl_count[i] != 0:
        output1.write("Number of Ca atoms with %d SiOAl (BO) = %d    %2.2lf\n" %(i,Ca_SiOAl_count[i],100.0*Ca_SiOAl_count[i]/total_atoms_of_type[4]))
    output1.write("\n")
    for i in xrange(10):
      if Ca_AlOAl_count[i] != 0:
        output1.write("Number of Ca atoms with %d AlOAl (BO) = %d    %2.2lf\n" %(i,Ca_AlOAl_count[i],100.0*Ca_AlOAl_count[i]/total_atoms_of_type[4]))
    output1.write("\n")
    for i in xrange(10):
      if Ca_Other_count[i] != 0:
        output1.write("Number of Ca atoms with %d Other than above Oxygen status = %d    %2.2lf\n" %(i,Ca_Other_count[i],100.0*Ca_Other_count[i]/total_atoms_of_type[4]))
    output1.write("\n\n")

    output1.write("#####################################################")
    output1.write("\n# BO and NBO of Mg info: Number\t %of Mg atoms\n")
    for i in xrange(10):
      if Mg_SiO_count[i] != 0:
        output1.write("Number of Mg atoms with %d SiO (NBO) = %d     %2.2lf\n" %(i,Mg_SiO_count[i],100.0*Mg_SiO_count[i]/total_atoms_of_type[5]))
    output1.write("\n")
    for i in xrange(10):
      if Mg_AlO_count[i] != 0:
        output1.write("Number of Mg atoms with %d AlO (NBO) = %d     %2.2lf\n" %(i,Mg_AlO_count[i],100.0*Mg_AlO_count[i]/total_atoms_of_type[5]))
    output1.write("\n")
    for i in xrange(10):
      if Mg_SiOSi_count[i] != 0:
        output1.write("Number of Mg atoms with %d SiOSi (BO) = %d     %2.2lf\n" %(i,Mg_SiOSi_count[i],100.0*Mg_SiOSi_count[i]/total_atoms_of_type[5]))
    output1.write("\n")
    for i in xrange(10):
      if Mg_SiOAl_count[i] != 0:
        output1.write("Number of Mg atoms with %d SiOAl (BO) = %d     %2.2lf\n" %(i,Mg_SiOAl_count[i],100.0*Mg_SiOAl_count[i]/total_atoms_of_type[5]))
    output1.write("\n")
    for i in xrange(10):
      if Mg_AlOAl_count[i] != 0:
        output1.write("Number of Mg atoms with %d AlOAl (BO) = %d     %2.2lf\n" %(i,Mg_AlOAl_count[i],100.0*Mg_AlOAl_count[i]/total_atoms_of_type[5]))
    output1.write("\n")
    for i in xrange(10):
      if Mg_Other_count[i] != 0:
        output1.write("Number of Mg atoms with %d Other than above Oxygen status = %d     %2.2lf\n" %(i,Mg_Other_count[i],100.0*Mg_Other_count[i]/total_atoms_of_type[5]))
    output1.write("\n\n")

    output1.write("#####################################################")
    output1.write("\n# BO and NBO of R info: Number\t %of R atoms\n")
    for i in xrange(10):
      if R_SiO_count[i] != 0:
        output1.write("Number of R atoms with %d SiO (NBO) = %d     %2.2lf\n" %(i,R_SiO_count[i],100.0*R_SiO_count[i]/(total_atoms_of_type[4]+total_atoms_of_type[5])))
    output1.write("\n")
    for i in xrange(10):
      if R_AlO_count[i] != 0:
        output1.write("Number of R atoms with %d AlO (NBO) = %d     %2.2lf\n" %(i,R_AlO_count[i],100.0*R_AlO_count[i]/(total_atoms_of_type[4]+total_atoms_of_type[5])))
    output1.write("\n")
    for i in xrange(10):
      if R_SiOSi_count[i] != 0:
        output1.write("Number of R atoms with %d SiOSi (BO) = %d     %2.2lf\n" %(i,R_SiOSi_count[i],100.0*R_SiOSi_count[i]/(total_atoms_of_type[4]+total_atoms_of_type[5])))
    output1.write("\n")
    for i in xrange(10):
      if R_SiOAl_count[i] != 0:
        output1.write("Number of R atoms with %d SiOAl (BO) = %d     %2.2lf\n" %(i,R_SiOAl_count[i],100.0*R_SiOAl_count[i]/(total_atoms_of_type[4]+total_atoms_of_type[5])))
    output1.write("\n")
    for i in xrange(10):
      if R_AlOAl_count[i] != 0:
        output1.write("Number of R atoms with %d AlOAl (BO) = %d     %2.2lf\n" %(i,R_AlOAl_count[i],100.0*R_AlOAl_count[i]/(total_atoms_of_type[4]+total_atoms_of_type[5])))
    output1.write("\n")
    for i in xrange(10):
      if R_Other_count[i] != 0:
        output1.write("Number of R atoms with %d Other than above Oxygen status = %d     %2.2lf\n" %(i,R_Other_count[i],100.0*R_Other_count[i]/(total_atoms_of_type[4]+total_atoms_of_type[5])))
    output1.write("\n\n")
 
    output1.write("###############################################\n")
    if total_atoms_of_type[4] != 0:
      output1.write("# Average type_of_oxygen per Ca\n")
      output1.write("Average SiNBO/Ca = %lf\n" %(1.0*count_SiNBO_per_Ca/total_atoms_of_type[4]))
      output1.write("Average AlNBO/Ca = %lf\n" %(1.0*count_AlNBO_per_Ca/total_atoms_of_type[4]))
      output1.write("Average SiOSi/Ca = %lf\n" %(1.0*count_SiOSi_per_Ca/total_atoms_of_type[4]))
      output1.write("Average SiOAl/Ca = %lf\n" %(1.0*count_SiOAl_per_Ca/total_atoms_of_type[4]))
      output1.write("Average AlOAl/Ca = %lf\n" %(1.0*count_AlOAl_per_Ca/total_atoms_of_type[4]))
      output1.write("Average Other_O/Ca = %lf\n" %(1.0*count_Other_per_Ca/total_atoms_of_type[4]))

    if total_atoms_of_type[5] != 0:
      output1.write("\n\n###############################################\n")
      output1.write("# Average type_of_oxygen per Mg\n")
      output1.write("Average SiNBO/Mg = %lf\n" %(1.0*count_SiNBO_per_Mg/total_atoms_of_type[5]))
      output1.write("Average AlNBO/Mg = %lf\n" %(1.0*count_AlNBO_per_Mg/total_atoms_of_type[5]))
      output1.write("Average SiOSi/Mg = %lf\n" %(1.0*count_SiOSi_per_Mg/total_atoms_of_type[5]))
      output1.write("Average SiOAl/Mg = %lf\n" %(1.0*count_SiOAl_per_Mg/total_atoms_of_type[5]))
      output1.write("Average AlOAl/Mg = %lf\n" %(1.0*count_AlOAl_per_Mg/total_atoms_of_type[5]))
      output1.write("Average Other_O/Mg = %lf\n" %(1.0*count_Other_per_Mg/total_atoms_of_type[5]))

    output1.close()
    output_R_1SiO=open(BASE_FILE+str("_R_1SiO.atoms"),'w')
    output_R_2SiO=open(BASE_FILE+str("_R_2SiO.atoms"),'w')
    output_R_3SiO=open(BASE_FILE+str("_R_3SiO.atoms"),'w')
    output_R_4SiO=open(BASE_FILE+str("_R_4SiO.atoms"),'w')
    output_R_5SiO=open(BASE_FILE+str("_R_5SiO.atoms"),'w')
    output_R_6SiO=open(BASE_FILE+str("_R_6SiO.atoms"),'w')
    write_imd_header4(output_R_1SiO,imd_box,atype,rc)
    write_imd_header4(output_R_2SiO,imd_box,atype,rc)
    write_imd_header4(output_R_3SiO,imd_box,atype,rc)
    write_imd_header4(output_R_4SiO,imd_box,atype,rc)
    write_imd_header4(output_R_5SiO,imd_box,atype,rc)
    write_imd_header4(output_R_6SiO,imd_box,atype,rc)
    for i in xrange(TotalAtoms):
      output_R_1SiO.write("%d  %d  %lf  %lf  %lf %d %d\n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i],SelectAtoms_R_1SiO[i]))
      output_R_2SiO.write("%d  %d  %lf  %lf  %lf %d %d\n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i],SelectAtoms_R_2SiO[i]))
      output_R_3SiO.write("%d  %d  %lf  %lf  %lf %d %d\n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i],SelectAtoms_R_3SiO[i]))
      output_R_4SiO.write("%d  %d  %lf  %lf  %lf %d %d\n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i],SelectAtoms_R_4SiO[i]))
      output_R_5SiO.write("%d  %d  %lf  %lf  %lf %d %d\n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i],SelectAtoms_R_5SiO[i]))
      output_R_6SiO.write("%d  %d  %lf  %lf  %lf %d %d\n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i],SelectAtoms_R_6SiO[i]))

    output_R_1SiO.close()
    output_R_2SiO.close()
    output_R_3SiO.close()
    output_R_4SiO.close()
    output_R_5SiO.close()
    output_R_6SiO.close()

  if Enable_phase2_calculations == "yes": 
    #######################################################################################################
    # working on bridging oxygens again
    # May be I am rewriting the part of existing code again on 16/12/2018
    SiOSi_count=[[0 for i in xrange(8)] for j in xrange(8)]	# Si(i)OSi(j)
    SiOAl_count=[[0 for i in xrange(8)] for j in xrange(8)]     # Si(i)OAl(j)
    AlOAl_count=[[0 for i in xrange(8)] for j in xrange(8)]     # Al(i)OAl(j)
    Ca_count_with_NBO = 0
    Mg_count_with_NBO = 0
    NBO_count=0
    NBO_count1=0
    SiNBO_count1=0
    AlNBO_count1=0
    BO_count=0
    output1=open(BASE_FILE+str("_BO_NBO_info_v03.data"),'w')
    output1.write("# BO and NBO information\n")
    for i in xrange(TotalAtoms):
      if atom_type[i] == 1:
        # Si-O-Si
        if atom_nnl_type[i].count(2) == 2:
          Si_count1=[0 for i1 in xrange(2)]
          Si_count2=0
          for j in nnl[i]:
            if atom_type[j] == 2 :
              Si_count1[Si_count2] = nnl_count[j]
              Si_count2=Si_count2+1
          if Si_count1[0] <= Si_count1[1] :
            SiOSi_count[Si_count1[0]][Si_count1[1]]=SiOSi_count[Si_count1[0]][Si_count1[1]]+1
          else :
            SiOSi_count[Si_count1[1]][Si_count1[0]]=SiOSi_count[Si_count1[1]][Si_count1[0]]+1

        # Si-O-Al
        if atom_nnl_type[i].count(2) == 1 and atom_nnl_type[i].count(3) == 1:
          Si_count1=0
          Al_count1=0
          for j in nnl[i]:
            if atom_type[j] == 2 :
              Si_count1 = nnl_count[j]
            elif atom_type[j] == 3 :
              Al_count1 = nnl_count[j]
          SiOAl_count[Si_count1][Al_count1]=SiOAl_count[Si_count1][Al_count1]+1

        # Al-O-Al
        if atom_nnl_type[i].count(3) == 2:
          Al_count1=[0 for i1 in xrange(2)]
          Al_count2=0
          for j in nnl[i]:
            if atom_type[j] == 3 :
              Al_count1[Al_count2] = nnl_count[j]
              Al_count2=Al_count2+1
          if Al_count1[0] <= Al_count1[1] :
            AlOAl_count[Al_count1[0]][Al_count1[1]]=AlOAl_count[Al_count1[0]][Al_count1[1]]+1
          else :
            AlOAl_count[Al_count1[1]][Al_count1[0]]=AlOAl_count[Al_count1[1]][Al_count1[0]]+1

        # Number of R atoms connected to NBO when the O is having  Si-O-(n)R or Al-O-(n)R, where n could be any number
        if (atom_nnl_type[i].count(2)+atom_nnl_type[i].count(3)) == 1:
          Ca_count_with_NBO = Ca_count_with_NBO + atom_nnl_type[i].count(4)
          Mg_count_with_NBO = Mg_count_with_NBO + atom_nnl_type[i].count(5)
          NBO_count1=NBO_count1+1
          if atom_nnl_type[i].count(2) == 1:
            SiNBO_count1=SiNBO_count1+1
          elif atom_nnl_type[i].count(3) == 1:
            AlNBO_count1=AlNBO_count1+1

    BO_count = BO_count + SiOSi_count[4][4] + SiOAl_count[4][4] + AlOAl_count[4][4]
    NBO_count = total_atoms_of_type[1]-BO_count
    output1.write("# Total Number of Oxygen Atoms                                          = %d  \n" %(total_atoms_of_type[1]))
    output1.write("# Total Number of Bridging Oxygen i.e Si4-O-Si4 + Si4-O-Al4 + Al4-O-Al4 = %d \t %2.2lf\n" %(BO_count,  100.0* BO_count/total_atoms_of_type[1]))
    output1.write("# Total Number of Non-Bridging Oxygen (including tri-clusters)          = %d \t %2.2lf\n" %(NBO_count, 100.0*NBO_count/total_atoms_of_type[1]))

    output1.write("# SiOSi \n")
    for i in xrange(8):
      for j in xrange(8):
        if SiOSi_count[i][j] != 0:
          output1.write("Si [%d] O Si [%d] = %d  \t %2.2lf \n" %(i,j,SiOSi_count[i][j],100.0*SiOSi_count[i][j]/total_atoms_of_type[1]))
    output1.write("\n# SiOAl \n")
    for i in xrange(8):
      for j in xrange(8):
        if SiOAl_count[i][j] != 0:
          output1.write("Si [%d] O Al [%d] = %d  \t %2.2lf \n" %(i,j,SiOAl_count[i][j],100.0*SiOAl_count[i][j]/total_atoms_of_type[1]))
    output1.write("\n# AlOAl \n")
    for i in xrange(8):
      for j in xrange(8):
        if AlOAl_count[i][j] != 0:
          output1.write("Al [%d] O Al [%d] = %d  \t %2.2lf \n" %(i,j,AlOAl_count[i][j],100.0*AlOAl_count[i][j]/total_atoms_of_type[1]))
    output1.write("# Number of NBO (excluding tri-cluster, i.e only Si-O or Al-O)    = %d   %2.2lf \n" %(NBO_count1,100.0*NBO_count1/total_atoms_of_type[1]))
    output1.write("# Number of Si-NBO (excluding tri-cluster, i.e only Si-O )        = %d   %2.2lf \n" %(SiNBO_count1,100.0*SiNBO_count1/total_atoms_of_type[1]))
    output1.write("# Number of Al-NBO (excluding tri-cluster, i.e only Al-O )        = %d   %2.2lf \n" %(AlNBO_count1,100.0*AlNBO_count1/total_atoms_of_type[1]))
    output1.write("# Number of Ca atoms connected with NBO (where NBO= Si-O or Al-O) = %d   %2.2lf \n" %(Ca_count_with_NBO,100.0*Ca_count_with_NBO/total_atoms_of_type[1]))
    output1.write("# Number of Mg atoms connected with NBO (where NBO= Si-O or Al-O) = %d   %2.2lf \n" %(Mg_count_with_NBO,100.0*Mg_count_with_NBO/total_atoms_of_type[1]))


