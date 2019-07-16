#!/usr/bin/python
"""
@Author : Sudheer Ganisetti
@Date   : Mo 25. Dez 18:34:13 CET 2017
This code is to compute nnl, nnlchk,  connectivity(Qn)  nd other properties of NAP glass (Na2O + Al2O3 + P2O5)
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
def write_imd_header4(output,box,atype1,r1,a1):
  #output
  if a1 == 0:
    output.write("#F A 1 1 3 1 0\n")
    output.write("#C number type x y z tot_nn \n")
  elif a1 == 1:
    output.write("#F A 1 1 3 3 0\n")
    output.write("#C number type x y z tot_nn     m_in_QmAlnP	n_in_QmAlnP\n")
  elif a1 == 2:
    output.write("#F A 1 1 3 3 0\n")
    output.write("#C number type x y z tot_nn     m_in_QmPnAl	n_in_QmPnAl\n")
  output.write("#X %lf 0.0 0.0 \n" %(box[0]))
  output.write("#Y 0.0 %lf 0.0 \n" %(box[1]))
  output.write("#Z 0.0 0.0 %lf \n" %(box[2]))
  output.write("## rcut of %s-O = %lf, %s-O = %lf, %s-O = %lf \n" %(str(atype1[3]),r1[3][1],str(atype1[5]),r1[5][1],str(atype1[6]),r1[6][1]))
  output.write("#E \n")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if __name__=="__main__":
  
  tot_argv=len(sys.argv)
  if tot_argv != 2:
     subprocess.call("sudheer_banner")
     print "************* S. Ganisetti *************"
     print "Error: usage is wrong"
     print "./this_program  ChkptFile"
     print "The program is prepared for a maximum of 6 atom types"
     print "By default we are assuming following atoms types and cutoffs"
     print "O(1), Si(2), Al(3), Ca(4), Na(5), P(6) => Si-O= 1.85 ; Al-O= 2.0 ; CaO= 2.9 ; MgO= 2.7"
     print "Please specify the cutoff radii of all required pair of atoms"
     print "If you want more atom types please modify MAX_ATOM_TYPES and also give corresponding atom pair cutoff radii"
     print "****************************************"
     sys.exit(0)

  # --------------------------------------------------------------------
  # The cutoff radii are hard coaded, please edit according to your need
  MAX_NEIGHBOURS=15
  given_O_type=1
  given_Si_type=2
  given_Al_type=3
  given_Ca_type=4
  given_Na_type=5
  given_P_type=6

  given_atom_types=[1,3,5,6]
  rc_AlO=2.5
  rc_NaO=3.0
  rc_PO=1.7
  MAX_ATOM_TYPES=6
  atype=["","O","","Al","","Na","P"]
  rc=[[0.0 for i in xrange(MAX_ATOM_TYPES+1)] for j in xrange(MAX_ATOM_TYPES+1)]
  #rc[2][1] = rc[1][2] = rc_SiO
  rc[3][1] = rc[1][3] = rc_AlO
  rc[5][1] = rc[1][5] = rc_NaO
  rc[6][1] = rc[1][6] = rc_PO
  rcut=0.0
  for i in range(MAX_ATOM_TYPES+1):
    for j in range(MAX_ATOM_TYPES+1):
       rcut=max(rcut,rc[i][j])
  # for debugging
  temp87=0
  # --------------------------------------------------------------------

  BASE_FILE=str(sys.argv[1])
  LAMMPS_DUMP_FILE=BASE_FILE+str('.dump')
  foutput1=BASE_FILE+str('.nnl')
  output1=open(foutput1,'w')
  output1.write("# nnl file: used cutoffs are r_AlO = %lf ; r_NaO = %lf ; r_PO = %lf\n" %(rc[1][3],rc[1][5],rc[1][6])) 

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
  write_imd_header4(output,imd_box,atype,rc,0)
  for i in range(TotalAtoms):
    output.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
  output.close()

  # **************************************************************************************
  # writing out individual cations and their anions
  for i in range(2,MAX_ATOM_TYPES+1):
    if i in given_atom_types:
      output=open(BASE_FILE+str("_")+str(atype[i])+str("-O.atoms"),'w')
      #write_imd_header1(output,imd_box,str(atype[i]),rc[i][1],str(atype[i]),rc[i][1])
      write_imd_header4(output,imd_box,atype,rc,0)
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

  # **************************************************************************************
  # writing out bridging and non-bridging oxygen
#  NBOS=[0 for i in xrange(TotalAtoms)]									# Non-Bridging Oxygen Status
#  OCa_count=[0 for i in xrange(TotalAtoms)]								# Number of Ca connected to each oxygen
#  OMg_count=[0 for i in xrange(TotalAtoms)]								# Number of Mg connected to each oxygen
#  BOCa_2Si_count=0											# Number of Ca connected to bridging oxygen and 2 Si
#  BOCa_2Al_count=0											# Number of Ca connected to bridging oxygen and 2 Al
#  BOCa_SiAl_count=0											# Number of Ca connected to bridging oxygen and 1 Si+  1Al
#  BOMg_2Si_count=0                                                                                      # Number of Mg connected to bridging oxygen and 2 Si
#  BOMg_2Al_count=0                                                                                      # Number of Mg connected to bridging oxygen and 2 Al
#  BOMg_SiAl_count=0                                                                                     # Number of Mg connected to bridging oxygen and 1 Si+  1Al
#  NBOCa_Si_count=0                                                                                      # Number of Ca connected to non-bridging oxygen and 2 Si
#  NBOCa_Al_count=0                                                                                      # Number of Ca connected to non-bridging oxygen and 2 Al
#  NBOMg_Si_count=0
#  NBOMg_Al_count=0
#  NBOSi_count	=NBOAl_count	=NBO3Si_count	=0
#  BOSiOSi_count	=BOSiOAl_count	=BOAlOAl_count	=0
#  NBOI=[[[[0 for i in xrange(8)] for j in xrange(8)] for k in xrange(8)] for l in xrange(8)]		# Non-Bridging Oxygen Information
#  NBOI1=[[[0 for i in xrange(8)] for j in xrange(8)] for k in xrange(8)]				# Non-bridging Oxygen Information
#  BOI=[[[[0 for i in xrange(8)] for j in xrange(8)] for k in xrange(8)] for l in xrange(8)]             # Bridging Oxygen Information
#  BOI1=[[[0 for i in xrange(8)] for j in xrange(8)] for k in xrange(8)]                                 # Bridging Oxygen Information
#  AlOAl_NBO_Count1=[0 for i in xrange(10)]
#  AlOAl_NBO_Count2=[0 for i in xrange(10)]
#  AlOAl_NBO_Count3=[0 for i in xrange(10)]
#  AlOAl_NBO_Count4=[0 for i in xrange(10)]
#  output1=open(BASE_FILE+str("_non-bridging_oxygen_info.data"),'w')
#  output2=open(BASE_FILE+str("_bridging_oxygen_info.data"),'w')
#  output3=open(BASE_FILE+str("_bridging_oxygen_additional_info.data"),'w')
#  output11=open(BASE_FILE+str("_non-bridging_oxygen.atoms"),'w')
#  output12=open(BASE_FILE+str("_non-bridging_oxygen_and_its_tetra.atoms"),'w')
#  output21=open(BASE_FILE+str("_bridging_oxygen.atoms"),'w')
#  output22=open(BASE_FILE+str("_bridging_oxygen_and_its_tetra.atoms"),'w')
#  output31=open(BASE_FILE+str("_bridging_oxygen_with2Si.atoms"),'w')
#  output32=open(BASE_FILE+str("_bridging_oxygen_with2Si_and_its_tetra.atoms"),'w')
#  output33=open(BASE_FILE+str("_bridging_oxygen_with2Si0Al0R.atoms"),'w')
#  output41=open(BASE_FILE+str("_bridging_oxygen_with2Al.atoms"),'w')
#  output42=open(BASE_FILE+str("_bridging_oxygen_with2Al_and_it_tetra.atoms"),'w')
#  output51=open(BASE_FILE+str("_bridging_oxygen_with1Si1Al.atoms"),'w')
#  output52=open(BASE_FILE+str("_bridging_oxygen_with1Si1Al_and_it_tetra.atoms"),'w')
#  output61=open(BASE_FILE+str("_non-bridging_oxygen_with0Si1AlnR.atoms"),'w')
#  output62=open(BASE_FILE+str("_non-bridging_oxygen_with0Si1AlnR_and_its_tetra.atoms"),'w')
#  write_imd_header3(output11,imd_box,atype,rc)
#  write_imd_header3(output12,imd_box,atype,rc)
#  write_imd_header3(output21,imd_box,atype,rc)
#  write_imd_header3(output22,imd_box,atype,rc)
#  write_imd_header3(output31,imd_box,atype,rc)
#  write_imd_header3(output32,imd_box,atype,rc)
#  write_imd_header3(output33,imd_box,atype,rc)
#  write_imd_header3(output41,imd_box,atype,rc)
#  write_imd_header3(output42,imd_box,atype,rc)
#  write_imd_header3(output51,imd_box,atype,rc)
#  write_imd_header3(output52,imd_box,atype,rc)
#  write_imd_header3(output61,imd_box,atype,rc)
#  write_imd_header3(output62,imd_box,atype,rc)
#  output1.write("# non-bridging_oxygen information \n")
#  output1.write("# Si	Al	Ca	Mg	Number_Of_Units \n")
#  output2.write("# bridging_oxygen information \n")
#  output2.write("# Si   Al      Ca      Mg      Number_Of_Units \n")
#  for i in range(TotalAtoms):
#    #if atom_type[i] == 1 and nnl_count[i] != 2 :    
#    if atom_type[i] == 1 :
#      neighbour_atom_types=[]
#      for j in nnl[i]:
#        neighbour_atom_types.append(atom_type[j])
#      #if (4 in set(neighbour_atom_types)) or (5 in set(neighbour_atom_types)):
#      #  output1.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
#      #else:
#      #  output2.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
#      number_of_Si=neighbour_atom_types.count(2)
#      number_of_Al=neighbour_atom_types.count(3)
#      number_of_Ca=neighbour_atom_types.count(4)
#      number_of_Mg=neighbour_atom_types.count(5)
#      if number_of_Si+number_of_Al < 2:
#        output11.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
#        output12.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
#        NBOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]=NBOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]+1
#        NBOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]=NBOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]+1
#        for j in nnl[i]:
#          output11.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
#          output12.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
#          if atom_type[j] == 2 or atom_type[j] == 3:
#            for k in nnl[j]:
#              if k != i :
#                output12.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
#        NBOS[i]=1
#        OCa_count[i]=number_of_Ca
#        OMg_count[i]=number_of_Mg
#        if number_of_Si == 1:
#          NBOCa_Si_count=NBOCa_Si_count+number_of_Ca
#          NBOMg_Si_count=NBOMg_Si_count+number_of_Mg
#          NBOSi_count=NBOSi_count+1
#        elif number_of_Al == 1:
#          NBOCa_Al_count=NBOCa_Al_count+number_of_Ca
#          NBOMg_Al_count=NBOMg_Al_count+number_of_Mg
#          NBOAl_count=NBOAl_count+1
#      elif number_of_Si == 2:
#        output21.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
#        output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
#        output31.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
#        output32.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
#        BOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]=BOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]+1
#        BOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]=BOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]+1
#        for j in nnl[i]:
#          output21.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
#          output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
#          output31.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
#          output32.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
#          if atom_type[j] == 2:
#            for k in nnl[j]:
#              if k != i :
#                output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
#                output32.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
#                if number_of_Al == 0 and number_of_Ca == 0 and number_of_Mg == 0:
#                  output33.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
#          if number_of_Al == 0 and number_of_Ca == 0 and number_of_Mg == 0:
#            output33.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
#        if number_of_Al == 0 and number_of_Ca == 0 and number_of_Mg == 0:
#          output33.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
#        OCa_count[i]=number_of_Ca
#        OMg_count[i]=number_of_Mg
#        BOCa_2Si_count=BOCa_2Si_count+number_of_Ca
#        BOMg_2Si_count=BOMg_2Si_count+number_of_Mg
#        BOSiOSi_count=BOSiOSi_count+1
#      elif number_of_Al == 2:
#        output21.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
#        output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
#        output41.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
#        output42.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
#        BOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]=BOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]+1
#        BOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]=BOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]+1
#        for j in nnl[i]:
#          output21.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
#          output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
#          output41.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
#          output42.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
#          if atom_type[j] == 3:
#            for k in nnl[j]:
#              if k != i :
#                output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
#                output42.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
#        OCa_count[i]=number_of_Ca
#        OMg_count[i]=number_of_Mg
#        BOCa_2Al_count=BOCa_2Al_count+number_of_Ca
#        BOMg_2Al_count=BOMg_2Al_count+number_of_Mg
#        BOAlOAl_count=BOAlOAl_count+1
#      elif number_of_Si == 1 and number_of_Al == 1:
#        output21.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
#        output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
#        output51.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
#        output52.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[i],atom_type[i],atom_posx[i],atom_posy[i],atom_posz[i],nnl_count[i]))
#        BOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]=BOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]+1
#        BOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]=BOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]+1
#        for j in nnl[i]:
#          output21.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
#          output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
#          output51.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
#          output52.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[j],atom_type[j],atom_posx[j],atom_posy[j],atom_posz[j],nnl_count[j]))
#          if atom_type[j] == 2 or atom_type[j] == 3:
#            for k in nnl[j]:
#              if k != i :
#                output22.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
#                output52.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
#        OCa_count[i]=number_of_Ca
#        OMg_count[i]=number_of_Mg
#        BOCa_SiAl_count = BOCa_SiAl_count + number_of_Ca
#        BOMg_SiAl_count = BOMg_SiAl_count + number_of_Mg
#        BOSiOAl_count=BOSiOAl_count+1
#      elif number_of_Si+number_of_Al > 2:
#        NBO3Si_count=NBO3Si_count+1
#        #if number_of_Al == 3 or number_of_Al=2 and :  
#        NBOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]=NBOI[number_of_Si][number_of_Al][number_of_Ca][number_of_Mg]+1
#        NBOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]=NBOI1[number_of_Si][number_of_Al][number_of_Ca+number_of_Mg]+1
 #       NBOS[i]=1
#
#  for i in range(8):
#    for j in range(8):
#      for k in range(8):
#        for l in range(8):
#          if NBOI[i][j][k][l] != 0:
#            output1.write("%d Si + %d Al + %d Ca + %d Mg = %d \n" %(i,j,k,l,NBOI[i][j][k][l]))
#          if BOI[i][j][k][l] != 0:
#            output2.write("%d Si + %d Al + %d Ca + %d Mg = %d \n" %(i,j,k,l,BOI[i][j][k][l]))
#  # computing total Ca and Mg connected to Non-Bridiging oxygen and Bridging Oxygen
#  sum1=0
#  sum2=0
#  sum3=0
#  sum4=0
#  sum11=0
#  sum31=0
#  for i in range(TotalAtoms):
#    # for Non bridging oxygen
#    if atom_type[i] == 1 and NBOS[i]==1 :
#      sum1=sum1+OCa_count[i]
#      sum2=sum2+OMg_count[i]
#      sum11=sum11+1
#    # For bridging oxygen
#    elif atom_type[i] == 1 and NBOS[i]==0 :
#      sum3=sum3+OCa_count[i]
#      sum4=sum4+OMg_count[i]
#      sum31=sum31+1
#
#  output1.write("\n\n\n\n***************************************\n")
#  output1.write("# Si	Al 	(R=either Ca or Mg)	Number_of_units\n")
#  output2.write("\n\n\n\n***************************************\n")
#  output2.write("# Si   Al      (R=either Ca or Mg)     Number_of_units\n")
#  for i in range(8):
#    for j in range(8):
#      for k in range(8):
#        if NBOI1[i][j][k] != 0:
#          output1.write("%d Si + %d Al + %d (R) = %d \n" %(i,j,k,NBOI1[i][j][k]))
#        if BOI1[i][j][k] != 0:
#          output2.write("%d Si + %d Al + %d (R) = %d \n" %(i,j,k,BOI1[i][j][k]))
#
#  output1.write("\n\nTotal number of Non-Bridging oxygen atoms 		= %d \n" %(sum11))
#  output1.write("Number of Ca Connected to Non-Bridging oxygen  	= %d \n" %(sum1))
#  output1.write("Number of Mg Connected to Non-Bridging oxygen		= %d \n" %(sum2))
#  output1.write("Number of Ca Connected to Non-Bridging oxygen (Si) 	= %d \n" %(NBOCa_Si_count))
#  output1.write("Number of Ca Connected to Non-Bridging oxygen (Al) 	= %d \n" %(NBOCa_Al_count))
#  output1.write("Number of Mg Connected to Non-Bridging oxygen (Si) 	= %d \n" %(NBOMg_Si_count))
#  output1.write("Number of Mg Connected to Non-Bridging oxygen (Al) 	= %d \n" %(NBOMg_Al_count))
#  output1.write("\n\nNon-Bridging Oxygen SiO 				= %d \n" %(NBOSi_count))
#  output1.write("Non-Bridging Oxygen AlO                            	= %d \n" %(NBOAl_count))
#  output1.write("Non-Bridging Oxygen 3(Si+Al)                           = %d \n" %(NBO3Si_count))

#  output2.write("\n\nTotal number of Bridging oxygen atoms		= %d \n" %(sum31))
#  output2.write("Number of Ca Connected to Bridging oxygen		= %d \n" %(sum3))
#  output2.write("Number of Mg Connected to Bridging oxygen 		= %d \n" %(sum4))
#  output2.write("Number of Ca Connected to Bridging oxygen (2Si)	= %d \n" %(BOCa_2Si_count))
#  output2.write("Number of Ca Connected to Bridging oxygen (1Si+1Al)	= %d \n" %(BOCa_SiAl_count))
#  output2.write("Number of Ca Connected to Bridging oxygen (2Al)	= %d \n" %(BOCa_2Al_count))
#  output2.write("Number of Mg Connected to Bridging oxygen (2Si)	= %d \n" %(BOMg_2Si_count))
#  output2.write("Number of Mg Connected to Bridging oxygen (1Si+1Al)	= %d \n" %(BOMg_SiAl_count))
#  output2.write("Number of Mg Connected to Bridging oxygen (2Al) 	= %d \n" %(BOMg_2Al_count))
#  output2.write("\n\nBridging Oxygen SiOSi                            	= %d \n" %(BOSiOSi_count))
#  output2.write("Bridging Oxygen SiOAl					= %d \n" %(BOSiOAl_count))
#  output2.write("Bridging Oxygen AlOAl					= %d \n" %(BOAlOAl_count))

  """
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

  """
  # **************************************************************************************
  # finding QAl(mAlPlusP)(nP) (Q0, Q1, Q2, Q3, Q4)
  QAl=[[0 for i in xrange(MAX_ATOM_TYPES+1)] for j in xrange(MAX_ATOM_TYPES+1)]						# QAl[Al+P][P]=Total units
  QAl_atoms=[[[-1 for k in xrange(50000)] for j in xrange(MAX_ATOM_TYPES+1)] for i in xrange(MAX_ATOM_TYPES+1)]		# QAl_atoms[Al+P][P][count]=atom_id
  QAl_Pstatus=[-1 for i in xrange(TotalAtoms)]
  QAl_Alstatus=[-1 for i in xrange(TotalAtoms)]
  QAl_P_list=[-1 for i in xrange(TotalAtoms)]
  SelectAtoms1=[-1 for i in xrange(TotalAtoms)]
  SelectAtoms2=[-1 for i in xrange(TotalAtoms)]
  SelectAtoms3=[-1 for i in xrange(TotalAtoms)]
  SelectAtoms4=[-1 for i in xrange(TotalAtoms)]

  output1=open(BASE_FILE+str("_QAl_info.data"), 'w')
  output2=open(BASE_FILE+str("_QAl_additional_info.data"), 'w')
  output1.write("# QAl information\n")
  output2.write("# Special Al which has 4 neighboring oxygen atoms but the neighbouring oxygen atoms has bonded with more than 4 cations\n")
  # looping all atoms
  for C1 in range(TotalAtoms):
    # go only Al atoms with 4 coordination
    if atom_type[C1] == given_Al_type and nnl_count[C1] == 4:
      count_Al=0
      count_P=0
      # first neighbours of Al (i.e O)
      for A1 in nnl[C1]:
        # second neighbours of Al (i.e cations)
        for C2 in nnl[A1]:
          # go only to other cations and also which is 4 coordinated
          if C1 != C2 and nnl_count[C2] == 4:
              # increase the count if the other cation is Si
              if atom_type[C2] == given_Al_type:
                count_Al=count_Al+1
              # increase the count if the other cation is Al
              elif atom_type[C2] == given_P_type:
                count_P=count_P+1
                QAl_P_list[C2]=1
                SelectAtoms3[C2]=1
                SelectAtoms4[C2]=1
                for O2_P in nnl[C2]:
                  SelectAtoms3[O2_P] = 1
                  SelectAtoms4[O2_P] = 1

      if count_Al < 5:
        QAl_atoms[count_Al][count_P][QAl[count_Al][count_P]]=C1
        QAl_Alstatus[C1]=count_Al
        QAl_Pstatus[C1]=count_P
        QAl[count_Al][count_P] = QAl[count_Al][count_P]+1
      else:
        output2.write("%d \n" %(C1))

  for i in range(5):
    for j in range(5):
        #if j<=i:
        if QAl[i][j] != 0:
          output1.write("QAl%d (%dP) = %d \n" %(i,j,QAl[i][j]))
        if j == 0:
          for k in QAl_atoms[i][j]:
            if k != -1:
              SelectAtoms1[k]=1
              for l in nnl[k]:
                SelectAtoms1[l]=1
        else:
          for k in QAl_atoms[i][j]:
            if k != -1:
              SelectAtoms2[k]=1
              SelectAtoms4[k]=1
              for l in nnl[k]:
                SelectAtoms2[l]=1
                SelectAtoms4[l]=1

  output3=open(BASE_FILE+str("_QAl0P_Al_tets.atoms"),'w')
  output4=open(BASE_FILE+str("_QAlmP_Al_tets.atoms"),'w')
  output5=open(BASE_FILE+str("_QAlmP_P_tets.atoms"),'w')
  output6=open(BASE_FILE+str("_QAlmP_Al_and_P_tets.atoms"),'w')
  write_imd_header4(output3,imd_box,atype,rc,1)
  write_imd_header4(output4,imd_box,atype,rc,1)
  write_imd_header4(output5,imd_box,atype,rc,1)
  write_imd_header4(output6,imd_box,atype,rc,1)
  for k in range(0,TotalAtoms):
    if SelectAtoms1[k] == 1:
      #print "%d \n" %(k)
      output3.write("%d  %d  %lf  %lf  %lf %d %d %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k],QAl_Alstatus[k],QAl_Pstatus[k]))
  for k in range(TotalAtoms):
    if SelectAtoms2[k] == 1:
      output4.write("%d  %d  %lf  %lf  %lf %d %d %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k],QAl_Alstatus[k],QAl_Pstatus[k]))
  for k in range(TotalAtoms):
    if SelectAtoms3[k] == 1:
      output5.write("%d  %d  %lf  %lf  %lf %d %d %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k],QAl_Alstatus[k],QAl_Pstatus[k]))
  for k in range(TotalAtoms):
    if SelectAtoms4[k] == 1:
      output6.write("%d  %d  %lf  %lf  %lf %d %d %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k],QAl_Alstatus[k],QAl_Pstatus[k]))

  output1.close()
  output2.close()
  output3.close()
  output4.close()
  output5.close()
  output6.close()

  # **************************************************************************************
  # finding QP(mAlPlusP)(nAl) (Q0, Q1, Q2, Q3, Q4)
  QP=[[0 for i in xrange(MAX_ATOM_TYPES+1)] for j in xrange(MAX_ATOM_TYPES+1)]                                         # QP[Al+P][Al]=Total units
  QP_atoms=[[[-1 for k in xrange(50000)] for j in xrange(MAX_ATOM_TYPES+1)] for i in xrange(MAX_ATOM_TYPES+1)]          # QP_atoms[P+Al][Al][count]=atom_id
  QP_Alstatus=[-1 for i in xrange(TotalAtoms)]
  QP_Pstatus=[-1 for i in xrange(TotalAtoms)]
  QP_Al_list=[-1 for i in xrange(TotalAtoms)]
  SelectAtoms1=[-1 for i in xrange(TotalAtoms)]
  SelectAtoms2=[-1 for i in xrange(TotalAtoms)]
  SelectAtoms3=[-1 for i in xrange(TotalAtoms)]
  SelectAtoms4=[-1 for i in xrange(TotalAtoms)]

  output1=open(BASE_FILE+str("_QP_info.data"), 'w')
  output2=open(BASE_FILE+str("_QP_additional_info.data"), 'w')
  output1.write("# QP information\n")
  output2.write("# Special P which has 4 neighboring oxygen atoms but the neighbouring oxygen atoms has bonded with more than 4 cations\n")
  # looping all atoms
  for C1 in range(TotalAtoms):
    # go only P atoms with 4 coordination
    if atom_type[C1] == given_P_type and nnl_count[C1] == 4:
      count_P=0
      count_Al=0
      # first neighbours of P (i.e O)
      for A1 in nnl[C1]:
        # second neighbours of P (i.e cations)
        for C2 in nnl[A1]:
          # go only to other cations and also which is 4 coordinated
          if C1 != C2:
              # increase the count if the other cation is Si
              if atom_type[C2] == given_P_type and nnl_count[C2] == 4:
                count_P=count_P+1
              # increase the count if the other cation is Al
              #elif atom_type[C2] == given_Al_type and nnl_count[C2] == 4:
              elif atom_type[C2] == given_Al_type:
                count_Al=count_Al+1
                QP_Al_list[C2]=1
                SelectAtoms3[C2]=1
                SelectAtoms4[C2]=1
                for O2_Al in nnl[C2]:
                  SelectAtoms3[O2_Al] = 1
                  SelectAtoms4[O2_Al] = 1

      if count_P + count_Al < 5:
        QP_atoms[count_P][count_Al][QP[count_P][count_Al]]=C1
        QP_Pstatus[C1]=count_P
        QP_Alstatus[C1]=count_Al
        QP[count_P][count_Al] = QP[count_P][count_Al]+1
      else:
        output2.write("%d \n" %(C1))

  for i in range(5):
    for j in range(5):
        #if j<=i:
        if QP[i][j] != 0:
          output1.write("QP%d (%dAl) = %d \n" %(i,j,QP[i][j]))
        if j == 0:
          for k in QP_atoms[i][j]:
            if k != -1:
              SelectAtoms1[k]=1
              for l in nnl[k]:
                SelectAtoms1[l]=1
        else:
          for k in QP_atoms[i][j]:
            if k != -1:
              SelectAtoms2[k]=1
              SelectAtoms4[k]=1
              for l in nnl[k]:
                SelectAtoms2[l]=1
                SelectAtoms4[l]=1

  output3=open(BASE_FILE+str("_QP0Al_P_tets.atoms"),'w')
  output4=open(BASE_FILE+str("_QPmAl_P_tets.atoms"),'w')
  output5=open(BASE_FILE+str("_QPmAl_Al_tets.atoms"),'w')
  output6=open(BASE_FILE+str("_QPmAl_P_and_Al_tets.atoms"),'w')

  write_imd_header4(output3,imd_box,atype,rc,2)
  write_imd_header4(output4,imd_box,atype,rc,2)
  write_imd_header4(output5,imd_box,atype,rc,2)
  write_imd_header4(output6,imd_box,atype,rc,2)
  for k in range(TotalAtoms):
    if SelectAtoms1[k] == 1:
      output3.write("%d  %d  %lf  %lf  %lf %d %d %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k],QP_Pstatus[k],QP_Alstatus[k]))
  for k in range(TotalAtoms):
    if SelectAtoms2[k] == 1:
      output4.write("%d  %d  %lf  %lf  %lf %d %d %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k],QP_Pstatus[k],QP_Alstatus[k]))
  for k in range(TotalAtoms):
    if SelectAtoms3[k] == 1:
      output5.write("%d  %d  %lf  %lf  %lf %d %d %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k],QP_Pstatus[k],QP_Alstatus[k]))
  for k in range(TotalAtoms):
    if SelectAtoms4[k] == 1:
      output6.write("%d  %d  %lf  %lf  %lf %d %d %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k],QP_Pstatus[k],QP_Alstatus[k]))

  output1.close()
  output2.close()
  output3.close()
  output4.close()
  output5.close()
  output6.close()
  # **************************************************************************************
  # finding Al5-O-X
  QAl5_AllNeighbours=[[0 for i in xrange(MAX_NEIGHBOURS)] for j in xrange(MAX_ATOM_TYPES+1)]
  QAl5=[[0 for i in xrange(6)] for j in xrange(6)]
  O3_count=0
  for C1 in range(TotalAtoms):
    if atom_type[C1] == given_Al_type and nnl_count[C1] == 5:
      QAl5_Al_count=0
      QAl5_P_count=0
      for A1 in nnl[C1]:
        for C2 in nnl[A1]:
          if C1 != C2 and nnl_count[C2] == 4:
            if atom_type[C2] == given_Al_type:
              QAl5_Al_count=QAl5_Al_count+1
            elif atom_type[C2] == given_P_type:
              QAl5_P_count=QAl5_P_count+1
          if C1 != C2 and (atom_type[C2] == given_Al_type or atom_type[C2] == given_P_type):
            QAl5_AllNeighbours[atom_type[C2]][nnl_count[C2]]=QAl5_AllNeighbours[atom_type[C2]][nnl_count[C2]]+1
      if QAl5_Al_count + QAl5_P_count < 6:
        QAl5[QAl5_Al_count][QAl5_P_count]=QAl5[QAl5_Al_count][QAl5_P_count]+1
      else:
        O3_count=O3_count+1

  output1=open(BASE_FILE+str("_Al5_info.data"),'w')
  output1.write("# Al5 with m Al4 and n P4\n")
  for i in range(5):
    for j in range(5):
      if QAl5[i][j] != 0:
        output1.write("QAl5 (%d Al4) (%d P4)  = %d\n" %(i,j,QAl5[i][j]))
  output1.write("Total O3 count associated with Al5 = %d (However, it needs to be checked)\n" %(O3_count))
  output1.write("\n\n# Al5 with all possible coordinations\n")
  for i in range(MAX_ATOM_TYPES+1):
    for j in range(MAX_NEIGHBOURS):
      if QAl5_AllNeighbours[i][j] != 0:
        output1.write("Al5 -O- %s%d  = %d\n" %(str(atype[i]),j,QAl5_AllNeighbours[i][j]))
  output1.close()

  # **************************************************************************************
  # finding Al(m)-O-Al(n) bonds
  AlOAl_atom1_list=[[[-1 for k in xrange(TotalAtoms)] for j in xrange(MAX_NEIGHBOURS)] for i in xrange(MAX_NEIGHBOURS)]
  AlOAl_atom2_list=[[[-1 for k in xrange(TotalAtoms)] for j in xrange(MAX_NEIGHBOURS)] for i in xrange(MAX_NEIGHBOURS)]
  AlOAl=[[ 0 for j in xrange(MAX_NEIGHBOURS)] for i in xrange(MAX_NEIGHBOURS)]
  SelectAtoms1=[[[-1 for k in xrange(TotalAtoms)] for j in xrange(MAX_NEIGHBOURS)] for i in xrange(MAX_NEIGHBOURS)]	# collect Al(m)-O-Al(n) no tetra hedra atoms 
  SelectAtoms2=[-1 for i in xrange(TotalAtoms)]		# Collect Al-O-Al irrespective of state of Al
  SelectAtoms3=[-1 for i in xrange(TotalAtoms)]		# Select all tetrahedra atoms which are invovled in Al(4)-O-Al(4)
  SelectAtoms4=[-1 for i in xrange(TotalAtoms)]		# Select all tetrahedra atoms which are invovled in Al(m)-O-Al(n)
  SelectAtoms5=[[[-1 for k in xrange(TotalAtoms)] for j in xrange(MAX_NEIGHBOURS)] for i in xrange(MAX_NEIGHBOURS)]	# collect Al,O and Na atoms where Na is attached to Al(m)-O-Al(n)
  SelectAtoms6=[[[-1 for k in xrange(TotalAtoms)] for j in xrange(MAX_NEIGHBOURS)] for i in xrange(MAX_NEIGHBOURS)]	# collect Al,O and Na atoms where Na is attached to Al(m)-O-Al(n) + Al tetra
  SelectAtoms7=[-1 for i in xrange(TotalAtoms)]		# collect Al,O and Na atoms where Na is attached to Al-O-Al irrespective of Al state		
  SelectAtoms8=[-1 for i in xrange(TotalAtoms)]		# collect Al,O and Na atoms where Na is attached to Al-O-Al irrespective of Al state and also Al tetras

  for C1 in range(TotalAtoms):
    if atom_type[C1] == given_Al_type:
      C1_coord=nnl_count[C1]
      for A1 in nnl[C1]:
        for C2 in nnl[A1]:
          if C1 != C2 and atom_type[C2] == given_Al_type:
            C2_coord=nnl_count[C2]
            SelectAtoms2[C1]=1
            SelectAtoms2[C2]=1
            SelectAtoms2[A1]=1
            SelectAtoms4[C1]=1
            SelectAtoms4[C2]=1
            for A2 in nnl[C1]:
              SelectAtoms4[A2]=1
            for A2 in nnl[C2]:
              SelectAtoms4[A2]=1

            if C1_coord <= C2_coord:
              AlOAl_atom1_list[C1_coord][C2_coord][AlOAl[C1_coord][C2_coord]] = C1
              AlOAl_atom2_list[C1_coord][C2_coord][AlOAl[C1_coord][C2_coord]] = C2
              AlOAl[C1_coord][C2_coord]=AlOAl[C1_coord][C2_coord]+1
              SelectAtoms1[C1_coord][C2_coord][C1]=1
              SelectAtoms1[C1_coord][C2_coord][C2]=1
              SelectAtoms1[C1_coord][C2_coord][A1]=1
              Selected_Na_Atom_Numbers=[]
              for C22 in nnl[A1]:
                if atom_type[C22] == given_Na_type:
                  Selected_Na_Atom_Numbers.append(C22)
              if len(Selected_Na_Atom_Numbers) != 0:
                SelectAtoms5[C1_coord][C2_coord][C1]=1
                SelectAtoms5[C1_coord][C2_coord][C2]=1
                SelectAtoms5[C1_coord][C2_coord][A1]=1
                SelectAtoms6[C1_coord][C2_coord][C1]=1
                SelectAtoms6[C1_coord][C2_coord][C2]=1
                SelectAtoms6[C1_coord][C2_coord][A1]=1
                for A2 in nnl[C1]:
                  SelectAtoms6[C1_coord][C2_coord][A2]=1
                for A2 in nnl[C2]:
                  SelectAtoms6[C1_coord][C2_coord][A2]=1
                for CNa in Selected_Na_Atom_Numbers:
                  SelectAtoms5[C1_coord][C2_coord][CNa]=1
                  SelectAtoms6[C1_coord][C2_coord][CNa]=1
                  
            elif C1_coord > C2_coord:
              AlOAl_atom1_list[C2_coord][C1_coord][AlOAl[C2_coord][C1_coord]] = C2
              AlOAl_atom2_list[C2_coord][C1_coord][AlOAl[C2_coord][C1_coord]] = C1
              AlOAl[C2_coord][C1_coord]=AlOAl[C2_coord][C1_coord]+1
              SelectAtoms1[C2_coord][C1_coord][C2]=1
              SelectAtoms1[C2_coord][C1_coord][C1]=1
              SelectAtoms1[C2_coord][C1_coord][A1]=1
              Selected_Na_Atom_Numbers=[]
              for C22 in nnl[A1]:
                if atom_type[C22] == given_Na_type:
                  Selected_Na_Atom_Numbers.append(C22)
              if len(Selected_Na_Atom_Numbers) != 0:
                SelectAtoms5[C2_coord][C1_coord][C1]=1
                SelectAtoms5[C2_coord][C1_coord][C2]=1
                SelectAtoms5[C2_coord][C1_coord][A1]=1
                SelectAtoms6[C2_coord][C1_coord][C1]=1
                SelectAtoms6[C2_coord][C1_coord][C2]=1
                SelectAtoms6[C2_coord][C1_coord][A1]=1
                for A2 in nnl[C1]:
                  SelectAtoms6[C2_coord][C1_coord][A2]=1
                for A2 in nnl[C2]:
                  SelectAtoms6[C2_coord][C1_coord][A2]=1
                for CNa in Selected_Na_Atom_Numbers:
                  SelectAtoms5[C2_coord][C1_coord][CNa]=1
                  SelectAtoms6[C2_coord][C1_coord][CNa]=1
  # writing Al(m)-O-Al(n)
  for i in range(MAX_NEIGHBOURS):
    for j in range(MAX_NEIGHBOURS):
      if i<=j and AlOAl[i][j] != 0:
        output1=open(BASE_FILE+str("_Al")+str(i)+str("Al")+str(j)+str(".imd"),'w')
        write_imd_header4(output1,imd_box,atype,rc,0)
        for k in range(TotalAtoms):
          if SelectAtoms1[i][j][k] == 1:
            output1.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
        output1.close()

  # output1 = writing out Al(m)-O-Al(n) and also the Na atoms in these network 
  # output2 = writing out Al(m)-O-Al(n) and also the Na atoms in these network and also Al tetrahedra included
  for i in range(MAX_NEIGHBOURS):
    for j in range(MAX_NEIGHBOURS):
      if i<=j and AlOAl[i][j] != 0:
        output1=open(BASE_FILE+str("_Na_atoms_in_Al")+str(i)+str("Al")+str(j)+str("_network.imd"),'w')
        output2=open(BASE_FILE+str("_Na_atoms_in_Al")+str(i)+str("Al")+str(j)+str("_network_and_Al_tets.imd"),'w')
        write_imd_header4(output1,imd_box,atype,rc,0)
        write_imd_header4(output2,imd_box,atype,rc,0)
        for k in range(TotalAtoms):
          if SelectAtoms5[i][j][k] == 1:
            output1.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
            SelectAtoms7[k] =1
          if SelectAtoms6[i][j][k] == 1:
            output2.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
            SelectAtoms8[k] =1
        output1.close()
        output2.close()
  #output3 = writing out Al-O-Al irrespective of Al state 
  #output4 = writing out Al-O-Al irrespective of Al state and also the Na atoms in these network and also Al tetrahedra included
  output3=open(BASE_FILE+str("_Na_atoms_in_AlOAl_network.imd"),'w')
  output4=open(BASE_FILE+str("_Na_atoms_in_AlOAl_network_and_Al_tetra.imd"),'w')
  write_imd_header4(output3,imd_box,atype,rc,0)
  write_imd_header4(output4,imd_box,atype,rc,0)
  for k in range(TotalAtoms):
    if SelectAtoms7[k] == 1:
       output3.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
    if SelectAtoms8[k] == 1:
       output4.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
  output3.close()
  output4.close()

  output1=open(BASE_FILE+str("_AlmOAln.imd"),'w')
  write_imd_header4(output1,imd_box,atype,rc,0)
  for k in range(TotalAtoms):
    if SelectAtoms2[k] == 1:
      output1.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
  output1.close()

  for k in range(TotalAtoms):
    if SelectAtoms1[4][4][k] == 1 and atom_type[k] == given_Al_type:
      SelectAtoms3[k] = 1
      for l in nnl[k]:
        SelectAtoms3[l] = 1

  output1=open(BASE_FILE+str("_Al4OAl4_tets.imd"),'w')
  output2=open(BASE_FILE+str("_AlmOAln_tets.imd"),'w')
  write_imd_header4(output1,imd_box,atype,rc,0)
  write_imd_header4(output2,imd_box,atype,rc,0)
  for k in range(TotalAtoms):
    if SelectAtoms3[k] == 1:
      output1.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
    if SelectAtoms4[k] == 1:
      output2.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[k],atom_type[k],atom_posx[k],atom_posy[k],atom_posz[k],nnl_count[k]))
  output1.close()

  # **************************************************************************************
  # finding Na neighbours information
  NaOAl=[[ 0 for j in xrange(MAX_NEIGHBOURS)] for i in xrange(MAX_NEIGHBOURS)]
  for C1 in range(TotalAtoms):
    if atom_type[C1] == given_Na_type:
      NaOAl_coord_count=[0 for i in xrange(MAX_NEIGHBOURS)]
      for A1 in nnl[C1]:
        for C2 in nnl[A1]:
          if C1 != C2 and atom_type[C2] == given_Al_type:
            NaOAl_coord_count[nnl_count[C2]]=NaOAl_coord_count[nnl_count[C2]]+1
      for i in range(MAX_NEIGHBOURS):
        NaOAl[nnl_count[C1]][i]=NaOAl[nnl_count[C1]][i]+NaOAl_coord_count[i]

  output1=open(BASE_FILE+str("_NaOAl_neighbours_info.data"),'w')
  output1.write("# Na_coordination \t Al_coordination \t TotalNumberOf_Na_units")
  for i in range(MAX_NEIGHBOURS):
    for j in range(MAX_NEIGHBOURS):
      if NaOAl[i][j] != 0:
        output1.write("%d  %d  %d \n" %(i,j,NaOAl[i][j]))
  output1.close()

  # **************************************************************************************
  # finding O neighbours information
  O_coord2=[[[ 0 for k in xrange(MAX_NEIGHBOURS)] for j in xrange(MAX_ATOM_TYPES+1)] for k in xrange(MAX_NEIGHBOURS)]
  for A1 in range(TotalAtoms):
    if atom_type[A1] == given_O_type:
      O_coord1=[[ 0 for k in xrange(MAX_NEIGHBOURS)] for j in xrange(MAX_ATOM_TYPES+1)]
      for C1 in nnl[A1]:
        #if atom_type[C1] == given_Al_type:
        #  OAl_coord[nnl_count[C1]] = OAl_coord[nnl_count[C1]]+1
        #elif atom_type[C1] == given_Na_type:
        #  ONa_coord[nnl_count[C1]] = ONa_coord[nnl_count[C1]]+1
        #elif atom_type[C1] == given_P_type:
        #  OP_coord[nnl_count[C1]] = OP_coord[nnl_count[C1]]+1
        O_coord1[atom_type[C1]][nnl_count[C1]]=O_coord1[atom_type[C1]][nnl_count[C1]]+1
      for i in range(MAX_ATOM_TYPES+1):
        for j in range(MAX_NEIGHBOURS):
          O_coord2[nnl_count[A1]][i][j]=O_coord2[nnl_count[A1]][i][j]+O_coord1[i][j]

  output1=open(BASE_FILE+str("_O_neighbours_info.data"),'w')
  output1.write("# O_coord \t C1_atom_type \t Number_of_C1_atoms \t Total_Such_O_units\n")
  for i in range(MAX_NEIGHBOURS):
    for j in range(MAX_ATOM_TYPES+1):
      for k in range(MAX_NEIGHBOURS):
        if O_coord2[i][j][k] != 0:
          output1.write("%d  %d  %d  %d \n" %(i,j,k,O_coord2[i][j][k]))
  output1.close()

  # **************************************************************************************
  # finding Na information within Al(m)-O-Al(n)
  # AlOAl_count1[Number_of_Al_4_fold][Number_of_Al_5_fold][Number_of_Al_6_fold][Number_of_Na]
  AlOAl_count1=[[[[[0 for i1 in xrange(8)] for i2 in xrange(8)] for i3 in xrange(8)] for i4 in xrange(8)] for i5 in xrange(8)]
  for A1 in range(TotalAtoms):
    if atom_type[A1] == given_O_type:
      count_Al=[0 for i in xrange(MAX_ATOM_TYPES+1)]
      count_P=0
      count_Na=0
      for C1 in nnl[A1]:
        if atom_type[C1] == given_Al_type:
          count_Al[nnl_count[C1]]=count_Al[nnl_count[C1]]+1
        #elif atom_type[C1] == given_P_type:
        #  count_P[nnl_count[C1]]=count_P[nnl_count[C1]]+1
        elif atom_type[C1] == given_P_type:
          count_P=count_P+1
        elif atom_type[C1] == given_Na_type:
          count_Na=count_Na+1
      AlOAl_count1[count_Al[4]][count_Al[5]][count_Al[6]][count_P][count_Na]=AlOAl_count1[count_Al[4]][count_Al[5]][count_P][count_Al[6]][count_Na]+1

  output1=open(BASE_FILE+str("_Na_info_in_AlOAl_And_AlOP.data"),'w')
  output1.write("# NumberOf_Al  (Al4  Al5  Al6)  Number_of_P  Number_of_Na  Number_of_oxygen\n")
  for i1 in range(8):
    for i2 in range(8):
      for i3 in range(8):
        for i4 in range(8):
          for i5 in xrange(8):
            if AlOAl_count1[i1][i2][i3][i4][i5] != 0 :
              output1.write("%d ( %d  %d  %d) (Al) \t + %d (P) \t + %d (Na) \t = %d \n" %(i1+i2+i3,i1,i2,i3,i4,i5,AlOAl_count1[i1][i2][i3][i4][i5]))
  output1.close()

  """
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
  """
  # **************************************************************************************
  # coordination of each atom
  Coord=[[0 for j in range(MAX_NEIGHBOURS+1)] for i in range(MAX_ATOM_TYPES+1)]		# Coord[atom_type][coordination]
  for i in range(TotalAtoms):
    Coord[atom_type[i]][nnl_count[i]]=Coord[atom_type[i]][nnl_count[i]]+1
  output_coord1=open(BASE_FILE+str("_average_coord"),'w')

  MAX_LENGTH_OF_COORD=0
  for i in range(1,MAX_ATOM_TYPES+1):
    if i in given_atom_types:
      avg_coord1=0
      avg_coord2=0
      output_coord=open(BASE_FILE+str("_")+str(atype[i])+str("_coord"),'w')
      for j in range(0,MAX_NEIGHBOURS+1):
        output_coord.write("%d %d \n" %(j,Coord[i][j]))
        MAX_LENGTH_OF_COORD=max(MAX_LENGTH_OF_COORD,Coord[i][j])
        avg_coord1=avg_coord1+(j*Coord[i][j])
        avg_coord2=avg_coord2+Coord[i][j]
      output_coord.close()
      output_coord1.write("%s %lf \n" %(atype[i],1.0*avg_coord1/avg_coord2))
  output_coord1.close()

  Coordination=[[[0 for k in range(MAX_LENGTH_OF_COORD)] for j in range(MAX_NEIGHBOURS+1)] for i in range(MAX_ATOM_TYPES+1)] 	# Coordination[atom_type][coordination][count] = atom_id
  Coordination_count=[[0 for j in range(MAX_NEIGHBOURS+1)] for i in range(MAX_ATOM_TYPES+1)]					# Coordination_count[atom_type][coordination]
  for i in range(TotalAtoms):
    Coordination[atom_type[i]][nnl_count[i]][Coordination_count[atom_type[i]][nnl_count[i]]]=i
    Coordination_count[atom_type[i]][nnl_count[i]]=Coordination_count[atom_type[i]][nnl_count[i]]+1

  for i in range(1,MAX_ATOM_TYPES+1):
    if i in given_atom_types:
      for j in range(MAX_NEIGHBOURS+1):
        if Coord[i][j] !=0:
          output=open(BASE_FILE+str("_")+str(atype[i])+str("_")+str(j)+str("_coord.atomsid"),'w')
          output1=open(BASE_FILE+str("_")+str(atype[i])+str("_")+str(j)+str("_coord.imd"),'w')
          write_imd_header4(output1,imd_box,atype,rc,0)
          for k in range(Coordination_count[i][j]):
            l=Coordination[i][j][k]
            output.write("%d  " %(atom_global_id[l]))
            output1.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[l],atom_type[l],atom_posx[l],atom_posy[l],atom_posz[l],nnl_count[l]))
            for m in nnl[l]:
              if m != -1:
                output.write("%d  " %(atom_global_id[m]))
                output1.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_global_id[m],atom_type[m],atom_posx[m],atom_posy[m],atom_posz[m],nnl_count[m]))
            output.write("\n")
          output.close
          output1.close

  """
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
  Ca_NBRS=[[0 for i in xrange(8)] for j in xrange(8)]
  Mg_NBRS=[[0 for i in xrange(8)] for j in xrange(8)]
  Ca_NBRS1=[[[[[0 for i1 in xrange(8)] for i2 in xrange(8)] for i3 in xrange(8)] for i4 in xrange(8)] for i5 in xrange(8)]
  Mg_NBRS1=[[[[[0 for i1 in xrange(8)] for i2 in xrange(8)] for i3 in xrange(8)] for i4 in xrange(8)] for i5 in xrange(8)]
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

  """
  
