#!/usr/bin/python
"""
@Author : Sudheer Ganisetti
@Date   : Sa 1. Sep 23:43:13 CEST 2018
Mean Square Displacement
"""
import numpy as np
import sys
import math as mt
import itertools as it
import subprocess
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if __name__=="__main__":
  
  tot_argv=len(sys.argv)
  if tot_argv != 3:
     subprocess.call("sudheer_banner")
     print "************* S. Ganisetti *************"
     print "Error: usage is wrong"
     print "./this_program  Ref_ChkptFile.dump   Cur_Chkpt.dump"
     print "Compute Mean Square Displacement"
     print "The considered elements and associated atom types are as follows:"
     print "O(1), Si(2), Al(3), Ca(4), Mg(5), P(6), Na(7), K(8) and Sr(9) "
     print "****************************************"
     sys.exit(0)

  # --------------------------------------------------------------------
  # The cutoff radii are hard coaded, please edit according to your need
  MAX_NEIGHBOURS=15
  MAX_ATOM_TYPES=9
  REF_FILE_NAME1=str(sys.argv[1])
  CUR_FILE_NAME1=str(sys.argv[2])
  #print "%s" %(REF_FILE_NAME1)
  ref_atoms_data1=np.loadtxt(REF_FILE_NAME1,dtype='int,int,float,float,float,float,float,float,float',skiprows=9)
  cur_atoms_data1=np.loadtxt(CUR_FILE_NAME1,dtype='int,int,float,float,float,float,float,float,float',skiprows=9)
  cur_atoms_data2=open(CUR_FILE_NAME1,'r')

  ref_TotalAtoms=len(ref_atoms_data1)
  cur_TotalAtoms=len(cur_atoms_data1)
  ref_TotalProperties=len(ref_atoms_data1[0])

  if ref_TotalAtoms == cur_TotalAtoms:
    TotalAtoms=ref_TotalAtoms
  else:
    print "Error: Total number of atoms in reference configuration is different than the current configuration"
    exit(0)

  ref_atoms_data2=[]
  for i in range(ref_TotalAtoms):
    for j in ref_atoms_data1[i]:
      ref_atoms_data2.append(j)
  ref_atoms_data2=np.array(ref_atoms_data2)
  ref_atoms_data3=np.reshape(ref_atoms_data2,(ref_TotalAtoms,ref_TotalProperties))
  MAX_ATOM_NUMBER=int(max(ref_atoms_data3.T[0]))
  
  #atom_status=[0 for i in xrange(TotalAtoms)]
  ref_atom_id		=[-1 for i in xrange(MAX_ATOM_NUMBER+2)]
  ref_atom_type		=[-1 for i in xrange(MAX_ATOM_NUMBER+2)]
  ref_atom_charge	=[0.0 for i in xrange(MAX_ATOM_NUMBER+2)]
  ref_atom_posx		=[-1.0 for i in xrange(MAX_ATOM_NUMBER+2)]
  ref_atom_posy		=[-1.0 for i in xrange(MAX_ATOM_NUMBER+2)]
  ref_atom_posz		=[-1.0 for i in xrange(MAX_ATOM_NUMBER+2)]
  cur_atom_id		=[-1 for i in xrange(MAX_ATOM_NUMBER+2)]
  cur_atom_type		=[-1 for i in xrange(MAX_ATOM_NUMBER+2)]
  cur_atom_charge	=[0.0 for i in xrange(MAX_ATOM_NUMBER+2)]
  cur_atom_posx		=[-1.0 for i in xrange(MAX_ATOM_NUMBER+2)]
  cur_atom_posy		=[-1.0 for i in xrange(MAX_ATOM_NUMBER+2)]
  cur_atom_posz		=[-1.0 for i in xrange(MAX_ATOM_NUMBER+2)]
  dx2			=[0.0 for i in xrange(MAX_ATOM_NUMBER+2)]
  dy2			=[0.0 for i in xrange(MAX_ATOM_NUMBER+2)]
  dz2			=[0.0 for i in xrange(MAX_ATOM_NUMBER+2)]
  dr2			=[0.0 for i in xrange(MAX_ATOM_NUMBER+2)]
  sum_dx2_O  = sum_dy2_O  = sum_dz2_O  = sum_dr2_O  =0.0
  sum_dx2_Si = sum_dy2_Si = sum_dz2_Si = sum_dr2_Si =0.0
  sum_dx2_Al = sum_dy2_Al = sum_dz2_Al = sum_dr2_Al =0.0
  sum_dx2_Ca = sum_dy2_Ca = sum_dz2_Ca = sum_dr2_Ca =0.0
  sum_dx2_Mg = sum_dy2_Mg = sum_dz2_Mg = sum_dr2_Mg =0.0
  sum_dx2_P  = sum_dy2_P  = sum_dz2_P  = sum_dr2_P  =0.0
  sum_dx2_Na = sum_dy2_Na = sum_dz2_Na = sum_dr2_Na =0.0
  sum_dx2_K  = sum_dy2_K  = sum_dz2_K  = sum_dr2_K  =0.0
  sum_dx2_Sr = sum_dy2_Sr = sum_dz2_Sr = sum_dr2_Sr =0.0
  Total_O    = Total_Si   = Total_Al   = Total_Ca   =0
  Total_Mg   = Total_P    = Total_Na   = Total_K    =0
  Total_Sr   = 0

  # Reading box dimensions
  count=1
  for line in cur_atoms_data2:
    data=line.strip().split()
    if count < 10:
      if count == 4:
        if cur_TotalAtoms != int(data[0]):
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

  # Reading reference atoms position
  for i in range(ref_TotalAtoms):
    j			=int(ref_atoms_data3[i][0])
    ref_atom_id[j]	=j
    ref_atom_type[j]	=int(ref_atoms_data3[i][1])
    ref_atom_charge[j]	=ref_atoms_data3[i][2]
    ref_atom_posx[j]	=ref_atoms_data3[i][3]
    ref_atom_posy[j]	=ref_atoms_data3[i][4]
    ref_atom_posz[j]	=ref_atoms_data3[i][5]

  # Reading current atoms position
  for i in range(cur_TotalAtoms):
    j			=int(cur_atoms_data1[i][0])
    cur_atom_id[j]	=j
    cur_atom_type[j]	=int(cur_atoms_data1[i][1])
    cur_atom_charge[j]	=cur_atoms_data1[i][2]
    cur_atom_posx[j]	=cur_atoms_data1[i][3]
    cur_atom_posy[j]	=cur_atoms_data1[i][4]
    cur_atom_posz[j]	=cur_atoms_data1[i][5]

  # compute MSD
  for i in ref_atom_id:
    dx=abs(cur_atom_posx[i]-ref_atom_posx[i])
    if dx > box_xx/2.0:
       dx=dx-box_xx
    dy=abs(cur_atom_posy[i]-ref_atom_posy[i])
    if dy > box_yy/2.0:
       dy=dy-box_yy
    dz=abs(cur_atom_posz[i]-ref_atom_posz[i])
    if dz > box_zz/2.0:
       dz=dz-box_zz
    dx2[i]=dx*dx
    dy2[i]=dy*dy
    dz2[i]=dz*dz
    dr2[i]=dx2[i]+dy2[i]+dz2[i]

  # Sum MSD based on type of atom
  for i in ref_atom_id:
    if ref_atom_type[i] == 1:
       sum_dx2_O = sum_dx2_O + dx2[i]
       sum_dy2_O = sum_dy2_O + dy2[i]
       sum_dz2_O = sum_dz2_O + dz2[i]
       sum_dr2_O = sum_dr2_O + dr2[i]
       Total_O=Total_O+1
    elif ref_atom_type[i] == 2:
       sum_dx2_Si = sum_dx2_Si + dx2[i]
       sum_dy2_Si = sum_dy2_Si + dy2[i]
       sum_dz2_Si = sum_dz2_Si + dz2[i]
       sum_dr2_Si = sum_dr2_Si + dr2[i]
       Total_Si=Total_Si+1
    elif ref_atom_type[i] == 3:
       sum_dx2_Al = sum_dx2_Al + dx2[i]
       sum_dy2_Al = sum_dy2_Al + dy2[i]
       sum_dz2_Al = sum_dz2_Al + dz2[i]
       sum_dr2_Al = sum_dr2_Al + dr2[i]
       Total_Al=Total_Al+1
    elif ref_atom_type[i] == 4:
       sum_dx2_Ca = sum_dx2_Ca + dx2[i]
       sum_dy2_Ca = sum_dy2_Ca + dy2[i]
       sum_dz2_Ca = sum_dz2_Ca + dz2[i]
       sum_dr2_Ca = sum_dr2_Ca + dr2[i]
       Total_Ca=Total_Ca+1
    elif ref_atom_type[i] == 5:
       sum_dx2_Mg = sum_dx2_Mg + dx2[i]
       sum_dy2_Mg = sum_dy2_Mg + dy2[i]
       sum_dz2_Mg = sum_dz2_Mg + dz2[i]
       sum_dr2_Mg = sum_dr2_Mg + dr2[i]
       Total_Mg=Total_Mg+1
    elif ref_atom_type[i] == 6:
       sum_dx2_P = sum_dx2_P + dx2[i]
       sum_dy2_P = sum_dy2_P + dy2[i]
       sum_dz2_P = sum_dz2_P + dz2[i]
       sum_dr2_P = sum_dr2_P + dr2[i]
       Total_P   = Total_P   + 1
    elif ref_atom_type[i] == 7:
       sum_dx2_Na = sum_dx2_Na + dx2[i]
       sum_dy2_Na = sum_dy2_Na + dy2[i]
       sum_dz2_Na = sum_dz2_Na + dz2[i]
       sum_dr2_Na = sum_dr2_Na + dr2[i]
       Total_Na   = Total_Na   + 1
    elif ref_atom_type[i] == 8:
       sum_dx2_K = sum_dx2_K + dx2[i]
       sum_dy2_K = sum_dy2_K + dy2[i]
       sum_dz2_K = sum_dz2_K + dz2[i]
       sum_dr2_K = sum_dr2_K + dr2[i]
       Total_K   = Total_K   + 1
    elif ref_atom_type[i] == 9:
       sum_dx2_Sr = sum_dx2_Sr + dx2[i]
       sum_dy2_Sr = sum_dy2_Sr + dy2[i]
       sum_dz2_Sr = sum_dz2_Sr + dz2[i]
       sum_dr2_Sr = sum_dr2_Sr + dr2[i]
       Total_Sr   = Total_Sr   + 1

  # write output
  if Total_O != 0:
    output_MSD_O  =open("z01_MSD_O.data",'w')
    #output_MSD_O.write("#dx2  dy2  dz2  MSD\n")
    output_MSD_O.write("%lf  %lf  %lf  %lf\n" %(sum_dx2_O/Total_O,sum_dy2_O/Total_O,sum_dz2_O/Total_O,sum_dr2_O/Total_O))
    output_MSD_O.close()
  if Total_Si != 0:
    output_MSD_Si =open("z02_MSD_Si.data",'w')
    #output_MSD_Si.write("#dx2  dy2  dz2  MSD\n")
    output_MSD_Si.write("%lf  %lf  %lf  %lf\n" %(sum_dx2_Si/Total_Si,sum_dy2_Si/Total_Si,sum_dz2_Si/Total_Si,sum_dr2_Si/Total_Si))
    output_MSD_Si.close()
  if Total_Al != 0:
    output_MSD_Al =open("z03_MSD_Al.data",'w')
    #output_MSD_Al.write("#dx2  dy2  dz2  MSD\n")
    output_MSD_Al.write("%lf  %lf  %lf  %lf\n" %(sum_dx2_Al/Total_Al,sum_dy2_Al/Total_Al,sum_dz2_Al/Total_Al,sum_dr2_Al/Total_Al))
    output_MSD_Al.close()
  if Total_Ca != 0:
    output_MSD_Ca =open("z04_MSD_Ca.data",'w')
    #output_MSD_Ca.write("#dx2  dy2  dz2  MSD\n")
    output_MSD_Ca.write("%lf  %lf  %lf  %lf\n" %(sum_dx2_Ca/Total_Ca,sum_dy2_Ca/Total_Ca,sum_dz2_Ca/Total_Ca,sum_dr2_Ca/Total_Ca))
    output_MSD_Ca.close()
  if Total_Mg != 0:
    output_MSD_Mg =open("z05_MSD_Mg.data",'w')
    #output_MSD_Mg.write("#dx2  dy2  dz2  MSD\n")
    output_MSD_Mg.write("%lf  %lf  %lf  %lf\n" %(sum_dx2_Mg/Total_Mg,sum_dy2_Mg/Total_Mg,sum_dz2_Mg/Total_Mg,sum_dr2_Mg/Total_Mg))
    output_MSD_Mg.close()
  if Total_P != 0:
    output_MSD_P =open("z06_MSD_P.data",'w')
    #output_MSD_P.write("#dx2  dy2  dz2  MSD\n")
    output_MSD_P.write("%lf  %lf  %lf  %lf\n" %(sum_dx2_P/Total_P,sum_dy2_P/Total_P,sum_dz2_P/Total_P,sum_dr2_P/Total_P))
    output_MSD_P.close()
  if Total_Na != 0:
    output_MSD_Na =open("z07_MSD_Na.data",'w')
    #output_MSD_Na.write("#dx2  dy2  dz2  MSD\n")
    output_MSD_Na.write("%lf  %lf  %lf  %lf\n" %(sum_dx2_Na/Total_Na,sum_dy2_Na/Total_Na,sum_dz2_Na/Total_Na,sum_dr2_Na/Total_Na))
    output_MSD_Na.close()
  if Total_K != 0:
    output_MSD_K =open("z08_MSD_K.data",'w')
    #output_MSD_K.write("#dx2  dy2  dz2  MSD\n")
    output_MSD_K.write("%lf  %lf  %lf  %lf\n" %(sum_dx2_K/Total_K,sum_dy2_K/Total_K,sum_dz2_K/Total_K,sum_dr2_K/Total_K))
    output_MSD_K.close()
  if Total_Sr != 0:
    output_MSD_Sr =open("z09_MSD_Sr.data",'w')
    #output_MSD_Sr.write("#dx2  dy2  dz2  MSD\n")
    output_MSD_Sr.write("%lf  %lf  %lf  %lf\n" %(sum_dx2_Sr/Total_Sr,sum_dy2_Sr/Total_Sr,sum_dz2_Sr/Total_Sr,sum_dr2_Sr/Total_Sr))
    output_MSD_Sr.close()




