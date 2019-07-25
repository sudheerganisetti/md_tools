
"""
@Author : Sudheer Ganisetti
@Date   : Mo 25. Dez 18:34:13 CET 2017
@         Do 20. Jun 20:34:10 CEST 2019
This code is to compute nnl, nnlchk,  connectivity(Qn) and A-O-A triplets for the sample (SiO2 + Al2O3 + P2O5 + Na2O + CaO + MgO + CaF2)
"""

import numpy as np
import sys
import math as mt
import itertools as it
import datetime
import ganisetti_tools

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if __name__=="__main__":

  # read command line first
  cmd=ganisetti_tools.read_command_line(sys.argv)
  if cmd.error != "none":
    sys.exit(0)
  atom_type_sym2num    = cmd.atom_type_sym2num
  atom_type_num2sym    = cmd.atom_type_num2sym
  bond_length_sym2num  = cmd.bond_length_sys2num
  bond_length_num2num  = cmd.bond_length_num2num
  rc                   = cmd.rc
  given_anions_sym2num = cmd.given_anions_sym2num
  given_cations_sym2num= cmd.given_cations_sym2num
  given_formers_sym2num= cmd.given_formers_sym2num
  given_modifiers_sym2num = cmd.given_modifiers_sym2num

  # main loop starts here
  BASE_FILE=str(sys.argv[1])
  LAMMPS_DUMP_FILE=BASE_FILE+str('.dump')
  MAX_NEIGHBOURS=8
  MAX_ATOM_TYPES=8

  # **************************************************************************************
  # get the atoms information
  config1	= ganisetti_tools.get_atoms_info_from_lammps(LAMMPS_DUMP_FILE)
  MAX_ATOM_TYPE=max(config1.type.values())

  # **************************************************************************************
  # compute total atoms of each atom type
  total_atoms_of_type_sym={}
  tmp=config1.type.values()
  tmp=list(tmp)
  for i in atom_type_sym2num.keys():
    total_atoms_of_type_sym[i]=tmp.count(atom_type_sym2num[i])

  # **************************************************************************************
  # computing nnl
  config1_nnl   = ganisetti_tools.compute_nnl(config1,rc,atom_type_num2sym)

  output_nnl=open(BASE_FILE+str('.nnl'),'w')
  output_nnl.write("# nnl file; used cutoffs are: ")
  for i in atom_type_num2sym.keys():
    if atom_type_sym2num['O'] != i:
      output_nnl.write("%s-O = %.2lf ; " %(atom_type_num2sym[i],rc[i][atom_type_sym2num['O']]))
  output_nnl.write("\n")

  for i in config1.id:
    output_nnl.write("%d\t" %(i))
    for j in config1_nnl.nnl[i]:
      output_nnl.write("%d\t" %(j))
    output_nnl.write("\n")

  # **************************************************************************************
  # writing out chkpt with nnl status
  output_nnl_imd=open(BASE_FILE+str("_withnnlstatus.atoms"),'w')
  ganisetti_tools.write_imd_header(output_nnl_imd,config1.box,rc,atom_type_sym2num)
  for i in config1.id:
    #output_nnl_imd.write("%d  %d  %lf  %lf  %lf %d \n" %(i,config1.type[i],config1.posx[i],config1.posy[i],config1.posz[i],config1_nnl.nnl_count[i]))
    ganisetti_tools.write_imd_atom(output_nnl_imd,i,config1,config1_nnl)
  output_nnl_imd.close()

  # **************************************************************************************
  # writing out individual cations and their anions
  for i in atom_type_sym2num.keys():
    if i != "O":
      output=open(BASE_FILE+str("_")+str(i)+str("-O.atoms"),'w')
      ganisetti_tools.write_imd_header(output,config1.box,rc,atom_type_sym2num)
      SelectAtoms1={}
      for j in config1.id:
        SelectAtoms1[j]=0
      for j in config1.id:
        if i == atom_type_num2sym[config1.type[j]]:
          tmp={j:1}
          SelectAtoms1.update(tmp)
          for k in config1_nnl.nnl[j]:
            tmp={k:1}
            SelectAtoms1.update(tmp)
      for j in config1.id:
        if SelectAtoms1[j] == 1:
          ganisetti_tools.write_imd_atom(output,j,config1,config1_nnl)
      output.close()


  # **************************************************************************************
  # compute coordination of each atom and average atom type
  config1_coord=ganisetti_tools.compute_coordination(config1,config1_nnl,atom_type_num2sym)
  for i in atom_type_sym2num.keys():
    output1=open(BASE_FILE+str("_")+str(i)+str("_coord"),'w')
    output1.write("# %s : coord  number_of_atoms  %% of atoms\n" %(str(i)))
    for j in range(config1_nnl.max_nnl_each_atom_type_sym[i]+1):
      temp1=config1_coord.individual_coord_sym[(i,j)]
      output1.write("%d\t%d\t\t%.2lf\n" %(j,temp1,float(temp1*100.0/total_atoms_of_type_sym[i])))
    output1.close()
  output1=open(BASE_FILE+str("_average_coord"),'w')
  output1.write("# AtomType_sym  AtomType_num   Number_of_Atoms AverageCoordination\n")
  for i in atom_type_sym2num.keys():
    output1.write("%s\t%d\t%d\t%.2lf\n" %(i,atom_type_sym2num[i],total_atoms_of_type_sym[i],config1_coord.avg_coord_sym[i]))
  output1.close()


  # **************************************************************************************
  # compute environment of each atom
  config1_env=ganisetti_tools.compute_each_atom_environment(config1,config1_nnl,atom_type_sym2num,atom_type_num2sym)
  '''
  for i in config1.id:
      for j in given_anions_sym2num.values():
        if config1.type[i] == j:
          for k in config1_env.env_atomid[i]:
  '''


  # **************************************************************************************
  # compute Q of each network former
  '''
  # This is deprecated
  for i in given_formers_sym2num.keys():
    output1=open(BASE_FILE+str("_Q")+str(i)+str(".data"),'w')
    output1.write("# Q(n) speciation of %s\n" %(i))
    output1.write("# n  number_of_units\n")
    config1_Q=ganisetti_tools.compute_Q(cmd,config1,config1_nnl,config1_env,i)
    for j in range(5):
      output1.write("%d %d\t%.2lf \n" %(j,config1_Q.Q[j],config1_Q.Q[j]*100.0/total_atoms_of_type_sym[i]))
    output1.close()
  '''

  config1_Q=ganisetti_tools.compute_Q(cmd,config1,config1_nnl,config1_env)
  for i in given_formers_sym2num.keys():
    output1=open(BASE_FILE+str("_Q")+str(i)+str(".data"),'w')
    output1.write("# Q(n) speciation of %s\n" %(i))
    output1.write("# n  number_of_units\n")
    for j in range(5):
      output1.write("%d %d\t%.2lf \n" % (j, config1_Q.Q_summary[(i,j)], config1_Q.Q_summary[(i,j)] * 100.0 / total_atoms_of_type_sym[i]))
      #print(config1_Q.Q_non4CoordFormers_list[i,j])
    output1.close()
  #for i in config1.id:
  #  if config1.type[i] in given_formers_sym2num.values():
  #    print(config1.type[i],config1_Q.Q_status_atomid2num[i])

  # **************************************************************************************
  # compute triplets
  config1_triplets=ganisetti_tools.compute_triplets(cmd,config1,config1_nnl)
  #for i in atom_type_sym2num.keys():
  #  print(i,config1_triplets.total_triplets_sym2count[i])

  output1=open(BASE_FILE+str("_triplets.data"),'w')
  output1.write("# Triplets based on anions\n")
  for i in atom_type_sym2num.keys():                                         # This is B in A-B-C
    if i in given_anions_sym2num.keys():
      temp1=config1_triplets.total_triplets_sym2count[i]
      output1.write("\n# Total triplets with centered with %s = %d\n" %(i,temp1))
      for j in atom_type_sym2num.keys():                                       # This is A in A-B-C
        for k in atom_type_sym2num.keys():                                     # This is C in A-B-C
          for m in range(config1_nnl.max_nnl_each_atom_type_sym[j] + 1):
            for n in range(config1_nnl.max_nnl_each_atom_type_sym[k] + 1):
              temp2=config1_triplets.triplets_AmBCn2count[(j,m,i,k,n)]
              if temp2 != 0:
                output1.write("%s %d - %s - %s %d\t = %d\t  %.2lf\n" %(j,m,i,k,n,temp2,temp2*100.0/temp1))
                #print("%s %d - %s - %s %d\t = %d\t  %.2lf\n" %(j,m,i,k,n,temp2,temp2*100.0/temp1))

  # **************************************************************************************
  # compute  ions distribution
  # this is deprecated
  '''
  sum1=0
  terminal_ions_sym_list=list(given_formers_sym2num.keys())
  config1_O_distribution_among_formers = ganisetti_tools.compute_ion_distribution(cmd, config1, config1_nnl,"O",list(given_formers_sym2num.keys()))
  for A in terminal_ions_sym_list:
    for C in terminal_ions_sym_list:
      for m in range(config1_nnl.max_nnl_each_atom_type_sym[A] + 1):
        for n in range(config1_nnl.max_nnl_each_atom_type_sym[C] + 1):
          temp1=config1_O_distribution_among_formers.triplets_AmBCn2count[(A,m,"O",C,n)]
          if temp1 != 0:
            temp2=temp1*100.0/config1_O_distribution_among_formers.total_triplets
            print("%s %d - O - %s %d\t=  %4d\t\t%.2lf " %(A,m,C,n,temp1,temp2))
            sum1=sum1+temp1
  print("total_o_atoms = %d ; sum of above triplets = %d" %(total_atoms_of_type_sym["O"],sum1))
  '''

  # **************************************************************************************
  # compute  anions distribution
  config1_anions_distribution=ganisetti_tools.compute_anions_distribution(cmd, config1, config1_nnl)

  print("total BA              = %d" %(config1_anions_distribution.total_anions_of_BA_4CoordFormer_sym2count['O']))
  print("total NBA             = %d" %(config1_anions_distribution.total_anions_of_NBA_former_sym2count['O']+config1_anions_distribution.total_anions_of_BA_non4CoordFormer_sym2count['O']))
  print("total triclusters     = %d" %(config1_anions_distribution.total_anions_of_tri_clusters_sym2count['O']))
  print("anions with zero NF   = %d" %(config1_anions_distribution.total_anions_with_zero_formers_sym2count['O']))
  print("any other type anions = %d" %(config1_anions_distribution.total_anions_of_any_other_type_sym2count['O']))
  set_a={}
  set_c={}
  count=0
  for i in given_formers_sym2num.keys():
    temp1={count:i}
    set_a.update(temp1)
    set_c.update(temp1)
    count=count+1
  for B in given_anions_sym2num.keys():
    for i in range(count):
      A=set_a[i]
      for j in range(i,count):
        C=set_c[j]
        for m in range(config1_nnl.max_nnl_each_atom_type_sym[A] + 1):
          for n in range(config1_nnl.max_nnl_each_atom_type_sym[C] + 1):
            temp1=config1_anions_distribution.triplets_AmBCn2count[(A,m,"O",C,n)]
            if temp1 != 0:
              print("%s %d %s %s %d   = %d " %(A,m,B,C,n,temp1))
        print("\n")
      print("\n")
