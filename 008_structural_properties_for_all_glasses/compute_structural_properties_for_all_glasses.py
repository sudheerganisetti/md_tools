#!/home/sudheer/bin/python
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
  error1=0
  error2=0
  debug = "yes"
  tot_argv=len(sys.argv)
  CREDBG    = '\33[41m'
  CREDBGEND = '\x1b[0m'
  while error1 == 0:
    # checking for minimum arguments
    if tot_argv < 8:
      error1=1
      error1_statement="Error info: a minimum of 7 arguments should be passed to the program"
      break
    if (tot_argv%2) !=0:
      error1=1
      error1_statement="Error info: odd number of arguments are expecting to be passed to the program"
      break
    commandline=[]
    recognized_argvs_idex=[]
    for i in sys.argv:
      commandline.append(i)
    commandline=list(commandline)
    atom_type={}		# store all given atom types				# atom_type['Si'] = 1
    bond_length={} 		# store all given bond lengths				# bond_length['Si'] = 2.2
    atom_type_sym2num={}	# store all given atom types based on chemical symbol 	# atom_type_sym2num['Si'] = 1
    atom_type_num2sym={}        # store all given atom types based on atom type		# atom_type_sym2num['1'] = Si
    bond_length_sym2num={}	# store all given bond lengths based on chemical symbol	# bond_length_sym2num['Si'] = 2.2 
    bond_length_num2num={}	# store all given bond lengths based on atom type	# bond_length_num2num['1'] = 2.2
    
    given_all_atom_types=[]
    # search for atom types
    if "-O" in commandline:
      O_index=commandline.index("-O")
      given_O_type=int(commandline[O_index+1])
      atom_type['O']			=given_O_type
      atom_type_sym2num['O']		=given_O_type
      atom_type_num2sym[given_O_type]	='O'
      given_all_atom_types.append(given_O_type)
      recognized_argvs_idex.append(O_index)
      recognized_argvs_idex.append(O_index+1)
    if "-Si" in commandline:
      Si_index=commandline.index("-Si")
      given_Si_type=int(commandline[Si_index+1])
      atom_type['Si']			= given_Si_type
      atom_type_sym2num['Si']		= given_Si_type
      atom_type_num2sym[given_Si_type]	= 'Si'
      given_all_atom_types.append(given_Si_type)
      recognized_argvs_idex.append(Si_index)
      recognized_argvs_idex.append(Si_index+1)
    if "-Al" in commandline:
      Al_index=commandline.index("-Al")
      given_Al_type=int(commandline[Al_index+1])
      atom_type['Al']			=given_Al_type
      atom_type_sym2num['Al']           =given_Al_type
      atom_type_num2sym[given_Al_type]  ='Al'
      given_all_atom_types.append(given_Al_type)
      recognized_argvs_idex.append(Al_index)
      recognized_argvs_idex.append(Al_index+1)
    if "-P" in commandline:
      P_index=commandline.index("-P")
      given_P_type=int(commandline[P_index+1])
      atom_type['P']			=given_P_type
      atom_type_sym2num['P']            =given_P_type
      atom_type_num2sym[given_P_type]   ='P'
      given_all_atom_types.append(given_P_type)
      recognized_argvs_idex.append(P_index)
      recognized_argvs_idex.append(P_index+1)
    if "-Ca" in commandline:
      Ca_index=commandline.index("-Ca")
      given_Ca_type=int(commandline[Ca_index+1])
      atom_type['Ca']			=given_Ca_type
      atom_type_sym2num['Ca']           =given_Ca_type
      atom_type_num2sym[given_Ca_type]  ='Ca'
      given_all_atom_types.append(given_Ca_type)
      recognized_argvs_idex.append(Ca_index)
      recognized_argvs_idex.append(Ca_index+1)
    if "-Mg" in commandline:
      Mg_index=commandline.index("-Mg")
      given_Mg_type=int(commandline[Mg_index+1])
      atom_type['Mg']			=given_Mg_type
      atom_type_sym2num['Mg']           =given_Mg_type
      atom_type_num2sym[given_Mg_type]  ='Mg'
      given_all_atom_types.append(given_Mg_type)
      recognized_argvs_idex.append(Mg_index)
      recognized_argvs_idex.append(Mg_index+1)
    if "-Sr" in commandline:
      Sr_index=commandline.index("-Sr")
      given_Sr_type=int(commandline[Sr_index+1])
      atom_type['Sr']			=given_Sr_type
      atom_type_sym2num['Sr']           =given_Sr_type
      atom_type_num2sym[given_Sr_type]  ='Sr'
      given_all_atom_types.append(given_Sr_type)
      recognized_argvs_idex.append(Sr_index)
      recognized_argvs_idex.append(Sr_index+1)
    if "-Na" in commandline:
      Na_index=commandline.index("-Na")
      given_Na_type=int(commandline[Na_index+1])
      atom_type['Na']			=given_Na_type
      atom_type_sym2num['Na']           =given_Na_type
      atom_type_num2sym[given_Na_type]  ='Na'
      given_all_atom_types.append(given_Na_type)
      recognized_argvs_idex.append(Na_index)
      recognized_argvs_idex.append(Na_index+1)
    if "-F" in commandline:
      F_index=commandline.index("-F")
      given_F_type=int(commandline[F_index+1])
      atom_type['F']                   =given_F_type
      atom_type_sym2num['F']           =given_F_type
      atom_type_num2sym[given_Na_type]  ='F'
      given_all_atom_types.append(given_F_type)
      recognized_argvs_idex.append(F_index)
      recognized_argvs_idex.append(F_index+1)

    # O atom type must be given
    if "O" not in atom_type_sym2num.keys():
      error1=1
      error1_statement="Error info: Oxygen atom type is not given"
      break
    # initilize a cutoff radii array
    rc=[[0.0 for i in range(max(given_all_atom_types)+1)] for j in range(max(given_all_atom_types)+1)]
    # search for bond lengths
    if "-SiO" in commandline:
      if "Si" not in atom_type_sym2num:
        error1=1
        error1_statement="Error info: Si atom type is not given but SiO bond length is given"
        break
      tmp_index=commandline.index("-SiO")
      bond_length['SiO']					= float(commandline[tmp_index+1])
      bond_length_sym2num['Si']					= float(commandline[tmp_index+1])
      bond_length_num2num[atom_type_sym2num['Si']]		= float(commandline[tmp_index+1])
      rc[atom_type_sym2num['Si']][atom_type_sym2num['O']]	= float(commandline[tmp_index+1])
      rc[atom_type_sym2num['O']][atom_type_sym2num['Si']]       = float(commandline[tmp_index+1])
      recognized_argvs_idex.append(tmp_index)
      recognized_argvs_idex.append(tmp_index+1)
    if "-AlO" in commandline:
      if "Al" not in atom_type_sym2num:
        error1=1
        error1_statement="Error info: Al atom type is not given but AlO bond length is given"
        break
      tmp_index=commandline.index("-AlO")
      bond_length['AlO']                                	= float(commandline[tmp_index+1])
      bond_length_sym2num['Al']					= float(commandline[tmp_index+1])
      bond_length_num2num[atom_type_sym2num['Al']]      	= float(commandline[tmp_index+1])
      rc[atom_type_sym2num['Al']][atom_type_sym2num['O']]       = float(commandline[tmp_index+1])
      rc[atom_type_sym2num['O']][atom_type_sym2num['Al']]       = float(commandline[tmp_index+1])
      recognized_argvs_idex.append(tmp_index)
      recognized_argvs_idex.append(tmp_index+1)
    if "-PO" in commandline:
      if "P" not in atom_type_sym2num:
        error1=1
        error1_statement="Error info: P atom type is not given but PO bond length is given"
        break
      tmp_index=commandline.index("-PO")
      bond_length['PO']						= float(commandline[tmp_index+1])
      bond_length_sym2num['P']					= float(commandline[tmp_index+1])
      bond_length_num2num[atom_type_sym2num['P']]		= float(commandline[tmp_index+1])
      rc[atom_type_sym2num['P']][atom_type_sym2num['O']]       	= float(commandline[tmp_index+1])
      rc[atom_type_sym2num['O']][atom_type_sym2num['P']]       	= float(commandline[tmp_index+1])
      recognized_argvs_idex.append(tmp_index)
      recognized_argvs_idex.append(tmp_index+1)
    if "-CaO" in commandline:
      if "Ca" not in atom_type_sym2num:
        error1=1
        error1_statement="Error info: Ca atom type is not given but CaO bond length is given"
        break
      tmp_index=commandline.index("-CaO")
      bond_length['CaO']                                	= float(commandline[tmp_index+1])
      bond_length_sym2num['Ca']					= float(commandline[tmp_index+1])
      bond_length_num2num[atom_type_sym2num['Ca']]      	= float(commandline[tmp_index+1])
      rc[atom_type_sym2num['Ca']][atom_type_sym2num['O']]       = float(commandline[tmp_index+1])
      rc[atom_type_sym2num['O']][atom_type_sym2num['Ca']]       = float(commandline[tmp_index+1])
      recognized_argvs_idex.append(tmp_index)
      recognized_argvs_idex.append(tmp_index+1)
    if "-MgO" in commandline:
      if "Mg" not in atom_type_sym2num:
        error1=1
        error1_statement="Error info: Mg atom type is not given but MgO bond length is given"
        break
      tmp_index=commandline.index("-MgO")
      bond_length['MgO']                                	= float(commandline[tmp_index+1])
      bond_length_sym2num['Mg']                        		= float(commandline[tmp_index+1])
      bond_length_num2num[atom_type_sym2num['Mg']]      	= float(commandline[tmp_index+1])
      rc[atom_type_sym2num['Mg']][atom_type_sym2num['O']]       = float(commandline[tmp_index+1])
      rc[atom_type_sym2num['O']][atom_type_sym2num['Mg']]       = float(commandline[tmp_index+1])
      recognized_argvs_idex.append(tmp_index)
      recognized_argvs_idex.append(tmp_index+1)
    if "-NaO" in commandline:
      if "Na" not in atom_type_sym2num:
        error1=1
        error1_statement="Error info: Na atom type is not given but NaO bond length is given"
        break
      tmp_index=commandline.index("-NaO")
      bond_length['NaO']                                	= float(commandline[tmp_index+1])
      bond_length_sym2num['Na']                        		= float(commandline[tmp_index+1])
      bond_length_num2num[atom_type_sym2num['Na']]      	= float(commandline[tmp_index+1])
      rc[atom_type_sym2num['Na']][atom_type_sym2num['O']]       = float(commandline[tmp_index+1])
      rc[atom_type_sym2num['O']][atom_type_sym2num['Na']]       = float(commandline[tmp_index+1])
      recognized_argvs_idex.append(tmp_index)
      recognized_argvs_idex.append(tmp_index+1)
    if "-SrO" in commandline:
      if "Sr" not in atom_type_sym2num:
        error1=1
        error1_statement="Error info: Sr atom type is not given but SrO bond length is given"
        break
      tmp_index=commandline.index("-SrO")
      bond_length['SrO']                                	= float(commandline[tmp_index+1])
      bond_length_sym2num['Sr']                        		= float(commandline[tmp_index+1])
      bond_length_num2num[atom_type_sym2num['Sr']]      	= float(commandline[tmp_index+1])
      rc[atom_type_sym2num['Sr']][atom_type_sym2num['O']]       = float(commandline[tmp_index+1])
      rc[atom_type_sym2num['O']][atom_type_sym2num['Sr']]       = float(commandline[tmp_index+1])
      recognized_argvs_idex.append(tmp_index)
      recognized_argvs_idex.append(tmp_index+1)


    if "-CaF" in commandline:
      if "Ca" not in atom_type_sym2num:
        error1=1
        error1_statement="Error info: Ca atom type is not given but CaF bond length is given"
        break
      if "F" not in atom_type_sym2num:
        error1=1
        error1_statement="Error info: F atom type is not given but CaF bond length is given"
        break
      tmp_index=commandline.index("-CaF")
      bond_length['CaF']                                        = float(commandline[tmp_index+1])
      bond_length_sym2num['Ca']                                 = float(commandline[tmp_index+1])
      bond_length_num2num[atom_type_sym2num['Ca']]              = float(commandline[tmp_index+1])
      rc[atom_type_sym2num['Ca']][atom_type_sym2num['F']]       = float(commandline[tmp_index+1])
      rc[atom_type_sym2num['F']][atom_type_sym2num['Ca']]       = float(commandline[tmp_index+1])
      recognized_argvs_idex.append(tmp_index)
      recognized_argvs_idex.append(tmp_index+1)
    
    # checking for bond lengths of all given atom types
    for i in atom_type_sym2num:
      if i != "O":
        if i not in bond_length_sym2num:
          error1=1
          error1_statement="Error info: " + str(i) +" atom type is given but " + str(i) + "O bond length is not given"
          break
    # checking for any unrecognized arguments passed into command line
    for i in range(2,len(commandline)):
      if i not in recognized_argvs_idex:
        error1=1
        error1_statement="Error info: " + str(commandline[i]) + " is an unrecognozed argument passed into the command line"
        break

    break # terminate the main while loop

  if debug == "no":
      for i in atom_type_sym2num:
        print(i,atom_type[i],atom_type_sym2num[i])
        if i != "O":
          print(bond_length_sym2num[i],bond_length_num2num[atom_type_sym2num[i]])
      print("\n\n")
      for i in bond_length_sym2num:
        print(i,bond_length[i+str("O")],bond_length_sym2num[i])
    
  if error1 != 0:
     ganisetti_tools.banner()
     print("************************************** S. Ganisetti **************************************")
     print("Error: usage is wrong")
     print("./this_program  ChkptFile -O 1 -Si 2 -Al 3 -Ca 4 -Mg 5 -SiO 2.0 -AlO 2.0 -CaO 3.0 -MgO 3.0")
     print("The program is prepared for chemical compositions: SiO2, Al2O3, P2O5, Na2O, CaO and MgO")
     print("Please specify both atom types and the cutoff radii of all required pair of atoms")
     print("%s %s %s" %(CREDBG,str(error1_statement),CREDBGEND))
     print("******************************************************************************************")
     sys.exit(0)

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
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
  for i in given_all_atom_types:
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
  ganisetti_tools.write_imd_header(output_nnl_imd,config1.box,rc,atom_type_sym2num,atom_type_num2sym)
  for i in config1.id:
    #output_nnl_imd.write("%d  %d  %lf  %lf  %lf %d \n" %(i,config1.type[i],config1.posx[i],config1.posy[i],config1.posz[i],config1_nnl.nnl_count[i]))
    ganisetti_tools.write_imd_atom(output_nnl_imd,i,config1,config1_nnl)
  output_nnl_imd.close()

  # **************************************************************************************
  # writing out individual cations and their anions
  for i in atom_type_sym2num.keys():
    if i != "O":
      output=open(BASE_FILE+str("_")+str(i)+str("-O.atoms"),'w')
      ganisetti_tools.write_imd_header(output,config1.box,rc,atom_type_sym2num,atom_type_num2sym)
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
  # compute environment of each atom
  config1_env=ganisetti_tools.compute_each_atom_environment(config1,config1_nnl,atom_type_sym2num)

  # **************************************************************************************
  # compute coordination of each atom
  config1_coord=ganisetti_tools.compute_coordination(config1,config1_nnl,atom_type_num2sym)
  output1=open(BASE_FILE+str("_average_coord"),'w')
  output1.write("# AtomType_sym  AtomType_num   Number_of_Atoms AverageCoordination\n")
  for i in atom_type_sym2num.keys():
    output1.write("%s\t%d\t%d\t%.2lf\n" %(i,atom_type_sym2num[i],total_atoms_of_type_sym[i],config1_coord.coord_sym[i]))
  output1.close()

