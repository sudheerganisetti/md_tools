
"""
@Author : Sudheer Ganisetti
@Date   : Mo 25. Dez 18:34:13 CET 2017
@         Do 20. Jun 20:34:10 CEST 2019
This code is to compute bond statistics
"""

import sys
import ganisetti_tools
import matplotlib as mpl
import matplotlib.pyplot as plt

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if __name__=="__main__":

  cmd1=ganisetti_tools.read_updated_command_line(sys.argv)
  # read parameter file
  read_parameter_file=ganisetti_tools.read_parameter_file(sys.argv[1])
  ganisetti_tools.print_error_message(read_parameter_file.error_status,read_parameter_file.error_messages)
  cmd=ganisetti_tools.read_command_line(read_parameter_file.all_arguments)

  if cmd.error != "none":
    sys.exit(0)

  # get the version details
  ganisetti_tools.get_ganisetti_tools_version()

  atom_type_sym2num    = cmd.atom_type_sym2num
  atom_type_num2sym    = cmd.atom_type_num2sym
  bond_length_sym2num  = cmd.bond_length_sys2num
  bond_length_num2num  = cmd.bond_length_num2num
  rc                   = cmd.rc
  given_anions_sym2num = cmd.given_anions_sym2num
  given_cations_sym2num= cmd.given_cations_sym2num
  given_formers_sym2num= cmd.given_formers_sym2num
  given_modifiers_sym2num = cmd.given_modifiers_sym2num

  # controllers for some specific glasses
  temp1=atom_type_sym2num.keys()
  BG_glass="no"
  if "Si" in temp1 and "P" in temp1 and "Na" in temp1 and "Ca" in temp1:
    BG_glass="yes"
  NAPS_glass="no"
  if "P" in temp1 and "Al" in temp1 and "Na" in temp1:
    NAPS_glass="yes"
  CMAS_glass="no"
  if "Ca" in temp1 and "Mg" in temp1 and "Al" in temp1 and "Si" in temp1:
    CMAS_glass="yes"

  # main loop starts here
  BASE_FILE=str(read_parameter_file.cur_file)
  # **************************************************************************************
  # get the atoms information
  ref_config    = ganisetti_tools.get_atoms_info_from_lammps(read_parameter_file.input_lammps_dump_file)
  cur_config    = ganisetti_tools.get_atoms_info_from_lammps(read_parameter_file.cur_lammps_dump_file)
  MAX_ATOM_TYPE=max(ref_config.type.values())

  # **************************************************************************************
  # compute total atoms of each atom type
  total_atoms_of_type_sym={}
  tmp=ref_config.type.values()
  tmp=list(tmp)
  for i in atom_type_sym2num.keys():
    total_atoms_of_type_sym[i]=tmp.count(atom_type_sym2num[i])

  # **************************************************************************************
  # computing nnl
  ref_config_nnl   = ganisetti_tools.compute_nnl(ref_config,rc,atom_type_num2sym)
  cur_config_nnl   = ganisetti_tools.compute_nnl(cur_config,rc,atom_type_num2sym)

  # **************************************************************************************
  # compute anions distribution
  ref_config_anions_distribution=ganisetti_tools.compute_anions_distribution(cmd, ref_config, ref_config_nnl)

  # **************************************************************************************
  # compute bond statistics
  bond_statistics = ganisetti_tools.compute_bond_statistics(given_cations_sym2num,ref_config,ref_config_nnl,cur_config_nnl)

  # **************************************************************************************
  # write bond statistics
  output1=open(BASE_FILE+str('_bond_statistics.data'),'w')
  output1.write("total bonds                  = %d\n" %(bond_statistics.total_bonds))
  output1.write("total broken bonds and %%    = %d \t %2.2lf\n" %(bond_statistics.total_broken_bonds,bond_statistics.total_broken_bonds/bond_statistics.total_bonds*100.0))
  output1.write("total new bonds and %%       = %d \t %2.2lf\n" %(bond_statistics.total_new_bonds,bond_statistics.total_new_bonds/bond_statistics.total_bonds*100.0))
  output1.write("total survived bonds and %%  = %d \t %2.2lf\n" %(bond_statistics.total_survived_bonds,bond_statistics.total_survived_bonds/bond_statistics.total_bonds*100.0))
  output1.write("\n\n")

# For NAPS glass
  if NAPS_glass == "yes":
    Na_type=atom_type_sym2num["Na"]
    Na_id__total_BA_neighbours  ={}
    Na_id__total_NBA_neighbours ={}
    Na_id__BA_neighbours_list   ={}
    Na_id__NBA_neighbours_list  ={}
    BA_and_NBA_neighbours_of_Na__Na_list ={}
    # initialization 
    for BA_neighbours in range(0,11):
      for NBA_neighbours in range(0,11):
        temp1={(BA_neighbours,NBA_neighbours):[]}
        BA_and_NBA_neighbours_of_Na__Na_list.update(temp1)

    # loop over each Na atom
    for i in ref_config.id:
      if ref_config.type[i] == Na_type:
        BA_count_for_i=0
        NBA_count_for_i=0
        temp_Na_id__BA_neighbours_list=[]
        temp_Na_id__NBA_neighbours_list=[]
        for j in ref_config_nnl.nnl[i]: # j=anion
          if j in ref_config_anions_distribution.BA_4CoordFormer_id2list.keys():
            BA_count_for_i=BA_count_for_i+1
            temp_Na_id__BA_neighbours_list.append(j)
          else :
            NBA_count_for_i=NBA_count_for_i+1
            temp_Na_id__NBA_neighbours_list.append(j)
        temp1={i:BA_count_for_i}
        temp2={i:NBA_count_for_i}
        Na_id__total_BA_neighbours.update(temp1)
        Na_id__total_NBA_neighbours.update(temp2)
        temp1={i:temp_Na_id__BA_neighbours_list}
        temp2={i:temp_Na_id__NBA_neighbours_list}
        Na_id__BA_neighbours_list.update(temp1)
        Na_id__NBA_neighbours_list.update(temp2)
        temp1=BA_and_NBA_neighbours_of_Na__Na_list[(BA_count_for_i,NBA_count_for_i)]
        temp1.append(i)
        temp2={(BA_count_for_i,NBA_count_for_i):temp1}
        BA_and_NBA_neighbours_of_Na__Na_list.update(temp2)

    output1.write("# BO_neighbours\tNBO_neighbours\tNumber_Of_Na_atoms\tNa-O_bonds\tBroken_Bonds\tSwitched_Bonds\tNew_Bonds\tSurvived_Bonds\n")
    for BA_neighbours in range(11):
      for NBA_neighbours in range(11):
        temp1=0
        temp2=0
        temp3=0
        temp4=0
        temp5=0
        for i in BA_and_NBA_neighbours_of_Na__Na_list[(BA_neighbours,NBA_neighbours)]:
          temp1=temp1 + len(bond_statistics.broken_bonds__id2list[i])
          temp2=temp2 + bond_statistics.switched_bonds__id2count[i]
          temp3=temp3 + len(bond_statistics.new_bonds__id2list[i])
          temp4=temp4 + len(bond_statistics.survived_bonds__id2list[i])
          temp5=temp5 + ref_config_nnl.nnl_count[i]
        output1.write("\t%d\t\t%d\t\t%d\t\t%d\t\t    %d\t\t   %d        \t%d   \t\t%d\n" %(BA_neighbours,NBA_neighbours,len(BA_and_NBA_neighbours_of_Na__Na_list[(BA_neighbours,NBA_neighbours)]),temp5,temp1,temp2,temp3,temp4))
    
    #for i in ref_config.id:
      #if ref_config.type[i] == Na_type:
        #print(len(bond_statistics.broken_bonds__id2list[i]),len(bond_statistics.new_bonds__id2list[i]),len(bond_statistics.survived_bonds__id2list[i]))
        #print(i, bond_statistics.broken_bonds__id2list[i])
        #print(i, bond_statistics.new_bonds__id2list[i])
        #print(i, ref_config_nnl.nnl[i])

  output1.close()




