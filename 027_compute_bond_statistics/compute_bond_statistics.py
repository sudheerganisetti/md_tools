
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
  output1=open(BASE_FILE+str('_bond_statistics_v01.data'),'w')
  output2=open(BASE_FILE+str('_bond_statistics_v02.data'),'w')
  output3=open(BASE_FILE+str('_bond_statistics_v03.data'),'w')
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

    # print the bond statistics of Na atoms based on the number of NBO and BO neighbors
    output1.write("# BO_neigh  NBO_neigh  tot_Na  tot_NaO  Brok_B  Switchd_B  New_B  Survd_B\n")
    output2.write("# BO_neigh  NBO_neigh  tot_Na_atoms  Na-O_bonds  tot_Brok_B  BO_Brok_B  ")
    output2.write("NBO_Brok_B  tot_New_B  BO_New_B  NBO_New_B  tot_Survd_B  BO_Survd_B  NBO_Survived_B\n")
    output3.write("#tot_BB  tot_BO_BB  tot_NBO_BB  tot_NB  tot_BO_NB  tot_NBO_NB  tot_SB  tot_BO_SB  tot_NBO_SB\n")
    tot_BO_BB   = 0
    tot_NBO_BB  = 0
    tot_BO_NB   = 0
    tot_NBO_NB  = 0
    tot_BO_SB   = 0
    tot_NBO_SB  = 0
    tot_NaO_BB  = 0
    tot_NaO_NB  = 0
    tot_NaO_SB  = 0

    for BA_neighbours in range(11):
      for NBA_neighbours in range(11):
        tot_BB=0
        tot_switchB=0
        tot_NB=0
        tot_SB=0
        tot_NaO_bonds=0
        BA_count_of_new_bonds       = 0
        NBA_count_of_new_bonds      = 0
        BA_count_of_broken_bonds    = 0
        NBA_count_of_broken_bonds   = 0
        BA_count_of_survived_bonds  = 0
        NBA_count_of_survived_bonds = 0
        for i in BA_and_NBA_neighbours_of_Na__Na_list[(BA_neighbours,NBA_neighbours)]:
          tot_BB          = tot_BB        + len(bond_statistics.broken_bonds__id2list[i])
          tot_switchB     = tot_switchB   + bond_statistics.switched_bonds__id2count[i]
          tot_NB          = tot_NB        + len(bond_statistics.new_bonds__id2list[i])
          tot_SB          = tot_SB        + len(bond_statistics.survived_bonds__id2list[i])
          tot_NaO_bonds   = tot_NaO_bonds + ref_config_nnl.nnl_count[i]
          
          for j in bond_statistics.new_bonds__id2list[i]:
            if j in ref_config_anions_distribution.BA_4CoordFormer_id2list.keys():
              BA_count_of_new_bonds       = BA_count_of_new_bonds + 1
              tot_BO_NB                   = tot_BO_NB + 1
            else: 
              NBA_count_of_new_bonds      = NBA_count_of_new_bonds +1
              tot_NBO_NB                  = tot_NBO_NB + 1
          for j in bond_statistics.broken_bonds__id2list[i]:
            if j in ref_config_anions_distribution.BA_4CoordFormer_id2list.keys():
              BA_count_of_broken_bonds    = BA_count_of_broken_bonds + 1
              tot_BO_BB                   = tot_BO_BB + 1
            else: 
              NBA_count_of_broken_bonds   = NBA_count_of_broken_bonds +1
              tot_NBO_BB                  = tot_NBO_BB + 1
          for j in bond_statistics.survived_bonds__id2list[i]:
            if j in ref_config_anions_distribution.BA_4CoordFormer_id2list.keys():
              BA_count_of_survived_bonds  = BA_count_of_survived_bonds + 1
              tot_BO_SB                   = tot_BO_SB + 1
            else: 
              NBA_count_of_survived_bonds = NBA_count_of_survived_bonds +1
              tot_NBO_SB                  = tot_NBO_SB + 1
        tot_NaO_BB = tot_NaO_BB + tot_BB 
        tot_NaO_NB = tot_NaO_NB + tot_NB 
        tot_NaO_SB = tot_NaO_SB + tot_SB 

        output1.write("\t%d \t %d \t" %(BA_neighbours,NBA_neighbours))
        output1.write("%d \t" %(len(BA_and_NBA_neighbours_of_Na__Na_list[(BA_neighbours,NBA_neighbours)])))
        output1.write("%d \t %d \t   %d \t   %d \t   %d \n" %(tot_NaO_bonds,tot_BB,tot_switchB,tot_NB,tot_SB))
        output2.write("\t%d\t%d\t" %(BA_neighbours,NBA_neighbours))
        output2.write("%d \t\t" %(len(BA_and_NBA_neighbours_of_Na__Na_list[(BA_neighbours,NBA_neighbours)])))
        output2.write("%d \t   %d \t\t %d \t %d \t\t" %(tot_NaO_bonds,tot_BB,BA_count_of_broken_bonds,NBA_count_of_broken_bonds))
        output2.write("%d \t %d \t   %d \t\t" %(tot_NB,BA_count_of_new_bonds,NBA_count_of_new_bonds))
        output2.write("%d \t \t%d \t %d \n" %(tot_SB,BA_count_of_survived_bonds,NBA_count_of_survived_bonds))
    output3.write("  %d \t    %d \t %d \t"  %(tot_NaO_BB,tot_BO_BB,tot_NBO_BB))
    output3.write(" %d \t %d \t\t %d \t" %(tot_NaO_NB,tot_BO_NB,tot_NBO_NB))
    output3.write(" %d \t   %d \t\t %d \n" %(tot_NaO_SB,tot_BO_SB,tot_NBO_SB))

    output5=open(BASE_FILE+str('_detailed_BB.data'),'w')
    output6=open(BASE_FILE+str('_detailed_NB.data'),'w')
    output7=open(BASE_FILE+str('_detailed_SB.data'),'w')
    BA_and_NBA_speciation_BB={}
    BA_and_NBA_speciation_NB={}
    BA_and_NBA_speciation_SB={}
    BA_and_NBA_speciation_BB1={}
    BA_and_NBA_speciation_NB1={}
    BA_and_NBA_speciation_SB1={}
    for i in ref_config.id:
      if ref_config.type[i] == Na_type:
        # collect detailed broken bonds
        for j in bond_statistics.broken_bonds__id2list[i]:
          temp1=( ("Si",ref_config_nnl.nnl_type_sym[j].count("Si")),("Al",ref_config_nnl.nnl_type_sym[j].count("Al")),("P",ref_config_nnl.nnl_type_sym[j].count("P")) )
          try:
            temp2={temp1:BA_and_NBA_speciation_BB[temp1]+1}
          except:
            temp2={temp1:1}
          BA_and_NBA_speciation_BB.update(temp2)

        # collect detailed new bonds
        for j in bond_statistics.new_bonds__id2list[i]:
          temp1=( ("Si",ref_config_nnl.nnl_type_sym[j].count("Si")),("Al",ref_config_nnl.nnl_type_sym[j].count("Al")),("P",ref_config_nnl.nnl_type_sym[j].count("P")) )
          try:
            temp2={temp1:BA_and_NBA_speciation_NB[temp1]+1}
          except:
            temp2={temp1:1}
          BA_and_NBA_speciation_NB.update(temp2)

        # collect detailed survived bonds
        for j in bond_statistics.survived_bonds__id2list[i]:
          temp1=( ("Si",ref_config_nnl.nnl_type_sym[j].count("Si")),("Al",ref_config_nnl.nnl_type_sym[j].count("Al")),("P",ref_config_nnl.nnl_type_sym[j].count("P")) )
          try:
            temp2={temp1:BA_and_NBA_speciation_SB[temp1]+1}
          except:
            temp2={temp1:1}
          BA_and_NBA_speciation_SB.update(temp2)

    output5.write("# num_Si\tnum_Al\tnum_P\ttot_BB\n")
    BA_and_NBA_speciation_BB1={k: v for k, v in sorted(BA_and_NBA_speciation_BB.items(), key=lambda item: item[1])}
    for i in BA_and_NBA_speciation_BB1.keys():
      for j in i:
        output5.write("%d\t" %(j[1]))
      output5.write("\t%d\n" %(BA_and_NBA_speciation_BB1[i]))
        
    output6.write("# num_Si\tnum_Al\tnum_P\ttot_NB\n")
    BA_and_NBA_speciation_NB1={k: v for k, v in sorted(BA_and_NBA_speciation_NB.items(), key=lambda item: item[1])}
    for i in BA_and_NBA_speciation_NB1.keys():
      for j in i:
        output6.write("%d\t" %(j[1]))
      output6.write("\t%d\n" %(BA_and_NBA_speciation_NB1[i]))

    output7.write("# num_Si\tnum_Al\tnum_P\ttot_SB\n")
    BA_and_NBA_speciation_SB1={k: v for k, v in sorted(BA_and_NBA_speciation_SB.items(), key=lambda item: item[1])}
    for i in BA_and_NBA_speciation_SB1.keys():
      for j in i:
        output7.write("%d\t" %(j[1]))
      output7.write("\t%d\n" %(BA_and_NBA_speciation_SB1[i]))
      
    output5.close()
    output6.close()
    output7.close()

    for i in range(11):
      for j in range(11):
        if len(BA_and_NBA_neighbours_of_Na__Na_list[(i,j)]) != 0:
          output4=open("List_of_Na_atoms_with___BO_"+str(i)+"___NBO_"+str(j)+".data",'w')
          output4.write("# BO %d \t\t NBO %d\n" %(i,j))
          for k in BA_and_NBA_neighbours_of_Na__Na_list[(i,j)]:
            output4.write("%d \n" %(k))
          output4.close()

  output1.close()
  output2.close()
  output3.close()

