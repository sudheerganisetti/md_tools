
"""
@Author : Sudheer Ganisetti
@Date   : Mo 25. Dez 18:34:13 CET 2017
@         Do 20. Jun 20:34:10 CEST 2019
This code is to compute nnl, nnlchk,  connectivity(Qn) and A-O-A triplets for the sample (SiO2 + Al2O3 + P2O5 + Na2O + CaO + MgO + CaF2)
"""

import sys
import ganisetti_tools
import matplotlib as mpl
import matplotlib.pyplot as plt

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if __name__=="__main__":

  # read command line first
  cmd=ganisetti_tools.read_command_line(sys.argv)
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
  # config1_nnl.nnl = {0:[1,2,3,4]}
  # config1_nnl.nnl_count = {0:4}
  # config1_nnl.nnl_type_num = {0:[1,1,1,1]}
  # config1_nnl.nnl_type_sym = {0:[Si,Si,Si,Si]}
  # config1_nnl.nnl_each_pair_distance={0:[1.58,1.62,1.59,1.60]}
  # config1_nnl.max_nnl_each_atom_type_sym={O:4,Si:4,Al:6,Ca:8}

  # **************************************************************************************
  # writing out chkpt with nnl status
  output1=open(BASE_FILE+str("_withnnlstatus.atoms"),'w')
  ganisetti_tools.write_imd_header(output1,config1.box,rc,atom_type_sym2num)
  for i in config1.id:
    #output_nnl_imd.write("%d  %d  %lf  %lf  %lf %d \n" %(i,config1.type[i],config1.posx[i],config1.posy[i],config1.posz[i],config1_nnl.nnl_count[i]))
    ganisetti_tools.write_imd_atom(output1,i,config1,config1_nnl)
  output1.close()

  # **************************************************************************************
  # writing out individual cations and their anions
  for i in atom_type_sym2num.keys():
    if i not in given_anions_sym2num.keys():
      output1=open(BASE_FILE+str("_")+str(i)+str("-anions.atoms"),'w')
      ganisetti_tools.write_imd_header(output1,config1.box,rc,atom_type_sym2num)
      SelectAtoms1={}
      for j in config1.id:
        tmp={j:0}
        SelectAtoms1.update(tmp)
      for j in config1.id:
        if i == atom_type_num2sym[config1.type[j]]:
          tmp={j:1}
          SelectAtoms1.update(tmp)
          for k in config1_nnl.nnl[j]:
            tmp={k:1}
            SelectAtoms1.update(tmp)
      for j in config1.id:
        if SelectAtoms1[j] == 1:
          ganisetti_tools.write_imd_atom(output1,j,config1,config1_nnl)
      output1.close()

  # for Bio Glass; SiO2 + P2O5 + Na2O + CaO
  if BG_glass == "yes":
    output1=open(BASE_FILE+str("_")+str("NaCa")+str("-anions.atoms"),'w')
    ganisetti_tools.write_imd_header(output1,config1.box,rc,atom_type_sym2num)
    SelectAtoms1={}
    for j in config1.id:
      tmp={j:0}
      SelectAtoms1.update(tmp)
    for j in config1.id:
      if atom_type_num2sym[config1.type[j]] == "Na" or atom_type_num2sym[config1.type[j]] == "Ca":
        tmp={j:1}
        SelectAtoms1.update(tmp)
        for k in config1_nnl.nnl[j]:
          tmp={k:1}
          SelectAtoms1.update(tmp)
    for j in config1.id:
      if SelectAtoms1[j] == 1:
        ganisetti_tools.write_imd_atom(output1,j,config1,config1_nnl)
    output1.close()

    output1=open(BASE_FILE+str("_")+str("F")+str("-cations.atoms"),'w')
    ganisetti_tools.write_imd_header(output1,config1.box,rc,atom_type_sym2num)
    SelectAtoms1={}
    for j in config1.id:
      tmp={j:0}
      SelectAtoms1.update(tmp)
    for j in config1.id:
      if atom_type_num2sym[config1.type[j]] == "F":
        tmp={j:1}
        SelectAtoms1.update(tmp)
        for k in config1_nnl.nnl[j]:
          tmp={k:1}
          SelectAtoms1.update(tmp)
    for j in config1.id:
      if SelectAtoms1[j] == 1:
        ganisetti_tools.write_imd_atom(output1,j,config1,config1_nnl)
    output1.close()

  # **************************************************************************************
  # compute coordination of each atom and average atom type
  config1_coord=ganisetti_tools.compute_coordination(config1,config1_nnl,atom_type_num2sym)
  for i in atom_type_sym2num.keys():
    output1=open(BASE_FILE+str("_")+str(i)+str("_coord"),'w')
    output1.write("# %s : coord  number_of_atoms  %% of atoms\n" %(str(i)))
    for j in range(config1_nnl.max_nnl_each_atom_type_sym[i]+1):
      temp1=config1_coord.individual_coord_sym[(i,j)]
      output1.write("%d\t%d\t\t%.2lf\n" %(j,temp1,float(temp1*100.0/total_atoms_of_type_sym[i])))
    output1.write("\n\n# Average coord = %.2lf\n" %(config1_coord.avg_coord_sym[i]))
    output1.write("\n\n# Coordination of %s with each atom type\n" %(i))
    for j in atom_type_sym2num.keys():
      if config1_coord.avg_coord_of_each_type[(i,j)] != 0:
        output1.write("%s   %.2lf\n" %(j,config1_coord.avg_coord_of_each_type[(i,j)]))
    output1.close()
  output1=open(BASE_FILE+str("_average_coord"),'w')
  output1.write("# AtomType_sym  AtomType_num   Number_of_Atoms AverageCoordination\n")
  for i in atom_type_sym2num.keys():
    output1.write("%s\t%d\t%d\t%.2lf\n" %(i,atom_type_sym2num[i],total_atoms_of_type_sym[i],config1_coord.avg_coord_sym[i]))
  output1.close()
  # config1_coord.avg_coord_sym={Si:4,Al:4.13,P:4}
  # config1_coord.avg_coord_num={2:4.0,3:4.13,4:4.0}
  # config1_coord.avg_coord_of_each_type={(Si,O):4,(P,O):4}
  # config1_coord.individual_coord_sym={(Si,4):200,(Al,4):150,(Al,5):10}

  # **************************************************************************************
  # compute environment of each atom
  config1_env=ganisetti_tools.compute_each_atom_environment(config1,config1_nnl,atom_type_sym2num,atom_type_num2sym)
  
  for i in given_anions_sym2num.keys():
    output1=open(BASE_FILE+str("_")+str(i)+str("_local_env.data"),'w')
    output1.write("# %s Environment \t\t\t Number ( %s O atoms )\t\t" %(str(i),str("%")))
    for j in given_cations_sym2num.keys():
      output1.write("%s %s-%s bonds \t" %(str("%"),j,i))
    output1.write("\n")
    for j in config1_env.all_possible_local_env_and_counts.keys():
      if j[0] == i:
        temp1=0
        temp2=[]
        number_of_tabs_for_formatting_in_the_output_file=0
        for k in j:
          if temp1 > 0:
            if k[1] != 0:
              temp2.append(k)
          else:
            temp1=temp1+1
        for k in temp2:
          output1.write("%s " %(str(k)))
          number_of_tabs_for_formatting_in_the_output_file=number_of_tabs_for_formatting_in_the_output_file+1
        temp3=config1_env.all_possible_local_env_and_counts[j]
        for num_tabs in range(5-number_of_tabs_for_formatting_in_the_output_file):
          output1.write("\t")
        # writing the percentage of anions
        output1.write(" => \t%d \t%.2lf \t\t" %(temp3,temp3*100.0/total_atoms_of_type_sym[i]))
        # writing the percentage of bonds for each type of atom
        for temp_cations_in_all_possible_local_env in j:
            if temp_cations_in_all_possible_local_env[0] not in given_anions_sym2num.keys():
              #if temp_all_cations == temp_cations_in_all_possible_local_env[0]:
              temp4=temp_cations_in_all_possible_local_env[0]
              temp5=temp_cations_in_all_possible_local_env[1]
              output1.write("%.2lf \t\t" %(temp3*100.0*temp5/(total_atoms_of_type_sym[temp4]*config1_coord.avg_coord_sym[temp4])))
        output1.write("\n")

    output1.close()
  #config1_env.all_possible_local_env_and_counts={(O,(Si,1),(Al,1)):154} # total O which is in Si-O-Al are 154
  #config1_env.env_atomid[0]={('Si',4):1,('Al',4):1}

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
    output1.write("# where n = number of bridging anions\n")
    output1.write("# n  number_of_units %_of_units\n")
    for j in range(5):
      output1.write("%d\t%d\t\t%.2lf \n" % (j, config1_Q.Q_summary[(i,j)], config1_Q.Q_summary[(i,j)] * 100.0 / total_atoms_of_type_sym[i]))
      #print(config1_Q.Q_non4CoordFormers_list[i,j])
    output1.close()
  #for i in config1.id:
  #  if config1.type[i] in given_formers_sym2num.values():
  #    print(config1.type[i],config1_Q.Q_status_atomid2num[i])
  #config1_Q.Q_summary={(Si,1):10,(Si,2):20,(Si,3):100,(Si,4):150} # For Si: Q1=10, Q2=20, Q3=100, Q4=150
  #config1_Q.Q_status_atomid2num={20:2,21:-1,22:1} # atom20:type=Si,BridgingAnions=2;atom21:type=O
  #config1_Q.Q_4CoordFormers_list={(Si,1):1,2,3,4} # list of Si with Q1=1,2,3,4
  #config1_Q.Q_non4CoordFormers_list={(Si,1):55,83} # list of Si with Q1 but one of the other Si is having non 4 coord = 55, 83

  # **************************************************************************************
  # compute triplets
  config1_triplets=ganisetti_tools.compute_triplets(cmd,config1,config1_nnl)

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
  #config1_triplets.triplets_AmBCn2count={(Si,4,O,Si,4):150} # Si4-O-Si4 = 150 triplets
  #config1_triplets.total_triplets_sym2count=total triplets

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
  output1=open(BASE_FILE+str("_BA_NBA_info1.data"),'w')
  output1.write("# Bridging and Non-Bridging Aninons distribution among the network formers\n")

  NBA_count = {}
  for i in given_anions_sym2num.keys():
    for j in given_formers_sym2num.keys():
      temp1 = {(i, j): 0}
      NBA_count.update(temp1)
  for i in config1_anions_distribution.NBA_former_id2list.keys():
    i_type = atom_type_num2sym[config1.type[i]]
    j = config1_anions_distribution.NBA_former_id2list[i]
    for k in j:
      k_type = atom_type_num2sym[config1.type[k]]
      temp1 = {(i_type, k_type): NBA_count[(i_type, k_type)] + 1}
      NBA_count.update(temp1)

  for k in given_anions_sym2num.keys():
    output1.write("#########################\n")
    output1.write("# Anion %s distribution \n" %(k))
    output1.write("#########################\n")
    output1.write("total %s               = %d\n" %(k,total_atoms_of_type_sym[k]))
    output1.write("total BA              = %d\n" %(config1_anions_distribution.total_anions_of_BA_4CoordFormer_sym2count[k]))
    output1.write("total NBA             = %d\n" %(config1_anions_distribution.total_anions_of_NBA_former_sym2count[k]+config1_anions_distribution.total_anions_of_BA_non4CoordFormer_sym2count[k]))
    output1.write("total triclusters     = %d\n" %(config1_anions_distribution.total_anions_of_tri_clusters_sym2count[k]))
    output1.write("anions with zero NF   = %d\n" %(config1_anions_distribution.total_anions_with_zero_formers_sym2count[k]))
    output1.write("any other type anions = %d\n" %(config1_anions_distribution.total_anions_of_any_other_type_sym2count[k]))
    output1.write("\n# Oxygen distribution in triplets \n")
    set_a={}
    set_c={}
    count=0
    for i in given_formers_sym2num.keys():
      temp1={count:i}
      set_a.update(temp1)
      set_c.update(temp1)
      count=count+1
    go_to_next_line = "no"
    for i in range(count):
      A=set_a[i]
      for j in range(i,count):
        C=set_c[j]
        for m in range(config1_nnl.max_nnl_each_atom_type_sym[A] + 1):
          for n in range(config1_nnl.max_nnl_each_atom_type_sym[C] + 1):
            temp1=config1_anions_distribution.triplets_AmBCn2count[(A,m,k,C,n)]
            if temp1 != 0:
              output1.write("%s %d %s %s %d   = %d \n" %(A,m,k,C,n,temp1))
              go_to_next_line="yes"
        if go_to_next_line == "yes":
          output1.write("\n")
          go_to_next_line="no"
      if go_to_next_line == "yes":
        output1.write("\n")
        go_to_next_line = "no"

    output1.write("\n# Non-bridging anions distribution\n")
    for j in given_formers_sym2num.keys():
      if NBA_count[(k,j)] != 0:
        output1.write("%s - %s %d\n" %(j,k,NBA_count[(k,j)]))
    output1.write("\n\n")

  # for Bio-Glass; SiO2 + P2O5 + Na2O + CaO
  if BG_glass == "yes":
    Si_typ=atom_type_sym2num["Si"]
    P_typ=atom_type_sym2num["P"]
    SelectAtoms1={}
    for i in config1.id:
      if Si_typ in config1_nnl.nnl_type_num[i] and P_typ in config1_nnl.nnl_type_num[i]:
        temp1={i:i}
        SelectAtoms1.update(temp1)
        for j in config1_nnl.nnl[i]:
          if Si_typ == config1.type[j] or P_typ == config1.type[j]:
            temp1={j:j}
            SelectAtoms1.update(temp1)
    output1=open(BASE_FILE+str("_SiOP.atoms"),'w')
    ganisetti_tools.write_imd_header(output1,config1.box,rc,atom_type_sym2num)
    for i in SelectAtoms1.keys():
      ganisetti_tools.write_imd_atom(output1,i,config1,config1_nnl)
    output1.close()

  # For NAPS glass
  if NAPS_glass == "yes":
    Na_type=atom_type_sym2num["Na"]
    # Na speciation
    total_BA_count=0
    total_NBA_count=0
    Na_BO_pair_distance=0.0
    Na_NBO_pair_distance=0.0
    Na_BA={}
    Na_NBA={}
    Na_descretization_in_BA  = {}
    Na_descretization_in_NBA = {}
    for Si_count in range(10):
      for Al_count in range(10):
        for P_count in range(10):
          temp1={(Si_count,Al_count,P_count):0}
          Na_descretization_in_BA.update(temp1)
          Na_descretization_in_NBA.update(temp1)
    for i in config1.id:
      if config1.type[i] == Na_type:
        BA_count=0
        NBA_count=0
        for j in config1_nnl.nnl[i]: # j=anion
          if j in config1_anions_distribution.BA_4CoordFormer_id2list.keys():
            BA_count=BA_count+1
            if ('Si',4) in config1_env.env_atomid[j]:
              Si_count=config1_env.env_atomid[j][('Si',4)]
            else:
              Si_count=0
            if ('Al',4) in config1_env.env_atomid[j]:
              Al_count=config1_env.env_atomid[j][('Al',4)]
            else:
              Al_count=0
            if ('P',4) in config1_env.env_atomid[j]:
              P_count=config1_env.env_atomid[j][('P',4)]
            else:
              P_count=0
            temp1={(Si_count,Al_count,P_count):Na_descretization_in_BA[(Si_count,Al_count,P_count)]+1}
            Na_descretization_in_BA.update(temp1)
          else :
            NBA_count=NBA_count+1
            Si_count=0
            Al_count=0
            P_count=0
            Si_count=config1_nnl.nnl_type_sym[j].count('Si')
            Al_count=config1_nnl.nnl_type_sym[j].count('Al')
            P_count=config1_nnl.nnl_type_sym[j].count('P')
            temp1={(Si_count,Al_count,P_count):Na_descretization_in_NBA[(Si_count,Al_count,P_count)]+1}
            Na_descretization_in_NBA.update(temp1)
        temp1={i:BA_count}
        temp2={i:NBA_count}
        Na_BA.update(temp1)
        Na_NBA.update(temp2)
        total_BA_count=total_BA_count+BA_count
        total_NBA_count=total_NBA_count+NBA_count
        # computing average bond length of Na-BO and Na-NBO
        count=0
        for j in config1_nnl.nnl[i]:
          if j in config1_anions_distribution.BA_4CoordFormer_id2list.keys():
      	    Na_BO_pair_distance=Na_BO_pair_distance+config1_nnl.nnl_each_pair_distance[i][count]
          else:
      	    Na_NBO_pair_distance=Na_NBO_pair_distance+config1_nnl.nnl_each_pair_distance[i][count]
          count=count+1
    output1=open(BASE_FILE+str("_BA_and_NBA_around_Na_atoms1.data"),'w')
    output1.write("# BA and NBA attached to Na\n")
    output1.write("# Na-Anion(type)\tnum_of_bonds\t%_of_Na\n")
    output1.write("BA\t\t\t= %d\t\t%.2f\n" %(total_BA_count,total_BA_count*100.0/(total_BA_count+total_NBA_count)))
    output1.write("NBA\t\t\t= %d\t\t%.2f\n" %(total_NBA_count,total_NBA_count*100.0/(total_BA_count+total_NBA_count)))
    output1.write("\n\n")
    output1.write("Average bond length of Na-BA  = %lf\n" %(Na_BO_pair_distance/total_BA_count))
    output1.write("Average bond length of Na-NBA = %lf\n" %(Na_NBO_pair_distance/total_NBA_count))
    output1.write("\n\n")
    output1.write("\n\n# Detailed speciation\n")
    output1.write("# Na-Anion(type)\tnum_of_bonds\t%_of_Na\n")
    for Si_count in range(3):
      for Al_count in range(3):
        for P_count in range(3):
          temp3=Na_descretization_in_BA[(Si_count,Al_count,P_count)]
          if temp3 !=0:
            temp1=[0 for i in range(2)]
            temp2=0
            if Si_count == 1:
              temp1[temp2]="Si"
              temp2=temp2+1
            elif Si_count == 2:
              temp1[0]="Si"
              temp1[1]="Si"
            if Al_count == 1:
              temp1[temp2]="Al"
              temp2=temp2+1
            elif Al_count == 2:
              temp1[0]="Al"
              temp1[1]="Al"
            if P_count == 1:
              temp1[temp2]="P"
              temp2=temp2+1
            elif P_count == 2:
              temp1[0] = "P"
              temp1[1] = "P"
            output1.write("Na in %s - BA - %s\t= %d\t\t%.2f\n" %(str(temp1[0]),str(temp1[1]),temp3,temp3*100.0/(total_BA_count+total_NBA_count)))
    other_type=0
    output1.write("\n")
    for Si_count in range(10):
      for Al_count in range(10):
        for P_count in range(10):
    #for i in Na_descretization_in_NBA.keys():
    #  Si_count=i[0]
    #  Al_count=i[1]
    #  P_count=i[2]
          temp3=Na_descretization_in_NBA[(Si_count,Al_count,P_count)]
          if temp3 != 0:
            if Si_count+Al_count+P_count == 1:
              if Si_count == 1:
                output1.write("Na in Si - NBA\t\t= %d\t\t%.2f\n" % (temp3, temp3 * 100.0 / (total_BA_count + total_NBA_count)))
              elif Al_count == 1:
                output1.write("Na in Al - NBA\t\t= %d\t\t%.2f\n" % (temp3, temp3 * 100.0 / (total_BA_count + total_NBA_count)))
              elif P_count == 1:
                output1.write("Na in P  - NBA\t\t= %d\t\t%.2f\n" % (temp3, temp3 * 100.0 / (total_BA_count + total_NBA_count)))
            else:
              other_type=other_type+temp3
    output1.write("Na in tri-clusters or over coordinated formers \t= %d\t%.2f\n" % (other_type, other_type * 100.0 / (total_BA_count + total_NBA_count)))
    output1.write("\n\n\n# atoms id \t Number_Of_BA\tNumber_of_NBA\n")
    for i in config1.id:
      if config1.type[i] == Na_type:
        output1.write("%d  => %d\t%d\n" %(i,Na_BA[i],Na_NBA[i]))
    output1.close()

    # For contour density plot of Na in BA and NBA
    Na_BA_and_NBA={}
    output1=open(BASE_FILE+str("_BA_and_NBA_around_Na_atoms2.data"),'w')
    output1.write("#BO\tNBO\tNumber_Of_Na_atoms\n")
    #for i in range(max(Na_BA.values())+1):
    #  for j in range(max(Na_NBA.values())+1):
    for i in range(11):
      for j in range(11):
        temp1={(i,j):0}
        Na_BA_and_NBA.update(temp1)

    for i in Na_BA.keys():
      temp1={(Na_BA[i],Na_NBA[i]):Na_BA_and_NBA[(Na_BA[i],Na_NBA[i])]+1}
      Na_BA_and_NBA.update(temp1)
    # preparing data for plot
    x=[]
    y=[]
    z=[]
    for i in Na_BA_and_NBA.keys():
      x.append(i[0])
      y.append(i[1])
      z.append(0.5*Na_BA_and_NBA[i]**2)
      output1.write("%d\t%d\t%d\n" %(i[0],i[1],Na_BA_and_NBA[i]))
    # preparing fig
    fig1=plt.figure()
    ax1= fig1.add_subplot(111)
    plt.scatter(x,y,s=z,alpha=0.5)
    ax1.set_xlabel("BO")
    ax1.set_ylabel("NBO")
    plt.savefig('Fig01_Na_in_BO_and_NBO_density_plot.png',bbox_inches='tight',dpi=600)
    output1.close()
    
    # Seperate the Na and BO, and Na and NBO atoms
    output1=open(BASE_FILE+str("_BA_and_NBA_around_Na__Na_BO_atoms.atoms"),'w')
    output2=open(BASE_FILE+str("_BA_and_NBA_around_Na__Na_NBO_atoms.atoms"),'w')
    output3=open(BASE_FILE+str("_BA_and_NBA_around_Na__BO_atoms_no_Na.atoms"),'w')
    ganisetti_tools.write_imd_header(output1,config1.box,rc,atom_type_sym2num)
    ganisetti_tools.write_imd_header(output2,config1.box,rc,atom_type_sym2num)
    ganisetti_tools.write_imd_header(output3,config1.box,rc,atom_type_sym2num)
    Na_in_BO__Na_O_atoms=[]
    Na_in_NBO__Na_O_atoms=[]
    BO_without_Na__O_atoms=[]
    for anion_id in config1_anions_distribution.modifiers_id2list.keys():
      if anion_id in config1_anions_distribution.BA_4CoordFormer_id2list.keys():
      	if len(config1_anions_distribution.modifiers_id2list[anion_id]) > 0: 	# the anion_id has atleast one modifier
      	  Na_in_BO__Na_O_atoms.append(anion_id)
      	  for j in config1_anions_distribution.modifiers_id2list[anion_id]:
      	    Na_in_BO__Na_O_atoms.append(j)
      	else:																	# the anion_id does not have any modifier
      	  BO_without_Na__O_atoms.append(anion_id)
      else :
      	Na_in_NBO__Na_O_atoms.append(anion_id)
      	for j in config1_anions_distribution.modifiers_id2list[anion_id]:
      	  Na_in_NBO__Na_O_atoms.append(j)
    Na_in_BO__Na_O_atoms=list(set(Na_in_BO__Na_O_atoms))
    Na_in_NBO__Na_O_atoms=list(set(Na_in_NBO__Na_O_atoms))
    BO_without_Na__O_atoms=list(set(BO_without_Na__O_atoms))
    for i in Na_in_BO__Na_O_atoms:
      ganisetti_tools.write_imd_atom(output1,i,config1,config1_nnl)
    for i in Na_in_NBO__Na_O_atoms:
      ganisetti_tools.write_imd_atom(output2,i,config1,config1_nnl)
    for i in BO_without_Na__O_atoms:
      ganisetti_tools.write_imd_atom(output3,i,config1,config1_nnl)
    output1.close()
    output2.close()
    output3.close()

    # All Na atoms and their neighbors, if neighbor O is BO assign type 11 or if neighbor O is NBO assign type 12 
    output1=open(BASE_FILE+str("_BA_and_NBA_around_Na__Na_BO_NBO_atoms.atoms"),'w')
    output2=open(BASE_FILE+str("_BO_and_Formers.atoms"),'w')
    ganisetti_tools.write_imd_header(output1,config1.box,rc,atom_type_sym2num)
    ganisetti_tools.write_imd_header(output2,config1.box,rc,atom_type_sym2num)
    SelectAtoms1={}
    for atom_id in config1.id:
      if config1.type[atom_id] == atom_type_sym2num["Na"]:
      	output1.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_id,config1.type[atom_id],config1.posx[atom_id],config1.posy[atom_id],config1.posz[atom_id],config1_nnl.nnl_count[atom_id]))
      	for Na_neighbour in config1_nnl.nnl[atom_id]:
      	  if Na_neighbour in config1_anions_distribution.BA_4CoordFormer_id2list.keys():
      	  	output1.write("%d  %d  %lf  %lf  %lf %d \n" %(Na_neighbour,11,config1.posx[Na_neighbour],config1.posy[Na_neighbour],config1.posz[Na_neighbour],config1_nnl.nnl_count[Na_neighbour]))
      	  else :
      	  	output1.write("%d  %d  %lf  %lf  %lf %d \n" %(Na_neighbour,12,config1.posx[Na_neighbour],config1.posy[Na_neighbour],config1.posz[Na_neighbour],config1_nnl.nnl_count[Na_neighbour]))
      #elif config1.type[atom_id] == atom_type_sym2num["Si"] or config1.type[atom_id] == atom_type_sym2num["Al"] or config1.type[atom_id] == atom_type_sym2num["P"]:
      #	temp1={atom_id:atom_id}
      #	SelectAtoms1.update(temp1)
      #	for i in config1_nnl.nnl[atom_id]:
      #	  temp1={i:i}
      #	  SelectAtoms1.update(temp1)
    #for i in SelectAtoms1.keys():
    #  ganisetti_tools.write_imd_atom(output2,i,config1,config1_nnl)
    #print(config1_anions_distribution.BA_4CoordFormer_id2list.keys())
    for anion_id in config1_anions_distribution.BA_4CoordFormer_id2list.keys():
      temp1={anion_id:anion_id}
      SelectAtoms1.update(temp1)
      for anion_neighbor in config1_anions_distribution.BA_4CoordFormer_id2list[anion_id]:
      	temp1={anion_neighbor:anion_neighbor}
      	SelectAtoms1.update(temp1)
    for i in SelectAtoms1.keys():
      ganisetti_tools.write_imd_atom(output2,i,config1,config1_nnl)
    output1.close()
    output2.close()
    
    # Q(Si) distribution
    output1=open(BASE_FILE+str("_Q")+str("Si.data"),'a+')
    output1.write("\n\n")
    output1.write("# Q Discretization\n")
    output1.write("# Qa(bAl)(cP); where a= total_network_cations; b=num_of_Al; c=num_of_P\n")
    output1.write("# a\tb\tc\tnum_of_units\t%_of_units\n")
    temp_Q=[[[0 for i in range(5)] for j in range(5)] for k in range(5)]
    for i in config1.id:
      if atom_type_num2sym[config1.type[i]] == "Si" and config1_nnl.nnl_count[i] == 4:
        Q_count = config1_Q.Q_status_atomid2num[i]
        Al_count=0
        P_count=0
        for j in config1_nnl.nnl[i]:
          if ('Al',4) in config1_env.env_atomid[j].keys():
            if config1_env.env_atomid[j][('Al',4)] == 1:
              Al_count = Al_count+1
          if ('P',4) in config1_env.env_atomid[j].keys():
            if config1_env.env_atomid[j][('P', 4)] == 1:
              P_count = P_count+1
        temp_Q[Q_count][Al_count][P_count]=temp_Q[Q_count][Al_count][P_count]+1
    for i in range(5):
      for j in range(5):
        for k in range(5):
          temp1=temp_Q[i][j][k]
          if temp1 != 0:
            output1.write("%d\t%d\t%d\t%d\t%.2f\n" %(i,j,k,temp1,temp1*100.0/total_atoms_of_type_sym["Si"]))

    # Q(P) distribtuion
    output1=open(BASE_FILE+str("_Q")+str("P.data"),'a+')
    output1.write("\n\n# Q Discretization\n")
    temp1_Q=[[0 for i in range(5)] for j in range(5)]
    temp2_Q=[[[0 for i in range(5)] for j in range(5)] for k in range(5)]
    temp3_Q=[0 for i in range(5)]
    for i in config1.id:
      if atom_type_num2sym[config1.type[i]] == "P" and config1_nnl.nnl_count[i] == 4:
        P_count=0
        Al_count=0
        Si_count=0
        for j in config1_nnl.nnl[i]:
          P_count=P_count+config1_env.env_atomid[j][('P',4)]-1
          if ('Al', 4) in config1_env.env_atomid[j].keys():
            if config1_env.env_atomid[j][('Al', 4)] == 1:
              Al_count = Al_count + 1
          if ('Si', 4) in config1_env.env_atomid[j].keys():
            if config1_env.env_atomid[j][('Si', 4)] == 1:
              Si_count = Si_count + 1
        temp3_Q[P_count]=temp3_Q[P_count]+1
        temp1_Q[P_count][Al_count]=temp1_Q[P_count][Al_count]+1
        temp2_Q[P_count][Al_count][Si_count]=temp2_Q[P_count][Al_count][Si_count]+1
    output1.write("# Q(aP); where a= num_of_P\n")
    output1.write("# a\tnum_of_units\t%_of_units\n")
    for i in range(5):
      temp3=temp3_Q[i]
      if temp3 !=0:
        output1.write("%d\t%d\t\t%.2f\n" %(i,temp3,temp3*100.0/total_atoms_of_type_sym["P"]))
    output1.write("\n\n# Q(aP)(bAl); where a= num_of_P; b=num_of_Al\n")
    output1.write("# a\tb\tnum_of_units\t%_of_units\n")
    for i in range(5):
      for j in range(5):
        temp1=temp1_Q[i][j]
        if temp1 != 0:
          output1.write("%d\t%d\t%d\t\t%.2f\n" %(i,j,temp1,temp1*100.0/total_atoms_of_type_sym["P"]))
    output1.write("\n\n# Q(aP)(bAl)(cSi); where a= num_of_P; b=num_of_Al; c=num_of_Si\n")
    output1.write("# a\tb\tc\tnum_of_units\t%_of_units\n")
    for i in range(5):
      for j in range(5):
        for k in range(5):
          temp2=temp2_Q[i][j][k]
          if temp2 != 0:
            output1.write("%d\t%d\t%d\t%d\t\t%.2f\n" %(i,j,k,temp2,temp2*100.0/total_atoms_of_type_sym["P"]))

  # For CMAS glass
  if CMAS_glass == "yes":
    #temp5={}
    #for i in range(5):
    #  temp1={i:0}
    #  temp5.update(temp1)
    #for i in config1.id:
    #  if config1.type[i] == atom_type_sym2num["Al"] and config1_nnl.nnl_count[i] == 5:      # i = 5 coord Al atoms
    #    for j in config1_nnl.nnl[i]:                                                        # j = Oxygen atoms
    #      # print(config1_env.env_atomid[j])
    #      temp1={}
    #      temp3=0
    #      for k in config1_env.env_atomid[j].keys():
    #        if list(k)[0] in given_formers_sym2num.keys():
    #          temp2={k:config1_env.env_atomid[j][k]}
    #          temp1.update(temp2)
    #          temp3=temp3+config1_env.env_atomid[j][k]
    #      temp4={temp3:temp5[temp3]+1}
    #      temp5.update(temp4)
    #      print(temp1)
    #for i in range(4):
    #  print(i,temp5[i])
    output1 = open(BASE_FILE + str("_") + str("Al") + str("_coord"), 'a+')
    output1.write("\n# 5 coordinated Al that are associated with tri-clusters and their percentage among total 5 coordinated Al\n")
    temp1=0
    for i in config1.id:
      if config1.type[i] == atom_type_sym2num["Al"] and config1_nnl.nnl_count[i] == 5:  # i = 5 coord Al atoms
        temp2=0
        for j in config1_nnl.nnl[i]:                                                    # j = Oxygen atoms
          if j in config1_anions_distribution.tri_cluster_former_id2list.keys():
            temp2=1
            # print(config1_env.env_atomid[j])
        if temp2 == 1:
          temp1 = temp1 + 1
          # print(config1_env.env_atomid[i])
    output1.write("%d\t%.2lf" %(temp1,temp1*100.0/config1_coord.individual_coord_sym[("Al",5)]))
    output1.close()


