
"""
@Author : Sudheer Ganisetti
@Date   : Mo 25. Dez 18:34:13 CET 2017
@         Do 20. Jun 20:34:10 CEST 2019
This code is to compute nnl, nnlchk,  connectivity(Qn) and A-O-A triplets for the sample (SiO2 + Al2O3 + P2O5 + Na2O + CaO + MgO + CaF2)
"""

import sys
import ganisetti_tools

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
    output1.write("# %s Environment \t\t\t Number ( Percentage )\n" %(str(i)))
    for j in config1_env.all_possible_local_env_and_counts.keys():
      if j[0] == i:
        temp1=0
        temp2=[]
        for k in j:
          if temp1 > 0:
            if k[1] != 0:
              temp2.append(k)
          else:
            temp1=temp1+1
        for k in temp2:
          output1.write("%s " %(str(k)))
        temp3=config1_env.all_possible_local_env_and_counts[j]
        output1.write("\t\t\t=> %d ( %.2lf )\n" %(temp3,temp3*100.0/total_atoms_of_type_sym[i]))
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
    output1.write("# n  number_of_units\n")
    for j in range(5):
      output1.write("%d %d\t%.2lf \n" % (j, config1_Q.Q_summary[(i,j)], config1_Q.Q_summary[(i,j)] * 100.0 / total_atoms_of_type_sym[i]))
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

  temp1=given_formers_sym2num.keys()

  #if "Na"  in temp1 and "Al" in temp1 and "P" in temp1 and "Si" in temp1:
  #  print(given_formers_sym2num.keys())