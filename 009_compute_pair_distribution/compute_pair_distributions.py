
"""
@Author : Sudheer Ganisetti
@Date   : Mo 22. Jul 17:08:50 CEST 2019
This code computes the distribution of cation-anion pairs
"""

import numpy as np
import sys
import math as mt
import itertools as it
import datetime
import ganisetti_tools

if __name__=="__main__":

  # read command line first
  cmd = ganisetti_tools.read_command_line(sys.argv)
  if cmd.error != "none":
    sys.exit(0)
  atom_type_sym2num = cmd.atom_type_sym2num
  atom_type_num2sym = cmd.atom_type_num2sym
  bond_length_sym2num = cmd.bond_length_sys2num
  bond_length_num2num = cmd.bond_length_num2num
  rc = cmd.rc
  given_anions_sym2num = cmd.given_anions_sym2num
  given_cations_sym2num = cmd.given_cations_sym2num
  given_formers_sym2num = cmd.given_formers_sym2num
  given_modifiers_sym2num = cmd.given_modifiers_sym2num

  # main loop starts here
  BASE_FILE = str(sys.argv[1])
  LAMMPS_DUMP_FILE = BASE_FILE + str('.dump')

  # **************************************************************************************
  # get the atoms information
  config1 = ganisetti_tools.get_atoms_info_from_lammps(LAMMPS_DUMP_FILE)
  rcut=7.0
  total_bins=100
  config1_pdf=ganisetti_tools.compute_pair_distribution(config1,cmd,rcut,total_bins)

  for i in given_cations_sym2num.keys():
    for j in given_anions_sym2num.keys():
      output1=open(BASE_FILE+str("_")+str(i)+str(j)+str(".pair"),'w')
      output1.write("# S.No  r    g(r) \n")
      for k in range(total_bins):
        output1.write("%d     %lf\t\t%lf\n" %(k,float(k*rcut/total_bins),config1_pdf.gr[(i,j,k)]))
      output1.close()

