#! /usr/bin/python
"""
Author    : Sudheer Ganisetti
Date      : Fr 15. Dez 15:02:18 CET 2017
            Compute density
"""
import numpy as np
import math
import sys
import random as rd
import subprocess
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if __name__=="__main__":

  # Checking for command line and correct usage  
  error=0
  error1=0
  error2=0
  error3=0
  error4=0
  tot_argv=len(sys.argv)
  if tot_argv < 3 or (tot_argv%2) !=0:
    error1=1
    error=1
  
  if error==0:
    commandline=[]
    for i in sys.argv:
      commandline.append(i)
    commandline=list(commandline)
    given_all_atoms=[]
    given_all_atom_types=[]
    given_O_type=-1
    given_Si_type=-1
    given_Al_type=-1
    given_Ca_type=-1
    given_Mg_type=-1
    given_Sr_type=-1
    if "-O" in commandline:
      O_index=commandline.index("-O")
      given_O_type=int(commandline[O_index+1])
      given_all_atoms.append("-O")
      given_all_atom_types.append(given_O_type)
    if "-Si" in commandline:
      Si_index=commandline.index("-Si")
      given_Si_type=int(commandline[Si_index+1])
      given_all_atoms.append("-Si")
      given_all_atom_types.append(given_Si_type)
    if "-Al" in commandline:
      Al_index=commandline.index("-Al")
      given_Al_type=int(commandline[Al_index+1])
      given_all_atoms.append("-Al")
      given_all_atom_types.append(given_Al_type)
    if "-Ca" in commandline:
      Ca_index=commandline.index("-Ca")
      given_Ca_type=int(commandline[Ca_index+1])
      given_all_atoms.append("-Ca")
      given_all_atom_types.append(given_Ca_type)
    if "-Mg" in commandline:
      Mg_index=commandline.index("-Mg")
      given_Mg_type=int(commandline[Mg_index+1])
      given_all_atoms.append("-Mg")
      given_all_atom_types.append(given_Mg_type)
    if "-Sr" in commandline:
      Sr_index=commandline.index("-Sr")
      given_Sr_type=int(commandline[Sr_index+1])
      given_all_atoms.append("-Sr")
      given_all_atom_types.append(given_Sr_type)

    atoms_we_have=set(["-O","-Si","-Al","-Ca","-Mg","-Sr"])
    given_all_atoms=set(given_all_atoms)
    if len(given_all_atoms-atoms_we_have) !=0:
      error2=1
      error=1
    if len(given_all_atoms) != (len(commandline)-2)/2:
      error3=1
      error=1

  if error ==0:
    File1=str(sys.argv[1])
    AtoCM=1e-8					# Angstrom to cm conversion factor
    CMtoA=1e8					# cm to Angstrom conversion factor
    AtomsPerMole=6.022140857*1e23 		# Avagadro's number
    factor=AtomsPerMole/(CMtoA*CMtoA*CMtoA)	# to convert from Moles and A into g/cm3
    m_Si =28.0855			# g/mole
    m_O  =15.999			# g/mole
    m_Al =26.981539		# g/mole
    m_Ca =40.078			# g/mole
    m_Mg =24.305			# g/mole   
    m_Sr =87.62			# g/mole
    MAX_TYPE=max(given_all_atom_types)
    mass_of_each_atom_of_type=[0.0 for i in range(0,MAX_TYPE+1)]
    if given_O_type != -1:
      mass_of_each_atom_of_type[given_O_type]=m_O
    if given_Si_type != -1:
      mass_of_each_atom_of_type[given_Si_type]=m_Si
    if given_Al_type != -1:
      mass_of_each_atom_of_type[given_Al_type]=m_Al
    if given_Ca_type != -1:
      mass_of_each_atom_of_type[given_Ca_type]=m_Ca
    if given_Mg_type != -1:
      mass_of_each_atom_of_type[given_Mg_type]=m_Mg
    if given_Sr_type != -1:
      mass_of_each_atom_of_type[given_Sr_type]=m_Sr

    data1=open(File1,'r')
    # Reading box header
    line_count=0
    for line in data1:
      if line_count == 5:
        data2=line.strip().split()
        xlo=float(data2[0])
        xhi=float(data2[1])
      if line_count == 6:
        data2=line.strip().split()
        ylo=float(data2[0])
        yhi=float(data2[1])
      if line_count == 7:
        data2=line.strip().split()
        zlo=float(data2[0])
        zhi=float(data2[1])
      line_count=line_count+1

    box_xx_length=xhi-xlo
    box_yy_length=yhi-ylo
    box_zz_length=zhi-zlo
    Volume=box_xx_length*box_yy_length*box_zz_length

    # Reading data for atoms type
    data3=np.loadtxt(File1,skiprows=9)
    # if atom types given in the command line are not equal to the atom types in the file
    if set(data3.T[1]) != set(given_all_atom_types):
      error4=1
      error=1

  if error ==0:
    number_of_atoms_of_type=[0 for i in xrange(0,MAX_TYPE+1)]
    for i in data3.T[1]:
      number_of_atoms_of_type[int(i)]=number_of_atoms_of_type[int(i)]+1

    #mass_of_type=[0.0 for i in xrange(0,MAX_TYPES)]
    total_mass=0.0
    for j in range(0,MAX_TYPE+1):
      temp1=number_of_atoms_of_type[j]*mass_of_each_atom_of_type[j]
      total_mass=total_mass+temp1

    # compute density
    density=total_mass/Volume/factor
    print density

    # sanitary check
    # If the composition contains SiO2 + Al2O3 + CaO + MgO + SrO
    temp1=0
    if given_Si_type != -1:
      temp1=temp1+number_of_atoms_of_type[given_Si_type]*2.0			# SiO2
    if given_Al_type != -1:
      temp1=temp1+number_of_atoms_of_type[given_Al_type]*3.0/2.0		# Al2O3
    if given_Ca_type != -1:
      temp1=temp1+number_of_atoms_of_type[given_Ca_type]			# CaO
    if given_Mg_type != -1:
      temp1=temp1+number_of_atoms_of_type[given_Mg_type]			# MgO
    if given_Sr_type != -1:
      temp1=temp1+number_of_atoms_of_type[given_Sr_type]			# SrO

    if temp1 != number_of_atoms_of_type[given_O_type]:
      print "--------------------------------------------------------------------------------------------------------------------------------------"
      print "| WARNING: It is assumed that the sample contains few of the following commponents                                                   |"
      print "|          SiO2 + Al2O3 + CaO + MgO + SrO                                                                                            |"
      print "|          if that is the case, the anions from the file are not equal to the anions calculated from the cations using above formula |"
      print "|          Please check it again! If you know what you are doing dont bother about this warning!                                     |"
      print "--------------------------------------------------------------------------------------------------------------------------------------"

  if error == 1 or error1 == 1 or error2 == 1 or error3 == 1 or error4 == 1:
     subprocess.call("sudheer_banner")
  if error ==1:
    print ""
    print "************* S. Ganisetti *************"
    print "Computing the density of lammps data file"
    print "correct usage: ./this_program  lammps.dump -O 1 -Si 2 -Al 3 -Ca 4 -Mg 5 ..."
    print "check carefully that the atom types you have in the file is same as the atom types you are "
    print "paasing to the program through the command line above"
    print "The data is available for O,Si,Al,Ca,Mg,Sr"
    print "Error could be due to the following reason(s)"
    if error1==1:
      print "  arguments passing to the program through command line are not enough"
    if error2==1:
      print "  some of the atoms you have passed through the command line are not available in the database"
    if error3==1:
      print "  some of the atoms you have passed through the command line might be wrong"
    if error4==1:
      print "  additional atom types are found in the file which you did not pass through the command line"
    print "****************************************"
    print ""
    sys.exit(0)

