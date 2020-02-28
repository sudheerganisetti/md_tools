
"""
@Author : Sudheer Ganisetti
@Date   : Mo 25. Dez 18:34:13 CET 2017
@         Do 20. Jun 20:34:10 CEST 2019
@         Fr 28. Feb 16:05:21 CET 2020
This code is to compute the local density of given atoms based on voxel descritization
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

  # controllers for some specific glasses
  temp1=atom_type_sym2num.keys()
  BG_glass="no"
  if "Si" in temp1 and "P" in temp1 and "Na" in temp1 and "Ca" in temp1:
    BG_glass="yes"
  NAPS_glass="no"
  if "P" in temp1 and "Al" in temp1 and "Na" in temp1:
    NAPS_glass="yes"

  RAS_glass="no"
  if "Si" in temp1 and "Al" in temp1:
    if "Mg" in temp1 or "Ca" in temp1 or "Sr" in temp1:
      RAS_glass="yes"


  # main loop starts here
  BASE_FILE=str(sys.argv[1])
  LAMMPS_DUMP_FILE=BASE_FILE+str('.dump')
  MAX_NEIGHBOURS=8
  MAX_ATOM_TYPES=8

  # **************************************************************************************
  # get the atoms information
  config1	= ganisetti_tools.get_atoms_info_from_lammps(LAMMPS_DUMP_FILE)

# compute voxel based atoms density
all_density_required_atoms=[]
for i in config1.id:
  if config1.type[i] == atom_type_sym2num["P"]:
  	all_density_required_atoms.append(i)

config1_atoms_density=ganisetti_tools.compute_atoms_density(config1,0.8,3,all_density_required_atoms)
output1=open("temporary.txt",'w')
ganisetti_tools.write_imd_header_custom_property(output1,config1.box,"local_atoms_density")
"""
This is depricated
output1.write("#F A 1 1 3 1 0\n")
output1.write("#C number type x y z voxel_based_atoms_density \n")
output1.write("#X %lf 0.0 0.0 \n" %(config1.box[0][1]-config1.box[0][0]))
output1.write("#Y 0.0 %lf 0.0 \n" %(config1.box[1][1]-config1.box[1][0]))
output1.write("#Z 0.0 0.0 %lf \n" %(config1.box[2][1]-config1.box[2][0]))
output1.write("## atoms density sudheer ganisetti's own idea\n")
output1.write("#E \n")
"""
count=1
for i,j,k in config1_atoms_density.atoms_density.keys():
  output1.write("%d 1 %.4lf %.4lf %.4lf %d\n" %(count,float(i),float(j),float(k),config1_atoms_density.atoms_density[(i,j,k)]))
  count=count+1
output1.close()

