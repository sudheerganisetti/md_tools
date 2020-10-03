"""
@Author : Sudheer Ganisetti
@Date   : Mo 25. Dez 18:34:13 CET 2017
@         Do 20. Jun 20:34:10 CEST 2019
@         Fr 28. Feb 16:05:21 CET 2020
@         Mi 30. Sep 23:23:46 CEST 2020
This code is to compute the clustered channels of given atoms based on voxel descritization
"""

import sys
import ganisetti_tools
import time
import datetime

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if __name__=="__main__":

  # read command line first
  cmd=ganisetti_tools.read_command_line(sys.argv)
  if cmd.error != "none":
    sys.exit(0)
  # get the version details
  ganisetti_tools.get_ganisetti_tools_version()
  start_time=time.time()
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
  IMD_FILE=BASE_FILE+str('.imd')
  MAX_NEIGHBOURS=8
  MAX_ATOM_TYPES=8

  # controlling parameter
  file_format               = "lammps"
  voxel_length_parameter    = 0.5 # for example Si-O bond length is 2.2 Ang, and you want to keep 5 voxels within the Si-O bond 
                                  # then the voxel_length_parameter you have to provide is 2.2/5 = 0.44
  voxel_smoothing_parameter = 6   # if you want to search for neighboring atoms within a separation distance of 2.7 Ang 
                                  # then the voxel_smoothing_parameter you have to provide is int(2.7/voxel_length_parameter) = int(6.1364) = 6

  # **************************************************************************************
  # get the atoms information
  if file_format == "imd":
    config1 = ganisetti_tools.get_atoms_info_from_imd(IMD_FILE)
  if file_format == "lammps":
    config1 = ganisetti_tools.get_atoms_info_from_lammps(LAMMPS_DUMP_FILE)


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
  
  #all_clustering_required_atoms=[1759,586,383,524,237,2649,363,465,93,1868]
  all_clustering_required_atoms=[]
  for i in config1.id:
    if config1.type[i] == atom_type_sym2num["Na"]:
      all_clustering_required_atoms.append(i)
      for j in config1_nnl.nnl[i]:
        all_clustering_required_atoms.append(j)
  
  print("total Na atoms = %d " %(len(all_clustering_required_atoms)))
  
  config1_Na_channels=ganisetti_tools.compute_clustered_channels(config1,config1_nnl.nnl,0.8,5,all_clustering_required_atoms)
  
  output1=open(BASE_FILE+str("_clusters.atoms"),'w')
  ganisetti_tools.write_imd_header_custom_property(output1,config1.box,"cluster_id")
  count=1
  
  for i,j,k in config1_Na_channels.cluster_cells_position_to_id.keys():
    output1.write("%d 1 %.4lf %.4lf %.4lf %d\n" %(count,float(i),float(j),float(k),config1_Na_channels.cluster_cells_position_to_id[(i,j,k)]))
    count=count+1
  output1.close()
  
end_time=time.time()
print("Total excution time = %s" %(str(datetime.timedelta(seconds=end_time-start_time))))
