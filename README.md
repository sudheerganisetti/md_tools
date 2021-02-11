
          #***********************************************************************************# 
          #                                                                                   # 
          #    \\\|///      //////  /     /  /////    /     /  ///////  ///////  /////        # 
          #    (  °  )      /       /     /  /    /   /     /  /        /        /    /       # 
          #   ( ^   ^ )     /       /     /  /     /  /     /  /        /        /     /      # 
          #  {  °   °  }    //////  /     /  /     /  ///////  //////   //////   //////       # 
          #   (  \^/  )          /  /     /  /     /  /     /  /        /        /    /       # 
          #    (  Ö  )           /  /     /  /    /   /     /  /        /        /     /      # 
          #     (_Ä_)       //////  ///////  /////    /     /  ///////  ///////  /     /      # 
          #                                                                                   # 
          #***********************************************************************************# 

The following figure is prepared with the data obtained by post-processing tools written under the framework 'md_tools'. The main intention for presenting the figure is to show only the visual capabilities rather than discussing the scientific findings. Therefore, you have to excuse me for not providing complete details of the figure. The tools list and their short description are presented following the figure.

![all_figures_v02](https://user-images.githubusercontent.com/7426968/107627483-a58b6100-6c5f-11eb-85e7-b12abc114505.png)
Figure description:
a) CMAS glass (11096 atoms) with the assignment of various colors to various atom types
b) Diffusion path of Na in NAP glass
c) Distribution of Li and K in Vanadium containing glasses
d) Brding Oxygen (BO) rich network region of NAPS glass
e) Non-Bridging Oxygen (NBO) rich Channel region of NAPS glass
f) Demonstrating the role of Na in the presence of Si
g) Separating the channel and network regions based on BO and NBO for NAPS glass
h) Identification of clusters in CMAS glass
i) Computing the presence of Na at various oxgen sites (NAPS glasses)
j) Percentage of Al-O-Al bonds in CAS glasses 
k) Pair distribution functions of CMAS glass
l) Distribution of O at various sites in NAPS glasses
m) Percentages of O in Al-O-V and V-O-V in Vanadium containing glasses
n) Extraction of triclusters (CMAS glass)
o) Influence of F on the network structure of Bio-glasses
(The figures are generated with [ovito](https://www.ovito.org) and [gnuplot](http://www.gnuplot.info/) software packages)

#                                     Author: Sudheer Ganisetti

The following tools are very useful for analysing various properties of amorphous materials (glasses) with given atomic positions and chemical atom types.
The main intension of these programs is to post-process the output files obtained by the atomistic simulations (LAMMPS and IMD programs).

01) 001_ganisetti_tools_module			: a module with a lot of classes and functions for computing various properties  
02) 002_compute_glass_density 			: computes the density of a glass  
03) 003_CMAS_glass_structural_properties	: computes the structural properties of CMAS glass  
04) 004_NAPS_glass_structural_properties	: computes the structural properties of NAPS glass  
05) 005_SiO2_Al2O3_RO_TernaryPhaseDiagram	: plots the ternary phase diagram of SiO2, Al2O3 and RO  
06) 006_Compute_MeanSquareDisplacement		: computes the mean square displacement of each atom type  
07) 007_StructureGeneration			: generate the structure of desired composition by placing the atoms randomly  
08) 008_structural_properties_for_all_glasses	: a generalized programm for computing structural properties of any kind of glass  
09) 009_compute_pair_distribution		: computes pair distribution function of any glass composition (still need to normalize)
10) 010_IMD_projection_of_pair_distribution	: computes projection of pair distribution function on x,y, and z-axis
11) 011_IMD_projection_of_angle_distribution    : computes azimuthal and polar angles from the angle distributions
12) 012_compute_broken_bonds			: computes broken bonds when given two files contains same number of atoms
13) 013_change_atom_type_of_selected_atoms	: changes the atom type of selected atoms
14) 014_ReadingAtomNumbersFromAFileAndCollectingThem : three different codes exist for reading atoms from a file and collect only these atoms from: (1) imd chkpt file (2) lammps dump file (3) nnl file
15) 015_compute_common_atoms_of_two_data_files	: get the common atoms of two data files (for example two ring data files)
16) 016_separating_each_tetrahedron		: separate each tetrahedron
17) 017_extracting_atoms_of_a_ring_and_its_neighbours : extract the atoms of a given ring and also their neighbors from a given chkpt file
18) 018_script_to_extract_edge_corner_tets_info : extract the edge and corner sharing tetrahedra information from a rings_code output file
19) 019_making_SLICES_compute_Guttman_rings	: make the sample into slices and compute the ring statistics with Gutman rings criterion using rings_code
20) 020_mean_square_displacement		: compute mean square displacement
21) 021_bond_orientations_alpha_beta		: compute orientation of each bond in terms of alpha and beta in spherical coordinates
22) 022_compute_local_atomic_density_based_on_voxels : the local atomic density has been calculated for the given atoms from a lammps dum file which can be useful to analyze the channel regions formed by the given atoms 
23) 023_compute_the_diffusion_path_of_each_given_atom_during_MSD: Computes the diffusion path of each given atom during MSD
24) 024_extract_given_list_of_atoms_from_the_given_chkpt_file : Extract the list of atoms provided by the atom numbers in a text file from the given chkpt file
25) 025_checking_for_network_glitches_using_ping  : A small code to report the network glitches
26) 026_compute_channels_clustering_based_on_voxels: Computes the channel regions provided by the type of atoms that the channel contains. For example, extract the channels containing the network modifier atoms 
27) 027_compute_bond_statistics                   : Computes broken bonds, switched bonds and new bonds during the MSD
28) 028_ContourPlot_of_ModifierCation_Connected_to_BO_and_NBO: Generates a nice illustration (contour plot) of percentage of Na atoms accommodated at various sites based on the bridging oxygen (BO) and non-bridging oxygen (NBO).

NOTE: IMD input or output file formats will not work for all the programs, however, this can be solved in near future.

