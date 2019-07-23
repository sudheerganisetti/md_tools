#                                     Author: Sudheer Ganisetti


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


The following tools are very useful for analysing various properties of amorphous materials (glasses) with given atomic positions and chemical atom types.
The main intension of these programs is to post-process the output files of the atomistic simulations (LAMMPS and IMD programs).

1) 001_ganisetti_tools_module			: a module with a lot of classes and functions for computing various properties

2) 002_compute_glass_density 			: computes the density of a glass

3) 003_CMAS_glass_structural_properties		: computes the structural properties of CMAS glass

4) 004_NAPS_glass_structural_properties		: computes the structural properties of NAPS glass

5) 005_SiO2_Al2O3_RO_TernaryPhaseDiagram	: plots the ternary phase diagram of SiO2, Al2O3 and RO

6) 006_Compute_MeanSquareDisplacement		: computes the mean square displacement of each atom type

7) 007_StructureGeneration			: generate the structure of desired composition by placing the atoms randomly 

8) 008_structural_properties_for_all_glasses	: a generalized programm for computing structural properties of any kind of glass

9) 009_compute_pair_distribution		: computes pair distribution function of any glass composition (still need to normalize)

NOTE: IMD input or output file formats will not work for all the programs, however, this can be solved in near future.

