# Sudheer Ganisetti  

/* *******************************************************************************************************************
 * Sudheer Ganisetti, Mo 21. Mär 21:58:04 CET 2016
 * A small piece of code to calculate the Mean Square Displacement (MSD)
 * How I am calculating the MSD
 * 1) Take reference sample and make it into voxels (of given size i.e VoxelSizeX,VoxelSizeY, VoxelSizeZ)
 * 2) Take the Current sample
 * 3) Deal with periodic boundary conditions
 * 4) calculate SD of each atom
 * 5) MSD = Average all SD's in a reference voxel
 * 6) Output SD and MSD to CurChkpt
 *********************************************************************************************************************/

gcc making_into_small_boxes_to_measure_MeanSquareDisplacement.c  -o MeanSquareDisplacement_WithPBC -lm


