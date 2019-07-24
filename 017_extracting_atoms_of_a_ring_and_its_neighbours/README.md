#Sudheer Ganisetti  
Mon Mar  7 15:34:55 CET 2016  
*************************************  
The atom numbers from the RingsDataList are fake atom numbers i.e based on the atom position in xyz file (but they are not the real atom numbers of IMD)  
So what I am doing here is  
1) Read the fake identification numbers from the RingsDataList  
2) with this we can find what is the position (i.e line number) in the xyz file  
3) with this we can find what the position (i.e line number) in IMD file  
4) get the atom number from that position (i.e from that line)  

In this script, you can control number of lines (i.e number of rings) need to reed from the RingsDataFile and also Size of the Ring  

./collecting_atom_numbers_for_chhkpt_files_based_on_rings_list_from_rings_code.bsh RingsDataList Chkptfile.chkpt  

We can run a single bash script to do all these jobs because for each task I have written a code separately
./run_little_pieces_of_codes_together_and_collect_the_atoms_in_chkpt_format.bsh RingsDataList Chkptfile.chkpt nnlfile.nnl

A better script can be written if time is given, unfortunately, this is not in my high priority list of tasks
