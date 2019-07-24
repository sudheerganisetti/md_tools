# Sudheer Ganisetti  
Tue Nov  3 09:31:14 CET 2015 

compile  
********  
gcc ChangingAtomTypeOfSelectedAtomsInAChkpt.c -o ChangingAtomTypeOfSelectedAtoms  
or   
gcc ChangingAtomTypeOfSelectedAtomsInAChkpt.c -lm -o ChangingAtomTypeOfSelectedAtoms  


run  
****  
ChangingAtomTypeOfSelectedAtoms ChkptFile.chkpt SelectedAtoms.data  


Info  
****  
You have to provide two files along with it....  
1) A chkpt file, where you want to change the atom type  
2) A .nnl or .txt or .dat any file with atom ids of the atoms you want to change their id  
SO, if the atom is Si the old type is 0 now changes to 2 for the selected atoms  
    if the atom is O the old type is 1 now changes to 3 for the selected atoms  

