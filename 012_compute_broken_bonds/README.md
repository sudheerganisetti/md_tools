# Sudheer Ganisetti

gcc compute_nnl_chkpt.c -o nnl_chkpt  
usage: ./nnl_diff_chkpt chkptfile 1 chkptfile2 nnlfile1 nnlfile 2  


output:  
you will get a ".nnlchkpt" for "chkpt2" with including additional information of   
tot_nn         = total number of nearest neighbours for each atom  
broken_bonds   = bonds present in chkpt 1 but not in chkpt 2  
new_bonds      = bonds present in chkpt 2 but not in chjpt 1  
switched bonds = no of bonds broken in chkpt 1 and at the same time formed new bonds in chkpt 2  


compute_nnl_chkpt_AdditionalColumnBondStatus.c  
according to Erik => bond status will be last one and which means  
 * 0 - no changes in bond  
 * 1 - new bonds were formed  
 * 2 - no of bonds broken = no of bonds formed  
 * 3 - bonds were broken  

05/10/2015  
when i change my code for getting bond status in the above order (Acc Erik's suggestion), I forget to enable the bond_switched column  
so bond status is getting correctly but switched_bonds are not written into the nnlchkpt file  
Now, I enable it  
Everything should work fine!  
So far, everything is working fine  

15/10/2015  
The switch bond logic is not correct to count switched bonds !!!!!!  
Till now I am calculating the switch bonds by considering as follows  
If a bond is broken between atom A and atom B and formed a new bond between atom A and atom C then I am counting as switched bonds for atom A = 1  
But this logic is misleading in some cases  
So trying to implementing it in a way to avoid this misleading in the next folder  
(However, the broken bond, new bond information is correct !!!!!)  

02/12/2015  
Now I am increasing the neighbors list checking count from 6 neighbours to 8 neighbours  
i.e If the number of nearest neighbours are more in the new configuration ==> New Bond (It doesnt matter, whether bonds are switched or not)  
    If the number of nearest neighbours are less in the new configuration ==> Broken Bond (It doesnt matter, whether bonds are switched or not)  
    If the number of nearest neighbours are same but neighbouring atoms are not same ==> Switched Bond  
    
    
    broken_bonds    = number of broken bonds    (it does not depend on the bond_status)  
    new_bonds       = number of new bonds       (it does not depend on the bond_status)  
    switched_bonds  = number of switched bonds  (it does not depend on the bond_status)  
    bond_status     = I didnt change  
     * 0 - no changes in bond  
     * 1 - new bonds were formed  
     * 2 - no of bonds broken = no of bonds formed  
     * 3 - bonds were broken  


