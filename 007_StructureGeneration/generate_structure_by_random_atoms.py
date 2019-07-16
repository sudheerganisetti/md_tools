#! /usr/bin/python
"""
Author    : Sudheer Ganisetti
Date      : Mi 13. Dez 15:05:30 CET 2017
            Generate random atoms
"""
import numpy as np
import math
import sys
import random as rd
import subprocess
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if __name__=="__main__":
  
  tot_argv=len(sys.argv)
  if tot_argv != 2:
     #subprocess.call("sudheer_banner")
     print "************* S. Ganisetti *************"
     print "Error: usage"
     print "./this_program  rcut"
     print "Generate random atoms"
     print "****************************************"
     sys.exit(0)

  rcut=float(sys.argv[1])  
  AtoCM=1e-8			# Angstrom to cm conversion factor
  CMtoA=1e8			# cm to Angstrom conversion factor
  AtomsPerMole=6.022140857*1e23 # Avagadro's number
  m_Si =28.0855			# g/mole
  m_O  =15.999			# g/mole
  m_Al =26.981539		# g/mole
  m_Ca =40.078			# g/mole
  m_Mg =24.305			# g/mole
  m_Sr =87.62			# g/mole
  m_P  =30.973762		# g/mole
  m_B  =10.811	      		# g/mole
  m_Na =22.98976		# g/mole
  m_Li =6.941			# g/mole
     
  one_atom_mass_Si =1.0*m_Si/AtomsPerMole	# grams
  one_atom_mass_O  =1.0*m_O/AtomsPerMole	# grams
  one_atom_mass_Al =1.0*m_Al/AtomsPerMole	# grams
  one_atom_mass_Ca =1.0*m_Ca/AtomsPerMole	# grams
  one_atom_mass_Mg =1.0*m_Mg/AtomsPerMole	# grams
  one_atom_mass_Sr =1.0*m_Sr/AtomsPerMole	# grams
  one_atom_mass_P  =1.0*m_P/AtomsPerMole	# grams
  one_atom_mass_B  =1.0*m_B/AtomsPerMole	# grams
  one_atom_mass_Na =1.0*m_Na/AtomsPerMole	# grams
  one_atom_mass_Li =1.0*m_Li/AtomsPerMole	# grams
  
  SiO2_unit_mass   =1.0*one_atom_mass_Si + 2.0*one_atom_mass_O		# grams
  Al2O3_unit_mass  =2.0*one_atom_mass_Al + 3.0*one_atom_mass_O		# grams
  CaO_unit_mass    =1.0*one_atom_mass_Ca + 1.0*one_atom_mass_O		# grams
  MgO_unit_mass    =1.0*one_atom_mass_Mg + 1.0*one_atom_mass_O		# grams
  SrO_unit_mass    =1.0*one_atom_mass_Sr + 1.0*one_atom_mass_O	        # grams
  P2O5_unit_mass   =2.0*one_atom_mass_P  + 5.0*one_atom_mass_O		# grams
  B2O3_unit_mass   =2.0*one_atom_mass_B  + 3.0*one_atom_mass_O		# grams
  Na2O_unit_mass   =2.0*one_atom_mass_Na + 1.0*one_atom_mass_O		# grams
  Li2O_unit_mass   =2.0*one_atom_mass_Li + 1.0*one_atom_mass_O        	# grams

  # Controlling parameters to compute box length based on moles of composition and density
  density	=2.56          			# g/cm3  # calculated this density for 2000 K roughly calculated from Courtail&Dingwell1999
  SiO2_moles	=0
  Al2O3_moles	=0
  CaO_moles	=0
  MgO_moles	=0
  SrO_moles	=0
  P2O5_moles	=0
  B2O3_moles	=0
  Na2O_moles	=0
  Li2O_moles	=0
  # give atom types in the order as you like
  # However, when using the Pedone2006 potential the parameters are coded by giving the following order
  # O=1 Si=2 Al=3 Ca=4 Mg=5 P=6 Na=7 K=8 Sr=9 Li=10
  given_O_atom_type	=1
  given_Si_atom_type	=2
  given_Al_atom_type	=3
  given_Ca_atom_type	=4
  given_Mg_atom_type	=5
  given_Sr_atom_type	=9
  given_P_atom_type	=6
  given_B_atom_type	=-1
  given_Na_atom_type	=7
  given_Li_atom_type	=10

  given_O_atom_charge     =-1.2
  given_Si_atom_charge    =2.4
  given_Al_atom_charge    =1.8
  given_Ca_atom_charge    =1.2
  given_Mg_atom_charge    =1.2
  given_Sr_atom_charge    =1.2
  given_P_atom_charge     =3.0
  given_B_atom_charge     =0.0
  given_Na_atom_charge    =0.6
  given_Li_atom_charge    =0.6

  max_aec	=9	#buffer for maximum atoms in each cell
   
  mass_SiO2  = SiO2_unit_mass	* SiO2_moles	# grams
  mass_Al2O3 = Al2O3_unit_mass	* Al2O3_moles	# grams
  mass_CaO   = CaO_unit_mass	* CaO_moles	# grams
  mass_MgO   = MgO_unit_mass	* MgO_moles	# grams
  mass_SrO   = SrO_unit_mass	* SrO_moles	# grams
  mass_P2O5  = P2O5_unit_mass	* P2O5_moles	# grams
  mass_B2O3  = B2O3_unit_mass	* B2O3_moles	# grams
  mass_Na2O  = Na2O_unit_mass	* Na2O_moles	# grams
  mass_Li2O  = Li2O_unit_mass   * Li2O_moles    # grams
  total_mass =mass_SiO2+mass_Al2O3+mass_CaO+mass_MgO+mass_SrO +mass_P2O5+mass_B2O3+mass_Na2O+mass_Li2O	# grams
  volume=total_mass*CMtoA*CMtoA*CMtoA/density	# A^3
  box_xx=np.cbrt(volume)
  box_yy=box_xx
  box_zz=box_xx
  
  total_atoms=SiO2_moles*3+Al2O3_moles*5+CaO_moles*2+MgO_moles*2+SrO_moles*2+P2O5_moles*7+B2O3_moles*5+Na2O_moles*3+Li2O_moles*3
  atom_id		=[-1 for i in xrange(0,total_atoms)]
  atom_type_xyz		=["O" for i in xrange(0,total_atoms)]
  atom_type_lammps	=[1 for i in xrange(0,total_atoms)]
  atom_charge		=[given_O_atom_charge for i in xrange(0,total_atoms)]
  reserved_atom		=[0 for i in xrange(0,total_atoms)]
  atom_posx		=[-1.0 for i in xrange(0,total_atoms)]
  atom_posy		=[-1.0 for i in xrange(0,total_atoms)]
  atom_posz		=[-1.0 for i in xrange(0,total_atoms)]

  # Now generate atoms randomly with in the cubic cell of box_xx
  # to keep the atoms not close to each other, check the distance of generated atom with neighbouring atoms
  # and avoid if the generated atom is close to any other atom with a distance less than rcut radius
  # devide the box into voxels (cells) based on rcut
  no_cells_X=int(box_xx/rcut)
  no_cells_Y=int(box_yy/rcut)
  no_cells_Z=int(box_zz/rcut)

  cell_width_X=box_xx/no_cells_X
  cell_width_Y=box_yy/no_cells_Y
  cell_width_Z=box_zz/no_cells_Z

  atoms_list_in_cell=[[[[-1 for i in xrange(0,max_aec)] for j in xrange(0,no_cells_Z)] for k in xrange(0,no_cells_Y)] for l in xrange(0,no_cells_X)]
  atoms_count_in_cell=[[[0 for j in xrange(0,no_cells_Z)] for k in xrange(0,no_cells_Y)] for l in xrange(0,no_cells_X)]
  atoms_count=0
  while atoms_count < total_atoms:
    rand_x=rd.uniform(0.0,box_xx)
    rand_y=rd.uniform(0.0,box_yy)
    rand_z=rd.uniform(0.0,box_zz)
    px_cell=int(rand_x/cell_width_X)
    py_cell=int(rand_y/cell_width_Y)
    pz_cell=int(rand_z/cell_width_Z)
    accept_the_atom="yes"
    for i in (-1,0,1):
      for j in (-1,0,1):
        for k in (-1,0,1):
          new_px_cell=px_cell+i
          new_py_cell=py_cell+j
          new_pz_cell=pz_cell+k
          # consider periodic boundary conditions
          if new_px_cell < 0:
            new_px_cell=no_cells_X-1
          if new_py_cell < 0:
            new_py_cell=no_cells_Y-1
          if new_pz_cell < 0:
            new_pz_cell=no_cells_Z-1

          if new_px_cell > no_cells_X-1:
            new_px_cell=0
          if new_py_cell > no_cells_Y-1:
            new_py_cell=0
          if new_pz_cell > no_cells_Z-1:
            new_pz_cell=0
          for l in atoms_list_in_cell[new_px_cell][new_py_cell][new_pz_cell]:
            if l != -1:
              dx=atom_posx[l]-rand_x
              dy=atom_posy[l]-rand_y
              dz=atom_posz[l]-rand_z
              r=math.sqrt(dx*dx+dy*dy+dz*dz)
              if r < rcut:
                accept_the_atom="no"

    if accept_the_atom == "yes":
      atom_id[atoms_count]=atoms_count
      atom_posx[atoms_count]=rand_x
      atom_posy[atoms_count]=rand_y
      atom_posz[atoms_count]=rand_z
      temp1=atoms_count_in_cell[px_cell][py_cell][pz_cell]
      atoms_list_in_cell[px_cell][py_cell][pz_cell][temp1]=atoms_count
      atoms_count_in_cell[px_cell][py_cell][pz_cell]=atoms_count_in_cell[px_cell][py_cell][pz_cell]+1
      atoms_count=atoms_count+1

    #atom_posx.append(rand_x)
    #atom_posy.append(rand_y)
    #atom_posz.append(rand_z)
    
  # Randomly pic atoms and give name as Si,Al,Ca,Mg,Sr,P,B,Na and Li
  count=0
  # Si
  while count < SiO2_moles:
    rand=rd.randint(0,total_atoms-1)
    if reserved_atom[rand] == 0:
      atom_type_xyz[rand]	="Si"
      atom_type_lammps[rand]	=given_Si_atom_type
      atom_charge[rand]		=given_Si_atom_charge
      reserved_atom[rand]	=1
      count=count+1
  count=0
  # Al
  while count < 2*Al2O3_moles:
    rand=rd.randint(0,total_atoms-1)
    if reserved_atom[rand] == 0:
      atom_type_xyz[rand]	="Al"
      atom_type_lammps[rand]	=given_Al_atom_type
      atom_charge[rand]		=given_Al_atom_charge
      reserved_atom[rand]	=1
      count=count+1
  count=0
  # Ca
  while count < CaO_moles:
    rand=rd.randint(0,total_atoms-1)
    if reserved_atom[rand] == 0:
      atom_type_xyz[rand]	="Ca"
      atom_type_lammps[rand]	=given_Ca_atom_type
      atom_charge[rand]		=given_Ca_atom_charge
      reserved_atom[rand]	=1
      count=count+1
  count=0
  # Mg
  while count < MgO_moles:
    rand=rd.randint(0,total_atoms-1)
    if reserved_atom[rand] == 0:
      atom_type_xyz[rand]	="Mg"
      atom_type_lammps[rand]	=given_Mg_atom_type
      atom_charge[rand]		=given_Mg_atom_charge
      reserved_atom[rand]	=1
      count=count+1
  count=0
  # Sr
  while count < SrO_moles:
    rand=rd.randint(0,total_atoms-1)
    if reserved_atom[rand] == 0:
      atom_type_xyz[rand]	="Sr"
      atom_type_lammps[rand]	=given_Sr_atom_type
      atom_charge[rand]		=given_Sr_atom_charge
      reserved_atom[rand]	=1
      count=count+1
  count=0
  # P
  while count < 2*P2O5_moles:
    rand=rd.randint(0,total_atoms-1)
    if reserved_atom[rand] == 0:
      atom_type_xyz[rand]	="P"
      atom_type_lammps[rand]	=given_P_atom_type
      atom_charge[rand]		=given_P_atom_charge
      reserved_atom[rand]	=1
      count=count+1
  count=0
  # B
  while count < 2*B2O3_moles:
    rand=rd.randint(0,total_atoms-1)
    if reserved_atom[rand] == 0:
      atom_type_xyz[rand]	="B"
      atom_type_lammps[rand]	=given_B_atom_type
      atom_charge[rand]		=given_B_atom_charge
      reserved_atom[rand]	=1
      count=count+1
  count=0
  # Na
  while count < 2*Na2O_moles:
    rand=rd.randint(0,total_atoms-1)
    if reserved_atom[rand] == 0:
      atom_type_xyz[rand]	="Na"
      atom_type_lammps[rand]	=given_Na_atom_type
      atom_charge[rand]		=given_Na_atom_charge
      reserved_atom[rand]	=1
      count=count+1
  count=0
  # Li
  while count < 2*Li2O_moles:
    rand=rd.randint(0,total_atoms-1)
    if reserved_atom[rand] == 0:
      atom_type_xyz[rand]	="Li"
      atom_type_lammps[rand]	=given_Li_atom_type
      atom_charge[rand]		=given_Li_atom_charge
      reserved_atom[rand]	=1
      count=count+1

  # write atoms in xyz format to a file
  Output1=open("RandomlyDistributedAtoms.xyz",'w')
  Output1.write("%d\n" %(total_atoms))
  Output1.write("# %d SiO2 + %d Al2O3 + %d CaO + %d MgO +%d SrO + %d P2O5 + %d B2O3 + %d Na2O + %d Li2O \n" %(SiO2_moles,Al2O3_moles,CaO_moles,MgO_moles,SrO_moles,P2O5_moles,B2O3_moles,Na2O_moles,Li2O_moles))
  for i in range(0,total_atoms):
      Output1.write("%s  %lf  %lf  %lf\n" %(atom_type_xyz[i],atom_posx[i],atom_posy[i],atom_posz[i]))

  total_atom_types=len(set(atom_type_lammps))

  # write atoms in lammps format
  Output2=open("RandomlyDistributedAtoms.lammps",'w')
  #Output2.write("# SudheerGanisetti %d SiO2 + %d Al2O3 + %d CaO + %d MgO +%d SrO + %d P2O5 + %d B2O3 + %d Na2O,id  type  charge  x  y  z  vx  vy  vz\n\n" %(SiO2_moles,Al2O3_moles,CaO_moles,MgO_moles,SrO_moles,P2O5_moles,B2O3_moles,Na2O_moles))
  Output2.write("# Sudheer, S.No  id      charge  x       y       z       vx      vy      vz\n\n")
  Output2.write("%d atoms\n" %(total_atoms))
  Output2.write("%d atom types\n" %(total_atom_types))
  Output2.write("0.0 %lf xlo xhi\n" %(box_xx))
  Output2.write("0.0 %lf ylo yhi\n" %(box_yy))
  Output2.write("0.0 %lf zlo zhi\n\n" %(box_zz))
  Output2.write("Atoms\n\n")
  for i in range(0,total_atoms):
    Output2.write("%d %d %.01lf %lf %lf %lf 0 0 0\n" %(i+1,atom_type_lammps[i],atom_charge[i],atom_posx[i],atom_posy[i],atom_posz[i]))

  Output1.close()
  Output2.close()
