#! /usr/bin/python
"""
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* This module is written by Sudheer Ganisetti
* These tools are useful for extracting and analyzing several glass properties 
* upon providing the MD simulations data
* usage: import ganisetti_tools
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
"""
import numpy as np
import sys
import itertools

#class banner:
def sudheer_banner():
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Sudheer banner
  * This is a simple and little fun thing created to display SUDHEER on the screen
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  print("\n\n")
  print("    #***********************************************************************************# ")
  print("    #                                                                                   # ")
  print("    #     \\\\|//       //////  /     /  /////    /     /  ///////  ///////  /////        #")
  print("    #    (  °  )      /       /     /  /    /   /     /  /        /        /    /       # ")
  print("    #   ( ^   ^ )     /       /     /  /     /  /     /  /        /        /     /      # ")
  print("    #  {  °   °  }    //////  /     /  /     /  ///////  //////   //////   //////       # ")
  print("    #   (  \^/  )          /  /     /  /     /  /     /  /        /        /    /       # ")
  print("    #    (  Ö  )           /  /     /  /    /   /     /  /        /        /     /      # ")
  print("    #     (_Ä_)       //////  ///////  /////    /     /  ///////  ///////  /     /      # ")
  print("    #                                                                                   # ")
  print("    #***********************************************************************************# ")
  print("\n\n")


class get_atoms_info_from_imd:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Reading IMD file
  * This class helps to read the imd format file
  * and store the data into an object
  * usage: config1=ganisetti_tools.get_atoms_info_from_imd(01.imd)
  * 
  * output:     config1.totalatoms         = TotalAtoms
  *             config1.box_xx             = box_xx
  *             config1.box_yy             = box_yy
  *             config1.box_zz             = box_zz
  *             config1.box                = np.array([[box_xx_lo,box_xx_hi],[box_yy_lo,box_yy_hi],[box_zz_lo,box_zz_hi]])
  *             config1.id[id]             =id
  *             config1.type[id]           =int(i[1])
  *             config1.mass[id]           =float(i[2])
  *             config1.posx[id]           =float(i[3])
  *             config1.posy[id]           =float(i[4])
  *             config1.posz[id]           =float(i[5])
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,IMD_DUMP_FILE):
    atoms_data1=np.loadtxt(IMD_DUMP_FILE)
    atoms_data2=open(IMD_DUMP_FILE,'r')
    count=1
    for each_line in atoms_data2:
      if count < 10:
        data=each_line.strip().split()
        if count == 1:
          if data[0] != "#F":
            print("#########################################################")
            print("#                                                       #")
            print("#  ERROR: It seems the given file is not in imd format  #")
            print("#                                                       #")
            print("#########################################################")
        if data[0] == "#X":
          box_xx = float(data[1])
          box_xy = float(data[2])
          box_xz = float(data[3])
        if data[0] == "#Y":
          box_yx = float(data[1])
          box_yy = float(data[2])
          box_yz = float(data[3])
        if data[0] == "#Z":
          box_zx = float(data[1])
          box_zy = float(data[2])
          box_zz = float(data[3])
      count = count + 1
   
    self.box_xx          = box_xx + min(box_yx,box_zx)
    self.box_yy          = box_yy + min(box_xy,box_zy)
    self.box_zz          = box_zz + min(box_xz,box_yz)
    self.box=np.array([[box_xx,box_xy,box_xz],[box_yx,box_yy,box_yz],[box_zx,box_zy,box_zz]])
    self.totalatoms      = len(atoms_data1)

    self.id     = {}
    self.type   = {}
    self.mass   = {}
    self.posx   = {}
    self.posy   = {}
    self.posz   = {}
    for i in atoms_data1:
      tmp_id            = int(i[0])
      self.id[tmp_id]   = tmp_id
      self.type[tmp_id] = int(i[1])
      self.mass[tmp_id] = float(i[2])
      self.posx[tmp_id] = float(i[3])
      self.posy[tmp_id] = float(i[4])
      self.posz[tmp_id] = float(i[5])


class get_atoms_info_from_lammps:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Reading LAMMPS DUMP file
  * This class helps to read the imd format file
  * and store the data into an object
  *
  * usage: config1=ganisetti_tools.get_atoms_info_from_lammps(01.dump)
  *
  * output:	config1.totalatoms      	= TotalAtoms
  *		config1.totalproperties 	= TotalProperties
  *		config1.box_xx          	= box_xx
  *		config1.box_yy          	= box_yy
  *		config1.box_zz          	= box_zz
  *		config1.box	     	= np.array([[box_xx_lo,box_xx_hi],[box_yy_lo,box_yy_hi],[box_zz_lo,box_zz_hi]])
  *    		config1.id[id]     	=id
  *    		config1.type[id]   	=int(i[1])
  *    		config1.charge[id] 	=float(i[2])
  *    		config1.posx[id]   	=float(i[3])
  *    		config1.posy[id]   	=float(i[4])
  *    		config1.posz[id]   	=float(i[5])
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,LAMMPS_DUMP_FILE):
    atoms_data1=np.loadtxt(LAMMPS_DUMP_FILE,dtype='int,int,float,float,float,float,float,float,float',skiprows=9)
    atoms_data2=open(LAMMPS_DUMP_FILE,'r')
    TotalAtoms=len(atoms_data1)
    TotalProperties=len(atoms_data1[0])
    # Reading box dimensions
    count=1
    for each_line in atoms_data2:
      if count < 10:
        data=each_line.strip().split()
        if count == 4:
          if TotalAtoms != int(data[0]):
            print("######################################################################################################################")
            print("#                                                                                                                    #")
            print("#  WARNING!: The total number of atoms mentioned in the file is not equal to the total number of atoms in the file   #")
            print("#                                                                                                                    #")
            print("######################################################################################################################")
        elif count == 6:
          box_xx_lo=float(data[0])
          box_xx_hi=float(data[1])
          box_xx=box_xx_hi-box_xx_lo
        elif count == 7:
          box_yy_lo=float(data[0])
          box_yy_hi=float(data[1])
          box_yy=box_yy_hi-box_yy_lo
        elif count == 8:
          box_zz_lo=float(data[0])
          box_zz_hi=float(data[1])
          box_zz=box_zz_hi-box_zz_lo
      count=count+1
    self.totalatoms      = TotalAtoms
    self.totalproperties = TotalProperties
    self.box_xx          = box_xx
    self.box_yy          = box_yy
    self.box_zz          = box_zz
    self.box=np.array([[box_xx_lo,box_xx_hi],[box_yy_lo,box_yy_hi],[box_zz_lo,box_zz_hi]])

    self.id     = {}
    self.type   = {}
    self.charge = {}
    self.posx   = {}
    self.posy   = {}
    self.posz   = {}
    for i in atoms_data1:
      tmp_id              =int(i[0])
      self.id[tmp_id]     =tmp_id
      self.type[tmp_id]   =int(i[1])
      self.charge[tmp_id] =float(i[2])
      self.posx[tmp_id]   =float(i[3])
      self.posy[tmp_id]   =float(i[4])
      self.posz[tmp_id]   =float(i[5])


class compute_nnl:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Compute Nearest Neighbouring List (nnl)
  * This class helps to compute the  nnl of given atoms
  * and store the data into an object
  *
  * usage: config1=ganisetti_tools.compute_nnl(config,rc,atom_type_num2sym)
  * where config		= ganisetti_tools.get_atoms_info_from_lammps(01.dump)
  *       rc[0][0]		= 2.2
  *       atom_type_num2sym	= {"1":"O", "2":"Si", "3":"Al", etc,.}
  *
  * output: config1.nnl
  *         config1.nnl_count
  *         config1.nnl_type_num
  *         config1.nnl_type_sym
  *         config1.nnl_each_pair_distance
  *         config1.max_nnl_each_atom_type_sym
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,config,rc,atom_type_num2sym):

    MAX_NEIGHBOURS=15
    TotalAtoms=config.totalatoms
    rcut=0.0
    for i in range(len(rc)):
      for j in rc[i]:
        rcut=float(max(j,rcut))
    
    # Devide the box into voxels
    MaxCells_X=int(config.box_xx/rcut)
    MaxCells_Y=int(config.box_yy/rcut)
    MaxCells_Z=int(config.box_zz/rcut)

    # LEC= Length of Each Cell
    LEC_X=float(config.box_xx/MaxCells_X)
    LEC_Y=float(config.box_yy/MaxCells_Y)
    LEC_Z=float(config.box_zz/MaxCells_Z)

    # cell_atoms_count[X][Y][Z]
    cell_atoms_count=[[[0 for i in range(MaxCells_Z)] for j in range(MaxCells_Y)] for k in range(MaxCells_X)]
    # cell_atoms_list[X][Y][Z][200]
    cell_atoms_list=[[[[0 for i in range(200)] for j in range(MaxCells_Z)] for k in range(MaxCells_Y)] for l in range(MaxCells_X)]

    # define parameters to store nnl related information
    self.nnl			            = {}
    self.nnl_count		            = {}
    self.nnl_type_num		        = {}
    self.nnl_type_sym               = {}
    self.nnl_each_pair_distance     = {}
    self.max_nnl_each_atom_type_sym = {}
    MAX_ATOM_NUMBER=max(config.id)
    nnl=[[-1 for j in range(MAX_NEIGHBOURS)] for i in range(MAX_ATOM_NUMBER+1)]
    nnl_count=[0 for i in range(MAX_ATOM_NUMBER+1)]
    #storing atom types of each neighbour in the nnl_types for getting things easy to work with
    #nnl_type=[[-1 for j in range(MAX_NEIGHBOURS)] for i in range(MAX_ATOM_NUMBER+1)]
    #store each pair distance to use them in later
    nnl_each_pair_distance=[[-1 for j in range(MAX_NEIGHBOURS)] for i in range(MAX_ATOM_NUMBER+1)]
    max1={}
    for i in atom_type_num2sym.keys():
      temp1={i:0}
      max1.update(temp1)

    # store atom positions locally
    atom_posx=config.posx
    atom_posy=config.posy
    atom_posz=config.posz
    atom_type=config.type

    # move the box to center and correspondingly atoms if the corner of the box is not at center
    for i in config.id:
      atom_posx[i]=config.posx[i]-config.box[0][0]
      atom_posy[i]=config.posy[i]-config.box[1][0]
      atom_posz[i]=config.posz[i]-config.box[2][0]
      # If the atoms are outside of the box then bring them into other side according to PBC
      if atom_posx[i] >= config.box_xx:
        atom_posx[i]=atom_posx[i]-config.box_xx
      if atom_posy[i] >= config.box_yy:
        atom_posy[i]=atom_posy[i]-config.box_yy
      if atom_posz[i] >= config.box_zz:
        atom_posz[i]=atom_posz[i]-config.box_zz
      #assign the atoms into respective cell number  #assign the atoms into respective cell number
      cellX=int(atom_posx[i]/LEC_X)
      cellY=int(atom_posy[i]/LEC_Y)
      cellZ=int(atom_posz[i]/LEC_Z)
      cell_atoms_list[cellX][cellY][cellZ][cell_atoms_count[cellX][cellY][cellZ]]=i
      cell_atoms_count[cellX][cellY][cellZ]=cell_atoms_count[cellX][cellY][cellZ]+1

    # Main loop to compute nnl and other information
    # (i,j,k) => selecting atom which is belongs to one cell
    for i in range(MaxCells_X):
      for j in range(MaxCells_Y):
        for k in range(MaxCells_Z):
          for atom1 in range(cell_atoms_count[i][j][k]):
            atom1_id=cell_atoms_list[i][j][k][atom1]
            pos_x1=atom_posx[atom1_id]
            pos_y1=atom_posy[atom1_id]
            pos_z1=atom_posz[atom1_id]

            # (l,m,n) => selecting neighbour cell 
            for l in np.arange(0,2,1):
              new_i=i+l
              an=new_i
              pbc_x=0.0
              # -----------------------------------------------------
              # Dealing with periodic boundary condition along x-axis
              if new_i < 0:
                # But the conditions will not lead to excute this case
                new_i=MaxCells_X-1
                print("WARNING!!! There is somethong wrong with the pbc in nnl, please check it well !!!\n")
              if new_i > MaxCells_X-1:
                new_i=0
                pbc_x=config.box_xx
              # -----------------------------------------------------
              for m in np.arange(-1*l,2,1):
                new_j=j+m
                bn=new_j
                pbc_y=0.0
                # -----------------------------------------------------
                # Dealing with periodic boundary condition along y-axis
                if new_j < 0:
                  new_j=MaxCells_Y-1
                  pbc_y=-1.0*config.box_yy
                if new_j > MaxCells_Y-1:
                  new_j=0
                  pbc_y=config.box_yy
                # -----------------------------------------------------
                n1=-1
                if l==m:
                  n1=-1*l
                for n in np.arange(n1,2,1):
                  new_k=k+n
                  cn=new_k
                  pbc_z=0.0
                  # -----------------------------------------------------
                  # Dealing with periodic boundary conditions along z-axis
                  if new_k < 0:
                    new_k=MaxCells_Z-1
                    pbc_z=-1.0*config.box_zz
                  if new_k > MaxCells_Z-1:
                    new_k = 0
                    pbc_z=config.box_zz
                  # debugging
                  #if i==0 and j==0 and k==0:
                  #  print ("%d %d %d %d %d %d   atom1_id %d") %(l,m,n,an,bn,cn,atom1_id)

                  # -----------------------------------------------------
                  # This is a trick to avoid double counting the interactions if the two atoms are in the same cell
                  temp=0
                  if i == new_i and j == new_j and k == new_k:
                    temp=atom1+1
                  # selecting the neighbour atom
                  for atom2 in range(temp,cell_atoms_count[new_i][new_j][new_k]):
                    atom2_id=cell_atoms_list[new_i][new_j][new_k][atom2]
                    pos_x2=atom_posx[atom2_id] + pbc_x
                    pos_y2=atom_posy[atom2_id] + pbc_y
                    pos_z2=atom_posz[atom2_id] + pbc_z

                    r=np.linalg.norm(np.array([pos_x1,pos_y1,pos_z1])-np.array([pos_x2,pos_y2,pos_z2]))
                    if r <= rc[atom_type[atom1_id]][atom_type[atom2_id]] and r > 0.2:
                      nnl[atom1_id][nnl_count[atom1_id]] = atom2_id
                      nnl[atom2_id][nnl_count[atom2_id]] = atom1_id
                      nnl_each_pair_distance[atom1_id][nnl_count[atom1_id]]=r
                      nnl_each_pair_distance[atom2_id][nnl_count[atom2_id]]=r
                      nnl_count[atom1_id]=nnl_count[atom1_id]+1
                      nnl_count[atom2_id]=nnl_count[atom2_id]+1                  

    # removing -1 from nnl[i]
    for i in range(MAX_ATOM_NUMBER+1):
      if i in config.id:
        temp1=[]
        temp2=[]
        temp3=[]
        temp4=[]
        for j in range(len(nnl[i])):
          if nnl[i][j] !=-1:
            temp1.append(nnl[i][j])
            temp2.append(atom_type[nnl[i][j]])
            temp3.append(nnl_each_pair_distance[i][j])
            temp4.append(atom_type_num2sym[atom_type[nnl[i][j]]])
        temp5={config.type[i]:int(max(max1[config.type[i]],nnl_count[i]))}
        max1.update(temp5)
        self.nnl[i]			            = temp1
        self.nnl_type_num[i]		    = temp2
        self.nnl_type_sym[i]		    = temp4
        self.nnl_each_pair_distance[i]	= temp3
        self.nnl_count[i]		        = nnl_count[i]

    # update maximum neighbors of each atom type
    temp1={}
    for i in atom_type_num2sym.keys():
      temp2={atom_type_num2sym[i]:max1[i]}
      temp1.update(temp2)
    self.max_nnl_each_atom_type_sym=temp1


def write_imd_header(output,box,rc,atom_type_sym2num):
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Writing the IMD header into given output file
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  output.write("#F A 1 1 3 1 0\n")
  output.write("#C number type x y z tot_nn \n")
  output.write("#X %lf 0.0 0.0 \n" %(box[0][1]-box[0][0]))
  output.write("#Y 0.0 %lf 0.0 \n" %(box[1][1]-box[1][0]))
  output.write("#Z 0.0 0.0 %lf \n" %(box[2][1]-box[2][0]))
  output.write("## rcut of ")
  for i in atom_type_sym2num.keys():
    if i != "O":
      output.write("%s-O = %.2lf ; " %(i,rc[atom_type_sym2num[i]][atom_type_sym2num['O']]))
  output.write("\n")
  output.write("#E \n")

def write_imd_atom(output,atom_id,config,config_nnl):
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Writing the IMD header into given output file
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  output.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_id,config.type[atom_id],config.posx[atom_id],config.posy[atom_id],config.posz[atom_id],config_nnl.nnl_count[atom_id]))



class compute_each_atom_environment:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Compute each atom's environment
  * This class helps to compute the environment of each atom
  * and store the data into an object
  *
  * usage: environment1=ganisetti_tools.compute_each_atom_environment(config,config_nnl,atom_type_sym2num,atom_type_num2sym)
  * where config                = ganisetti_tools.get_atoms_info_from_lammps(01.dump)
  *       config_nnl            = ganisetti_tools.compute_nnl(config,rc,atom_type_num2sym)
  *       atom_type_sym2num     = {"O":"1", "Si":"2", "Al":"3", etc,.}
  *       atom_type_num2sym     = {1:"O", 2:"Si", 3:"Al", etc,.}
  *
  * output: deprecated => config1.environment_atomnum2sym[0] = {'O':0,'Si':2,'Al':1}
  *         deprecated => config1.env_atomnum2nnlsymandcoord[0] = (('Si',4),('Al',4))
  *         deprecated => config1.env_atomnum2countofnnlsymandcoord[0] = (1,1)
  *         config1.env_atomid[0]={('Si',4):1,('Al',4):1}  # which means Si[4]-O-Al[4]
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,config,config_nnl,atom_type_sym2num,atom_type_num2sym):
    #self.env_atomnum2sym={}
    #self.env_atomnum2nnlsymandcoord={}
    #self.env_atomnum2countofnnlsymandcoord={}
    #self.env_total={}
    #for i in config.id:
    #  temp_each_atom_type_count={}
      #temp1={'id':i}
      #temp_each_atom_type_count.update(temp1)
    #  for j in atom_type_sym2num.keys():
    #    temp1={j:config_nnl.nnl_type_sym[i].count(j)}
    #    if temp1 != 0:
    #      temp_each_atom_type_count.update(temp1)
    #  self.env_atomnum2sym[i] = temp_each_atom_type_count
    self.env_atomid={}
    for i in config.id:
      temp_local_env = {}
      for j in config_nnl.nnl[i]:
        temp1 = (atom_type_num2sym[config.type[j]], config_nnl.nnl_count[j])
        if temp1 in temp_local_env:
          temp2 = {temp1: temp_local_env[temp1] + 1}
          temp_local_env.update(temp2)
        else:
          temp2 = {temp1: 1}
          temp_local_env.update(temp2)
      #self.env_atomnum2nnlsymandcoord[i] = tuple(temp_local_env.keys())
      #self.env_atomnum2countofnnlsymandcoord[i] = tuple(temp_local_env.values())
      self.env_atomid[i]=temp_local_env

class compute_coordination:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Compute each atom's type coordination
  * This class helps to compute the coordination of each atom type
  * and store the data into an object
  *
  * usage: environment1=ganisetti_tools.compute_coordination(config,config_nnl,atom_type_num2sym)
  * where config                = ganisetti_tools.get_atoms_info_from_lammps(01.dump)
  *       config_nnl            = ganisetti_tools.compute_nnl(config,rc,atom_type_num2sym)
  *       atom_type_num2sym     = {"1":"O", "2":"Si", "3":"Al", etc,.}
  *
  * output: config1.avg_coord_sym  # config1.avg_coord_sym['Si'] = 4
  *         config1.avg_coord_num  # config1.avg_coord_num[atom_type_sym2num['Si']] = 4
  *         config1.individual_coord_sym # config1.individual_coord_sym[('Si',4)] = 100 atoms
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,config,config_nnl,atom_type_num2sym):
    self.avg_coord_num          = {}
    self.avg_coord_sym          = {}
    temp1_count                 = {}
    temp2_count                 = {}
    atoms_count                 = {}
    total_coord                 = {}
    nnl_individual_coord_count  = {}
    for i in atom_type_num2sym.keys():
      temp1_count={i:0}
      atoms_count.update(temp1_count)
      total_coord.update(temp1_count)
      for j in range(config_nnl.max_nnl_each_atom_type_sym[atom_type_num2sym[i]]+1):
        temp2_count={(atom_type_num2sym[i],j):0}
        nnl_individual_coord_count.update(temp2_count)

    for i in config.id:
      temp1_count={config.type[i]:atoms_count[config.type[i]]+1}
      atoms_count.update(temp1_count)
      temp2_count={config.type[i]:total_coord[config.type[i]]+config_nnl.nnl_count[i]}
      total_coord.update(temp2_count)
      temp3=atom_type_num2sym[config.type[i]]
      temp4_count={(temp3,config_nnl.nnl_count[i]):nnl_individual_coord_count[(temp3,config_nnl.nnl_count[i])]+1}
      nnl_individual_coord_count.update(temp4_count)
    for i in atom_type_num2sym.keys():
      self.avg_coord_num[i]=float(total_coord[i]/atoms_count[i])
      self.avg_coord_sym[atom_type_num2sym[i]]=float(total_coord[i]/atoms_count[i])
    self.individual_coord_sym=nnl_individual_coord_count


class read_command_line_deprecated:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Reading command line and arguments
  * This class helps to read the command line arguments
  * and store the data into relavent parameters
  * usage: cmd=ganisetti_tools.read_command_line(sys.argv)
  *
  * output: cmd.atom_type            = stores given atom types                           # atom_type['Si'] = 1
  *         cmd.atom_type_sym2num    = stores given atom types based on chemical symbol  # atom_type_sym2num['Si'] = 1
  *         cmd.atom_type_num2sym    = stores given atom types based on atom type		 # atom_type_sym2num['1'] = Si
  *         cmd.bond_length          = stores given bond lengths	                     # bond_length['Si'] = 2.2
  *         cmd.bond_length_sys2num  = stores given bond lengths based on chemical symbol# bond_length_sym2num['Si']=2.2
  *         cmd.bond_length_num2num  = stores given bond lengths based on atom type	     # bond_length_num2num['1'] =2.2
  *         cmd.rc                   = stores the given cutoff radii#rc[atom_type_sym2num['Si']][atom_type_sym2num['O']]
  *         cmd.given_all_atom_types = stores the given atom types # given_all_atom_types = 1,2,3,4, etc.
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,sys_argv):
    error1 = 0
    error2 = 0
    debug = "yes"
    tot_argv = len(sys_argv)
    CREDBG = '\33[41m'
    CREDBGEND = '\x1b[0m'
    while error1 == 0:
      # checking for minimum arguments
      if tot_argv < 8:
        error1 = 1
        error1_statement = "Error info: a minimum of 7 arguments should be passed to the program"
        break
      if (tot_argv % 2) != 0:
        error1 = 1
        error1_statement = "Error info: odd number of arguments are expecting to be passed to the program"
        break
      commandline = []
      recognized_argvs_idex = []
      for i in sys_argv:
        commandline.append(i)
      commandline = list(commandline)
      atom_type           = {}  # store all given atom types				            # atom_type['Si'] = 1
      bond_length         = {}  # store all given bond lengths				            # bond_length['Si'] = 2.2
      atom_type_sym2num   = {}  # store all given atom types based on chemical symbol 	# atom_type_sym2num['Si'] = 1
      atom_type_num2sym   = {}  # store all given atom types based on atom type		    # atom_type_sym2num['1'] = Si
      bond_length_sym2num = {}  # store all given bond lengths based on chemical symbol	# bond_length_sym2num['Si']=2.2
      bond_length_num2num = {}  # store all given bond lengths based on atom type	    # bond_length_num2num['1'] = 2.2

      given_all_atom_types = []
      # search for atom types
      if "-O" in commandline:
        O_index                           = commandline.index("-O")
        given_O_type                      = int(commandline[O_index + 1])
        atom_type['O']                    = given_O_type
        atom_type_sym2num['O']            = given_O_type
        atom_type_num2sym[given_O_type]   = 'O'
        given_all_atom_types.append(given_O_type)
        recognized_argvs_idex.append(O_index)
        recognized_argvs_idex.append(O_index + 1)
      if "-Si" in commandline:
        Si_index                          = commandline.index("-Si")
        given_Si_type                     = int(commandline[Si_index + 1])
        atom_type['Si']                   = given_Si_type
        atom_type_sym2num['Si']           = given_Si_type
        atom_type_num2sym[given_Si_type]  = 'Si'
        given_all_atom_types.append(given_Si_type)
        recognized_argvs_idex.append(Si_index)
        recognized_argvs_idex.append(Si_index + 1)
      if "-Al" in commandline:
        Al_index                          = commandline.index("-Al")
        given_Al_type                     = int(commandline[Al_index + 1])
        atom_type['Al']                   = given_Al_type
        atom_type_sym2num['Al']           = given_Al_type
        atom_type_num2sym[given_Al_type]  = 'Al'
        given_all_atom_types.append(given_Al_type)
        recognized_argvs_idex.append(Al_index)
        recognized_argvs_idex.append(Al_index + 1)
      if "-P" in commandline:
        P_index                           = commandline.index("-P")
        given_P_type                      = int(commandline[P_index + 1])
        atom_type['P']                    = given_P_type
        atom_type_sym2num['P']            = given_P_type
        atom_type_num2sym[given_P_type]   = 'P'
        given_all_atom_types.append(given_P_type)
        recognized_argvs_idex.append(P_index)
        recognized_argvs_idex.append(P_index + 1)
      if "-Ca" in commandline:
        Ca_index                          = commandline.index("-Ca")
        given_Ca_type                     = int(commandline[Ca_index + 1])
        atom_type['Ca']                   = given_Ca_type
        atom_type_sym2num['Ca']           = given_Ca_type
        atom_type_num2sym[given_Ca_type]  = 'Ca'
        given_all_atom_types.append(given_Ca_type)
        recognized_argvs_idex.append(Ca_index)
        recognized_argvs_idex.append(Ca_index + 1)
      if "-Mg" in commandline:
        Mg_index                          = commandline.index("-Mg")
        given_Mg_type                     = int(commandline[Mg_index + 1])
        atom_type['Mg']                   = given_Mg_type
        atom_type_sym2num['Mg']           = given_Mg_type
        atom_type_num2sym[given_Mg_type]  = 'Mg'
        given_all_atom_types.append(given_Mg_type)
        recognized_argvs_idex.append(Mg_index)
        recognized_argvs_idex.append(Mg_index + 1)
      if "-Sr" in commandline:
        Sr_index                          = commandline.index("-Sr")
        given_Sr_type                     = int(commandline[Sr_index + 1])
        atom_type['Sr']                   = given_Sr_type
        atom_type_sym2num['Sr']           = given_Sr_type
        atom_type_num2sym[given_Sr_type]  = 'Sr'
        given_all_atom_types.append(given_Sr_type)
        recognized_argvs_idex.append(Sr_index)
        recognized_argvs_idex.append(Sr_index + 1)
      if "-Na" in commandline:
        Na_index                          = commandline.index("-Na")
        given_Na_type                     = int(commandline[Na_index + 1])
        atom_type['Na']                   = given_Na_type
        atom_type_sym2num['Na']           = given_Na_type
        atom_type_num2sym[given_Na_type]  = 'Na'
        given_all_atom_types.append(given_Na_type)
        recognized_argvs_idex.append(Na_index)
        recognized_argvs_idex.append(Na_index + 1)
      if "-F" in commandline:
        F_index                           = commandline.index("-F")
        given_F_type                      = int(commandline[F_index + 1])
        atom_type['F']                    = given_F_type
        atom_type_sym2num['F']            = given_F_type
        atom_type_num2sym[given_Na_type]  = 'F'
        given_all_atom_types.append(given_F_type)
        recognized_argvs_idex.append(F_index)
        recognized_argvs_idex.append(F_index + 1)

      # O atom type must be given
      if "O" not in atom_type_sym2num.keys():
        error1 = 1
        error1_statement = "Error info: Oxygen atom type is not given"
        break
      # initilize a cutoff radii array
      rc = [[0.0 for i in range(max(given_all_atom_types) + 1)] for j in range(max(given_all_atom_types) + 1)]
      # search for bond lengths
      if "-SiO" in commandline:
        if "Si" not in atom_type_sym2num:
          error1 = 1
          error1_statement = "Error info: Si atom type is not given but SiO bond length is given"
          break
        tmp_index                                           = commandline.index("-SiO")
        bond_length['SiO']                                  = float(commandline[tmp_index + 1])
        bond_length_sym2num['Si']                           = float(commandline[tmp_index + 1])
        bond_length_num2num[atom_type_sym2num['Si']]        = float(commandline[tmp_index + 1])
        rc[atom_type_sym2num['Si']][atom_type_sym2num['O']] = float(commandline[tmp_index + 1])
        rc[atom_type_sym2num['O']][atom_type_sym2num['Si']] = float(commandline[tmp_index + 1])
        recognized_argvs_idex.append(tmp_index)
        recognized_argvs_idex.append(tmp_index + 1)
      if "-AlO" in commandline:
        if "Al" not in atom_type_sym2num:
          error1 = 1
          error1_statement = "Error info: Al atom type is not given but AlO bond length is given"
          break
        tmp_index                                           = commandline.index("-AlO")
        bond_length['AlO']                                  = float(commandline[tmp_index + 1])
        bond_length_sym2num['Al']                           = float(commandline[tmp_index + 1])
        bond_length_num2num[atom_type_sym2num['Al']]        = float(commandline[tmp_index + 1])
        rc[atom_type_sym2num['Al']][atom_type_sym2num['O']] = float(commandline[tmp_index + 1])
        rc[atom_type_sym2num['O']][atom_type_sym2num['Al']] = float(commandline[tmp_index + 1])
        recognized_argvs_idex.append(tmp_index)
        recognized_argvs_idex.append(tmp_index + 1)
      if "-PO" in commandline:
        if "P" not in atom_type_sym2num:
          error1 = 1
          error1_statement = "Error info: P atom type is not given but PO bond length is given"
          break
        tmp_index                                           = commandline.index("-PO")
        bond_length['PO']                                   = float(commandline[tmp_index + 1])
        bond_length_sym2num['P']                            = float(commandline[tmp_index + 1])
        bond_length_num2num[atom_type_sym2num['P']]         = float(commandline[tmp_index + 1])
        rc[atom_type_sym2num['P']][atom_type_sym2num['O']]  = float(commandline[tmp_index + 1])
        rc[atom_type_sym2num['O']][atom_type_sym2num['P']]  = float(commandline[tmp_index + 1])
        recognized_argvs_idex.append(tmp_index)
        recognized_argvs_idex.append(tmp_index + 1)
      if "-CaO" in commandline:
        if "Ca" not in atom_type_sym2num:
          error1 = 1
          error1_statement = "Error info: Ca atom type is not given but CaO bond length is given"
          break
        tmp_index                                           = commandline.index("-CaO")
        bond_length['CaO']                                  = float(commandline[tmp_index + 1])
        bond_length_sym2num['Ca']                           = float(commandline[tmp_index + 1])
        bond_length_num2num[atom_type_sym2num['Ca']]        = float(commandline[tmp_index + 1])
        rc[atom_type_sym2num['Ca']][atom_type_sym2num['O']] = float(commandline[tmp_index + 1])
        rc[atom_type_sym2num['O']][atom_type_sym2num['Ca']] = float(commandline[tmp_index + 1])
        recognized_argvs_idex.append(tmp_index)
        recognized_argvs_idex.append(tmp_index + 1)
      if "-MgO" in commandline:
        if "Mg" not in atom_type_sym2num:
          error1 = 1
          error1_statement = "Error info: Mg atom type is not given but MgO bond length is given"
          break
        tmp_index                                           = commandline.index("-MgO")
        bond_length['MgO']                                  = float(commandline[tmp_index + 1])
        bond_length_sym2num['Mg']                           = float(commandline[tmp_index + 1])
        bond_length_num2num[atom_type_sym2num['Mg']]        = float(commandline[tmp_index + 1])
        rc[atom_type_sym2num['Mg']][atom_type_sym2num['O']] = float(commandline[tmp_index + 1])
        rc[atom_type_sym2num['O']][atom_type_sym2num['Mg']] = float(commandline[tmp_index + 1])
        recognized_argvs_idex.append(tmp_index)
        recognized_argvs_idex.append(tmp_index + 1)
      if "-NaO" in commandline:
        if "Na" not in atom_type_sym2num:
          error1 = 1
          error1_statement = "Error info: Na atom type is not given but NaO bond length is given"
          break
        tmp_index                                           = commandline.index("-NaO")
        bond_length['NaO']                                  = float(commandline[tmp_index + 1])
        bond_length_sym2num['Na']                           = float(commandline[tmp_index + 1])
        bond_length_num2num[atom_type_sym2num['Na']]        = float(commandline[tmp_index + 1])
        rc[atom_type_sym2num['Na']][atom_type_sym2num['O']] = float(commandline[tmp_index + 1])
        rc[atom_type_sym2num['O']][atom_type_sym2num['Na']] = float(commandline[tmp_index + 1])
        recognized_argvs_idex.append(tmp_index)
        recognized_argvs_idex.append(tmp_index + 1)
      if "-SrO" in commandline:
        if "Sr" not in atom_type_sym2num:
          error1 = 1
          error1_statement = "Error info: Sr atom type is not given but SrO bond length is given"
          break
        tmp_index                                           = commandline.index("-SrO")
        bond_length['SrO']                                  = float(commandline[tmp_index + 1])
        bond_length_sym2num['Sr']                           = float(commandline[tmp_index + 1])
        bond_length_num2num[atom_type_sym2num['Sr']]        = float(commandline[tmp_index + 1])
        rc[atom_type_sym2num['Sr']][atom_type_sym2num['O']] = float(commandline[tmp_index + 1])
        rc[atom_type_sym2num['O']][atom_type_sym2num['Sr']] = float(commandline[tmp_index + 1])
        recognized_argvs_idex.append(tmp_index)
        recognized_argvs_idex.append(tmp_index + 1)

      if "-CaF" in commandline:
        if "Ca" not in atom_type_sym2num:
          error1 = 1
          error1_statement = "Error info: Ca atom type is not given but CaF bond length is given"
          break
        if "F" not in atom_type_sym2num:
          error1 = 1
          error1_statement = "Error info: F atom type is not given but CaF bond length is given"
          break
        tmp_index                                           = commandline.index("-CaF")
        bond_length['CaF']                                  = float(commandline[tmp_index + 1])
        bond_length_sym2num['Ca']                           = float(commandline[tmp_index + 1])
        bond_length_num2num[atom_type_sym2num['Ca']]        = float(commandline[tmp_index + 1])
        rc[atom_type_sym2num['Ca']][atom_type_sym2num['F']] = float(commandline[tmp_index + 1])
        rc[atom_type_sym2num['F']][atom_type_sym2num['Ca']] = float(commandline[tmp_index + 1])
        recognized_argvs_idex.append(tmp_index)
        recognized_argvs_idex.append(tmp_index + 1)

      # checking for bond lengths of all given atom types
      for i in atom_type_sym2num:
        if i != "O":
          if i not in bond_length_sym2num:
            error1 = 1
            error1_statement = "Error info: " + str(i) + " atom type is given but " + str(
              i) + "O bond length is not given"
            break
      # checking for any unrecognized arguments passed into command line
      for i in range(2, len(commandline)):
        if i not in recognized_argvs_idex:
          error1 = 1
          error1_statement = "Error info: " + str(
            commandline[i]) + " is an unrecognozed argument passed into the command line"
          break

      break  # terminate the main while loop

    if debug == "no":
      for i in atom_type_sym2num:
        print(i, atom_type[i], atom_type_sym2num[i])
        if i != "O":
          print(bond_length_sym2num[i], bond_length_num2num[atom_type_sym2num[i]])
      print("\n\n")
      for i in bond_length_sym2num:
        print(i, bond_length[i + str("O")], bond_length_sym2num[i])

    if error1 != 0:
      ganisetti_tools.banner()
      print("************************************** S. Ganisetti **************************************")
      print("Error: usage is wrong")
      print("./this_program  ChkptFile -O 1 -Si 2 -Al 3 -Ca 4 -Mg 5 -SiO 2.0 -AlO 2.0 -CaO 3.0 -MgO 3.0")
      print("The program is prepared for chemical compositions: SiO2, Al2O3, P2O5, Na2O, CaO and MgO")
      print("Please specify both atom types and the cutoff radii of all required pair of atoms")
      print("%s %s %s" % (CREDBG, str(error1_statement), CREDBGEND))
      print("******************************************************************************************")
      sys.exit(0)

    self.atom_type            = atom_type
    self.atom_type_sym2num    = atom_type_sym2num
    self.atom_type_num2sym    = atom_type_num2sym
    self.bond_length          = bond_length
    self.bond_length_sys2num  = bond_length_sym2num
    self.bond_length_num2num  = bond_length_num2num
    self.rc                   = rc
    self.given_all_atom_types = given_all_atom_types


class read_command_line:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Reading command line and arguments
  * This class helps to read the command line arguments
  * and stores the data into relavent parameters
  * usage: cmd=ganisetti_tools.read_command_line(sys.argv)
  *
  * output: cmd.atom_type_sym2num    = stores given atom types based on chemical symbol  # atom_type_sym2num['Si'] = 1
  *         cmd.atom_type_num2sym    = stores given atom types based on atom type		 # atom_type_sym2num['1'] = Si
  *         cmd.bond_length_sys2num  = stores given bond lengths based on chemical symbol# bond_length_sym2num['Si']=2.2
  *         cmd.bond_length_num2num  = stores given bond lengths based on atom type	     # bond_length_num2num['1'] =2.2
  *         cmd.rc                   = stores the given cutoff radii#rc[atom_type_sym2num['Si']][atom_type_sym2num['O']]
  *         cmd.error                = if "none" is returned then the command line is successfully passed
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,sys_argv):
    error1_statement = "none"
    self.error = "none"
    tot_argv = len(sys_argv)
    CREDBG = '\33[31m' # \33[41m for red background
    CREDBGEND = '\x1b[0m'
    known_cations = ['Si','Al','P','Na','Ca','Mg','Sr'] # The cations list should be updated if you want to add a new cation
    known_anions  = ['O','F']                      # The anions list should be updated if you want to add a new anion
    known_all_atom_types = known_cations + known_anions
    known_network_formers=['Si','Al','P']
    known_modifiers      =['Na','Ca','Mg','Sr']

    # while is using to handle the error rather iterating the loop
    while error1_statement == "none":
      # checking for minimum arguments
      if tot_argv < 8:
        error1_statement = "Error info: a minimum of 7 arguments should be passed to the program"
        break
      if (tot_argv % 2) != 0:
        error1_statement = "Error info: odd number of arguments are expecting to be passed to the program"
        break
      commandline = []
      recognized_argvs_index = []
      for i in sys_argv:
        commandline.append(i)
      commandline         = list(commandline)

      # define the parameters
      # atom_type           = {} # store all given atom types				             # atom_type['Si'] = 1
      # bond_length         = {} # store all given bond lengths				             # bond_length['Si'] = 2.2
      atom_type_sym2num     = {} # store all given atom types based on chemical symbol 	 # atom_type_sym2num['Si'] = 1
      atom_type_num2sym     = {} # store all given atom types based on atom type		 # atom_type_sym2num['1'] = Si
      bond_length_sym2num   = {} # store all given bond lengths based on chemical symbol # bond_length_sym2num['Si']=2.2
      bond_length_num2num   = {} # store all given bond lengths based on atom type	     # bond_length_num2num['1']=2.2
      given_anions_sym2num  = {} # store all given anions
      given_cations_sym2num = {} # store all given cations
      #given_all_atom_types = [] # store the given atom types as a list                # given_all_atom_type =1,2, etc.
      given_formers_sym2num = {} # store all given network formers
      given_modifiers_sym2num = {} # store all given network modifiers

      # search for atom type
      for known_atom_type in known_all_atom_types:
        known_atom_type_string = str("-")+known_atom_type
        if known_atom_type_string in commandline:
          given_atom_type_sym                     =  known_atom_type
          known_atom_type_index                   = int(commandline.index(known_atom_type_string))
          given_atom_type_num                     = int(commandline[known_atom_type_index + 1])
          atom_type_sym2num[given_atom_type_sym]  = given_atom_type_num
          atom_type_num2sym[given_atom_type_num]  = given_atom_type_sym
          recognized_argvs_index.append(known_atom_type_index)
          recognized_argvs_index.append(known_atom_type_index+1)
          if given_atom_type_sym in known_cations:
            given_cations_sym2num[given_atom_type_sym]=atom_type_sym2num[given_atom_type_sym]
          if given_atom_type_sym in known_anions:
            given_anions_sym2num[given_atom_type_sym]=atom_type_sym2num[given_atom_type_sym]
          if given_atom_type_sym in known_network_formers:
            given_formers_sym2num[given_atom_type_sym]=atom_type_sym2num[given_atom_type_sym]
          if given_atom_type_sym in known_modifiers:
            given_modifiers_sym2num[given_atom_type_sym]=atom_type_sym2num[given_atom_type_sym]

      # at least one anion should be given
      if len(given_anions_sym2num.keys()) < 1 or len(given_cations_sym2num.keys()) < 1:
        error1_statement="Error info: at least one anion and one cation should be given"
        break

      # initilize a 2d array for cutoff radii
      rc = [[0.0 for i in range(max(atom_type_num2sym.keys()) + 1)] for j in range(max(atom_type_num2sym.keys()) + 1)]

      # search for atom pairs
      known_atom_pairs=[(x,y) for x in known_cations for y in known_anions]
      for known_atom_pair in known_atom_pairs:
        kap1=known_atom_pair[0]
        kap2=known_atom_pair[1]
        known_atom_pair_string=str("-")+kap1+kap2
        if known_atom_pair_string in commandline:
          known_atom_pair_index = int(commandline.index(known_atom_pair_string))
          bond_length_sym2num[(kap1,kap2)] = float(commandline[known_atom_pair_index + 1])
          bond_length_num2num[(atom_type_sym2num[kap1],atom_type_sym2num[kap2])]=float(commandline[known_atom_pair_index+1])
          rc[atom_type_sym2num[kap1]][atom_type_sym2num[kap2]] = float(commandline[known_atom_pair_index + 1])
          rc[atom_type_sym2num[kap2]][atom_type_sym2num[kap1]] = float(commandline[known_atom_pair_index + 1])
          recognized_argvs_index.append(known_atom_pair_index)
          recognized_argvs_index.append(known_atom_pair_index+1)
          # checking if bond lenth of an atom is given but one of the atom types is not given
          if not kap1 in atom_type_sym2num.keys():
            error1_statement="Error info: "+str(kap1)+" atom type is not given but "+str(known_atom_pair)+" bond length is given"
            break
          if not kap2 in atom_type_sym2num.keys():
            error1_statement="Error info: "+str(kap2)+" atom type is not given but "+str(known_atom_pair)+" bond length is given"
            break

      # checking if one of the atom types is given but bond length is not given
      for i in given_cations_sym2num.keys():
        for j in given_anions_sym2num.keys():
          if str("-")+i+j not in commandline:
            error1_statement=str(i)+" and "+str(j)+" atoms are given but "+str(i)+str(j)+" bondlength is not given"
            break

      # checking for any unrecognized arguments passed into the command line
      for i in range(2, len(commandline)):
        if i not in recognized_argvs_index:
          error1_statement="Error info: an unrecognized argument "+str(commandline[i])+" is passed into the command line"
          break
      break  # breaking the main while loop

    if error1_statement != "none":
      sudheer_banner()
      print("************************************** S. Ganisetti **************************************")
      print("Error: usage is wrong")
      print("correct usage is: ./this_program  ChkptFile -O 1 -Si 2 -Al 3 -Ca 4 -Mg 5 -SiO 2.0 -AlO 2.0 -CaO 3.0 -MgO 3.0")
      #print("The program is prepared for chemical compositions: SiO2, Al2O3, P2O5, Na2O, CaO and MgO")
      print("This program is prepared for ",end=" ")
      for i in range(len(known_all_atom_types)-2):
        print(known_all_atom_types[i],end=", ")
      print("%s and %s atom types" %(known_all_atom_types[i+1],known_all_atom_types[i+2]))
      print("Please specify all atom types and the cutoff radii of all required atom pairs\n")
      print("additional error information: ")
      print("%s %s %s" % (CREDBG, str(error1_statement), CREDBGEND))
      print("******************************************************************************************")
      self.error = "error"

    if error1_statement == "none":
      self.atom_type_sym2num    = atom_type_sym2num
      self.atom_type_num2sym    = atom_type_num2sym
      self.bond_length_sys2num  = bond_length_sym2num
      self.bond_length_num2num  = bond_length_num2num
      self.rc                   = rc
      self.given_anions_sym2num = given_anions_sym2num
      self.given_cations_sym2num= given_cations_sym2num
      self.given_formers_sym2num= given_formers_sym2num
      self.given_modifiers_sym2num = given_modifiers_sym2num
    # self.atom_type            = atom_type
    # self.bond_length          = bond_length
    # self.given_all_atom_types = given_all_atom_types


class compute_pair_distribution:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Computing pair distribution
  * This class helps to compute the pair distribution
  * and stores the data into relavent parameters
  * usage: config1_pdf=ganisetti_tools.compute_pair_distribution(config,cmd,rcut,total_bins)
  *
  * input:  config1	= ganisetti_tools.get_atoms_info_from_lammps(LAMMPS_DUMP_FILE)
  *         cmd     = ganisetti_tools.read_command_line(sys.argv)
  *
  * output: config1_pdf.gr[('Si','O',bin_number)]=distribution
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,config,cmd,rcut,total_bins):

    atom_type_sym2num=cmd.atom_type_sym2num
    atom_type_num2sym=cmd.atom_type_num2sym
    given_anions=cmd.given_anions_sym2num
    given_cations=cmd.given_cations_sym2num
    bin_size=float(1.0*rcut/total_bins)

    self.gr={}
    gr={}
    temp1=[(i,j,k) for i in given_cations for j in given_anions for k in range(int(total_bins))]
    for i in temp1:
      temp2={i:0}
      gr.update(temp2)

    # Devide the box into voxels
    MaxCells_X = int(config.box_xx / rcut)
    MaxCells_Y = int(config.box_yy / rcut)
    MaxCells_Z = int(config.box_zz / rcut)

    # LEC= Length of Each Cell
    LEC_X = float(config.box_xx / MaxCells_X)
    LEC_Y = float(config.box_yy / MaxCells_Y)
    LEC_Z = float(config.box_zz / MaxCells_Z)

    # cell_atoms_count[X][Y][Z]
    cell_atoms_count = [[[0 for i in range(MaxCells_Z)] for j in range(MaxCells_Y)] for k in range(MaxCells_X)]
    # cell_atoms_list[X][Y][Z][200]
    cell_atoms_list = [[[[0 for i in range(200)] for j in range(MaxCells_Z)] for k in range(MaxCells_Y)] for l in
                       range(MaxCells_X)]

    # store atom positions locally
    atom_posx = config.posx
    atom_posy = config.posy
    atom_posz = config.posz
    atom_type = config.type

    # move the box to center and correspondingly atoms if the corner of the box is not at center
    for i in config.id:
      atom_posx[i] = config.posx[i] - config.box[0][0]
      atom_posy[i] = config.posy[i] - config.box[1][0]
      atom_posz[i] = config.posz[i] - config.box[2][0]
      # If the atoms are outside of the box then bring them into other side according to PBC
      if atom_posx[i] >= config.box_xx:
        atom_posx[i] = atom_posx[i] - config.box_xx
      if atom_posy[i] >= config.box_yy:
        atom_posy[i] = atom_posy[i] - config.box_yy
      if atom_posz[i] >= config.box_zz:
        atom_posz[i] = atom_posz[i] - config.box_zz
      # assign the atoms into respective cell number  #assign the atoms into respective cell number
      cellX = int(atom_posx[i] / LEC_X)
      cellY = int(atom_posy[i] / LEC_Y)
      cellZ = int(atom_posz[i] / LEC_Z)
      cell_atoms_list[cellX][cellY][cellZ][cell_atoms_count[cellX][cellY][cellZ]] = i
      cell_atoms_count[cellX][cellY][cellZ] = cell_atoms_count[cellX][cellY][cellZ] + 1

    # Main loop to compute pair distribution and other information
    # (i,j,k) => selecting atom which is belongs to one cell
    for i in range(MaxCells_X):
      for j in range(MaxCells_Y):
        for k in range(MaxCells_Z):
          for atom1 in range(cell_atoms_count[i][j][k]):
            atom1_id = cell_atoms_list[i][j][k][atom1]
            pos_x1 = atom_posx[atom1_id]
            pos_y1 = atom_posy[atom1_id]
            pos_z1 = atom_posz[atom1_id]

            # (l,m,n) => selecting neighbour cell
            for l in np.arange(0, 2, 1):
              new_i = i + l
              an = new_i
              pbc_x = 0.0
              # -----------------------------------------------------
              # Dealing with periodic boundary condition along x-axis
              if new_i < 0:
                # But the conditions will not lead to excute this case
                new_i = MaxCells_X - 1
                print("WARNING!!! There is somethong wrong with the pbc in nnl, please check it well !!!\n")
              if new_i > MaxCells_X - 1:
                new_i = 0
                pbc_x = config.box_xx
              # -----------------------------------------------------
              for m in np.arange(-1 * l, 2, 1):
                new_j = j + m
                bn = new_j
                pbc_y = 0.0
                # -----------------------------------------------------
                # Dealing with periodic boundary condition along y-axis
                if new_j < 0:
                  new_j = MaxCells_Y - 1
                  pbc_y = -1.0 * config.box_yy
                if new_j > MaxCells_Y - 1:
                  new_j = 0
                  pbc_y = config.box_yy
                # -----------------------------------------------------
                n1 = -1
                if l == m:
                  n1 = -1 * l
                for n in np.arange(n1, 2, 1):
                  new_k = k + n
                  cn = new_k
                  pbc_z = 0.0
                  # -----------------------------------------------------
                  # Dealing with periodic boundary conditions along z-axis
                  if new_k < 0:
                    new_k = MaxCells_Z - 1
                    pbc_z = -1.0 * config.box_zz
                  if new_k > MaxCells_Z - 1:
                    new_k = 0
                    pbc_z = config.box_zz
                  # debugging
                  # if i==0 and j==0 and k==0:
                  #  print ("%d %d %d %d %d %d   atom1_id %d") %(l,m,n,an,bn,cn,atom1_id)

                  # -----------------------------------------------------
                  # This is a trick to avoid double counting the interactions if the two atoms are in the same cell
                  temp = 0
                  if i == new_i and j == new_j and k == new_k:
                    temp = atom1 + 1
                  # selecting the neighbour atom
                  for atom2 in range(temp, cell_atoms_count[new_i][new_j][new_k]):
                    atom2_id = cell_atoms_list[new_i][new_j][new_k][atom2]
                    temp_atm1=atom_type_num2sym[atom_type[atom1_id]]
                    temp_atm2=atom_type_num2sym[atom_type[atom2_id]]
                    if ( temp_atm1 in given_anions and temp_atm2 in given_cations) or (temp_atm1 in given_cations and temp_atm2 in given_anions):
                      pos_x2 = atom_posx[atom2_id] + pbc_x
                      pos_y2 = atom_posy[atom2_id] + pbc_y
                      pos_z2 = atom_posz[atom2_id] + pbc_z

                      r = np.linalg.norm(np.array([pos_x1, pos_y1, pos_z1]) - np.array([pos_x2, pos_y2, pos_z2]))
                      if r < rcut and r > 0.2: # rcut value is excluded
                        temp_bin=int(1.0*r/bin_size)
                        if temp_atm1 in given_cations:
                          temp3={(temp_atm1,temp_atm2,temp_bin):gr[(temp_atm1,temp_atm2,temp_bin)]+1}
                        else :
                          temp3 = {(temp_atm2, temp_atm1, temp_bin): gr[(temp_atm2, temp_atm1, temp_bin)] + 1}
                        gr.update(temp3)
    self.gr=gr

class compute_Q_deprecated:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Computing Qn speciation
  *
  * usage: config1_QSi=ganisetti_tools.compute_Q(cmd,config,config_nnl,config_env,'Si')
  *
  * input:  config1	= ganisetti_tools.get_atoms_info_from_lammps(LAMMPS_DUMP_FILE)
  *         cmd     = ganisetti_tools.read_command_line(sys.argv)
  *
  * output: config1_QSi.Q = list(Q[0],Q[1],Q[2],Q[3] and Q[4])
  *         config1_QSi.Q_status_atomid2num = {0:-1,1;-1,2:1,3:-1} 0,1,3 are other atom tpyes, 2 is requred type
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,cmd,config,config_nnl,config_env,T):
    self.Q_status_atomid2num={}
    atom_type_sym2num = cmd.atom_type_sym2num
    atom_type_num2sym = cmd.atom_type_num2sym
    given_anions_sym2num = cmd.given_anions_sym2num
    given_cations_sym2num = cmd.given_cations_sym2num
    given_formers_sym2num = cmd.given_formers_sym2num
    given_modifiers_sym2num = cmd.given_modifiers_sym2num
    Q=[0 for i in range(5)]

    Q_status_atomid2num={}
    Q_neighboring_formers_list_atomid2list={}
    for i in config.id:
      temp1={i:-1}
      Q_status_atomid2num.update(temp1)

    for i in config.id:
      if atom_type_num2sym[config.type[i]] == T and config_nnl.nnl_count[i] == 4:
        BA_count=0
        NBA_count=0
        for j in config_nnl.nnl[i]: # j=O1, O2, O3, O4
          Q_Former4Coord       = 0
          Q_FormerButNot4Coord = 0
          Q_NonFormerAnyCoord  = 0

          for k in config_env.env_atomid[j].keys(): # k=[(Si,4),(Al,4),(Ca,5)]
            if k[0] in given_formers_sym2num.keys():
              if k[1] == 4:
                Q_Former4Coord = Q_Former4Coord + config_env.env_atomid[j][k]
              else :
                Q_FormerButNot4Coord = Q_FormerButNot4Coord + config_env.env_atomid[j][k] # counts (Si,3),(Al,5)
                # this is nothing but Si4-O-Al5
            else:
              Q_NonFormerAnyCoord = Q_NonFormerAnyCoord + 1
          if Q_Former4Coord == 2:
            BA_count = BA_count + 1
          else:
            NBA_count = NBA_count + 1
        Q[BA_count] = Q[BA_count] + 1
        temp1={i:BA_count}
        Q_status_atomid2num.update(temp1)
        '''
          for l in given_formers_sym2num.keys(): # l=Si,Al,P
            for m in range(config_nnl.max_nnl_each_atom_type_sym[l]+1): # m = 0 to 5
              if m == 4:
                # Q main
                if (l,m) in config_env.env_atomid[j].keys(): # counts only (Si,4) or (Al,4) or (P,4):1
                  Q_Former4Coord=Q_Former4Coord+config_env.env_atomid[j][(l,m)]
              else:
                if 
                Q_FormerButNot4Coord=Q_FormerButNot4Coord+config_env.env_atomid[j][(l,m)] # counts (Si,3), (Al,5), etc.
          Q_NonFormerAnyCoord=len(config_env.env_atomid[j].keys())-Q_FormerButNot4Coord-Q_Former4Coord # counts (Ca,5)
          if Q_Former4Coord == 2:
            BA_count=BA_count+1
          else:
            NBA_count=NBA_count+1
        Q[BA_count]=Q[BA_count]+1
        '''
    self.Q=Q
    self.Q_status_atomid2num

class compute_Q:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Computing Qn speciation
  *
  * usage: config1_QSi=ganisetti_tools.compute_Q(cmd,config,config_nnl,config_env,'Si')
  *
  * input:  config1	= ganisetti_tools.get_atoms_info_from_lammps(LAMMPS_DUMP_FILE)
  *         cmd     = ganisetti_tools.read_command_line(sys.argv)
  *
  * output: config1_QSi.Q_summary = {(Si,0):10,(Si,1):22, etc.}
  *         config1_QSi.Q_status_atomid2num = {0:-1,1;-1,2:1,3:-1} 0,1,3 are other atom tpyes, 2 is requred type
  *         config1_QSi.Q_4CoordFormers_list = {(Si,0):[1,2,3,etc.],(Si,1):[1,2,3,etc.]}
  *         config1_QSi.Q_non4CoordFormers_list = {(Si,1):[1,2,3, etc],(Al,1)}
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self, cmd, config, config_nnl, config_env):
    self.Q_status_atomid2num      = {}
    self.Q_summary                = {}
    self.Q_4CoordFormers_list     = {}
    self.Q_non4CoordFormers_list  = {}
    atom_type_sym2num = cmd.atom_type_sym2num
    atom_type_num2sym = cmd.atom_type_num2sym
    given_anions_sym2num = cmd.given_anions_sym2num
    given_cations_sym2num = cmd.given_cations_sym2num
    given_formers_sym2num = cmd.given_formers_sym2num
    given_modifiers_sym2num = cmd.given_modifiers_sym2num

    Q_summary               = {}
    Q_status_atomid2num     = {}
    Q_4CoordFormers_list    = {}
    Q_non4CoordFormers_list = {}
    for i in given_formers_sym2num.keys():
      for j in range(5):
        temp1 = {(i,j):0}
        Q_summary.update(temp1)
        temp1 = {(i,j):[]}
        Q_4CoordFormers_list.update(temp1)
        Q_non4CoordFormers_list.update(temp1)
    for i in config.id:
      temp1={i:-1}
      Q_status_atomid2num.update(temp1)

    for i in config.id:
      if config.type[i] in given_formers_sym2num.values() and config_nnl.nnl_count[i] == 4: # the main Q cation with 4 coordination
        BridgingAnions_count            = 0   # for atom i
        NonBridgingAnions_count         = 0   # for atom i
        Triplet_NonBridgingAnions_count = 0   # for atom i
        NetworkFormersList_4coord       = []  # for atom i
        NetworkFormersList_non4coord    = []  # for atom i
        for j in config_nnl.nnl[i]: # j = Q neighbouring anions: O1,O2,O3,O4
          Anion_4CoordinatedNetworkFormers_count    = 0 # for atom j
          Anion_non4CoordinatedNetworkFormers_count = 0 # for atom j
          for k in config_nnl.nnl[j]: # k =Q neighboring cations Si,Si (including the Q cation)
            if config.type[k] in given_formers_sym2num.values():
              if config_nnl.nnl_count[k] == 4:
                Anion_4CoordinatedNetworkFormers_count=Anion_4CoordinatedNetworkFormers_count+1
                NetworkFormersList_4coord.append(k)   # This is for Si[4]-O-Si[4], Si[4]-O-Al[4], etc
              else :
                Anion_non4CoordinatedNetworkFormers_count=Anion_non4CoordinatedNetworkFormers_count+1
                NetworkFormersList_non4coord.append(k) # This is for Si[4]-O-Si[4], Si[4]-O-Al[5], etc
          if Anion_4CoordinatedNetworkFormers_count == 2:
            BridgingAnions_count            = BridgingAnions_count+1
          elif Anion_4CoordinatedNetworkFormers_count == 1:
            NonBridgingAnions_count         = NonBridgingAnions_count+1
          elif Anion_4CoordinatedNetworkFormers_count == 3:
            Triplet_NonBridgingAnions_count = Triplet_NonBridgingAnions_count+1
        temp1={(atom_type_num2sym[config.type[i]],BridgingAnions_count):Q_summary[(atom_type_num2sym[config.type[i]],BridgingAnions_count)]+1}
        Q_summary.update(temp1)
        self.Q_summary=Q_summary

        temp1={i:BridgingAnions_count}
        Q_status_atomid2num.update(temp1)
        self.Q_status_atomid2num=Q_status_atomid2num

        NetworkFormersList_4coord=NetworkFormersList_4coord+Q_4CoordFormers_list[(atom_type_num2sym[config.type[i]],BridgingAnions_count)]
        NetworkFormersList_4coord=list(dict.fromkeys(NetworkFormersList_4coord))
        temp1 = {(atom_type_num2sym[config.type[i]],BridgingAnions_count):NetworkFormersList_4coord}
        Q_4CoordFormers_list.update(temp1)
        self.Q_4CoordFormers_list=Q_4CoordFormers_list

        NetworkFormersList_non4coord = NetworkFormersList_non4coord + Q_non4CoordFormers_list[(atom_type_num2sym[config.type[i]],BridgingAnions_count)]
        NetworkFormersList_non4coord = list(dict.fromkeys(NetworkFormersList_non4coord))
        temp1 = {(atom_type_num2sym[config.type[i]], BridgingAnions_count): NetworkFormersList_non4coord}
        Q_non4CoordFormers_list.update(temp1)
        self.Q_non4CoordFormers_list=Q_non4CoordFormers_list
'''
class compute_anions_distribution:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Computing anions speciation
  *
  * usage: config1_anions_dist=ganisetti_tools.compute_anions_distribution(cmd,config,config_nnl,config_env,'Si')
  *
  * input:  config1	= ganisetti_tools.get_atoms_info_from_lammps(LAMMPS_DUMP_FILE)
  *         cmd     = ganisetti_tools.read_command_line(sys.argv)
  *
  * output: config1_QSi.Q_summary = {(Si,0):10,(Si,1):22, etc.}
  *         config1_QSi.Q_status_atomid2num = {0:-1,1;-1,2:1,3:-1} 0,1,3 are other atom tpyes, 2 is requred type
  *         config1_QSi.Q_4CoordFormers_list = {(Si,0):[1,2,3,etc.],(Si,1):[1,2,3,etc.]}
  *         config1_QSi.Q_non4CoordFormers_list = {(Si,1):[1,2,3, etc],(Al,1)}
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self, cmd, config, config_nnl, config_env):
    self.Q_status_atomid2num      = {}
    self.Q_summary                = {}
    self.Q_4CoordFormers_list     = {}
    self.Q_non4CoordFormers_list  = {}
    atom_type_sym2num = cmd.atom_type_sym2num
    atom_type_num2sym = cmd.atom_type_num2sym
    given_anions_sym2num = cmd.given_anions_sym2num
    given_cations_sym2num = cmd.given_cations_sym2num
    given_formers_sym2num = cmd.given_formers_sym2num
    given_modifiers_sym2num = cmd.given_modifiers_sym2num

    bonds_count={}
    total_bonds_count={}
    total_atom_types={}
    total_atom_types.update(given_anions_sym2num)
    total_atom_types.update(given_cations_sym2num)
    for i in total_atom_types:
      temp1={i:0}
      total_bonds_count.update(temp1)

    for i in config.id:
      for j in config_nnl.nnl[i]:
        temp1={atom_type_num2sym[config.type[j]]:bonds_count[atom_type_num2sym[config.type[j]]]+1}
      if config.type[i] in given_anions_sym2num.values(): # i = anion
        for j in config_nnl.nnl[i]: # j = neighbouring cations
          temp1={atom_type_num2sym[config.type[j]]:bonds_count[atom_type_num2sym[config.type[j]]]+1}
          bonds_count.update(temp1)
'''
class compute_triplets:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Computing triplets
  *
  * usage: config1_triplets=ganisetti_tools.compute_triplets(cmd,config,config_nnl)
  *
  * input:  config1	= ganisetti_tools.get_atoms_info_from_lammps(LAMMPS_DUMP_FILE)
  *         cmd     = ganisetti_tools.read_command_line(sys.argv)
  *
  * output: config1_triplets.triplets_AmBCn2count = {(A,m,B,C,n):count, etc.} # A-B-C
  *         config1_triplets.total_triplets_sym2count= {O:100,Si:120,etc.} O, Si are cetral atoms
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self, cmd, config, config_nnl):
    self.triplets_AmBCn2count={}
    self.total_triplets_sym2count={}
    atom_type_sym2num = cmd.atom_type_sym2num
    atom_type_num2sym = cmd.atom_type_num2sym
    given_anions_sym2num = cmd.given_anions_sym2num
    given_cations_sym2num = cmd.given_cations_sym2num
    given_formers_sym2num = cmd.given_formers_sym2num
    given_modifiers_sym2num = cmd.given_modifiers_sym2num

    triplets_count  = {}
    total_triplets  = {}

    for i in atom_type_sym2num.keys():       # This is B in A-B-C
      for j in atom_type_sym2num.keys():     # This is A in A-B-C
        for k in atom_type_sym2num.keys():   # This is C in A-B-C
          for m in range(config_nnl.max_nnl_each_atom_type_sym[j] + 1):
            for n in range(config_nnl.max_nnl_each_atom_type_sym[k] + 1):
              temp1={(j,m,i,k,n):0}
              triplets_count.update(temp1)
      temp2={i:0}
      total_triplets.update(temp2)

    for i in config.id:
      if config_nnl.nnl_count[i] > 1:
        B = atom_type_num2sym[config.type[i]] # B is the middle atom in A-B-C
        for j in list( itertools.combinations(config_nnl.nnl[i], 2) ): # j=all combinations of (A and C) in A-B-C
          A = atom_type_num2sym[config.type[j[0]]]
          C = atom_type_num2sym[config.type[j[1]]]
          m=config_nnl.nnl_count[j[0]]
          n=config_nnl.nnl_count[j[1]]

          temp1={(A,m,B,C,n):triplets_count[(A,m,B,C,n)]+1}
          triplets_count.update(temp1)

          temp2={B:total_triplets[B]+1}
          total_triplets.update(temp2)
    self.triplets_AmBCn2count=triplets_count
    self.total_triplets_sym2count=total_triplets


class compute_ion_distribution:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Computing given ion distribution among the given atoms
  *
  * usage: config1_O_dist_among_formers=ganisetti_tools.compute_ion_distribution(cmd,config,config_nnl,O,given_formers_list)
  *
  * input:  config1	= ganisetti_tools.get_atoms_info_from_lammps(LAMMPS_DUMP_FILE)
  *         cmd     = ganisetti_tools.read_command_line(sys.argv)
  *
  * output: config1_triplets.triplets_AmBCn2count = {(A,m,B,C,n):count, etc.} # A-B-C
  *         config1_triplets.total_triplets_sym2count= {O:100,Si:120,etc.} O, Si are cetral atoms
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self, cmd, config, config_nnl,central_ion_sym,terminal_ions_sym_list):
    self.triplets_AmBCn2count={}
    self.bonds_in_triplets_sym2count={}
    atom_type_num2sym = cmd.atom_type_num2sym

    bonds_count={}
    for A in terminal_ions_sym_list:
      temp_bonds_count={A:0}
      bonds_count.update(temp_bonds_count)

    triplets_count={}
    for A in terminal_ions_sym_list:
      for m in range(config_nnl.max_nnl_each_atom_type_sym[A] + 1):
        for C in terminal_ions_sym_list:
          for n in range(config_nnl.max_nnl_each_atom_type_sym[C] + 1):
            temp_triplets_count={(A,m,central_ion_sym,C,n):0}
            triplets_count.update(temp_triplets_count)

    total_triplets=0
    for i in config.id:
      B=atom_type_num2sym[config.type[i]]
      if  B == central_ion_sym:
        # counting the triplets in which B is at the centre while A and C are at the two terminals i.e, A-B-C
        if config_nnl.nnl_count[i] > 1:
          for j in list( itertools.combinations(config_nnl.nnl[i], 2) ): # j=all combinations of (A and C) in A-B-C
            A = atom_type_num2sym[config.type[j[0]]]
            C = atom_type_num2sym[config.type[j[1]]]
            if A in terminal_ions_sym_list and C in terminal_ions_sym_list:
              m = config_nnl.nnl_count[j[0]]
              n = config_nnl.nnl_count[j[1]]
              temp_bonds_count = {A: bonds_count[A] + 1}
              bonds_count.update(temp_bonds_count)
              temp_bonds_count = {C: bonds_count[C] + 1}
              bonds_count.update(temp_bonds_count)
              if A == B and m == n:
                temp_triplets_count = {(A, m, B, C, n): triplets_count[(A, m, B, C, n)] + 1}
                triplets_count.update(temp_triplets_count)
              else :
                temp_triplets_count = {(A, m, B, C, n): triplets_count[(A, m, B, C, n)] + 1}
                triplets_count.update(temp_triplets_count)
                temp_triplets_count = {(C, n, B, A, m): triplets_count[(C, n, B, A, m)] + 1}
                triplets_count.update(temp_triplets_count)
              total_triplets=total_triplets+1

      self.bonds_in_triplets_sym2count  = bonds_count
      self.triplets_AmBCn2count         = triplets_count
      self.total_triplets               = total_triplets
