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

class banner:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Sudheer banner
  * This is a simple and little fun thing created to display SUDHEER on the screen
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self):
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
    self.nnl			= {}
    self.nnl_count		= {}
    self.nnl_type_num		= {}
    self.nnl_type_sym           = {}
    self.nnl_each_pair_distance = {}
    MAX_ATOM_NUMBER=max(config.id)
    nnl=[[-1 for j in range(MAX_NEIGHBOURS)] for i in range(MAX_ATOM_NUMBER+1)]
    nnl_count=[0 for i in range(MAX_ATOM_NUMBER+1)]
    #storing atom types of each neighbour in the nnl_types for getting things easy to work with
    #nnl_type=[[-1 for j in range(MAX_NEIGHBOURS)] for i in range(MAX_ATOM_NUMBER+1)]
    #store each pair distance to use them in later
    nnl_each_pair_distance=[[-1 for j in range(MAX_NEIGHBOURS)] for i in range(MAX_ATOM_NUMBER+1)]

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
        self.nnl[i]			= temp1
        self.nnl_type_num[i]		= temp2
        self.nnl_type_sym[i]		= temp4
        self.nnl_each_pair_distance[i]	= temp3
        self.nnl_count[i]		= nnl_count[i]

def write_imd_header(output,box,rc,atom_type_sym2num,atom_type_num2sym):
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
  * usage: environment1=ganisetti_tools.compute_each_atom_environment(config,config_nnl,atom_type_sym2num)
  * where config                = ganisetti_tools.get_atoms_info_from_lammps(01.dump)
  *       config_nnl            = ganisetti_tools.compute_nnl(config,rc,atom_type_num2sym)
  *       atom_type_sym2num     = {"O":"1", "Si":"2", "Al":"3", etc,.}
  *
  * output: config1.environment_sym
  *       : config1.environment_num
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,config,config_nnl,atom_type_sym2num):
    self.environment_sym={}
    for i in config.id:
      temp_each_atom_type_count={}
      #temp1={'id':i}
      #temp_each_atom_type_count.update(temp1)
      for j in atom_type_sym2num.keys():
        temp1={j:config_nnl.nnl_type_sym[i].count(j)}
        temp_each_atom_type_count.update(temp1)
      self.environment_sym[i] = temp_each_atom_type_count

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
  * output: config1.coord_sym
  *       : config1.coord_num
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,config,config_nnl,atom_type_num2sym):
    self.coord_num={}
    self.coord_sym={}
    temp1_count={}
    temp2_count={}
    atoms_count={}
    total_coord={}
    for i in atom_type_num2sym.keys():
      temp1_count={i:0}
      atoms_count.update(temp1_count)
      total_coord.update(temp1_count)
    for i in config.id:
      temp1_count={config.type[i]:atoms_count[config.type[i]]+1}
      atoms_count.update(temp1_count)
      temp2_count={config.type[i]:total_coord[config.type[i]]+config_nnl.nnl_count[i]}
      total_coord.update(temp2_count)
    for i in atom_type_num2sym.keys():
      self.coord_num[i]=float(total_coord[i]/atoms_count[i])
      self.coord_sym[atom_type_num2sym[i]]=float(total_coord[i]/atoms_count[i])

