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