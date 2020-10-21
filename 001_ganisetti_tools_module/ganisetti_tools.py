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
import datetime
import time
import git
import os

__author__  = "Sudheer Ganisetti"
__version__ = "73"
__email__   = "sudheerganisetti@gmail.com"
__status__  = "under preparation..."

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

def get_ganisetti_tools_version():
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * The ganisetti_tools version is defined based on the git commit of ganisetti_tools
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  time_now              = datetime.datetime.now()
  output1 = open("version_and_date_info.txt", 'w')
  output1.write("Author : Sudheer Ganisetti\n")
  output1.write("Date   : %s\n" % (time_now))
  output1.write("command used:\n")
  for i in sys.argv:
    output1.write("%s  " % (str(i)))
  output1.write("\n\ninformation of module ganisetti_tools\n")

  if os.path.isdir("/data/ganisetti/TOOLS/github_repos/md_tools"):
    git_repo              = git.Repo("/data/ganisetti/TOOLS/github_repos/md_tools")
    git_log               = git_repo.heads.master.log()
    git_last_commit       = git_log[-1]
    git_last_commit_time  = datetime.datetime.fromtimestamp(int(git_last_commit[3][0])).strftime("%a %B %d %Y %I:%M:%S")
    git_commits_number    = len(git_log)

    output1.write("version number   : %d\n" %(git_commits_number))
    output1.write("version time     : %s\n" %(str(git_last_commit_time)))
    output1.write("version stamp    : %s\n" %(str(git_last_commit[1])))
    output1.write("version last %s\n" %(git_last_commit[4]))
  else :
    output1.write("git version details are not avialable\n")
  output1.close()

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
  *             config1.imd_box            = np.array([[box_xx,box_yy,box_zz],[box_yx,box_yy,box_yz],[box_zx,box_zy,box_zz]])
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
            break
          else:
            imd_header_properties1=data
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
        if data[0] == "#C":
          imd_header_properties2=data
      count = count + 1
  
    '''
    # The following loop is to generalize the reading imd file but did not complete
    imd_header_index2property={}
    imd_header_property2index={}
    for i in range(1,len(imd_header_properties2)):
      temp1={i-1:imd_header_properties2[i]}
      imd_header_index2property.update(temp1)
      temp1={imd_header_properties2[i]:i-1}
      imd_header_property2index.update(temp1)
    if imd_header_index2property[imd_header_property2index[]] != "number":
      warning_alarm=1
      warnin_message="WARNING: please make sure the first column of input file is 'number' "
    if imd_header_index2property[1] != "mass":
    
    imd_header_properties3={} 
    imd_header_properties4={}
    imd_header_properties5={}
      #temp1=imd_header_properties1_summation[i-1]+int(imd_header_properties1[i])
      #imd_header_properties1_summation.append(temp1)
    #id_index              = 0  # the first column of imd file should always be atom id 
    #type_index            = 1  # the second column of imd file should always be atom type
    #if imd_header_properties1[4] == "1":
    #  mass_index          = count1 + 1
    #  if imd_header_properties2[3] != "mass":
    #    print("#############################################################################################################")
    #    print("#                                                                                                           #")
    #    print("#  WARNING: It seems the third column is not mass but considered it as mass so becareful with the analysis  #")
    #    print("#                                                                                                           #")
    #    print("#############################################################################################################")
    #posx_index            = 
    #posy_index            = 
    #posz_index            = 
    #extra_property1_index =
    #extra_property2_index =
    #extra_property3_index =
    #extra_property4_index =
    '''
    box_xx_lo = min(box_xx,box_yx,box_zx)
    box_yy_lo = min(box_xy,box_yy,box_zy)
    box_zz_lo = min(box_xz,box_yz,box_zz)
    box_xx_hi = max(box_xx,box_yx,box_zx)
    box_yy_hi = max(box_xy,box_yy,box_zy)
    box_zz_hi = max(box_xz,box_yz,box_zz)
    self.box_xx          = box_xx + min(box_yx,box_zx)
    self.box_yy          = box_yy + min(box_xy,box_zy)
    self.box_zz          = box_zz + min(box_xz,box_yz)
    self.imd_box         = np.array([[box_xx,box_xy,box_xz],[box_yx,box_yy,box_yz],[box_zx,box_zy,box_zz]])
    self.box             = np.array([[box_xx_lo,box_xx_hi],[box_yy_lo,box_yy_hi],[box_zz_lo,box_zz_hi]])
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
  *		config1.box	     	        = np.array([[box_xx_lo,box_xx_hi],[box_yy_lo,box_yy_hi],[box_zz_lo,box_zz_hi]])
  *   config1.id[id]     	      = id
  *   config1.type[id]   	      = int(i[1])
  *   config1.charge[id] 	      = float(i[2])
  *   config1.posx[id]   	      = float(i[3])
  *   config1.posy[id]   	      = float(i[4])
  *   config1.posz[id]   	      = float(i[5])
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
  * output: config1.nnl                           = {0:[1,2,3,4]}
  *         config1.nnl_count                     = {0:4}
  *         config1.nnl_type_num                  = {0:[1,1,1,1]}
  *         config1.nnl_type_sym                  = {0:[Si,Si,Si,Si]}
  *         config1.nnl_each_pair_distance        = {0:[1.58,1.62,1.59,1.60]}
  *         config1.max_nnl_each_atom_type_sym    = {O:4,Si:4,Al:6,Ca:8}
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

def write_imd_header_basic(output,box):
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Writing basic version of IMD header into given output file
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  output.write("#F A 1 1 3 0 0\n")
  output.write("#C number type x y z \n")
  output.write("#X %lf 0.0 0.0 \n" %(box[0][1]-box[0][0]))
  output.write("#Y 0.0 %lf 0.0 \n" %(box[1][1]-box[1][0]))
  output.write("#Z 0.0 0.0 %lf \n" %(box[2][1]-box[2][0]))
  output.write("#E \n")


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

def write_imd_header_custom_property(output,box,property_name):
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Writing the IMD header into given output file with a custom property name
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  output.write("#F A 1 1 3 1 0\n")
  output.write("#C number type x y z %s \n" %(str(property_name)))
  output.write("#X %lf 0.0 0.0 \n" %(box[0][1]-box[0][0]))
  output.write("#Y 0.0 %lf 0.0 \n" %(box[1][1]-box[1][0]))
  output.write("#Z 0.0 0.0 %lf \n" %(box[2][1]-box[2][0]))
  output.write("#E \n")

def write_imd_atom(output,atom_id,config,config_nnl):
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Writing the IMD atom data into a given output file
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  output.write("%d  %d  %lf  %lf  %lf %d \n" %(atom_id,config.type[atom_id],config.posx[atom_id],config.posy[atom_id],config.posz[atom_id],config_nnl.nnl_count[atom_id]))

def write_imd_atom_custom_property(output,atom_id,config,prop):
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Writing the IMD atom data into a given output file with custom property 
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  output.write("%d  %d  %lf  %lf  %lf %s \n" %(atom_id,config.type[atom_id],config.posx[atom_id],config.posy[atom_id],config.posz[atom_id],str(prop)))

def compute_pair_distance(i,j,config):
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Compute distance between two atoms
  *
  * usage: pair_distance = ganisetti_tools.compute_pair_distance(i,j,config1)
  * 
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  rx = np.linalg.norm(config.posx[i]-config.posx[j])
  ry = np.linalg.norm(config.posy[i]-config.posy[j])
  rz = np.linalg.norm(config.posz[i]-config.posz[j])
  box_xl=(config.box[0][1]-config.box[0][0])
  box_yl=(config.box[1][1]-config.box[1][0])
  box_zl=(config.box[2][1]-config.box[2][0])
  if rx > box_xl/2.0:
    rx = box_xl-rx
  if ry > box_yl/2.0:
    ry=box_yl-ry
  if rz > box_zl/2.0:
    rz=box_zl-rz
  return float(np.linalg.norm([rx,ry,rz]))

class compute_each_atom_environment:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Compute each atom's environment
  * This class helps to compute the environment of each atom
  * and store the data into an object
  *
  * usage: efg_env=ganisetti_tools.compute_each_atom_environment(config,config_nnl,atom_type_sym2num,atom_type_num2sym)
  * where config                = ganisetti_tools.get_atoms_info_from_lammps(01.dump)
  *       config_nnl            = ganisetti_tools.compute_nnl(config,rc,atom_type_num2sym)
  *       atom_type_sym2num     = {"O":"1", "Si":"2", "Al":"3", etc,.}
  *       atom_type_num2sym     = {1:"O", 2:"Si", 3:"Al", etc,.}
  *
  * output: deprecated => config1.environment_atomnum2sym[0] = {'O':0,'Si':2,'Al':1}
  *         deprecated => config1.env_atomnum2nnlsymandcoord[0] = (('Si',4),('Al',4))
  *         deprecated => config1.env_atomnum2countofnnlsymandcoord[0] = (1,1)
  *         cfg_env.env_atomid[0]={('Si',4):1,('Al',4):1}  # which means Si[4]-O-Al[4]
  *         cfg_env.all_possible_local_env_and_counts={(O,(Si,2),(Al,2)):154} # the number of O which has 2Si+2Al = 154
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
    self.all_possible_local_env_and_counts ={}

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

    local_env_count = {}
    # Since creating a complete matrix to store the environment count is very memory expensive
    # only the necessary elements of the matrix is created
    # for this we need to run the loop two times, first time to initilized it and second time for counting it
    for i in config.id:
      local_env = [atom_type_num2sym[config.type[i]]]
      for j in atom_type_sym2num.keys():
        local_env.append((j, list(config_nnl.nnl_type_sym[i]).count(j)))
      temp1 = {tuple(local_env): 0}
      local_env_count.update(temp1)

    for i in config.id:
      local_env = [atom_type_num2sym[config.type[i]]]
      for j in atom_type_sym2num.keys():
        local_env.append((j, list(config_nnl.nnl_type_sym[i]).count(j)))
      temp1 = {tuple(local_env): local_env_count[tuple(local_env)] + 1}
      local_env_count.update(temp1)
    self.all_possible_local_env_and_counts=local_env_count
    # In the future, I would like to rewrite this as follows
    # all_possible_local_env_and_counts_new={O:{((Si,1),(Al,1)):50,((Si,1),(Al,1),(Na,1)):200}}

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
  *         config1.avg_coord_of_each_type={(Si,O):3.0,(Si,F):1.0} total_coord of Si =4 of which O is 3 and F is 1
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,config,config_nnl,atom_type_num2sym):
    self.avg_coord_num          = {}
    self.avg_coord_sym          = {}
    self.avg_coord_of_each_type = {}
    atoms_count                 = {}
    total_coord                 = {}
    nnl_individual_coord_count  = {}
    avg_coord_of_each_type      = {}
    for i in atom_type_num2sym.keys():
      i_sym = atom_type_num2sym[i]
      temp1_count={i:0}
      atoms_count.update(temp1_count)
      total_coord.update(temp1_count)
      for j in range(config_nnl.max_nnl_each_atom_type_sym[atom_type_num2sym[i]]+1):
        temp2_count={(atom_type_num2sym[i],j):0}
        nnl_individual_coord_count.update(temp2_count)
      for j_sym in atom_type_num2sym.values():
        temp_avg_coord_of_each_type={(i_sym,j_sym):0}
        avg_coord_of_each_type.update(temp_avg_coord_of_each_type)

    # Main loop
    for i in config.id:
      temp1_count={config.type[i]:atoms_count[config.type[i]]+1}
      atoms_count.update(temp1_count)
      temp2_count={config.type[i]:total_coord[config.type[i]]+config_nnl.nnl_count[i]}
      total_coord.update(temp2_count)
      i_sym=atom_type_num2sym[config.type[i]]
      temp4_count={(i_sym,config_nnl.nnl_count[i]):nnl_individual_coord_count[(i_sym,config_nnl.nnl_count[i])]+1}
      nnl_individual_coord_count.update(temp4_count)
      for j in config_nnl.nnl[i]:
        j_sym=atom_type_num2sym[config.type[j]]
        temp_avg_coord_of_each_type={(i_sym,j_sym):avg_coord_of_each_type[(i_sym,j_sym)]+1}
        avg_coord_of_each_type.update(temp_avg_coord_of_each_type)

    for i in atom_type_num2sym.keys():
      self.avg_coord_num[i]=float(total_coord[i]/atoms_count[i])
      self.avg_coord_sym[atom_type_num2sym[i]]=float(total_coord[i]/atoms_count[i])
      i_sym=atom_type_num2sym[i]
      for j_sym in atom_type_num2sym.values():
        temp_avg_coord_of_each_type={(i_sym,j_sym):float(avg_coord_of_each_type[(i_sym,j_sym)]/atoms_count[i])}
        avg_coord_of_each_type.update(temp_avg_coord_of_each_type)

    self.individual_coord_sym=nnl_individual_coord_count
    self.avg_coord_of_each_type=avg_coord_of_each_type


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
    known_cations = ['Si','Al','P','Na','Ca','Mg','Sr','Li','K','V'] # The cations list should be updated if you want to add a new cation
    known_anions  = ['O','F']                      # The anions list should be updated if you want to add a new anion
    known_all_atom_types = known_cations + known_anions
    known_network_formers=['Si','Al','P','V']
    known_modifiers      =['Na','Ca','Mg','Sr','Li','K']

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
      print("correct usage is: ./this_program  parameter_file")
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
    self.Q_status_atomid2num=Q_status_atomid2num

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
  *         config1_QSi.Q_status_atomid2num = {0:-1,1:-1,2:1,3:-1} 0,1,3 are other atom tpyes, 2 is requred type
  *         config1_QSi.Q_4CoordFormers_list = {(Si,0):[1,2,3,etc.],(Si,1):[1,2,3,etc.]}
  *         config1_QSi.Q_non4CoordFormers_list = {(Si,1):[1,2,3, etc],(Al,1)}
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self, cmd, config, config_nnl, config_env):
    self.Q_status_atomid2num            = {}
    self.Q_summary                      = {}
    self.Q_4CoordFormers_list           = {}
    self.Q_non4CoordFormers_list        = {}
    self.Q_4CoordFormers_id2sym_list    = {}
    self.Q_non4CoordFormers_id2sym_list = {}
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
        NetworkFormersList_4coord_sym   = []  # for atom i
        NetworkFormersList_non4coord_sym= []  # for atom i
        for j in config_nnl.nnl[i]: # j = Q neighbouring anions: O1,O2,O3,O4
          Anion_4CoordinatedNetworkFormers_count    = 0 # for atom j
          Anion_non4CoordinatedNetworkFormers_count = 0 # for atom j
          for k in config_nnl.nnl[j]: # k = Q neighboring cations Si,Si (including the Q cation)
            if config.type[k] in given_formers_sym2num.values():
              if config_nnl.nnl_count[k] == 4:
                Anion_4CoordinatedNetworkFormers_count=Anion_4CoordinatedNetworkFormers_count+1
                NetworkFormersList_4coord.append(k)   # This is for Si[4]-O-Si[4], Si[4]-O-Al[4], etc
                if k != i:
                  NetworkFormersList_4coord_sym.append(config.type[k])
              else :
                Anion_non4CoordinatedNetworkFormers_count=Anion_non4CoordinatedNetworkFormers_count+1
                NetworkFormersList_non4coord.append(k) # This is for Si[4]-O-Si[5], Si[4]-O-Al[5], etc
                if k != i:
                  NetworkFormersList_non4coord_sym.append(config.type[k])
          if Anion_4CoordinatedNetworkFormers_count == 2:
            BridgingAnions_count            = BridgingAnions_count+1
          elif Anion_4CoordinatedNetworkFormers_count == 1:
            NonBridgingAnions_count         = NonBridgingAnions_count+1
          elif Anion_4CoordinatedNetworkFormers_count == 3:
            Triplet_NonBridgingAnions_count = Triplet_NonBridgingAnions_count+1
        temp1={(atom_type_num2sym[config.type[i]],BridgingAnions_count):Q_summary[(atom_type_num2sym[config.type[i]],BridgingAnions_count)]+1}
        Q_summary.update(temp1)

        temp1={i:BridgingAnions_count}
        Q_status_atomid2num.update(temp1)

        NetworkFormersList_4coord=NetworkFormersList_4coord + Q_4CoordFormers_list[(atom_type_num2sym[config.type[i]],BridgingAnions_count)]
        NetworkFormersList_4coord=list(dict.fromkeys(NetworkFormersList_4coord))
        temp1 = {(atom_type_num2sym[config.type[i]],BridgingAnions_count):NetworkFormersList_4coord}
        Q_4CoordFormers_list.update(temp1)

        NetworkFormersList_non4coord = NetworkFormersList_non4coord + Q_non4CoordFormers_list[(atom_type_num2sym[config.type[i]],BridgingAnions_count)]
        NetworkFormersList_non4coord = list(dict.fromkeys(NetworkFormersList_non4coord))
        temp1 = {(atom_type_num2sym[config.type[i]], BridgingAnions_count): NetworkFormersList_non4coord}
        Q_non4CoordFormers_list.update(temp1)
        self.Q_4CoordFormers_id2sym_list[i]     = NetworkFormersList_4coord_sym
        self.Q_non4CoordFormers_id2sym_list[i]  = NetworkFormersList_non4coord_sym
    self.Q_summary = Q_summary
    self.Q_status_atomid2num=Q_status_atomid2num
    self.Q_4CoordFormers_list = Q_4CoordFormers_list
    self.Q_non4CoordFormers_list=Q_non4CoordFormers_list

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


class compute_ion_distribution_deprecated:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Computing given ion distribution among the given atoms
  *
  * usage: config1_O_dist_among_formers=ganisetti_tools.compute_ion_distribution(cmd,config,config_nnl,O,given_formers_list)
  *
  * input:  config1	= ganisetti_tools.get_atoms_info_from_lammps(LAMMPS_DUMP_FILE)
  *         cmd     = ganisetti_tools.read_command_line(sys.argv)
  *
  * output: config1_O_dist_among_formers.triplets_AmBCn2count = {(A,m,B,C,n):count, etc.} # A-B-C
  *         config1_O_dist_among_formers.total_triplets= total number of triplets with O at centre
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


class compute_anions_distribution:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Computing all anions distribution among the given network formers
  *
  * usage: config1_ad=ganisetti_tools.compute_anions_distribution(cmd,config,config_nnl,O,given_formers_list)
  *
  * input:  config1	= ganisetti_tools.get_atoms_info_from_lammps(LAMMPS_DUMP_FILE)
  *         cmd     = ganisetti_tools.read_command_line(sys.argv)
  *
  * output: config1_ad.triplets_AmBCn2count = {(A,m,B,C,n):count, etc.} # gives the count of B atoms in as A-B-C
  *         config1_ad.total_triplets= total number of triplets with O at centre
  *         config1_ad.BA_4CoordFormer_id2list = {'O':[1,2,5,10,etc.]}  # 1 = Si[4], 2 = Si[4] => Si[4]-O-Si[4]
  *         config1_ad.BA_non4CoordFormer_id2list = {'O':[3,6,29,etc.]} # 3 = Al[5], 6 = Si[4] => Al[5]-O-Si[4]
  *         config1_ad.modifiers_id2list = {'O':[23,45,60,etc]} # if O is BA then these are charge compensators
  *                                                             # if O is NBA then these are network modifiers
  *         config1_ad.tri_cluster_former_id2list = {O:[21,27,29]} # tricluster atoms
  *         config1_ad.NBA_former_id2list = {'O':[33,37,56]} # these are network formers in Si-NBO, Al-NBO, P-NBO etc
  *         config1_ad.anions_with_zero_formers_id2list ={'O':[1,11,2]} these are list of anions connected to no former
  *         config1_ad.anions_of_any_other_type_id2list ={'O':[22,13,66]} = Any_other_type_of_anion
  *
  *         config1_ad.total_anions_of_BA_4CoordFormer_sym2count    ={'O':20} total anions in Si[4]-O-Si[4]
  *         config1_ad.total_anions_of_BA_non4CoordFormer_sym2count ={'O':31} total_anions in Al[5]-O-Si[4]
  *         config1_ad.total_anions_of_tri_clusters_sym2count       ={'O':24} total anions in triclusters
  *         config1_ad.total_anions_of_NBA_former_sym2count         ={'O':54} total anions in Si-NBO,Al-NBO, P-NBO
  *         config1_ad.total_anions_with_zero_formers_sym2count     ={'O':76} total anions_connected_to_no_formers
  *         config1_ad.total_anions_of_any_other_type_sym2count     ={'O':93} total of Any_other_type_anions
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self, cmd, config, config_nnl):
    self.triplets_AmBCn2count={}
    #self.bonds_in_triplets_sym2count={}

    atom_type_sym2num         = cmd.atom_type_sym2num
    atom_type_num2sym         = cmd.atom_type_num2sym
    given_anions_sym2num      = cmd.given_anions_sym2num
    given_cations_sym2num     = cmd.given_cations_sym2num
    given_formers_sym2num     = cmd.given_formers_sym2num
    given_modifiers_sym2num   = cmd.given_modifiers_sym2num

    '''
    bonds_count={}
    for A in terminal_ions_sym_list:
      temp_bonds_count={A:0}
      bonds_count.update(temp_bonds_count)
    '''
    triplets_count={}
    for B in given_anions_sym2num.keys():
      for A in given_formers_sym2num.keys():
        for m in range(config_nnl.max_nnl_each_atom_type_sym[A] + 1):
          for C in given_formers_sym2num.keys():
            for n in range(config_nnl.max_nnl_each_atom_type_sym[C] + 1):
              temp_triplets_count={(A,m,B,C,n):0}
              triplets_count.update(temp_triplets_count)

    BA_4Coord_former                    = {}
    BA_non4Coord_former                 = {}
    modifiers                           = {}
    tri_cluster_former                  = {}
    NBA_former                          = {}
    anions_connected_to_no_former       = {}
    Any_other_type_anion                = {}
    BA_4Coord_former_anion_count        = {}
    BA_non4Coord_former_anion_count     = {}
    NBA_former_anion_count              = {}
    tri_cluster_former_anion_count      = {}
    anions_connected_to_no_former_count = {}
    Any_other_type_anion_count          = {}
    for i in given_anions_sym2num.keys():
      temp1={i:0}
      BA_4Coord_former_anion_count.update(temp1)
      BA_non4Coord_former_anion_count.update(temp1)
      NBA_former_anion_count.update(temp1)
      tri_cluster_former_anion_count.update(temp1)
      anions_connected_to_no_former_count.update(temp1)
      Any_other_type_anion_count.update(temp1)

    temp1_debug=0
    for i in config.id:
      B=atom_type_num2sym[config.type[i]]
      if B in given_anions_sym2num.keys(): # i = B = anion
        temp1_debug = temp1_debug +1
        temp1_4Coord_former     = []
        temp1_non4Coord_former  = []
        temp1_modifier          = []
        temp3_4Coord_former     = {i:[]}
        temp3_non4Coord_former  = {i:[]}
        temp3_modifier          = {i:[]}
        temp4_4Coord_former     = 0
        temp4_non4Coord_former  = 0
        temp4_modifier          = 0
        for j in config_nnl.nnl[i]: # j = neighbouring cations [23,48,88,24]
          if config.type[j] in given_formers_sym2num.values():
            if config_nnl.nnl_count[j] == 4:
              temp1_4Coord_former.append(j)
              temp2_4Coord_former={i:temp1_4Coord_former}
              temp3_4Coord_former.update(temp2_4Coord_former)
              temp4_4Coord_former=temp4_4Coord_former+1
            else :
              temp1_non4Coord_former.append(j)
              temp2_non4Coord_former={i:temp1_non4Coord_former}
              temp3_non4Coord_former.update(temp2_non4Coord_former)
              temp4_non4Coord_former=temp4_non4Coord_former+1
          else:
            temp1_modifier.append(j)
            temp2_modifier={i:temp1_modifier}
            temp3_modifier.update(temp2_modifier)
            temp4_modifier=temp4_modifier+1

        if temp4_4Coord_former == 2 and temp4_non4Coord_former == 0:
          BA_4Coord_former.update(temp3_4Coord_former)
          temp1={B:BA_4Coord_former_anion_count[B]+1}
          BA_4Coord_former_anion_count.update(temp1)
        elif temp4_4Coord_former + temp4_non4Coord_former == 2:
          temp5=[]
          temp6=temp3_4Coord_former[i]
          for j in temp6:
            temp5.append(j)
          temp6=temp3_non4Coord_former[i]
          for j in temp6:
            temp5.append(j)
          temp6={i:temp5}
          BA_non4Coord_former.update(temp6)
          temp1={B:BA_non4Coord_former_anion_count[B]+1}
          BA_non4Coord_former_anion_count.update(temp1)
        elif temp4_4Coord_former + temp4_non4Coord_former == 3:
          temp5=[]
          temp6=temp3_4Coord_former[i]
          for j in temp6:
            temp5.append(j)
          temp6=temp3_non4Coord_former[i]
          for j in temp6:
            temp5.append(j)
          temp6={i:temp5}
          tri_cluster_former.update(temp6)
          temp1={B:tri_cluster_former_anion_count[B]+1}
          tri_cluster_former_anion_count.update(temp1)
        elif temp4_4Coord_former + temp4_non4Coord_former == 1:
          temp5=[]
          temp6=temp3_4Coord_former[i]
          for j in temp6:
            temp5.append(j)
          temp6=temp3_non4Coord_former[i]
          for j in temp6:
            temp5.append(j)
          temp6={i:temp5}
          NBA_former.update(temp6)
          temp1={B:NBA_former_anion_count[B]+1}
          NBA_former_anion_count.update(temp1)
        elif temp4_4Coord_former + temp4_non4Coord_former == 0:
          temp6={i:i}
          anions_connected_to_no_former.update(temp6)
          temp1={B:anions_connected_to_no_former_count[B]+1}
          anions_connected_to_no_former_count.update(temp1)
        else :
          temp6={i:i}
          Any_other_type_anion.update(temp6)
          temp1={B:Any_other_type_anion_count[B]+1}
          Any_other_type_anion_count.update(temp1)
        modifiers.update(temp3_modifier)

        # the following loop is to compute triplets i.e Si[4]-O-Al[4], Si[4]-O-Al[5], etc.
        # Do not include triclusters
        if temp4_4Coord_former + temp4_non4Coord_former != 3:
          for j in list(itertools.combinations(config_nnl.nnl[i], 2)):  # j=all combinations of (A and C) in A-B-C
            A = atom_type_num2sym[config.type[j[0]]]
            C = atom_type_num2sym[config.type[j[1]]]
            if A in given_formers_sym2num.keys() and C in given_formers_sym2num.keys():
              m = config_nnl.nnl_count[j[0]]
              n = config_nnl.nnl_count[j[1]]
              #temp_bonds_count = {A: bonds_count[A] + 1}
              #bonds_count.update(temp_bonds_count)
              #temp_bonds_count = {C: bonds_count[C] + 1}
              #bonds_count.update(temp_bonds_count)
              if A == C and m == n:
                temp_triplets_count = {(A, m, B, C, n): triplets_count[(A, m, B, C, n)] + 1}
                triplets_count.update(temp_triplets_count)
              else:
                temp_triplets_count = {(A, m, B, C, n): triplets_count[(A, m, B, C, n)] + 1}
                triplets_count.update(temp_triplets_count)
                temp_triplets_count = {(C, n, B, A, m): triplets_count[(C, n, B, A, m)] + 1}
                triplets_count.update(temp_triplets_count)
              #total_triplets = total_triplets + 1

    self.BA_4CoordFormer_id2list                      = BA_4Coord_former
    self.BA_non4CoordFormer_id2list                   = BA_non4Coord_former
    self.modifiers_id2list                            = modifiers
    self.tri_cluster_former_id2list                   = tri_cluster_former
    self.NBA_former_id2list                           = NBA_former
    self.anions_with_zero_formers_id2list             = anions_connected_to_no_former
    self.anions_of_any_other_type_id2list             = Any_other_type_anion

    self.total_anions_of_BA_4CoordFormer_sym2count    = BA_4Coord_former_anion_count
    self.total_anions_of_BA_non4CoordFormer_sym2count = BA_non4Coord_former_anion_count
    self.total_anions_of_tri_clusters_sym2count       = tri_cluster_former_anion_count
    self.total_anions_of_NBA_former_sym2count         = NBA_former_anion_count
    self.total_anions_with_zero_formers_sym2count     = anions_connected_to_no_former_count
    self.total_anions_of_any_other_type_sym2count     = Any_other_type_anion_count

    self.triplets_AmBCn2count                         = triplets_count


class compute_atoms_density:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * This is my own idea: -- Sudheer
  * Compute the atoms density at each point in space based on the surrounding atoms
  *
  * usage: config1=ganisetti_tools.compute_atoms_density(config,voxel_length,voxels_count_for_smoothing,given_atom_ids,atom_type_num2sym)
  * where config    = ganisetti_tools.get_atoms_info_from_lammps(01.dump)
  *       voxel_length = 2.5/5
  *       voxels_count_for_smoothing = 3
  *       atom_type_num2sym = {"1":"O", "2":"Si", "3":"Al", etc,.}
  *
  * output: config1.density
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,config,rcut,vox_smooth,given_atoms):


    # Devide the box into voxels
    MaxCells_X=int(config.box_xx/rcut)
    MaxCells_Y=int(config.box_yy/rcut)
    MaxCells_Z=int(config.box_zz/rcut)

    # LEC= Length of Each Cell
    LEC_X=float(config.box_xx/MaxCells_X)
    LEC_Y=float(config.box_yy/MaxCells_Y)
    LEC_Z=float(config.box_zz/MaxCells_Z)

    # define parameters
    self.atoms_density={}
    temp_atoms_density={}
    temp_factor1=1.73205080756888*vox_smooth  # sqrt(3)*vox_smooth
    # store atom positions locally
    atom_posx=config.posx
    atom_posy=config.posy
    atom_posz=config.posz
    # move the box to center and correspondingly atoms if the corner of the box is not at center
    for i in given_atoms:
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
      #assign the atoms into respective cell number  #assign the atoms into respective cell number; cell numbers starts with 0
      cellX=int(atom_posx[i]/LEC_X)
      cellY=int(atom_posy[i]/LEC_Y)
      cellZ=int(atom_posz[i]/LEC_Z)

      # add 1 to the voxel that belong to the atom 'i' and also the surrounding 'vox_smooth' number of voxels 
      for ix in range(cellX-vox_smooth,cellX+vox_smooth+1):
        raw_ix=ix-cellX
        # Dealing with periodic boundary conditions
        if ix < 0:
          ix = MaxCells_X+ix
        if ix >= MaxCells_X:
          ix = ix - MaxCells_X

        for iy in range(cellY-vox_smooth,cellY+vox_smooth+1):
          raw_iy=iy-cellY
          # Dealing with periodic boundary conditions
          if iy < 0:
            iy = MaxCells_Y + iy
          if iy >= MaxCells_Y:
            iy = iy - MaxCells_Y

          for iz in range(cellZ-vox_smooth,cellZ+vox_smooth+1):
            raw_iz=iz-cellZ
            # Dealing with periodic boundary conditions
            if iz < 0:
              iz = MaxCells_Z + iz
            if iz >= MaxCells_Z:
              iz = iz - MaxCells_Z

            temp1=(ix,iy,iz) in temp_atoms_density.keys()
            if temp1 == False:
              temp2={(ix,iy,iz):0.0}
              temp_atoms_density.update(temp2)
            temp_factor2=pow(raw_ix,2)+pow(raw_iy,2)+pow(raw_iz,2)
            if temp_factor2 == 0:
              temp_factor3=temp_factor1
            else :
              temp_factor3=temp_factor1/np.sqrt(temp_factor2)

            temp3={(ix,iy,iz):temp_atoms_density[(ix,iy,iz)]+temp_factor3}
            temp_atoms_density.update(temp3)
            
    temp2_atoms_density={}
    for i,j,k in temp_atoms_density.keys():
      temp1={(i*LEC_X,j*LEC_Y,k*LEC_Z):temp_atoms_density[(i,j,k)]}
      temp2_atoms_density.update(temp1)

    self.atoms_density = temp2_atoms_density
          
class compute_clustered_channels_based_on_voxels:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * This is my own idea: -- Sudheer
  * Compute the channels by clustering (grouping) the given atoms
  *
  * usage: config1=ganisetti_tools.compute_clustered_channels(config,voxel_length,voxels_count_for_smoothing,given_atom_ids)
  * where config    = ganisetti_tools.get_atoms_info_from_lammps(01.dump)
  *       voxel_length = 2.5/5
  *       voxels_count_for_smoothing = 3
  *       atom_type_num2sym = {"1":"O", "2":"Si", "3":"Al", etc,.}
  *
  * output: config1.channels
  * 
  * IT IS NOT WORKING AS I EXPECTED STILL NEED TO WORK ON IT
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,config,rcut,vox_smooth,given_atoms):

    self.vox_smooth=vox_smooth
    # Devide the box into voxels
    self.MaxCells_X=int(config.box_xx/rcut)
    self.MaxCells_Y=int(config.box_yy/rcut)
    self.MaxCells_Z=int(config.box_zz/rcut)

    # LEC= Length of Each Cell
    self.LEC_X=float(config.box_xx/self.MaxCells_X)
    self.LEC_Y=float(config.box_yy/self.MaxCells_Y)
    self.LEC_Z=float(config.box_zz/self.MaxCells_Z)

    # define parameters
    cluster_count = 1
    all_neighbouring_cells_of_the_given_atoms=[]
    cluster_id = {}
    self.clusters_position_to_id={}
    self.clusters_id_to_position={}
    active_cell_neighbours={}

    # store atom positions locally
    atom_posx=config.posx
    atom_posy=config.posy
    atom_posz=config.posz
    # move the box to center and correspondingly atoms if the corner of the box is not at center
    for i in given_atoms:
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
      #assign the atoms into respective cell number  #assign the atoms into respective cell number; cell numbers starts with 0
      self.cellX=int(atom_posx[i]/self.LEC_X)
      self.cellY=int(atom_posy[i]/self.LEC_Y)
      self.cellZ=int(atom_posz[i]/self.LEC_Z)

      # compute neighboring cells of each atom, these cells are also referred as active cells
      # store them in a list, i.e, all_neighbouring_cells_of_the_given_atoms_original
      # to avoid working with large data of emty cells
      neighbours1=self.get_neighbour_cells()
      for i in neighbours1:
      	if i not in all_neighbouring_cells_of_the_given_atoms:
          all_neighbouring_cells_of_the_given_atoms.append(i)
    print("total active cells = %d " %(len(all_neighbouring_cells_of_the_given_atoms)))

    # initializing a dictinary for storing the active cell neighbors of all active cells 
    for i in all_neighbouring_cells_of_the_given_atoms:
      temp1={i:[]}
      active_cell_neighbours.update(temp1)
    start_time1=time.time()

    # for each active cell store the neighbors of only the active cells
    for i in all_neighbouring_cells_of_the_given_atoms:   # i = (ix,iy,iz) = cell position
      self.cellX=i[0]
      self.cellY=i[1]
      self.cellZ=i[2]
      neighbours1=self.get_neighbour_cells()
      neighbours1.remove(i)                                   # remove the own cell
      for j in neighbours1:
        if j in all_neighbouring_cells_of_the_given_atoms:
          temp1=active_cell_neighbours[i]
          temp1.append(j)
          temp2=list(set(temp1))
          temp3={i:temp2}
          active_cell_neighbours.update(temp3)		#active_cell_neighbors[(x1,y1,z1)]=[(x2,y2,z2),(x3,y3,z3)]
    end_time1=time.time()
    print("storing active cells is finished; time taken = %s " %(str(datetime.timedelta(seconds=end_time1-start_time1))))

    # main loop to find the clusters
    all_unclustered_cells=all_neighbouring_cells_of_the_given_atoms
    while len(all_unclustered_cells) != 0:
      collected_cells_of_the_cluster={}
      temp1={all_unclustered_cells[0]:1}
      collected_cells_of_the_cluster.update(temp1)
      while len(collected_cells_of_the_cluster) !=0 :
      	the_cell_of_the_cluster=list(collected_cells_of_the_cluster.keys())[0]
      	temp1={the_cell_of_the_cluster:cluster_count}
      	cluster_id.update(temp1)
      	all_unclustered_cells.remove(the_cell_of_the_cluster)
      	collected_cells_of_the_cluster.pop(the_cell_of_the_cluster)
      	for the_neighbour_of_the_cell_of_the_cluster in active_cell_neighbours[the_cell_of_the_cluster]:
      	  temp1={the_neighbour_of_the_cell_of_the_cluster:1}
      	  collected_cells_of_the_cluster.update(temp1)
      	temp1=active_cell_neighbours[the_cell_of_the_cluster]
      	temp1.append(the_cell_of_the_cluster)
      	for pair in list(itertools.combinations(temp1,2)):
      	  try:
      	    active_cell_neighbours[pair[0]].remove(pair[1])
      	    active_cell_neighbours[pair[1]].remove(pair[0])
      	  except:
      	    pass
      cluster_count=cluster_count+1

    '''
    while len(all_unclustered_cells) != 0:						# this loop repeats for each cluster
      first_cell_of_the_cluster=all_unclustered_cells[0]
      temp1={first_cell_of_the_cluster:cluster_count}
      cluster_id.update(temp1)
      all_unclustered_cells.remove(first_cell_of_the_cluster)	# remove the cell from the list of all unclustered cells if the cell is given with a channel id
      cluster1=[]
      for first_cell_neighbour in active_cell_neighbours[first_cell_of_the_cluster]:
      	active_cell_neighbours[first_cell_neighbour].remove(first_cell_of_the_cluster)
        cluster1.update(first_cell_neighbour)
      while len(cluster1) != 0:									# this loop repeats for each cell of a cluster
      	next_cell_of_the_cluster=cluster1[0]
      	temp1={next_cell_of_the_cluster:cluster_count}
      	cluster_id.update(temp1)
      	all_unclustered_cells.remove(next_cell_of_the_cluster)
      	cluster1.remove(next_cell_of_the_cluster)
      	for next_cell_neighbour in active_cell_neighbours[next_cell_of_the_cluster]:
      	  active_cell_neighbours[next_cell_neighbour].remove(next_cell_of_the_cluster)
      	  cluster1.update(next_cell_neighbour)
    '''
    '''
    while len(all_unclustered_cells) != 0:
      collected_cells_of_the_cluster=[]
      collected_cells_of_the_cluster.append(all_unclustered_cells[0])
      while len(collected_cells_of_the_cluster) !=0 :
      	the_cell_of_the_cluster=collected_cells_of_the_cluster[0]
      	temp1={the_cell_of_the_cluster:cluster_count}
      	cluster_id.update(temp1)
      	print(the_cell_of_the_cluster)
      	all_unclustered_cells.remove(the_cell_of_the_cluster)
      	collected_cells_of_the_cluster.remove(the_cell_of_the_cluster)
      	for the_neighbour_of_the_cell_of_the_cluster in active_cell_neighbours[the_cell_of_the_cluster]:
      	  active_cell_neighbours[the_neighbour_of_the_cell_of_the_cluster].remove(the_cell_of_the_cluster)
      	  collected_cells_of_the_cluster.append(the_neighbour_of_the_cell_of_the_cluster)
      cluster_count=cluster_count+1
    '''

    '''
    for i in all_neighbouring_cells_of_the_given_atoms:		# i = (ix,iy,iz) = cell position
      temp1={i:cluster_count}
      cluster_id.update(temp1)
      all_neighbouring_cells_of_the_given_atoms.remove(i)   # remove the entry from the list

      cluster1=[]
      #cells_that_are_already_accounted_in_cluster1=[]
      cluster1.append(i)
      for m in active_cell_neighbours[i]:
        active_cell_neighbours[m].remove(i)
      #finished_searching_cells_of_one_complete_cluster="no"
      #while finished_searching_cells_of_one_complete_cluster == "no":
      for j in cluster1:
          #self.cellX=j[0]
          #self.cellY=j[1]
          #self.cellZ=j[2]
          #neighbour_cells_of_cell_j=self.get_neighbour_cells()	# the output is a list = [(x1,y1,z1),(x2,y2,z2),(x3,y3,z3)]
          #cells_that_are_already_accounted_in_cluster1.append(j)
          neighbour_cells_of_cell_j=active_cell_neighbours[j]
          for k in neighbour_cells_of_cell_j:
            print("==> %s" %(str(k)))
            if k in all_neighbouring_cells_of_the_given_atoms:
              print("==>  ==> %s" %(str(k)))
              temp1={k:cluster_count}
              cluster_id.update(temp1)
              all_neighbouring_cells_of_the_given_atoms.remove(k)   # remove the entry from the list
              #if k not in cells_that_are_already_accounted_in_cluster1:
              for m in active_cell_neighbours[k]:
                active_cell_neighbours[m].remove(k)
              cluster1.append(k)
          print("j = %s" %(str(j)))
      #finished_searching_cells_of_one_complete_cluster="yes"
      print("finished cluster: %d" %(cluster_count))
      cluster_count =cluster_count+1
    '''
    for i in range(cluster_count-1):
      temp1={i+1:[]}
      self.clusters_id_to_position.update(temp1)

    for i in cluster_id.keys():
      temp1=(i[0]*self.LEC_X, i[1]*self.LEC_Y, i[2]*self.LEC_Z)		# i=cell id; i*LEC become cell position
      temp2={temp1:cluster_id[i]}
      self.clusters_position_to_id.update(temp2)
      temp2=self.clusters_id_to_position[cluster_id[i]]
      temp2.append(temp1)
      temp3={cluster_id[i]:temp2}
      self.clusters_id_to_position.update(temp3)

  def get_neighbour_cells(self):
      all_neighbour_cells=[]
      # get all surrounding 'vox_smooth' number of voxels 
      for ix in range(self.cellX-self.vox_smooth,self.cellX+self.vox_smooth+1):
        raw_ix=ix-self.cellX
        temp1=ix-self.cellX
        # Dealing with periodic boundary conditions
        if ix < 0:
          ix = self.MaxCells_X+ix
        if ix >= self.MaxCells_X:
          ix = ix - self.MaxCells_X

        for iy in range(self.cellY-self.vox_smooth,self.cellY+self.vox_smooth+1):
          raw_iy=iy-self.cellY
          temp2=iy-self.cellY
          # Dealing with periodic boundary conditions
          if iy < 0:
            iy = self.MaxCells_Y + iy
          if iy >= self.MaxCells_Y:
            iy = iy - self.MaxCells_Y

          for iz in range(self.cellZ-self.vox_smooth,self.cellZ+self.vox_smooth+1):
            raw_iz=iz-self.cellZ
            temp3=iz-self.cellZ
            # Dealing with periodic boundary conditions
            if iz < 0:
              iz = self.MaxCells_Z + iz
            if iz >= self.MaxCells_Z:
              iz = iz - self.MaxCells_Z
            #print(temp1*self.LEC_X,temp2*self.LEC_Y,temp3*self.LEC_Z,pow(temp1*self.LEC_X,2)+pow(temp2*self.LEC_Y,2)+pow(temp3*self.LEC_Z,2),self.LEC_X*self.LEC_Y*self.LEC*pow(self.rcut,3))
            #print(temp1,temp2,temp3)
            if (pow(temp1,2)+pow(temp2,2)+pow(temp3,2)) < pow(self.vox_smooth,2):
              all_neighbour_cells.append((ix,iy,iz))
      return all_neighbour_cells


class compute_clustered_channels:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * This is my own idea: -- Sudheer
  * Compute the channels by clustering (grouping) the given atoms
  *
  * usage: config1=ganisetti_tools.compute_clustered_channels(config,config_nnl,voxel_length,voxels_count_for_smoothing,given_atom_ids)
  * where config    = ganisetti_tools.get_atoms_info_from_lammps(01.dump)
  *       voxel_length = 2.5/5
  *       voxels_count_for_smoothing = 3
  *       atom_type_num2sym = {"1":"O", "2":"Si", "3":"Al", etc,.}
  *
  * output: config1.channels
  * 
  * HURRAY!, THIS IS WORKING AS I EXPECTED
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,config,config_nnl,rcut,vox_smooth,given_atoms):

    # some initializations
    all_unclustered_atoms={}
    cluster_id={}
    cluster_count=1
    clusters_id_to_position={}
    clusters_position_to_id={}
    config_nnl_of_given_atoms={}
    self.cluster_cells_id_to_position={}
    self.cluster_cells_position_to_id={}
    self.cluster_atoms_id_to_cluster_id={}

    for i in given_atoms:
      temp1={i:1}
      all_unclustered_atoms.update(temp1)
      temp1=[]
      for j in config_nnl[i]:
        if j in given_atoms:
          temp1.append(j)
      temp2={i:temp1}
      config_nnl_of_given_atoms.update(temp2)

    # main loop for identifying the clusters
    while len(all_unclustered_atoms) != 0:
      collected_atoms_of_the_cluster={}
      temp1={list(all_unclustered_atoms.keys())[0]:1}
      collected_atoms_of_the_cluster.update(temp1)
      while len(collected_atoms_of_the_cluster) != 0:
        current_atom_of_the_cluster=list(collected_atoms_of_the_cluster.keys())[0]
        temp1={current_atom_of_the_cluster:cluster_count}
        cluster_id.update(temp1)
        all_unclustered_atoms.pop(current_atom_of_the_cluster)
        collected_atoms_of_the_cluster.pop(current_atom_of_the_cluster)
        for neighbour_of_current_atom_of_the_cluster in config_nnl_of_given_atoms[current_atom_of_the_cluster]:
      	  temp1={neighbour_of_current_atom_of_the_cluster:1}
      	  collected_atoms_of_the_cluster.update(temp1)
        temp1=config_nnl[current_atom_of_the_cluster]
        temp1.append(current_atom_of_the_cluster)
        for pair in list(itertools.combinations(temp1,2)):
          try:
            config_nnl_of_given_atoms[pair[0]].remove(pair[1])
          except:
      	    pass
          try:
            config_nnl_of_given_atoms[pair[1]].remove(pair[0])
          except:
            pass
      cluster_count=cluster_count+1
    self.cluster_atoms_id_to_cluster_id=cluster_id
    
    # now descritize the sample into voxels so that it looks nicer for visualization
    self.vox_smooth=vox_smooth
    # Devide the box into voxels
    self.MaxCells_X=int(config.box_xx/rcut)
    self.MaxCells_Y=int(config.box_yy/rcut)
    self.MaxCells_Z=int(config.box_zz/rcut)

    # LEC= Length of Each Cell
    self.LEC_X=float(config.box_xx/self.MaxCells_X)
    self.LEC_Y=float(config.box_yy/self.MaxCells_Y)
    self.LEC_Z=float(config.box_zz/self.MaxCells_Z)

    # store atom positions locally
    atom_posx=config.posx
    atom_posy=config.posy
    atom_posz=config.posz
    # move the box to center and correspondingly atoms if the corner of the box is not at center
    for i in given_atoms:
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
      #assign the atoms into respective cell number  #assign the atoms into respective cell number; cell numbers starts with 0
      self.cellX=int(atom_posx[i]/self.LEC_X)
      self.cellY=int(atom_posy[i]/self.LEC_Y)
      self.cellZ=int(atom_posz[i]/self.LEC_Z)

      # compute neighboring cells of each atom
      neighbours1=self.get_neighbour_cells()
      for j in neighbours1:
        temp1={(j[0]*self.LEC_X,j[1]*self.LEC_Y,j[2]*self.LEC_Z):cluster_id[i]}
        self.cluster_cells_position_to_id.update(temp1)


  def get_neighbour_cells(self):
      all_neighbour_cells=[]
      # get all surrounding 'vox_smooth' number of voxels 
      for ix in range(self.cellX-self.vox_smooth,self.cellX+self.vox_smooth+1):
        raw_ix=ix-self.cellX
        temp1=ix-self.cellX
        # Dealing with periodic boundary conditions
        if ix < 0:
          ix = self.MaxCells_X+ix
        if ix >= self.MaxCells_X:
          ix = ix - self.MaxCells_X

        for iy in range(self.cellY-self.vox_smooth,self.cellY+self.vox_smooth+1):
          raw_iy=iy-self.cellY
          temp2=iy-self.cellY
          # Dealing with periodic boundary conditions
          if iy < 0:
            iy = self.MaxCells_Y + iy
          if iy >= self.MaxCells_Y:
            iy = iy - self.MaxCells_Y

          for iz in range(self.cellZ-self.vox_smooth,self.cellZ+self.vox_smooth+1):
            raw_iz=iz-self.cellZ
            temp3=iz-self.cellZ
            # Dealing with periodic boundary conditions
            if iz < 0:
              iz = self.MaxCells_Z + iz
            if iz >= self.MaxCells_Z:
              iz = iz - self.MaxCells_Z
            #print(temp1*self.LEC_X,temp2*self.LEC_Y,temp3*self.LEC_Z,pow(temp1*self.LEC_X,2)+pow(temp2*self.LEC_Y,2)+pow(temp3*self.LEC_Z,2),self.LEC_X*self.LEC_Y*self.LEC*pow(self.rcut,3))
            #print(temp1,temp2,temp3)
            if (pow(temp1,2)+pow(temp2,2)+pow(temp3,2)) < pow(self.vox_smooth,2):
              all_neighbour_cells.append((ix,iy,iz))
      return all_neighbour_cells


class read_parameter_file:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Since the complexity of the programme is increasing as the number of properties 
  * that are calculated are increasing, the parameters which are passing to the 
  * programme through the command line are now feed through a parameter file
  *
  * This class helps to read the parameter file
  * usage: 
  *   rpf = ganisetti_tools.read_parameter_file(parameter_file_name)
  *
  * output:
  *   rpf.input_file_format="imd" or "dump"
  *   rpf.imd_chkpt_file="chkpt_file.chkpt"
  *   rpf.lammps_dump_file="chkpt_file.dump"
  *   rpf.error_status = "yes" or "no"
  *   rpf.error_messages = ["error_message1", "error_message2", etc. ]
  *   rpf.all_arguments = ["python_script", "chkpt", "-O", "1", -Si", "2", "-SiO", "2.0", etc.]
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,parameter_file1):
    parameter_file=open(parameter_file1,'r')
    input_file_format   = "none"
    error_status        = False
    error_message       = ""
    read_bs_parameters  = False
    cur_file_exists     = False
    all_arguments       = []
    known_atoms         = []
    known_atom_pairs    = []
    error_messages      = []
    all_arguments.append(sys.argv[0])

    for each_line in parameter_file:
      data=each_line.split()
      if len(data) != 0 and data[0] != "#":
        # read known atoms
        if data[0] == "known_atoms":
          for i in range(2,len(data)):
            known_atoms.append(data[i])
          known_atom_pairs1=list( itertools.combinations(known_atoms, 2))
          for i,j in known_atom_pairs1:
            temp1=str(i)+str(j)
            known_atom_pairs.append(temp1)
            temp1=str(j)+str(i)
            known_atom_pairs.append(temp1)

        # read input file format
        if data[0] == "input_file_format":
          if len(data) == 3:
            input_file_format=data[2]
            self.input_file_format=input_file_format
          else:
            error_status=True
            error_messages.append("error in parameter file => 'input_file_format' is missing!")

        # read input file name
        if data[0] == "input":
          if input_file_format == "none":
            error_status= True
            error_messages.append("error in parameter file => please specify 'input_file_format' before 'input'")
          if len(data) == 3:
            all_arguments.append(data[2])
            self.input_file=data[2]
            #input_file=dump[2]+str(".")+str(input_file_format)
            if input_file_format == "imd":
              try:
                imd_chkpt_file=data[2]+str(".chkpt")
                temp1=open(imd_chkpt_file,'r')
                self.input_imd_chkpt_file=imd_chkpt_file
              except:
                error_status = True
                error_messages.append("error in parameter file => the input '.chkpt' file is missing!")
            elif input_file_format == "dump":
              try:
                lammps_dump_file=data[2]+str(".dump")
                temp1=open(lammps_dump_file,'r')
                self.input_lammps_dump_file=lammps_dump_file
              except:
                error_status = True
                error_messages.append("error in parameter file => the input '.dump' file is missing!")
          else:
            error_status = True
            error_messages.append("error in parameter file => 'input' is missing!")

        # read atom types and cutoff distances
        if data[0] in known_atoms:
          if len(data) == 3:
            all_arguments.append(str("-")+str(data[0]))
            all_arguments.append(data[2])
        if data[0] in known_atom_pairs:
          if len(data) == 3:
            all_arguments.append(str("-")+str(data[0]))
            all_arguments.append(data[2])

        # read keywords for computing bond statistics
        if data[0] == "compute_bs" and (data[2]).upper() == 'TRUE':
          read_bs_parameters = True
        if read_bs_parameters == True:
          if data[0] == "bs_cur_file":
            if len(data) == 3:
              cur_file=data[2]
              self.cur_file=cur_file
              try :
                if input_file_format == "imd":
                  cur_file1=cur_file+str(".chkpt")
                  open(cur_file1,'r')
                  cur_file_exists=True
                  self.cur_imd_chkpt_file=cur_file1
                elif input_file_format == "dump":
                  cur_file1=cur_file+str(".dump")
                  open(cur_file1,'r')
                  cur_file_exists=True
                  self.cur_lammps_dump_file=cur_file1
                else:
                  error_status=True
                  error_messages.append("error in parameter file => 'input_file_format' is missing")
              except:
                error_status=True
                error_messages.append("error in parameter file => the file %s does not exist!" %(cur_file1))
            else:
              error_status=True
              error_messages.append("error in parameter file => 'bs_cur_file' is missing!")
    if read_bs_parameters == True and cur_file_exists == False:
      error_status=True
      error_messages.append("error in parameter file => 'bs_cur_file' is missing!\n   bs_cur_file is needed inorder to compute bond statistics\n \
  if you do not want to compute bond statistics then use compute_bs = False in the parameter file")

    self.all_arguments = all_arguments
    self.error_status = error_status
    self.error_messages = error_messages


class generate_parameter_file_template:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Since the complexity of the programme is increasing as the number of properties 
  * that are calculated are increasing, the parameters which are passing to the 
  * programme through the command line are now feed through a parameter file
  *
  * This class helps to generate template of parameter file
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self):
    time_now              = datetime.datetime.now()
    output1=open("parameter_file_template.param", 'w')
    output1.write("# Author : Sudheer Ganisetti\n")
    output1.write("# Date   : %s \n" %(time_now))
    output1.write("# Info   : the parameter file can be used to pass the arguments to the python program upon calling 'ganisetti_tools' module\n")
    output1.write("           for computing a variety of properties of various glasses prepared with the molecular dynamics simulations\n\n")
    output1.write("# all known atoms\n")
    output1.write("known_atoms = O F Si Al P Na Ca Mg Sr Li K V\n\n")
    output1.write("# Known file formats = imd, dump\n")
    output1.write("input_file_format = dump\n\n")
    output1.write("# chkpt file\n")
    output1.write("input = chkpt_file\n\n")
    output1.write("# Atom Types\nO = 1\nSi = 2\nAl = 3\nP  = 6\nNa = 7\n\n")
    output1.write("# Cutoff Distances\n")
    output1.write("SiO = 2.0\nAlO = 2.40\nPO  = 2.00\nNaO = 3.15\n\n")
    output1.write("# bond_statistics\n# compute_bs = True or False\n")
    output1.write("# bs_ref_file = the 'input' file ; bs_cur_file = current file\n")
    output1.write("compute_bs = True\n")
    output1.write("bs_cur_file = current_chkpt_file\n\n")

class read_updated_command_line:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * Since the complexity of the programme is increasing as the number of properties 
  * that are calculated are increasing, the parameters which are passing to the 
  * programme through the command line are now feed through a parameter file
  * 
  * This class helps to read the command line and check the existance of parameter file
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,sys_argv):
  	error=False
  	if len(sys_argv) == 1 or sys_argv[1] == "-h" or sys_argv[1] == "--help":
  		error=True
  	elif sys_argv[1] == "-g":
  		generate_parameter_file_template()
  		print("The parameter file template is successfully generated")
  		sys.exit()
  	else:
  		try :
  			temp=open(sys_argv[1],'r')
  		except:
  			error=True
  	if error == True:
  		sudheer_banner()
  		CREDBG = '\33[31m' # \33[41m for red background
  		CREDBGEND = '\x1b[0m'
  		print("************************************** S. Ganisetti **************************************")
  		print("Error: usage is wrong")
  		print("\nThe correct usage is:")
  		print("%s   python3  %s  parameter_file  %s" % (CREDBG, str(sys_argv[0]), CREDBGEND))
  		print("\nIf you do not have the parameter_file then generate the template using the following command:")
  		print("%s   python3  %s  -g  %s " % (CREDBG, str(sys_argv[0]), CREDBGEND))
  		print("******************************************************************************************")
  		sys.exit()


def print_error_message(error_status,error_messages):
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * error_messages = ["error_message1", "error_message2", etc. ] 
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  if error_status == True:
    sudheer_banner()
    CREDBG = '\33[31m' # \33[41m for red background
    CREDBGEND = '\x1b[0m'
    print("************************************** S. Ganisetti **************************************\n")
    for i in error_messages:
      print("%s   %s  %s" % (CREDBG, str(i), CREDBGEND))
      print("\n******************************************************************************************")
      sys.exit()

class compute_bond_statistics:
  """
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * This class helps to compute the bond statistics 
  * 
  *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  """
  def __init__(self,given_cations_sym2num,ref_config,ref_config_nnl,cur_config_nnl):
    survived_bonds          = {}
    new_bonds               = {}
    broken_bonds            = {}
    bond_status             = {}
    switched_bonds__id2count= {}
    total_bonds             = 0
    total_new_bonds         = 0
    total_broken_bonds      = 0
    total_survived_bonds    = 0
    total_switched_bonds    = 0

    for i in ref_config.id:
      if ref_config.type[i] in given_cations_sym2num.values():
        temp1={i:[]}
        survived_bonds.update(temp1)
        new_bonds.update(temp1)
        broken_bonds.update(temp1)

        temp_broken_bonds  = list(set(ref_config_nnl.nnl[i]) - set(cur_config_nnl.nnl[i]))
        temp_new_bonds     = list(set(cur_config_nnl.nnl[i]) - set(ref_config_nnl.nnl[i]))
        temp_survived_bonds= list(set(cur_config_nnl.nnl[i]) & set(ref_config_nnl.nnl[i]))

        # for bond_status = survive
        if len(temp_broken_bonds) == 0 and len(temp_new_bonds) == 0:
          temp1={i:"survived_bond"}
          bond_status.update(temp1)
          for j in temp_survived_bonds:
            temp1=survived_bonds[i]
            temp1.append(j)
            temp2={i:temp1}
            survived_bonds.update(temp2)
            total_bonds=total_bonds+1
            total_survived_bonds=total_survived_bonds+1
          temp1={i:0}
          switched_bonds__id2count.update(temp1)
        # for bond_status = new
        elif ref_config_nnl.nnl_count[i] < cur_config_nnl.nnl_count[i]:
          temp1={i:"new_bond"}
          bond_status.update(temp1)
          for j in temp_broken_bonds:
            temp1=broken_bonds[i]
            temp1.append(j)
            temp2={i:temp1}
            broken_bonds.update(temp2)
            total_bonds=total_bonds+1
            total_broken_bonds=total_broken_bonds+1
          for j in temp_new_bonds:
            temp1=new_bonds[i]
            temp1.append(j)
            temp2={i:temp1}
            new_bonds.update(temp2)
            total_bonds=total_bonds+1
            total_new_bonds=total_new_bonds+1
          for j in temp_survived_bonds:
            temp1=survived_bonds[i]
            temp1.append(j)
            temp2={i:temp1}
            survived_bonds.update(temp2)
            total_bonds=total_bonds+1
            total_survived_bonds=total_survived_bonds+1
          total_switched_bonds = total_switched_bonds + len(temp_broken_bonds)
          temp1={i:len(temp_broken_bonds)}
          switched_bonds__id2count.update(temp1)
        # for bond_status = broken
        elif ref_config_nnl.nnl_count[i] > cur_config_nnl.nnl_count[i]:
          temp1={i:"broken_bond"}
          bond_status.update(temp1)
          for j in temp_broken_bonds:
            temp1=broken_bonds[i]
            temp1.append(j)
            temp2={i:temp1}
            broken_bonds.update(temp2)
            total_bonds=total_bonds+1
            total_broken_bonds=total_broken_bonds+1
          for j in temp_new_bonds:
            temp1=new_bonds[i]
            temp1.append(j)
            temp2={i:temp1}
            new_bonds.update(temp2)
            total_bonds=total_bonds+1
            total_new_bonds=total_new_bonds+1
          for j in temp_survived_bonds:
            temp1=survived_bonds[i]
            temp1.append(j)
            temp2={i:temp1}
            survived_bonds.update(temp2)
            total_bonds=total_bonds+1
            total_survived_bonds=total_survived_bonds+1
          total_switched_bonds = total_switched_bonds + len(temp_new_bonds)
          temp1={i:len(temp_new_bonds)}
          switched_bonds__id2count.update(temp1)
        # for bond_status = switched
        elif ref_config_nnl.nnl_count[i] == cur_config_nnl.nnl_count[i] and len(temp_broken_bonds) == len(temp_new_bonds):
          temp1={i:"switched_bond"}
          bond_status.update(temp1)
          for j in temp_broken_bonds:
            temp1=broken_bonds[i]
            temp1.append(j)
            temp2={i:temp1}
            broken_bonds.update(temp2)
            total_bonds=total_bonds+1
            total_broken_bonds=total_broken_bonds+1
          for j in temp_new_bonds:
            temp1=new_bonds[i]
            temp1.append(j)
            temp2={i:temp1}
            new_bonds.update(temp2)
            total_bonds=total_bonds+1
            total_new_bonds=total_new_bonds+1
          for j in temp_survived_bonds:
            temp1=survived_bonds[i]
            temp1.append(j)
            temp2={i:temp1}
            survived_bonds.update(temp2)
            total_bonds=total_bonds+1
            total_survived_bonds=total_survived_bonds+1
          total_switched_bonds = total_switched_bonds + len(temp_new_bonds)
          temp1={i:len(temp_new_bonds)}
          switched_bonds__id2count.update(temp1)
    self.total_bonds              = total_bonds
    self.total_broken_bonds       = total_broken_bonds
    self.total_new_bonds          = total_new_bonds
    self.total_switched_bonds     = total_switched_bonds
    self.total_survived_bonds     = total_survived_bonds
    self.switched_bonds__id2count = switched_bonds__id2count
    self.bond_status__id2text     = bond_status
    self.new_bonds__id2list       = new_bonds
    self.broken_bonds__id2list    = broken_bonds
    self.survived_bonds__id2list  = survived_bonds

