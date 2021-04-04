#!/usr/bin/python
"""
@Author : S. Ganisetti
@Date   : Sa 6. Aug 13:41:41 CEST 2016
@modified: Fr 2. Apr 17:42:34 CEST 2021
The code is to convert LAMMPS dump file to xyz file
"""
import numpy as np
import math
import sys

def call_error_message():
     print("************* S. Ganisetti *************")
     print("Error: usage is wrong")
     print("./this_program lammps.dump O 1 Si 2 Al 3 etc.")
     print("This programm coonverts lammps dump file to xyz")
     print("****************************************")


""" **************** Main Function **************** """
if __name__=="__main__":
  chemical_elements = []
  chemical_elements += ['H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne'] # 1 to 10
  chemical_elements += ['Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca'] # 11 to 20
  chemical_elements += ['Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'] # 21 to 30
  chemical_elements += ['Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr'] # 31 to 40
  chemical_elements += ['Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'] # 41 to 50
  chemical_elements += ['Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd'] # 51 to 60
  chemical_elements += ['Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb'] # 61 to 70
  chemical_elements += ['Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'] # 71 to 80
  chemical_elements += ['Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th'] # 81 to 90
  chemical_elements += ['Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm'] # 91 to 100
  chemical_elements += ['Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds'] # 101 to 110
  chemical_elements += ['Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts',  'Og']            # 111 to 118

  dump_file=str(sys.argv[1])
  try:
    dump1=np.loadtxt(dump_file,dtype='int,int,float,float,float,float,float,float,float',skiprows=9)
  except:
    call_error_message()
    print("Please check the existance of file %s in the given path" (sys.argv[1]))
    exit()

  atom_types={}
  for i in range(2,len(sys.argv)):
    atom_string = sys.argv[i]
    if i%2 == 0:
        temp1 = str(atom_string)
        if temp1 not in chemical_elements:
            call_error_message()
            print("%s is probably not a chemical element, please check it!" (temp1))
            exit()
    else :
        temp2={}
        try:
            temp2={int(atom_string):temp1}
            atom_types.update(temp2)
        except TypeError:
            call_error_message()
            print("%s is not a valid atom type for %s, please check!" (atom_string,temp1))
            exit()
        except:
            call_error_message()
            print("The atom number %d is given more than once, please check on it!" (int(atom_string)))
            exit()

  xyz_file=str(dump_file[:-4])+str('xyz')
  output1=open(xyz_file,'w')

  # reading lammps input file
  data2={}
  atom_ids = []
  for i in range(len(dump1)):
    atom_id=dump1[i][0]
    try:
        atom_type=atom_types[dump1[i][1]]
    except KeyError:
        call_error_message()
        print("atom type " +str(dump1[i][1])+ " is not found")
        exit()
    posx=dump1[i][3]
    posy=dump1[i][4]
    posz=dump1[i][5]

    temp1={atom_id:[atom_type, posx, posy, posz]}
    data2.update(temp1)
    atom_ids.append(atom_id)

  output1.write(str(len(atom_ids))+"\n")
  output1.write("# ")
  for i in sys.argv:
    output1.write(str(i) + " ")
  output1.write("\n")


  atom_ids = np.array(atom_ids)
  atom_ids.sort()
  for i in atom_ids:
    for j in data2[i]:
        output1.write(str(j)+" ")
    output1.write("\n")
