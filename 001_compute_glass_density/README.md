This piece of code is used to compute the density of a glass with given atom positions and atom types in a lammps dump file format

usage ./compute_glass_density.py  lammps.dump -O 1 -Si 2 -Al 3 -Ca 4 -Mg 5 ...

get help by just typing ./compute_glass_density.py

At the moment the data base contains a combination of following chemical compositions:
SiO2 + Al2O3 + CaO + MgO + SrO + Na2O + P2O5

