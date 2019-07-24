# Author: Sudheer Ganisetti

imd_util is modified to get the projection of RDF on each of the axes

usage:
compile: make pair
run    : pair -e 5.0 -n 500 -p parameter.param

cat parameter.param
coordname       input_imd.chkpt
outfiles        output.pair

box_from_header 1
pbc_dirs        1 1 1
ntypes          2
nbl_margin      0.5
nbl_size        1.2

