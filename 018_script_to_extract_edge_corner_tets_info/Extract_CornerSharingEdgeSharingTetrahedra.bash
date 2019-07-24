#!/bin/bash

echo "# 1)ChkptNo  2)Edge-Sharing_O  3)Edge-Sharing_O%  4)Corner-Sharing_O  5)Corner-Sharing_O%  6)Number_of_tetrahedra" > Edge_Corner_Sharing_Tetrahedra.data 
echo "# 7)Corner-sharing_tetrahedra  8)edge-sharing_tetrahedra  9)Corner-sharing_tetrahedra%  10)edge-sharing_tetrahedra%" >> Edge_Corner_Sharing_Tetrahedra.data
for i in {00000..00203}
do

cd $i
cd bonds
echo $i > temp0.txt
grep "Edge-Sharing  :" bond-prop.dat | awk '{print $3,$5}' >> temp1.txt
grep "Corner-Sharing:" bond-prop.dat | awk '{print $2,$4}' > temp2.txt
grep "Number of tetrahedra:" bond-prop.dat | awk '{print $4}' > temp3.txt
grep "Corner-sharing tetrahedra:" bond-prop.dat | awk '{print $3,$7}' > temp4.txt
grep "Which gives Td(CS):" bond-prop.dat | awk '{print $4,$8}' > temp5.txt

paste temp0.txt temp1.txt temp2.txt temp3.txt temp4.txt temp5.txt > temp6.txt
cat temp6.txt >> ../../Edge_Corner_Sharing_Tetrahedra.data

rm -rf temp0.txt temp1.txt temp2.txt temp3.txt temp4.txt temp5.txt temp6.txt

cd ..
cd ..

done
