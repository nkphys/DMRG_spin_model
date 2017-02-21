#!/bin/bash -x
length=$1
file_2p="sites_SzSz_$length.txt"
file_1p="sites_Sz_$length.txt"
file_4p="sites_4point_$length.txt"
rm -rf $file_2p
rm -rf $file_1p
rm -rf $file_4p

#2_point
j=0;
let l_1=length-1;
while [ $j -lt $l_1 ];
do
let i=j+1;
while [ $i -lt $length ];
do 
printf "$j $i\n" >> $file_2p; 
let i=i+1;
done
let j=j+1;
done

#4_point
j=0;
let l_3=length-3;
let l_2=length-2;
while [ $j -lt $l_3 ];
do
let j_=j+1;
let i=j+2;
while [ $i -lt $l_1 ];
do
let i_=i+1;
printf "$j $j_ $i $i_\n" >> $file_4p;
let i=i+1;
done
let j=j+1;
done

#1_point
j=0;
while [ $j -lt $length ];
do
printf "$j\n" >> $file_1p;
let j=j+1;
done
