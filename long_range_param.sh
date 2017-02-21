#!/bin/bash -x
Length=$1
J1_=$2
J2_=$3
alpha_=$4
mag_=$5


alpha_by4=$(bc <<< "scale=10;$alpha_/4")
J1_by2=$(bc <<< "scale=10;$J1_/2")
J2_by2=$(bc <<< "scale=10;$J2_/2")

J1_zz=$J1_
J2_zz=$J2_
J1_pm=$(bc <<< "scale=10;$J1_by2+$alpha_by4")
J1_mp=$(bc <<< "scale=10;$J1_by2+$alpha_by4")
J2_pm=$J2_by2
J2_mp=$J2_by2
J1_pp=$alpha_by4
J1_mm=$alpha_by4
J2_pp="0.0"
J2_mm="0.0"


file_Jzz=Jzz_longrange_$Length.txt
file_Jpm=Jpm_longrange_$Length.txt
file_Jmp=Jmp_longrange_$Length.txt
file_Jpp=Jpp_longrange_$Length.txt
file_Jmm=Jmm_longrange_$Length.txt
file_H=H_file_$Length.txt

printf "Input params beaing created for L=$Length system \n"
rm $file_Jzz;
rm $file_Jpm;
rm $file_Jmp;
rm $file_Jpp;
rm $file_Jmm;
rm $file_H;

H_mag="$mag_";
J_zz=('0.0');
J_pm=('0.0');
J_mp=('0.0');
J_pp=('0.0');
J_mm=('0.0');
H=("$H_mag");
counter=1
while [ $counter -lt $Length ];
do
J_zz=("${J_zz[@]}" "0.0")
J_pm=("${J_pm[@]}" "0.0")
J_mp=("${J_mp[@]}" "0.0")
J_pp=("${J_pp[@]}" "0.0")
J_mm=("${J_mm[@]}" "0.0")
H=("${H[@]}" "$H_mag")

let counter=counter+1
done
#printf ${#J_zz[@]}
#--------------RIGHT HERE THE CONNECTIONS VALUES like A[d];where d=|i-j|------
J_zz[1]="$J1_zz"
J_zz[2]="$J2_zz"
J_pm[1]="$J1_pm"
J_pm[2]="$J2_pm"
J_mp[1]="$J1_mp"
J_mp[2]="$J2_mp"
J_pp[1]="$J1_pp"
J_pp[2]="$J2_pp"
J_mm[1]="$J1_mm"
J_mm[2]="$J2_mm"
#----------------CONNECTIONS VALUES WRITTEN-----------------------------------

#Jzz
site1=0
while [ $site1 -lt $Length ];
do
site2=0
while [ $site2 -lt $Length ];
do
if [ $site1 -lt $site2 ]; then
dis=$(echo "$site2 -$site1" | bc)
fi
if [ $site2 -lt $site1 ]; then
dis=$(echo "$site1 -$site2" | bc)
fi
if [ "$site2" -eq "$site1" ]; then
dis=0
fi
printf '%.10g ' ${J_zz[$dis]}  >> $file_Jzz
printf '%.10g ' ${J_pm[$dis]} >> $file_Jpm
printf '%.10g ' ${J_mp[$dis]} >> $file_Jmp
printf '%.10g ' ${J_pp[$dis]} >> $file_Jpp
printf '%.10g ' ${J_mm[$dis]} >> $file_Jmm
#printf "$dis "
let site2=site2+1
done
printf "\n" >> $file_Jzz
printf "\n" >> $file_Jpm
printf "\n" >> $file_Jmp
printf "\n" >> $file_Jpp
printf "\n" >> $file_Jmm
printf '%.10g ' ${H[$site1]} >> $file_H
printf "\n" >> $file_H
let site1=site1+1
done 

