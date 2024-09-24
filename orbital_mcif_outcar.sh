#!/bin/bash

no_of_species=$(awk ' {if ( NR == 6 )  print NF }' POSCAR)        # To find the type and total no.of species or elements
#echo "The no. of species is : "$no_of_species

# Creating the cif file from POSCAR

scl=$(awk 'NR==2 {print $1}' POSCAR)                      # To find the scaling parameter in POSCAR
#echo "Scale in POSCAR= $scl"
echo -e "\n_chemical_name_common\t '$(awk 'NR==1 {print $0}' POSCAR)' " > orbital_POSCAR.mcif
echo -e "_cell_length_a\t\t $(awk -v m=$scl 'NR==3 {print  m*sqrt($1^2+$2^2+$3^2)}' POSCAR|bc -l)">>orbital_POSCAR.mcif
echo -e "_cell_length_b\t\t $(awk -v m=$scl 'NR==4 {print  m*sqrt($1^2+$2^2+$3^2)}' POSCAR|bc -l)">>orbital_POSCAR.mcif
echo -e "_cell_length_c\t\t $(awk -v m=$scl 'NR==5 {print  m*sqrt($1^2+$2^2+$3^2)}' POSCAR|bc -l)">>orbital_POSCAR.mcif

a1=$(awk 'NR==3{print $1}' POSCAR)
a2=$(awk 'NR==3{print $2}' POSCAR)
a3=$(awk 'NR==3{print $3}' POSCAR)
b1=$(awk 'NR==4{print $1}' POSCAR)
b2=$(awk 'NR==4{print $2}' POSCAR)
b3=$(awk 'NR==4{print $3}' POSCAR)
c1=$(awk 'NR==5{print $1}' POSCAR)
c2=$(awk 'NR==5{print $2}' POSCAR)
c3=$(awk 'NR==5{print $3}' POSCAR)

#echo $a1 $a2 $a3 $b1 $b2 $b3 $c1 $c2 $c3

i_alpha=$(echo "($b1*$c1+$b2*$c2+$b3*$c3)/(sqrt(($b1)^2+($b2)^2+($b3)^2)*sqrt(($c1)^2+($c2)^2+($c3)^2))" | bc -l)

if [ "${i_alpha}" == 0  ]
          then
          alpha=90.0
      else
          alpha=$(echo "scale=5;180+a(sqrt(1-($i_alpha)^2)/$i_alpha)*180/3.14159265359" | bc -l)
fi

i_beta=$(echo "($a1*$c1+$a2*$c2+$a3*$c3)/(sqrt(($a1)^2+($a2)^2+($a3)^2)*sqrt(($c1)^2+($c2)^2+($c3)^2))" | bc -l)
if [ "${i_beta}" == 0  ]
          then
          beta=90.0
      else
          beta=$(echo "scale=5;180+a(sqrt(1-($i_beta)^2)/$i_beta)*180/3.14159265359" | bc -l)
fi

i_gamma=$(echo "($b1*$a1+$b2*$a2+$b3*$a3)/(sqrt(($b1)^2+($b2)^2+($b3)^2)*sqrt(($a1)^2+($a2)^2+($a3)^2))" | bc -l)
if [ "${i_gamma}" == 0  ]
          then
          gamma=90.0
      else
          gamma=$(echo "scale=5;180+a(sqrt(1-($i_gamma)^2)/$i_gamma)*180/3.14159265359" | bc -l)
fi

echo -e "_cell_angle_alpha\t\t $alpha">> orbital_POSCAR.mcif
echo -e "_cell_angle_beta\t\t $beta">> orbital_POSCAR.mcif
echo -e "_cell_angle_gamma\t\t $gamma">> orbital_POSCAR.mcif

echo -e "\nloop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z">>orbital_POSCAR.mcif

NION=0
for (( i=1; i<=$no_of_species; i++ ))
do
      spec[$i]=$(awk -v j="$i" '{ if ( NR == 6 ) print $j }' POSCAR)
      tot_atm_spec[$i]=$(awk -v j="$i" '{ if ( NR == 7 ) print $j }' POSCAR)

           for (( w=1; w<=${tot_atm_spec[$i]}; w++  ))
           do
                   echo -e "${spec[$i]}$w\t $(awk -v x=$(($NION+$w+8)) 'NR==x {print $0}' POSCAR)">>orbital_POSCAR.mcif

           done

      st_no[$i]=$(($NION+1))
      echo "Species $i : " ${spec[$i]} ${tot_atm_spec[$i]} " Starts from atom no. " ${st_no[$i]}
      NION=$(($NION+${tot_atm_spec[$i]}))

done

#echo "NION= " $NION
offst1=$((($NION+8)*3))            # offset or the region to look for in OUTCAR
#offst2=$(($offst1-$NION-9))
#echo "Offset1  " $offst1
#echo "Offset2  " $offst2



# To put mcif information in orbital_POSCAR.mcif

echo -e -n '\nloop_
_atom_site_moment.label
_atom_site_moment.crystalaxis_x
_atom_site_moment.crystalaxis_y
_atom_site_moment.crystalaxis_z'>>  orbital_POSCAR.mcif


tot_mag_ion=0        # initializing variables
inpt="Q"

while [ "$inpt" != "0" ]                  #finding which species is magnetic
do
	tot_mag_ion=$(($tot_mag_ion+1))
	read -p "Enter the magnetic species in correct order (enter 0 after you are done): " mag_ion[$tot_mag_ion]
        inpt=${mag_ion[$tot_mag_ion]}
#        echo "You entered: $inpt "
	for (( j=1; j<=$no_of_species; j++ ))                                        #finding species no. of magnetic ions
	do
             if [ "${spec[$j]}" == "${mag_ion[$tot_mag_ion]}"  ]
	     then
	     mag_ion_no[$tot_mag_ion]=$j
#	     echo "Species no. :  ${mag_ion_no[$tot_mag_ion]} "
	     fi
	done
done

tot_mag_ion=$(($tot_mag_ion-1))                     #to correct for the 0 entered by user


grep -A$offst1 "orbital moment (x" OUTCAR  >  temp.txt   # to create a file with required magnetization of last iteration

last_col=$(grep -w -m1 "1 " temp.txt |awk ' END{print NF}')  # to ascertain the no. of orbitals i.e. till d or f
#echo -e "\n Last column " $last_col

grep -A4 "from SAXIS" OUTCAR > temp2.txt

for (( p=1; p<=$tot_mag_ion; p++ ))
do
	for (( q=1; q<=${tot_atm_spec[${mag_ion_no[$p]}]} ; q++ ))
	do
	atm_no=$((${st_no[${mag_ion_no[$p]}]}+$q-1))       # current atom no and also correct for extra 1 from q
       echo -n -e "\n"${mag_ion[$p]}""$q " Current atom no. "$atm_no " "
	echo -n -e "\n"${mag_ion[$p]}""$q>> orbital_POSCAR.mcif



		for (( r=1; r<=3; r++  ))
	        do

		m_axis[$r]=$(grep -w "$atm_no " temp.txt |awk -v j=$r -v k=$last_col 'NR==j {print $k}')
          #    echo -n -e "\t"${m_axis[$r]}

	        done


	        for (( r=1; r<=3; r++  ))
		do
		m=0
		v=1
        		for (( u=1; u<=3; u++ ))
        		do

        		matrx[$r,$u]=$(awk -v a=$(($r+2)) -v b=$v 'NR==a {print $b}' temp2.txt|bc -l)
                m=$(echo "$m+$(echo "${m_axis[$u]}*${matrx[$r,$u]}"|bc -l)"|bc -l)
	            v=$(($v+2))

        		done

	        echo -n -e "\t"
                printf '%2.3f' $m
#        	echo -n -e "\t"$m>>orbital_POSCAR.mcif
                echo -n -e "\t" >> orbital_POSCAR.mcif
                printf '%2.3f' $m >> orbital_POSCAR.mcif
		done
        echo -n -e "\t |L| : $(tail -1 orbital_POSCAR.mcif| awk '{print sqrt($2^2+$3^2+$4^2)}')"
        echo -n -e "\t Theta L in x-z Plane: $(tail -1 orbital_POSCAR.mcif| awk '{print atan2 ($4,$2) * 180 / 3.14159265359}'|bc -l)"
        done


done

echo -e "\n"
rm temp.txt --verbose
rm temp2.txt --verbose
echo -e "\n The Total Moment therefore is: \n"
paste POSCAR.mcif orbital_POSCAR.mcif | tail -$(($q-1)) | awk '{print $1,$2+$6,$3+$7,$4+$8}' | awk '{print $0, "   |J| ", sqrt($2^2+$3^2+$4^2) , "    Theta J in x-z Plane ", atan2($4,$2)*180/3.14159265359}'
echo -e "\n----End of program---- "
