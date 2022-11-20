#!/bin/bash

## FILE TO READ
file=log.waveDyMFoam

## CREATE HEADER, WRITE IT ON 'RB_DATA.dat'
printf 'Time \t CORx \t CORy \t CORz \t COMx \t COMy \t COMz \t Rxx \t Rxy \t Rxz \t Ryx \t Ryy \t Ryz \t Rzx \t Rzy \t Rzz \t vx \t vy \t vz \t ax \t ay \t az \t wx \t wy \t wz \t Lx \t Ly \t Lz \t Tx \t Ty \t Tz \n' > RB_DATA.dat

## SEARCH FOR THE DESIRED KEYWORDS
## SAVE THE DATA ON 'my_output'
my_output="$(awk '
	/Time/{TIME=$3; next}
	/Centre of rotation/{CORx=$4; CORy=$5; CORz=$6; next} 
	/Centre of mass/{COMx=$4; COMy=$5; COMz=$6; next}
	/Orientation/{Rxx=$2; Rxy=$3; Rxz=$4; Ryx=$5; Ryy=$6; Ryz=$7; Rzx=$8; Rzy=$9; Rzz=$10; next}
	/Linear velocity/{vx=$3; vy=$4; vz=$5; next}
	/Linear acceleration/{ax=$3; ay=$4; az=$5; next}
	/Angular velocity/{wx=$3; wy=$4; wz=$5; next}
	/Angular momentum/{Lx=$3; Ly=$4; Lz=$5; next}
	/Torque/{Tx=$2; Ty=$3; Tz=$4;

	printf("%8.6f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	TIME,CORx,CORy,CORz,COMx,COMy,COMz,Rxx,Rxy,Rxz,Ryx,Ryy,Ryz,Rzx,Rzy,Rzz,vx,vy,vz,ax,ay,az,wx,wy,wz,Lx,Ly,Lz,Tx,Ty,Tz);}' $file)"

## PRINT 'my_output' 
echo "$my_output" >> temp.dat

## REMOVE PARENTHESES FROM 'my_output' 
sed -e 's/(//g' -e 's/)//g' temp.dat >> RB_DATA.dat
rm temp.dat