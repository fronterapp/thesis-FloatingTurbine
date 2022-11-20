#!/bin/bash

## FILE TO READ
file=log.waveDyMFoam

## CREATE HEADER, WRITE IT ON 'RB_DATA.dat'
printf 'Time \t F1x \t F1y \t F1z \t F2x \t F2y \t F2z \t F3x \t F3y \t F3z \n' > MOORING_DATA.dat

## SEARCH FOR THE DESIRED KEYWORDS
## SAVE THE DATA ON 'my_output'
my_output="$(awk '
	/Time/{TIME=$3; next}
	/Restraint catenaryLine1/{F1x=$6; F1y=$7; F1z=$8; next} 
	/Restraint catenaryLine2/{F2x=$6; F2y=$7; F2z=$8; next} 
	/Restraint catenaryLine3/{F3x=$6; F3y=$7; F3z=$8; 
	printf("%8.6f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	TIME,F1x,F1y,F1z,F2x,F2y,F2z,F3x,F3y,F3z);}' $file)"

## PRINT 'my_output' 
echo "$my_output" >> temp.dat

## REMOVE PARENTHESES FROM 'my_output' 
sed -e 's/(//g' -e 's/)//g' temp.dat >> MOORING_DATA.dat
rm temp.dat