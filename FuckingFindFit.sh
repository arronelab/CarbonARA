#!/bin/bash

for i in {1..3}

do
   echo " Run number : $i "

   ./predictStructure $1/human_SMARCAL1Saxs.dat $1/ $1/coordinates1.dat none $1/varyingSectionSecondary1.dat 1 none none 0.02 0.25 3 $1/fittys/fitmolecule$i $1/fittys/scatter$i.dat $1/mixtureFile.dat 1
done

