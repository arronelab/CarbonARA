#!/bin/bash
# $1 is the molecule name, $2 the fit name
ScatterFile=$1/$1Saxs.dat
fileLocs=$1/
initialCoordsFile=$1/coordinates.dat
noStructures=1
pairedPredictions=none
fixedsections=$1/varyingSectionSecondary.dat 
crystalSymmetry=none
withinMonomerHydroCover=none
betweenMonomerHydroCover=none
kmin=0.02;
kmax=0.025;
maxNoFitSteps=10

mkdir $1/$2


for i in {1..5}
do
    ./predictStructure $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $betweenMonomerHydroCover $kmin $kmax $maxNoFitSteps $1/$2/mol$i $1/$2/scatter$i.dat $1/mixtureFile.dat
done
 
