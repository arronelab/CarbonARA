#!/bin/bash
ScatterFile=Example/$1/$1Saxs.dat
fileLocs=Example/$1/
initialCoordsFile=Example/$1/coordinates.dat
noStructures=1
pairedPredictions=none
fixedsections=none
crystalSymmetry=none
withinMonomerHydroCover=none
betweenMonomerHydroCover=none
kmin=0.01;
kmax=0.1;
maxNoFitSteps=10

if [ ! -d Example/$1/$2 ]; then
  mkdir -p Example/$1/$2;
fi

for i in {1..1}
do
    predictStructure $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $betweenMonomerHydroCover $kmin $kmax $maxNoFitSteps Example/$1/$2/mol$i Example/$1/$2/scatter$i.dat Example/$1/mixtureFile.dat
done
 
