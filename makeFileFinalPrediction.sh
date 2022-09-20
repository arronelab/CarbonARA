#/bin/sh
rm finalSrc/*.o

g++ -c -O3 -std=gnu++14 -o src/point.o  src/point.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/optimizationPoint.o src/optimizationPoint.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/polyHelix.o src/polyHelix.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/randomMolGen.o src/randomMolGen.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/ktlMoleculeRandom.o src/ktlMoleculeRandom.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/experimentalData.o src/experimentalData.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/hydrationShellRandom.o src/hydrationShellRandom.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/binaryFind.o src/binaryFind.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/getSections.o src/getSections.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/mutualWind2.o src/mutualWind2.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/localWrithe.o src/localWrithe.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/moleculeFitAndState.o src/moleculeFitAndState.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/mainPredictionFinal.o src/mainPredictionFinal.cpp -fopenmp
g++ -O3 -std=gnu++14 -o predictStructure src/point.o src/optimizationPoint.o src/polyHelix.o src/randomMolGen.o src/ktlMoleculeRandom.o src/experimentalData.o src/hydrationShellRandom.o src/binaryFind.o src/mutualWind2.o src/localWrithe.o src/moleculeFitAndState.o src/mainPredictionFinal.o -fopenmp