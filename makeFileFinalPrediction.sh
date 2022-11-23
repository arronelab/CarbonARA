#/bin/sh
rm src/*.o

g++ -c -O3 -std=gnu++14 -o src/point.o  src/point.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/polyHelix.o src/polyHelix.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/randomMolGen.o src/randomMolGen.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/ktlMoleculeRandom.o src/ktlMoleculeRandom.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/experimentalData.o src/experimentalData.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/hydrationShellRandom.o src/hydrationShellRandom.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/writhe.o src/writhe.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/mainPredictionFinal.o src/writhe.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/moleculeFitAndState.o src/moleculeFitAndState.cpp -fopenmp
g++ -c -O3 -std=gnu++14 -o src/mainPredictionFinal.o src/mainPredictionFinal.cpp -fopenmp

g++ -c -O3 -std=gnu++14 -o src/Flexible_generator.o src/Flexible_generator.cpp -fopenmp

g++ -O3 -std=gnu++14 -o predictStructure src/point.o src/polyHelix.o src/randomMolGen.o src/ktlMoleculeRandom.o src/experimentalData.o src/hydrationShellRandom.o src/writhe.o src/moleculeFitAndState.o src/mainPredictionFinal.o -fopenmp

g++ -O3 -std=gnu++14 -o generate_structure src/point.o src/polyHelix.o src/randomMolGen.o src/ktlMoleculeRandom.o src/experimentalData.o src/hydrationShellRandom.o src/writhe.o src/moleculeFitAndState.o src/Flexible_generator.o -fopenmp
