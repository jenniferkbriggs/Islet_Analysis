#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=16:00:00
#SBATCH --partition=amilan
#SBATCH --output=out_%j.out

module purge
module load gcc/14.2.0
module load boost/1.86.0

#Execute

rm potential.txt sodium.txt potassium.txt calcium.txt caer.txt atp.txt adp.txt IRP.txt PP.txt DP.txt FIP.txt RIP.txt Cap.txt time.txt noise.txt
rm spheres Islet RandomGenerator Beta

#compile sphere and generate sphere input file
g++ neighbor.C spheres.C box.C sphere.C event.C heap.C read_input.C -o spheres 
./spheres sphereInput

#create islet
g++ -I /projects/$USER/boost/boost_1_58_0/ -fopenmp GenerateSphere.cpp -o Islet
./Islet isletInput

#compile and general random numbers
g++ -std=c++0x -I /projects/$USER/boost/boost_1_58_0/ -fopenmp RandomVars.cpp -o RandomGenerator
./RandomGenerator

#run main code for solving differential equations
g++ -I /projects/$USER/boost/boost_1_58_0/ -fopenmp MainFile.cpp -o Beta
./Beta 0B 
