How to use MATLAB Version of GADIHOM
DOS REIS F. 
04.04.2022
Introduction
This is the MATLAB version of GADIHOM, a Genetic Algorithm for Inverse HOMogenization1. 
The code is composed with modules :
1. GADIHOM.m : main function module of the genetic algorithm
2. Homogenization.m : asymptotic homogenization module
3. Mesh.m : mesh module 
4. Compliance.m : Compute the compliance tensor from mechanics moduli
5. Save_matrix.m : save a matrix in csv format
6. Fcost.m : cost function for the genetic algorithm
7. mechanic_moduli.m : extract mechanic moduli from homogenized compliance tensor 
How to use 
*  modify the data values in a ExampleUsei.m code to call the main function GADIHOM(…)

* Run 
* Results are stored in various csv file :
o Tb.csv : width’s beam values of best lattice found. See the appendix in paper to know the store order used
o Mechanic_homogenized.csv : homogenized mechanic moduli values for best lattice. Contain : 
o Other files containing the topology of the lattice
1 Not all the features are included. A translated fastest and enhanced C code was available.
---------------

------------------------------------------------------------

---------------

------------------------------------------------------------

