# Sex_chromosome_degeneration

The different programs simulate N diploid individuals with discrete generations. Males are XY and females are XX (see article for more detailed description of the programs).

master: standard model, with one pair of male-specific and female-specific trans acting modifiers for each gene on the sex chromosome, and asymmetric function for stabilizing selection on level of gene expression

symmetric: same as standard model, with symmetric function for stabilizing selection on level of gene expression

single_trans: only one pair of male-specific and female-specific trans acting modifiers for all genes on the sex chromosome

The programs use the Mersenne Twister random number generator, and the file MersenneTwister.h (available from http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/MersenneTwister.h) must be included in the program before compilation.

The file « sexchr.h » provides function prototypes and global variables, in particular the structure « ind » that represents an individual. The main function (in "main.cpp") calls functions (defined in "files.cpp") that read parameters from an input file ("parameters.txt"), write parameters in the file "results.txt" and then call the function "recursion" that performs the simulation. The file "parameters.txt" provides parameter values in the order indicated on the first line of the file; each line corresponding to a parameter set must start with * (several sets of parameters can be added, the program will run one after the other). The file « fitness.cpp" contains the functions needed to compute the fitness of males and females, and « rec.cpp » the functions to construct recombinant chromosomes. The file « record.cpp » contains functions to measure different statistics from the population.

The parameter « Rep » specifies the number of replicate simulations to be performed for each parameter set. For each replicate, an output file is produced with the parameter values in the file name, giving the different statistics measured at different generations (averaged over all genes if output=1, or for each gene if output=0). A file called « summary…txt » is also produced, showing the averages over replicates (this file is updated after each replicate, and is filled only when output=1).

Using gcc, the programs can by compiled using the command
g++ -O3 *.cpp
