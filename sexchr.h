// Header file: definitions of global variables, function prototypes

#ifndef SEXCHR_H
#define SEXCHR_H

#include <vector>
#include <algorithm>
#include <iostream>
#include "MersenneTwister.h"
using namespace std;

// Global variables:

#define inputFile "parameters.txt"     // names of input
#define outputFile "results.txt"		// and output files

// "ind": represents an individual

struct ind
{
    double * gene; // allelic values at loci affecting to fitness
    double * cis; // values at cis loci controlling regulation
    double * transm; // trans values for each trans male locus
    double * transf; // trans values for each trans female locus
};

// Function prototypes:

void openInFile();
void openOutFile();
void writeParameters(int Nv, double sigtrans, double sigcis, int nbSv, int NbGenv, int NbPrelimv, int pasv, double s, double s_max, double I, double U_g, double U_c, double U_tm,double U_tf, double Rgv, double Rgc, int outputv);
bool readFile(int &Nr, double &sigtransr, double &sigcisr, int &nbSr, int &NbGenr, int &NbPrelimr, int &pasr, double &sr, double &s_maxr, double &Ir, double &U_gr, double &U_cr, double &U_tm,double &U_tf, double &Rgr, double &Rcr, int &Rep, int &outputr);
void recursion(int Nv, double sigtrans, double sigcis, int nbSv, int NbGenv, int NbPrelimv, int pasv, double s, double s_max, double I, double U_g, double U_c, double U_tm, double U_tf, double Rg, double Rc, int Rep, int output, double** allAverages);
double gammln(const double xx);
double poisdev(const double xm);
double gasdev();
double binldev(const double pp, const int n);
void rec(ind &offspring, ind &father, ind &mother, double Rg, double Rc, double ML, int nS, int off_sex);
void recInit(ind &offspring, ind &father, ind &mother, double Rg, double Rc, double ML, int nS, int off_sex);
double Wmale(ind &parent, double h0, double s_max, double Q0, double I, int nbS);
double Wfemale(ind &parent, double h0, double s_max, double Q0, double I, int nbS);
void record_output(ind * pop, double * Wtot, double ** measures, double * popAve, int nbSv, int Nmales, int Nv, double h0);
void record_averages(ind * pop, double * Wtot, double * popAve, int nbSv, int Nmales, int Nv, double expdom, double smax);

#endif
