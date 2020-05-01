// Functions to open input and output files,
// read parameter values from input file and
// write them in output file.

#include "sexchr.h"
#include <iostream>
#include <fstream>
using namespace std;

extern FILE * inFile;
extern FILE * outFile;

// opens input file:

void openInFile()
{
	inFile = fopen(inputFile,"r");
}


// opens output file:

void openOutFile()
{
	outFile = fopen(outputFile,"a");
}


// reads parameter values from input file,
// returns 1 if end of input file, else returns 0:

bool readFile(int &Nr, double &sigtransr, double &sigcisr, int &nbSr, int &NbGenr, int &NbPrelimr, int &pasr, double &sr, double &s_maxr, double &Ir, double &U_gr, double &U_cr, double &U_tm,double &U_tf, double &Rgr, double &Rcr, int &Rep, int &outputr)
{
	int x;
	bool term;
	do {x = fgetc(inFile);} while (!((x == '*') || (x == EOF)));
		// each parameter set must start with *
	if (x == EOF)
	{
		term = true;
	}
	else
	{
		fscanf(inFile,"%d ",&Nr);
		fscanf(inFile,"%lf ",&sigtransr);
		fscanf(inFile,"%lf ",&sigcisr);
		fscanf(inFile,"%d ",&nbSr);
		fscanf(inFile,"%d ",&NbGenr);
        fscanf(inFile,"%d ",&NbPrelimr);
		fscanf(inFile,"%d ",&pasr);
        fscanf(inFile,"%lf ",&sr);
        fscanf(inFile,"%lf ",&s_maxr);
        fscanf(inFile,"%lf ",&Ir);
        fscanf(inFile,"%lf ",&U_gr);
        fscanf(inFile,"%lf ",&U_cr);
	    fscanf(inFile,"%lf ",&U_tm);
	    fscanf(inFile,"%lf ",&U_tf);
	    fscanf(inFile,"%lf ",&Rgr);
        fscanf(inFile,"%lf ",&Rcr);
        fscanf(inFile,"%d ",&Rep);
        fscanf(inFile,"%d ",&outputr);

		term = false;
	}
	return term;
}


// writes parameter values in output file:

void writeParameters(int Nv, double sigtrans, double sigcis, int nbSv, int NbGenv, int NbPrelimv, int pasv, double s, double s_max, double I, double U_g, double U_c, double U_tm,double U_tf, double Rgv, double Rcv, int outputv)
{
	fprintf(outFile,"\n_________________________________________\n");
	fprintf(outFile,"\nN = %d", Nv);
	fprintf(outFile,", sigtrans = %g", sigtrans);
	fprintf(outFile,", sigcis = %g", sigcis);
	fprintf(outFile,", nbS = %d", nbSv);
	fprintf(outFile,", generations = %d", NbGenv);
    fprintf(outFile,", NbPrelim = %d", NbPrelimv);
	fprintf(outFile,", pas = %d", pasv);
    fprintf(outFile,"\ns = %g", s);
    fprintf(outFile,", s_max = %g", s_max);
    fprintf(outFile,", s = %g", I);
    fprintf(outFile,", U_g = %g", U_g);
    fprintf(outFile,", U_c = %g", U_c);
    fprintf(outFile,", U_tm = %g", U_tm);
    fprintf(outFile,", U_tf = %g", U_tf);
    fprintf(outFile,", Rg = %g", Rgv);
    fprintf(outFile,", Rc = %g", Rcv);
    fprintf(outFile,", output = %d", outputv);
}
