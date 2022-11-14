#include "sexchr.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <climits>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

extern MTRand rnd; // random number generator
extern FILE * outFile; // output file


/*----------------------------------------------------------
Function recursion: iterates the life cycle.
Parameters are:
Nv: population size
sigtrans: std deviation of distribution of mutational effects on trans regulators
sigcis: std deviation of distribution of mutational effects on cis regulators
nbSv: number of fitness loci per sex chromosome; equal to number of cis regulators
NbGenv: number of generations
NbPrelimv: number of preliminary generations with recombination between sex chromosomes
pasv: time interval between measurements
s: mean deleterious effect of mutations in fitness loci
smax: fitness effect of gene knock-out (maximal deleterious effect)
I: intensity of stabilizing selection on gene expression
U_g: mutation rate per fitness locus
U_c: mutation rate per cis regulator
U_tm: mutation rate per trans regulator in males
U_tf: mutation rate per trans regulator in females
Rg: genetic distance between genes
Rc: genetic distance between a cis-regulator and its regulated gene !!! RECOMBINATION FUNCTION ONLY VALID FOR Rc < Rg !!!!
output: 0 if all loci recorded
-----------------------------------------------------------*/



void recursion(int Nv, double sigtrans, double sigcis, int nbSv, int NbGenv, int NbPrelimv, int pasv, double s, double s_max, double I, double U_g, double U_c, double U_tm, double U_tf, double Rg, double Rc, int Rep, int output, double** allAverages)
{
    // variables:

	int i, j, k, nm, gen, mut, p1, p2, indiv, site, chrom, off_sex, Nmales_1, Nfemales_1, NjuvM, NjuvF;
	double w, wmmax, wfmax, h, dg, dc, dt, varm, moyc;

	double h0 = 0.25;
    double Q0 = 2;
    int N_1 = Nv - 1;
    int nbSm1 = nbSv - 1;
    int twonbSm1 = 2 * nbSv - 1;
    int Nmales = Nv / 2;
    double MLength = Rg * (nbSv - 1) + Rc;
    double UgTot = 2*Nv*U_g*nbSv;
    double UcTot = 2*Nv*U_c*nbSv;
    double UtmTot = 2*Nv*U_tm*nbSv;
    double UtfTot = 2*Nv*U_tf*nbSv;
    double expdom = -log(h0) / log(2.0);

    // creates result file (with parameter values in file name):

    char fileName[256];
    stringstream nomF;
    nomF << "N" << Nv << "_nbS" << nbSv << "_NbPrelim" << NbPrelimv << "_s" << s << "_I" << I << "_sigtrans" << sigtrans << "_sigcis" << sigcis << "_Ug" << U_g << "_Uc" << U_c << "_Utm" << U_tm << "_Utf" << U_tf << "_Rg" << Rg << "_Rc" << Rc << "_Rep" << Rep << ".txt";
    nomF >> fileName; // Writing "nomF" to nomFichier file
    ofstream fout;
    fout.open(fileName);

	// population of N individuals:
	ind * pop = new ind [Nv];
	ind * temp = new ind [Nv]; // temp for saving current generation during offspring production
    ind * pc;

    // array for measures from population: for each locus: sbarX, sbarY, hY, hXinact, attractXmale, attractXAct, attractXInact, attractY, dosY, dosXmale, dosXAct, dosXInact

    double ** measures;
    double * popAverages;

    // will hold measures for each locus: sbarX, sbarY, hY, attractX, cisx, cisy, transmale, transfemale (see record.cpp):

    measures = new double *[nbSv];
    for(i = 0; i < nbSv; i++)
        measures[i] = new double[8];

    // will hold averages over the whole population: Wbarmale, WbarFemale,
    // averages over loci of sbarX, sbarY, hY, attractX, cisx, cisy, transmale, transfemale,
    // number of loci half-dead, dead, half-silenced, silenced on Y (see record.cpp):

    if (output == 0)
        popAverages = new double [2];  // only Wbarmale, WbarFemale
    else
        popAverages = new double [14];

	// fitnesses:

	double * Wtot = new double [Nv];

	// for time length measure:

	time_t debut, fin; // "time_t" declares a time-calendar type; Debut-start, Fin-End
	struct tm *ptr; // sets "ptr" as a * for tm
	debut = time(0); // sets debut time to calendar time (not important)

    // initialization:

	for (i = 0; i < Nv; i++) // Loop through [Nv] times iterating over i (over each individual in the population)
    {
        // GENES
        pop[i].gene = new double [2*nbSv]; // two chromosomes
        temp[i].gene = new double [2*nbSv];

        // CIS regulators
        pop[i].cis = new double [2*nbSv]; // two chromosomes
        temp[i].cis = new double [2*nbSv];

        // TRANS regulators Effect in Male
        pop[i].transm = new double [2*nbSv]; // two chromosomes, trans reg's  corresponding to each gene
        temp[i].transm = new double [2*nbSv];

        // TRANS regulators Effect in female
        pop[i].transf = new double [2*nbSv]; // two chromosomes, trans reg's  corresponding to each gene
        temp[i].transf = new double [2*nbSv];


        for (k = 0; k < nbSv; k++)
        {
            pop[i].gene[k] = 1.0; // initialize all genes with 1 (WT/good) allele for Chromosome #1
            pop[i].gene[nbSv+k] = 1.0; // same for Chromosome #2
            pop[i].cis[k] = 1; // initialize all cis regulators
            pop[i].cis[nbSv+k] = 1; // same for Chromosome #2
            pop[i].transm[k] = 1; // initialize all trans male-effect regulators
            pop[i].transm[nbSv+k] = 1; // same for Chromosome #2
            pop[i].transf[k] = 1; // initialize all trans female-effect regulators
            pop[i].transf[nbSv+k] = 1; // same for Chromosome #2
         }

    }

  	// generations:
	for (gen = 0; gen <= NbGenv; gen++) // Iterates over all generations
	{
        // Draw mutations:

        // Gene mutations
        mut = int(poisdev(UgTot)); // 2NuL number of expected mutations in total population

        for (nm = 0; nm < mut; nm++) // iterates over number of muts
        {
            indiv = rnd.randInt(N_1); // selects random individual; "rnd" in "MersenneTwister.h"
            site = rnd.randInt(twonbSm1); // selects random loci

            // draw mutational effect on fitness from exp distribution
            dg = s * (-log(rnd.rand())); // Draw from exp dist of mean s, lambda = 1/s

            // For a set s_max limit
            if (pop[indiv].gene[site] * (1 - dg) >= 1-s_max)
                pop[indiv].gene[site] *= (1 - dg);
            else
                pop[indiv].gene[site] = 1-s_max; // sets floor to 1-s_max
        }

        // Cis mutations
        mut = int(poisdev(UcTot)); // 2NuL number of expected mutations in total population
        for (nm = 0; nm < mut; nm++)
        {
            indiv = rnd.randInt(N_1);
            site = rnd.randInt(twonbSm1);

            // draw mutational effect on fitness from Gaussian distribution
            dc = sigcis * gasdev();
            pop[indiv].cis[site] += dc;
        }

        // Trans mutations male
        mut = int(poisdev(UtmTot)); // 2NuL number of expected mutations in total population
        for (nm = 0; nm < mut; nm++)
        {
            indiv = rnd.randInt(N_1);
            site = rnd.randInt(twonbSm1);

            // draw mutational effect on fitness from Gaussian distribution
            dc = sigtrans * gasdev();
            pop[indiv].transm[site] += dc;
        }

         // Trans mutations female
        mut = int(poisdev(UtfTot)); // 2NuL number of expected mutations in total population
        for (nm = 0; nm < mut; nm++)
        {
            indiv = rnd.randInt(N_1);
            site = rnd.randInt(twonbSm1);

            // draw mutational effect on fitness from Gaussian distribution
            dc = sigtrans * gasdev();
            pop[indiv].transf[site] += dc;
        }

        // fitnesses:

        wmmax = 0;
        for (i = 0; i < Nmales; i++) // Loops through each male in the population
        {
            w = Wmale(pop[i], expdom, s_max, Q0, I, nbSv);
            Wtot[i] = w;
            if (wmmax < w)
                wmmax = w;
        }

        wfmax = 0;
        for (i = Nmales; i < Nv; i++) // Loops through each female in the population
        {
            w = Wfemale(pop[i], expdom, s_max, Q0, I, nbSv);
            Wtot[i] = w;
            if (wfmax < w)
                wfmax = w;
        }

        Nmales_1 = Nmales-1;
        Nfemales_1 = Nv-Nmales-1;
        NjuvM = 0;
        NjuvF = 0;


        // scaling cis and trans effects at the end of burnin generations

        if (gen == NbPrelimv)
        {
            record_output(pop, Wtot, measures, popAverages, nbSv, Nmales, Nv, expdom);

            for (i = 0; i < Nv; i++)
                for (k = 0; k < nbSv; k++)
                {
                    moyc = ((2.0*(Nv-Nmales)+Nmales)*measures[k][4] + Nmales*measures[k][5]) / (2.0*Nv);
                    pop[i].cis[k] = pop[i].cis[k] / moyc;
                    pop[i].cis[nbSv + k] = pop[i].cis[nbSv + k] / moyc;
                    pop[i].transm[k] = pop[i].transm[k] / measures[k][6];
                    pop[i].transm[nbSv + k] = pop[i].transm[nbSv + k] / measures[k][6];
                    pop[i].transf[k] = pop[i].transf[k] / measures[k][7];
                    pop[i].transf[nbSv + k] = pop[i].transf[nbSv + k] / measures[k][7];
                }
        }

        // measures phenotypic moments and writes in result file every "pasv" generations:

        if (gen % pasv == 0) // "gen" is the generational index (iterable) and once it reaches a multiple of the time interval per generation "pasv", write results
        {
            if (output == 0)  // records quantities for each locus individually (and writes in result file)
            {
                record_output(pop, Wtot, measures, popAverages, nbSv, Nmales, Nv, expdom);

                fout << gen << " " << popAverages[0] << " " << popAverages[1] << " " ;
                for (j = 0; j < nbSv; j++)
                    for (i = 0; i < 7; i++)
                        fout << measures[j][i] << " ";
                fout << endl;
            }
            else  // records only averages across loci (writes in result file, and update allAverages table averaging over replicates):
            {
                record_averages(pop, Wtot, popAverages, nbSv, Nmales, Nv, expdom, s_max);

                fout << gen << " ";
                for (i = 0; i < 14; i++)
                {
                    fout << popAverages[i] << " ";
                    allAverages[gen/pasv][i] = ((allAverages[gen/pasv][i])*Rep + popAverages[i]) / (Rep+1);
                }
                fout << endl;
            }
        }

        // next generation (meiosis):

		for (j = 0; j < Nv; j++)
		{
            do{
				p1 = rnd.randInt(Nmales_1); // selects the father

			} while (rnd.rand() > (Wtot[p1] / wmmax));

            do{
                p2 = rnd.randInt(Nfemales_1); // selects the mother

            } while (rnd.rand() > (Wtot[Nmales + p2] / wfmax));

            // Sex determination:

            if (rnd.rand() < 0.5) // offspring is male
                off_sex = 0;
            else
                off_sex = 1;

            // recombination:

            if (gen < NbPrelimv)
            {
                if (off_sex == 0)
                {
                    recInit(temp[NjuvM], pop[p1], pop[Nmales + p2], Rg, Rc, MLength, nbSv, off_sex);
                    NjuvM += 1;
                }
                else
                {
                    recInit(temp[Nv-1-NjuvF], pop[p1], pop[Nmales + p2], Rg, Rc, MLength, nbSv, off_sex);
                    NjuvF += 1;
                }
            }
            else
            {
                if (off_sex == 0)
                {
                    rec(temp[NjuvM], pop[p1], pop[Nmales + p2], Rg, Rc, MLength, nbSv, off_sex);
                    NjuvM += 1;
                }
                else
                {
                    rec(temp[Nv-1-NjuvF], pop[p1], pop[Nmales + p2], Rg, Rc, MLength, nbSv, off_sex);
                    NjuvF += 1;
                }
            }
		}

        // update population:

        pc = pop;
        pop = temp;
        temp = pc;

        Nmales = NjuvM;
	}

    fin = time(0); // NOT IMPORTANT (C function for getting the calendar time to write to the results file)

    // writes in output file:

    fprintf(outFile, "\n\nResults in file ");
    // fprintf(outFile, fileName );
    fprintf(outFile, "\n");

    // time length:

    int temps = int(difftime(fin, debut));
    fprintf(outFile, "\n%d generations took %d hour(s) %d minute(s) %d seconds\n", NbGenv, temps / 3600, (temps % 3600) / 60, temps % 60);

    // date and time:

    ptr=localtime(&fin); //  & = address of; localtime(&fin) gives the corresponding time zone time for "fin"
    // fprintf(outFile, asctime(ptr)); // ERIC: prints asctime(ptr) to the FILE stream with pointer fichierS; asctime(time) gives "time" in human-readable format

	// FREES MEMORY AFTER RUNNING
    for (i = 0; i < Nv; i++)
    {
        // delete [] "pointer" deletes/frees the array memory for the given pointer
        delete [] pop[i].gene;
        delete [] temp[i].gene;
        delete [] pop[i].transm;
        delete [] temp[i].transm;
        delete [] pop[i].transf;
        delete [] temp[i].transf;
        delete [] pop[i].cis;
        delete [] temp[i].cis;
    }
	delete [] pop;
    delete [] temp;
    delete [] Wtot;
    for(i = 0; i < nbSv; i++)
        delete [] measures[i];
    delete [] measures;
    delete [] popAverages;
}
