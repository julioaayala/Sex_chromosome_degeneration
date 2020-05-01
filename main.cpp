// main() function

#include "sexchr.h"
#include "MersenneTwister.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;

// random number generator (MTRand class defined in MersenneTwister.h)

MTRand rnd;

// pointers to input and output files:

FILE * inFile;
FILE * outFile;

int main()
{
	// variables:
    
    int repet,i, j, Nt, nbS, NbGen, NbPrelim, pas, Rep, output;
	double sigtrans,sigcis, s, s_max, I, U_g, U_c, U_tm, U_tf, Rg, Rc;

	// opens input and output files:

    bool fin;
    openInFile();
    openOutFile();
    fin = false;
    ofstream fout;

    // reads parameter values from input file, and writes them in output file:

    do
    {
        // reads parameters from input file:
        
        fin = readFile(Nt, sigtrans, sigcis, nbS, NbGen, NbPrelim, pas, s, s_max, I, U_g, U_c, U_tm, U_tf, Rg, Rc, Rep, output);
        
        if (!fin)
        {
            // allAverages table will contain averages (over replicates) of the 14 different quantities measured during the simulations,
            // at the different time steps (every "pas" generations). Filled by recursion function, only in output=1 mode:
            
            double** allAverages = new double *[(NbGen / pas) + 1];
            for(i = 0; i < (NbGen / pas) + 1; i++)
                allAverages[i] = new double[14];
            
            // initialization:

            for(i = 0; i < (NbGen / pas) + 1; i++)
                for (j = 0; j < 14; j++)
                    allAverages[i][j] = 0;
            
            // writes parameter values in output file:

            writeParameters(Nt, sigtrans,sigcis, nbS, NbGen, NbPrelim, pas, s, s_max, I, U_g, U_c, U_tm,U_tf, Rg, Rc, output);

            // runs the simulation (with Rep replicates):
            
            for (repet = 0; repet < Rep; repet++)
            {
                recursion(Nt, sigtrans,sigcis, nbS, NbGen, NbPrelim, pas, s, s_max, I, U_g, U_c, U_tm, U_tf, Rg, Rc, repet, output, allAverages);
                
                // writes current content of "allAverages" table in output file:

                char fileName[256];
                stringstream nomF;
                nomF << "summary_N" << Nt << "_nbS" << nbS << "_NbPrelim" << NbPrelim << "_s" << s << "_I" << I << "_sigtrans" << sigtrans << "_sigcis" << sigcis << "_Ug" << U_g << "_Uc" << U_c << "_Utm" << U_tm << "_Utf" << U_tf << "_Rg" << Rg << "_Rc" << Rc <<".txt"; // results file naming convention, parameter values in file name
                nomF >> fileName; // Writing "nomF" to nomFichier file
                fout.open(fileName);
                fout << repet+1 << endl;
                for(i = 0; i < (NbGen / pas) + 1; i++)
                {
                    fout << i*pas << " ";
                    for (j = 0; j < 14; j++)
                        fout << allAverages[i][j] << " ";
                    fout << endl;
                }
                fout.close();
            }

            for(i = 0; i < (NbGen/pas + 1); i++)
                delete [] allAverages[i];
            delete [] allAverages;
        }
    } while (!fin);

	// closes files:

	fclose(inFile);
	fclose(outFile);

	return 0;
}
