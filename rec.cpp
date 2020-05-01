#include "sexchr.h"
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

extern MTRand rnd;

// recInit function: generates offspring from father and mother, burn-in generations (before recombination arrest):

void recInit(ind &offspring, ind &father, ind &mother, double Rg, double Rc, double ML, int nS, int off_sex)
{
    vector<double> posCoG;
    int i, j, locus, nbCo, ChrInit, strand;

    // paternally inherited gamete:
    
    double Rsex = Rg; // Recombination rate between the first cis reg and the sex locus. Only matters in males where the sex locus is heterozygous. Set to Rg to avoid adding a parameter
    nbCo = poisdev(ML+Rsex); // number of cross-overs
    for (j = 0; j < nbCo; j++)
        posCoG.push_back(rnd.randExc(ML+Rsex)); // position of cross-over j
    sort(posCoG.begin(), posCoG.end());

    if (off_sex == 0) // if offspring is a male, it inherits the sex-determining locus from the proto-Y
        ChrInit = 0;
    else
        ChrInit = 1;

    // CIS REG:
    if (Rc!=0.5)  //do this unless you have the Rc=0.5 mode, which enforces free recombination of cis regulators
   {
       locus = 0; //
       for (j = 0; j < nbCo; j++)
       {
           strand = ((j+ChrInit)%2)*nS;
           while ((locus < nS) && ((locus * Rg + Rsex) < posCoG[j])) // while locus is on the left of cross-over j; note that cis regulator i is at position Rg*i
           {
               offspring.cis[locus] = father.cis[strand + locus]; // copy from chrInit if nb of Co is even; copy from other chromosome otherwise; note that (j+chrInit)%2 equals 0 when chrInit = 0 and j is even or chrInit = 1 and j is odd
               locus++;
           }
       }
       strand = ((nbCo+ChrInit)%2)*nS;
       while (locus < nS) // cis regulators on the right of the last crossover; if total number of crossovers is even, all positions are taken from chrInit
       {
           offspring.cis[locus] = father.cis[strand + locus];
           locus++;
       }
   }

    // GENES:
    locus = 0;
    for (j = 0; j < nbCo; j++)
    {
        strand = ((j+ChrInit)%2)*nS;
        while ((locus < nS) && ((locus * Rg + Rc + Rsex) < posCoG[j])) // note that position of gene i is i*Rg + Rc
        {
            offspring.gene[locus] = father.gene[strand + locus]; // copy from chrInit if nb of Co is even; copy from other chromosome otherwise; note that (j+chrInit)%2 equals 0 when chrInit = 0 and j is even or chrInit = 1 and j is odd
            locus++;
        }
    }
    strand = ((nbCo+ChrInit)%2)*nS;
    while (locus < nS)
    {
        offspring.gene[locus] = father.gene[strand + locus];
        locus++;
    }

    // maternally inherited gamete:

    nbCo = poisdev(ML); // number of cross-overs
    for (j = 0; j < nbCo; j++)
        posCoG.push_back(rnd.randExc(ML)); // position of cross-over j
    sort(posCoG.begin(), posCoG.end());

    ChrInit = rnd.randInt(1); // which X chromosome contributes at position 0

    // CIS REG:
    if (Rc!=0.5) //do this unless you have the Rc=0.5 mode, which enforces free recombination of cis regulators
    {
        locus = 0;
        for (j = 0; j < nbCo; j++)
        {
            strand = ((j+ChrInit)%2)*nS;
            while ((locus < nS) && ((locus * Rg) < posCoG[j])) // while locus is on the left of cross-over j; note that cis regulator i is at position Rg*i
            {
                offspring.cis[nS+locus] = mother.cis[strand + locus]; // copy from chrInit if nb of Co is even; copy from other chromosome otherwise; note that (j+chrInit)%2 equals 0 when chrInit = 0 and j is even or chrInit = 1 and j is odd
                locus++;
            }
        }
        strand = ((nbCo+ChrInit)%2)*nS;
        while (locus < nS) // cis regulators on the right of the last crossover; if total number of crossovers is even, all positions are taken from chrInit
        {
            offspring.cis[nS+locus] = mother.cis[strand + locus];
            locus++;
        }
    }

    // GENES:
    locus = 0;
    for (j = 0; j < nbCo; j++)
    {
        strand = ((j+ChrInit)%2)*nS;
        while ((locus < nS) && ((locus * Rg + Rc) < posCoG[j])) // note that position of gene i is i*Rg + Rc
        {
            offspring.gene[nS+locus] = mother.gene[strand + locus]; // copy from chrInit if nb of Co is even; copy from other chromosome otherwise; note that (j+chrInit)%2 equals 0 when chrInit = 0 and j is even or chrInit = 1 and j is odd
            locus++;
        }
    }
    strand = ((nbCo+ChrInit)%2)*nS;
    while (locus < nS)
    {
        offspring.gene[nS+locus] = mother.gene[strand + locus];
        locus++;
    }

    // CIS REGULATORS in FREE REC MODE

    if (Rc == 0.5)
    {
        for (i = 0; i < nS; i++)
        {
            // paternally inherited cis modifiers:
            ChrInit = rnd.randInt(1);
            offspring.cis[i] = father.cis[ChrInit*nS+i];

            // maternally inherited cis modifiers:
            ChrInit = rnd.randInt(1);
            offspring.cis[i+nS] = mother.cis[ChrInit*nS+i];
        }
    }

    // TRANS REGULATORS

    for (i = 0; i < nS; i++)
    {
        // paternally inherited trans male modifiers:
        ChrInit = rnd.randInt(1);
        offspring.transm[i] = father.transm[ChrInit*nS+i];

        // maternally inherited trans modifiers:
        ChrInit = rnd.randInt(1);
        offspring.transm[i+nS] = mother.transm[ChrInit*nS+i];
    }

    for (i = 0; i < nS; i++)
    {
        // paternally inherited trans female modifiers:
        ChrInit = rnd.randInt(1);
        offspring.transf[i] = father.transf[ChrInit*nS+i];

        // maternally inherited trans modifiers:
        ChrInit = rnd.randInt(1);
        offspring.transf[i+nS] = mother.transf[ChrInit*nS+i];
    }
}




// rec function: generates offspring from father and mother, after recombination arrest:

void rec(ind &offspring, ind &father, ind &mother, double Rg, double Rc, double ML, int nS, int off_sex)
{
    vector<double> posCoG;
    int i, j, locus, nbCo, ChrInit, strand;

    if (off_sex == 0) // if offspring is a male
    {
        for (i = 0; i < nS; i++)
        {
            offspring.gene[i] = father.gene[i]; // father transmits his Y chromosome
            offspring.cis[i] = father.cis[i];
        }
    }
    else
    {
        for (i = 0; i < nS; i++)
        {
            offspring.gene[i] = father.gene[nS+i]; // father transmits his X chromosome
            offspring.cis[i] = father.cis[nS+i];
        }
    }

    nbCo = poisdev(ML); // number of cross-overs
    for (j = 0; j < nbCo; j++)
        posCoG.push_back(rnd.randExc(ML)); // position of cross-over j
    sort(posCoG.begin(), posCoG.end());

    ChrInit = rnd.randInt(1); // which X chromosome contributes at position 0

    // CIS REG:
    if (Rc!=0.5)
    {
        locus = 0;
        for (j = 0; j < nbCo; j++)
        {
            strand = ((j+ChrInit)%2)*nS;
            while ((locus < nS) && ((locus * Rg) < posCoG[j])) // while locus is on the left of cross-over j; note that cis regulator i is at position Rg*i
            {
                offspring.cis[nS+locus] = mother.cis[strand + locus]; // copy from chrInit if nb of Co is even; copy from other chromosome otherwise; note that (j+chrInit)%2 equals 0 when chrInit = 0 and j is even or chrInit = 1 and j is odd
                locus++;
            }
        }
        strand = ((nbCo+ChrInit)%2)*nS;
        while (locus < nS) // cis regulators on the right of the last crossover; if total number of crossovers is even, all positions are taken from chrInit
        {
            offspring.cis[nS+locus] = mother.cis[strand + locus];
            locus++;
        }
    }

    // GENES:
    locus = 0;
    for (j = 0; j < nbCo; j++)
    {
        strand = ((j+ChrInit)%2)*nS;
        while ((locus < nS) && ((locus * Rg + Rc) < posCoG[j])) // note that position of gene i is i*Rg + Rc
        {
            offspring.gene[nS+locus] = mother.gene[strand + locus]; // copy from chrInit if nb of Co is even; copy from other chromosome otherwise; note that (j+chrInit)%2 equals 0 when chrInit = 0 and j is even or chrInit = 1 and j is odd
            locus++;
        }
    }
    strand = ((nbCo+ChrInit)%2)*nS;
    while (locus < nS)
    {
        offspring.gene[nS+locus] = mother.gene[strand + locus];
        locus++;
    }

    // CIS REGULATORS
    if (Rc==0.5)
    {
        for (i = 0; i < nS; i++)
        {
            // maternally inherited cis regulators:
            ChrInit = rnd.randInt(1);
            offspring.cis[i+nS] = mother.cis[ChrInit*nS+i];
        }
    }

    // TRANS REGULATORS

    for (i = 0; i < nS; i++)
    {
        // paternally inherited trans male modifiers:
        ChrInit = rnd.randInt(1);
        offspring.transm[i] = father.transm[ChrInit*nS+i];

        // maternally inherited trans male modifiers:
        ChrInit = rnd.randInt(1);
        offspring.transm[i+nS] = mother.transm[ChrInit*nS+i];
    }
    for (i = 0; i < nS; i++)
    {
        // paternally inherited trans female modifiers:
        ChrInit = rnd.randInt(1);
        offspring.transf[i] = father.transf[ChrInit*nS+i];

        // maternally inherited trans female modifiers:
        ChrInit = rnd.randInt(1);
        offspring.transf[i+nS] = mother.transf[ChrInit*nS+i];
    }
}
