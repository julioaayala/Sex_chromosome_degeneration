#include "sexchr.h"
#include <cmath>
using namespace std;

void record_output(ind * pop, double * Wtot, double ** measures, double * popAve, int nbSv, int Nmales, int Nv, double expdom)
{
    // measures from population: for each locus: sbarX, sbarY, hY, attractX, cisx, cisy, transm, transf

    double WbarMales, WbarFemales, y, transmale,transfemale, c1, c2;
    int i, j;
    int Nfemales = Nv - Nmales;
    int twoNfemale = 2*(Nv - Nmales);

    for (j = 0; j < nbSv; j++)
        for (i = 0; i < 8; i++)
            measures[j][i] = 0;

    WbarMales = 0;
    for (i = 0; i < Nmales; i++)
    {
        WbarMales += Wtot[i];


        for (j = 0; j < nbSv; j++)
        {
            transmale = (((pop[i].transm[j] > 0) ? pop[i].transm[j] : 0)+((pop[i].transm[nbSv+j] > 0) ? pop[i].transm[nbSv+j] : 0))/2.0;
            transfemale = (((pop[i].transf[j] > 0) ? pop[i].transf[j] : 0)+((pop[i].transf[nbSv+j] > 0) ? pop[i].transf[nbSv+j] : 0))/2.0;
            c1 = ((pop[i].cis[j] > 0) ? pop[i].cis[j] : 0);
            c2 = ((pop[i].cis[nbSv+j] > 0) ? pop[i].cis[nbSv+j] : 0);

            measures[j][0] += 1 - pop[i].gene[nbSv+j];
            measures[j][1] += 1 - pop[i].gene[j];
            y = c1 / (c1 + c2);
            measures[j][2] += pow(y, expdom);
            measures[j][4] += c2;
            measures[j][5] += c1;
            measures[j][6] += transmale;
            measures[j][7] += transfemale;

        }
    }
    WbarMales /= Nmales;

    WbarFemales = 0;
    for (i = Nmales; i < Nv; i++)
    {
        WbarFemales += Wtot[i];


        for (j = 0; j < nbSv; j++)
        {
            transmale = (((pop[i].transm[j] > 0) ? pop[i].transm[j] : 0)+((pop[i].transm[nbSv+j] > 0) ? pop[i].transm[nbSv+j] : 0))/2.0;
            transfemale = (((pop[i].transf[j] > 0) ? pop[i].transf[j] : 0)+((pop[i].transf[nbSv+j] > 0) ? pop[i].transf[nbSv+j] : 0))/2.0;
            c1 = ((pop[i].cis[j] > 0) ? pop[i].cis[j] : 0);
            c2 = ((pop[i].cis[nbSv+j] > 0) ? pop[i].cis[nbSv+j] : 0);

            measures[j][0] += 2 - pop[i].gene[j] - pop[i].gene[nbSv+j];
            measures[j][3] += (c1 + c2) * transfemale;
            measures[j][4] += (c1 + c2) ;
            measures[j][6] += transmale;
            measures[j][7] += transfemale;

        }

    }
    for (j = 0; j < nbSv; j++)
    {
        measures[j][0] /= (Nmales + twoNfemale); // sbarX
        measures[j][1] /= Nmales; // sbarY
        measures[j][2] /= Nmales; // hY
        measures[j][3] /= twoNfemale; // attractX
        measures[j][4] /= (Nmales + twoNfemale); // cisX
        measures[j][5] /= Nmales; // cisY
        measures[j][6] /= (Nmales+Nfemales); // transmale
        measures[j][7] /= (Nmales+Nfemales); // transfemale

    }
    WbarFemales /= Nfemales;

    popAve[0] = WbarMales;
    popAve[1] = WbarFemales;

}

void record_averages(ind * pop, double * Wtot, double * popAve, int nbSv, int Nmales, int Nv, double expdom, double smax)
{
    // measures from population: Wbarmale, WbarFemale, sbarXbar, sbarYbar, hYbar, attractXbar, cisxbar, cisybar, transmbar, transfbar,
    // number of loci half-dead, dead, half-silenced, silenced on Y

    double WbarMales, WbarFemales, y, transmale,transfemale, c1, c2;
    int i, j;
    int twoNfemale = 2*(Nv - Nmales);
    double halfminW = (1-smax/2.0);
    double minW = 1-smax;
    int cmptyM = 0;


    for (i = 0; i < 14; i++)
        popAve[i] = 0;

    WbarMales = 0;
    for (i = 0; i < Nmales; i++)
    {
        WbarMales += Wtot[i];


        for (j = 0; j < nbSv; j++)
        {
            transmale = (((pop[i].transm[j] > 0) ? pop[i].transm[j] : 0)+((pop[i].transm[nbSv+j] > 0) ? pop[i].transm[nbSv+j] : 0))/2.0;
            transfemale = (((pop[i].transf[j] > 0) ? pop[i].transf[j] : 0)+((pop[i].transf[nbSv+j] > 0) ? pop[i].transf[nbSv+j] : 0))/2.0;
            c1 = ((pop[i].cis[j] > 0) ? pop[i].cis[j] : 0);
            c2 = ((pop[i].cis[nbSv+j] > 0) ? pop[i].cis[nbSv+j] : 0);

            popAve[2] += 1 - pop[i].gene[nbSv+j]; //sbarxbar
            popAve[3] += 1 - pop[i].gene[j];       //sbarybar
            if (c1 + c2 > 0)
            {
                y = c1 / (c1 + c2);
                popAve[4] += pow(y, expdom);    //hybar
                cmptyM++;
            }
            popAve[6] += c2;  //cisxbar
            popAve[7] += c1;  //cisybar
            popAve[8] += transmale;  //transm
            popAve[9] += transfemale;  //transm
            if (pop[i].gene[j] < halfminW)
                popAve[10] += 1;
            if (pop[i].gene[j] == minW)
                popAve[11] += 1;
            if (y < 0.25)
                popAve[12] += 1;
            if (y < 0.01)
                popAve[13] += 1;
        }
    }
    WbarMales /= Nmales;

    WbarFemales = 0;
    for (i = Nmales; i < Nv; i++)
    {
        WbarFemales += Wtot[i];


        for (j = 0; j < nbSv; j++)
        {
            transmale = (((pop[i].transm[j] > 0) ? pop[i].transm[j] : 0)+((pop[i].transm[nbSv+j] > 0) ? pop[i].transm[nbSv+j] : 0))/2.0;
            transfemale = (((pop[i].transf[j] > 0) ? pop[i].transf[j] : 0)+((pop[i].transf[nbSv+j] > 0) ? pop[i].transf[nbSv+j] : 0))/2.0;
            c1 = ((pop[i].cis[j] > 0) ? pop[i].cis[j] : 0);
            c2 = ((pop[i].cis[nbSv+j] > 0) ? pop[i].cis[nbSv+j] : 0);

            popAve[2] += 2 - pop[i].gene[j] - pop[i].gene[nbSv+j];
            popAve[5] +=  (c1 + c2) * transfemale ;
            popAve[6] +=  (c1 + c2)  ;
            popAve[8] += transmale;  //transm
            popAve[9] += transfemale;  //transm

        }

    }
    WbarFemales /= (Nv - Nmales);

    popAve[0] = WbarMales;
    popAve[1] = WbarFemales;
    popAve[2] /= (nbSv *(Nmales + twoNfemale)); // sbarXbar;
    popAve[3] /= (nbSv * Nmales); // sbarYbar ;
    popAve[4] /= cmptyM;
    popAve[5] /= (nbSv * twoNfemale); //attractXbar
    popAve[6] /= (nbSv *(Nmales + twoNfemale)); //cisxbar
    popAve[7] /= (nbSv *Nmales); // cisybar
    popAve[8] /= (nbSv * Nv); // transbarmale
    popAve[9] /= (nbSv * Nv); // transbarfemale
    popAve[10] /= (nbSv * Nmales);
    popAve[11] /= (nbSv * Nmales);
    popAve[12] /= (nbSv * Nmales);
    popAve[13] /= (nbSv * Nmales);
}



