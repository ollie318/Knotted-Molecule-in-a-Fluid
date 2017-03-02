#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Main.h"
#include "random.h"
#include "KnotAnalysis.h"

#define Boltzmann 1.38064852E-23
#define pi 3.1415926535897932
#define AvogadroNum 6.02E23

//MAIN PROG

int main(int argc, char *argv[]){

    CONSTANTS c;
    VEC* PositionArrayOld;
    VEC* PositionArrayNew;
    VEC* frames;
    VEC* FENEArray;
    VEC * BrownianArray;
    VEC * PotentialArray;

    char* paramfile = NULL;

    if(argc != 2){
      die("Please use: ./a.out <paramfile>\n", __LINE__, __FILE__);
    }
    else paramfile = argv[1];

    initialise(&c, &PositionArrayOld, &PositionArrayNew, &frames, &FENEArray, &BrownianArray, &PotentialArray, paramfile);

    int loopcount;
    for(loopcount = 0; loopcount < c.maxIters; loopcount++){

        updateFrames(c, loopcount, frames, PositionArrayOld);
        timestep(c, PositionArrayOld, PositionArrayNew, FENEArray, BrownianArray, PotentialArray);

        for(int j = 1; j < c.N; j ++){
            PositionArrayOld[j] = PositionArrayNew[j];
        }
    }

    writeVTF(c, frames);
    writeKnotAnalysis(c, frames);
    finalise(&c, &PositionArrayOld, &PositionArrayNew, &frames, &FENEArray, &BrownianArray, &PotentialArray);

    return 0;
}

int initialise(CONSTANTS* c, VEC** PositionArrayOld, VEC** PositionArrayNew, VEC** frames, VEC** FENEArray, VEC** BrownianArray, VEC** PotentialArray, const char* paramfile){
    // array constants
    FILE* params;
    params = fopen(paramfile, "r");

    if(params == NULL){
      die("File not found", __LINE__, __FILE__);
    }

    int readvalue;

    readvalue = fscanf(params, "%d\n", &(c->N));
    if(readvalue != 1) die("Could not read N", __LINE__, __FILE__);

    readvalue = fscanf(params, "%d\n", &(c->maxIters));
    if(readvalue != 1) die("Could not read number of iterations", __LINE__, __FILE__);

    readvalue = fscanf(params, "%lf\n", &(c->BeadRadi));
    if(readvalue != 1) die("Could not read bead radii", __LINE__, __FILE__);

    readvalue = fscanf(params, "%lf\n", &(c->FluidViscos));
    if(readvalue != 1) die("Could not read fluid viscosity", __LINE__, __FILE__);

    readvalue = fscanf(params, "%lf\n", &(c->FlowVel));
    if(readvalue != 1) die("Could not read flow velocity", __LINE__, __FILE__);

    readvalue = fscanf(params, "%lf\n", &(c->h));
    if(readvalue != 1) die("Could not read timestep", __LINE__, __FILE__);

    readvalue = fscanf(params, "%lf\n", &(c->T));
    if(readvalue != 1) die("Could not read temperature", __LINE__, __FILE__);

    c->D = (Boltzmann * c->T) / (6 * pi * c->FluidViscos * c->BeadRadi);    //Diffusion coefficient
    c->eta = 6*pi*c->FluidViscos*c->BeadRadi;

    //Spring coefficient calcs
    readvalue = fscanf(params, "%lf\n", &(c->N_k));
    if(readvalue != 1) die("Could not read N_k", __LINE__, __FILE__);

    readvalue = fscanf(params, "%lf\n", &(c->b_k));
    if(readvalue != 1) die("Could not read b_k", __LINE__, __FILE__);

    c->N_ks = c->N_k / (c->N-1);
    c->L_s = c->N_ks * c->b_k;
    c->H = (3*Boltzmann*c->T) / (c->L_s * c->b_k);                             //Taken from Simons paper, values for polystyrene not DNA

    //c.m = 1.9927E-26;
    c->m = 0.104 / AvogadroNum;                                          //Bead mass for styrene
    c->Q_0 = c->N_ks * c->b_k;

    c->MaxExtension = c->L_s;
    c->a = -475896;
    c->b = &(c->a);

    (*PositionArrayOld) = (VEC*) malloc(sizeof(VEC) * c->N);
    if(PositionArrayOld == NULL) die("cannot allocate memory for PositionArrayOld", __LINE__, __FILE__);

    (*PositionArrayNew) = (VEC*) malloc(sizeof(VEC) * c->N);
    if(PositionArrayNew == NULL) die("cannot allocate memory for PositionArrayNew", __LINE__, __FILE__);

    (*frames) = (VEC*) malloc(sizeof(VEC) * c->N * c->maxIters);
    if(frames == NULL) die("cannot allocate memory for frames", __LINE__, __FILE__);

    (*FENEArray) = (VEC*) malloc(sizeof(VEC) * c->N);
    if(FENEArray == NULL) die("cannot allocate memory for FENEArray", __LINE__, __FILE__);

    (*BrownianArray) = (VEC *) malloc(sizeof(VEC) * c->N);
    if(BrownianArray == NULL) die("cannot allocate memory for BrownianArray", __LINE__, __FILE__);

    (*PotentialArray) = (VEC *) malloc(sizeof(VEC) * c->N);
    if(PotentialArray == NULL) die("cannot allocate memory for PotentialArray", __LINE__, __FILE__);

    (*PositionArrayOld)[0].xcoord = 0;                                         //Initialising pos1 to 0,0,0
    (*PositionArrayOld)[0].ycoord = 0;
    (*PositionArrayOld)[0].zcoord = 0;

    CalcKnotPos(*c, *PositionArrayOld);

    return EXIT_SUCCESS;
}

int CalcKnotPos(CONSTANTS c, VEC* PositionArrayOld){

    FILE *knot;
    knot = fopen("8_19_32beads.txt", "r");

    if(knot == NULL) die("cannot find knot coordinate file", __LINE__, __FILE__);

    for(int i = 0; i < c.N; i++){
        double Testxcoord = 0.0, Testycoord = 0.0, Testzcoord = 0.0;
        int beadnumber;

        if((i >= 0 && i < 15) || (i > 45 && i < c.N )){
          if(i >= 0 && i < 15){
            PositionArrayOld[m].xcoord = i * (c.Q_0 * 0.8);
            PositionArrayOld[m].ycoord = 0.0;
            PositionArrayOld[m].zcoord = 0.0;
          }
          else{
            PositionArrayOld[m].xcoord = (i - 31) * (c.Q_0 * 0.8);
            PositionArrayOld[m].ycoord = 0.0;
            PositionArrayOld[m].zcoord = 0.0;
          }
        }

        else{
            fscanf(knot, "%d\t%lf\t%lf\t%lf", &beadnumber, &Testxcoord, &Testycoord, &Testzcoord);
            PositionArrayOld[m].xcoord = 15 * (c.Q_0 * 0.8) + Testxcoord;
            PositionArrayOld[m].ycoord = Testycoord;
            PositionArrayOld[m].zcoord = Testzcoord;
        }
    }

    return EXIT_SUCCESS;

}

int updateFrames(CONSTANTS c, int CurrentFrame, VEC* frames, VEC* positions){
    #pragma omp parallel for
    for (int i = 0; i < c.N; i++) {
        frames[CurrentFrame*c.N + i] = positions[m];
    }

    return EXIT_SUCCESS;
}

int timestep(CONSTANTS c, VEC* PositionArrayOld, VEC* PositionArrayNew, VEC* FENEArray, VEC* BrownianArray, VEC* PotentialArray){

    PositionArrayNew[0].xcoord = 0;
    PositionArrayNew[0].ycoord = 0;
    PositionArrayNew[0].zcoord = 0;

    // #pragma omp parallel for
    int i;
    for(i = 1; i < c.N; i ++) {

        if(i == c.N-1){
          FENEArray[m] = FENEForce(PositionArrayOld[i-1], PositionArrayOld[m], PositionArrayOld[m], c);
        }
        else FENEArray[m] = FENEForce(PositionArrayOld[i-1], PositionArrayOld[m], PositionArrayOld[i+1], c);
    }

    // #pragma omp parallel for
    int j;
    for(j = 1; j < c.N; j ++) {
        BrownianArray[j] = Brownian(c);
    }

    // #pragma omp parallel for
    int k;
    for(k = 1; k < c.N; k ++) {
        PotentialArray[k] = potential(c, PositionArrayOld, k);
    }

    // #pragma omp parallel for
    int m;
    for(m = 1; m < c.N; m ++) {
        PositionArrayNew[m].xcoord = PositionArrayOld[m].xcoord + c.h*((FENEArray[m].xcoord + BrownianArray[m].xcoord + PotentialArray[m].xcoord)/c.eta);

        PositionArrayNew[m].ycoord = PositionArrayOld[m].ycoord + c.h*((FENEArray[m].ycoord + BrownianArray[m].ycoord + PotentialArray[m].ycoord)/c.eta);

        PositionArrayNew[m].zcoord = PositionArrayOld[m].zcoord + c.h*((FENEArray[m].zcoord + BrownianArray[m].zcoord + PotentialArray[m].zcoord)/c.eta);
    }
    return EXIT_SUCCESS;
}

VEC FENEForce(VEC nMinusOnePos, VEC nPos, VEC nPosPlusOne, CONSTANTS c){
    VEC FENEForces;

    double TotalSep1 = sqrt(pow(nPosPlusOne.xcoord - nPos.xcoord, 2) + pow(nPosPlusOne.ycoord - nPos.ycoord, 2) + pow(nPosPlusOne.zcoord - nPos.zcoord, 2));

    double TotalSep2 = sqrt(pow(nPos.xcoord - nMinusOnePos.xcoord, 2) + pow(nPos.ycoord - nMinusOnePos.ycoord, 2) + pow(nPos.zcoord - nMinusOnePos.zcoord, 2));

    double Q_x1 = nPosPlusOne.xcoord - nPos.xcoord;
    double x1 = (c.H * Q_x1) / (1 - ((TotalSep1*TotalSep1)/ (c.MaxExtension * c.MaxExtension)));                             //Q_x1 is bond length for x & right, x2 left, y in y-dir etc

    double Q_x2 = nPos.xcoord - nMinusOnePos.xcoord;
    double x2 = (c.H * Q_x2) / (1 - ((TotalSep2*TotalSep2)/ (c.MaxExtension * c.MaxExtension)));                           //Values usually between -10 & 10 -- -297 on one run?

    double Q_y1 = nPosPlusOne.ycoord - nPos.ycoord;
    double y1 = (c.H * Q_y1) / (1 - ((TotalSep1*TotalSep1)/ (c.MaxExtension * c.MaxExtension)));

    double Q_y2 = nPos.ycoord - nMinusOnePos.ycoord;
    double y2 = (c.H * Q_y2) / (1 - ((TotalSep2*TotalSep2)/ (c.MaxExtension * c.MaxExtension)));

    double Q_z1 = nPosPlusOne.zcoord - nPos.zcoord;
    double z1 = (c.H * Q_z1) / (1 - ((TotalSep1*TotalSep1)/ (c.MaxExtension * c.MaxExtension)));

    double Q_z2 = nPos.zcoord - nMinusOnePos.zcoord;
    double z2 = (c.H * Q_z2) / (1 - ((TotalSep2*TotalSep2)/ (c.MaxExtension * c.MaxExtension)));

    FENEForces.xcoord = x1 - x2;
    FENEForces.ycoord = y1 - y2;
    FENEForces.zcoord = z1 - z2;

    return FENEForces;
}

VEC  Brownian(CONSTANTS c){
    VEC BrownianForces;

    BrownianForces.xcoord = sqrt((6 * Boltzmann * c.T*c.eta)/c.h) * GenGaussRand(c);            //random number from 1 to -1, Gaussian distribution
    BrownianForces.ycoord = sqrt((6 * Boltzmann * c.T*c.eta)/c.h) * GenGaussRand(c);
    BrownianForces.zcoord = sqrt((6 * Boltzmann * c.T*c.eta)/c.h) * GenGaussRand(c);            //On the scale E-2

    return BrownianForces;
}

VEC potential(CONSTANTS c, VEC* PositionArrayOld, int i){
    double sepX, sepY, sepZ, TotalSep, epsilon, sigma, potX, potY, potZ;
    VEC  pot;

    sigma = 2 * c.BeadRadi;                                           //r where attraction/repulsion changes
    epsilon = 1.0;                                                    //Depth of the weakly attractive well for atom

    pot.xcoord = 0.0;
    pot.ycoord = 0.0;
    pot.zcoord = 0.0;

    #pragma omp parallel for
    for(int j = 0; j < c.N; j++)
    {
        if(j == i){
          potX = 0.0;
          potY = 0.0;
          potZ = 0.0;
        }

        else{
        sepX = PositionArrayOld[m].xcoord - PositionArrayOld[j].xcoord;
        sepY = PositionArrayOld[m].ycoord - PositionArrayOld[j].ycoord;
        sepZ = PositionArrayOld[m].zcoord - PositionArrayOld[j].zcoord;
        TotalSep = sqrt( sepX*sepX + sepY*sepY + sepZ*sepZ );

        if(TotalSep != 0){
          potX = (5*(epsilon) * sepX * (pow(sigma, 5)/pow(TotalSep, 7)))/AvogadroNum;
          potY = (5*(epsilon) * sepY * (pow(sigma, 5)/pow(TotalSep, 7)))/AvogadroNum;
          potZ = (5*(epsilon) * sepZ * (pow(sigma, 5)/pow(TotalSep, 7)))/AvogadroNum;
        }

        else die("Molecule separation is 0", __LINE__, __FILE__);

        pot.xcoord += potX;
        pot.ycoord += potY;
        pot.zcoord += potZ;

      }

    }
    return pot;
}

int writeVTF(CONSTANTS c, VEC* frames){

    FILE* File_BeadPos;
    File_BeadPos = fopen("File_BeadPos.vtf", "w");
    if(File_BeadPos == NULL) die("Bead coordinate file could not be opened", __LINE__, __FILE__);

    fprintf(File_BeadPos, "atom 0:%d\tradius 1.0\tname S\n" , c.N);

    int k;
    for(k = 0; k < c.N-1; k++){
        fprintf(File_BeadPos, "bond %d:%d\n", k, k+1);
    }

    int j;								//j represents number of one bead
    for(j = 0; j < c.N*c.maxIters; j++)
    {
        if (j%c.N == 0)
            fprintf(File_BeadPos, "\ntimestep\n");
            fprintf(File_BeadPos, "%.14lf\t%.14lf\t%.14lf\n", 10E6 * frames[j].xcoord, 10E6 * frames[j].ycoord, 10E6 * frames[j].zcoord);		//Writes value in um, micrometres
    }

    fclose(File_BeadPos);

    return EXIT_SUCCESS;
}

int writeKnotAnalysis(CONSTANTS c, VEC* frames){
    FILE* KnotAnalysis;
    KnotAnalysis = fopen("KnotAnalysis.txt", "w");

    if(KnotAnalysis == NULL) die("Knot Analysis file could not be opened", __LINE__, __FILE__);


    double* chain = (double*) malloc( sizeof(double) * 3 * c.N );

    for(int j = 0; j<c.N*c.maxIters; j += 10*c.N){ // taking every 10th frame
      for(int i = 0; i < c.N; i++ ) { // taking each bead for that frame
        VEC p = frames[ j + i ];
        chain[3*i] = p.xcoord;
        chain[3*i + 1] = p.ycoord;
        chain[3*i + 2] = p.zcoord;
      }
    }

      fprintf(KnotAnalysis, "Start\tEnd\tPosition\n" );

    // jKN* PolymerKnot;
    // PolymerKnot = jKN_alloc(chain, c.N);
    // KnotScan(PolymerKnot);
    // if (PolymerKnot->state>0){
    //     double kstart, kend, kpos;
		// 		kstart = floor(0.5+PolymerKnot->start)+1;
		// 		kend = floor(0.5+PolymerKnot->end)+1;
		// 		kpos = floor(0.5+PolymerKnot->position)+1;
    //     fprintf(KnotAnalysis, "%lf\t%lf\t%lf\n", kstart, kend, kpos);
    // }
    // else fprintf(KnotAnalysis, "state = 0");

    fclose(KnotAnalysis);
    return EXIT_SUCCESS;
}

double GenGaussRand(CONSTANTS c){
    double OutputGauss_1;

    double input_1 = ran2(c.b);
    double input_2 = ran2(c.b);

    OutputGauss_1 = sqrt(-2 * log(input_1) ) * cos(2 * pi * input_2);          //If using a standard Gaussian

    return OutputGauss_1;
}

int finalise(CONSTANTS* c, VEC** PositionArrayOld, VEC** PositionArrayNew, VEC** frames, VEC** FENEArray, VEC** BrownianArray, VEC** PotentialArray){
    free(*PositionArrayOld);
    *PositionArrayOld = NULL;
    free(*PositionArrayOld);
    *PositionArrayOld = NULL;
    free(*frames);
    *frames = NULL;
    free(*FENEArray);
    *FENEArray = NULL;
    free(*BrownianArray);
    *BrownianArray = NULL;
    free(*PotentialArray);
    *PotentialArray = NULL;

    return EXIT_SUCCESS;
}

void die(const char* message, const int line, const char* file){
    fprintf(stderr, "%s. \nLine: %d\t File: %s\n", message, line, file);
    fflush(stderr);
    exit(EXIT_FAILURE);
}
