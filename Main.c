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
    POSITION* PositionArrayOld;
    POSITION* PositionArrayNew;
    POSITION* frames;
    char* paramfile = NULL;

    if(argc != 2){
      die("Please use: ./a.out <paramfile>\n", __LINE__, __FILE__);
    }
    else paramfile = argv[1];

    initialise(&c, &PositionArrayOld, &PositionArrayNew, &frames, paramfile);

    int loopcount;
    for(loopcount = 0; loopcount < c.maxIters; loopcount++){

        updateFrames(c, loopcount, frames, PositionArrayOld);
        timestep(c, PositionArrayOld, PositionArrayNew);
        int j;
        for(j = 1; j < c.N; j ++){
            PositionArrayOld[j] = PositionArrayNew[j];
        }
    }

    writeVTF(c, frames);
    writeKnotAnalysis(c, frames);
    finalise(&c, &PositionArrayOld, &PositionArrayNew, &frames);

    return 0;
}

int initialise(CONSTANTS* c, POSITION** PositionArrayOld, POSITION** PositionArrayNew, POSITION** frames, const char* paramfile){
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

    (*PositionArrayOld) = (POSITION*) malloc(sizeof(POSITION) * c->N);
    if(PositionArrayOld == NULL) die("cannot allocate memory for PositionArrayOld", __LINE__, __FILE__);

    (*PositionArrayNew) = (POSITION*) malloc(sizeof(POSITION) * c->N);
    if(PositionArrayNew == NULL) die("cannot allocate memory for PositionArrayNew", __LINE__, __FILE__);

    (*frames) = (POSITION*) malloc(sizeof(POSITION) * c->N * c->maxIters);
    if(frames == NULL) die("cannot allocate memory for frames", __LINE__, __FILE__);

    (*PositionArrayOld)[0].xPos = 0;                                         //Initialising pos1 to 0,0,0
    (*PositionArrayOld)[0].yPos = 0;
    (*PositionArrayOld)[0].zPos = 0;

    CalcKnotPos(*c, *PositionArrayOld);

    return EXIT_SUCCESS;
}

int CalcKnotPos(CONSTANTS c, POSITION* PositionArrayOld){

    FILE *knot;
    knot = fopen("8_19_32beads.txt", "r");

    if(knot == NULL) die("cannot find knot coordinate file", __LINE__, __FILE__);


    #pragma omp parallel for
    for(int i = 0; i < c.N; i++){
        double TestxPos = 0.0, TestyPos = 0.0, TestzPos = 0.0;
        int beadnumber;

        if((i >= 0 && i < 15) || (i > 45 && i < c.N )){
          if(i >= 0 && i < 15){
            PositionArrayOld[i].xPos = i * (c.Q_0 * 0.8);
            PositionArrayOld[i].yPos = 0.0;
            PositionArrayOld[i].zPos = 0.0;
          }
          else{
            PositionArrayOld[i].xPos = (i - 31) * (c.Q_0 * 0.8);
            PositionArrayOld[i].yPos = 0.0;
            PositionArrayOld[i].zPos = 0.0;
          }
        }

        else{
            fscanf(knot, "%d\t%lf\t%lf\t%lf", &beadnumber, &TestxPos, &TestyPos, &TestzPos);
            PositionArrayOld[i].xPos = 15 * (c.Q_0 * 0.8) + TestxPos;
            PositionArrayOld[i].yPos = TestyPos;
            PositionArrayOld[i].zPos = TestzPos;
        }
    }

    return EXIT_SUCCESS;

}

int updateFrames(CONSTANTS c, int CurrentFrame, POSITION* frames, POSITION* positions){
    #pragma omp parallel for
    for (int i = 0; i < c.N; i++) {
        frames[CurrentFrame*c.N + i] = positions[i];
    }

    return EXIT_SUCCESS;
}

int timestep(CONSTANTS c, POSITION* PositionArrayOld, POSITION* PositionArrayNew){

    int i;

    PositionArrayNew[0].xPos = 0;
    PositionArrayNew[0].yPos = 0;
    PositionArrayNew[0].zPos = 0;

    for(i = 1; i < c.N-1; i ++) {
        PositionArrayNew[i] = Forces(PositionArrayOld[i-1], PositionArrayOld[i], PositionArrayOld[i+1], PositionArrayNew[i], PositionArrayOld, c, i);

    }

    PositionArrayNew[c.N-1] = ForcesLast(PositionArrayOld[c.N-2], PositionArrayOld[c.N-1], PositionArrayOld[c.N-1], PositionArrayNew[i], PositionArrayOld, c, i);

    return EXIT_SUCCESS;
}

POSITION Forces(POSITION nMinusOnePos, POSITION nPosOld, POSITION nPosPlusOne, POSITION nPosNew, POSITION* PositionArrayOld,  CONSTANTS c, int i){

    BROWNIAN BrownianForces = Brownian(c);
    FENE FENEForces = FENEForce(nMinusOnePos, nPosOld, nPosPlusOne, c);
    POTENTIAL pot = potential(c, PositionArrayOld, i);

    nPosNew.xPos = nPosOld.xPos + c.h*((FENEForces.FENE_x1 - FENEForces.FENE_x2 + BrownianForces.BrownianForce_x + pot.potentialX)/c.eta);
    // printf("%.12lf\n", pot.potentialX);

    nPosNew.yPos = nPosOld.yPos + c.h*((FENEForces.FENE_y1 - FENEForces.FENE_y2 + BrownianForces.BrownianForce_y + pot.potentialY)/c.eta);

    nPosNew.zPos = nPosOld.zPos + c.h*((FENEForces.FENE_z1 - FENEForces.FENE_z2 + BrownianForces.BrownianForce_z + pot.potentialZ)/c.eta);
    return nPosNew;
}

POSITION ForcesLast(POSITION nMinusOnePos, POSITION nPosOld, POSITION nPosPlusOne, POSITION nPosNew, POSITION* PositionArrayOld, CONSTANTS c, int i){

    BROWNIAN BrownianForces = Brownian(c);
    FENE FENEForces = FENEForce(nMinusOnePos, nPosOld, nPosPlusOne, c);
    POTENTIAL pot = potential(c, PositionArrayOld, i);

    nPosNew.xPos = nPosOld.xPos + (c.h*((FENEForces.FENE_x1 - FENEForces.FENE_x2)/c.eta + (BrownianForces.BrownianForce_x/c.eta) + (pot.potentialX/c.eta)));

    nPosNew.yPos = nPosOld.yPos + (c.h*((FENEForces.FENE_y1 - FENEForces.FENE_y2)/c.eta + (BrownianForces.BrownianForce_y/c.eta) + (pot.potentialY/c.eta)));

    nPosNew.zPos = nPosOld.zPos + (c.h*((FENEForces.FENE_z1 - FENEForces.FENE_z2)/c.eta + (BrownianForces.BrownianForce_z/c.eta) + (pot.potentialZ/c.eta)));

    return nPosNew;
}

FENE FENEForce(POSITION nMinusOnePos, POSITION nPos, POSITION nPosPlusOne, CONSTANTS c){
    FENE FENEForces;

    double TotalSep1 = sqrt(pow(nPosPlusOne.xPos - nPos.xPos, 2) + pow(nPosPlusOne.yPos - nPos.yPos, 2) + pow(nPosPlusOne.zPos - nPos.zPos, 2));

    double TotalSep2 = sqrt(pow(nPos.xPos - nMinusOnePos.xPos, 2) + pow(nPos.yPos - nMinusOnePos.yPos, 2) + pow(nPos.zPos - nMinusOnePos.zPos, 2));

    double Q_x1 = nPosPlusOne.xPos - nPos.xPos;
    FENEForces.FENE_x1 = (c.H * Q_x1) / (1 - ((TotalSep1*TotalSep1)/ (c.MaxExtension * c.MaxExtension)));                             //Q_x1 is bond length for x & right, x2 left, y in y-dir etc

    double Q_x2 = nPos.xPos - nMinusOnePos.xPos;
    FENEForces.FENE_x2 = (c.H * Q_x2) / (1 - ((TotalSep2*TotalSep2)/ (c.MaxExtension * c.MaxExtension)));                           //Values usually between -10 & 10 -- -297 on one run?

    double Q_y1 = nPosPlusOne.yPos - nPos.yPos;
    FENEForces.FENE_y1 = (c.H * Q_y1) / (1 - ((TotalSep1*TotalSep1)/ (c.MaxExtension * c.MaxExtension)));

    double Q_y2 = nPos.yPos - nMinusOnePos.yPos;
    FENEForces.FENE_y2 = (c.H * Q_y2) / (1 - ((TotalSep2*TotalSep2)/ (c.MaxExtension * c.MaxExtension)));

    double Q_z1 = nPosPlusOne.zPos - nPos.zPos;
    FENEForces.FENE_z1 = (c.H * Q_z1) / (1 - ((TotalSep1*TotalSep1)/ (c.MaxExtension * c.MaxExtension)));

    double Q_z2 = nPos.zPos - nMinusOnePos.zPos;
    FENEForces.FENE_z2 = (c.H * Q_z2) / (1 - ((TotalSep2*TotalSep2)/ (c.MaxExtension * c.MaxExtension)));

    return FENEForces;
}

BROWNIAN Brownian(CONSTANTS c){
    BROWNIAN BrownianForces;

    BrownianForces.BrownianForce_x = sqrt((6 * Boltzmann * c.T*c.eta)/c.h) * GenGaussRand(c);            //random number from 1 to -1, Gaussian distribution
    BrownianForces.BrownianForce_y = sqrt((6 * Boltzmann * c.T*c.eta)/c.h) * GenGaussRand(c);
    BrownianForces.BrownianForce_z = sqrt((6 * Boltzmann * c.T*c.eta)/c.h) * GenGaussRand(c);            //On the scale E-2

    return BrownianForces;
}

double GenGaussRand(CONSTANTS c){
    double OutputGauss_1;

    double input_1 = ran2(c.b);
    double input_2 = ran2(c.b);

    OutputGauss_1 = sqrt(-2 * log(input_1) ) * cos(2 * pi * input_2);          //If using a standard Gaussian

    return OutputGauss_1;
}

POTENTIAL potential(CONSTANTS c, POSITION* PositionArrayOld, int i){

    double sepX, sepY, sepZ, TotalSep, epsilon, sigma, potX, potY, potZ;
    POTENTIAL pot;

    sigma = 2 * c.BeadRadi;                                           //r where attraction/repulsion changes
    epsilon = 1.0;                                                    //Depth of the weakly attractive well for atom

    pot.potentialX = 0.0;
    pot.potentialY = 0.0;
    pot.potentialZ = 0.0;

    // #pragma omp parallel for
    for(int j = 0; j < c.N; j++)
    {

      sepX = PositionArrayOld[i].xPos - PositionArrayOld[j].xPos;
      sepY = PositionArrayOld[i].yPos - PositionArrayOld[j].yPos;
      sepZ = PositionArrayOld[i].zPos - PositionArrayOld[j].zPos;
      TotalSep = sqrt( sepX*sepX + sepY*sepY + sepZ*sepZ );

      potX = (5*(epsilon) * sepX * (pow(sigma, 5)/pow(TotalSep, 7)))/AvogadroNum;
      potY = (5*(epsilon) * sepY * (pow(sigma, 5)/pow(TotalSep, 7)))/AvogadroNum;
      potZ = (5*(epsilon) * sepZ * (pow(sigma, 5)/pow(TotalSep, 7)))/AvogadroNum;

      if(j == i){
        potX = 0.0;
        potY = 0.0;
        potZ = 0.0;
      }

      pot.potentialX += potX;
      pot.potentialY += potY;
      pot.potentialZ += potZ;



    }
    return pot;
}

int writeVTF(CONSTANTS c, POSITION* frames){

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
            fprintf(File_BeadPos, "%.14lf\t%.14lf\t%.14lf\n", 10E6 * frames[j].xPos, 10E6 * frames[j].yPos, 10E6 * frames[j].zPos);		//Writes value in um, micrometres
    }

    fclose(File_BeadPos);

    return EXIT_SUCCESS;
}

int writeKnotAnalysis(CONSTANTS c, POSITION* frames){
    FILE* KnotAnalysis;
    KnotAnalysis = fopen("KnotAnalysis.txt", "w");

    if(KnotAnalysis == NULL) die("Knot Analysis file could not be opened", __LINE__, __FILE__);


    double* chain = (double*) malloc( sizeof(double) * 3 * c.N );

    #pragma omp parallel for
    for(int j = 0; j<c.N*c.maxIters; j += 10*c.N){ // taking every 10th frame
      for(int i = 0; i < c.N; i++ ) { // taking each bead for that frame
        POSITION p = frames[ j + i ];
        chain[3*i] = p.xPos;
        chain[3*i + 1] = p.yPos;
        chain[3*i + 2] = p.zPos;
      }
    }

      fprintf(KnotAnalysis, "Start\tEnd\tPosition\n" );

    jKN* PolymerKnot;
    PolymerKnot = jKN_alloc(chain, c.N);
    KnotScan(PolymerKnot);
    if (PolymerKnot->state>0){
        double kstart, kend, kpos;
				kstart = floor(0.5+PolymerKnot->start)+1;
				kend = floor(0.5+PolymerKnot->end)+1;
				kpos = floor(0.5+PolymerKnot->position)+1;
        fprintf(KnotAnalysis, "%lf\t%lf\t%lf\n", kstart, kend, kpos);
    }
    else fprintf(KnotAnalysis, "state = 0");

    fclose(KnotAnalysis);
    return EXIT_SUCCESS;
}

int finalise(CONSTANTS* c, POSITION** PositionArrayOld, POSITION** PositionArrayNew, POSITION** frames){
    free(*PositionArrayOld);
    *PositionArrayOld = NULL;
    free(*PositionArrayOld);
    *PositionArrayOld = NULL;
    free(*frames);
    *frames = NULL;

    return EXIT_SUCCESS;
}

void die(const char* message, const int line, const char* file){
    fprintf(stderr, "%s. \nLine: %d\t File: %s\n", message, line, file);
    fflush(stderr);
    exit(EXIT_FAILURE);
}
