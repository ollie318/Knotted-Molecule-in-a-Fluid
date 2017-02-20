#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Main.h"
#include "random.h"

#define Boltzmann 1.38064852E-23
#define pi 3.1415926535897932
#define AvogadroNum 6.02E23

//MAIN PROG

int main(void)
{
    CONSTANTS c;
    POSITION* PositionArrayOld;
    POSITION* PositionArrayNew;
    POSITION* frames;

    initialise(&c, &PositionArrayOld, &PositionArrayNew, &frames);

    int loopcount;
    for(loopcount = 0; loopcount < c.maxIters; loopcount++){

        updateFrames(c, loopcount, frames, PositionArrayOld);
        timestep(c, PositionArrayOld, PositionArrayNew);
        int j;
        for(j = 1; j < c.N; j ++){
            PositionArrayOld[j] = PositionArrayNew[j];
        }
    }

    writeValues(c, frames);
    finalise(&c, &PositionArrayOld, &PositionArrayNew, &frames);

    return 0;
}

int initialise(CONSTANTS* c, POSITION** PositionArrayOld, POSITION** PositionArrayNew, POSITION** frames){
    // array constants
    c->N = 60;
    c->maxIters = 10000;

    c->BeadRadi = 100E-9;                                                //diameter of polystyrene
    c->FluidViscos = 2;                                             //Fluid viscosity
    c->FlowVel = 13.0;                                                    //Fluid velocity, NOT relative velocity as needed for Stokes Law
    c->h = 0.001;                                                          //Time step
    c->T = 298;                                                          //Temperature in Kelvin
    c->D = (Boltzmann * c->T) / (6 * pi * c->FluidViscos * c->BeadRadi);    //Diffusion coefficient
    c->eta = 6*pi*c->FluidViscos*c->BeadRadi;

    //Spring coefficient calcs
    c->N_k = 2626;
    c->b_k = 1.8E-9;
    c->N_ks = c->N_k / (c->N-1);
    c->L_s = c->N_ks * c->b_k;
    c->H = (3*Boltzmann*c->T) / (c->L_s * c->b_k);                              //Taken from Simons paper, values for polystyrene not DNA

    //c.m = 1.9927E-26;
    c->m = 0.104 / AvogadroNum;                                          //Bead mass for styrene
    c->Q_0 = c->N_ks * c->b_k;

    c->MaxExtension = 5 * c->Q_0;

    (*PositionArrayOld) = (POSITION*) malloc(sizeof(POSITION) * c->N);
    (*PositionArrayNew) = (POSITION*) malloc(sizeof(POSITION) * c->N);
    (*frames) = (POSITION*) malloc(sizeof(POSITION) * c->N * c->maxIters);

    (*PositionArrayOld)[0].xPos = 0;                                         //Initialising pos1 to 0,0,0
    (*PositionArrayOld)[0].yPos = 0;
    (*PositionArrayOld)[0].zPos = 0;

    (*PositionArrayOld)[1].xPos = (*PositionArrayOld)[0].xPos + (0.8 * c->Q_0);
    (*PositionArrayOld)[1].yPos = (*PositionArrayOld)[0].yPos + (0.8 * c->Q_0);
    (*PositionArrayOld)[1].zPos = (*PositionArrayOld)[0].zPos + (0.8 * c->Q_0);

    CalcKnotPos(*c, &PositionArrayOld);

    return EXIT_SUCCESS;
}

int CalcKnotPos(CONSTANTS c, POSITION* PositionArrayOld){

    FILE *knot;
    knot = fopen("8_19_32beads.txt", "r");
    int beadnumber;
    double TestxPos = 0.0, TestyPos = 0.0, TestzPos = 0.0;

    if(knot != NULL){
      int i = 0;
      for(i = 0; i < c.N; i++){
          if((i >= 0 && i < 15) || (i > 47 && i < c.N )){
            if(i >= 0 && i < 15){
              PositionArrayOld[i].xPos = i * (c.Q_0 * 0.8);
              PositionArrayOld[i].yPos = 0.0;
              PositionArrayOld[i].zPos = 0.0;
            }
            else{
              PositionArrayOld[i].xPos = (i - 32) * (c.Q_0 * 0.8);
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
    else {
      printf("Error opening file\n");
      exit(0);
    }
}

int updateFrames(CONSTANTS c, int CurrentFrame, POSITION* frames, POSITION* positions){
    int i;
    for (i = 0; i < c.N; i++) {
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

    nPosNew.xPos = nPosOld.xPos + (c.h*((FENEForces.FENE_x1 - FENEForces.FENE_x2)/c.eta  + (BrownianForces.BrownianForce_x/c.eta) + (pot.potentialX/c.eta)));
    // printf("%.12lf\n", (c.h*((FENEForces.FENE_x1 - FENEForces.FENE_x2)/c.eta )));

    nPosNew.yPos = nPosOld.yPos + (c.h*((FENEForces.FENE_y1 - FENEForces.FENE_y2)/c.eta  + (BrownianForces.BrownianForce_y/c.eta) + (pot.potentialY/c.eta)));

    nPosNew.zPos = nPosOld.zPos + (c.h*((FENEForces.FENE_z1 - FENEForces.FENE_z2)/c.eta  + (BrownianForces.BrownianForce_z/c.eta) + (pot.potentialZ/c.eta)));
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

    double Q_x1 = nPosPlusOne.xPos - nPos.xPos;
    FENEForces.FENE_x1 = (c.H * Q_x1) / (1 - ((Q_x1 * Q_x1)/ (c.MaxExtension * c.MaxExtension)));                             //Q_x1 is bond length for x & right, x2 left, y in y-dir etc

    double Q_x2 = nPos.xPos - nMinusOnePos.xPos;
    FENEForces.FENE_x2 = (c.H * Q_x2) / (1 - ((Q_x2 * Q_x2)/ (c.MaxExtension * c.MaxExtension)));                           //Values usually between -10 & 10 -- -297 on one run?

    double Q_y1 = nPosPlusOne.yPos - nPos.yPos;
    FENEForces.FENE_y1 = (c.H * Q_y1) / (1 - ((Q_y1 * Q_y1)/ (c.MaxExtension * c.MaxExtension)));

    double Q_y2 = nPos.yPos - nMinusOnePos.yPos;
    FENEForces.FENE_y2 = (c.H * Q_y2) / (1 - ((Q_y2 * Q_y2)/ (c.MaxExtension * c.MaxExtension)));

    double Q_z1 = nPosPlusOne.zPos - nPos.zPos;
    FENEForces.FENE_z1 = (c.H * Q_z1) / (1 - ((Q_z1 * Q_z1)/ (c.MaxExtension * c.MaxExtension)));

    double Q_z2 = nPos.zPos - nMinusOnePos.zPos;
    FENEForces.FENE_z2 = (c.H * Q_z2) / (1 - ((Q_z2 * Q_z2)/ (c.MaxExtension * c.MaxExtension)));

    return FENEForces;
}

BROWNIAN Brownian(CONSTANTS c){
    BROWNIAN BrownianForces;

    BrownianForces.BrownianForce_x = sqrt((6 * Boltzmann * c.T*c.eta)/c.h) * GenGaussRand();            //random number from 1 to -1, Gaussian distribution
    BrownianForces.BrownianForce_y = sqrt((6 * Boltzmann * c.T*c.eta)/c.h) * GenGaussRand();
    BrownianForces.BrownianForce_z = sqrt((6 * Boltzmann * c.T*c.eta)/c.h) * GenGaussRand();            //On the scale E-2

    return BrownianForces;
}

double GenGaussRand(){
    double OutputGauss_1;

    double input_1 = ran2(long *idum);
    double input_2 = ran2(long *idum);

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

    int j;
    for(j = 0; j < c.N; j++)
    {
        if(j == i){
          potX = 0.0;
          potY = 0.0;
          potZ = 0.0;
        }

        else{
        sepX = PositionArrayOld[i].xPos - PositionArrayOld[j].xPos;
        sepY = PositionArrayOld[i].yPos - PositionArrayOld[j].yPos;
        sepZ = PositionArrayOld[i].zPos - PositionArrayOld[j].zPos;
        TotalSep = sqrt( sepX*sepX + sepY*sepY + sepZ*sepZ );

        potX = (4*(epsilon) * sepX * (pow(sigma, 4)/pow(TotalSep, 6)))/AvogadroNum;
        potY = (4*(epsilon) * sepY * (pow(sigma, 4)/pow(TotalSep, 6)))/AvogadroNum;
        potZ = (4*(epsilon) * sepZ * (pow(sigma, 4)/pow(TotalSep, 6)))/AvogadroNum;

        pot.potentialX += potX;
        pot.potentialY += potY;
        pot.potentialZ += potZ;

        }
    }
    return pot;
}

int writeValues(CONSTANTS c, POSITION* frames){

    FILE *File_BeadPos;
    File_BeadPos = fopen("File_BeadPos.vtf", "w");

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
