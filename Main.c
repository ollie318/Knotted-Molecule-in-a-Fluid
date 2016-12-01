#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Boltzmann 1.38064852E-23
#define pi 3.1415926535897932
#define N 20

//STRUCTURE/s

typedef struct
{
    double xPos, yPos, zPos;
} POSITION;

typedef struct
{
    double theta, phi;
} ANGLES;

typedef struct
{
    double Gauss_1, Gauss_2;
} TWO_GAUSS;

typedef struct
{
    double FENE_x1, FENE_x2, FENE_y1, FENE_y2, FENE_z1, FENE_z2;
} FENE;

typedef struct
{
    double BrownianForce_x, BrownianForce_y, BrownianForce_z;
} BROWNIAN;

typedef struct
{
    double BeadRadi, FluidViscos, h, T, D, FlowVel, H, m;
} CONSTANTS;

//FUNCTION PROTOTYPES

POSITION CalcNextBallPos(POSITION nMinusOnePos, POSITION nPos, CONSTANTS c);

double GenRandDouble(double minDoub, double maxDoub);                   //Note: not truly random, testing purposes only!

ANGLES CalcNextAngles(CONSTANTS c);

TWO_GAUSS BoxMullerTrans (CONSTANTS c, double input_1, double input_2);

FENE FENEForce(POSITION nMinusOnePos, POSITION nPos, POSITION Pos, CONSTANTS c);

double DragForce (CONSTANTS c);

BROWNIAN Brownian(CONSTANTS c);

void printFile(FILE *File_BeadPos, POSITION *PositionArray);

void update(POSITION nPos, POSITION nMinusOnePos, POSITION *nPosPlusOne, CONSTANTS c);


//MAIN PROG

int main()
{
    
    //double const radius = 1.0;
    
    
    //INITIALISE VARIABLES & CONSTANTS
    
    CONSTANTS c;
    c.BeadRadi = 0.002;                                         //diameter of a DNA helix
    c.FluidViscos = 1;
    c.FlowVel = 12;
    c.h = 0.01;
    c.T = 298;                                                          //Temperature in Kelvin
    c.D = (Boltzmann * c.T) / (6 * pi * c.FluidViscos * c.BeadRadi);    //Diffusion coefficient
    c.H = (3*Boltzmann*c.T)/(1.8E-9 * 4.7E-12);                     //Taken from Simons paper, values for polystyrene not DNA
    c.m = 1.9927E-26;
    POSITION PositionArray[N];
    double radius = 154E-12;

    FILE *File_BeadPos;
    
    File_BeadPos = fopen("File_BeadPos.vtf", "w");
    
    PositionArray[0].xPos = 0;                  //Initialising pos1 to 0,0,0
    PositionArray[0].yPos = 0;
    PositionArray[0].zPos = 0;
    
    PositionArray[1].xPos = PositionArray[0].xPos + radius;
    PositionArray[1].yPos = PositionArray[0].yPos + radius;
    PositionArray[1].zPos = PositionArray[0].zPos + radius;
    
    int i;
    for (i = 0; i <= N; i++)
    {
        PositionArray[i] = CalcNextBallPos(PositionArray[i-2], PositionArray[i-1], c);
    }
    
    fprintf(File_BeadPos, "atom 0:%d\tradius 1.0\tname S\n" ,  N);
    
    
    int loopcount;
    int noOfRuns;
    double t = 0;
    
    noOfRuns = 1000;
    
    for(loopcount = 0; loopcount < noOfRuns; loopcount++){
        t += c.h;
        int i =1;
        for(i = 1; i<=N; i++){
            update(PositionArray[i-1], PositionArray[i+1], &PositionArray[i], c);
        }
        printFile(File_BeadPos, PositionArray);
    }
    
    fclose(File_BeadPos);
    
    return 0;
}

//FUNCTIONS

POSITION CalcNextBallPos(POSITION nMinusOnePos, POSITION nPos, CONSTANTS c)
{
    POSITION nPosPlusOne;
    /*
     nPosPlusOne.xPos = 2*nPos.xPos - nMinusOnePos.xPos;    //Minus contribution from angles, just use rad?
     nPosPlusOne.yPos = 2*nPos.yPos - nMinusOnePos.yPos;
     nPosPlusOne.zPos = 2*nPos.zPos - nMinusOnePos.zPos;
     */
    
    
    ANGLES nAngles = CalcNextAngles(c);
    
    nPosPlusOne.xPos = nPos.xPos + sin(nAngles.theta) * cos(nAngles.phi);
    nPosPlusOne.yPos = nPos.yPos + sin(nAngles.theta) * sin(nAngles.phi);      //This line appears to have the error
    nPosPlusOne.zPos = nPos.zPos + cos(nAngles.theta);
    
    return nPosPlusOne;
}


double GenRandDouble(double minDoub, double maxDoub)
{
    double randDoub;
    double fraction = rand() / (RAND_MAX + 1.0);
    randDoub = minDoub + (maxDoub - minDoub) * fraction;
    return randDoub;
}


ANGLES CalcNextAngles(CONSTANTS c)
{
    ANGLES NewAngles;
    NewAngles.theta = GenRandDouble(-pi/4, pi/4);
    NewAngles.phi = GenRandDouble(-pi, pi);
    
    return NewAngles;
}


TWO_GAUSS BoxMullerTrans (CONSTANTS c, double input_1, double input_2)
{
    TWO_GAUSS OutputGauss;
    
    OutputGauss.Gauss_1 = sqrt(-2 * log(input_1) ) * cos(2 * pi * input_2);          //If using a standard Gaussian
    OutputGauss.Gauss_1 = OutputGauss.Gauss_1 * sqrt(2 * c.D * c.h);                     //Multiply by standard deviation and add mean (0) for our Gaussian
    
    OutputGauss.Gauss_2 = sqrt(-2 * log(input_1) ) * sin(2 * pi * input_2);
    OutputGauss.Gauss_2 = OutputGauss.Gauss_2 * sqrt(2 * c.D * c.h);
    
    return OutputGauss;
}


FENE FENEForce(POSITION nMinusOnePos, POSITION nPos, POSITION nPosPlusOne, CONSTANTS c)
{
    FENE FENEForces;
    double Q_0 = 236.34E-9;
    
    FENEForces.FENE_x1 = (c.H * (nPosPlusOne.xPos - nPos.xPos)) / (1 - pow(nPosPlusOne.xPos - nPos.xPos, 2) / Q_0);
    FENEForces.FENE_x2 = (c.H * (nPos.xPos - nMinusOnePos.xPos)) / (1 - pow(nPos.xPos - nMinusOnePos.xPos, 2) / pow(Q_0, 2));
    FENEForces.FENE_y1 = (c.H * (nPosPlusOne.yPos - nPos.yPos)) / (1 - pow(nPosPlusOne.yPos - nPos.yPos, 2) / Q_0);
    FENEForces.FENE_y2 = (c.H * (nPos.yPos - nMinusOnePos.yPos)) / (1 - pow(nPos.yPos - nMinusOnePos.yPos, 2) / pow(Q_0, 2));
    FENEForces.FENE_z1 = (c.H * (nPosPlusOne.zPos - nPos.zPos)) / (1 - pow(nPosPlusOne.zPos - nPos.zPos, 2) / Q_0);
    FENEForces.FENE_z2 = (c.H * (nPos.zPos - nMinusOnePos.zPos)) / (1 - pow(nPos.zPos - nMinusOnePos.zPos, 2) / pow(Q_0, 2));
    return FENEForces;
}


double DragForce (CONSTANTS c)
{
    double StokeForce;
    StokeForce = - 6 * pi * c.FluidViscos * c.FlowVel * c.BeadRadi;
    return StokeForce;
}

BROWNIAN Brownian(CONSTANTS c){
    BROWNIAN BrownianForces;
    
    BrownianForces.BrownianForce_x = sqrt(6 * c.D / c.h) * GenRandDouble(-1, 1);            //random number from 1 to -1
    BrownianForces.BrownianForce_y = sqrt(6 * c.D / c.h) * GenRandDouble(-1, 1);
    BrownianForces.BrownianForce_z = sqrt(6 * c.D / c.h) * GenRandDouble(-1, 1);
    
    return BrownianForces;
}

void printFile(FILE *File_BeadPos, POSITION *PositionArray){
    fprintf(File_BeadPos, "\ntimestep\n");
    int j;
    for(j = 0; j <= N; j++)
    {
        fprintf(File_BeadPos, "%lf\t%lf\t%lf\n", PositionArray[j].xPos, PositionArray[j].yPos, PositionArray[j].zPos);
    }
}

void update(POSITION nPos, POSITION nMinusOnePos, POSITION* nPosPlusOne, CONSTANTS c){
    BROWNIAN BrownianForces = Brownian(c);
    FENE FENEForces = FENEForce(nMinusOnePos, nPos, *nPosPlusOne, c);
    
    nPosPlusOne -> xPos += c.h*(FENEForces.FENE_x1 + BrownianForces.BrownianForce_x);
    
    nPosPlusOne -> yPos += c.h*(FENEForces.FENE_y1 + BrownianForces.BrownianForce_y);
    
    nPosPlusOne -> zPos += c.h*(FENEForces.FENE_z1 + BrownianForces.BrownianForce_z + DragForce(c));
}









































