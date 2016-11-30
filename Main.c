#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Boltzmann 1.38064852E-23
#define pi 3.1415926535897932

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
    double BeadRadi, FluidViscos, h, T, D, FlowVel, H;
    const int N;
} CONSTANTS;

//FUNCTION PROTOTYPES

POSITION CalcNextBallPos(POSITION nMinusTwoPos, POSITION nMinusOnePos);

double GenRandDouble(double minDoub, double maxDoub);                   //Note: not truly random, testing purposes only!

ANGLES CalcNextAngles(CONSTANTS c);

TWO_GAUSS BoxMullerTrans (CONSTANTS c, double input_1, double input_2);

FENE FENEForce(POSITION nMinusTwoPos, POSITION nMinusOnePos, POSITION Pos, CONSTANTS c);

double DragForce (CONSTANTS c);

BROWNIAN Brownian(CONSTANTS c);

void printFile(FILE *File_BeadPos, int N);

void update(POSITION nMinusOnePos, POSITION nMinusTwoPos, POSITION nPos, CONSTANTS c);


//MAIN PROG

int main()
{

    //double const radius = 1.0;
    

    //INITIALISE VARIABLES & CONSTANTS
    
    CONSTANTS c;
    c.BeadRadi = 0.000000002;       									//diameter of a DNA helix
    c.FluidViscos = 1;
    c.h = 0.01;        
    c.T = 298;                                                          //Temperature in Kelvin
    c.D = (Boltzmann * c.T) / (6 * pi * c.FluidViscos * c.BeadRadi);    //Diffusion coefficient
    c.H = (3*Boltzmann*c.T)/(1.8 E-9 * 4.7 E-12);       				//Taken from Simons paper, values for polystyrene not DNA
    c.N = 20;
    POSITION PositionArray[c.N];
    double radius = 1;
    
    FILE *File_BeadPos;

    File_BeadPos = fopen("File_BeadPos.txt", "w");

    PositionArray[0].xPos = 0;                  //Initialising pos1 to 0,0,0
    PositionArray[0].yPos = 0;
    PositionArray[0].zPos = 0;
    
    PositionArray[1].xPos = PositionArray[0].xPos + radius;
    PositionArray[1].yPos = PositionArray[0].yPos + radius;
    PositionArray[1].zPos = PositionArray[0].zPos + radius;

    for (int i = 0; i <= c.N; ++i)
    {
        PositionArray[i] = CalcNextBallPos(PositionArray[i-2], PositionArray[i-1]);
    }

    fprintf(File_BeadPos, "atom 0:%lf\tradius 1.0\tname S\n" ,  N);

    for(double t =0; t<100; t+=c.h){
        printFile(File_BeadPos, N);
        update(nMinusOnePos, nMinusTwoPos, nPos, c);
    }

    fclose(File_BeadPos);
    
    return 0;
}

//FUNCTIONS

POSITION CalcNextBallPos(POSITION nMinusTwoPos, POSITION nMinusOnePos)
{
    POSITION nPos;
    /*
    nPos.xPos = 2*nMinusOnePos.xPos - nMinusTwoPos.xPos;    //Minus contribution from angles, just use rad?
    nPos.yPos = 2*nMinusOnePos.yPos - nMinusTwoPos.yPos;
    nPos.zPos = 2*nMinusOnePos.zPos - nMinusTwoPos.zPos;
    */


    ANGLES nAngles = CalcNextAngles(c);

    nPos.xPos = nMinusOnePos.xPos + sin(nAngles.theta) * cos(nAngles.phi);
    nPos.yPos = nMinusOnePos.yPos + sin(nAngles.theta) * sin(nAngles.phi);      //This line appears to have the error
    nPos.zPos = nMinusOnePos.zPos + cos(nAngles.theta);

    return nPos;
}


double GenRandDouble(double minDoub, double maxDoub)
{
    double randDoub;
    double fraction = rand() / RAND_MAX;
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

    OutputGauss.Gauss_1 = sqrt(-2 * ln(input_1) ) * cos(2 * pi * input_2);          //If using a standard Gaussian
    OutputGauss.Gauss_1 = OutputGauss.Gauss_1 * sqrt(2 c.D * c.h) + 0;                     //Multiply by standard deviation and add mean (0) for our Gaussian

    OutputGauss.Gauss_2 = sqrt(-2 * ln(input_1) ) * sin(2 * pi * input_2);
    OutputGauss.Gauss_2 = OutputGauss.Gauss_2 * sqrt(2 c.D * c.h) + 0;

    return OutputGauss;
}


FENE FENEForce(POSITION nMinusTwoPos, POSITION nMinusOnePos, POSITION nPos, CONSTANTS c)
{
    FENE FENEForces;
    double Q_0 = 236.34E-9;

    FENEForces.FENE_x1 = (c.H * (nPos.xPos - nMinusOnePos.xPos)) / (1 - pow(nPos.xPos - nMinusOnePos.xPos, 2) / Q_0);
    FENEForces.FENE_x2 = (c.H * (nMinusOnePos.xPos - nMinusTwoPos.xPos)) / (1 - pow(nMinusOnePos.xPos - nMinusTwoPos.xPos, 2) / pow(Q_0, 2);
    FENEForces.FENE_y1 = (c.H * (nPos.yPos - nMinusOnePos.yPos)) / (1 - pow(nPos.yPos - nMinusOnePos.yPos, 2) / Q_0);
    FENEForces.FENE_y2 = (c.H * (nMinusOnePos.yPos - nMinusTwoPos.yPos)) / (1 - pow(nMinusOnePos.yPos - nMinusTwoPos.yPos, 2) / pow(Q_0, 2);
    FENEForces.FENE_z1 = (c.H * (nPos.zPos - nMinusOnePos.zPos)) / (1 - pow(nPos.zPos - nMinusOnePos.zPos, 2) / Q_0);
    FENEForces.FENE_z2 = (c.H * (nMinusOnePos.zPos - nMinusTwoPos.zPos)) / (1 - pow(nMinusOnePos.zPos - nMinusTwoPos.zPos, 2) / pow(Q_0, 2);

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

    BrownianForces.BrownianForce_x = sqrt(6 * c.D / c.h) * GenRandDouble(-1, 1);			//random number from 1 to -1
    BrownianForces.BrownianForce_y = sqrt(6 * c.D / c.h) * GenRandDouble(-1, 1);
    BrownianForces.BrownianForce_z = sqrt(6 * c.D / c.h) * GenRandDouble(-1, 1);

    return BrownianForces;
}

void printFile(FILE* File_BeadPos, int N){
    fprintf(File_BeadPos, "timestep\n");
    for(int i = 0; i <= N; i++)
    {        
        fprintf(File_BeadPos, "%lf\t%lf\t%lf\n\n", PositionArray[i].xPos, PositionArray[i].yPos, PositionArray[i].zPos);
    }   
}

void update(POSITION nMinusOnePos, POSITION nMinusTwoPos, POSITION nPos, CONSTANTS c){
    FENEForce(POSITION nMinusTwoPos, POSITION nMinusOnePos, POSITION nPos, CONSTANTS c);
    DragForce(CONSTANTS c);
    Brownian(CONSTANTS c);

}









































