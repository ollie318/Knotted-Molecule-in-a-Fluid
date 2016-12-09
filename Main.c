#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Boltzmann 1.38064852E-23
#define pi 3.1415926535897932
#define N 20
#define AvogadroNum 6.02E23

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
    double BeadRadi, FluidViscos, h, T, D, FlowVel, H, m, Q_0, N_k, N_ks, b_k, L_s;
} CONSTANTS;


//FUNCTION PROTOTYPES

POSITION CalcNextBallPos(POSITION nMinusOnePos, ANGLES LastAngles, CONSTANTS c);

double GenRandDouble(double minDoub, double maxDoub);                   //Note: not truly random, testing purposes only!

ANGLES CalcNextAngles(CONSTANTS c);

TWO_GAUSS BoxMullerTrans (CONSTANTS c, double input_1, double input_2);

FENE FENEForce(POSITION nMinusOnePos, POSITION nPos, POSITION nPosPlusOne, CONSTANTS c);

double DragForce (CONSTANTS c);

BROWNIAN Brownian(CONSTANTS c);

void printFile(FILE *File_BeadPos, POSITION *PositionArrayNew);

POSITION update(POSITION nMinusOnePos, POSITION nPosOld, POSITION nPosPlusOne, POSITION nPosNew, CONSTANTS c);


//MAIN PROG

int main()
{
    
    
    //INITIALISE VARIABLES & CONSTANTS
    
    
    
    FILE *File_BeadPos;
    
    File_BeadPos = fopen("File_BeadPos.vtf", "w");
    
    CONSTANTS c;
    c.BeadRadi = 100E-9;                                                //diameter of polystyrene
    c.FluidViscos = 8.9E-4;                                             //Fluid viscosity
    c.FlowVel = 1.2;                                                    //Fluid velocity, NOT relative velocity as needed for Stokes Law
    c.h = 0.1;                                                          //Time step?
    c.T = 298;                                                          //Temperature in Kelvin
    c.D = (Boltzmann * c.T) / (6 * pi * c.FluidViscos * c.BeadRadi);    //Diffusion coefficient
    
    //Spring coefficient calcs
    c.N_k = 2626;
    c.b_k = 1.8E-9;
    c.N_ks = c.N_k/N;
    c.L_s = c.N_ks*c.b_k;
    c.H = (3*Boltzmann*c.T)/(c.L_s*c.b_k);                              //Taken from Simons paper, values for polystyrene not DNA
    
    //c.m = 1.9927E-26;
    c.m = 0.104 / AvogadroNum;                                          //Bead mass for styrene
    c.Q_0 = 236.34E-9;                                                  //Average bond length
    POSITION PositionArrayOld[N];
    POSITION PositionArrayNew[N];
    
    PositionArrayOld[0].xPos = 0;                                          //Initialising pos1 to 0,0,0
    PositionArrayOld[0].yPos = 0;
    PositionArrayOld[0].zPos = 0;
    
    PositionArrayOld[1].xPos = PositionArrayOld[0].xPos + (0.8 * c.Q_0);
    PositionArrayOld[1].yPos = PositionArrayOld[0].yPos;
    PositionArrayOld[1].zPos = PositionArrayOld[0].zPos;
    
    int i;
    ANGLES LastBondAngles;
    LastBondAngles.phi = 0;
    LastBondAngles.theta = pi/2;                                                             //Writes initial positions
    for (i = 2; i <= N; i++)                                            //i starts at 2 because 0 & 1 already set
    {
        PositionArrayOld[i] = CalcNextBallPos(PositionArrayOld[i-1], LastBondAngles, c);
    }
    
    
    //WRITES TEXT NEEDED FOR VMD
    fprintf(File_BeadPos, "atom 0:%d\tradius 1.0\tname S\n" ,  N);
    
    int k;
    for(k = 0; k < N; k++){
        fprintf(File_BeadPos, "bond %d:%d\n", k, k+1);
    }
    
    
    printFile(File_BeadPos, PositionArrayOld);
    
    //LOOP FOR TIMESTEPPING
    int loopcount;
    int noOfRuns;
    double t = 0;
    
    noOfRuns = 10000;
    
    int i_dash;
    
    for(loopcount = 0; loopcount < noOfRuns; loopcount++){
        t += c.h;
        
        for(i_dash = 1; i_dash <= N; i_dash ++){
            
            PositionArrayNew[i_dash] = update(PositionArrayOld[i_dash-1], PositionArrayOld[i_dash], PositionArrayOld[i_dash+1], PositionArrayNew[i_dash], c);
        }
        
        printFile(File_BeadPos, PositionArrayNew);
        
        int j;
        for(j = 1; j <= N; j ++){
            PositionArrayOld[j] = PositionArrayNew[j];
        }

    }
    
    fclose(File_BeadPos);
    
    return 0;
}

//FUNCTIONS

POSITION CalcNextBallPos(POSITION nMinusOnePos, ANGLES LastAngles, CONSTANTS c)
{
    POSITION nPos;
    /*
     nPosPlusOne.xPos = 2*nPos.xPos - nMinusOnePos.xPos;    //Minus contribution from angles, just use rad?
     nPosPlusOne.yPos = 2*nPos.yPos - nMinusOnePos.yPos;
     nPosPlusOne.zPos = 2*nPos.zPos - nMinusOnePos.zPos;
     */
    
    
    ANGLES nAngles = CalcNextAngles(c);
    
    nPos.xPos = nMinusOnePos.xPos + ( c.Q_0 * (sin(LastAngles.theta + nAngles.theta) * cos(LastAngles.phi + nAngles.phi)) );
    nPos.yPos = nMinusOnePos.yPos + ( c.Q_0 * (sin(LastAngles.theta + nAngles.theta) * sin(LastAngles.phi + nAngles.phi)) );
    nPos.zPos = nMinusOnePos.zPos + ( c.Q_0 * (cos(LastAngles.theta + nAngles.theta)) );
    
    return nPos;
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
    OutputGauss.Gauss_1 = OutputGauss.Gauss_1 * sqrt(2 * c.D * c.h);                 //Multiply by standard deviation and add mean (0) for our Gaussian
    
    OutputGauss.Gauss_2 = sqrt(-2 * log(input_1) ) * sin(2 * pi * input_2);
    OutputGauss.Gauss_2 = OutputGauss.Gauss_2 * sqrt(2 * c.D * c.h);
    
    return OutputGauss;
}


FENE FENEForce(POSITION nMinusOnePos, POSITION nPos, POSITION nPosPlusOne, CONSTANTS c)
{
    FENE FENEForces;
    
    double Q_x1 = nPosPlusOne.xPos - nPos.xPos;                                                             //Qs are current spring length, non-equilibrium
    FENEForces.FENE_x1 = (c.H * Q_x1) / (1 - ( pow(Q_x1, 2) / pow(c.Q_0, 2) ));                             //Q_x1 is bond length for x & right, x2 left, y in y-dir etc
    
    double Q_x2 = nPos.xPos - nMinusOnePos.xPos;
    FENEForces.FENE_x2 = (c.H * (Q_x2)) / (1 - ( pow(Q_x2, 2) / pow(c.Q_0, 2) ));                           //Values usually between -10 & 10 -- -297 on one run?
    
    double Q_y1 = nPosPlusOne.yPos - nPos.yPos;
    FENEForces.FENE_y1 = (c.H * Q_y1) / (1 - (pow( Q_y1, 2 ) / pow (c.Q_0, 2)));
    
    double Q_y2 = nPos.yPos - nMinusOnePos.yPos;
    FENEForces.FENE_y2 = (c.H * Q_y2) / (1 - (pow( Q_y2, 2 ) / pow (c.Q_0, 2)));
    
    double Q_z1 = nPosPlusOne.zPos - nPos.zPos;
    FENEForces.FENE_z1 = (c.H * Q_z1) / (1 - (pow( Q_z1, 2 ) / pow (c.Q_0, 2)));
    
    double Q_z2 = nPos.zPos - nMinusOnePos.zPos;
    FENEForces.FENE_z2 = (c.H * Q_z2) / (1 - (pow( Q_z2, 2 ) / pow (c.Q_0, 2)));
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
    BrownianForces.BrownianForce_z = sqrt(6 * c.D / c.h) * GenRandDouble(-1, 1);            //On the scale E-2
    
    return BrownianForces;
}

void printFile(FILE *File_BeadPos, POSITION *PositionArrayNew){
    fprintf(File_BeadPos, "\ntimestep\n");
    int j;
    for(j = 0; j <= N; j++)
    {
        fprintf(File_BeadPos, "%.14lf\t%.14lf\t%.14lf\n", PositionArrayNew[j].xPos, PositionArrayNew[j].yPos, PositionArrayNew[j].zPos);
    }
}

POSITION update(POSITION nMinusOnePos, POSITION nPosOld, POSITION nPosPlusOne, POSITION nPosNew, CONSTANTS c){
    
    BROWNIAN BrownianForces = Brownian(c);
    FENE FENEForces = FENEForce(nMinusOnePos, nPosOld, nPosPlusOne, c);

    double StokeDragForce = DragForce(c);
    
    nPosNew.xPos = nPosOld.xPos + c.h*( FENEForces.FENE_x1 - FENEForces.FENE_x2 + BrownianForces.BrownianForce_x);            //x INCREASES by FENE acc (a = F/m) and Brownian acc
    
    nPosNew.yPos = nPosOld.yPos + c.h*( FENEForces.FENE_y1 - FENEForces.FENE_y2 + BrownianForces.BrownianForce_y);
    
    nPosNew.zPos = nPosOld.zPos + c.h*( FENEForces.FENE_z1 - FENEForces.FENE_z2 + BrownianForces.BrownianForce_z + StokeDragForce);
    
    return nPosNew;
}









































