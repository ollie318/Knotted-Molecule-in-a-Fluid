#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Main.h"
// #include "mtwist.h"

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

int updateFrames(CONSTANTS c, int CurrentFrame, POSITION* frames, POSITION* positions)
{
    int i;
    for (i = 0; i < c.N; i++) {
        frames[CurrentFrame*c.N + i] = positions[i];
    }

    return EXIT_SUCCESS;
}

int timestep(CONSTANTS c, POSITION* PositionArrayOld, POSITION* PositionArrayNew)
{

    int i;

    PositionArrayNew[0].xPos = 0;
    PositionArrayNew[0].yPos = 0;
    PositionArrayNew[0].zPos = 0;

    for(i = 1; i < c.N-1; i ++) {
        PositionArrayNew[i] = Forces(PositionArrayOld[i-1], PositionArrayOld[i], PositionArrayOld[i+1], PositionArrayNew[i], PositionArrayNew, c, i);

    }

    PositionArrayNew[c.N-1] = ForcesLast(PositionArrayOld[c.N-2], PositionArrayOld[c.N-1], PositionArrayOld[c.N-1], PositionArrayNew[i], PositionArrayNew, c, i);

    return EXIT_SUCCESS;
}

int initialise(CONSTANTS* c, POSITION** PositionArrayOld, POSITION** PositionArrayNew, POSITION** frames)
{
    // array constants
    c->N = 45;
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

    ANGLES LastBondAngles;
    LastBondAngles.phi = 0;
    LastBondAngles.theta = pi/2;

    int i;
    for (i = 2; i < c->N; i++){
//         CalcKnotPos(*c, &((*PositionArrayOld)[i]), (*PositionArrayOld)[i-1], i);
        CalcNextBallPos((*PositionArrayOld)[i-1], &((*PositionArrayOld)[i]), LastBondAngles, *c);
    }

    return 0;
}

int CalcNextBallPos(POSITION nMinusOnePos, POSITION* nPos, ANGLES LastAngles, CONSTANTS c)
{

    nPos -> xPos = nMinusOnePos.xPos + ( c.Q_0 * 0.8 );
    nPos -> yPos = nMinusOnePos.yPos + ( c.Q_0 * 0.8 );
    nPos -> zPos = nMinusOnePos.zPos + ( c.Q_0 * 0.8 );

    return EXIT_SUCCESS;
}

ANGLES CalcNextAngles(CONSTANTS c)
{
    ANGLES NewAngles;
    NewAngles.theta = GenRandDouble(-pi/4, pi/4);
    NewAngles.phi = GenRandDouble(-pi, pi);

    return NewAngles;
}

double GenRandDouble(double minDoub, double maxDoub)
{
    double randDoub;
    double fraction = rand() / (RAND_MAX + 1.0);
    randDoub = minDoub + (maxDoub - minDoub) * fraction;
    return randDoub;
}


// double GenRandDoubleMT(CONSTANTS c)
// {
//   seedMT();
//
// 	float x1, x2, w, y1;
// 	static float y2;
// 	static int use_last = 0;
//
// 	if (use_last)		        /* use value from previous call */
// 	{
// 		y1 = y2;
// 		use_last = 0;
// 	}
// 	else
// 	{
// 		do {
// 			x1 = 2.0 * randomMT() - 1.0;
// 			x2 = 2.0 * randomMT() - 1.0;
// 			w = x1 * x1 + x2 * x2;
// 		} while ( w >= 1.0 );
//
// 		w = sqrt( (-2.0 * log( w ) ) / w );
// 		y1 = x1 * w;
// 		y2 = x2 * w;
// 		use_last = 1;
// 	}
//
// 	return( y1 * sqrt(2 * c.D * c.h));
// }


TWO_GAUSS BoxMullerTrans (CONSTANTS c, double input_1, double input_2)
{
    TWO_GAUSS OutputGauss;

    OutputGauss.Gauss_1 = sqrt(-2 * log(input_1) ) * cos(2 * pi * input_2);          //If using a standard Gaussian
    OutputGauss.Gauss_1 = OutputGauss.Gauss_1 * sqrt(2 * c.D * c.h);                 //Multiply by standard deviation and add mean (0) for our Gaussian

    OutputGauss.Gauss_2 =sqrt(-2 * log(input_1) ) * sin(2 * pi * input_2);
    OutputGauss.Gauss_2 = OutputGauss.Gauss_2 * sqrt(2 * c.D * c.h);

    return OutputGauss;
}


FENE FENEForce(POSITION nMinusOnePos, POSITION nPos, POSITION nPosPlusOne, CONSTANTS c)
{
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

    BrownianForces.BrownianForce_x = sqrt((6 * Boltzmann * c.T*c.eta)/c.h) * GenRandDouble(-1, 1);            //random number from 1 to -1
    BrownianForces.BrownianForce_y = sqrt((6 * Boltzmann * c.T*c.eta)/c.h) * GenRandDouble(-1, 1);			  //will need Gaussian dist number -1 to 1
    BrownianForces.BrownianForce_z = sqrt((6 * Boltzmann * c.T*c.eta)/c.h) * GenRandDouble(-1, 1);            //On the scale E-2

    return BrownianForces;
}

POSITION Forces(POSITION nMinusOnePos, POSITION nPosOld, POSITION nPosPlusOne, POSITION nPosNew, POSITION* PositionArrayNew,  CONSTANTS c, int i){

    BROWNIAN BrownianForces = Brownian(c);
    FENE FENEForces = FENEForce(nMinusOnePos, nPosOld, nPosPlusOne, c);
    POTENTIAL pot = potential(c, PositionArrayNew, i);

    nPosNew.xPos = nPosOld.xPos + (c.h*((FENEForces.FENE_x1 - FENEForces.FENE_x2)/c.eta  + (BrownianForces.BrownianForce_x/c.eta) + pot.potentialX/c.eta));
    // printf("%.12lf\n", (c.h*((FENEForces.FENE_x1 - FENEForces.FENE_x2)/c.eta )));

    nPosNew.yPos = nPosOld.yPos + (c.h*((FENEForces.FENE_y1 - FENEForces.FENE_y2)/c.eta  + (BrownianForces.BrownianForce_y/c.eta) + pot.potentialX/c.eta));

    nPosNew.zPos = nPosOld.zPos + (c.h*((FENEForces.FENE_z1 - FENEForces.FENE_z2)/c.eta  + (BrownianForces.BrownianForce_z/c.eta) + pot.potentialX/c.eta));
    return nPosNew;
}

POSITION ForcesLast(POSITION nMinusOnePos, POSITION nPosOld, POSITION nPosPlusOne, POSITION nPosNew, POSITION* PositionArrayNew, CONSTANTS c, int i){

    BROWNIAN BrownianForces = Brownian(c);
    FENE FENEForces = FENEForce(nMinusOnePos, nPosOld, nPosPlusOne, c);
    POTENTIAL pot = potential(c, PositionArrayNew, i);

    nPosNew.xPos = nPosOld.xPos + (c.h*((FENEForces.FENE_x1 - FENEForces.FENE_x2)/c.eta + (BrownianForces.BrownianForce_x/c.eta) + pot.potentialX/c.eta));

    nPosNew.yPos = nPosOld.yPos + (c.h*((FENEForces.FENE_y1 - FENEForces.FENE_y2)/c.eta + (BrownianForces.BrownianForce_y/c.eta) + pot.potentialX/c.eta));

    nPosNew.zPos = nPosOld.zPos + (c.h*((FENEForces.FENE_z1 - FENEForces.FENE_z2)/c.eta + (BrownianForces.BrownianForce_z/c.eta) + pot.potentialX/c.eta));

    return nPosNew;
}

POTENTIAL potential(CONSTANTS c, POSITION* PositionArrayNew, int i){
    double sepX, sepY, sepZ, epsilon, sigma, potX, potY, potZ;
    POTENTIAL pot;

    sigma = c.BeadRadi;
    epsilon = 1;

    pot.potentialX = 0.0;
    pot.potentialY = 0.0;
    pot.potentialZ = 0.0;

    int j;
    for(j = 0; j < c.N; j++)
    {

        sepX = PositionArrayNew[i].xPos - PositionArrayNew[j].xPos;
        sepY = PositionArrayNew[i].yPos - PositionArrayNew[j].yPos;
        sepZ = PositionArrayNew[i].zPos - PositionArrayNew[j].zPos;

        potX = -5*(epsilon/sigma) * pow(sigma, 5)/pow(sepX, 6);
        potY = -5*(epsilon/sigma) * pow(sigma, 5)/pow(sepY, 6);
        potZ = -5*(epsilon/sigma) * pow(sigma, 5)/pow(sepZ, 6);

        pot.potentialX += potX;
        pot.potentialY += potY;
        pot.potentialZ += potZ;
    }

    return pot;
}

int CalcKnotPos(CONSTANTS c, POSITION* PositionArrayNew, POSITION nMinusOnePos, int i){

    double phi_knot = 0;
    double p, q, r, TestxPos = 0.0, TestyPos = 0.0, TestzPos = 0.0 ;					//variables in knot eq, wiki
    p = 4.0;
    q = 3.0;							//For 8_19 knot

    double bondlength = 0.0;

    phi_knot = 0;

    if((i > 1 && i < 15) || (i > 30 && i < c.N )){
        PositionArrayNew[i].xPos = nMinusOnePos.xPos + (c.Q_0 * 0.8);
        PositionArrayNew[i].yPos = nMinusOnePos.yPos + (c.Q_0 * 0.8);
        PositionArrayNew[i].zPos = nMinusOnePos.zPos + (c.Q_0 * 0.8);
    }

    else{
        while(bondlength < c.Q_0*0.8){

            phi_knot += 0.01 ;
            r = cos(q * phi_knot) + 2;

            TestxPos = (r * cos(p * phi_knot))/3000000;
            TestyPos = (r * sin(p * phi_knot))/30000000;
            TestzPos = (- sin(q * phi_knot))/30000000;

            bondlength = sqrt(pow(TestxPos - PositionArrayNew[i-1].xPos, 2) + pow(TestyPos - PositionArrayNew[i-1].yPos, 2) + pow(TestzPos - PositionArrayNew[i-1].zPos, 2));

        }

        PositionArrayNew[i].xPos = nMinusOnePos.xPos + (c.Q_0 * 0.8);
        PositionArrayNew[i].yPos = TestyPos;
        PositionArrayNew[i].zPos = TestzPos;
    }
    return EXIT_SUCCESS;
}

int writeValues(CONSTANTS c, POSITION* frames)
{

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

int finalise(CONSTANTS* c, POSITION** PositionArrayOld, POSITION** PositionArrayNew, POSITION** frames)
{
    free(*PositionArrayOld);
    *PositionArrayOld = NULL;
    free(*PositionArrayOld);
    *PositionArrayOld = NULL;
    free(*frames);
    *frames = NULL;

    return EXIT_SUCCESS;
}
