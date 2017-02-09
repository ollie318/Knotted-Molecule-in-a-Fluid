#ifndef MAIN_H
#define MAIN_H

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
    double BeadRadi, FluidViscos, h, T, D, FlowVel, H, m, Q_0;
    double N_k, N_ks, b_k, L_s, eta, MaxExtension;
    int N, maxIters;
} CONSTANTS;

typedef struct
{
    double potentialX, potentialY, potentialZ;
}POTENTIAL;

//FUNCTION PROTOTYPES

int CalcNextBallPos(POSITION nMinusOnePos, POSITION* nPos, ANGLES LastAngles, CONSTANTS c);

double GenRandDouble(double minDoub, double maxDoub);

double GenRandDoubleMT(CONSTANTS c);                   //Note: not truly random, testing purposes only!

ANGLES CalcNextAngles(CONSTANTS c);

int CalcKnotPos(CONSTANTS c, POSITION* PositionArrayNew, POSITION nMinusOnePos, int i);

TWO_GAUSS BoxMullerTrans (CONSTANTS c, double input_1, double input_2);

FENE FENEForce(POSITION nMinusOnePos, POSITION nPos, POSITION nPosPlusOne, CONSTANTS c);

double DragForce (CONSTANTS c);

BROWNIAN Brownian(CONSTANTS c);

POSITION Forces(POSITION nMinusOnePos, POSITION nPosOld, POSITION nPosPlusOne, POSITION nPosNew, POSITION* PositionArrayOld,  CONSTANTS c, int i);		//Gives new pos from sum forces

POSITION ForcesLast(POSITION nMinusOnePos, POSITION nPosOld, POSITION nPosPlusOne, POSITION nPosNew, POSITION* PositionArrayOld, CONSTANTS c, int i);


int initialise(CONSTANTS* c, POSITION** PositionArrayOld, POSITION** PositionArrayNew, POSITION** frames);			//Sets constants, allocates memory for array of pointers

int timestep(CONSTANTS c, POSITION* PositionArrayOld, POSITION* PositionArrayNew);									//Loop that calls Forces for each bead, will inc potential

int updateFrames(CONSTANTS c, int CurrentFrame, POSITION* frames, POSITION* positions);								//Writes current set of positions into frames, frames written to file at end

int writeValues(CONSTANTS c, POSITION* frames);																		//Writes all bead positions to file after arrays finished

int finalise(CONSTANTS* c, POSITION** PositionArrayOld, POSITION** PositionArrayNew, POSITION** frames);			//frees memory

POTENTIAL potential(CONSTANTS c, POSITION* PositionArrayOld, int i);																//Lennard-Jones Potential


#endif
