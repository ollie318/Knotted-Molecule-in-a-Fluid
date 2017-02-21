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
    long a;
    long *b;
} CONSTANTS;

typedef struct
{
    double potentialX, potentialY, potentialZ;
}POTENTIAL;

//FUNCTION PROTOTYPES

int initialise(CONSTANTS* c, POSITION** PositionArrayOld, POSITION** PositionArrayNew, POSITION** frames, const char* paramfile);			//Sets constants, allocates memory for array of pointers

int CalcKnotPos(CONSTANTS c, POSITION* PositionArrayOld);

int updateFrames(CONSTANTS c, int CurrentFrame, POSITION* frames, POSITION* positions);								//Writes current set of positions into frames, frames written to file at end

int timestep(CONSTANTS c, POSITION* PositionArrayOld, POSITION* PositionArrayNew);									//Loop that calls Forces for each bead, will inc potential

POSITION Forces(POSITION nMinusOnePos, POSITION nPosOld, POSITION nPosPlusOne, POSITION nPosNew, POSITION* PositionArrayOld,  CONSTANTS c, int i);		//Gives new pos from sum forces

POSITION ForcesLast(POSITION nMinusOnePos, POSITION nPosOld, POSITION nPosPlusOne, POSITION nPosNew, POSITION* PositionArrayOld, CONSTANTS c, int i);

FENE FENEForce(POSITION nMinusOnePos, POSITION nPos, POSITION nPosPlusOne, CONSTANTS c);

BROWNIAN Brownian(CONSTANTS c);

double GenGaussRand (CONSTANTS c);

POTENTIAL potential(CONSTANTS c, POSITION* PositionArrayOld, int i);																//Lennard-Jones Potential

int writeValues(CONSTANTS c, POSITION* frames);																		//Writes all bead positions to file after arrays finished

int finalise(CONSTANTS* c, POSITION** PositionArrayOld, POSITION** PositionArrayNew, POSITION** frames);			//frees memory

void die(const char* message, const int line, const char* file);

#endif
