#ifndef MAIN_H
#define MAIN_H

typedef struct
{
    double xcoord, ycoord, zcoord;
} VEC;

typedef struct
{
    double BeadRadi, FluidViscos, h, T, D, FlowVel, H, m, Q_0;
    double N_k, N_ks, b_k, L_s, eta, MaxExtension, StokesFric;
    int N, maxIters;
    long a;
    long *b;
} CONSTANTS;


//FUNCTION PROTOTYPES

int initialise(CONSTANTS* c, VEC** PositionArrayOld, VEC** PositionArrayNew, VEC** frames, VEC** VECArray, VEC** BrownianArray, VEC** PotentialArray, const char* paramfile);			//Sets constants, allocates memory for array of pointers

int CalcKnotPos(CONSTANTS c, VEC* PositionArrayOld);

int updateFrames(CONSTANTS c, int CurrentFrame, VEC* frames, VEC* positions);								//Writes current set of positions into frames, frames written to file at end

int timestep(CONSTANTS c, VEC* PositionArrayOld, VEC* PositionArrayNew, VEC* VECArray, VEC* BrownianArray, VEC* PotentialArray, long* seed);						//Loop that calls Forces for each bead, will inc potential

VEC FENEForce(VEC nPos, VEC nPosPlusOne, CONSTANTS c);

int initialiseseed(int numseeds, long** seed, CONSTANTS c);

VEC Brownian(long* seed, CONSTANTS c, int tid);

double GenGaussRand (long* seed, CONSTANTS c, int tid);

VEC  potential(CONSTANTS c, VEC* PositionArrayOld, VEC* PotentialArray, int i);																//Lennard-Jones Potential

int writeVTF(CONSTANTS c, VEC* frames);																		//Writes all bead positions to file after arrays finished

int writeKnotAnalysis(CONSTANTS c, VEC* frames);

int finalise(CONSTANTS* c, VEC** PositionArrayOld, VEC** PositionArrayNew, VEC** frames, VEC** VECArray, VEC** BrownianArray, VEC** PotentialArray);			//frees memory

void die(const char* message, const int line, const char* file);

#endif
