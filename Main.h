#ifndef MAIN_H
#define MAIN_H

typedef struct
{
    double xcoord, ycoord, zcoord;
} VEC;

typedef struct
{
    double BeadRadi, FluidViscos, h, T, D, FlowVel, H, m, Q_0, PipeRad;
    double N_k, N_ks, b_k, L_s, eta, MaxExtension;
    int N, maxIters;
} CONSTANTS;

//FUNCTION PROTOTYPES

int initialise(CONSTANTS* c, VEC** PositionArrayOld, VEC** PositionArrayNew, VEC** frames, VEC** VECArray, VEC** BrownianArray, VEC** PotentialArray, long** seed, const char* paramfile);			//Sets constants, allocates memory for array of pointers

int CalcKnotPos(CONSTANTS c, VEC* PositionArrayOld);

int updateFrames(CONSTANTS c, int CurrentFrame, VEC* frames, VEC* positions);								//Writes current set of positions into frames, frames written to file at end

int timestep(CONSTANTS c, VEC* PositionArrayOld, VEC* PositionArrayNew, VEC* VECArray, VEC* BrownianArray, VEC* PotentialArray, gsl_rng** seed);						//Loop that calls Forces for each bead, will inc potential

VEC FENEForce(VEC nPos, VEC nPosPlusOne, CONSTANTS c);

VEC Brownian(CONSTANTS c, gsl_rng* seednum);

double StokesFlow(CONSTANTS c, VEC OldPos);

VEC  potential(CONSTANTS c, VEC* PositionArrayOld, VEC* PotentialArray, int i);																//Lennard-Jones Potential

int writeVTF(CONSTANTS c, VEC* frames, char* buffer);																		//Writes all bead positions to file after arrays finished

int writeKnotAnalysis(CONSTANTS c, VEC* frames, char* buffer);

int finalise(CONSTANTS* c, VEC** PositionArrayOld, VEC** PositionArrayNew, VEC** frames, VEC** VECArray, VEC** BrownianArray, VEC** PotentialArray);			//frees memory

double GenGaussRand (CONSTANTS c, gsl_rng* seednum);

void die(const char* message, const int line, const char* file);

#endif
