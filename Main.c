#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Boltzmann 1.38064852E-23
#define pi 3.1415926535897932
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
    double BeadRadi, FluidViscos, h, T, D, FlowVel, H, m, Q_0;
    double N_k, N_ks, b_k, L_s, eta;
    int N, maxIters;
} CONSTANTS;

//FUNCTION PROTOTYPES

int CalcNextBallPos(POSITION nMinusOnePos, POSITION* nPos, ANGLES LastAngles, CONSTANTS c);

double GenRandDouble(double minDoub, double maxDoub);                   //Note: not truly random, testing purposes only!

ANGLES CalcNextAngles(CONSTANTS c);

TWO_GAUSS BoxMullerTrans (CONSTANTS c, double input_1, double input_2);

FENE FENEForce(POSITION nMinusOnePos, POSITION nPos, POSITION nPosPlusOne, CONSTANTS c);

double DragForce (CONSTANTS c);

BROWNIAN Brownian(CONSTANTS c);

POSITION Forces(POSITION nMinusOnePos, POSITION nPosOld, POSITION nPosPlusOne, POSITION nPosNew, CONSTANTS c);		//Gives new pos from sum forces

POSITION ForcesLast(POSITION nMinusOnePos, POSITION nPosOld, POSITION nPosPlusOne, POSITION nPosNew, CONSTANTS c);


int initialise(CONSTANTS* c, POSITION** PositionArrayOld, POSITION** PositionArrayNew, POSITION** frames);			//Sets constants, allocates memory for array of pointers

int timestep(CONSTANTS c, POSITION* PositionArrayOld, POSITION* PositionArrayNew);									//Loop that calls Forces for each bead, will inc collision

int updateFrames(CONSTANTS c, int CurrentFrame, POSITION* frames, POSITION* positions);								//Writes current set of positions into frames

int writeValues(CONSTANTS c, POSITION* frames);																		//Writes all bead positions to file after arrays finished

int finalise(CONSTANTS* c, POSITION** PositionArrayOld, POSITION** PositionArrayNew, POSITION** frames);			//frees memory

int collision(CONSTANTS c, POSITION* PositionArrayNew);																//checks for overlap & nudges away if needed

//MAIN PROG

int main()
{
  CONSTANTS c;
  POSITION* PositionArrayOld;
  POSITION* PositionArrayNew;
  POSITION* frames;

  initialise(&c, &PositionArrayOld, &PositionArrayNew, &frames);

  int loopcount;
  for(loopcount = 0; loopcount < c.maxIters; loopcount++){

      timestep(c, PositionArrayOld, PositionArrayNew);
      updateFrames(c, loopcount, frames, PositionArrayNew);
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
      PositionArrayNew[i] = Forces(PositionArrayOld[i-1], PositionArrayOld[i], PositionArrayOld[i+1], PositionArrayNew[i], c);
  }

  PositionArrayNew[c.N-1] = ForcesLast(PositionArrayOld[c.N-2], PositionArrayOld[c.N-1], PositionArrayOld[c.N-1], PositionArrayNew[c.N-1], c);
  // collision(c, PositionArrayNew);

  return EXIT_SUCCESS;
}

int initialise(CONSTANTS* c, POSITION** PositionArrayOld, POSITION** PositionArrayNew, POSITION** frames)
{
  // array constants
  c->N = 20;
  c->maxIters = 10000;


  c->BeadRadi = 100E-9;                                                //diameter of polystyrene
  c->FluidViscos = 30;                                             //Fluid viscosity
  c->FlowVel = 13.0;                                                    //Fluid velocity, NOT relative velocity as needed for Stokes Law
  c->h = 0.001;                                                          //Time step?
  c->T = 298;                                                          //Temperature in Kelvin
  c->D = (Boltzmann * c->T) / (6 * pi * c->FluidViscos * c->BeadRadi);    //Diffusion coefficient
  c->eta = 6*pi*c->FluidViscos*c->BeadRadi;

  //Spring coefficient calcs
  c->N_k = 2626;
  c->b_k = 1.8E-9;
  c->N_ks = c->N_k/c->N;
  c->L_s = c->N_ks*c->b_k;
  c->H = (3*Boltzmann*c->T)/(c->L_s*c->b_k);                              //Taken from Simons paper, values for polystyrene not DNA

  //c.m = 1.9927E-26;
  c->m = 0.104 / AvogadroNum;                                          //Bead mass for styrene
  c->Q_0 = c->N_ks * c->b_k;


  (*PositionArrayOld) = (POSITION*) malloc(sizeof(POSITION) * c->N);
  (*PositionArrayNew) = (POSITION*) malloc(sizeof(POSITION) * c->N);
  (*frames) = (POSITION*) malloc(sizeof(POSITION) * c->N * c->maxIters);

  (*PositionArrayOld)[0].xPos = 0;                                         //Initialising pos1 to 0,0,0
  (*PositionArrayOld)[0].yPos = 0;
  (*PositionArrayOld)[0].zPos = 0;

  (*PositionArrayOld)[1].xPos = (*PositionArrayOld)[0].xPos + (0.8 * c->Q_0);
  (*PositionArrayOld)[1].yPos = (*PositionArrayOld)[0].yPos;
  (*PositionArrayOld)[1].zPos = (*PositionArrayOld)[0].zPos;

  ANGLES LastBondAngles;
  LastBondAngles.phi = 0;
  LastBondAngles.theta = pi/2;

  int i;
  for (i = 2; i < c->N; i++){
      CalcNextBallPos((*PositionArrayOld)[i-1], &((*PositionArrayOld)[i]), LastBondAngles, *c);
  }

  return 0;
}

int CalcNextBallPos(POSITION nMinusOnePos, POSITION* nPos, ANGLES LastAngles, CONSTANTS c)
{
  ANGLES nAngles = CalcNextAngles(c);

  nPos -> xPos = nMinusOnePos.xPos + ( c.Q_0 * (sin(LastAngles.theta + nAngles.theta) * cos(LastAngles.phi + nAngles.phi)) );
  nPos -> yPos = nMinusOnePos.yPos + ( c.Q_0 * (sin(LastAngles.theta + nAngles.theta) * sin(LastAngles.phi + nAngles.phi)) );
  nPos -> zPos = nMinusOnePos.zPos + ( c.Q_0 * (cos(LastAngles.theta + nAngles.theta)) );

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
//Qs are current spring length, non-equilibrium
  FENEForces.FENE_x1 = (c.H * Q_x1) / (1 - ( pow(Q_x1, 2) / pow(c.Q_0, 2)));                             //Q_x1 is bond length for x & right, x2 left, y in y-dir etc
  double Q_x2 = nPos.xPos - nMinusOnePos.xPos;
  FENEForces.FENE_x2 = (c.H * Q_x2) / (1 - ( pow(Q_x2, 2) / pow(c.Q_0, 2)));                           //Values usually between -10 & 10 -- -297 on one run?

  double Q_y1 = nPosPlusOne.yPos - nPos.yPos;
  FENEForces.FENE_y1 = (c.H * Q_y1) / (1 - (pow(Q_y1, 2) / pow (c.Q_0, 2)));

  double Q_y2 = nPos.yPos - nMinusOnePos.yPos;
  FENEForces.FENE_y2 = (c.H * Q_y2) / (1 - (pow(Q_y2, 2) / pow (c.Q_0, 2)));

  double Q_z1 = nPosPlusOne.zPos - nPos.zPos;
  FENEForces.FENE_z1 = (c.H * Q_z1) / (1 - (pow(Q_z1, 2) / pow (c.Q_0, 2)));

  double Q_z2 = nPos.zPos - nMinusOnePos.zPos;
  FENEForces.FENE_z2 = (c.H * Q_z2) / (1 - (pow(Q_z2, 2) / pow (c.Q_0, 2)));

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

  BrownianForces.BrownianForce_x = sqrt((6 * Boltzmann * c.T)/(6 * pi * c.FluidViscos * c.BeadRadi * c.h)) * GenRandDouble(-1, 1);            //random number from 1 to -1
  BrownianForces.BrownianForce_y = sqrt((6 * Boltzmann * c.T)/(6 * pi * c.FluidViscos * c.BeadRadi * c.h)) * GenRandDouble(-1, 1);			  //will need Gaussian dist number -1 to 1
  BrownianForces.BrownianForce_z = sqrt((6 * Boltzmann * c.T)/(6 * pi * c.FluidViscos * c.BeadRadi * c.h)) * GenRandDouble(-1, 1);            //On the scale E-2

  return BrownianForces;
}

POSITION Forces(POSITION nMinusOnePos, POSITION nPosOld, POSITION nPosPlusOne, POSITION nPosNew, CONSTANTS c){

  BROWNIAN BrownianForces = Brownian(c);
  FENE FENEForces = FENEForce(nMinusOnePos, nPosOld, nPosPlusOne, c);

  double StokeDragForce = DragForce(c);

  nPosNew.xPos = nPosOld.xPos + (c.h*((FENEForces.FENE_x1 - FENEForces.FENE_x2)/c.eta  + BrownianForces.BrownianForce_x));
  // printf("%.12lf\n", (c.eta * c.FlowVel));

  nPosNew.yPos = nPosOld.yPos + (c.h*((FENEForces.FENE_x1 - FENEForces.FENE_x2)/c.eta  + BrownianForces.BrownianForce_y));

  nPosNew.zPos = nPosOld.zPos + (c.h*((FENEForces.FENE_x1 - FENEForces.FENE_x2)/c.eta  + BrownianForces.BrownianForce_z));
  return nPosNew;
}

POSITION ForcesLast(POSITION nMinusOnePos, POSITION nPosOld, POSITION nPosPlusOne, POSITION nPosNew, CONSTANTS c){

  BROWNIAN BrownianForces = Brownian(c);
  FENE FENEForces = FENEForce(nMinusOnePos, nPosOld, nPosPlusOne, c);

  double StokeDragForce = DragForce(c);

  nPosNew.xPos = nPosOld.xPos + (c.h*((FENEForces.FENE_x1 - FENEForces.FENE_x2)/c.eta + BrownianForces.BrownianForce_x));            //x INCREASES by FENE acc (a = F/m) and Brownian acc

  nPosNew.yPos = nPosOld.yPos + (c.h*((FENEForces.FENE_x1 - FENEForces.FENE_x2)/c.eta + BrownianForces.BrownianForce_y));

  nPosNew.zPos = nPosOld.zPos + (c.h*((FENEForces.FENE_x1 - FENEForces.FENE_x2)/c.eta + BrownianForces.BrownianForce_z));

  return nPosNew;
}

int collision(CONSTANTS c, POSITION* PositionArrayNew){
  double separation;

  int i, j;
  for(i = 0; i < c.N; i++)
  {
    for(j = 0; j < c.N; j++)
    {
      if(i == j)
        separation = 0;	//separation is a double so will not hold exactly 0
      else
        separation = sqrt(pow(PositionArrayNew[i].xPos - PositionArrayNew[j].xPos, 2) + pow(PositionArrayNew[i].yPos - PositionArrayNew[j].yPos, 2) + pow(PositionArrayNew[i].zPos - PositionArrayNew[j].zPos, 2));

      if(separation <= 2*c.BeadRadi && i != j)
      {
        POSITION unit;
        unit.xPos = (PositionArrayNew[i].xPos - PositionArrayNew[j].xPos) / separation;
        unit.yPos = (PositionArrayNew[i].yPos - PositionArrayNew[j].yPos) / separation;
        unit.zPos = (PositionArrayNew[i].zPos - PositionArrayNew[j].zPos) / separation;

        double nudge = (separation - 2*c.BeadRadi) / 2.0;

        PositionArrayNew[i].xPos += nudge * unit.xPos;
        PositionArrayNew[i].yPos += nudge * unit.yPos;
        PositionArrayNew[i].zPos += nudge * unit.zPos;

        PositionArrayNew[j].xPos -= nudge * unit.xPos;
        PositionArrayNew[j].yPos -= nudge * unit.yPos;
        PositionArrayNew[j].zPos -= nudge * unit.zPos;

      }
    }
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
