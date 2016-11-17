#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//STRUCTURE/s

typedef struct
{
    double xPos, yPos, zPos;
} POSITION;

typedef struct
{
    double theta, psi;
} ANGLES;

//FUNCTION PROTOTYPES

POSITION CalcNextBallPos(POSITION nMinusTwoPos, POSITION nMinusOnePos);

double GenRandDouble(double minDoub, double maxDoub);                   //Note: not truly random, testing purposes only!

ANGLES CalcNextAngles();

double FENEForce(POSITION nMinusTwoPos, POSITION nMinusOnePos, POSITION Pos);


//MAIN PROG

int main()
{
    int const N = 5;

    double const radius = 1.0;

    POSITION PositionArray[N];

    PositionArray[0].xPos = 0;					//Initialising pos1 to 0,0,0
    PositionArray[0].yPos = 0;
    PositionArray[0].zPos = 0;

    PositionArray[1].xPos = PositionArray[0].xPos + radius;
    PositionArray[1].yPos = PositionArray[0].yPos + radius;
    PositionArray[1].zPos = PositionArray[0].zPos + radius;

    int loopcount;

    for(loopcount = 2; loopcount <= N; loopcount++)
    {
        PositionArray[loopcount] = CalcNextBallPos(PositionArray[loopcount-2], PositionArray[loopcount-1]);
    }
    return 0;
}

//FUNCTIONS

POSITION CalcNextBallPos(POSITION nMinusTwoPos, POSITION nMinusOnePos)
{
    POSITION nPos;
    nPos.xPos = 2*nMinusOnePos.xPos - nMinusTwoPos.xPos;    //Minus contribution from angles, just use rad?
    nPos.yPos = 2*nMinusOnePos.yPos - nMinusTwoPos.yPos;
    nPos.zPos = 2*nMinusOnePos.zPos - nMinusTwoPos.zPos;


    /*
    ANGLES nAngles = CalcNextAngles();

    nPos.xPos = nMinusOnePos.xPos + sin(nAngles.theta) * cos(nAngles.phi);
    nPos.yPos = nMinusOnePos.yPos + sin(nAngles.theta) * sin(nAngles.phi);
    nPos.zPos = nMinusOnePos.zPos + cos(theta);
    */


    return nPos;
}


double GenRandDouble(double minDoub, double maxDoub)
{
    double randDoub;
    double fraction = rand() / RAND_MAX;
    randDoub = minDoub + (maxDoub - minDoub) * fraction;
    return randDoub;
}

ANGLES CalcNextAngles()
{
    ANGLES NewAngles;
    NewAngles.theta = GenRandDouble(-pi/4, pi/4);
    NewAngles.phi = GenRandDouble(-pi, pi);

    return NewAngles;
}

double FENEForce(POSITION nMinusTwoPos, POSITION nMinusOnePos, POSITION Pos)
{

    double FENEForcex1, FENEForcex2, FENEForcey1, FENEForcey2, FENEForcez1, FENEForcez2, Q1, Q2, Q_0, H;

    H = ;
    Q_0 = ;
    FENEForcex1 = (H*(Pos.xPos-nMinusOnePos.xPos))/(1 - pow(Pos.xPos-nMinusOnePos.xPos, 2)/Q_0);
    FENEForcex2 = (H*(nMinusOnePos.xPos-nMinusTwoPos.xPos))/(1 - pow(nMinusOnePos.xPos-nMinusTwoPos.xPos, 2)/pow(Q_0, 2);
    FENEForcey1 = (H*(Pos.yPos-nMinusOnePos.yPos))/(1 - pow(Pos.yPos-nMinusOnePos.yPos, 2)/Q_0);
    FENEForcey2 = (H*(nMinusOnePos.yPos-nMinusTwoPos.yPos))/(1 - pow(nMinusOnePos.yPos-nMinusTwoPos.yPos, 2)/pow(Q_0, 2);
    FENEForcez1 = (H*(Pos.zPos-nMinusOnePos.zPos))/(1 - pow(Pos.zPos-nMinusOnePos.zPos, 2)/Q_0);
    FENEForcez2 = (H*(nMinusOnePos.zPos-nMinusTwoPos.zPos))/(1 - pow(nMinusOnePos.zPos-nMinusTwoPos.zPos, 2)/pow(Q_0, 2);

}






























