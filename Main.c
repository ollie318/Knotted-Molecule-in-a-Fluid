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

double DragForce (double BeadRadi, double FluidViscos, double FlowVel);

//MAIN PROG

int main()
{
    int const N = 10;

    double const radius = 1.0;

    POSITION PositionArray[N];
    
    FILE *File_BeadPos;

    File_BeadPos = fopen("File_BeadPos.txt", "w");

    PositionArray[0].xPos = 0;                  //Initialising pos1 to 0,0,0
    PositionArray[0].yPos = 0;
    PositionArray[0].zPos = 0;

    fprintf(File_BeadPos, "%lf%c%lf%c%lf\n", PositionArray[0].xPos, 9, PositionArray[0].yPos, 9, PositionArray[0].zPos);
    
    PositionArray[1].xPos = PositionArray[0].xPos + radius;
    PositionArray[1].yPos = PositionArray[0].yPos + radius;
    PositionArray[1].zPos = PositionArray[0].zPos + radius;

    fprintf(File_BeadPos, "%lf%c%lf%c%lf\n", PositionArray[1].xPos, 9, PositionArray[1].yPos, 9, PositionArray[1].zPos);
    
    int loopcount;

    for(loopcount = 2; loopcount <= N; loopcount++)
    {
        PositionArray[loopcount] = CalcNextBallPos(PositionArray[loopcount-2], PositionArray[loopcount-1]);
        
        fprintf(File_BeadPos, "%lf%c%lf%c%lf\n", PositionArray[loopcount].xPos, 9, PositionArray[loopcount].yPos, 9, PositionArray[loopcount].zPos);
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


    ANGLES nAngles = CalcNextAngles();

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


ANGLES CalcNextAngles()
{
    ANGLES NewAngles;
    NewAngles.theta = GenRandDouble(-pi/4, pi/4);
    NewAngles.phi = GenRandDouble(-pi, pi);

    return NewAngles;
}


double FENEForce(POSITION nMinusTwoPos, POSITION nMinusOnePos, POSITION nPos)
{

    double FENEForce_x1, FENEForce_x2, FENEForce_y1, FENEForce_y2, FENEForce_z1, FENEForce_z2, Q1, Q2, Q_0, H;

    H = ;
    Q_0 = ;
    FENEForce_x1 = (H* (nPos.xPos - nMinusOnePos.xPos)) / (1 - pow(nPos.xPos - nMinusOnePos.xPos, 2) / Q_0);
    FENEForce_x2 = (H* (nMinusOnePos.xPos - nMinusTwoPos.xPos)) / (1 - pow(nMinusOnePos.xPos - nMinusTwoPos.xPos, 2) / pow(Q_0, 2);
    FENEForce_y1 = (H* (nPos.yPos - nMinusOnePos.yPos)) / (1 - pow(nPos.yPos - nMinusOnePos.yPos, 2) / Q_0);
    FENEForce_y2 = (H* (nMinusOnePos.yPos - nMinusTwoPos.yPos)) / (1 - pow(nMinusOnePos.yPos - nMinusTwoPos.yPos, 2) / pow(Q_0, 2);
    FENEForce_z1 = (H* (nPos.zPos - nMinusOnePos.zPos)) / (1 - pow(nPos.zPos - nMinusOnePos.zPos, 2) / Q_0);
    FENEForce_z2 = (H* (nMinusOnePos.zPos - nMinusTwoPos.zPos)) / (1 - pow(nMinusOnePos.zPos - nMinusTwoPos.zPos, 2) / pow(Q_0, 2);

}


double DragForce (double BeadRadi, double FluidViscos, double FlowVel)
{
    double StokeForce;

    StokeForce = - 6 * M_PI * FluidViscos * FlowVel * BeadRadi;

    return StokeForce;
}