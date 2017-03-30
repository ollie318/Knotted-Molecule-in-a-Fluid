#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

int main(){
    for(int v = 15; v <= 200; v += 5){
        char buffer1[200];
        sprintf(buffer1, "cd Output_Files/Flow_Velocity/Vel%d\nqsub KnotVel%d\ncd ../../..", v, v);
        system(buffer1);
    }

    for(double eta = 1E-3; eta <= 1; eta+=eta){
        char buffer2[200];
        sprintf(buffer2, "cd Output_Files/Flow_Viscosity/Vel_100_Vis_%.3lf\nqsub KnotVis_%.3lf\ncd ../../..", eta, eta);
        system(buffer2);
    }

    for(double eta = 1E-3; eta <= 1; eta+=eta){
        char buffer3[200];
        sprintf(buffer3, "cd Output_Files/Flow_Viscosity/Vel_200_Vis_%.3lf\nqsub KnotVis_%.3lf\ncd ../../..", eta, eta);
        system(buffer3);
    }

    for(int i = 60; i <= 150; i+=5){
        char buffer4[200];
        sprintf(buffer4, "cd Output_Files/Chain_Length/Bead_num_%d\nqsub Knotlength_%d\ncd ../../..", i, i);
        system(buffer4);
    }
}
