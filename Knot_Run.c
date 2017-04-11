#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

int main(){
    for(int v = 15; v <= 200; v += 5){
        char buffer1[200];
        sprintf(buffer1, "cd Output_Files2/Flow_Velocity/Vel%d\nqsub KnotVel%d\ncd ../../..", v, v);
        system(buffer1);
    }

    for(double eta = 1E-3; eta < 1E-1; eta+=1E-3){
        char buffer2[200];
        sprintf(buffer2, "cd Output_Files2/Flow_Viscosity/Vel_100_Vis_%.3lf\nqsub KnotVis_%.3lf\ncd ../../..", eta, eta);
        system(buffer2);
    }


}
