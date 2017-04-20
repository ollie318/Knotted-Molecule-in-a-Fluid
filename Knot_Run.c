#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

int main(){

    for(double eta = 1E-3; eta < 1E-1; eta+=1E-3){
        char buffer2[200];
        sprintf(buffer2, "cd Output_Files2/Flow_Viscosity/Vel_20_Vis_%.3lf\nqsub KnotVis_%.3lf\ncd ../../..", eta, eta);
        system(buffer2);
    }


}
