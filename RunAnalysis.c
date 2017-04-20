#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

int main(){
    for(int v=15; v<=200; v+=5){
        char buffer1[200];
        sprintf(buffer1, "cd Flow_Velocity/Vel%d\n../../Analysis", v);
        system(buffer1);
    }

    for(double eta = 1E-3; eta<1; eta+=eta){
        char buffer2[200];
        sprintf(buffer2, "cd Flow_Viscosity/Vel_100_Vis_%.3lf\n../../Analysis", eta);
        system(buffer2);
    }

    for(int n=60; n<150; n+=5){
        char buffer2[200];
        sprintf(buffer2, "cd Chain_Length/Bead_num_%d\n../../Analysis", n);
        system(buffer2);
    }
}
