#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

int main(){
    for(int v = 15; v <= 200; v += 5){
        char buffer1[80], buffer2[80], buffer3[80];
        sprintf(buffer1, "cd Output_Files/Flow_Velocity/Vel%d\nqsub KnotVel%d\ncd ../../..", v, v);
        system(buffer1);

    }

}
