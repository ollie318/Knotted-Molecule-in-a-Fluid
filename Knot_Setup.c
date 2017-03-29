#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>


#define NUMBEADS 60
#define NUMKNOTS 1
#define KNOTTYPE "8_19"

int main()
{

    FILE *Params;
    FILE *Script;
    char root[120], params[120], script[120];
    int v,k;
    double eta, j;

    mkdir("Output_Files", 0777);
    mkdir("Output_Files/Flow_Viscosity", 0777);
    mkdir("Output_Files/Flow_Velocity", 0777);
    mkdir("Output_Files/Chain_Length", 0777);

    for(v = 15; v <= 200; v += 5){
        sprintf(params, "Output_Files/Flow_Velocity/Vel%d", v);
        mkdir(params, 0777);

        sprintf(params, "Output_Files/Flow_Velocity/Vel%d/8_19Params.txt", v);
        Params = fopen(params, "w");

        fprintf(Params, "60\n");
        fprintf(Params, "10000000\n");
        fprintf(Params, "24E-9\n");
        fprintf(Params, "1E-2\n");
        fprintf(Params, "%d\n", v);
        fprintf(Params, "0.0000001\n");
        fprintf(Params, "298\n");
        fprintf(Params, "2626\n");
        fprintf(Params, "1.8E-9\n");

        fclose(Params);

        sprintf(script, "Output_Files/Flow_Velocity/Vel%d/KnotVel%d", v, v);
        Script = fopen(script, "w");

        fprintf(Script,"#!/bin/bash\n");
        fprintf(Script,"#PBS -l nodes=1:ppn=4,walltime=04:00:00\n");
        fprintf(Script,"cd $PBS_O_WORKDIR\n");
        fprintf(Script,"../../../a.out 8_19Params.txt");

        fclose(Script);
    }

    for(eta = 1E-3; eta < 1; eta+=eta){
        sprintf(params, "Output_Files/Flow_Viscosity/Vel_200_Vis_%.3lf", eta);
        mkdir(params, 0777);

        sprintf(params, "Output_Files/Flow_Viscosity/Vel_200_Vis_%.3lf/8_19Params.txt", eta);
        Params = fopen(params, "w");

        fprintf(Params, "60\n");
        fprintf(Params, "10000000\n");
        fprintf(Params, "24E-9\n");
        fprintf(Params, "%lf\n", eta);
        fprintf(Params, "200\n");
        fprintf(Params, "0.0000001\n");
        fprintf(Params, "298\n");
        fprintf(Params, "2626\n");
        fprintf(Params, "1.8E-9\n");

        fclose(Params);

        sprintf(script, "Output_Files/Flow_Viscosity/Vel_200_Vis_%.3lf/KnotVis_%.3lf", eta, eta);
        Script = fopen(script, "w");

        fprintf(Script,"#!/bin/bash\n");
        fprintf(Script,"#PBS -l nodes=1:ppn=4,walltime=04:00:00\n");
        fprintf(Script,"cd $PBS_O_WORKDIR\n");
        fprintf(Script,"../../../a.out 8_19Params.txt");

        fclose(Script);

    }

    for(eta = 1E-3; eta < 1; eta+=eta){
        sprintf(params, "Output_Files/Flow_Viscosity/Vel_100_Vis_%.3lf", eta);
        mkdir(params, 0777);

        sprintf(params, "Output_Files/Flow_Viscosity/Vel_100_Vis_%.3lf/8_19Params.txt", eta);
        Params = fopen(params, "w");

        fprintf(Params, "60\n");
        fprintf(Params, "10000000\n");
        fprintf(Params, "24E-9\n");
        fprintf(Params, "%lf\n", eta);
        fprintf(Params, "100\n");
        fprintf(Params, "0.0000001\n");
        fprintf(Params, "298\n");
        fprintf(Params, "2626\n");
        fprintf(Params, "1.8E-9\n");

        fclose(Params);

        sprintf(script, "Output_Files/Flow_Viscosity/Vel_100_Vis_%.3lf/KnotVis_%.3lf", eta, eta);
        Script = fopen(script, "w");

        fprintf(Script,"#!/bin/bash\n");
        fprintf(Script,"#PBS -l nodes=1:ppn=4,walltime=04:00:00\n");
        fprintf(Script,"cd $PBS_O_WORKDIR\n");
        fprintf(Script,"../../../a.out 8_19Params.txt");

        fclose(Script);

    }

    for(int i = 60; i <= 150; i+=5){
        sprintf(params, "Output_Files/Chain_Length/Bead_num_%d", i);
        mkdir(params, 0777);

        sprintf(params, "Output_Files/Chain_Length/Bead_num_%d/8_19Params.txt", i);
        Params = fopen(params, "w");

        j = (2626/60) * i;
        k = (j+0.5);

        fprintf(Params, "%d\n", i);
        fprintf(Params, "10000000\n");
        fprintf(Params, "24E-9\n");
        fprintf(Params, "1E-3\n");
        fprintf(Params, "100\n");
        fprintf(Params, "0.0000001\n");
        fprintf(Params, "298\n");
        fprintf(Params, "%d\n", k);
        fprintf(Params, "1.8E-9\n");

        fclose(Params);

        sprintf(script, "Output_Files/Chain_Length/Bead_num_%d/Knotlength_%d", i, i);
        Script = fopen(script, "w");

        fprintf(Script,"#!/bin/bash\n");
        fprintf(Script,"#PBS -l nodes=1:ppn=4,walltime=04:00:00\n");
        fprintf(Script,"cd $PBS_O_WORKDIR\n");
        fprintf(Script,"../../../a.out 8_19Params.txt");

        fclose(Script);

    }

}
