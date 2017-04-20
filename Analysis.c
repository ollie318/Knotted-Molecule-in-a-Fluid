#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include "Length_analysis.h"
#include "Time_analysis.h"


int main ()
{
  DIR *dp;
  struct dirent *ep;
  char filename[70];
  int chainlength, iterations, vel;
  double beaddiam, vis;
  FILE *Params;
  int found = 1;

  dp = opendir ("./");
  if (dp != NULL)
    {
      while ((ep = readdir (dp))){
        strncpy(filename, ep->d_name, 70);

        if(filename[0] == 56){
            Params = fopen(filename, "r");
            if(Params == NULL)printf("File not opened\n" );
            fscanf(Params, "%d\n%d\n%lf\n%lf\n%d", &chainlength, &iterations, &beaddiam, &vis, &vel);
            fclose(Params);
            found = 0;
        }
    }
      (void) closedir (dp);
    }
  else
    perror ("Couldn't open the directory");

    // dp = opendir ("./");
    // if (dp != NULL)
    //   {
    //     while ((ep = readdir (dp))){
    //       strncpy(filename, ep->d_name, 70);
    //       int letter = filename[0];
    //       if(filename[0] == 66 && found == 0)Lengthanalysis(filename, chainlength);
    //     }
    //     (void) closedir (dp);
    //   }
    // else
    //   perror ("Couldn't open the directory");

    double x[20] = {0.0};
    double totx = 0.0;
    int file_num = 0;

    FILE* Results_File;
    Results_File = fopen("../UnknotTime.txt", "a");

      dp = opendir ("./");
      if (dp != NULL)
        {
          while ((ep = readdir (dp))){
            strncpy(filename, ep->d_name, 70);

            if(filename[0] == 75 && filename[4] == 65){
                x[file_num] = Timeanalysis(filename);
                totx+=x[file_num];
                file_num++;
            }
        }
          (void) closedir (dp);
        }
      else perror ("Couldn't open the directory");

    double mean = (double) totx/(file_num);
    double meandiffsq = 0.0;

    for(int i=0; i<file_num; i++){
        meandiffsq += pow(x[i] - mean, 2);
    }

    double sd = sqrt(meandiffsq/(file_num-1));
    fprintf(Results_File, "%d\t%lf\t%.5lf\t%.5lf\n", vel, vis, (double) totx/file_num, sd);
    fclose(Results_File);

  return 0;
}
