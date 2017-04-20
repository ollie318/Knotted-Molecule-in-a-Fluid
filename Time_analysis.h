#ifndef TIME_ANALYSIS_MAIN_H
#define TIME_ANALYSIS_MAIN_H

double Timeanalysis(char *test)
{
    long file_size;

    FILE* Knot_File;
    char buffer1[250];
    sprintf(buffer1, "%s", test);
    Knot_File = fopen(buffer1, "r");

    fseek(Knot_File, 0, SEEK_END);
    file_size = ftell(Knot_File);
    fseek(Knot_File, 0, SEEK_SET);

    char *myArray;
    myArray = malloc(file_size + 1);    //Extra byte to zero-terminate

    if(NULL == Knot_File)
    {
        printf("Null file\n");
        return 0;
    }

    fread(myArray, 1, file_size, Knot_File);
    myArray[file_size] = '\0';

    long line_num = 0;       //Num of lines passed, e.g. to get to line 5 set while(line_num != 4) then lines passed is 4, line is 5

    double time_unknot;
    long loopcount;
    long file_sizeMINUS_50 = file_size - 50;
    for(loopcount = 0; loopcount < file_sizeMINUS_50;)
    {
        if(myArray[loopcount] == '\n' )
        {
            line_num++;
        }

        loopcount++;
        if(myArray[loopcount] == 115 && loopcount>1000){
            time_unknot = line_num * 0.00005;
            break;
        }
    }

    free(myArray);

    return time_unknot;
}


#endif
