#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

int main()
{
    const int const_chain_length = 60; //+2 intermediate lines
    const int const_line_repeat = 62;
    long file_size;

    FILE* Coord_File;
    Coord_File = fopen("BeadPos_Sat_Apr__1_21_59_13_2017_.vtf","r");

    FILE* Results_File;
    Results_File = fopen("Lengths_Apr_1_21_59.txt", "w");

    fseek(Coord_File, 0, SEEK_END);
    file_size = ftell(Coord_File);
    fseek(Coord_File, 0, SEEK_SET);

    char *myArray;
    myArray = malloc(file_size + 1);    //Extra byte to zero-terminate

    if((NULL == Coord_File) || (NULL == Results_File))
    {
        printf("Null file\n");
        return 0;
    }

    fread(myArray, 1, file_size, Coord_File);
    myArray[file_size] = '\0';
    fclose(Coord_File);

    double dbl_x, dbl_y, dbl_z;                         //When converted to doubles
    double end_to_end_length;                           //Will hold sqrt(x^2 + y^2 + z^2), end-to-end-length
    char buff;                                          //Holds
    char *ptrChar = &buff;

    int line_num = 0;       //Num of lines passed, e.g. to get to line 5 set while(line_num != 4) then lines passed is 4, line is 5
    long loopcount = 0;      //Num of character

    int chain_start = 62;
    int chain_end = chain_start + const_chain_length - 1;     //last bead in chain

    for(loopcount = 0; loopcount < (file_size-50);)
    {
        while(line_num != chain_end)
        {
            if(myArray[loopcount] == '\n')
            {
                line_num++;
            }
            loopcount++;
        }

        dbl_x = strtod(&myArray[loopcount], &ptrChar);
        dbl_y = strtod(ptrChar+1, &ptrChar);            //ptrChar points to first char not in dbl_x
        dbl_z = strtod(ptrChar+1, NULL);
        end_to_end_length = sqrt(dbl_x*dbl_x + dbl_y*dbl_y + dbl_z*dbl_z);

        fprintf(Results_File, "%.13f", end_to_end_length);
        fprintf(Results_File, "\n");

        chain_end += const_line_repeat;
    }

    fclose(Results_File);
    free(myArray);

    return 0;
}
