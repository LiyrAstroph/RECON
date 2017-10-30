#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

/*!
 * read parameter set from parameter file.
 */
void read_parset()
{
  if(thistask == roottask)
  {
    #define MAXTAGS 300
    #define DOUBLE 1
    #define STRING 2
    #define INT 3

    FILE *fparam;
    int i, j, nt;
    char str[200], buf1[200], buf2[200], buf3[200];
    int id[MAXTAGS];
    void *addr[MAXTAGS];
    char tag[MAXTAGS][50];

    nt = 0;
    strcpy(tag[nt], "FileDir");
    addr[nt] = &parset.file_dir;
    id[nt++] = STRING;

    strcpy(tag[nt], "FileName");
    addr[nt] = &parset.file_name;
    id[nt++] = STRING;

    strcpy(tag[nt], "V");
    addr[nt] = &parset.V;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "W");
    addr[nt] = &parset.W;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "PSDType");
    addr[nt] = &parset.psd_type;
    id[nt++] = INT;

    strcpy(tag[nt], "PSDArg");
    addr[nt] = &parset.str_psd_arg;
    id[nt++] = STRING;

    strcpy(tag[nt], "FileSim");
    addr[nt] = &parset.file_sim;
    id[nt++] = STRING;

    strcpy(tag[nt], "ND_Sim");
    addr[nt] = &parset.nd_sim;
    id[nt++] = INT;

    strcpy(tag[nt], "DT");
    addr[nt] = &parset.DT;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "FBad");
    addr[nt] = &parset.fbad;
    id[nt++] = DOUBLE;
    
    char fname[200];
    sprintf(fname, "%s", parset.param_file);
    
    fparam = fopen(fname, "r");
    if(fparam == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }

    while(!feof(fparam))
    {
      sprintf(str,"empty");

      fgets(str, 200, fparam);
      if(sscanf(str, "%s%s%s", buf1, buf2, buf3)<2)
        continue;
      if(buf1[0]=='#')
        continue;
      for(i=0, j=-1; i<nt; i++)
        if(strcmp(buf1, tag[i]) == 0)
        {
          j = i;
          tag[i][0] = 0;
          //printf("%s %s\n", buf1, buf2);
          break;
        }
      if(j >=0)
      {
        switch(id[j])
        {
          case DOUBLE:
            *((double *) addr[j]) = atof(buf2);
            break;
          case STRING:
            strcpy(addr[j], buf2);
            break;
          case INT:
            *((int *)addr[j]) = (int) atof(buf2);
            break;
        }
      }
      else
      {
        fprintf(stderr, "# Error in file %s: Tag '%s' is not allowed or multiple defined.\n", 
                      parset.param_file, buf1);
        exit(0);
      }
    }
    fclose(fparam);
  }
  
  MPI_Bcast(&parset, sizeof(parset), MPI_BYTE, roottask, MPI_COMM_WORLD);
  return;
}
