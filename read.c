/*
 * RECON Copyright (C) 2018 Yan-Rong Li
 * A package for measuring spectral power and reconstructing time series in AGN.
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

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

    strcpy(tag[nt], "FreqLimit");
    addr[nt] = &parset.freq_limit;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "PSDType");
    addr[nt] = &parset.psd_type;
    id[nt++] = INT;

    strcpy(tag[nt], "PSDArg");
    addr[nt] = &parset.str_psd_arg;
    id[nt++] = STRING;

    strcpy(tag[nt], "Domain");
    addr[nt] = &parset.flag_domain;
    id[nt++] = INT;

    strcpy(tag[nt], "FlagEndMatch");
    addr[nt] = &parset.flag_endmatch;
    id[nt++] = INT;

    strcpy(tag[nt], "FlagWhiteNoise");
    addr[nt] = &parset.flag_whitenoise;
    id[nt++] = INT;

    strcpy(tag[nt], "FlagSaveOutput");
    addr[nt] = &parset.flag_saveoutput;
    id[nt++] = INT;

    strcpy(tag[nt], "PSDPeriodModel");
    addr[nt] = &parset.psdperiod_model;
    id[nt++] = INT;

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

    strcpy(tag[nt], "FErr");
    addr[nt] = &parset.ferr;
    id[nt++] = DOUBLE;
    
    char fname[200];
    sprintf(fname, "%s", parset.param_file);
    
    fparam = fopen(fname, "r");
    if(fparam == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }

    //default values
    parset.flag_endmatch = 0;
    parset.flag_whitenoise = 0;

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

    if(parset.V < 1.0)
    {
      printf("Incorrect V=%f.\n V should be larger than 1.\n", parset.V);
      exit(0);
    }

    if(parset.W < 1.0)
    {
      printf("Incorrect W=%f.\n V should be larger than 1.\n", parset.V);
      exit(0);
    }

    if(parset.freq_limit <=0.0)
    {
      printf("Incorrect freq_limit=%f\n freq_limit should be positive.\n", parset.freq_limit);
      exit(0);
    }

    if(parset.ferr <=0.0)
    {
      printf("#Incorrect ferr=%f\n. ferr should be positive.\n", parset.ferr);
      exit(0);
    }

    if(recon_flag_psd == 1)
    {
      parset.flag_endmatch = 0;
    }
  }
  
  MPI_Bcast(&parset, sizeof(parset), MPI_BYTE, roottask, MPI_COMM_WORLD);
  return;
}

int read_data(char *fname, int n, double *t, double *f, double *e)
{
  FILE *fp;
  int i;
  char buf[200];
  double slope;

  fp = fopen(fname, "r");
  if(fp==NULL)
  {
    printf("Cannot open file %s.\n", fname);
    exit(0);
  }

  for(i=0; i<n; i++)
  {
    fgets(buf, 200, fp);
    sscanf(buf, "%lf %lf %lf\n", &t[i], &f[i], &e[i]);
  }
  fclose(fp);

  // end matching 
  if(parset.flag_endmatch == 1)
  {
    slope = (f[n-1] - f[0])/(t[n-1] - t[0]);
    for(i=0; i<n; i++)
    {
      f[i] -= (slope*(t[i] - t[0]));
    }
  }
  
  return 0;
}


void read_sim_arg()
{
  switch(parset.psd_type)
  {
    case 0: // single power-law
      sscanf(parset.str_psd_arg, "%lf:%lf:%lf", &parset.psd_arg[0], &parset.psd_arg[1], &parset.psd_arg[2]);
        
      if(parset.psd_arg[0] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 1st PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[0] == 0.0)
      {
        parset.psd_arg[0] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[0] = log(parset.psd_arg[0]);
      }
        
      if(parset.psd_arg[2] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 3rd PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[2] == 0.0)
      {
        parset.psd_arg[2] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[2] = log(parset.psd_arg[2]);
      }
    
    break;
  
    case 1: // damped random walk
      sscanf(parset.str_psd_arg, "%lf:%lf:%lf", &parset.psd_arg[0], &parset.psd_arg[1], &parset.psd_arg[2]);
 
      if(parset.psd_arg[0] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 1st PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[0] == 0.0)
      {
        parset.psd_arg[0] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[0] = log(parset.psd_arg[0]);
      }

      if(parset.psd_arg[1] <=0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 2nd PSDArg.\n");
          exit(0);
        }
      }
      else
      {
        parset.psd_arg[1] = log(parset.psd_arg[1]);
      }

      if(parset.psd_arg[2] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 3rd PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[2] == 0.0)
      {
        parset.psd_arg[2] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[2] = log(parset.psd_arg[2]);
      }

      break;
      
    case 2: // bending power-law
      sscanf(parset.str_psd_arg, "%lf:%lf:%lf:%lf:%lf", &parset.psd_arg[0], &parset.psd_arg[1], &parset.psd_arg[2],
                                                            &parset.psd_arg[3], &parset.psd_arg[4]);
      if(parset.psd_arg[0] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 1st PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[0] == 0.0)
      {
        parset.psd_arg[0] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[0] = log(parset.psd_arg[0]);
      }

      if(parset.psd_arg[3] <=0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 4th PSDArg.\n");
          exit(0);
        }
      }
      else
      {
        parset.psd_arg[3] = log(parset.psd_arg[3]);
      }

      if(parset.psd_arg[4] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 5th PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[4] == 0.0)
      {
        parset.psd_arg[4] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[4] = log(parset.psd_arg[4]);
      }

      break;

    case 3: // single power-law and periodic
      sscanf(parset.str_psd_arg, "%lf:%lf:%lf:%lf:%lf:%lf", &parset.psd_arg[0], &parset.psd_arg[1], &parset.psd_arg[2],
                &parset.psd_arg[3], &parset.psd_arg[4], &parset.psd_arg[5]);

      if(parset.psd_arg[0] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 1st PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[0] == 0.0)
      {
        parset.psd_arg[0] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[0] = log(parset.psd_arg[0]);
      }

      if(parset.psd_arg[2] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 3rd PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[2] == 0.0)
      {
        parset.psd_arg[2] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[2] = log(parset.psd_arg[2]);
      }

      if(parset.psd_arg[3] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 4th PSDArg.\n");
          exit(0);
         }
      }
      else if(parset.psd_arg[3] == 0.0)
      {
        parset.psd_arg[3] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[3] = log(parset.psd_arg[3]);
      }
        
      if(parset.psd_arg[4] <= 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 5th PSDArg.\n");
           exit(0);
        }
      }
      else
      {
        parset.psd_arg[4] = log(parset.psd_arg[4]);
      }

      if(parset.psd_arg[5] <= 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 5th PSDArg.\n");
           exit(0);
        }
      }
      else
      {
        parset.psd_arg[5] = log(parset.psd_arg[5]);
      }
      break;
    
    case 4:  // drw + periodic
      sscanf(parset.str_psd_arg, "%lf:%lf:%lf:%lf:%lf:%lf", &parset.psd_arg[0], &parset.psd_arg[1], &parset.psd_arg[2], 
                                                &parset.psd_arg[3], &parset.psd_arg[4], &parset.psd_arg[5]);

      if(parset.psd_arg[0] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 1st PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[0] == 0.0)
      {
        parset.psd_arg[0] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[0] = log(parset.psd_arg[0]);
      }

      if(parset.psd_arg[1] <=0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 2nd PSDArg.\n");
          exit(0);
        }
      }
      else
      {
        parset.psd_arg[1] = log(parset.psd_arg[1]);
      }

      if(parset.psd_arg[2] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 3rd PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[2] == 0.0)
      {
        parset.psd_arg[2] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[2] = log(parset.psd_arg[2]);
      }

      if(parset.psd_arg[3] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 4th PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[3] == 0.0)
      {
        parset.psd_arg[3] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[3] = log(parset.psd_arg[3]);
      }
      
      if(parset.psd_arg[4] <= 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 5th PSDArg.\n");
          exit(0);
        }
      }
      else
      {
        parset.psd_arg[4] = log(parset.psd_arg[4]);
      }

      if(parset.psd_arg[5] <= 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 5th PSDArg.\n");
          exit(0);
        }
      }
      else
      {
        parset.psd_arg[5] = log(parset.psd_arg[5]);
      }
      

      break;

    case 5:   // bending power-law + periodic 
      sscanf(parset.str_psd_arg, "%lf:%lf:%lf:%lf:%lf:%lf:%lf:%lf", &parset.psd_arg[0], &parset.psd_arg[1], &parset.psd_arg[2],
              &parset.psd_arg[3], &parset.psd_arg[4], &parset.psd_arg[5], &parset.psd_arg[6], &parset.psd_arg[7]);

      if(parset.psd_arg[0] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 1st PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[0] == 0.0)
      {
        parset.psd_arg[0] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[0] = log(parset.psd_arg[0]);
      }

      if(parset.psd_arg[3] <=0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 4th PSDArg.\n");
          exit(0);
        }
      }
      else
      {
        parset.psd_arg[3] = log(parset.psd_arg[3]);
      }

      if(parset.psd_arg[4] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 5th PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[4] == 0.0)
      {
        parset.psd_arg[4] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[4] = log(parset.psd_arg[4]);
      }

      if(parset.psd_arg[5] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 6th PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[5] == 0.0)
      {
        parset.psd_arg[5] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[5] = log(parset.psd_arg[5]);
      }
      
      if(parset.psd_arg[6] <= 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 7th PSDArg.\n");
          exit(0);
        }
      }
      else
      {
        parset.psd_arg[6] = log(parset.psd_arg[6]);
      }

      if(parset.psd_arg[7] <= 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 8th PSDArg.\n");
          exit(0);
        }
      }
      else
      {
        parset.psd_arg[7] = log(parset.psd_arg[7]);
      }
      
      break;

    default:
      parset.psd_arg[0] = log(1.0e0);
      parset.psd_arg[1] = 1.5;
      parset.psd_arg[2] = -DBL_MAX;
      break;
  }
}
