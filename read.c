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
#include <ctype.h>
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
    int i, j, k, nt;
    char str[256], buf1[256], buf2[256], buf3[256];
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

    strcpy(tag[nt], "FlagDomain");
    addr[nt] = &parset.flag_domain;
    id[nt++] = INT;

    strcpy(tag[nt], "PSDModel");
    addr[nt] = &parset.psd_model;
    id[nt++] = STRING;

    strcpy(tag[nt], "PSDType");
    addr[nt] = &parset.psd_type_str;
    id[nt++] = STRING;
    
    strcpy(tag[nt], "PSDPeriodModel");
    addr[nt] = &parset.psdperiod_model;
    id[nt++] = STRING;

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
    
    char fname[256];
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
    parset.flag_saveoutput = 0;
    parset.psd_model_enum = simple;
    parset.psd_type = 0;
    parset.harmonic_term_num=1;
    parset.carma_p = 1;
    parset.carma_q = 0;
    parset.psdperiod_enum = none;
    parset.flag_domain = 0;

    while(!feof(fparam))
    {
      sprintf(str,"empty");

      fgets(str, 256, fparam);
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
            // convert to lower case
            for(k=0; k<strlen(buf2); k++)
            {
              buf2[k] = tolower(buf2[k]);
            }
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

    // check domain
    if(parset.flag_domain != 0)
    {
      printf("# Currently only support fequency domain.\n Set FlagDomain to be 0.\n");
      exit(0);
    }

    // check psd model
    if(strcmp(parset.psd_model, "simple") ==0 )
    {
      parset.psd_model_enum = simple;
      parset.psd_type = atoi(parset.psd_type_str);
    }
    else if(strcmp(parset.psd_model, "harmonic") ==0 )
    {
      parset.psd_model_enum = harmonic;
      parset.harmonic_term_num = atoi(parset.psd_type_str);
    }
    else if(strcmp(parset.psd_model, "carma") ==0 )
    {
      parset.psd_model_enum = carma;
      sscanf(parset.psd_type_str, "%d:%d", &parset.carma_p, &parset.carma_q);
    }
    else
    {
      printf("# Incorrect PSDModel=%s.\nPSDModel should [simple, harmonic, carma].\n", parset.psd_model);
      exit(0);
    }
    
    // check psd type
    if(parset.psd_model_enum == simple)
    {
      if(parset.psd_type > 2 || parset.psd_type <0)
      {
        printf("# Incorrect PSDModel=%d.\nPSDModel should lie in the range [0-2].\n", parset.psd_type);
        exit(0);
      }
    }
    else if(parset.psd_model_enum == harmonic)
    {
      if(parset.harmonic_term_num > 5 || parset.harmonic_term_num < 1)
      {
        printf("# Incorrect PSDModel=%d.\nPSDModel should lie in the range [1-5].\n", parset.harmonic_term_num);
        exit(0);
      }

      //periodic PSD only valid for "simple" case
      strcpy(parset.psdperiod_model, "none");
      printf("# Reset PSDPeriodModel to be none.\n");

    }
    else if(parset.psd_model_enum == carma)
    {
      if(parset.carma_p <= parset.carma_q)
      {
        printf("# CARMA order p should be larger than q.\n");
        exit(0);
      }
      if(parset.carma_p > 10 || parset.carma_p <1)
      {
        printf("# Incorrect CARMA order p=%d.\np should lie in the range [1-10].\n", parset.carma_p);
        exit(0);
      }
      if(parset.carma_q > 9 || parset.carma_p < 0)
      {
        printf("# Incorrect CARMA order q=%d.\np should lie in the range [0-9].\n", parset.carma_q);
        exit(0);
      }
      //periodic PSD only valid for "simple" case
      strcpy(parset.psdperiod_model, "none");
      printf("# Reset PSDPeriodModel to be none.\n");
    }
    else
    {
      //periodic PSD only valid for "simple" case
      strcpy(parset.psdperiod_model, "none");
      printf("# Reset PSDPeriodModel to be none.\n");
    }
    
    // check periodic PSD model
    if(strcmp(parset.psdperiod_model, "none") ==0 )
    {
      parset.psdperiod_enum = none;
    }
    else if(strcmp(parset.psdperiod_model, "delta") ==0 )
    {
      parset.psdperiod_enum = delta;
    }
    else if(strcmp(parset.psdperiod_model, "gaussian") == 0)
    {
      parset.psdperiod_enum = gaussian;
    }
    else if(strcmp(parset.psdperiod_model, "lorentzian") == 0)
    {
      parset.psdperiod_enum = lorentzian;
    }
    else
    {
      printf("Incorrect PSDPeriodModel=%s.\nPSDPeriodModel should [none, delta, gaussian, lorentzian].\n", parset.psdperiod_model);
      exit(0);
    }
    
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

    if(parset.ferr < 0.0)
    {
      printf("#Incorrect ferr=%f. FErr should not be negative.\n", parset.ferr);
      exit(0);
    }

    if(recon_flag_cal_psd == 1)
    {
      parset.slope_endmatch = 0.0;
      parset.flag_endmatch = 0;
    }

  }
  MPI_Bcast(&parset, sizeof(parset), MPI_BYTE, roottask, MPI_COMM_WORLD);

  return;
}
/*!
 * read data 
 */
int read_data(char *fname, int n, double *t, double *f, double *e)
{
  FILE *fp;
  int i;
  char buf[256];

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
  parset.slope_endmatch = 0.0; 
  if(parset.flag_endmatch == 1)
  {
    parset.slope_endmatch = (f[n-1] - f[0])/(t[n-1] - t[0]);
    for(i=0; i<n; i++)
    {
      f[i] -= (parset.slope_endmatch*(t[i] - t[0]));
    }
  }
  
  return 0;
}

/*
 * deal with input PSD argument.
 */
void read_sim_arg()
{
  int i;
  char *str, *pstr;
  str = parset.str_psd_arg;

  //read parameter values
  for(i=0; i<num_params_psd_tot; i++)
  {
    if(sscanf(str, "%lf", &parset.psd_arg[i]) < 1)
    {
      if(thistask == roottask)
      {
        printf("# PSDArg lacks inputs. Only %d values specified.\n", i);
      }
      exit(0);
    }
    
    if(thistask == roottask)
      printf("%d: %f\n", i, parset.psd_arg[i]);
    
    if(i >= num_params_psd_tot-1)
      break;

    pstr = strchr(str, ':');
    if(pstr == NULL)
    {
      if(thistask == roottask)
      {
        printf("# PSDArg lacks inputs. Only %d values specified.\n", i+1);
      }
      exit(0);
    }
    str = pstr+1;
  }
  
  if(parset.psd_model_enum == simple)
  {
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

    switch(parset.psd_type)
    {
      case 0: // single power-law               
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
   
      default:
        parset.psd_arg[0] = log(1.0e0);
        parset.psd_arg[1] = 1.5;
        parset.psd_arg[2] = -DBL_MAX;
        break;
    }
    
    if(parset.psdperiod_enum > delta)
    {      
      i = parset.num_params_psd+0;
      if(parset.psd_arg[i] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 4th PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[i] == 0.0)
      {
        parset.psd_arg[i] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[i] = log(parset.psd_arg[i]);
      }
       
      i = parset.num_params_psd+1;   
      if(parset.psd_arg[i] <= 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 5th PSDArg.\n");
          exit(0);
        }
      }
      else
      {
        parset.psd_arg[i] = log(parset.psd_arg[i]);
      }
   
      i = parset.num_params_psd+2;
      if(parset.psd_arg[i] <= 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 6th PSDArg.\n");
          exit(0);
        }
      }
      else
      {
        parset.psd_arg[i] = log(parset.psd_arg[i]);
      }
    }
    else if(parset.psdperiod_enum > none)
    {
      i = parset.num_params_psd+0;
      if(parset.psd_arg[i] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 4th PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[i] == 0.0)
      {
        parset.psd_arg[i] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[i] = log(parset.psd_arg[i]);
      }
       
      i = parset.num_params_psd+1;   
      if(parset.psd_arg[i] <= 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 5th PSDArg.\n");
          exit(0);
        }
      }
      else
      {
        parset.psd_arg[i] = log(parset.psd_arg[i]);
      }
   
      i = parset.num_params_psd+2;
      if(parset.psd_arg[i] < 0.0 || parset.psd_arg[i] > 1.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 6th PSDArg.\n");
          exit(0);
        }
      }
    }
  }

  // harmonic PSD
  if(parset.psd_model_enum == harmonic)
  {
    for(i=0; i<2+3*(parset.harmonic_term_num-1)+1; i++)
    {
      if(parset.psd_arg[i] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 1st PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[i] == 0.0)
      {
        parset.psd_arg[i] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[i] = log(parset.psd_arg[i]);
      }
    }
  }
 

  //carma PSD
  if(parset.psd_model_enum == carma)
  {
    for(i=0; i<parset.carma_p + parset.carma_q +1 +1; i++)
    {
      if(parset.psd_arg[i] < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect 1st PSDArg.\n");
          exit(0);
        }
      }
      else if(parset.psd_arg[i] == 0.0)
      {
        parset.psd_arg[i] = -DBL_MAX;
      }
      else
      {
        parset.psd_arg[i] = log(parset.psd_arg[i]);
      }
    }
  }

  return;
}


/*!
 * get number of particles from the option file.
 */
void get_num_particles(char *fname)
{
  FILE *fp;
  char buf[256], buf1[256];
  fp = fopen(fname, "r");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }

  buf[0]='#';
  while(buf[0]=='#')
  {
    fgets(buf, 256, fp);
    if(sscanf(buf, "%s", buf1) < 1)  // a blank line
    {
      buf[0] = '#';
    }
  }
  sscanf(buf, "%d", &num_particles);
  fclose(fp);
}

/*!
 * get maximum number of levels from the option file.
 */
void get_max_num_levels(char *fname)
{
  FILE *fp;
  char buf[256], buf1[256];
  fp = fopen(fname, "r");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }

  buf[0]='#';
  while(buf[0]=='#')
  {
    fgets(buf, 256, fp);
    if(sscanf(buf, "%s", buf1) < 1)  // a blank line
    {
      buf[0] = '#';
    }
  }
  fgets(buf, 256, fp);
  fgets(buf, 256, fp);
  fgets(buf, 256, fp);
  fgets(buf, 256, fp);
  sscanf(buf, "%d", &max_num_levels);
  fclose(fp);
}