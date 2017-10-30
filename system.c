/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file system.c
 *  \brief get time used by the program.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <dirent.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <signal.h>
#include <string.h>
#include <mpi.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

/*!
 * This function provides the clock time of the system.
 */
double second()
{
#ifdef WALLCLOCK
	return MPI_Wtime();
#else
	return ((double) clock()) /CLOCKS_PER_SEC;
#endif
}

/*!
 *  This function calculates the time difference.
 */
double timediff(double t0, double t1)
{
	double dt;
	dt = t1 -t0;
	if(dt<0)
	{
#ifdef WALLCLOCK
	dt =0;
#else   
	dt = t1 + pow(2, sizeof(clock_t)*8) / CLOCKS_PER_SEC - t0;  /* 1 byte = 8 bits */
#endif
  }
	return dt;
}

/*!
 *  This function converts clock time into hms.
 */
void get_hms(double dt, int *h, int *m, double *s)
{
	*h = (int) floor(dt/3600);
	*m = (int) floor((dt - (*h)*3600)/60);
  *s = dt - (*h)*3600 - (*m)*60;
  return;
}


int remove_restart_file()
{
  struct dirent *rdf;
  DIR *rd;
  char dstr[200];

  sprintf(dstr, "%s/data/", parset.file_dir);
  rd = opendir(dstr);

  if(rd == NULL)
  {
    printf("Cannot open directory %s.\n", dstr);
    exit(0);
  }

  while ((rdf = readdir(rd)) != NULL) 
  {
    if(strncmp(rdf->d_name, "restart_dnest.txt_", strlen("restart_dnest.txt_"))==0)
    {
      sprintf(dstr, "%s/data/%s", parset.file_dir, rdf->d_name);

      if(remove(dstr)!=0)
      {
        printf("Cannot remove file %s.\n", dstr);
      }
    }
  }
  closedir(rd);
  return 0;
}