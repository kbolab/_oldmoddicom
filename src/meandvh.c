/* Function for calculatinge mean and median dvh */
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <time.h>


void meanmediandvh(double *dvh, int *Nrow, int *Nboot, double *meanDvh, double *sampledvh, int *Nhist) {
	int n, m, p; /* counters */
  int dvhindex;  /* index of histograms */  
  /* initialize random seed */
  srand ( time(NULL) );  
	for (n=0; n < *Nboot; n++) { /* iteration among bootstrap dvh */
    /* creation of resampled histogram series */    
		for (m=0; m < *Nhist; m++) {  /* iteration among number of histograms */
      dvhindex = rand() % *Nhist; /* create the index of the sampled dvh */
      for (p=0; p < *Nrow; p++) { /* start creating sampled dvh */
        sampledvh[m * *Nrow + p] = dvh[dvhindex * *Nrow + p]; /* set element in the dvh */        
      }      
		}  /* once created the sample series of dvhs it's time for calculating mean and median */
    for (m=0; m < *Nhist; m++) {
      for (p=0; p < *Nrow; p++) {
        meanDvh[n * *Nrow + p] += sampledvh[m * *Nrow + p];
      }            
    }
    for (p=0; p < *Nrow; p++) {
      meanDvh[n * *Nrow + p] = meanDvh[n * *Nrow + p] / *Nhist;
    }
	}  
}

