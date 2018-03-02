#ifndef GLOBALS_H
#define GLOBALS_H

#define EPS 0.000001

#include <stdio.h>
#include <stdlib.h>

int current_frame_particle_number;  // number of particles
int box_type; //if the system in non-cubic or NPT, get box size info from a datafile
int FRAMES; // frames to read from input xmol file
int SAMPLEFREQ; // frequency at which to take frames from the xmol file
int max_particle_number; // The number of particles in the largest XYZ frame

extern int num_cluster_types;  // The number of items in the cluster names array
extern int cluster_size[]; // A list of the number of particles in each cluster type
extern char* cluster_names[];  // A list of strings of cluster names
extern int* do_cluster_list[];  // A list of pointers to the do_clusts variables
extern int* num_cluster_list[];  // A list of pointers to the nclusts variables
extern char** raw_list[];  // A list of pointers to the "s" raw storage variables
extern int*** cluster_list[];  // A list of pointers to the "hc" cluster storage variables


struct xyz_info {
    int total_frames;
    int data_width;
    long *num_particles;
    long *frame_offsets;
};


char fInputParamsName[1000];    // name of parameters file and r... coordinates file and memsize file
char *fXmolName, *fBoxSizeName; //Name of xyz file, name of file which contains info on box
int box_offsets[1000];    // Offsets of each line in the box file
double *x, *y, *z;  // positions in x y and z directions of N particles
int *particle_type; // particle type
double sidex, sidey, sidez, half_sidex, half_sidey, half_sidez;
double tiltxy,tiltxz,tiltyz;


double rcutAA,rcutAA2,rcutAB,rcutAB2,rcutBB,rcutBB2;    // diameters of AB and BB interactions for binary interactions
double fc;  // Voronoi adjustment parameter
int Vor;    // 0 use simple bond length method Get_Bonds(), 1 use Voronoi method Get_Bonds_With_Voronoi()
int PBCs;   // 0 do not impliment periodic boundary conditions, 1 implement periodic boundary conditions
int nB; // max number of bonds per particle
int USELIST;    // 0  do not use cell list, 1 use cell list
int doWriteBonds;   // write bonds files out

int doWriteClus;    // write out indices of each detected cluster
int doWriteRaw; // write raw_*** cluster xmol files out
int do11AcenXyz; // write centres of 11A
int do13AcenXyz; // write centres of 13A
int doWritePopPerFrame; // write pop_per_frame file

int initNoStatic;   // initial size of static cluster arrays
int incrStatic; // when full, increment static cluster arrays by this amount
int initNoClustPerPart; // initial size of clusters per part arrays
int incrClustPerPart;   // when full, increment cluster per part arrays by this amount

int PRINTINFO; // print running information about progress

int *cnb; // Current Number of Bonds for particles {1,...,N}
int **bNums;    // list of particles (indices j) bonded to particle at index i
double **bondlengths;   // length of bonds in the bond network and squared
int maxnb; // max number of bonds to one particle
int correctedBonds; // max number of bonds to one particle


int n_cells_x, n_cells_y, n_cells_z, n_cells_total;   // number of cells per box length, total number of cells
int *head;   // head of cell array
int *llist; // linked list array
int *map; // list of neighbouring cells for cell i
double inv_cell_len_x, inv_cell_len_y, inv_cell_len_z;

// Whether to perform analysis of this type of cluster
int dosp3, dosp3a, dosp3b, dosp3c;
int dosp4, dosp4a, dosp4b, dosp4c;
int dosp5, dosp5a, dosp5b, dosp5c;
int do6Z, do7K;
int do8A, do8B, do8K;
int do9A, do9B, do9K;
int do10A, do10B, do10K, do10W;
int do11A, do11B, do11C, do11E, do11F, do11W;
int do12A, do12B, do12D, do12E, do12K;
int do13A, do13B, do13K;
int doFCC, doHCP, doBCC9, doBCC15;

// number of clusters of particlar type in current frame
int nsp3a, nsp3b, nsp3c;
int nsp4a, nsp4b, nsp4c;
int nsp5a, nsp5b, nsp5c;
int n6Z, n7K;
int n8A, n8B, n8K;
int n9A, n9B, n9K;
int n10A, n10B, n10K, n10W;
int n11A, n11B, n11C, n11E, n11F, n11W;
int n12A, n12B, n12D, n12E, n12K;
int n13A, n13B, n13K;
int nFCC, nHCP, nBCC_9, nBCC_15;

// max size of cluster storage arrays in dimension i
int msp3a, msp3b, msp3c;
int msp4a, msp4b, msp4c;
int msp5a, msp5b, msp5c;
int m6Z, m7K;
int m8A, m8B, m8K;
int m9A, m9B, m9K;
int m10A, m10B, m10K, m10W;
int m11A, m11B, m11C, m11E, m11F, m11W;
int m12A, m12B, m12D, m12E, m12K;
int m13A, m13B, m13K;
int mFCC, mHCP, mBCC_9, mBCC_15;

// cluster storage arrays (index i denotes number/identifier of cluster, j lists particles in cluster)
int **hcsp3a, **hcsp3b, **hcsp3c;
int **hcsp4a, **hcsp4b, **hcsp4c;
int **hcsp5a, **hcsp5b, **hcsp5c;
int **hc6Z, **hc7K;
int **hc8A, **hc8B, **hc8K;
int **hc9A, **hc9B, **hc9K;
int **hc10A, **hc10B, **hc10K, **hc10W;
int **hc11A, **hc11B, **hc11C, **hc11E, **hc11F, **hc11W;
int **hc12A, **hc12B, **hc12D, **hc12E, **hc12K;
int **hc13A, **hc13B, **hc13K;
int **hcFCC, **hcHCP, **hcBCC_9, **hcBCC_15;

// mem lists the clusters of that type each particle is in, index i is the particle index, j is the cluster id
// nmem lists the number of clusters of that type each particle is in, index i is the number of particles
// mmem lists the width of mem, the maximum number of clusters of the specified type associated with a single particle
int **mem_sp3b, *nmem_sp3b, mmem_sp3b;
int **mem_sp3c, *nmem_sp3c, mmem_sp3c;
int **mem_sp4b, *nmem_sp4b, mmem_sp4b;
int **mem_sp4c, *nmem_sp4c, mmem_sp4c;
int **mem_sp5b, *nmem_sp5b, mmem_sp5b;
int **mem_sp5c, *nmem_sp5c, mmem_sp5c;

// Raw lists of particle identity, output to RAW_clust files and reset each frame
char *ssp3a, *ssp3b, *ssp3c;
char *ssp4a, *ssp4b, *ssp4c;
char *ssp5a, *ssp5b, *ssp5c;
char *s6Z, *s7K;
char *s8A, *s8B, *s8K;
char *s9A, *s9B, *s9K;
char *s10A, *s10B, *s10K, *s10W;
char *s11A, *s11B, *s11C, *s11E, *s11F, *s11W;
char *s12A, *s12B, *s12D, *s12E, *s12K;
char *s13A, *s13B, *s13K;
char *sFCC, *sHCP, *sBCC_9, *sBCC_15;

// Lists of particle population of each cluster type in each frame, index i is the frame number,
// index j is the cluster type
double **pop_per_frame;

// The average population of each cluster type over all frames, index i is cluster type
double *mean_pop_per_frame;

// Gross number of clusters of the specified type accumulated over all frames
int *gross_clusters;

// Net number of clusters of the specified type accumulated over all frames
int *net_clusters;

// Total number of clusters of the specified type accumulated over all frames
int *total_clusters;

// Variable used in the counting of net clusters, index i is particle number
int *a6, *a7, *a8, *a9, *a10, *a11, *a12, *a13;

#endif