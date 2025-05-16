#ifndef GLOBALS_H
#define GLOBALS_H

#define EPS 0.000001

#include <stdio.h>
#include <stdlib.h>

int box_type;                    //if the system in non-cubic or NPT, get box size info from a datafile
int frames_to_analyse;           // frames to read from input xmol file
int num_cluster_types;           // The number of items in the cluster names array

extern int cluster_size[];        // A list of the number of particles in each cluster type
extern char* cluster_names[];     // A list of strings of cluster names
extern int* do_cluster_list[];    // A list of pointers to the do_clusts variables
extern int* num_cluster_list[];   // A list of pointers to the nclusts variables
extern char** raw_list[];         // A list of pointers to the "s" raw storage variables
extern int*** cluster_list[];     // A list of pointers to the "hc" cluster storage variables
extern int* cluster_list_width[]; // A list of pointers to the "m" storage width variables


struct xyz_info {
    long total_frames;
    long max_frames;
    long *num_particles;
    long *frame_offsets;
};

char *fXmolName, *fBoxSizeName; //Name of xyz file, name of file which contains info on box
long box_offsets[1000];          // Offsets of each line in the box file
double *x, *y, *z;              // positions in x y and z directions of N particles
int *particle_type;             // species of particle, index is particle number

double rcutAA,rcutAA2,rcutAB,rcutAB2,rcutBB,rcutBB2;    // diameters of AB and BB interactions for binary interactions
double min_cutAA, min_cutAA2;
double fc;                  // Voronoi adjustment parameter
int use_voronoi_bonds;      // 0 use simple bond length method build_bond_network(), 1 use Voronoi method Get_Bonds_With_Voronoi()
int PBCs;                   // 0 do not impliment periodic boundary conditions, 1 implement periodic boundary conditions
int max_num_bonds;          // max number of bonds per particle
int use_cell_list;          // 0  do not use cell list, 1 use cell list
int analyse_all_clusters;   // 0 Read clusters to analyse from file, 1 analyse all clusters

int doWriteBonds;           // write bonds files out
int doWriteClus;            // write out indices of each detected cluster
int doWriteRaw;             // write raw_*** cluster xmol files out
int do11AcenXyz;            // write centres of 11A
int do13AcenXyz;            // write centres of 13A
int eleven_A_number;        // The location of the 11A cluster in the cluster list
int thirteen_A_number;
int doWritePopPerFrame;     // write pop_per_frame file
int doWriteXYZ;             // Write clusters as XYZ file

int incrStatic;             // when full, increment cluster arrays by this amount

// Lists of particle population of each cluster type in each frame, index i is the frame number,
// index j is the cluster type
double **pop_per_frame;

// The average population of each cluster type over all frames, index i is cluster type
double *mean_pop_per_frame;

// Gross number of particles in the specified cluster type accumulated over all frames
int *num_gross_particles;

// Total number of clusters of the specified type accumulated over all frames
int *total_clusters;

// Per frame variables

// Box and bond variables

double sidex, sidey, sidez, half_sidex, half_sidey, half_sidez;
double tiltxy,tiltxz,tiltyz;
long particles_in_current_frame;

int *num_bonds;                                       // Current Number of Bonds for particles {1,...,N}
int **bond_list;                                      // list of particles (indices j) bonded to particle at index i
double **squared_bondlengths;                         // length of bonds in the bond network and squared
int maxnb;                                            // max number of bonds to one particle
int correctedBonds;                                   // bonds adjusted due to voronoi assymetry

int num_sort_columns;                                 // Number of columns to iterate over with quicksort

int n_cells_x, n_cells_y, n_cells_z, n_cells_total;   // number of cells per box length, total number of cells
int *head;                                            // head of cell array
int *linked_list;                                     // linked list array
double cell_len_x, cell_len_y, cell_len_z;

// Cluster variables

// whether to perform analysis of this type of cluster
int dosp3,   dosp3a,  dosp3b,   dosp3c;
int dosp4,   dosp4a,  dosp4b,   dosp4c;
int dosp5,   dosp5a,  dosp5b,   dosp5c;                                                     // 12 here but spN not in others so count 9
int do6A,    do6MW,   do6Z;                                                                 // 12
int do7K,    do7MW,   do7PAB,   do7T_a, do7T_s;                                             // 17
int do8A,    do8B,    do8K,     do8MW,  do8O,   do8PAA,  do8PAB,  do8PBB;                   // 25
int do9A,    do9B,    do9K,     do9MW,  do9PAA, do9PAB,  do9PBB,  do9S;                     // 33
int do10A,   do10B,   do10K,    do10MW, do10O,  do10PAA, do10PAB, do10PBB, do10S,   do10W;  // 43
int do11A,   do11B,   do11C,    do11E,  do11F,  do11MW,  do11O,   do11PAA, do11PAB;         // 52
int do11PBB, do11S,   do11SB,   do11W;                                                      // 56
int do12A,   do12B,   do12D,    do12E,  do12K,  do12MW,  do12O;                             // 63
int do12PAA, do12PAB, do12PBB,  do12S,  do12SB;                                             // 68
int do13A,   do13B,   do13K,    do13MW;                                                     // 72
int do13PAA, do13PAB, do13PBB,  do13S,  do13SB;                                             // 77
int do14O,   doFCC,   doHCP,    doBCC9;                                                     // 81

// number of clusters of particlar type in current frame
int nsp3a,  nsp3b,  nsp3c;
int nsp4a,  nsp4b,  nsp4c;
int nsp5a,  nsp5b,  nsp5c;                                                       // 9
int n6A,    n6MW,   n6Z;                                                         // 12
int n7K,    n7MW,   n7PAB,  n7T_a, n7T_s;                                        // 17
int n8A,    n8B,    n8K,    n8MW,  n8O,   n8PAA,  n8PAB,  n8PBB;                 // 25
int n9A,    n9B,    n9K,    n9MW,  n9PAA, n9PAB,  n9PBB,  n9S;                   // 33
int n10A,   n10B,   n10K,   n10MW, n10O,  n10PAA, n10PAB, n10PBB, n10S, n10W;    // 43
int n11A,   n11B,   n11C,   n11E,  n11F,  n11MW,  n11O,   n11PAA, n11PAB;        // 52
int n11PBB, n11S,   n11SB,  n11W;                                                // 56
int n12A,   n12B,   n12D,   n12E,  n12K,  n12MW,  n12O;                          // 63
int n12PAA, n12PAB, n12PBB, n12S,  n12SB;                                        // 68
int n13A,   n13B,   n13K,   n13MW;                                               // 72
int n13PAA, n13PAB, n13PBB, n13S,  n13SB;                                        // 77
int n14O, nFCC,   nHCP,   nBCC_9;                                                // 81

// max size of cluster storage arrays in dimension i
int msp3a,  msp3b,  msp3c;
int msp4a,  msp4b,  msp4c;
int msp5a,  msp5b,  msp5c;                                                       // 9
int m6A,    m6MW,   m6Z;                                                         // 12
int m7K,    m7MW,   m7PAB,  m7T_a, m7T_s;                                        // 17
int m8A,    m8B,    m8K,    m8MW,  m8O,   m8PAA,  m8PAB,  m8PBB;                 // 25
int m9A,    m9B,    m9K,    m9MW,  m9PAA, m9PAB,  m9PBB,  m9S;                   // 33
int m10A,   m10B,   m10K,   m10MW, m10O,  m10PAA, m10PAB, m10PBB, m10S,   m10W;  // 43
int m11A,   m11B,   m11C,   m11E,  m11F,  m11MW,  m11O,   m11PAA, m11PAB;        // 52
int m11PBB, m11S,   m11SB,  m11W;                                                // 56
int m12A,   m12B,   m12D,   m12E,  m12K,  m12MW, m12O;                           // 63
int m12PAA, m12PAB, m12PBB, m12S,  m12SB;                                        // 68
int m13A,   m13B,   m13K,   m13MW;                                               // 72
int m13PAA, m13PAB, m13PBB, m13S,  m13SB;                                        // 77
int m14O,   mFCC,   mHCP,   mBCC_9;                                              // 81

// cluster storage arrays (index i denotes number/identifier of cluster, j lists particles in cluster)
int **hcsp3a,  **hcsp3b,  **hcsp3c;
int **hcsp4a,  **hcsp4b,  **hcsp4c;
int **hcsp5a,  **hcsp5b,  **hcsp5c;                                                                          // 9
int **hc6A,    **hc6MW,   **hc6Z;                                                                            // 12
int **hc7K,    **hc7MW,   **hc7PAB,  **hc7T_a, **hc7T_s;                                                     // 17
int **hc8A,    **hc8B,    **hc8K,    **hc8MW,  **hc8O,   **hc8PAA,  **hc8PAB,  **hc8PBB;                     // 25
int **hc9A,    **hc9B,    **hc9K,    **hc9MW,  **hc9PAA, **hc9PAB,  **hc9PBB,  **hc9S;                       // 33
int **hc10A,   **hc10B,   **hc10K,   **hc10MW, **hc10O,  **hc10PAA, **hc10PAB, **hc10PBB, **hc10S, **hc10W;  // 43
int **hc11A,   **hc11B,   **hc11C,   **hc11E,  **hc11F,  **hc11MW,  **hc11O,   **hc11PAA, **hc11PAB;         // 52
int **hc11PBB, **hc11S,   **hc11SB,  **hc11W;                                                                // 56
int **hc12A,   **hc12B,   **hc12D,   **hc12E,  **hc12K,  **hc12MW,  **hc12O;                                 // 63
int **hc12PAA, **hc12PAB, **hc12PBB, **hc12S,  **hc12SB;                                                     // 68                                    // 68
int **hc13A,   **hc13B,   **hc13K,   **hc13MW;                                                               // 72
int **hc13PAA, **hc13PAB, **hc13PBB, **hc13S,  **hc13SB;                                                     // 77                                    // 68
int **hc14O,   **hcFCC,   **hcHCP,   **hcBCC_9;                                                              // 81

// Raw lists of particle identity, output to RAW_clust files and reset each frame
char *ssp3a,  *ssp3b,  *ssp3c;
char *ssp4a,  *ssp4b,  *ssp4c;
char *ssp5a,  *ssp5b,  *ssp5c;                                                             // 9
char *s6A,    *s6MW,   *s6Z;                                                               // 12
char *s7K,    *s7MW,   *s7PAB,  *s7T_a, *s7T_s;                                            // 17
char *s8A,    *s8B,    *s8K,    *s8MW,  *s8O,   *s8PAA,  *s8PAB,  *s8PBB;                  // 25
char *s9A,    *s9B,    *s9K,    *s9MW,  *s9PAA, *s9PAB,  *s9PBB,  *s9S;                    // 33
char *s10A,   *s10B,   *s10K,   *s10MW, *s10O,  *s10PAA, *s10PAB, *s10PBB,  *s10S, *s10W;  // 43
char *s11A,   *s11B,   *s11C,   *s11E,  *s11F,  *s11MW,  *s11O,   **s11PAA, **s11PAB;      // 52
char *s11PBB, *s11S,   *s11SB,  *s11W;                                                     // 56
char *s12A,   *s12B,   *s12D,   *s12E,  *s12K,  *s12MW,  *s12O;                            // 63
char *s12PAA, *s12PAB, *s12PBB, *s12S,  *s12SB;                                            // 68
char *s13A,   *s13B,   *s13K,   *s13MW;                                                    // 72
char *s13PAA, *s13PAB, *s13PBB, *s13S,  *s13SB;                                            // 77
char *s14O,   *sFCC,   *sHCP,   *sBCC_9;                                                   // 81

// mem lists the clusters of that type each particle is in, index i is the particle index, j is the cluster id
// nmem lists the number of clusters of that type each particle is in, index i is the number of particles
// mmem lists the width of mem, the maximum number of clusters of the specified type associated with a single particle (the largest value in nmem)
int **mem_sp3b, *nmem_sp3b, mmem_sp3b;
int **mem_sp3c, *nmem_sp3c, mmem_sp3c;
int **mem_sp4b, *nmem_sp4b, mmem_sp4b;
int **mem_sp4c, *nmem_sp4c, mmem_sp4c;
int **mem_sp5b, *nmem_sp5b, mmem_sp5b;
int **mem_sp5c, *nmem_sp5c, mmem_sp5c;

#endif
