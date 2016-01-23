#ifndef GLOBALS_H
#define GLOBALS_H

#define EPS 0.000001


#include <stdio.h>
#include <stdlib.h>
//// START: Global variables
int rank,size;
int N, NA;  // number of particles
int doBinary;   // if two particle species, look at correlations between them
int ISNOTCUBIC; //if the system in non-cubic or NPT, get box size info from a datafile
int FRAMES; // frames to read from input xmol file
int STARTFROM;  // start reading from this frame in the xmol file
int SAMPLEFREQ; // frequency at which to take frames from the xmol file
int TOTALFRAMES;
double TSTART;
double FRAMETSTEP;
double TFINAL;

char fInputParamsName[1000], fDynMemSizeName[1000], fgsblName[1000], fPotentialParamsName[1000];    // name of parameters file and r... coordinates file and memsize file
char *fXmolName, *fBoxSizeName; //Name of xyz file, name of file which contains info on box
double *x, *y, *z;  // positions in x y and z directions of N particles
int *rtype; // particle type
double side, halfSide;  // box side length
double sidex, sidey, sidez, halfSidex, halfSidey, halfSidez; //NPT_FIX
double tilt;


double rcutAA,rcutAA2,rcutAB,rcutAB2,rcutBB,rcutBB2;    // diameters of AB and BB interactions for binary interactions
double fc;  // Voronoi adjustment parameter
int Vor;    // 0 use simple bond length method Bonds_GetBonds(), 1 use Voronoi method Bonds_GetBondsV()
int PBCs;   // 0 do not impliment periodic boundary conditions, 1 implement periodic boundary conditions
int nB; // max number of bonds per particle
int USELIST;    // 0  do not use cell list, 1 use cell list
int doWriteBonds;   // write bonds files out

int doWriteClus;    // write out indices of each detected cluster
int doWriteRaw; // write raw_*** cluster xmol files out
int do11AcenXmol; // write centres of 11A   
int do13AcenXmol; // write centres of 13A   
int doWritePopPerFrame; // write pop_per_frame file
double binWidth;    // width of bins for bond length histogram
int doBLDistros;    // calculate bond length distributions
int doClusBLDistros;    // calculate bond length distributions for each cluster type
int doClusBLDeviation;  // look at cluster bond lengths relative to ground state
int donbDistros;        // do distribution of number of bonds to a particle
int doBondedCen;    // look at cluster bond lengths relative to ground state
int doClusComp; // look at composition of clusters in terms of A and B species

int doPotential;    // do potential energy calculations
int WHICHPOTENTIAL; // which potential to use

int doCoslovich;    // do Coslovich voronoi faces analysis

int initNoStatic;   // initial size of static cluster arrays
int incrStatic; // when full, increment static cluster arrays by this amount
int initNoClustPerPart; // initial size of clusters per part arrays
int incrClustPerPart;   // when full, increment cluster per part arrays by this amount

int doDynamics; // do dynamics analysis for clusters
int initNoLifetimes;    // initial number of lifetimes for dynamic clusters
int initNoDynamicClusters;  // initial number of dynamic clusters
int incrDynamicClusters;    // when full, increment dynamic cluster arrays by this amount
int doSubClusts; // write sub clusters in tcc files
double talpha; // alpha relaxtion time
int PRINTINFO; // print running information about progress


int *cnb; // Current Number of Bonds for particles {1,...,N}
int **bNums;    // list of particles (indices j) bonded to particle at index i
double **bondlengths;   // length of bonds in the bond network and squared
int maxnb; // max number of bonds to one particle
int correctedBonds; // max number of bonds to one particle

int BLDistroNoBins; // number of bins in bond network histogram

int M, M_pot, ncells, ncells_pot;   // number of cells per box length, total number of cells
int *head, *head_pot;   // head of cell array
int *llist, *llist_pot; // linked list array
int *map, *map_pot; // list of neighbouring cells for cell i
double cellSide, invcellSide, cellSide_pot, invcellSide_pot;

double sigma_AB, sigma_BB, epsilon_AB, epsilon_BB;
double sigma_AB6, sigma_BB6, quarteruTail_AA, quarteruTail_AB, quarteruTail_BB;
double rcut, rcut_AB, rcut_BB;
double rcut2, rcut_AB2, rcut_BB2;
double stoddardford_AA, stoddardford_AB, stoddardford_BB;
double *psigma;
char fSigmaName[1000];
double RHO0, mepsilon, mcut;
double uMorseTail, mcut2;
double KAPPA, yukepsilon, yukcut;
double uYukTail, yukcut2;
double ipl_exp, ipl_pre;
double half_ipl_exp;
double sigma_AB2, sigma_BB2;
double uTail_AA, uTail_AB, uTail_BB;
double cubic_a_AA, cubic_a_AB, cubic_a_BB;
double cubic_a_AA2, cubic_a_AB2, cubic_a_BB2;
double cubicA_AA, cubicA_AB, cubicA_BB;
double cubicB_AA, cubicB_AB, cubicB_BB;

double TWOPARTRHO, NOPOTENTIALBINS, WRITEPOTMIN, WRITEPOTMAX;
double *part_pot, potential, av_potential, av_pot_check;

char *s_s_0_2_8, *s_s_1_2_5_3, *s_s_1_2_5_2, *s_s_0_3_6, *s_s_0_0_12;
int nCos_s_0_2_8, nCos_s_1_2_5_3, nCos_s_1_2_5_2, nCos_s_0_3_6, nCos_s_0_0_12;
int np_s_0_2_8, np_s_1_2_5_3, np_s_1_2_5_2, np_s_0_3_6, np_s_0_0_12;
double *pop_per_frame_s_0_2_8, *pop_per_frame_s_1_2_5_3, *pop_per_frame_s_1_2_5_2, *pop_per_frame_s_0_3_6, *pop_per_frame_s_0_0_12;

char *s_b_0_2_8_4, *s_b_0_2_8_5, *s_b_0_3_6_6, *s_b_0_1_10_4;
int nCos_b_0_2_8_4, nCos_b_0_2_8_5, nCos_b_0_3_6_6, nCos_b_0_1_10_4;
int np_b_0_2_8_4, np_b_0_2_8_5, np_b_0_3_6_6, np_b_0_1_10_4;
double *pop_per_frame_b_0_2_8_4, *pop_per_frame_b_0_2_8_5, *pop_per_frame_b_0_3_6_6, *pop_per_frame_b_0_1_10_4;

int dosp3, dosp3a, dosp3b, dosp3c;
int dosp4, dosp4a, dosp4b, dosp4c;
int dosp5, dosp5a, dosp5b, dosp5c;
int do6Z, do7K, do8A, do8B, do8K, do9A, do9B, do9K, do10A, do10B, do10K, do10W;
int do11A, do11B, do11C, do11E, do11F, do11W, do12A, do12B, do12D, do12E, do12K;
int do13A, do13B, do13K, doFCC, doHCP, doBCC9, doBCC15;

int BLDistroNoSamples, BLDistroNoSamplesAA, BLDistroNoSamplesAB, BLDistroNoSamplesBB;   // number of samples in bond length histogram
int BLDistroNoSamplessp3, BLDistroNoSamplessp3a, BLDistroNoSamplessp3b, BLDistroNoSamplessp3c;  // number of samples in bond length histogram
int BLDistroNoSamplessp4, BLDistroNoSamplessp4a, BLDistroNoSamplessp4b, BLDistroNoSamplessp4c;  // number of samples in bond length histogram
int BLDistroNoSamples6A;    // number of samples in bond length histogram
int BLDistroNoSamplessp5, BLDistroNoSamplessp5a, BLDistroNoSamplessp5b, BLDistroNoSamplessp5c;  // number of samples in bond length histogram
int BLDistroNoSamples6Z, BLDistroNoSamples7K;   // number of samples in bond length histogram
int BLDistroNoSamples8A, BLDistroNoSamples8B, BLDistroNoSamples8K;  // number of samples in bond length histogram   
int BLDistroNoSamples9A, BLDistroNoSamples9B, BLDistroNoSamples9K;  // number of samples in bond length histogram
int BLDistroNoSamples10A, BLDistroNoSamples10B, BLDistroNoSamples10K, BLDistroNoSamples10W; // number of samples in bond length histogram
int BLDistroNoSamples11A, BLDistroNoSamples11B, BLDistroNoSamples11C, BLDistroNoSamples11E, BLDistroNoSamples11F, BLDistroNoSamples11W; // number of samples in bond length histogram
int BLDistroNoSamples12A, BLDistroNoSamples12B, BLDistroNoSamples12D, BLDistroNoSamples12E, BLDistroNoSamples12K;   // number of samples in bond length histogram
int BLDistroNoSamples13A, BLDistroNoSamples13B, BLDistroNoSamples13K;   // number of samples in bond length histogram
int BLDistroNoSamplesFCC, BLDistroNoSamplesHCP, BLDistroNoSamplesBCC_9, BLDistroNoSamplesBCC_15;    // number of samples in bond length histogram

int *BLDistro, *BLDistroAA,*BLDistroAB,*BLDistroBB; // number of bond length histogram
int *BLDistrosp3, *BLDistrosp3a, *BLDistrosp3b, *BLDistrosp3c;  // number of bond length histogram
int *BLDistrosp4, *BLDistrosp4a, *BLDistrosp4b, *BLDistrosp4c;  // number of bond length histogram
int *BLDistro6A;    // number of bond length histogram
int *BLDistrosp5, *BLDistrosp5a, *BLDistrosp5b, *BLDistrosp5c;  // number of bond length histogram
int *BLDistro6Z, *BLDistro7K;   // number of bond length histogram
int *BLDistro8A, *BLDistro8B, *BLDistro8K;      // number of bond length histogram
int *BLDistro9A, *BLDistro9B, *BLDistro9K;  // number of bond length histogram
int *BLDistro10A, *BLDistro10B, *BLDistro10K, *BLDistro10W; // number of bond length histogram
int *BLDistro11A, *BLDistro11B, *BLDistro11C, *BLDistro11E, *BLDistro11F, *BLDistro11W; // number of bond length histogram
int *BLDistro12A, *BLDistro12B, *BLDistro12D, *BLDistro12E, *BLDistro12K;   // number of bond length histogram
int *BLDistro13A, *BLDistro13B, *BLDistro13K;   // number of bond length histogram
int *BLDistroFCC, *BLDistroHCP, *BLDistroBCC_9, *BLDistroBCC_15;    // number of bond length histogram

double meanBL, meanBLAA, meanBLAB, meanBLBB;
double meanBLsp3, meanBLsp3a, meanBLsp3b, meanBLsp3c;   // number of sp3a/b/c respectively
double meanBLsp4, meanBLsp4a, meanBLsp4b, meanBLsp4c;   // number of sp4a/b/c respectively
double meanBL6A;
double meanBLsp5, meanBLsp5a, meanBLsp5b, meanBLsp5c;   // number of sp5a/b/c respectively
double meanBL6Z, meanBL7K;  // number of clusters of particlar type
double meanBL8A, meanBL8B, meanBL8K;    
double meanBL9A, meanBL9B, meanBL9K;
double meanBL10A, meanBL10B, meanBL10K, meanBL10W;
double meanBL11A, meanBL11B, meanBL11C, meanBL11E, meanBL11F, meanBL11W;    // number of clusters of particlar type
double meanBL12A, meanBL12B, meanBL12D, meanBL12E, meanBL12K;
double meanBL13A, meanBL13B, meanBL13K;
double meanBLFCC, meanBLHCP, meanBLBCC_9, meanBLBCC_15;

int *nsp3, *nsp3a, *nsp3b, *nsp3c;  // number of sp3a/b/c respectively
int *nsp4, *nsp4a, *nsp4b, *nsp4c;  // number of sp4a/b/c respectively
int *n6A;
int *nsp5, *nsp5a, *nsp5b, *nsp5c;  // number of sp5a/b/c respectively
int *n6Z, *n7K; // number of clusters of particlar type
int *n8A, *n8B, *n8K;   
int *n9A, *n9B, *n9K;
int *n10A, *n10B, *n10K, *n10W;
int *n11A, *n11B, *n11C, *n11E, *n11F, *n11W;   // number of clusters of particlar type
int *n12A, *n12B, *n12D, *n12E, *n12K;
int *n13A, *n13B, *n13K;
int *nFCC, *nHCP, *nBCC_9, *nBCC_15;


int nbDistroNoSamples, nbDistroNoSamplesAA, nbDistroNoSamplesAB, nbDistroNoSamplesBA, nbDistroNoSamplesBB; // max number of bonds to one particle
int *nbDistro, *nbDistroAA, *nbDistroAB, *nbDistroBA, *nbDistroBB; // max number of bonds to one particle
double meannb, meannbAA, meannbAB, meannbBA, meannbBB;
int maxto3, maxto4, maxto5; // max number of bonds to one particle
int totNclus; // total number of particles identified in clusters

int *nsp3c_spindlebonds, *nsp4c_spindlebonds, *n6A_spindlebonds, *nsp5c_spindlebonds;
int *nsp3_excess_spindles, *nsp4_excess_spindles, *nsp5_excess_spindles;    // number of sp3a/b/c _excess_spindlesed basic clusters

int msp3, msp3a, msp3b, msp3c;  // max size of **sp** arrays in dimension i
int msp4, msp4a, msp4b, msp4c, m6A; // max size of **sp** arrays in dimension i
int msp5, msp5a, msp5b, msp5c;  // max size of **sp** arrays in dimension i
int m6Z, m7K;   // max size of m** arrays in dimension i
int m8A, m8B, m8K;  // max size of m** arrays in dimension i
int m9A, m9B, m9K;  // max size of m** arrays in dimension i
int m10A, m10B, m10K, m10W; // max size of m** arrays in dimension i
int m11A, m11B, m11C, m11E, m11F, m11W; // max size of m** arrays in dimension i
int m12A, m12B, m12D, m12E, m12K;   // max size of m** arrays in dimension i
int m13A, m13B, m13K;   // max size of m** arrays in dimension i
int mFCC, mHCP, mBCC_9, mBCC_15;    // max size of **sp** arrays in dimension i

double gsbl_sp3, gsbl_sp3a, gsbl_sp3b, gsbl_sp3c;   // gsbl_ax size of **sp** arrays in digsbl_ension i
double gsbl_sp4, gsbl_sp4a, gsbl_sp4b, gsbl_sp4c;   // gsbl_ax size of **sp** arrays in digsbl_ension i
double gsbl_sp5, gsbl_sp5a, gsbl_sp5b, gsbl_sp5c;   // gsbl_ax size of **sp** arrays in digsbl_ension i
double gsbl_6A, gsbl_6Z, gsbl_7K;   // gsbl_ax size of gsbl_** arrays in digsbl_ension i
double gsbl_8A, gsbl_8B, gsbl_8K;   // gsbl_ax size of gsbl_** arrays in digsbl_ension i
double gsbl_9A, gsbl_9B, gsbl_9K;   // gsbl_ax size of gsbl_** arrays in digsbl_ension i
double gsbl_10A, gsbl_10B, gsbl_10K, gsbl_10W;  // gsbl_ax size of gsbl_** arrays in digsbl_ension i
double gsbl_11A, gsbl_11B, gsbl_11C, gsbl_11E, gsbl_11F, gsbl_11W;  // gsbl_ax size of gsbl_** arrays in digsbl_ension i
double gsbl_12A, gsbl_12B, gsbl_12D, gsbl_12E, gsbl_12K;    // gsbl_ax size of gsbl_** arrays in digsbl_ension i
double gsbl_13A, gsbl_13B, gsbl_13K;    // gsbl_ax size of gsbl_** arrays in digsbl_ension i
double gsbl_FCC, gsbl_HCP, gsbl_BCC_9, gsbl_BCC_15; // gsbl_ax size of **sp** arrays in digsbl_ension i

double mean_bl_mom_sp3, mean_bl_mom_sp3a, mean_bl_mom_sp3b, mean_bl_mom_sp3c;   // mean_bl_mom_ax size of **sp** arrays in dimean_bl_mom_ension i
double mean_bl_mom_sp4, mean_bl_mom_sp4a, mean_bl_mom_sp4b, mean_bl_mom_sp4c;   // mean_bl_mom_ax size of **sp** arrays in dimean_bl_mom_ension i
double mean_bl_mom_sp5, mean_bl_mom_sp5a, mean_bl_mom_sp5b, mean_bl_mom_sp5c;   // mean_bl_mom_ax size of **sp** arrays in dimean_bl_mom_ension i
double mean_bl_mom_6A, mean_bl_mom_6Z, mean_bl_mom_7K;  // mean_bl_mom_ax size of mean_bl_mom_** arrays in dimean_bl_mom_ension i
double mean_bl_mom_8A, mean_bl_mom_8B, mean_bl_mom_8K;  // mean_bl_mom_ax size of mean_bl_mom_** arrays in dimean_bl_mom_ension i
double mean_bl_mom_9A, mean_bl_mom_9B, mean_bl_mom_9K;  // mean_bl_mom_ax size of mean_bl_mom_** arrays in dimean_bl_mom_ension i
double mean_bl_mom_10A, mean_bl_mom_10B, mean_bl_mom_10K, mean_bl_mom_10W;  // mean_bl_mom_ax size of mean_bl_mom_** arrays in dimean_bl_mom_ension i
double mean_bl_mom_11A, mean_bl_mom_11B, mean_bl_mom_11C, mean_bl_mom_11E, mean_bl_mom_11F, mean_bl_mom_11W;    // mean_bl_mom_ax size of mean_bl_mom_** arrays in dimean_bl_mom_ension i
double mean_bl_mom_12A, mean_bl_mom_12B, mean_bl_mom_12D, mean_bl_mom_12E, mean_bl_mom_12K; // mean_bl_mom_ax size of mean_bl_mom_** arrays in dimean_bl_mom_ension i
double mean_bl_mom_13A, mean_bl_mom_13B, mean_bl_mom_13K;   // mean_bl_mom_ax size of mean_bl_mom_** arrays in dimean_bl_mom_ension i
double mean_bl_mom_FCC, mean_bl_mom_HCP, mean_bl_mom_BCC_9, mean_bl_mom_BCC_15; // mean_bl_mom_ax size of **sp** arrays in dimean_bl_mom_ension i

double *bl_mom_sp3, *bl_mom_sp3a, *bl_mom_sp3b, *bl_mom_sp3c;   // *bl_mom_ax size of **sp** arrays in di*bl_mom_ension i
double *bl_mom_sp4,*bl_mom_sp4a, *bl_mom_sp4b, *bl_mom_sp4c;    // *bl_mom_ax size of **sp** arrays in di*bl_mom_ension i
double *bl_mom_sp5, *bl_mom_sp5a, *bl_mom_sp5b, *bl_mom_sp5c;   // *bl_mom_ax size of **sp** arrays in di*bl_mom_ension i
double *bl_mom_6A, *bl_mom_6Z, *bl_mom_7K;  // *bl_mom_ax size of *bl_mom_** arrays in di*bl_mom_ension i
double *bl_mom_8A, *bl_mom_8B, *bl_mom_8K;  // *bl_mom_ax size of *bl_mom_** arrays in di*bl_mom_ension i
double *bl_mom_9A, *bl_mom_9B, *bl_mom_9K;  // *bl_mom_ax size of *bl_mom_** arrays in di*bl_mom_ension i
double *bl_mom_10A, *bl_mom_10B, *bl_mom_10K, *bl_mom_10W;  // *bl_mom_ax size of *bl_mom_** arrays in di*bl_mom_ension i
double *bl_mom_11A, *bl_mom_11B, *bl_mom_11C, *bl_mom_11E, *bl_mom_11F, *bl_mom_11W;    // *bl_mom_ax size of *bl_mom_** arrays in di*bl_mom_ension i
double *bl_mom_12A, *bl_mom_12B, *bl_mom_12D, *bl_mom_12E, *bl_mom_12K; // *bl_mom_ax size of *bl_mom_** arrays in di*bl_mom_ension i
double *bl_mom_13A, *bl_mom_13B, *bl_mom_13K;   // *bl_mom_ax size of *bl_mom_** arrays in di*bl_mom_ension i
double *bl_mom_FCC, *bl_mom_HCP, *bl_mom_BCC_9, *bl_mom_BCC_15; // *bl_mom_ax size of **sp** arrays in di*bl_mom_ension i

int **sp3a, **sp3b, **sp3c; // sp3a/b/c arrays (index i denotes number of cluster, j lists particles in cluster)
int **sp4a, **sp4b, **sp4c; // sp5a/b/c arrays (index i denotes number of cluster, j lists particles in cluster)
int **hc6A;
int **sp5a, **sp5b, **sp5c; // sp6a/b/c arrays (index i denotes number of cluster, j lists particles in cluster)
int **hc6Z, **hc7K; // cluster storage arrays (index i denotes number/identifier of cluster, j lists particles in cluster)
int **hc8A, **hc8B, **hc8K;
int **hc9A, **hc9B, **hc9K;
int **hc10A, **hc10B, **hc10K, **hc10W;
int **hc11A, **hc11B, **hc11C, **hc11E, **hc11F, **hc11W;  // cluster storage arrays (index i denotes number/identifier of cluster, j lists particles in cluster)
int **hc12A, **hc12B, **hc12D, **hc12E, **hc12K;
int **hc13A, **hc13B, **hc13K;
int **hcFCC, **hcHCP, **hcBCC_9, **hcBCC_15;

int **mem_sp3b, *nmem_sp3b, mmem_sp3b;
int **mem_sp3c, *nmem_sp3c, mmem_sp3c;
int **mem_sp4b, *nmem_sp4b, mmem_sp4b;
int **mem_sp4c, *nmem_sp4c, mmem_sp4c;
int **mem_sp5b, *nmem_sp5b, mmem_sp5b;
int **mem_sp5c, *nmem_sp5c, mmem_sp5c;

char *ssp3, *ssp3a, *ssp3b, *s5A;
char *ssp4, *ssp4a, *ssp4b, *s6A;
char *ssp5, *ssp5a, *ssp5b, *s7A;
char *s6Z, *s7K;
char *s8A, *s8B, *s8K;
char *s9A, *s9B, *s9K;
char *s10A, *s10B, *s10K, *s10W;
char *s11A, *s11B, *s11C, *s11E, *s11F, *s11W;
char *s12A, *s12B, *s12D, *s12E, *s12K;
char *s13A, *s13B, *s13K;
char*sFCC, *sHCP, *sBCC_9, *sBCC_15;

char *s9B_cen, *s9K_cen;
char *s10B_cen, *s10K_cen, *s10W_cen;
char *s11A_cen, *s11B_cen, *s11C_cen, *s11W_cen;
char *s12A_cen, *s12B_cen, *s12K_cen;
char *s13A_cen, *s13B_cen, *s13K_cen;
char *sFCC_cen, *sHCP_cen, *sBCC_9_cen, *sBCC_15_cen;

char *s9B_shell, *s9K_shell;
char *s10B_shell, *s10K_shell, *s10W_shell;
char *s11A_shell, *s11B_shell, *s11C_shell, *s11W_shell;
char *s12A_shell, *s12B_shell, *s12K_shell;
char *s13A_shell, *s13B_shell, *s13K_shell;
char *sFCC_shell, *sHCP_shell, *sBCC_9_shell, *sBCC_15_shell;

double *pop_per_frame_sp3, *pop_per_frame_sp3a, *pop_per_frame_sp3b, *pop_per_frame_sp3c;
double *pop_per_frame_sp4, *pop_per_frame_sp4a, *pop_per_frame_sp4b, *pop_per_frame_6A;
double *pop_per_frame_sp5, *pop_per_frame_sp5a, *pop_per_frame_sp5b, *pop_per_frame_sp5c;
double *pop_per_frame_6Z, *pop_per_frame_7K;
double *pop_per_frame_8A, *pop_per_frame_8B, *pop_per_frame_8K;
double *pop_per_frame_9A, *pop_per_frame_9B, *pop_per_frame_9K;
double *pop_per_frame_10A, *pop_per_frame_10B, *pop_per_frame_10K, *pop_per_frame_10W;
double *pop_per_frame_11A, *pop_per_frame_11B, *pop_per_frame_11C, *pop_per_frame_11E, *pop_per_frame_11F, *pop_per_frame_11W;
double *pop_per_frame_12A, *pop_per_frame_12B, *pop_per_frame_12D, *pop_per_frame_12E, *pop_per_frame_12K;
double *pop_per_frame_13A, *pop_per_frame_13B, *pop_per_frame_13K;
double *pop_per_frame_FCC, *pop_per_frame_HCP, *pop_per_frame_BCC_9, *pop_per_frame_BCC_15;

double mean_pop_per_frame_sp3, mean_pop_per_frame_sp3a, mean_pop_per_frame_sp3b, mean_pop_per_frame_sp3c;
double mean_pop_per_frame_sp4, mean_pop_per_frame_sp4a, mean_pop_per_frame_sp4b, mean_pop_per_frame_6A;
double mean_pop_per_frame_sp5, mean_pop_per_frame_sp5a, mean_pop_per_frame_sp5b, mean_pop_per_frame_sp5c;
double mean_pop_per_frame_6Z, mean_pop_per_frame_7K;
double mean_pop_per_frame_8A, mean_pop_per_frame_8B, mean_pop_per_frame_8K;
double mean_pop_per_frame_9A, mean_pop_per_frame_9B, mean_pop_per_frame_9K;
double mean_pop_per_frame_10A, mean_pop_per_frame_10B, mean_pop_per_frame_10K, mean_pop_per_frame_10W;
double mean_pop_per_frame_11A, mean_pop_per_frame_11B, mean_pop_per_frame_11C, mean_pop_per_frame_11E, mean_pop_per_frame_11F, mean_pop_per_frame_11W;
double mean_pop_per_frame_12A, mean_pop_per_frame_12B, mean_pop_per_frame_12D, mean_pop_per_frame_12E, mean_pop_per_frame_12K;
double mean_pop_per_frame_13A, mean_pop_per_frame_13B, mean_pop_per_frame_13K;
double mean_pop_per_frame_FCC, mean_pop_per_frame_HCP, mean_pop_per_frame_BCC_9, mean_pop_per_frame_BCC_15;

double av_pot_sp3, av_pot_sp3a, av_pot_sp3b, av_pot_sp3c;
double av_pot_sp4, av_pot_sp4a, av_pot_sp4b, av_pot_sp4c;
double av_pot_sp5, av_pot_sp5a, av_pot_sp5b, av_pot_sp5c;
double av_pot_6A, av_pot_6Z, av_pot_7K;
double av_pot_8A, av_pot_8B, av_pot_8K;
double av_pot_9A, av_pot_9B, av_pot_9K;
double av_pot_10A, av_pot_10B, av_pot_10K, av_pot_10W;
double av_pot_11A, av_pot_11B, av_pot_11C, av_pot_11E, av_pot_11F, av_pot_11W;
double av_pot_12A, av_pot_12B, av_pot_12D, av_pot_12E, av_pot_12K;
double av_pot_13A, av_pot_13B, av_pot_13K;
double av_pot_FCC, av_pot_HCP, av_pot_BCC_9, av_pot_BCC_15;

double av_pot_cen_9B, av_pot_cen_9K;
double av_pot_cen_10B, av_pot_cen_10K, av_pot_cen_10W;
double av_pot_cen_11A, av_pot_cen_11B, av_pot_cen_11C, av_pot_cen_11W;
double av_pot_cen_12A, av_pot_cen_12B, av_pot_cen_12K;
double av_pot_cen_13A, av_pot_cen_13B, av_pot_cen_13K;
double av_pot_cen_FCC, av_pot_cen_HCP, av_pot_cen_BCC_9, av_pot_cen_BCC_15;

int cnt_av_pot_cen_9B, cnt_av_pot_cen_9K;
int cnt_av_pot_cen_10B, cnt_av_pot_cen_10K, cnt_av_pot_cen_10W;
int cnt_av_pot_cen_11A, cnt_av_pot_cen_11B, cnt_av_pot_cen_11C, cnt_av_pot_cen_11W;
int cnt_av_pot_cen_12A, cnt_av_pot_cen_12B, cnt_av_pot_cen_12K;
int cnt_av_pot_cen_13A, cnt_av_pot_cen_13B, cnt_av_pot_cen_13K;
int cnt_av_pot_cen_FCC, cnt_av_pot_cen_HCP, cnt_av_pot_cen_BCC_9, cnt_av_pot_cen_BCC_15;

double av_pot_shell_9B, av_pot_shell_9K;
double av_pot_shell_10B, av_pot_shell_10K, av_pot_shell_10W;
double av_pot_shell_11A, av_pot_shell_11B, av_pot_shell_11C, av_pot_shell_11W;
double av_pot_shell_12A, av_pot_shell_12B, av_pot_shell_12K;
double av_pot_shell_13A, av_pot_shell_13B, av_pot_shell_13K;
double av_pot_shell_FCC, av_pot_shell_HCP, av_pot_shell_BCC_9, av_pot_shell_BCC_15;

int cnt_av_pot_shell_9B, cnt_av_pot_shell_9K;
int cnt_av_pot_shell_10B, cnt_av_pot_shell_10K, cnt_av_pot_shell_10W;
int cnt_av_pot_shell_11A, cnt_av_pot_shell_11B, cnt_av_pot_shell_11C, cnt_av_pot_shell_11W;
int cnt_av_pot_shell_12A, cnt_av_pot_shell_12B, cnt_av_pot_shell_12K;
int cnt_av_pot_shell_13A, cnt_av_pot_shell_13B, cnt_av_pot_shell_13K;
int cnt_av_pot_shell_FCC, cnt_av_pot_shell_HCP, cnt_av_pot_shell_BCC_9, cnt_av_pot_shell_BCC_15;

int *n_distro_bonded_to_cen_9B, *n_distro_bonded_to_cen_9K;
int *n_distro_bonded_to_cen_10B, *n_distro_bonded_to_cen_10K, *n_distro_bonded_to_cen_10W;
int *n_distro_bonded_to_cen_11A, *n_distro_bonded_to_cen_11B, *n_distro_bonded_to_cen_11C, *n_distro_bonded_to_cen_11W;
int *n_distro_bonded_to_cen_12A, *n_distro_bonded_to_cen_12B, *n_distro_bonded_to_cen_12K;
int *n_distro_bonded_to_cen_13A, *n_distro_bonded_to_cen_13B, *n_distro_bonded_to_cen_13K;
int *n_distro_bonded_to_cen_FCC, *n_distro_bonded_to_cen_HCP, *n_distro_bonded_to_cen_BCC_9, *n_distro_bonded_to_cen_BCC_15;

int n_bonded_to_cen_9B, n_bonded_to_cen_9K;
int n_bonded_to_cen_10B, n_bonded_to_cen_10K, n_bonded_to_cen_10W;
int n_bonded_to_cen_11A, n_bonded_to_cen_11B, n_bonded_to_cen_11C, n_bonded_to_cen_11W;
int n_bonded_to_cen_12A, n_bonded_to_cen_12B, n_bonded_to_cen_12K;
int n_bonded_to_cen_13A, n_bonded_to_cen_13B, n_bonded_to_cen_13K;
int n_bonded_to_cen_FCC, n_bonded_to_cen_HCP, n_bonded_to_cen_BCC_9, n_bonded_to_cen_BCC_15;

int n_distro_sp3[4], n_distro_sp3a[4], n_distro_sp3b[5], n_distro_sp3c[6];
int n_distro_sp4[5], n_distro_sp4a[5], n_distro_sp4b[6], n_distro_sp4c[7], n_distro_6A[7];
int n_distro_sp5[6], n_distro_sp5a[6], n_distro_sp5b[7], n_distro_sp5c[8];
int n_distro_6Z[7], n_distro_7K[8];
int n_distro_8A[9], n_distro_8B[9], n_distro_8K[9];
int n_distro_9A[10], n_distro_9B[10], n_distro_9K[10];
int n_distro_10A[11], n_distro_10B[11], n_distro_10K[11], n_distro_10W[11];
int n_distro_11A[12], n_distro_11B[12], n_distro_11C[12], n_distro_11E[12], n_distro_11F[12], n_distro_11W[12];
int n_distro_12A[13], n_distro_12B[13], n_distro_12D[13], n_distro_12E[13], n_distro_12K[13];
int n_distro_13A[14], n_distro_13B[14], n_distro_13K[14];
int n_distro_FCC[14], n_distro_HCP[14], n_distro_BCC_9[10], n_distro_BCC_15[16];

int n_distro_cen_9B[2], n_distro_cen_9K[2];
int n_distro_cen_10B[2], n_distro_cen_10K[2], n_distro_cen_10W[2];
int n_distro_cen_11A[2], n_distro_cen_11B[2], n_distro_cen_11C[2], n_distro_cen_11W[2];
int n_distro_cen_12A[2], n_distro_cen_12B[2], n_distro_cen_12K[2];
int n_distro_cen_13A[2], n_distro_cen_13B[2], n_distro_cen_13K[2];
int n_distro_cen_FCC[2], n_distro_cen_HCP[2], n_distro_cen_BCC_9[2], n_distro_cen_BCC_15[2];

int n_distro_shell_9B[9], n_distro_shell_9K[9];
int n_distro_shell_10B[10], n_distro_shell_10K[10], n_distro_shell_10W[10];
int n_distro_shell_11A[11], n_distro_shell_11B[11], n_distro_shell_11C[11], n_distro_shell_11W[11];
int n_distro_shell_12A[12], n_distro_shell_12B[12], n_distro_shell_12K[12];
int n_distro_shell_13A[13], n_distro_shell_13B[13], n_distro_shell_13K[13];
int n_distro_shell_FCC[13], n_distro_shell_HCP[13], n_distro_shell_BCC_9[9], n_distro_shell_BCC_15[15];

int nAsp3, nAsp3a, nAsp3b, nAsp3c;
int nAsp4, nAsp4a, nAsp4b, nAsp4c, nA6A;
int nAsp5, nAsp5a, nAsp5b, nAsp5c;
int nA6Z, nA7K;
int nA8A, nA8B, nA8K;
int nA9A, nA9B, nA9K;
int nA10A, nA10B, nA10K, nA10W;
int nA11A, nA11B, nA11C, nA11E, nA11F, nA11W;
int nA12A, nA12B, nA12D, nA12E, nA12K;
int nA13A, nA13B, nA13K;
int nAFCC, nAHCP, nABCC_9, nABCC_15;

int nA_cen_9B, nA_cen_9K;
int nA_cen_10B, nA_cen_10K, nA_cen_10W;
int nA_cen_11A, nA_cen_11B, nA_cen_11C, nA_cen_11W;
int nA_cen_12A, nA_cen_12B, nA_cen_12K;
int nA_cen_13A, nA_cen_13B, nA_cen_13K;
int nA_cen_FCC, nA_cen_HCP, nA_cen_BCC_9, nA_cen_BCC_15;

int nA_shell_9B, nA_shell_9K;
int nA_shell_10B, nA_shell_10K, nA_shell_10W;
int nA_shell_11A, nA_shell_11B, nA_shell_11C, nA_shell_11W;
int nA_shell_12A, nA_shell_12B, nA_shell_12K;
int nA_shell_13A, nA_shell_13B, nA_shell_13K;
int nA_shell_FCC, nA_shell_HCP, nA_shell_BCC_9, nA_shell_BCC_15;

int nBsp3, nBsp3a, nBsp3b, nBsp3c;
int nBsp4, nBsp4a, nBsp4b, nBsp4c, nB6A;
int nBsp5, nBsp5a, nBsp5b, nBsp5c;
int nB6Z, nB7K;
int nB8A, nB8B, nB8K;
int nB9A, nB9B, nB9K;
int nB10A, nB10B, nB10K, nB10W;
int nB11A, nB11B, nB11C, nB11E, nB11F, nB11W;
int nB12A, nB12B, nB12D, nB12E, nB12K;
int nB13A, nB13B, nB13K;
int nBFCC, nBHCP, nBBCC_9, nBBCC_15;

int nB_cen_9B, nB_cen_9K;
int nB_cen_10B, nB_cen_10K, nB_cen_10W;
int nB_cen_11A, nB_cen_11B, nB_cen_11C, nB_cen_11W;
int nB_cen_12A, nB_cen_12B, nB_cen_12K;
int nB_cen_13A, nB_cen_13B, nB_cen_13K;
int nB_cen_FCC, nB_cen_HCP, nB_cen_BCC_9, nB_cen_BCC_15;

int nB_shell_9B, nB_shell_9K;
int nB_shell_10B, nB_shell_10K, nB_shell_10W;
int nB_shell_11A, nB_shell_11B, nB_shell_11C, nB_shell_11W;
int nB_shell_12A, nB_shell_12B, nB_shell_12K;
int nB_shell_13A, nB_shell_13B, nB_shell_13K;
int nB_shell_FCC, nB_shell_HCP, nB_shell_BCC_9, nB_shell_BCC_15;

int ngsp3, ngsp3a, ngsp3b, ng5A;
int ngsp4, ngsp4a, ngsp4b, ng6A;
int ngsp5, ngsp5a, ngsp5b, ng7A;
int ng6Z, ng7K;
int ng8A, ng8B, ng8K;
int ng9A, ng9B, ng9K;
int ng10A, ng10B, ng10K, ng10W;
int ng11A, ng11B, ng11C, ng11E, ng11F, ng11W;
int ng12A, ng12B, ng12D, ng12E, ng12K;
int ng13A, ng13B, ng13K;
int ngFCC, ngHCP, ngBCC_9, ngBCC_15;

int nn5A;
int nn6A, nn6Z, nn7K;
int nn7A;
int nn8A, nn8B, nn8K;
int nn9A, nn9B, nn9K;
int nn10A, nn10B, nn10K, nn10W;
int nn11A, nn11B, nn11C, nn11E, nn11F, nn11W;   // number of particles in clusers net
int nn12A, nn12B, nn12D, nn12E, nn12K;
int nn13A, nn13B, nn13K;
int nnFCC, nnHCP, nnBCC_9, nnBCC_15;

int ncsp3, ncsp3a, ncsp3b, nc5A;
int ncsp4, ncsp4a, ncsp4b, ncsp4c, nc6A;
int ncsp5, ncsp5a, ncsp5b, nc7A;    // total number of clusers
int nc6Z, nc7K;
int nc8A, nc8B, nc8K;
int nc9A, nc9B, nc9K;
int nc10A, nc10B, nc10K, nc10W;
int nc11A, nc11B, nc11C, nc11E, nc11F, nc11W;   // total number of clusers
int nc12A, nc12B, nc12D, nc12E, nc12K;
int nc13A, nc13B, nc13K;
int ncFCC, ncHCP, ncBCC_9, ncBCC_15;    // total number of clusers

int ncsp3_excess_spindles, ncsp4_excess_spindles, ncsp5_excess_spindles;    // total number of _excess_spindlesed basic clusters
int ncsp3c_spindlebonds, ncsp4c_spindlebonds, nc6A_spindlebonds, ncsp5c_spindlebonds;   // total number of spindle bonds

int *a5, *a6, *a7, *a8, *a9, *a10, *a11, *a12, *a13, *a15;

FILE *wsp3, *wsp3a, *wsp3b, *w5A;
FILE *wsp4, *wsp4a, *wsp4b, *wsp4c, *w6A, *w6Z, *w7K;
FILE *wsp5, *wsp5a, *wsp5b, *w7A;
FILE *w8A, *w8B, *w8K;
FILE *w9A, *w9B, *w9K;
FILE *w10A, *w10B, *w10K, *w10W;
FILE *w11A, *w11B, *w11C, *w11E, *w11F, *w11W;
FILE *w12A, *w12B, *w12D, *w12E, *w12K;
FILE *w13A, *w13B, *w13K;
FILE *wFCC, *wHCP, *wBCC_9, *wBCC_15;

FILE *bondsout;
FILE *fPopPerFrame;

FILE *file_raw_sp3, *file_raw_sp3a, *file_raw_sp3b, *file_raw_5A;
FILE *file_raw_sp4, *file_raw_sp4a, *file_raw_sp4b, *file_raw_6A;
FILE *file_raw_sp5, *file_raw_sp5a, *file_raw_sp5b, *file_raw_7A;
FILE *file_raw_6Z, *file_raw_7K;
FILE *file_raw_8A, *file_raw_8B, *file_raw_8K;  
FILE *file_raw_9A, *file_raw_9B, *file_raw_9K;
FILE *file_raw_10A, *file_raw_10B, *file_raw_10K, *file_raw_10W;
FILE *file_raw_11A, *file_raw_11B, *file_raw_11C, *file_raw_11E, *file_raw_11F, *file_raw_11W;
FILE *file_raw_12A, *file_raw_12B, *file_raw_12D, *file_raw_12E, *file_raw_12K;
FILE *file_raw_13A, *file_raw_13B, *file_raw_13K;
FILE *file_raw_FCC, *file_raw_HCP, *file_raw_BCC_9, *file_raw_BCC_15;

FILE *file_11A_cen_xmol;
FILE *file_13A_cen_xmol;

// dynamic variables
int dyn_MaxLives; // maximum number of times I witness one cluster
int dyn_TotClus, dyn_TotEvents, dyn_Tot_a_Events;
int dyn_nsp3, dyn_nsp3a, dyn_nsp3b, dyn_nsp3c;  // nudyn_mber of sp3a/b/c respectively
int dyn_nsp4, dyn_nsp4a, dyn_nsp4b, dyn_n6A;    // nudyn_mber of sp4a/b 6A respectively
int dyn_nsp5, dyn_nsp5a, dyn_nsp5b, dyn_nsp5c;  // nudyn_mber of sp5a/b/c respectively
int dyn_n6Z, dyn_n7K;   // nudyn_mber of clusters of particlar type
int dyn_n8A, dyn_n8B, dyn_n8K;  
int dyn_n9A, dyn_n9B, dyn_n9K;
int dyn_n10A, dyn_n10B, dyn_n10K, dyn_n10W;
int dyn_n11A, dyn_n11B, dyn_n11C, dyn_n11E, dyn_n11F, dyn_n11W; // nudyn_mber of clusters of particlar type
int dyn_n12A, dyn_n12B, dyn_n12D, dyn_n12E, dyn_n12K;
int dyn_n13A, dyn_n13B, dyn_n13K;
int dyn_nFCC, dyn_nHCP, dyn_nBCC_9, dyn_nBCC_15;

int dyn_esp3, dyn_esp3a, dyn_esp3b, dyn_esp3c;  // nudyn_mber of sp3a/b/c respectively
int dyn_esp4, dyn_esp4a, dyn_esp4b, dyn_e6A;    // nudyn_mber of sp4a/b 6A respectively
int dyn_esp5, dyn_esp5a, dyn_esp5b, dyn_esp5c;  // nudyn_mber of sp5a/b/c respectively
int dyn_e6Z, dyn_e7K;   // nudyn_mber of clusters of particlar type
int dyn_e8A, dyn_e8B, dyn_e8K;  
int dyn_e9A, dyn_e9B, dyn_e9K;
int dyn_e10A, dyn_e10B, dyn_e10K, dyn_e10W;
int dyn_e11A, dyn_e11B, dyn_e11C, dyn_e11E, dyn_e11F, dyn_e11W; // nudyn_mber of clusters of particlar type
int dyn_e12A, dyn_e12B, dyn_e12D, dyn_e12E, dyn_e12K;
int dyn_e13A, dyn_e13B, dyn_e13K;
int dyn_eFCC, dyn_eHCP, dyn_eBCC_9, dyn_eBCC_15;

int dyn_asp3, dyn_asp3a, dyn_asp3b, dyn_asp3c;  // nudyn_mber of sp3a/b/c respectively
int dyn_asp4, dyn_asp4a, dyn_asp4b, dyn_a6A;    // nudyn_mber of sp4a/b 6A respectively
int dyn_asp5, dyn_asp5a, dyn_asp5b, dyn_asp5c;  // nudyn_mber of sp5a/b/c respectively
int dyn_a6Z, dyn_a7K;   // nudyn_mber of clusters of particlar type
int dyn_a8A, dyn_a8B, dyn_a8K;  
int dyn_a9A, dyn_a9B, dyn_a9K;
int dyn_a10A, dyn_a10B, dyn_a10K, dyn_a10W;
int dyn_a11A, dyn_a11B, dyn_a11C, dyn_a11E, dyn_a11F, dyn_a11W; // nudyn_mber of clusters of particlar type
int dyn_a12A, dyn_a12B, dyn_a12D, dyn_a12E, dyn_a12K;
int dyn_a13A, dyn_a13B, dyn_a13K;
int dyn_aFCC, dyn_aHCP, dyn_aBCC_9, dyn_aBCC_15;

int dyn_msp3, dyn_msp3a, dyn_msp3b, dyn_msp3c;  // dyn_max size of **sp** arrays in didyn_mension i
int dyn_msp4, dyn_msp4a, dyn_msp4b, dyn_m6A;    // dyn_max size of **sp** arrays in didyn_mension i
int dyn_msp5, dyn_msp5a, dyn_msp5b, dyn_msp5c;  // dyn_max size of **sp** arrays in didyn_mension i
int dyn_m6Z, dyn_m7K;   // dyn_max size of dyn_m** arrays in didyn_mension i
int dyn_m8A, dyn_m8B, dyn_m8K;  // dyn_max size of dyn_m** arrays in didyn_mension i
int dyn_m9A, dyn_m9B, dyn_m9K;  // dyn_max size of dyn_m** arrays in didyn_mension i
int dyn_m10A, dyn_m10B, dyn_m10K, dyn_m10W; // dyn_max size of dyn_m** arrays in didyn_mension i
int dyn_m11A, dyn_m11B, dyn_m11C, dyn_m11E, dyn_m11F, dyn_m11W; // dyn_max size of dyn_m** arrays in didyn_mension i
int dyn_m12A, dyn_m12B, dyn_m12D, dyn_m12E, dyn_m12K;   // dyn_max size of dyn_m** arrays in didyn_mension i
int dyn_m13A, dyn_m13B, dyn_m13K;   // max size of dyn_m** arrays in dimension i
int dyn_mFCC, dyn_mHCP, dyn_mBCC_9, dyn_mBCC_15;    // dyn_max size of **sp** arrays in didyn_mension i

double sum_hist_sp3, sum_hist_sp3a, sum_hist_sp3b, sum_hist_sp3c;   // sum_hist_ax size of **sp** arrays in disum_hist_ension i
double sum_hist_sp4, sum_hist_sp4a, sum_hist_sp4b, sum_hist_6A; // sum_hist_ax size of **sp** arrays in disum_hist_ension i
double sum_hist_sp5, sum_hist_sp5a, sum_hist_sp5b, sum_hist_sp5c;   // sum_hist_ax size of **sp** arrays in disum_hist_ension i
double sum_hist_6Z, sum_hist_7K;    // sum_hist_ax size of sum_hist_** arrays in disum_hist_ension i
double sum_hist_8A, sum_hist_8B, sum_hist_8K;   // sum_hist_ax size of sum_hist_** arrays in disum_hist_ension i
double sum_hist_9A, sum_hist_9B, sum_hist_9K;   // sum_hist_ax size of sum_hist_** arrays in disum_hist_ension i
double sum_hist_10A, sum_hist_10B, sum_hist_10K, sum_hist_10W;  // sum_hist_ax size of sum_hist_** arrays in disum_hist_ension i
double sum_hist_11A, sum_hist_11B, sum_hist_11C, sum_hist_11E, sum_hist_11F, sum_hist_11W;  // sum_hist_ax size of sum_hist_** arrays in disum_hist_ension i
double sum_hist_12A, sum_hist_12B, sum_hist_12D, sum_hist_12E, sum_hist_12K;    // sum_hist_ax size of sum_hist_** arrays in disum_hist_ension i
double sum_hist_13A, sum_hist_13B, sum_hist_13K;    // max size of sum_hist_** arrays in dimension i
double sum_hist_FCC, sum_hist_HCP, sum_hist_BCC_9, sum_hist_BCC_15; // sum_hist_ax size of **sp** arrays in disum_hist_ension i

double sum_max_hist_sp3, sum_max_hist_sp3a, sum_max_hist_sp3b, sum_max_hist_sp3c;   // sum_max_hist_ax size of **sp** arrays in disum_max_hist_ension i
double sum_max_hist_sp4, sum_max_hist_sp4a, sum_max_hist_sp4b, sum_max_hist_6A; // sum_max_hist_ax size of **sp** arrays in disum_max_hist_ension i
double sum_max_hist_sp5, sum_max_hist_sp5a, sum_max_hist_sp5b, sum_max_hist_sp5c;   // sum_max_hist_ax size of **sp** arrays in disum_max_hist_ension i
double sum_max_hist_6Z, sum_max_hist_7K;    // sum_max_hist_ax size of sum_max_hist_** arrays in disum_max_hist_ension i
double sum_max_hist_8A, sum_max_hist_8B, sum_max_hist_8K;   // sum_max_hist_ax size of sum_max_hist_** arrays in disum_max_hist_ension i
double sum_max_hist_9A, sum_max_hist_9B, sum_max_hist_9K;   // sum_max_hist_ax size of sum_max_hist_** arrays in disum_max_hist_ension i
double sum_max_hist_10A, sum_max_hist_10B, sum_max_hist_10K, sum_max_hist_10W;  // sum_max_hist_ax size of sum_max_hist_** arrays in disum_max_hist_ension i
double sum_max_hist_11A, sum_max_hist_11B, sum_max_hist_11C, sum_max_hist_11E, sum_max_hist_11F, sum_max_hist_11W;  // sum_max_hist_ax size of sum_max_hist_** arrays in disum_max_hist_ension i
double sum_max_hist_12A, sum_max_hist_12B, sum_max_hist_12D, sum_max_hist_12E, sum_max_hist_12K;    // sum_max_hist_ax size of sum_max_hist_** arrays in disum_max_hist_ension i
double sum_max_hist_13A, sum_max_hist_13B, sum_max_hist_13K;    // max size of sum_max_hist_** arrays in dimension i
double sum_max_hist_FCC, sum_max_hist_HCP, sum_max_hist_BCC_9, sum_max_hist_BCC_15; // sum_max_hist_ax size of **sp** arrays in disum_max_hist_ension i

int **dyn_sp3, **dyn_sp3a, **dyn_sp3b, **dyn_sp3c;  // sp3a/b/c arrays (index i denotes nudyn_mber of cluster, j lists particles in cluster)
int **dyn_sp4, **dyn_sp4a, **dyn_sp4b, **dyn_hc6A;  // sp5a/b/c arrays (index i denotes nudyn_mber of cluster, j lists particles in cluster)
int **dyn_sp5, **dyn_sp5a, **dyn_sp5b, **dyn_sp5c;  // sp6a/b/c arrays (index i denotes number of cluster, j lists particles in cluster)
int **dyn_hc6Z, **dyn_hc7K; // cluster storage arrays (index i denotes nudyn_mber/identifier of cluster, j lists particles in cluster)
int **dyn_hc8A, **dyn_hc8B, **dyn_hc8K;
int **dyn_hc9A, **dyn_hc9B, **dyn_hc9K;
int **dyn_hc10A, **dyn_hc10B, **dyn_hc10K, **dyn_hc10W;
int **dyn_hc11A, **dyn_hc11B, **dyn_hc11C, **dyn_hc11E, **dyn_hc11F, **dyn_hc11W;  // cluster storage arrays (index i denotes nudyn_mber/identifier of cluster, j lists particles in cluster)
int **dyn_hc12A, **dyn_hc12B, **dyn_hc12D, **dyn_hc12E, **dyn_hc12K;
int **dyn_hc13A, **dyn_hc13B, **dyn_hc13K;
int **dyn_hcFCC, **dyn_hcHCP, **dyn_hcBCC_9, **dyn_hcBCC_15;

int **dyn_lsp3, **dyn_lsp3a, **dyn_lsp3b, **dyn_lsp3c;  // sp3a/b/c arrays (index i denotes nudyn_lmber of cluster, j lists particles in cluster)
int **dyn_lsp4, **dyn_lsp4a, **dyn_lsp4b, **dyn_l6A;    // sp5a/b/c arrays (index i denotes nudyn_lmber of cluster, j lists particles in cluster)
int **dyn_lsp5, **dyn_lsp5a, **dyn_lsp5b, **dyn_lsp5c;  // sp6a/b/c arrays (index i denotes number of cluster, j lists particles in cluster)
int **dyn_l6Z, **dyn_l7K;   // cluster storage arrays (index i denotes nudyn_mber/identifier of cluster, j lists particles in cluster)
int **dyn_l8A, **dyn_l8B, **dyn_l8K;
int **dyn_l9A, **dyn_l9B, **dyn_l9K;
int **dyn_l10A, **dyn_l10B, **dyn_l10K, **dyn_l10W;
int **dyn_l11A, **dyn_l11B, **dyn_l11C, **dyn_l11E, **dyn_l11F, **dyn_l11W;  // cluster storage arrays (index i denotes nudyn_mber/identifier of cluster, j lists particles in cluster)
int **dyn_l12A, **dyn_l12B, **dyn_l12D, **dyn_l12E, **dyn_l12K;
int **dyn_l13A, **dyn_l13B, **dyn_l13K;
int **dyn_lFCC, **dyn_lHCP, **dyn_lBCC_9, **dyn_lBCC_15;

int **dyn_sub_6Z, **dyn_sub_7K; // cluster storage arrays (index i denotes nudyn_mber/identifier of cluster, j lists particles in cluster)
int **dyn_sub_8A, **dyn_sub_8B, **dyn_sub_8K;
int **dyn_sub_9A, **dyn_sub_9B, **dyn_sub_9K;
int **dyn_sub_10A, **dyn_sub_10B, **dyn_sub_10K, **dyn_sub_10W;
int **dyn_sub_11A, **dyn_sub_11B, **dyn_sub_11C, **dyn_sub_11E, **dyn_sub_11F, **dyn_sub_11W;  // cluster storage arrays (index i denotes nudyn_mber/identifier of cluster, j lists particles in cluster)
int **dyn_sub_12A, **dyn_sub_12B, **dyn_sub_12D, **dyn_sub_12E, **dyn_sub_12K;
int **dyn_sub_13B, **dyn_sub_13K;
int **dyn_sub_FCC, **dyn_sub_HCP, **dyn_sub_BCC_9, **dyn_sub_BCC_15;

int *dyn_up_sp3b, *dyn_up_sp3c;
int *dyn_up_sp4b, *dyn_up_sp4c;
int *dyn_up_sp5b, *dyn_up_sp5c;
int *dyn_up_9B, *dyn_up_9K, *dyn_up_10B, *dyn_up_11A, *dyn_up_11C, *dyn_up_11F;

int *dyn_hist_sp3, *dyn_hist_sp3a, *dyn_hist_sp3b, *dyn_hist_sp3c;
int *dyn_hist_sp4, *dyn_hist_sp4a, *dyn_hist_sp4b, *dyn_hist_6A;
int *dyn_hist_sp5, *dyn_hist_sp5a, *dyn_hist_sp5b, *dyn_hist_sp5c;
int *dyn_hist_6Z, *dyn_hist_7K;
int *dyn_hist_8A, *dyn_hist_8B, *dyn_hist_8K;
int *dyn_hist_9A, *dyn_hist_9B, *dyn_hist_9K;
int *dyn_hist_10A, *dyn_hist_10B, *dyn_hist_10K, *dyn_hist_10W;
int *dyn_hist_11A, *dyn_hist_11B, *dyn_hist_11C, *dyn_hist_11E, *dyn_hist_11F, *dyn_hist_11W;
int *dyn_hist_12A, *dyn_hist_12B, *dyn_hist_12D, *dyn_hist_12E, *dyn_hist_12K;
int *dyn_hist_13A, *dyn_hist_13B, *dyn_hist_13K;
int *dyn_hist_FCC, *dyn_hist_HCP, *dyn_hist_BCC_9, *dyn_hist_BCC_15;

double *dyn_norm_hist_sp3, *dyn_norm_hist_sp3a, *dyn_norm_hist_sp3b, *dyn_norm_hist_sp3c;
double *dyn_norm_hist_sp4, *dyn_norm_hist_sp4a, *dyn_norm_hist_sp4b, *dyn_norm_hist_6A;
double *dyn_norm_hist_sp5, *dyn_norm_hist_sp5a, *dyn_norm_hist_sp5b, *dyn_norm_hist_sp5c;
double *dyn_norm_hist_6Z, *dyn_norm_hist_7K;
double *dyn_norm_hist_8A, *dyn_norm_hist_8B, *dyn_norm_hist_8K;
double *dyn_norm_hist_9A, *dyn_norm_hist_9B, *dyn_norm_hist_9K;
double *dyn_norm_hist_10A, *dyn_norm_hist_10B, *dyn_norm_hist_10K, *dyn_norm_hist_10W;
double *dyn_norm_hist_11A, *dyn_norm_hist_11B, *dyn_norm_hist_11C, *dyn_norm_hist_11E, *dyn_norm_hist_11F, *dyn_norm_hist_11W;
double *dyn_norm_hist_12A, *dyn_norm_hist_12B, *dyn_norm_hist_12D, *dyn_norm_hist_12E, *dyn_norm_hist_12K;
double *dyn_norm_hist_13A, *dyn_norm_hist_13B, *dyn_norm_hist_13K;
double *dyn_norm_hist_FCC, *dyn_norm_hist_HCP, *dyn_norm_hist_BCC_9, *dyn_norm_hist_BCC_15;

double *dyn_norm_corrected_hist_sp3, *dyn_norm_corrected_hist_sp3a, *dyn_norm_corrected_hist_sp3b, *dyn_norm_corrected_hist_sp3c;
double *dyn_norm_corrected_hist_sp4, *dyn_norm_corrected_hist_sp4a, *dyn_norm_corrected_hist_sp4b, *dyn_norm_corrected_hist_6A;
double *dyn_norm_corrected_hist_sp5, *dyn_norm_corrected_hist_sp5a, *dyn_norm_corrected_hist_sp5b, *dyn_norm_corrected_hist_sp5c;
double *dyn_norm_corrected_hist_6Z, *dyn_norm_corrected_hist_7K;
double *dyn_norm_corrected_hist_8A, *dyn_norm_corrected_hist_8B, *dyn_norm_corrected_hist_8K;
double *dyn_norm_corrected_hist_9A, *dyn_norm_corrected_hist_9B, *dyn_norm_corrected_hist_9K;
double *dyn_norm_corrected_hist_10A, *dyn_norm_corrected_hist_10B, *dyn_norm_corrected_hist_10K, *dyn_norm_corrected_hist_10W;
double *dyn_norm_corrected_hist_11A, *dyn_norm_corrected_hist_11B, *dyn_norm_corrected_hist_11C, *dyn_norm_corrected_hist_11E, *dyn_norm_corrected_hist_11F, *dyn_norm_corrected_hist_11W;
double *dyn_norm_corrected_hist_12A, *dyn_norm_corrected_hist_12B, *dyn_norm_corrected_hist_12D, *dyn_norm_corrected_hist_12E, *dyn_norm_corrected_hist_12K;
double *dyn_norm_corrected_hist_13A, *dyn_norm_corrected_hist_13B, *dyn_norm_corrected_hist_13K;
double *dyn_norm_corrected_hist_FCC, *dyn_norm_corrected_hist_HCP, *dyn_norm_corrected_hist_BCC_9, *dyn_norm_corrected_hist_BCC_15;

double *dyn_decay_norm_hist_sp3, *dyn_decay_norm_hist_sp3a, *dyn_decay_norm_hist_sp3b, *dyn_decay_norm_hist_sp3c;
double *dyn_decay_norm_hist_sp4, *dyn_decay_norm_hist_sp4a, *dyn_decay_norm_hist_sp4b, *dyn_decay_norm_hist_6A;
double *dyn_decay_norm_hist_sp5, *dyn_decay_norm_hist_sp5a, *dyn_decay_norm_hist_sp5b, *dyn_decay_norm_hist_sp5c;
double *dyn_decay_norm_hist_6Z, *dyn_decay_norm_hist_7K;
double *dyn_decay_norm_hist_8A, *dyn_decay_norm_hist_8B, *dyn_decay_norm_hist_8K;
double *dyn_decay_norm_hist_9A, *dyn_decay_norm_hist_9B, *dyn_decay_norm_hist_9K;
double *dyn_decay_norm_hist_10A, *dyn_decay_norm_hist_10B, *dyn_decay_norm_hist_10K, *dyn_decay_norm_hist_10W;
double *dyn_decay_norm_hist_11A, *dyn_decay_norm_hist_11B, *dyn_decay_norm_hist_11C, *dyn_decay_norm_hist_11E, *dyn_decay_norm_hist_11F, *dyn_decay_norm_hist_11W;
double *dyn_decay_norm_hist_12A, *dyn_decay_norm_hist_12B, *dyn_decay_norm_hist_12D, *dyn_decay_norm_hist_12E, *dyn_decay_norm_hist_12K;
double *dyn_decay_norm_hist_13A, *dyn_decay_norm_hist_13B, *dyn_decay_norm_hist_13K;
double *dyn_decay_norm_hist_FCC, *dyn_decay_norm_hist_HCP, *dyn_decay_norm_hist_BCC_9, *dyn_decay_norm_hist_BCC_15;

double *dyn_decay_norm_corrected_hist_sp3, *dyn_decay_norm_corrected_hist_sp3a, *dyn_decay_norm_corrected_hist_sp3b, *dyn_decay_norm_corrected_hist_sp3c;
double *dyn_decay_norm_corrected_hist_sp4, *dyn_decay_norm_corrected_hist_sp4a, *dyn_decay_norm_corrected_hist_sp4b, *dyn_decay_norm_corrected_hist_6A;
double *dyn_decay_norm_corrected_hist_sp5, *dyn_decay_norm_corrected_hist_sp5a, *dyn_decay_norm_corrected_hist_sp5b, *dyn_decay_norm_corrected_hist_sp5c;
double *dyn_decay_norm_corrected_hist_6Z, *dyn_decay_norm_corrected_hist_7K;
double *dyn_decay_norm_corrected_hist_8A, *dyn_decay_norm_corrected_hist_8B, *dyn_decay_norm_corrected_hist_8K;
double *dyn_decay_norm_corrected_hist_9A, *dyn_decay_norm_corrected_hist_9B, *dyn_decay_norm_corrected_hist_9K;
double *dyn_decay_norm_corrected_hist_10A, *dyn_decay_norm_corrected_hist_10B, *dyn_decay_norm_corrected_hist_10K, *dyn_decay_norm_corrected_hist_10W;
double *dyn_decay_norm_corrected_hist_11A, *dyn_decay_norm_corrected_hist_11B, *dyn_decay_norm_corrected_hist_11C, *dyn_decay_norm_corrected_hist_11E, *dyn_decay_norm_corrected_hist_11F, *dyn_decay_norm_corrected_hist_11W;
double *dyn_decay_norm_corrected_hist_12A, *dyn_decay_norm_corrected_hist_12B, *dyn_decay_norm_corrected_hist_12D, *dyn_decay_norm_corrected_hist_12E, *dyn_decay_norm_corrected_hist_12K;
double *dyn_decay_norm_corrected_hist_13A, *dyn_decay_norm_corrected_hist_13B, *dyn_decay_norm_corrected_hist_13K;
double *dyn_decay_norm_corrected_hist_FCC, *dyn_decay_norm_corrected_hist_HCP, *dyn_decay_norm_corrected_hist_BCC_9, *dyn_decay_norm_corrected_hist_BCC_15;

int *dyn_total_hist_sp3, *dyn_total_hist_sp3a, *dyn_total_hist_sp3b, *dyn_total_hist_sp3c;
int *dyn_total_hist_sp4, *dyn_total_hist_sp4a, *dyn_total_hist_sp4b, *dyn_total_hist_6A;
int *dyn_total_hist_sp5, *dyn_total_hist_sp5a, *dyn_total_hist_sp5b, *dyn_total_hist_sp5c;
int *dyn_total_hist_6Z, *dyn_total_hist_7K;
int *dyn_total_hist_8A, *dyn_total_hist_8B, *dyn_total_hist_8K;
int *dyn_total_hist_9A, *dyn_total_hist_9B, *dyn_total_hist_9K;
int *dyn_total_hist_10A, *dyn_total_hist_10B, *dyn_total_hist_10K, *dyn_total_hist_10W;
int *dyn_total_hist_11A, *dyn_total_hist_11B, *dyn_total_hist_11C, *dyn_total_hist_11E, *dyn_total_hist_11F, *dyn_total_hist_11W;
int *dyn_total_hist_12A, *dyn_total_hist_12B, *dyn_total_hist_12D, *dyn_total_hist_12E, *dyn_total_hist_12K;
int *dyn_total_hist_13A, *dyn_total_hist_13B, *dyn_total_hist_13K;
int *dyn_total_hist_FCC, *dyn_total_hist_HCP, *dyn_total_hist_BCC_9, *dyn_total_hist_BCC_15;

double *dyn_norm_total_hist_sp3, *dyn_norm_total_hist_sp3a, *dyn_norm_total_hist_sp3b, *dyn_norm_total_hist_sp3c;
double *dyn_norm_total_hist_sp4, *dyn_norm_total_hist_sp4a, *dyn_norm_total_hist_sp4b, *dyn_norm_total_hist_6A;
double *dyn_norm_total_hist_sp5, *dyn_norm_total_hist_sp5a, *dyn_norm_total_hist_sp5b, *dyn_norm_total_hist_sp5c;
double *dyn_norm_total_hist_6Z, *dyn_norm_total_hist_7K;
double *dyn_norm_total_hist_8A, *dyn_norm_total_hist_8B, *dyn_norm_total_hist_8K;
double *dyn_norm_total_hist_9A, *dyn_norm_total_hist_9B, *dyn_norm_total_hist_9K;
double *dyn_norm_total_hist_10A, *dyn_norm_total_hist_10B, *dyn_norm_total_hist_10K, *dyn_norm_total_hist_10W;
double *dyn_norm_total_hist_11A, *dyn_norm_total_hist_11B, *dyn_norm_total_hist_11C, *dyn_norm_total_hist_11E, *dyn_norm_total_hist_11F, *dyn_norm_total_hist_11W;
double *dyn_norm_total_hist_12A, *dyn_norm_total_hist_12B, *dyn_norm_total_hist_12D, *dyn_norm_total_hist_12E, *dyn_norm_total_hist_12K;
double *dyn_norm_total_hist_13A, *dyn_norm_total_hist_13B, *dyn_norm_total_hist_13K;
double *dyn_norm_total_hist_FCC, *dyn_norm_total_hist_HCP, *dyn_norm_total_hist_BCC_9, *dyn_norm_total_hist_BCC_15;

double *dyn_decay_norm_total_hist_sp3, *dyn_decay_norm_total_hist_sp3a, *dyn_decay_norm_total_hist_sp3b, *dyn_decay_norm_total_hist_sp3c;
double *dyn_decay_norm_total_hist_sp4, *dyn_decay_norm_total_hist_sp4a, *dyn_decay_norm_total_hist_sp4b, *dyn_decay_norm_total_hist_6A;
double *dyn_decay_norm_total_hist_sp5, *dyn_decay_norm_total_hist_sp5a, *dyn_decay_norm_total_hist_sp5b, *dyn_decay_norm_total_hist_sp5c;
double *dyn_decay_norm_total_hist_6Z, *dyn_decay_norm_total_hist_7K;
double *dyn_decay_norm_total_hist_8A, *dyn_decay_norm_total_hist_8B, *dyn_decay_norm_total_hist_8K;
double *dyn_decay_norm_total_hist_9A, *dyn_decay_norm_total_hist_9B, *dyn_decay_norm_total_hist_9K;
double *dyn_decay_norm_total_hist_10A, *dyn_decay_norm_total_hist_10B, *dyn_decay_norm_total_hist_10K, *dyn_decay_norm_total_hist_10W;
double *dyn_decay_norm_total_hist_11A, *dyn_decay_norm_total_hist_11B, *dyn_decay_norm_total_hist_11C, *dyn_decay_norm_total_hist_11E, *dyn_decay_norm_total_hist_11F, *dyn_decay_norm_total_hist_11W;
double *dyn_decay_norm_total_hist_12A, *dyn_decay_norm_total_hist_12B, *dyn_decay_norm_total_hist_12D, *dyn_decay_norm_total_hist_12E, *dyn_decay_norm_total_hist_12K;
double *dyn_decay_norm_total_hist_13A, *dyn_decay_norm_total_hist_13B, *dyn_decay_norm_total_hist_13K;
double *dyn_decay_norm_total_hist_FCC, *dyn_decay_norm_total_hist_HCP, *dyn_decay_norm_total_hist_BCC_9, *dyn_decay_norm_total_hist_BCC_15;

int *dyn_max_hist_sp3, *dyn_max_hist_sp3a, *dyn_max_hist_sp3b, *dyn_max_hist_sp3c;
int *dyn_max_hist_sp4, *dyn_max_hist_sp4a, *dyn_max_hist_sp4b, *dyn_max_hist_6A;
int *dyn_max_hist_sp5, *dyn_max_hist_sp5a, *dyn_max_hist_sp5b, *dyn_max_hist_sp5c;
int *dyn_max_hist_6Z, *dyn_max_hist_7K;
int *dyn_max_hist_8A, *dyn_max_hist_8B, *dyn_max_hist_8K;
int *dyn_max_hist_9A, *dyn_max_hist_9B, *dyn_max_hist_9K;
int *dyn_max_hist_10A, *dyn_max_hist_10B, *dyn_max_hist_10K, *dyn_max_hist_10W;
int *dyn_max_hist_11A, *dyn_max_hist_11B, *dyn_max_hist_11C, *dyn_max_hist_11E, *dyn_max_hist_11F, *dyn_max_hist_11W;
int *dyn_max_hist_12A, *dyn_max_hist_12B, *dyn_max_hist_12D, *dyn_max_hist_12E, *dyn_max_hist_12K;
int *dyn_max_hist_13A, *dyn_max_hist_13B, *dyn_max_hist_13K;
int *dyn_max_hist_FCC, *dyn_max_hist_HCP, *dyn_max_hist_BCC_9, *dyn_max_hist_BCC_15;

double *dyn_norm_max_hist_sp3, *dyn_norm_max_hist_sp3a, *dyn_norm_max_hist_sp3b, *dyn_norm_max_hist_sp3c;
double *dyn_norm_max_hist_sp4, *dyn_norm_max_hist_sp4a, *dyn_norm_max_hist_sp4b, *dyn_norm_max_hist_6A;
double *dyn_norm_max_hist_sp5, *dyn_norm_max_hist_sp5a, *dyn_norm_max_hist_sp5b, *dyn_norm_max_hist_sp5c;
double *dyn_norm_max_hist_6Z, *dyn_norm_max_hist_7K;
double *dyn_norm_max_hist_8A, *dyn_norm_max_hist_8B, *dyn_norm_max_hist_8K;
double *dyn_norm_max_hist_9A, *dyn_norm_max_hist_9B, *dyn_norm_max_hist_9K;
double *dyn_norm_max_hist_10A, *dyn_norm_max_hist_10B, *dyn_norm_max_hist_10K, *dyn_norm_max_hist_10W;
double *dyn_norm_max_hist_11A, *dyn_norm_max_hist_11B, *dyn_norm_max_hist_11C, *dyn_norm_max_hist_11E, *dyn_norm_max_hist_11F, *dyn_norm_max_hist_11W;
double *dyn_norm_max_hist_12A, *dyn_norm_max_hist_12B, *dyn_norm_max_hist_12D, *dyn_norm_max_hist_12E, *dyn_norm_max_hist_12K;
double *dyn_norm_max_hist_13A, *dyn_norm_max_hist_13B, *dyn_norm_max_hist_13K;
double *dyn_norm_max_hist_FCC, *dyn_norm_max_hist_HCP, *dyn_norm_max_hist_BCC_9, *dyn_norm_max_hist_BCC_15;

double *dyn_norm_corrected_max_hist_sp3, *dyn_norm_corrected_max_hist_sp3a, *dyn_norm_corrected_max_hist_sp3b, *dyn_norm_corrected_max_hist_sp3c;
double *dyn_norm_corrected_max_hist_sp4, *dyn_norm_corrected_max_hist_sp4a, *dyn_norm_corrected_max_hist_sp4b, *dyn_norm_corrected_max_hist_6A;
double *dyn_norm_corrected_max_hist_sp5, *dyn_norm_corrected_max_hist_sp5a, *dyn_norm_corrected_max_hist_sp5b, *dyn_norm_corrected_max_hist_sp5c;
double *dyn_norm_corrected_max_hist_6Z, *dyn_norm_corrected_max_hist_7K;
double *dyn_norm_corrected_max_hist_8A, *dyn_norm_corrected_max_hist_8B, *dyn_norm_corrected_max_hist_8K;
double *dyn_norm_corrected_max_hist_9A, *dyn_norm_corrected_max_hist_9B, *dyn_norm_corrected_max_hist_9K;
double *dyn_norm_corrected_max_hist_10A, *dyn_norm_corrected_max_hist_10B, *dyn_norm_corrected_max_hist_10K, *dyn_norm_corrected_max_hist_10W;
double *dyn_norm_corrected_max_hist_11A, *dyn_norm_corrected_max_hist_11B, *dyn_norm_corrected_max_hist_11C, *dyn_norm_corrected_max_hist_11E, *dyn_norm_corrected_max_hist_11F, *dyn_norm_corrected_max_hist_11W;
double *dyn_norm_corrected_max_hist_12A, *dyn_norm_corrected_max_hist_12B, *dyn_norm_corrected_max_hist_12D, *dyn_norm_corrected_max_hist_12E, *dyn_norm_corrected_max_hist_12K;
double *dyn_norm_corrected_max_hist_13A, *dyn_norm_corrected_max_hist_13B, *dyn_norm_corrected_max_hist_13K;
double *dyn_norm_corrected_max_hist_FCC, *dyn_norm_corrected_max_hist_HCP, *dyn_norm_corrected_max_hist_BCC_9, *dyn_norm_corrected_max_hist_BCC_15;

double *dyn_decay_norm_max_hist_sp3, *dyn_decay_norm_max_hist_sp3a, *dyn_decay_norm_max_hist_sp3b, *dyn_decay_norm_max_hist_sp3c;
double *dyn_decay_norm_max_hist_sp4, *dyn_decay_norm_max_hist_sp4a, *dyn_decay_norm_max_hist_sp4b, *dyn_decay_norm_max_hist_6A;
double *dyn_decay_norm_max_hist_sp5, *dyn_decay_norm_max_hist_sp5a, *dyn_decay_norm_max_hist_sp5b, *dyn_decay_norm_max_hist_sp5c;
double *dyn_decay_norm_max_hist_6Z, *dyn_decay_norm_max_hist_7K;
double *dyn_decay_norm_max_hist_8A, *dyn_decay_norm_max_hist_8B, *dyn_decay_norm_max_hist_8K;
double *dyn_decay_norm_max_hist_9A, *dyn_decay_norm_max_hist_9B, *dyn_decay_norm_max_hist_9K;
double *dyn_decay_norm_max_hist_10A, *dyn_decay_norm_max_hist_10B, *dyn_decay_norm_max_hist_10K, *dyn_decay_norm_max_hist_10W;
double *dyn_decay_norm_max_hist_11A, *dyn_decay_norm_max_hist_11B, *dyn_decay_norm_max_hist_11C, *dyn_decay_norm_max_hist_11E, *dyn_decay_norm_max_hist_11F, *dyn_decay_norm_max_hist_11W;
double *dyn_decay_norm_max_hist_12A, *dyn_decay_norm_max_hist_12B, *dyn_decay_norm_max_hist_12D, *dyn_decay_norm_max_hist_12E, *dyn_decay_norm_max_hist_12K;
double *dyn_decay_norm_max_hist_13A, *dyn_decay_norm_max_hist_13B, *dyn_decay_norm_max_hist_13K;
double *dyn_decay_norm_max_hist_FCC, *dyn_decay_norm_max_hist_HCP, *dyn_decay_norm_max_hist_BCC_9, *dyn_decay_norm_max_hist_BCC_15;

double *dyn_decay_norm_corrected_max_hist_sp3, *dyn_decay_norm_corrected_max_hist_sp3a, *dyn_decay_norm_corrected_max_hist_sp3b, *dyn_decay_norm_corrected_max_hist_sp3c;
double *dyn_decay_norm_corrected_max_hist_sp4, *dyn_decay_norm_corrected_max_hist_sp4a, *dyn_decay_norm_corrected_max_hist_sp4b, *dyn_decay_norm_corrected_max_hist_6A;
double *dyn_decay_norm_corrected_max_hist_sp5, *dyn_decay_norm_corrected_max_hist_sp5a, *dyn_decay_norm_corrected_max_hist_sp5b, *dyn_decay_norm_corrected_max_hist_sp5c;
double *dyn_decay_norm_corrected_max_hist_6Z, *dyn_decay_norm_corrected_max_hist_7K;
double *dyn_decay_norm_corrected_max_hist_8A, *dyn_decay_norm_corrected_max_hist_8B, *dyn_decay_norm_corrected_max_hist_8K;
double *dyn_decay_norm_corrected_max_hist_9A, *dyn_decay_norm_corrected_max_hist_9B, *dyn_decay_norm_corrected_max_hist_9K;
double *dyn_decay_norm_corrected_max_hist_10A, *dyn_decay_norm_corrected_max_hist_10B, *dyn_decay_norm_corrected_max_hist_10K, *dyn_decay_norm_corrected_max_hist_10W;
double *dyn_decay_norm_corrected_max_hist_11A, *dyn_decay_norm_corrected_max_hist_11B, *dyn_decay_norm_corrected_max_hist_11C, *dyn_decay_norm_corrected_max_hist_11E, *dyn_decay_norm_corrected_max_hist_11F, *dyn_decay_norm_corrected_max_hist_11W;
double *dyn_decay_norm_corrected_max_hist_12A, *dyn_decay_norm_corrected_max_hist_12B, *dyn_decay_norm_corrected_max_hist_12D, *dyn_decay_norm_corrected_max_hist_12E, *dyn_decay_norm_corrected_max_hist_12K;
double *dyn_decay_norm_corrected_max_hist_13A, *dyn_decay_norm_corrected_max_hist_13B, *dyn_decay_norm_corrected_max_hist_13K;
double *dyn_decay_norm_corrected_max_hist_FCC, *dyn_decay_norm_corrected_max_hist_HCP, *dyn_decay_norm_corrected_max_hist_BCC_9, *dyn_decay_norm_corrected_max_hist_BCC_15;


#endif 