/* Ring statistics and cluster analysis program in C */
/* by Stephen R. Williams, RSC, The Australian National University, Australia */
/* + Alex Malins, Bristol Centre for Complexity Sciences, University of Bristol, United Kingdom */
/* Early 2010 */
/* Not for general consumption */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//// START: Global variables
FILE *fXmol, *fBonds;
char *fXmolName;
char fInputParamsName[1000], fBondsName[1000], fDynamicDatName[1000];	// name of parameters file and r... coordinates file and memsize file
double ***r;	// positions in x y and z directions of N particles
char *rtype; // particle type
int ***bNums;	// list of particles (indices j) bonded to particle at index i
double **bondlengths;	// length of bonds in the bond network and squared
int **cnb; // Current Number of Bonds for particles {1,...,N}

int doDynamicAnalysis, go_away_for, doSplits, doDynFileRewrite, doWriteNewClustsFile, steady_Nclus, doWriteLifetimeDistro, t_h, doDisplacementDistro, doWriteIn, doWriteCtsIn;
int doWriteMSDClusNonClusCts, max_regions, doRegions;
int PRINTINFO; // print running information about progress

int start_from, end_by, useable_frames;
int percent_start, percent_end, useable_percent; // what frame to start and end calculation of pop_per_frame

// from TCC inputparameters.in file
int ISNOTCUBIC;
char fBoxSizeName[1000]; //NPT stuff: name of file which contains info on box
int FRAMES, STARTFROM, SAMPLEFREQ, doDynamics, doSubClusts;
double talpha, fc; // alpha relaxtion time
int Vor, PBCs;	// 0 do not impliment periodic boundary conditions, 1 implement periodic boundary conditions
double rcutAA, rcutAB, rcutBB, binWidth;
int nB;	// max number of bonds per particle

// from .params file
int N, NA, doBinary, TOTALFRAMES;
double rho, TSTART, FRAMETSTEP, TFINAL;

double sidex, halfSidex;	// box side length
double sidey, halfSidey;	// box side length
double sidez, halfSidez;	// box side length

/* legacy outputting single cluster
int **trial_neighbours, *trial_no_neighbours, max_no_neighbours, *exists;

int output_clust, output_neighbours, output_all, output_bonds, clust_nB;
double radA, radB, reduce, sphere_size, bond_thickness; // alpha relaxtion time
int jmol_sphere_size, jmol_bond_thickness;*/

int **temp_parts, *temp_parts_A, *temp_parts_A_cen, *temp_parts_A_shell, **temp_events, **temp_sub, **temp_neighbours, **no_temp_neighbours, **temp_splits;
double **temp_lifetimes;
int calc_maxlength, *calc_histo;
int **the_clusts;
int **the_clusts_mixture;
int **clustered_parts;

double range2, thebins;
int nobins;

int *no_regions, no_regions_noPBCs;
int *used_part, **used_bond, **bond_type_x, **bond_type_y, **bond_type_z;
int *used_part_noPBCs, **used_bond_noPBCs;
int **part_move;
int **regions, *regions_noPBCs;
double **regions_Rg, **regions_mean_cluster_lifetime;
FILE *f_region_cluster_lifetimes;

int dyn_msp3, dyn_msp3a, dyn_msp3b, dyn_msp3c;
int dyn_msp4, dyn_msp4a, dyn_msp4b, dyn_m6A;
int dyn_msp5, dyn_msp5a, dyn_msp5b, dyn_msp5c;
int dyn_m6Z, dyn_m7K;
int dyn_m8A, dyn_m8B, dyn_m8K;
int dyn_m9A, dyn_m9B, dyn_m9K;
int dyn_m10A, dyn_m10B, dyn_m10K, dyn_m10W;
int dyn_m11A, dyn_m11B, dyn_m11C, dyn_m11E, dyn_m11F, dyn_m11W;
int dyn_m12A, dyn_m12B, dyn_m12D, dyn_m12E, dyn_m12K;
int dyn_m13A, dyn_m13B, dyn_m13K;
int dyn_mFCC, dyn_mHCP, dyn_mBCC_9, dyn_mBCC_15;

int mClust;

FILE *bondsin, *bondsout;

//// END: Global variables

//// START: routines
void Error(char *msg);
int Setup_GetFirstIntFromLine(FILE *stream);
double Setup_GetFirstDoubleFromLine(FILE *stream);
void Setup_ReadIniFile(char *filename);	// initialize RingStat routine
void Setup_ReadBox(FILE *); // read in box size file
void Setup_ReadXmol();	// output single float valued arrays in gopenmol xyz format
void ReadXmol(int e_frame, int f_frame);
void Setup_ReadBonds();
void Setup_InitVars();
void Setup_ResetVars();
void Setup_FreeVars();
void Input_13A(int clustSize);
void Take_Clust(int clustSize, char *strClustType, int takeClust);
void Take_Clust_Neighbours(int clustSize, char *strClustType, int takeClust);
void Take_All(int clustSize, char *strClustType);
int Bonds_BondCheck(int f, int i, int j);
void Write_ClustXmol_Parts(int *arr, int f, int elements_clust, int elements_neighbours, char *filename, FILE *xmolfile, FILE *vmdfile, FILE *jmolfile);
void Write_ClustXmol_Bonds(int *arr, int f, int elements_clust, int elements_neighbours, FILE *vmdfile, FILE *jmolfile);
void Write_WholeClustXmol(int *arr, int f, char *filename, FILE *xmolfile, FILE *vmdfile, FILE *jmolfile);
//// END: routines


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!BEGIN DOUBLE QUICKSORT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//  quickSort
//
//  This public-domain C implementation by Darel Rex Finley.
//
//  * Returns YES if sort was successful, or NO if the nested
//    pivots went too deep, in which case your array will have
//    been re-ordered, but probably not sorted correctly.
//
//  * This function assumes it is called with valid parameters.
//
//  * Example calls:
//    quickSort(&myArray[0],5); // sorts elements 0, 1, 2, 3, and 4
//    quickSort(&myArray[3],5); // sorts elements 3, 4, 5, 6, and 7

int doubleQuickSort(double *arr, int elements) {
	double piv;
	int beg[500], end[500], i, L, R;

	i=0;
	beg[0]=0;
	end[0]=elements;
	while (i>=0) {
		L=beg[i]; R=end[i]-1;
		if (L<R) {
			piv=arr[L];
			if (i==500-1) return 0;
			while (L<R) {
				while (arr[R]>=piv && L<R) R--;
				if (L<R) arr[L++]=arr[R];
				while (arr[L]<=piv && L<R) L++;
				if (L<R) arr[R--]=arr[L];
			}
			arr[L]=piv;
			beg[i+1]=L+1;
			end[i+1]=end[i];
			end[i++]=L;
		}
		else i--;
	}
	return 1;
}
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END QUICKSORT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
