/* Alex Malins - alex.malins@gmail.com */
/* TCC: A topological cluster classification code with temporal tracking of clusters. */
/* Not for general consumption */	

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tools.h"
#include "setup.h"

double Bonds_GetR2(int i, int j);	// get separation between particles i and j
double Bonds_GetR2_PBCs(int i, int j);	// get wrapped separation between particles i and j
int icell(int tix, int tiy, int tiz);
void links();
void Bonds_GetBonds(int f);	// Get bonds using simple lengths
int Bonds_BondCheck(int i, int j);	// Returns 1 if i & j are bonded; 0 otherwise

void Rings_gSP3(int f, int n0);	// get SP3/4/5 rings including particle n0
void Rings_gSP4(int f, int n0, int n1, int n2);	// {n0,n1,n2} is not an SP3 ring, is it an SP4 or SP5 ring
void Rings_gSP5(int f, int n0, int n1, int n2, int n3);	// {n0,n1,n2,n3} is not an SP4 ring, is it an SP5 ring
void Rings_aSP3(int f, int n0, int n1, int n2);	// Take {n0,n1,n2}, check SP3 ring and if so detect SP3a/b/c cluster
void Rings_aSP4(int f, int n0, int n1, int n2, int n3);	// Take {n0,n1,n2,n3}, check SP4 ring and if so detect SP4a/b/c cluster
void Rings_aSP5(int f, int n0, int n1, int n2, int n3, int n4);	// Take {n0,n1,n2,n3,n4}, check SP5 ring and if so detect SP5a/b/c cluster
void Rings_setSP3c(int f); // store cluster 5A D3h from Bonds_aSP3()
void Rings_setSP4c(int f);	// store cluster 6A Oh from Bonds_aSP4()
void Rings_setSP5c(int f);	// store cluster 7A D5h from Bonds_aSP5()

void Clusters_Get6Z_C2v(int f);	// Detect 6Z clusters
void Clusters_Get7K(int f);	// Detect 7K clusters
void Clusters_Get8A_D2d(int f);	// Detect 8A D2d clusters
void Clusters_Get8B_Cs(int f);	// Detect 8B Cs clusters
void Clusters_Get8K(int f);	// Detect 8K clusters
void Clusters_Get9A_D3h(int f);	// Detect 9A D3h clusters
void Clusters_Get9K_10K(int f);	// Detect 9K clusters
int Clusters_Get10K(int f, char *ach, char *ach_cen, char *ach_shell);	// Detect 10K clusters
void Clusters_Get10W(int f);	// Detect 10W clusters

void Clusters_Get9B_10B_11B_11E_12D(int f);	// Detect 9B 10B 11B 11E 12D clusters together
void Clusters_Get10B_C3v(int f, int i, int j, char *ach, char *ach_cen, char *ach_shell, char *ach1, char *ach1_cen, char *ach1_shell);	// Return 1 if 9B is also 10B cluster
int Clusters_Get11B_C2v(int f, char *ach, char *ach_cen, char *ach_shell);	// Detect 11B C2v clusters
void Clusters_Get11E_12D(int f, int i, int j, int sp1, int sp2i, int sp2j, char *ach1, char *ach2);	// Return 11Es from 9B
int Clusters_Get11W_Cs(int f, char *ach, char *ach_cen, char *ach_shell);	// Detect 11W C2s clusters, return 1 if 10B is an 11W
int Clusters_Get12D_D2d(int f, int i, int j, int k, int sp1, int sp2, char *ach1);	// Return 1 if 12B is also 11E


void Clusters_Get10A_C3v(int f);	// Detect 10A C3v clusters
void Clusters_Get11A_12K(int f);	// Detect 11A D4d clusters
int Clusters_Get12K(int f, int SP3_1, int SP3_2, int SP3_3, char *ach, char *ach_cen, char *ach_shell);	// Detect 12K clusters

void Clusters_Get11C_12A(int f);	// Detect 11C Cs & 12A C2v clusters
int Clusters_Get12A_C2v(int f, char *ach, char *ach_cen, char *ach_shell);	// Detect 12A C2v cluster

void Clusters_Get11F_12E_13K(int f);	// Detect 11F C2v & 12E D3h & 13K
int Clusters_Get12E_D3h(int f, int j, char *ach);	// Return 1 is 11F is also 12E
int Clusters_Get13K(int f, int sp3c_i, int sp3c_j, int the6A_i, char *ach, char *ach_cen, char *ach_shell);	// Detect 13K clusters


void Clusters_Get12B_13A(int f);	// Detect 12B & 13A D5h clusters together

void Clusters_Get13B_D5h(int f);	// Detect 13B D5h clusters, i.e. twisted icosahedra

void Clusters_GetFCC(int f);	// Detect 13 particle FCC clusters
void Clusters_GetHCP(int f);	// Detect 13 particle HCP clusters
void Clusters_GetBCC_9(int f);	// Detect 9 particle BCC clusters
void Clusters_GetBCC_15(int f);	// Detect 15 particle BCC clusters

void Stats_Init();	// initialize Stats routine
void Stats_Reset();	// reset Cluster routine variables
void Stats_FreeMem();	// free memory from stats variables
void Stats_Analyse();	// output Cluster statistics to file
void Stats_SetA();	// Set arrays to true if the ith particle is a member of any clusters with this or a larger number of particles

