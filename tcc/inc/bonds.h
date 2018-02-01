#ifndef BONDS_H
#define BONDS_H 

#include "globals.h"
#include "tools.h"
//// START: Bonds routines
double Bonds_GetR2(int i, int j) ;  // get separation between particles i and j

double Bonds_GetR2_PBCs(int i, int j) ; // get PBC wrapped separation between particles i and j

void Bonds_WriteBonds(int f) ;

void Bonds_TickBLDistro(double , int *, int *) ;

void Bonds_WriteBLDistro(char *, int *, int *, double *) ;

void Bonds_WriteBLDistroClust(char *, int *, int *,double *) ;

void Bonds_TallynbDistro();

void Bonds_WritenbDistro(char *, int *, int *, double *) ;

void Bonds_CheckSymmetric();

void Bonds_GetBonds(int ) ;    // Get bonds using simple lengths

void Bonds_GetBondsV(int )  ;  // Get bonds using Voronoi

void Bonds_GetBondsV_CellList(int ) ; // Get bonds using Voronoi


int Bonds_BondCheck(int , int ) ; // Returns 1 if i & j are bonded; 0 otherwise

int Bonds_cnb_j(int , int ) ; // Returns number k of if bNums[i][k]=j; -1 if i and j unbonded

#endif