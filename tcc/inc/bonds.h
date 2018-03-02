#ifndef BONDS_H
#define BONDS_H 

double Get_Interparticle_Distance(int i, int j) ;  // get separation between particles i and j

double Get_Interparticle_Distance_With_PBCs(int i, int j) ; // get PBC wrapped separation between particles i and j

void Are_All_Bonds_Symmetric();

void Get_Bonds();    // Get bonds using simple lengths

void Get_Simple_Bonds();

void Bonds_GetBondsV()  ;  // Get bonds using Voronoi

int Bonds_BondCheck(int , int ) ; // Returns 1 if i & j are bonded; 0 otherwise

#endif