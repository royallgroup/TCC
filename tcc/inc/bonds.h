#ifndef BONDS_H
#define BONDS_H 

double Get_Interparticle_Distance(int i, int j);

double Get_Interparticle_Distance_With_PBCs(int i, int j);

void Are_All_Bonds_Symmetric();

void Get_Bonds();

void Get_Simple_Bonds();

void Too_Many_Bonds(int particle_1, int particle_2);

void Add_New_Bond(int particle_1, int particle_2, double squared_distance);

void Check_For_Valid_Bond(int particle_1, int particle_2, double squared_distance);

void Check_Num_Bonds(int particle_1, int particle_2, double squared_distance);

void Get_Bonds_With_Voronoi();

int Bonds_BondCheck(int , int );

#endif