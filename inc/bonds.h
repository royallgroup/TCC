#ifndef BONDS_H
#define BONDS_H 

double Get_Interparticle_Distance(int i, int j);

void get_distance_components(int i, int j, double *dx, double *dy, double *dz);

void enforce_PBCs(double *dx, double *dy, double *dz);

void Are_All_Bonds_Symmetric();

void build_bond_network(int frame_number);

void Get_Simple_Bonds();

void too_many_bonds(int particle_1, int particle_2, const char *method_name);

void Add_New_Bond(int particle_1, int particle_2, double squared_distance);

void Check_For_Valid_Bond(int particle_1, int particle_2, double squared_distance);

void Check_Num_Bonds(int particle_1, int particle_2, double squared_distance);

void Get_Bonds_With_Voronoi();

int Bonds_BondCheck(int , int );

#endif