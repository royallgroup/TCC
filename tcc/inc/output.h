#ifndef TCC_OUTPUT_H
#define TCC_OUTPUT_H
#include <stdio.h>

void Write_Raw(int f);
void Write_Raw_Particle_Types(int f, FILE *thefile, const char *sarr);
void Bonds_WriteBonds(int f);
void Write_Cluster_Centers_xyz(int f, int cluster_type);
void Write_Cluster(int f);
void Write_Cluster_Compostions(int f, int n, int **hc, int clusSize, int cluster_number);
void Write_Pop_Per_Frame(int f);

#endif //TCC_OUTPUT_H
