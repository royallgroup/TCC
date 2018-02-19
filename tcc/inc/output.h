#ifndef TCC_OUTPUT_H
#define TCC_OUTPUT_H
#include <stdio.h>

void Write_Raw_Xmol(int f, FILE *thefile, const char *sarr);
void Write_Cluster_Centers_xyz(int f, int cluster_type);
void Write_Raw(int f);
void Write_Cluster_Xmol(int f, int n, int **hc, int clusSize, int cluster_number);
void Write_Cluster(int f);
void Write_Cluster_sp3(int f, int cluster_number);
void Write_Cluster_sp4(int f, int cluster_number);
void Write_Cluster_sp5(int f, int cluster_number);
void Write_Pop_Per_Frame(int f);

#endif //TCC_OUTPUT_H
