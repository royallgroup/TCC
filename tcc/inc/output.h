#ifndef TCC_OUTPUT_H
#define TCC_OUTPUT_H
#include <stdio.h>

void Write_Raw_Init();
void Write_Raw_Xmol(int f, FILE *thefile, const char *sarr);
void Write_11A_cen_xmol(int f);
void Write_13A_cen_xmol(int f);
void Write_Raw(int f);
void Write_Raw_Close();
void Write_Cluster_Init();
void Write_Cluster_Close();
void Write_Cluster_Xmol(int f, FILE *writeout, int n, int **hc, int clusSize);
void Write_Cluster(int f);
void Write_Cluster_sp3(int f, FILE *writeout);
void Write_Cluster_sp4(int f, FILE *writeout);
void Write_Cluster_sp5(int f, FILE *writeout);
void Write_Pop_Per_Frame(int f);

#endif //TCC_OUTPUT_H
