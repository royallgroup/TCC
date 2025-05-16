#ifndef TCC_13SB_H
#define TCC_13SB_H

void Clusters_Get13SB();
int check_5A_13SB(int spindle, int spindle2, int *clust_11SB);
int overlap_13SB(int *clust_11SB, int *clust_5A);
int overlap_13SB_sp3a(int *clust_11SB, int *clust_sp3a,int spindle);
int check_unique_13SB(int *clust11SB, int s1, int s2);
void add_13SB(int *clust11SB, int spindle1, int spindle2);
#endif