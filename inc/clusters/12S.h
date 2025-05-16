#ifndef TCC_12S_H
#define TCC_12S_H

void Clusters_Get12S();
void add_12S(int **clust_12S);
int check_unique_12S(int **clust_12S);
int common_ring_spindle(const int *clust1, const int *clust2, const int *clust3);
//void get_12S(const int *clust_7A, const int *clust_9B, int **clust_12S);
int common_spindle_ring(const int *clust1, const int *clust2, const int *clust3);
void get_12S2(int **clust_12S, const int *cluster1, const int *cluster2, const int *cluster3);
void get_12S(const int *first_clust_7A, const int *second_clust_7A, const int *third_clust_7A, int **clust_12S);
int overlap(const int *cluster1, const int *cluster2, int l1, int l2);
int overlap_4A_12O(int *clust_4A, int *clust_6A1, int *clust_6A2);
#endif