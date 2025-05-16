#ifndef TCC_11S_H
#define TCC_11S_H

void Clusters_Get11S();
void add_11S(int **clust);
int check_unique_11S(int **clust);
void get_11S(const int *clust_7A, const int *clust_9B, int **clust_11S);
int overlap_7A_9B(const int *clust_7A,const int *clust_9B);
int i_in_ring(const int *clust_7A, int i);
int i_in_spindle(const int *clust_7A, int i);
int common_spindle(const int *clust1, const int *clust2, int r1, int r2);
int common_ring(const int *clust1, const int *clust2, int r1, int r2);
#endif
