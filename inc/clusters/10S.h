#ifndef TCC_10S_H
#define TCC_10S_H

void Clusters_Get10S();
void add_10S(int **clust);
void get_new_rings(int *new_ring_p, const int* clust1, const int *clust2, int r1, int r2);
void get_10S(int **array, const int *clust1, const int *clust2);
#endif
