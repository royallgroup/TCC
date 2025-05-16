#ifndef TCC_12O_H
#define TCC_12O_H

void Clusters_Get12O();
int check_unique_12O(int **new_12O_cluster);
void add_12O(int **new_12O_cluster);
void get_12O(int **array, int *first_4A_cluster, int *second_4A_cluster, int *first_6A_cluster, int *second_6A_cluster);
int overlap_6A_ring(int *first_6A_cluster, int *second_6A_cluster);
int overlap_6A_4A(int *clust_6A, int *clust_4A);
int shared_4A(int *third_4A_cluster, int *fourth_4A_cluster, int *first_4A_cluster, int *second_4A_cluster);
int overlap_4A_4A(int *clust_4A1, int *clust_4A2);
int overlap_6A_6A_12O(int **array, int *clust_6A1, int *clust_6A2);
int overlap_6A_6A_4A_4A(int **array,int *clust_6A1, int *clust_4A1, int *clust_6A2, int *clust_4A2);
int overlap_4A_12O(int *clust_4A, int *clust_6A1, int *clust_6A2);
#endif
