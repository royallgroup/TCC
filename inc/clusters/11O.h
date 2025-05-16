#ifndef TCC_11O_H
#define TCC_11O_H

void Clusters_Get11O();
int check_unique_11O(int **arr_11O);
void add_11O(int **arr_11O);
int overlap_6A_4A_11O(int **array, int *clust_6A, int *clust_4A);
void get_11O(int ** arr_12S, int *first_6A_cluster, int *second_6A_cluster, int *first_4A_cluster, int *second_4A_cluster);
int overlap_6A_6A_11O(int **array, int *clust_6A1, int *clust_6A2);
#endif