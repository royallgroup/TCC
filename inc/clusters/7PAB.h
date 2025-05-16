#ifndef TCC_7PAB_H
#define TCC_7PAB_H

void Clusters_Get7PAB();
int nring_bonded_nbonded(int bonded_id, int nbonded_id, const int *clust_2);
int nring_in_cluster(const int *clust_1, const int *clust_2, int clust_1_size);
void add_7PAB(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int new_particle) ;
int get_new_particle_P2(const int *new_3_ring, int bonded_id, int nbonded_id);
int check_unique7PAB(const int *old_clust, int new_particle, int new_clust_size);
int check_5A(const int *arr_6Z, int new_part);
#endif
