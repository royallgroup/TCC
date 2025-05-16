#ifndef TCC_8PAA_H
#define TCC_8PAA_H

void Clusters_Get8PAA();
void add_8PAA(int bonded_id, int nbonded_id, int bonded_id_2, int nbonded_id_2, int spindle_id, int spindle_id_2, int final_part, int new_particle);
int get_new_particle_8PAA(const int *new_5A_cluster, int spindle_id); 
int crossover_part(const int *old_7T_clust, int new_part, int clust_1_size);
int crossover_spiral(const int *clust_1, const int *clust_2, int clust_1_size);
int check_unique8PAA(const int *old_7T, int new_particle, int new_cluster_size);
void get_final_part(int *final_parts,const int *old_7T, int bonded_id, int nbonded_id, int bonded_id_2, int nbonded_id_2, int spindle_id, int spindle_id_2, int clust_size, int nparts);
#endif
