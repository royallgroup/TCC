#ifndef TCC_8PBB_H
#define TCC_8PBB_H

void Clusters_Get8PBB();
void add_8PBB(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int new_particle) ;
int check_unique8PBB(const int *old_clust, int new_particle, int new_clust_size);
#endif