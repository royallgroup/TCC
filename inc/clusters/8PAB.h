#ifndef TCC_8PAB_H
#define TCC_8PAB_H

void Clusters_Get8PAB();
void add_8PAB(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int new_particle) ;
int check_unique8PAB(const int *old_clust, int new_particle, int new_clust_size);
#endif
