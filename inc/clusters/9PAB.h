#ifndef TCC_9PAB_H
#define TCC_9PAB_H

void Clusters_Get9PAB();
void add_9PAB(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int final_part2, int new_particle) ;
int check_unique9PAB(const int *old_clust, int new_particle, int new_clust_size);
#endif
