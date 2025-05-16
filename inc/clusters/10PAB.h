#ifndef TCC_10PAB_H
#define TCC_10PAB_H

void Clusters_Get10PAB();
void add_10PAB(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int final_part2, int final_part3, int new_particle) ;
int check_unique10PAB(const int *old_clust, int new_particle, int new_clust_size);
#endif
