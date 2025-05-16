#ifndef TCC_10PAA_H
#define TCC_10PAA_H

void Clusters_Get10PAA();
void add_10PAA(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int final_part2, int final_part3, int new_particle) ;
int check_unique10PAA(const int *old_clust, int new_particle, int new_clust_size);
#endif
