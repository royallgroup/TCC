#ifndef TCC_12PAA_H
#define TCC_12PAA_H

void Clusters_Get12PAA();
void add_12PAA(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int final_part2, int final_part3, int final_part4, int final_part5, int new_particle) ;
int check_unique12PAA(const int *old_clust, int new_particle, int new_clust_size);
#endif
