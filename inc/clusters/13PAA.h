#ifndef TCC_13PAA_H
#define TCC_13PAA_H

void Clusters_Get13PAA();
void add_13PAA(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int final_part2, int final_part3, int final_part4, int final_part5, int final_part6, int new_particle) ;
int check_unique13PAA(const int *old_clust, int new_particle, int new_clust_size);
#endif
