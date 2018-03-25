#ifndef TCC_11C_H
#define TCC_11C_H

void Clusters_Get11C();

int get_11C_spindle_particles(int *uncommon_spindle, int id_first_7A, int id_second7A, int *common_spindle);

int get_bonded_7A_ring_particles(int *ar, int id_first_7A, int id_second7A, int ncom);

void resize_hc11C();

void Cluster_Write_11C();

#endif