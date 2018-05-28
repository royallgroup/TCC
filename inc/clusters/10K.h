#ifndef TCC_10K_H
#define TCC_10K_H

void Clusters_Get10K();

int count_extra_particles(const int *first_9K_cluster, int first_9K_common_id, int *extra_particle_id);

void Cluster_Write_10K(int id_9K, int ep);

#endif
