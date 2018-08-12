#ifndef TCC_10W_H
#define TCC_10W_H

void Clusters_Get10W();

int get_shell_particle_ids(int *shell_ids, const int *neighbouring_sp5_ids);

void Cluster_Write_10W(int center_id, int *shell_ids);

#endif
