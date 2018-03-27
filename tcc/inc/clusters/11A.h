#ifndef TCC_11A_H
#define TCC_11A_H

void Clusters_Get11A();

int Check_unique_6A_rings(int first_6A_id, int second_6A_id);

int Check_6A_rings_bonded(int first_6A_id, int second_6A_id);

void Cluster_Write_11A(int first_6A_id, int second_6A_id, const int sother[], int scom);

#endif
