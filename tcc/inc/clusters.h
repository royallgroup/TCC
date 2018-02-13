#ifndef CLUSTERS_H
#define CLUSTERS_H

void Clusters_Get6Z_C2v(int f);
void Cluster_Write_6Z(int f);

void Clusters_Get7K(int f);
void Cluster_Write_7K(int f);

void Clusters_Get8A_D2d(int f);
void Cluster_Write_8A(int f);

void Clusters_Get8B_Cs(int f);
void Clusters_8B_loop(int f, int i, int clusSize, int s1, int s2);
void Cluster_Write_8B(int f);

void Clusters_Get8K(int f);
void Cluster_Write_8K(int f);

void Clusters_Get9A_D3h(int f);
void Cluster_Write_9A(int f);

void Clusters_Get9B_10B_11B_11E_12D(int f);
void Cluster_Write_9B(int f);

void Clusters_Get10B_C3v(int f, int i, int j);
void Cluster_Write_10B(int f);

int Clusters_Get11B_C2v(int f);
void Cluster_Write_11B(int f);

void Clusters_Get11W(int f);
int get_11W_extra_particle(int id_10B, int spindle_10B);
int is_particle_in_10B(int particle_id, int id_10B);
int is_particle_bonded_to_7As(int id_10B, int extra_particle);
void resize_hc11W(int f);
void populate_hc11W(int f, int id_10B, int extra_particle);
void populate_s11W(int f);

void Clusters_Get11E_12D(int f, int i, int j, int sp1, int sp2i, int sp2j);
void Clust_Write_11E(int f);

int Clusters_Get12D_D2d(int f, int j, int k, int sp1, int sp2);
void Cluster_Write_12D(int f);

void Clusters_Get9K(int f);
void Cluster_Write_9k(int f, const int trial[]);

void Clusters_Get10K(int f);
int is_particle_in_9K(int id_9k, int id_particle);
void Cluster_Write_10K(int f, int id_9K, int ep);

void Clusters_Get10A_C3v(int f);
void Cluster_Write_10A(int f);

void Clusters_Get10W(int f);
void Cluster_Write_10W(int f);

void Clusters_Get11A(int f) ;
int Check_unique_6A_rings(int first_6A_id, int second_6A_id);
int Check_6A_rings_bonded(int first_6A_id, int second_6A_id);
void Cluster_Write_11A(int f, int first_6A_id, int second_6A_id, const int sother[], int scom);

void Clusters_Get12K(int f);
void get_12K_ring_bonds(int ptr_11A, int (*sp3_rings)[3]);
void find_12K_cluster(int f, int ptr_11A, const int *sp3_ring);
int is_particle_in_11A(int ptr_11A, int particle_id);
void Cluster_Write_12K(int f, int ep, int id_11A);

void Clusters_Get11C(int f);
void resize_hc11C(int f);
void Cluster_Write_11C(int f);

int Clusters_Get12A(int f);
int get_12A_extra_particle(int id_11C);
int is_particle_in_11C(int particle_id, int id_11C);
int bond_check_12A_extra_particle(int id_11C, int ep);
void resize_hc12A(int f);
void populate_hc12A(int f, int id_11C, int ep);
void populate_s12A(int f);

void Clusters_Get11F_12E_13K(int f);
void Cluster_Write_11F(int f);

int Clusters_Get12E_D3h(int f, int j);
void Cluster_Write_12E(int f);

int Clusters_Get13K(int f, int sp3c_i, int sp3c_j, int the6A_i);
void Cluster_Write_13K(int f);

void Clusters_Get12B_13A(int f);
void Clust_Write_12B(int f);
void Clust_Write_13A(int f);

void Clusters_Get13B_D5h(int f);
void Cluster_Write_13B(int f);

void Clusters_GetFCC(int f);
void Cluster_Write_FCC(int f);

void Clusters_GetHCP(int f);
void Cluster_Write_HCP(int f, int i, int j, int j2, int k);

void Clusters_GetBCC_9(int f);
void Cluster_Write_BCC9(int f);

void Clusters_GetBCC_15(int f);
void Cluster_Write_BCC_15(int f, int clusSize);

#endif