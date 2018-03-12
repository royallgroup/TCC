#ifndef CLUSTERS_H
#define CLUSTERS_H

void Clusters_Get6Z();
void Cluster_Write_6Z();

void Clusters_Get7K();
void Cluster_Write_7K();

void Clusters_Get8A();
void Cluster_Write_8A();

void Clusters_Get8B();
void Clusters_8B_loop(int i, int clusSize, int s1, int s2);
void Cluster_Write_8B();

void Clusters_Get8K();
void Cluster_Write_8K();

void Clusters_Get9A();
void Cluster_Write_9A();

void Clusters_Get9B_10B_11B_11E_12D();
void Cluster_Write_9B();

void Clusters_Get10B(int i, int j);
void Cluster_Write_10B();

int Clusters_Get11B();
void Cluster_Write_11B();

void Clusters_Get11W();
int get_11W_extra_particle(int id_10B, int spindle_10B);
int is_particle_in_10B(int particle_id, int id_10B);
int is_particle_bonded_to_7As(int id_10B, int extra_particle);
void resize_hc11W();
void populate_hc11W(int id_10B, int extra_particle);
void populate_s11W();

void Clusters_Get11E_12D(int i, int j, int sp1, int sp2i, int sp2j);
void Clust_Write_11E();

int Clusters_Get12D(int j, int k, int sp1, int sp2);
void Cluster_Write_12D();

void Clusters_Get9K();
void Cluster_Write_9k(const int trial[]);

void Clusters_Get10K();
int is_particle_in_9K(int id_9k, int id_particle);
void Cluster_Write_10K(int id_9K, int ep);

void Clusters_Get10A();
void Cluster_Write_10A();

void Clusters_Get10W();
void Cluster_Write_10W();

void Clusters_Get11A();
int Check_unique_6A_rings(int first_6A_id, int second_6A_id);
int Check_6A_rings_bonded(int first_6A_id, int second_6A_id);
void Cluster_Write_11A(int first_6A_id, int second_6A_id, const int sother[], int scom);

void Clusters_Get12K();
void get_12K_ring_bonds(int ptr_11A, int (*sp3_rings)[3]);
void find_12K_cluster(int ptr_11A, const int *sp3_ring);
int is_particle_in_11A(int ptr_11A, int particle_id);
void Cluster_Write_12K(int ep, int id_11A);

void Clusters_Get11C();
int get_11C_spindle_particles(int *uncommon_spindle, int id_first_7A, int id_second7A, int *common_spindle);
int get_bonded_7A_ring_particles(int *ar, int id_first_7A, int id_second7A, int ncom);
void resize_hc11C();
void Cluster_Write_11C();

int Clusters_Get12A();
int get_12A_extra_particle(int id_11C);
int is_particle_in_11C(int particle_id, int id_11C);
int bond_check_12A_extra_particle(int id_11C, int ep);
void resize_hc12A();
void populate_hc12A(int id_11C, int ep);
void populate_s12A();

void Clusters_Get11F_12E_13K();
void Cluster_Write_11F();

int Clusters_Get12E(int j);
void Cluster_Write_12E();

int Clusters_Get13K(int sp3c_i, int sp3c_j, int the6A_i);
void Cluster_Write_13K();

void Clusters_Get12B_13A();
void Clust_Write_12B();
void Clust_Write_13A();

void Clusters_Get13B();
void Cluster_Write_13B();

void Clusters_GetFCC();
void Cluster_Write_FCC();

void Clusters_GetHCP();
void Cluster_Write_HCP(int i, int j, int j2, int k);

void Clusters_GetBCC_9();
void Cluster_Write_BCC9();

void Clusters_GetBCC_15();
void Cluster_Write_BCC_15(int clusSize);

#endif