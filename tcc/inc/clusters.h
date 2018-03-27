#ifndef CLUSTERS_H
#define CLUSTERS_H


void Clusters_Get9B_10B_11B_11E_12D();
void Cluster_Write_9B();

void Clusters_Get10B(int i, int j);
void Cluster_Write_10B();

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

int Clusters_Get13K(int sp3c_i, int sp3c_j, int the6A_i);
void Cluster_Write_13K();

void Clusters_Get12B_13A();
void Clust_Write_12B();
void Clust_Write_13A();

void Clusters_Get13B();
void Cluster_Write_13B();

#endif