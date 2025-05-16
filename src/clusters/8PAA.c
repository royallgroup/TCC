#include <globals.h>
#include <tools.h>
#include "8PAA.h"
#include "7T.h"
//!  An 8PAA is a 7T with an additional particle
/*!
* A 7T cluster contains two 6Z clusters at the top and bottom of the cluster. 
* An 8PAA takes these two 6Z clusters and tries to make a 7T from the cluster and an extra particle.
*  Cluster output: BBBBBBBB
*  Storage order: original 7T particles x 6, new 5A spindle)
*/

void Clusters_Get8PAA() {
    for (int old_7T_id = 0; old_7T_id < n7T_a; ++old_7T_id) {
    int *old_7T_cluster = hc7T_a[old_7T_id];
    int clust_6Z[7]; 
    clust_6Z[0]= old_7T_cluster[2];
    clust_6Z[1] = old_7T_cluster[3];
    clust_6Z[2]= old_7T_cluster[0];
    clust_6Z[3] = old_7T_cluster[1];
    clust_6Z[4] = old_7T_cluster[5];
    clust_6Z[5] = old_7T_cluster[4];
    clust_6Z[6] = old_7T_cluster[6];
        int spindle_id = clust_6Z[5];
        int spindle_id_2 = clust_6Z[4];
        int nbonded_id = clust_6Z[3];
        int bonded_id = clust_6Z[0]; 
        int nbonded_id_2 = clust_6Z[1];
        int bonded_id_2 = clust_6Z[2]; 
        int final_part1 = clust_6Z[6]; 
        for (int new_5A_pointer = 0; new_5A_pointer < nmem_sp3c[spindle_id]; ++new_5A_pointer) {
            //int new_5A_id = mem_sp3c[spindle_id][new_5A_pointer];
            int *new_5A_cluster = hcsp3c[new_5A_pointer];
            if (is_particle_spindle_of_5A(spindle_id_2, new_5A_cluster) == 1) {
                if(is_particle_ring_of_5A(bonded_id,new_5A_cluster) == 1){
                    if(is_particle_ring_of_5A(nbonded_id,new_5A_cluster) == 1){
                        if(is_particle_ring_of_5A(spindle_id,new_5A_cluster) == 1){
                            int new_particle_id = get_new_particle(new_5A_cluster, spindle_id_2);
                                int arr_6Z[6] = {bonded_id, nbonded_id_2, bonded_id_2, nbonded_id, spindle_id, spindle_id_2};
                                int bond_counter = check_ring_bonds(new_5A_cluster, arr_6Z);
                                    if (bond_counter == 25 ){
                                        if (check_unique8PAA(old_7T_cluster, new_particle_id,8) == 0){
                                        add_8PAA(bonded_id, nbonded_id_2, bonded_id_2, nbonded_id,spindle_id, spindle_id_2, final_part1, new_particle_id);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
} 


// Given a 6Z and an 8PAA, returns the number of common particles
int crossover_spiral(const int *clust_1, const int *clust_2, int clust_1_size){
    int s1 = clust_1_size;
    int cross_count =0; 
    for(int m=0; m<6; ++m){
        for(int n=0; n<s1; ++n){
        if(clust_1[n]==clust_2[m]){
            cross_count += 1;
        }
    }
    }
    return cross_count;
}

//Given an 7T and an extra particle, returns 1 if the extra particle is already in the cluster
int crossover_part(const int *clust_1, int particle, int clust_1_size){
    int s1 = clust_1_size;
    int cross_count =0; 
    for(int n=0; n<s1; ++n){
        if(clust_1[n]==particle){
            cross_count += 1;
        }
    }
    return cross_count;
}

void add_8PAA(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part, int new_particle) {
    int clusSize = 8;
    if (n8PAA == m8PAA) {
        hc8PAA = resize_2D_int(hc8PAA, m8PAA, m8PAA + incrStatic, clusSize, -1);
        m8PAA = m8PAA + incrStatic;
    }
    hc8PAA[n8PAA][0] = bonded_id;
    hc8PAA[n8PAA][1] = nbonded_id_2;
    hc8PAA[n8PAA][2] = bonded_id_2;
    hc8PAA[n8PAA][3] = nbonded_id;
    hc8PAA[n8PAA][4] = spindle_id;
    hc8PAA[n8PAA][5] = spindle_id_2;
    hc8PAA[n8PAA][6] = final_part;
    hc8PAA[n8PAA][7] = new_particle;

    for (int i = 0; i < 8; ++i) {
        s8PAA[hc8PAA[n8PAA][i]] = 'B';
    }
    ++n8PAA;
}

int check_ring_bonds_8PAA(const int *new_5A_cluster, const int *old_7T_cluster) {
    int bond_counter = 0;

    for (int ring_pointer = 0; ring_pointer < 3; ++ring_pointer) {
        if (new_5A_cluster[ring_pointer] == old_7T_cluster[0]) {
            if (new_5A_cluster[ring_pointer] == old_7T_cluster[1]) {
                if (new_5A_cluster[ring_pointer] == old_7T_cluster[6]) {
                   return 1; 
                }
            }
        
        }
        if (new_5A_cluster[ring_pointer] == old_7T_cluster[2]) {
            if (new_5A_cluster[ring_pointer] == old_7T_cluster[3]) {
                if (new_5A_cluster[ring_pointer] == old_7T_cluster[6]) {
                    return 1;
                }
            }  
        }
    }
    return 0;
}

int get_new_particle_8PAA(const int *new_5A_cluster, int spindle_id) {
    if (new_5A_cluster[3] == spindle_id) {
        //printf("5A %i\n", new_5A_cluster[3]);
        return new_5A_cluster[4];
    }
    else if (new_5A_cluster[4] == spindle_id) {
        return new_5A_cluster[3];
    }
    else {
        Error("New spindle not found.");
        return 0;
    }
}



int check_unique8PAA(const int *old_clust, int new_particle, int new_clust_size){
    int u;
    for (int old_8PAA_id = 0; old_8PAA_id < n8PAA; ++old_8PAA_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<7; ++r){
                if(hc8PAA[old_8PAA_id][q] == old_clust[r]){
                    u += 1;
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc8PAA[old_8PAA_id][p] == new_particle){
                    u += 1;
                }
        }
    if(u == new_clust_size){
        return 1;            
    }
    }
    return 0;
}
