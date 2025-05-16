#include <globals.h>
#include <tools.h>
#include "9PAA.h"
#include "8PAA.h"
#include "7T.h"
//!  A 9PAA is an 8PAA with an additional particle
/*!
* A 8PAA cluster contains two 6Z clusters at the top and bottom of the cluster. 
* A 9PAA takes these two 6Z clusters and tries to make a 7T from the cluster and an extra particle.
*  Cluster output: BBBBBBBBB
*  Storage order: original 8PAA particles x 8, new 5A spindle)
*/

    void Clusters_Get9PAA() {
    for (int old_8PAA_id = 0; old_8PAA_id < n8PAA; ++old_8PAA_id) {
        int *old_8PAA_cluster = hc8PAA[old_8PAA_id];
        int clust_6Z[8][2]; 
        clust_6Z[0][0] = old_8PAA_cluster[3];
        clust_6Z[1][0] = old_8PAA_cluster[2];
        clust_6Z[2][0] = old_8PAA_cluster[5];
        clust_6Z[3][0] = old_8PAA_cluster[7];
        clust_6Z[4][0] = old_8PAA_cluster[0];
        clust_6Z[5][0] = old_8PAA_cluster[4];
        clust_6Z[6][0] = old_8PAA_cluster[1];
        clust_6Z[7][0] = old_8PAA_cluster[6];
        
        clust_6Z[0][1] = old_8PAA_cluster[1];
        clust_6Z[1][1] = old_8PAA_cluster[0];
        clust_6Z[2][1] = old_8PAA_cluster[4];
        clust_6Z[3][1] = old_8PAA_cluster[6];
        clust_6Z[4][1] = old_8PAA_cluster[2];
        clust_6Z[5][1] = old_8PAA_cluster[5];
        clust_6Z[6][1] = old_8PAA_cluster[3];
        clust_6Z[7][1] = old_8PAA_cluster[7];

        for(int g = 0; g<2; ++g){
            int spindle_id_2 = clust_6Z[5][g];
            int spindle_id = clust_6Z[4][g];
            int nbonded_id = clust_6Z[3][g];
            int bonded_id = clust_6Z[0][g]; 
            int nbonded_id_2 = clust_6Z[1][g];
            int bonded_id_2 = clust_6Z[2][g]; 
            int final_part1 = clust_6Z[6][g]; 
            int final_part2 = clust_6Z[7][g]; 
            for (int new_5A_pointer = 0; new_5A_pointer < nmem_sp3c[spindle_id_2]; ++new_5A_pointer) {
                //int new_5A_id = mem_sp3c[spindle_id_2][new_5A_pointer];
                int *new_5A_cluster = hcsp3c[new_5A_pointer];
                if (is_particle_spindle_of_5A(spindle_id_2, new_5A_cluster) == 1) {
                    if(is_particle_ring_of_5A(bonded_id,new_5A_cluster) == 1){
                        if(is_particle_ring_of_5A(nbonded_id,new_5A_cluster) == 1){
                            if(is_particle_ring_of_5A(spindle_id,new_5A_cluster) == 1){
                                int new_particle_id = get_new_particle(new_5A_cluster, spindle_id_2);
                                    if(crossover_part(old_8PAA_cluster, new_particle_id, 9) == 0){
                                    int arr_6Z[6] = {bonded_id, nbonded_id_2, bonded_id_2, nbonded_id, spindle_id, spindle_id_2};
                                    int bond_counter = check_ring_bonds(new_5A_cluster, arr_6Z);
                                        if (bond_counter == 25 ){
                                            if (check_unique9PAA(old_8PAA_cluster, new_particle_id,9) == 0){
                                            add_9PAA(bonded_id, nbonded_id_2, bonded_id_2, nbonded_id,spindle_id, spindle_id_2, final_part1, final_part2, new_particle_id);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    } 
}


void add_9PAA(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int final_part2, int new_particle) {
    int clusSize = 9;
    if (n9PAA == m9PAA) {
        hc9PAA = resize_2D_int(hc9PAA, m9PAA, m9PAA + incrStatic, clusSize, -1);
        m9PAA = m9PAA + incrStatic;
    }
    hc9PAA[n9PAA][0] = bonded_id;
    hc9PAA[n9PAA][1] = nbonded_id_2;
    hc9PAA[n9PAA][2] = bonded_id_2;
    hc9PAA[n9PAA][3] = nbonded_id;
    hc9PAA[n9PAA][4] = spindle_id;
    hc9PAA[n9PAA][5] = spindle_id_2;
    hc9PAA[n9PAA][6] = final_part1;
    hc9PAA[n9PAA][7] = final_part2;
    hc9PAA[n9PAA][8] = new_particle;

    for (int i = 0; i < 9; ++i) {
        s9PAA[hc9PAA[n9PAA][i]] = 'B';
    }
    ++n9PAA;
}

int check_unique9PAA(const int *old_clust, int new_particle, int new_clust_size){
    int u;
    for (int old_9PAA_id = 0; old_9PAA_id < n9PAA; ++old_9PAA_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<8; ++r){
          
                if(hc9PAA[old_9PAA_id][q] == old_clust[r]){
                    u += 1;
             
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc9PAA[old_9PAA_id][p] == new_particle){
                    u += 1;
                }
        }

    if(u == new_clust_size){
        return 1;            
    }
    }
    return 0;
}
