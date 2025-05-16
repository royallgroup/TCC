#include <globals.h>
#include <tools.h>
#include "10PAB.h"
#include "7T.h"
#include "7PAB.h"
#include "8PAB.h"
#include "9PAA.h"
#include "12PAB.h"
#include "13PBB.h"
//!  A 13PBB cluster is made of a 12PAB cluster with one extra particle which forms a 3 ring
/*!
*  13PBB is made up of a 12PAB cluster and a sp3a cluster
*  The sp3a cluster shares two particles with the 13PBB cluster
*  One particle of the sp3a is a bonded spindle of the 12PAB, and one is a non- bonded spindle
*
*  Cluster output: BBBBBBBBBBBBB
*  Storage order: original 12PAB particles x 12, new sp3a particle)
*/

    void Clusters_Get13PBB() {
    for (int old_12PAB_id = 0; old_12PAB_id < n12PAB; ++old_12PAB_id) {
        int *old_12PAB_cluster = hc12PAB[old_12PAB_id];
        int clust_12PAB[12]; 
        
        clust_12PAB[0] = old_12PAB_cluster[9];
        clust_12PAB[1] = old_12PAB_cluster[1];
        clust_12PAB[2] = old_12PAB_cluster[6];
        clust_12PAB[3] = old_12PAB_cluster[10];
        clust_12PAB[4] = old_12PAB_cluster[8];
        clust_12PAB[5] = old_12PAB_cluster[7];
        clust_12PAB[6] = old_12PAB_cluster[2];
        clust_12PAB[7] = old_12PAB_cluster[5];
        clust_12PAB[8] = old_12PAB_cluster[4];
        clust_12PAB[9] = old_12PAB_cluster[0];
        clust_12PAB[10] = old_12PAB_cluster[3];
        clust_12PAB[11] = old_12PAB_cluster[11];
        int spindle_id_2 = clust_12PAB[5];
        int spindle_id = clust_12PAB[4];
        int nbonded_id = clust_12PAB[3];
        int bonded_id = clust_12PAB[0]; 
        int nbonded_id_2 = clust_12PAB[1];
        int bonded_id_2 = clust_12PAB[2]; 
        int final_part1 = clust_12PAB[6]; 
        int final_part2 = clust_12PAB[7]; 
        int final_part3 = clust_12PAB[8]; 
        int final_part4 = clust_12PAB[9]; 
        int final_part5 = clust_12PAB[10]; 
        int final_part6 = clust_12PAB[11]; 

        int arr_6Z[6] = {bonded_id, nbonded_id_2, bonded_id_2, nbonded_id, spindle_id, spindle_id_2};
        for (int new_3_id = 0; new_3_id < nsp3a; ++new_3_id ){ //iterate over all sp3a rings
            int *new_3_cluster = hcsp3a[new_3_id];
            int rinc = nring_in_cluster(old_12PAB_cluster, new_3_cluster, 12);
            if(nring_in_cluster(old_12PAB_cluster, new_3_cluster, 12) == 2){
                if(nring_bonded_nbonded(bonded_id, nbonded_id, new_3_cluster) == 2){ //two particles in ring are in 6Z
                   int bond_counter = check_ring_bonds(new_3_cluster, arr_6Z);
                    if (bond_counter == 9){
                        int new_part = get_new_particle_P2(new_3_cluster, bonded_id, nbonded_id);
                        if (check_unique13PBB(old_12PAB_cluster, new_part, 13) == 0){
                            if(check_5A(arr_6Z, new_part) ==0){
                                add_13PBB(bonded_id, nbonded_id_2, bonded_id_2, nbonded_id,spindle_id, spindle_id_2, final_part1, final_part2, final_part3, final_part4, final_part5, final_part6, new_part);
                                }
                            }   
                        }
                    }
                }
            }
    } 
} 


void add_13PBB(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int final_part2, int final_part3, int final_part4, int final_part5, int final_part6, int new_particle) {
    int clusSize = 13;
    if (n13PBB == m13PBB) {
        hc13PBB = resize_2D_int(hc13PBB, m13PBB, m13PBB + incrStatic, clusSize, -1);
        m13PBB = m13PBB + incrStatic;
    }
    hc13PBB[n13PBB][0] = bonded_id;
    hc13PBB[n13PBB][1] = nbonded_id_2;
    hc13PBB[n13PBB][2] = bonded_id_2;
    hc13PBB[n13PBB][3] = nbonded_id;
    hc13PBB[n13PBB][4] = spindle_id;
    hc13PBB[n13PBB][5] = spindle_id_2;
    hc13PBB[n13PBB][6] = final_part1;
    hc13PBB[n13PBB][7] = final_part2;
    hc13PBB[n13PBB][8] = final_part3;
    hc13PBB[n13PBB][9] = final_part4;
    hc13PBB[n13PBB][10] = final_part5;
    hc13PBB[n13PBB][11] = final_part6;
    hc13PBB[n13PBB][12] = new_particle;

    for (int i = 0; i < 12; ++i) {
       s13PBB[hc13PBB[n13PBB][i]] = 'B';
    }
    ++n13PBB;
}

int check_unique13PBB(const int *old_clust, int new_particle, int new_clust_size){
    int u;
    for (int old_13PBB_id = 0; old_13PBB_id < n13PBB; ++old_13PBB_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<12; ++r){
                if(hc13PBB[old_13PBB_id][q] == old_clust[r]){
                    u += 1;
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc13PBB[old_13PBB_id][p] == new_particle){
                    u += 1;
                }
        }
    if(u == new_clust_size){
        return 1;            
    }
    }

    for (int old_13PAB_id = 0; old_13PAB_id < n13PAB; ++old_13PAB_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<12; ++r){
                if(hc13PAB[old_13PAB_id][q] == old_clust[r]){
                    u += 1;
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc13PAB[old_13PAB_id][p] == new_particle){
                    u += 1;
                }
        }
    if(u == new_clust_size){
        return 1;            
    }
    }
    return 0;
}