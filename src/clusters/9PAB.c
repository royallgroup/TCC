#include <globals.h>
#include <tools.h>
#include "9PAB.h"
#include "8PAA.h"
#include "7PAB.h"
#include "7T.h"
//!  A 9PAB is an 8PAA with an additional particle
/*!
* A 8PAA cluster contains two 6Z clusters at the top and bottom of the cluster. 
* A 9PAB takes these two 6Z clusters and tries to make a 7T from the cluster and an extra particle.
*  Cluster output: BBBBBBBBB
*  Storage order: original 8PAA particles x 8, new 5A spindle)
*/

    void Clusters_Get9PAB() {
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
            int arr_6Z[6] = {bonded_id, nbonded_id_2, bonded_id_2, nbonded_id, spindle_id, spindle_id_2};
            for (int new_3_id = 0; new_3_id < nsp3a; ++new_3_id ){ //iterate over all sp3a rings
                int *new_3_cluster = hcsp3a[new_3_id];
                if(nring_in_cluster(old_8PAA_cluster, new_3_cluster, 8) == 2){
                    if(nring_bonded_nbonded(bonded_id, nbonded_id, new_3_cluster) == 2){ //two particles in ring are in 6Z
                    int bond_counter = check_ring_bonds(new_3_cluster, arr_6Z);
                        if (bond_counter == 9){
                            int new_part = get_new_particle_P2(new_3_cluster, bonded_id, nbonded_id);
                            if (check_unique9PAB(old_8PAA_cluster, new_part, 9) == 0){
                                if(check_5A(arr_6Z, new_part) == 0){
                                add_9PAB(bonded_id, nbonded_id_2, bonded_id_2, nbonded_id,spindle_id, spindle_id_2, final_part1, final_part2, new_part);
                                }
                                }
                            }
                        }
                    }
                }
        }
    } 
}


void add_9PAB(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int final_part2, int new_particle) {
    int clusSize = 9;
    //printf("new_particle %i\n", new_particle);
    if (n9PAB == m9PAB) {
        hc9PAB = resize_2D_int(hc9PAB, m9PAB, m9PAB + incrStatic, clusSize, -1);
        m9PAB = m9PAB + incrStatic;
    }
    hc9PAB[n9PAB][0] = bonded_id;
    hc9PAB[n9PAB][1] = nbonded_id_2;
    hc9PAB[n9PAB][2] = bonded_id_2;
    hc9PAB[n9PAB][3] = nbonded_id;
    hc9PAB[n9PAB][4] = spindle_id;
    hc9PAB[n9PAB][5] = spindle_id_2;
    hc9PAB[n9PAB][6] = final_part1;
    hc9PAB[n9PAB][7] = final_part2;
    hc9PAB[n9PAB][8] = new_particle;

    for (int i = 0; i < 9; ++i) {
        s9PAB[hc9PAB[n9PAB][i]] = 'B';
    }
    ++n9PAB;
}

int check_unique9PAB(const int *old_clust, int new_particle, int new_clust_size){
    int u;
    for (int old_9PAB_id = 0; old_9PAB_id < n9PAB; ++old_9PAB_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<8; ++r){
                //printf("%i %i\n", hc11P[old_11P_id][q],old_clust[r]);
                if(hc9PAB[old_9PAB_id][q] == old_clust[r]){
                    u += 1;
             
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc9PAB[old_9PAB_id][p] == new_particle){
                    u += 1;
                }
        }

    if(u == new_clust_size){
        return 1;            
    }
    }

    for (int old_9PAA_id = 0; old_9PAA_id < n9PAA; ++old_9PAA_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<8; ++r){
                //printf("%i %i\n", hc11P[old_11P_id][q],old_clust[r]);
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
