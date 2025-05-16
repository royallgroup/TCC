#include <globals.h>
#include <tools.h>
#include "7T.h"
#include "7PAB.h"
#include "9PAA.h"
#include "8PAB.h"
#include "12PAB.h"
#include "12PAB.h"
//!  A 12PAB is a 11PAA with an additional particle
/*!
* A 11PAA cluster contains two 6Z clusters at the top and bottom of the cluster. 
* A 12PAB takes these two 6Z clusters and tries to make a 7T from the cluster and an extra particle.
*  Cluster output: BBBBBBBBBB
*  Storage order: original 11PAA particles x 9, new 5A spindle)
*/

    void Clusters_Get12PAB() {
    for (int old_11PAA_id = 0; old_11PAA_id < n11PAA; ++old_11PAA_id) {
        int *old_11PAA_cluster = hc11PAA[old_11PAA_id];
        int clust_11PAA[11][2]; 
        clust_11PAA[0][0] = old_11PAA_cluster[3];
        clust_11PAA[1][0] = old_11PAA_cluster[2];
        clust_11PAA[2][0] = old_11PAA_cluster[5];
        clust_11PAA[3][0] = old_11PAA_cluster[10];
        clust_11PAA[4][0] = old_11PAA_cluster[0];
        clust_11PAA[5][0] = old_11PAA_cluster[4];
        clust_11PAA[6][0] = old_11PAA_cluster[1];
        clust_11PAA[7][0] = old_11PAA_cluster[6];
        clust_11PAA[8][0] = old_11PAA_cluster[7];
        clust_11PAA[9][0] = old_11PAA_cluster[8];
        clust_11PAA[10][0] = old_11PAA_cluster[9];
        
        clust_11PAA[0][1] = old_11PAA_cluster[8];
        clust_11PAA[1][1] = old_11PAA_cluster[2];
        clust_11PAA[2][1] = old_11PAA_cluster[1];
        clust_11PAA[3][1] = old_11PAA_cluster[9];
        clust_11PAA[4][1] = old_11PAA_cluster[7];
        clust_11PAA[5][1] = old_11PAA_cluster[6];
        clust_11PAA[6][1] = old_11PAA_cluster[5];
        clust_11PAA[7][1] = old_11PAA_cluster[4];
        clust_11PAA[8][1] = old_11PAA_cluster[0];
        clust_11PAA[9][1] = old_11PAA_cluster[3];
        clust_11PAA[10][1] = old_11PAA_cluster[10];

        for(int g = 0; g<2; ++g){
            int spindle_id_2 = clust_11PAA[5][g];
            int spindle_id = clust_11PAA[4][g];
            int nbonded_id = clust_11PAA[3][g];
            int bonded_id = clust_11PAA[0][g]; 
            int nbonded_id_2 = clust_11PAA[1][g];
            int bonded_id_2 = clust_11PAA[2][g]; 
            int final_part1 = clust_11PAA[6][g]; 
            int final_part2 = clust_11PAA[7][g]; 
            int final_part3 = clust_11PAA[8][g]; 
            int final_part4 = clust_11PAA[9][g];
            int final_part5 = clust_11PAA[10][g];
            int arr_6Z[6] = {bonded_id, nbonded_id_2, bonded_id_2, nbonded_id, spindle_id, spindle_id_2};
            for (int new_3_id = 0; new_3_id < nsp3a; ++new_3_id ){ //iterate over all sp3a rings
                int *new_3_cluster = hcsp3a[new_3_id];
                if(nring_in_cluster(old_11PAA_cluster, new_3_cluster, 11) == 2){
                    if(nring_bonded_nbonded(bonded_id, nbonded_id, new_3_cluster) == 2){ //two particles in ring are in 6Z
                    int bond_counter = check_ring_bonds(new_3_cluster, arr_6Z);
                        if (bond_counter == 9){
                            int new_part = get_new_particle_P2(new_3_cluster, bonded_id, nbonded_id);
                            if (check_unique12PAB(old_11PAA_cluster, new_part, 12) == 0){
                                if(check_5A(arr_6Z, new_part) == 0){
                                    add_12PAB(bonded_id, nbonded_id_2, bonded_id_2, nbonded_id,spindle_id, spindle_id_2, final_part1, final_part2, final_part3, final_part4, final_part5, new_part);
                                }
                                }
                            }
                        }
                    }
                }
            }
        } 
    } 

void add_12PAB(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int final_part2, int final_part3, int final_part4, int final_part5, int new_particle){
    int clusSize = 12;
    if (n12PAB == m12PAB) {
        hc12PAB = resize_2D_int(hc12PAB, m12PAB, m12PAB + incrStatic, clusSize, -1);
        m12PAB = m12PAB + incrStatic;
    }
    hc12PAB[n12PAB][0] = bonded_id;
    hc12PAB[n12PAB][1] = nbonded_id_2;
    hc12PAB[n12PAB][2] = bonded_id_2;
    hc12PAB[n12PAB][3] = nbonded_id;
    hc12PAB[n12PAB][4] = spindle_id;
    hc12PAB[n12PAB][5] = spindle_id_2;
    hc12PAB[n12PAB][6] = final_part1;
    hc12PAB[n12PAB][7] = final_part2;
    hc12PAB[n12PAB][8] = final_part3;
    hc12PAB[n12PAB][9] = final_part4;
    hc12PAB[n12PAB][10] = final_part5;
    hc12PAB[n12PAB][11] = new_particle;

    for (int i = 0; i < 12; ++i) {
        s12PAB[hc12PAB[n12PAB][i]] = 'B';
    }
    ++n12PAB;
}

int check_unique12PAB(const int *old_clust, int new_particle, int new_clust_size){
    int u;
    for (int old_12PAB_id = 0; old_12PAB_id < n12PAB; ++old_12PAB_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<11; ++r){
                if(hc12PAB[old_12PAB_id][q] == old_clust[r]){
                    u += 1;
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc12PAB[old_12PAB_id][p] == new_particle){
                    u += 1;
                }
        }
    if(u == new_clust_size){
        return 1;            
    }
    }

    for (int old_12PAA_id = 0; old_12PAA_id < n12PAA; ++old_12PAA_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<11; ++r){
                if(hc12PAA[old_12PAA_id][q] == old_clust[r]){
                    u += 1;
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc12PAA[old_12PAA_id][p] == new_particle){
                    u += 1;
                }
        }
    if(u == new_clust_size){
        return 1;            
    }
    }
    return 0;
}