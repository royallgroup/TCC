#include <globals.h>
#include <tools.h>
#include "11SB.h"
#include "11SB.h"

//!  A 11SB is a made of a 9B cluster and two 5A clusters
/*!
* One spindle of each 9B is a ring of the other 9B
* There is one common ring particle
*/

void Clusters_Get11SB() {
    for (int old_9B_id = 0; old_9B_id < n9B; ++old_9B_id) {
        int spindle1, spindle2;
        int *old_9B_cluster = hc9B[old_9B_id];
        for (int first_5A_id = 0; first_5A_id < nsp3c; ++first_5A_id) {
            int *first_5A_cluster = hcsp3c[first_5A_id];
            if (first_5A_cluster[3] == old_9B_cluster[8]){
                spindle1 = first_5A_cluster[4];
                }
            else if (first_5A_cluster[4] == old_9B_cluster[8]){
                spindle1 = first_5A_cluster[3];
            }
            if( first_5A_cluster[0] == old_9B_cluster[2] || first_5A_cluster[0] == old_9B_cluster[3] || first_5A_cluster[0] == old_9B_cluster[7]){
                if( first_5A_cluster[1] == old_9B_cluster[2] || first_5A_cluster[1] == old_9B_cluster[3] || first_5A_cluster[1] == old_9B_cluster[7]){
                    if( first_5A_cluster[2] == old_9B_cluster[2] || first_5A_cluster[2] == old_9B_cluster[3] || first_5A_cluster[2] == old_9B_cluster[7]){
                        for (int second_5A_id = 0; second_5A_id < nsp3c; ++second_5A_id) {
                            if (first_5A_id != second_5A_id){
                                int *second_5A_cluster = hcsp3c[second_5A_id];
                                if (second_5A_cluster[3] == old_9B_cluster[8]){
                                    spindle2 = second_5A_cluster[4];
                                }
                                else if (second_5A_cluster[4] == old_9B_cluster[8]){
                                    spindle2 = second_5A_cluster[3];
                                }
                                if(spindle1 != spindle2){
                                    if( second_5A_cluster[0] == old_9B_cluster[0] || second_5A_cluster[0] == old_9B_cluster[1] || second_5A_cluster[0] == old_9B_cluster[6]){
                                        if( second_5A_cluster[1] == old_9B_cluster[0] || second_5A_cluster[1] == old_9B_cluster[1] || second_5A_cluster[1] == old_9B_cluster[6]){
                                            if( second_5A_cluster[2] == old_9B_cluster[0] || second_5A_cluster[2] == old_9B_cluster[1] || second_5A_cluster[2] == old_9B_cluster[6]){
                                                if(check_unique_11SB(old_9B_cluster, spindle1, spindle2)== 0){
                                                    add_11SB(old_9B_cluster, spindle1, spindle2);
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
        }
    }

int check_unique_11SB(int *clust9B, int s1, int s2){
    int u;
    for (int old_11SB_id = 0; old_11SB_id < n11SB; ++old_11SB_id) {
        u = 0;
        for (int q = 0; q < 9; ++q){
            for (int r = 0; r < 11; ++r){
                if(hc11SB[old_11SB_id][r] == clust9B[q]){
                    u += 1;
                }
            }

        }
            for (int r = 0; r < 11; ++r){
                if(hc11SB[old_11SB_id][r] == s1 || hc11SB[old_11SB_id][r] == s2){
                    u += 1;
                }
            }
    if(u == 11){
        return 1;            
    }
    }
    return 0;
}

void add_11SB(int *clust9B, int s1, int s2) {
    int clusSize = 11;
    //printf("new_particle %i\n", new_particle);
    if (n11SB == m11SB) {
        hc11SB = resize_2D_int(hc11SB, m11SB, m11SB + incrStatic, clusSize, -1);
        m11SB = m11SB + incrStatic;
    }
    hc11SB[n11SB][0] = clust9B[0];
    hc11SB[n11SB][1] = clust9B[1];
    hc11SB[n11SB][2] = clust9B[2];
    hc11SB[n11SB][3] = clust9B[3];
    hc11SB[n11SB][4] = clust9B[4];
    hc11SB[n11SB][5] = clust9B[5];
    hc11SB[n11SB][6] = clust9B[6];
    hc11SB[n11SB][7] = clust9B[7];
    hc11SB[n11SB][8] = clust9B[8];
    hc11SB[n11SB][9] = s1;
    hc11SB[n11SB][10] = s2;

    for (int i = 0; i < 11; ++i) {
        s11SB[hc11SB[n11SB][i]] = 'B';
    }
    ++n11SB;
}