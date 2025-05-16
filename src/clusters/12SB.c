#include <globals.h>
#include <tools.h>
#include "12SB.h"
#include "11S.h"
//!  A 12SB is an 11SB and a 5A 
// the [8] particle in the 11SB is one spindle of the 5A
// the ring particles of the 5A are all in the [0:7] of the 11SB

void Clusters_Get12SB() {
    for (int first_11SB_id = 0; first_11SB_id < n11SB; ++first_11SB_id) {
        int *new_12SB = malloc(12 * sizeof(int));
        int *first_11SB_cluster = hc11SB[first_11SB_id];
        for (int first_5A_id = 0; first_5A_id < nsp3c; ++first_5A_id) {
            int *first_5A_cluster = hcsp3c[first_5A_id];
            if(first_11SB_cluster[8] == first_5A_cluster[3] || first_11SB_cluster[8] == first_5A_cluster[4]){
                if(first_11SB_cluster[8] != first_5A_cluster[3] || first_11SB_cluster[8] != first_5A_cluster[4]){
                    if(common_ring(first_5A_cluster, first_11SB_cluster, 3,11) == 3){
                        int final_part;
                        int count = 0;
                        int count2 = 0;
                        if(first_11SB_cluster[8] == first_5A_cluster[3]){
                            final_part = first_5A_cluster[4];
                        }
                        if(first_11SB_cluster[8] == first_5A_cluster[4]){
                            final_part = first_5A_cluster[3];
                        }
                        for(int i = 0; i < 11; ++i){
                            if(first_11SB_cluster[i] == final_part){
                                count += 1;
                            }
                        }
                        for(int j = 0; j < 3; ++j){
                            if(first_5A_cluster[j] == first_11SB_cluster[0] ||first_5A_cluster[j] == first_11SB_cluster[0]||
                            first_5A_cluster[j] == first_11SB_cluster[2]|| first_5A_cluster[j] == first_11SB_cluster[3]){
                                count2 += 1;
                            }
                        }
                        if(count == 0 && count2 == 1){
                            if(check_unique_12SB(first_11SB_cluster, final_part) == 0){
                                //printf("5A %i %i %i %i %i\n",first_5A_cluster[0],first_5A_cluster[1],first_5A_cluster[2],
                                //first_5A_cluster[3],first_5A_cluster[4]);
                                //printf("11SB %i %i %i %i %i %i %i %i %i %i %i\n",first_11SB_cluster[0],first_11SB_cluster[1],first_11SB_cluster[2],
                                //first_11SB_cluster[3],first_11SB_cluster[4],first_11SB_cluster[5],first_11SB_cluster[6],first_11SB_cluster[7],
                                //first_11SB_cluster[8],first_11SB_cluster[9],first_11SB_cluster[10]);
                                //printf("final part %i\n", final_part);
                                add_12SB(first_11SB_cluster, final_part);
                            }
                        }
                    }
                }
            }
        }
    }
}

void add_12SB(const int *old_11SB, int new_part) {
    int clusSize = 12;
    if (n12SB == m12SB) {
        hc12SB = resize_2D_int(hc12SB, m12SB, m12SB + incrStatic, clusSize, -1);
        m12SB = m12SB + incrStatic;
    }
    hc12SB[n12SB][0] = old_11SB[0];
    hc12SB[n12SB][1] = old_11SB[1];
    hc12SB[n12SB][2] = old_11SB[2];
    hc12SB[n12SB][3] = old_11SB[3];
    hc12SB[n12SB][4] = old_11SB[4];
    hc12SB[n12SB][5] = old_11SB[5];
    hc12SB[n12SB][6] = old_11SB[6];
    hc12SB[n12SB][7] = old_11SB[7];
    hc12SB[n12SB][8] = old_11SB[8];
    hc12SB[n12SB][9] = old_11SB[9];
    hc12SB[n12SB][10] = old_11SB[10];
    hc12SB[n12SB][11] = new_part;
    for (int i = 0; i < 12 ; ++i) {
        s12SB[hc12SB[n12SB][i]] = 'B';
    }
    ++n12SB;
}

int check_unique_12SB(const int *old_11SB, int new_part){
    int u;
    for (int old_12SB_id = 0; old_12SB_id < n12SB; ++old_12SB_id) {
        u = 0;
        for (int q = 0; q < 11; ++q){
            for (int r = 0; r < 12; ++r){
                if(hc12SB[old_12SB_id][r] == old_11SB[q]){
                    u += 1;
                }
            }

        }
    for (int l = 0; l < 12; ++l){
        if(hc12SB[old_12SB_id][l] == new_part){
            u += 1;
        }
    }       
    if(u == 12){
        return 1;            
    }
    }
    return 0;
}