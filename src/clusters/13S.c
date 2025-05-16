#include <globals.h>
#include <tools.h>
#include "13S.h"
#include "12SB.h"
//!  A 13S is an 12SB and a 5A 
// the [8] particle in the 12SB is one spindle of the 5A
// the ring particles of the 5A are all in the [0:7] of the 12SB

void Clusters_Get13S() {
    int final_part;
    for (int first_12SB_id = 0; first_12SB_id < n12SB; ++first_12SB_id) {
        int *first_12SB_cluster = hc12SB[first_12SB_id];
        for (int first_5A_id = 0; first_5A_id < nsp3c; ++first_5A_id) {
            int *first_5A_cluster = hcsp3c[first_5A_id];
            int count1 = 0;
            int count2 = 0;
            int count3 = 0;
            int count4 = 0;
            for(int i = 0; i < 3; ++i){
                if(first_5A_cluster[i] == first_12SB_cluster[0] || first_5A_cluster[i] == first_12SB_cluster[1] ||
                first_5A_cluster[i] == first_12SB_cluster[2] || first_5A_cluster[i] == first_12SB_cluster[3]){
                    count1 += 1;
                }
                if(first_5A_cluster[i] == first_12SB_cluster[9] || first_5A_cluster[i] == first_12SB_cluster[10]){
                    count2 += 1;
                }
                if(first_5A_cluster[i] == first_12SB_cluster[11]){
                    count3 += 1;
                }
            }
            if(count1 == 1 && count2 == 1 && count3 == 1){
                if(first_5A_cluster[3] == first_12SB_cluster[6] || first_5A_cluster[3] == first_12SB_cluster[7] ||
                first_5A_cluster[4] == first_12SB_cluster[6] || first_5A_cluster[4] == first_12SB_cluster[7]){
                    if(first_5A_cluster[3] == first_12SB_cluster[6] || first_5A_cluster[3] == first_12SB_cluster[7]){
                        final_part = first_5A_cluster[4];
                    }
                    if(first_5A_cluster[4] == first_12SB_cluster[6] || first_5A_cluster[4] == first_12SB_cluster[7]){
                        final_part = first_5A_cluster[3];
                    }
                    for(int j = 0; j < 12; ++j){
                        for(int k = 0; k < 5; ++k){
                            if(first_12SB_cluster[j] == first_5A_cluster[k]){
                                count4 += 1;
                            }
                        }
                    }
                    if(count4 == 4){
                        if(check_unique_13S(first_12SB_cluster, final_part) == 0){
                            add_13S(first_12SB_cluster, final_part);
                        }
                    }
                }
            }
        }
    }
}

void add_13S(const int *old_12SB, int new_part) {
    int clusSize = 13;
    if (n13S == m13S) {
        hc13S = resize_2D_int(hc13S, m13S, m13S + incrStatic, clusSize, -1);
        m13S = m13S + incrStatic;
    }
    hc13S[n13S][0] = old_12SB[0];
    hc13S[n13S][1] = old_12SB[1];
    hc13S[n13S][2] = old_12SB[2];
    hc13S[n13S][3] = old_12SB[3];
    hc13S[n13S][4] = old_12SB[4];
    hc13S[n13S][5] = old_12SB[5];
    hc13S[n13S][6] = old_12SB[6];
    hc13S[n13S][7] = old_12SB[7];
    hc13S[n13S][8] = old_12SB[8];
    hc13S[n13S][9] = old_12SB[9];
    hc13S[n13S][10] = old_12SB[10];
    hc13S[n13S][11] = old_12SB[11];
    hc13S[n13S][12] = new_part;
    for (int i = 0; i < 13; ++i) {
        s13S[hc13S[n13S][i]] = 'B';
    }
    ++n13S;
}

int check_unique_13S(const int *old_12SB, int new_part){
    int u;
    for (int old_13S_id = 0; old_13S_id < n13S; ++old_13S_id) {
        u = 0;
        for (int q = 0; q < 12; ++q){
            for (int r = 0; r < 13; ++r){
                if(hc13S[old_13S_id][r] == old_12SB[q]){
                    u += 1;
                }
            }

        }
    for (int l = 0; l < 13; ++l){
        if(hc13S[old_13S_id][l] == new_part){
            u += 1;
        }
    }       
    if(u == 13){
        return 1;            
    }
    }
    return 0;
}