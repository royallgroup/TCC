#include <globals.h>
#include <tools.h>
#include "11S.h"
#include "8PAB.h"
#include "7PAB.h"
#include "9S.h"
#include "7T.h"
//!  A 11S is a made of a 9B cluster intersecting with a 7A cluster
/*!
* 5 particles in the 7A are shared with the 9B
* The other spindle of each outer 7A is a ring if the central 7A: s3, s4
* There is one common ring cr1 shared between the three 7As
* each outer 7A shares one ring with the central one: cr2, cr3
* the outer 7As have two ring particles not shared by the other 7As: r1, r2, r3, r4
* storage order: s1, s2, s3, s4, cr1, cr2, cr3, r1, r2, r3, r4
*/

void Clusters_Get11S() {
    int new_part;
    for (int first_9B_id = 0; first_9B_id < n9B; ++first_9B_id) {

        int *first_9B_cluster = hc9B[first_9B_id];
            for (int first_7A_id = 0; first_7A_id < nsp5c; ++first_7A_id) {
                int *first_7A_cluster = hcsp5c[first_7A_id];
                if(overlap_7A_9B(first_7A_cluster,first_9B_cluster) == 5){
                if(i_in_ring(first_7A_cluster,first_9B_cluster[4]) == 1 || i_in_ring(first_7A_cluster,first_9B_cluster[5]) == 1){
                    if(i_in_ring(first_7A_cluster,first_9B_cluster[8]) == 1){
                        if(i_in_ring(first_7A_cluster,first_9B_cluster[0]) == 1 || i_in_ring(first_7A_cluster,first_9B_cluster[1]) == 1 ||
                        i_in_ring(first_7A_cluster,first_9B_cluster[2]) == 1 || i_in_ring(first_7A_cluster,first_9B_cluster[3]) == 1){
                            if(i_in_spindle(first_7A_cluster,first_9B_cluster[0]) == 1 || i_in_spindle(first_7A_cluster,first_9B_cluster[1]) == 1 ||
                            i_in_spindle(first_7A_cluster,first_9B_cluster[2]) == 1 || i_in_spindle(first_7A_cluster,first_9B_cluster[3]) == 1){
                                if(i_in_spindle(first_7A_cluster,first_9B_cluster[6]) == 1 || i_in_spindle(first_7A_cluster,first_9B_cluster[7]) == 1){
                                    int *new_11S = malloc(11 * sizeof(int));
                                    get_11S(first_7A_cluster, first_9B_cluster, &new_11S);
                                    if(check_unique_11S(&new_11S) == 0){
                                        add_11S(&new_11S);
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

int check_unique_11S(int **clust){
    int u;
    for (int old_11S_id = 0; old_11S_id < n11S; ++old_11S_id) {
        u = 0;
        for (int q = 0; q < 11; ++q){
            for (int r = 0; r < 11; ++r){
                if(hc11S[old_11S_id][r] == (*clust)[q]){
                    u += 1;
                }
            }

        }
    if(u == 11){
        return 1;            
    }
    }
    return 0;
}

void get_11S(const int *clust_7A, const int *clust_9B, int **clust_11S){
    (*clust_11S)[0] = clust_9B[0];
    (*clust_11S)[1] = clust_9B[1];
    (*clust_11S)[2] = clust_9B[2];
    (*clust_11S)[3] = clust_9B[3];
    (*clust_11S)[4] = clust_9B[4];
    (*clust_11S)[5] = clust_9B[5];
    (*clust_11S)[6] = clust_9B[6];
    (*clust_11S)[7] = clust_9B[7];
    (*clust_11S)[8] = clust_9B[8];
    int u = 9;
    int v;
    for(int i = 0; i < 7; ++i){
        v = 0;
        for(int j = 0; j < 9; ++j){
            if(clust_9B[j] == clust_7A[i]){
                v += 1;
            }
        }
        if(v == 0){
            (*clust_11S)[u] = clust_7A[i];
            u += 1;
        }
    }
    return;
    }

int overlap_7A_9B(const int *clust_7A,const int *clust_9B){
    int u = 0;
    for(int i = 0; i < 9; ++i){
        for(int j = 0; j < 7; ++j){
            if(clust_7A[j] == clust_9B[i]){
                u += 1;
            }
        }
    }
    return u;
}

int i_in_ring(const int *clust_7A, int i){
    for(int j = 0; j < 5; ++j){
        if (clust_7A[j] == i){
        return 1;
        }
    }
    return 0;
}

int i_in_spindle(const int *clust_7A, int i){
    for(int j = 5; j < 7; ++j){
        if (clust_7A[j] == i){
        return 1;
        }
    }
    return 0;
}


void add_11S(int **clust) {
    int clusSize = 11;
    if (n11S == m11S) {
        hc11S = resize_2D_int(hc11S, m11S, m11S + incrStatic, clusSize, -1);
        m11S = m11S + incrStatic;
    }
    hc11S[n11S][0] = (*clust)[0];
    hc11S[n11S][1] = (*clust)[1];
    hc11S[n11S][2] = (*clust)[2];
    hc11S[n11S][3] = (*clust)[3];
    hc11S[n11S][4] = (*clust)[4];
    hc11S[n11S][5] = (*clust)[5];
    hc11S[n11S][6] = (*clust)[6];
    hc11S[n11S][7] = (*clust)[7];
    hc11S[n11S][8] = (*clust)[8];
    hc11S[n11S][9] = (*clust)[9];
    hc11S[n11S][10] = (*clust)[10];
    for (int i = 0; i < 11 ; ++i) {
        s11S[hc11S[n11S][i]] = 'B';
    }
    ++n11S;
}

int common_spindle(const int *clust1, const int *clust2, int r1, int r2){
    int count = 0;
    for(int i = r1; i < r1 + 2; ++i){
        for(int j = r2; j < r2 + 2; ++j){
            if(clust1[i] == clust2[j]){
                count += 1;
            }
        }
    }
    return count;
}

int common_ring(const int *clust1, const int *clust2, int r1, int r2){
    int count = 0;
    for(int i = 0; i < r1; ++i){
        for(int j = 0; j < r2; ++j){
            if(clust1[i] == clust2[j]){
                count += 1;
            }
        }
    }
    return count;
}
