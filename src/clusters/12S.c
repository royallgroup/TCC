#include <globals.h>
#include <tools.h>
#include "12S.h"
#include "11S.h"
//!  A 12S is a made of three intersecting 7A clusters 1,2 3
// Cluster 1 and 2 share 3 particles
// Cluster 1 and 3 share 3 particles
// Cluster 2 and 3 share 4 particles
// Cluster 1 has a spindle which is a ring of both 2 and 3
// Cluster 1 has a ring which is a spindle of both 2 and 3
// 1 and 2 share 1 ring
// 1 and 3 share 1 ring
// 2 and 3 share 1 ring



void Clusters_Get12S() {
    for (int first_7A_id = 0; first_7A_id < nsp5c; ++first_7A_id) {
        int *new_12S = malloc(12 * sizeof(int));
        int *first_7A_cluster = hcsp5c[first_7A_id];
        for (int second_7A_id = 0; second_7A_id < nsp5c; ++second_7A_id) {
            int *second_7A_cluster = hcsp5c[second_7A_id];
            if(second_7A_id != first_7A_id){
                if(overlap(first_7A_cluster, second_7A_cluster,7,7) == 5){
                    for (int third_7A_id = 0; third_7A_id < nsp5c; ++third_7A_id) {
                        int *third_7A_cluster = hcsp5c[third_7A_id];
                        if(second_7A_id != third_7A_id && first_7A_id != third_7A_id){
                            if(overlap(first_7A_cluster, third_7A_cluster,7,7) == 3){
                                if(overlap(second_7A_cluster, third_7A_cluster,7,7) == 3){
                                    
                                    if(common_ring_spindle(first_7A_cluster, second_7A_cluster, third_7A_cluster) == 1){
                                        if(common_spindle_ring(first_7A_cluster, second_7A_cluster, third_7A_cluster) == 1){
                                        //if(common_ring(third_7A_cluster, second_7A_cluster, 5,5) == 1){
                                            if(common_ring(first_7A_cluster, second_7A_cluster, 5,5) == 2){
                                                //int *new_12S = malloc(12 * sizeof(int));
                                                get_12S2(&new_12S,first_7A_cluster,second_7A_cluster,third_7A_cluster);
                                                //get_12S(first_7A_cluster, second_7A_cluster, third_7A_cluster, &new_12S);
                                                if(check_unique_12S(&new_12S) == 0){
                                                    //printf("boop\n");
                                                    add_12S(&new_12S);
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

// cluster 1 has one ring particle which is a spindle of both 2 and 3
int common_ring_spindle(const int *clust1, const int *clust2, const int *clust3){
    int count = 0;
    for(int i = 0; i < 5; ++i){
        for(int j = 5; j < 7; ++j){
            for(int k = 5; k < 7; ++k){
            if(clust1[j] == clust3[i] && clust2[k] == clust3[i]){
                count += 1;
            }
        }
    }
    }
    return count;
}

int common_spindle_ring(const int *clust1, const int *clust2, const int *clust3){
    int count = 0;
    for(int i = 5; i < 7; ++i){
        for(int j = 0; j < 5; ++j){
            for(int k = 0; k < 5; ++k){
            if(clust1[j] == clust3[i] && clust2[k] == clust3[i]){
                count += 1;
            }
        }
    }
    }
    return count;
}

int overlap(const int *cluster1, const int *cluster2, int l1, int l2){
    int u = 0;
    for(int i = 0;i < l1; ++i){
        for(int j = 0;j < l2; ++j){
            if(cluster1[i] == cluster2[j]){
                u += 1;
            }
        }        
    }
    return u;
}

void get_12S2(int **clust_12S, const int *cluster1, const int *cluster2, const int *cluster3){
    (*clust_12S)[0] = cluster1[0];
    (*clust_12S)[1] = cluster1[1];
    (*clust_12S)[2] = cluster1[2];
    (*clust_12S)[3] = cluster1[3];
    (*clust_12S)[4] = cluster1[4];
    (*clust_12S)[5] = cluster1[5];
    (*clust_12S)[6] = cluster1[6];
    int idx = 0;
    int u;
    for(int i = 0; i < 7; ++i){
        u = 0;
        for(int j = 0; j < 7; ++j){
            if((*clust_12S)[j] == cluster2[i]){
                u += 1;
            }
        }   
        if(u == 0){
           (*clust_12S)[7+idx] =  cluster2[i];
           idx += 1;
        }    
    }

    for(int k = 0; k < 7; ++k){
        u = 0;
        for(int l = 0; l < 7+idx; ++l){
            if((*clust_12S)[l] == cluster3[k]){
                u += 1;
            }
        }   
        if(u == 0){
           (*clust_12S)[7+idx] =  cluster3[k];
           idx += 1;
        }    
    }
    return;
    }

void get_12S(const int *first_clust_7A, const int *second_clust_7A, const int *third_clust_7A, int **clust_12S){
    (*clust_12S)[0] = third_clust_7A[0];
    (*clust_12S)[1] = third_clust_7A[1];
    (*clust_12S)[2] = third_clust_7A[2];
    (*clust_12S)[3] = third_clust_7A[3];
    (*clust_12S)[4] = third_clust_7A[4];
    (*clust_12S)[5] = third_clust_7A[5];
    (*clust_12S)[6] = third_clust_7A[6];
    int u = 7;
    int v;
    for(int i = 0; i < 7; ++i){
        v = 0;
        for(int j = 0; j < 7; ++j){
            if((*clust_12S)[j] == first_clust_7A[i]){
                v += 1;
            }
        }
        if(v == 0){
            (*clust_12S)[u] = first_clust_7A[i];
            u += 1;
        }
    }
    for(int i = 0; i < 7; ++i){
        v = 0;
        for(int j = 0; j < 9; ++j){
            if((*clust_12S)[j] == second_clust_7A[i]){
                v += 1;
            }
        }
        if(v == 0){
            (*clust_12S)[u] = second_clust_7A[i];
            u += 1;
        }
    }

    return;
    }

int check_unique_12S(int **clust){
    int u;
    for (int old_12S_id = 0; old_12S_id < n12S; ++old_12S_id) {
        u = 0;
        for (int q = 0; q < 12; ++q){
            for (int r = 0; r < 12; ++r){
                if(hc12S[old_12S_id][r] == (*clust)[q]){
                    u += 1;
                }
            }

        }
    if(u >= 12){
        return 1;            
    }
    }
    return 0;
}

void add_12S(int **clust) {
    int clusSize = 12;
    if (n12S == m12S) {
        hc12S = resize_2D_int(hc12S, m12S, m12S + incrStatic, clusSize, -1);
        m12S = m12S + incrStatic;
    }
    hc12S[n12S][0] = (*clust)[0];
    hc12S[n12S][1] = (*clust)[1];
    hc12S[n12S][2] = (*clust)[2];
    hc12S[n12S][3] = (*clust)[3];
    hc12S[n12S][4] = (*clust)[4];
    hc12S[n12S][5] = (*clust)[5];
    hc12S[n12S][6] = (*clust)[6];
    hc12S[n12S][7] = (*clust)[7];
    hc12S[n12S][8] = (*clust)[8];
    hc12S[n12S][9] = (*clust)[9];
    hc12S[n12S][10] = (*clust)[10];
    hc12S[n12S][11] = (*clust)[11];
    for (int i = 0; i < 11 ; ++i) {
        s12S[hc12S[n12S][i]] = 'B';
    }
    ++n12S;
}

