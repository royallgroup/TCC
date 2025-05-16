#include <globals.h>
#include <tools.h>
#include "8O.h"
#include "12O.h"

//!  A 8O is a made of two intersecting 6A clusters and 2 4A clusters
// The two 6A clusters share 2 ring particles
// two rings and one spindle of each 6A cluster are in a 4A cluster
// No 4A cluster overlaps with these two 4A clusters


void Clusters_Get8O() {
    int *array = malloc(8 * sizeof(int));
    int u;
    for (int first_6A_id = 0; first_6A_id < nsp4c; ++first_6A_id) {
        int *first_6A_cluster = hcsp4c[first_6A_id];
        for(int first_4A_id = 0; first_4A_id < nsp3b; ++first_4A_id){
            int *first_4A_cluster = hcsp3b[first_4A_id];
            u = 0;
            for(int i = 0; i < 4; ++i){
                for(int j = 0; j < 6; ++j){
                    if(first_6A_cluster[j] == first_4A_cluster[i]){
                    u += 1;
                }
                }                
            }
            //printf("%i\n", u);
            if(u == 3){
                for(int second_4A_id = first_4A_id; second_4A_id < nsp3b; ++second_4A_id){
                    int *second_4A_cluster = hcsp3b[second_4A_id];
                    //printf("%i %i %i\n", first_6A_id, first_4A_id, second_4A_id);
                    if(overlap_4A_4A(first_4A_cluster, second_4A_cluster) == 1){
                        if(overlap_6A_4A_8O(&array, first_6A_cluster, first_4A_cluster,second_4A_cluster) == 1){
                            if(check_unique_8O(&array) == 0){
                                add_8O(&array);
                            }
                        }
                    }
                }
            }
        }
    }
}

int overlap_6A_4A_8O(int **array,const int *clust_6A,const int *clust_4A1,const int *clust_4A2){
    int u = 0;
    int v;
    for(int i = 0; i < 4; ++i){
        for(int j = 0; j < 4; ++j){
            if(clust_4A1[i] == clust_4A2[j]){
                u += 1;
            }
        }    
    }
    if(u == 0){
        u = 0;
        for(int i = 0; i < 4; ++i){
            v = 0;
            for(int j = 0; j < 6; ++j){
                if(clust_4A1[i] == clust_6A[j]){
                    (*array)[u] = clust_4A1[i];
                    u += 1;
                    v += 1;
                }
            }    
            if(v == 0){
                (*array)[6] = clust_4A1[i];
            }
        }
        if(u == 3){
            u = 0;
            for(int i = 0; i < 4; ++i){
                v = 0;
                 for(int j = 0; j < 6; ++j){
                    if(clust_4A2[i] == clust_6A[j]){
                        (*array)[u+3] = clust_4A2[i];
                        u += 1;
                        v += 1;
                    }
                }
                if(v == 0){
                (*array)[7] = clust_4A2[i];
                }    
            }
        }
        if(u == 3){
            return 1;
        }
    }
    return 0;
}


int check_unique_8O(int **array){
    int u;
    for (int old_8O_id = 0; old_8O_id < n8O; ++old_8O_id) {
        u = 0;
        for (int r = 0; r < 8; ++r){
            for (int q = 0; q < 8; ++q){
                if(hc8O[old_8O_id][r] == (*array)[q]){
                    u += 1;
                }
            }

        }
    
    if(u == 8){
        return 1;           
    }
    }
    return 0;
}

void add_8O(int **array) {
    int clusSize = 8;
    if (n8O == m8O) {
        hc8O = resize_2D_int(hc8O, m8O, m8O + incrStatic, clusSize, -1);
        m8O = m8O + incrStatic;
    }
    hc8O[n8O][0] = (*array)[0];
    hc8O[n8O][1] = (*array)[1];
    hc8O[n8O][2] = (*array)[2];
    hc8O[n8O][3] = (*array)[3];
    hc8O[n8O][4] = (*array)[4];
    hc8O[n8O][5] = (*array)[5];
    hc8O[n8O][6] = (*array)[6];
    hc8O[n8O][7] = (*array)[7];

    for (int i = 0; i < 8; ++i) {
        s8O[hc8O[n8O][i]] = 'B';
    }
    ++n8O;
}
