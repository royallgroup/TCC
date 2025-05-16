#include <globals.h>
#include <tools.h>
#include "11O.h"
#include "12O.h"

//!  A 11O is a made of two intersecting 6A clusters and 2 4A clusters
// The two 6A clusters share 2 ring particles
// two rings and one spindle of each 6A cluster are in a 4A cluster
// No 4A cluster overlaps with these two 4A clusters


void Clusters_Get11O() {

    int i,j,u,v;
    int *array = malloc(11 * sizeof(int));
    //get 6a
    for (int first_6A_id = 0; first_6A_id < nsp4c; ++first_6A_id) {
        int *first_6A_cluster = hcsp4c[first_6A_id];
        for (int second_6A_id = 0; second_6A_id < nsp4c; ++second_6A_id) {
            int *second_6A_cluster = hcsp4c[second_6A_id];
            //6a share three particles
            if(overlap_6A_6A_11O(&array, first_6A_cluster, second_6A_cluster)== 1){
                for(int first_4A_id = 0; first_4A_id < nsp3b; ++first_4A_id){
                    int *first_4A_cluster = hcsp3b[first_4A_id];
                    u = 0;
                    v = 0;
                    //first 6a shares 3 particles with 1st 4a, none of which are shares with the 2nd 6a
                    for(i = 0; i < 4; ++i){
                        for(j = 0; j < 6; ++j){
                                if(first_4A_cluster[i] == second_6A_cluster[j]){
                                    u +=1;
                                }
                            if(first_4A_cluster[i] == first_6A_cluster[j]){
                                v +=1;
                            }
                        }
                    }
                    if(u == 0 && v == 3){
                        //second 4a
                        for(int second_4A_id = 0; second_4A_id < nsp3b; ++second_4A_id){
                        int *second_4A_cluster = hcsp3b[second_4A_id];
                        if(first_4A_id != second_4A_id){
                            
                            if(overlap_4A_4A(first_4A_cluster, second_4A_cluster) == 1){
                                u = 0;
                                v = 0;
                                if(u ==0 && v == 0){
                                    u = 0;
                                    v = 0;
                                    for(i = 0; i < 4; ++i){
                                        for(j = 0; j < 6; ++j){
                                                if(second_4A_cluster[i] == first_6A_cluster[j]){
                                                    u +=1;
                                                }
                                            if(second_4A_cluster[i] == second_6A_cluster[j]){
                                                v +=1;
                                            }
                                        }
                                    }
                                    if(u == 0 && v == 3){
                                        array[3] = first_4A_cluster[0];
                                        array[4] = first_4A_cluster[1];
                                        array[5] = first_4A_cluster[2];
                                        array[6] = first_4A_cluster[3];

                                        array[7] = second_4A_cluster[0];
                                        array[8] = second_4A_cluster[1];
                                        array[9] = second_4A_cluster[2];
                                        array[10] = second_4A_cluster[3];
                                        if(check_unique_11O(&array) == 0){                                
                                            add_11O(&array);
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

int overlap_6A_4A_11O(int **array, int *clust_4A1, int *clust_4A2){
    int a = 0;
    int b = 0;
    int c = 0;
    int d = 0;
    for(int i = 3; i < 6; ++i){
        for(int j = 0; j < 4; ++j){
            if((*array)[i] == clust_4A1[j]){
                a += 1;
            }
            if((*array)[i] == clust_4A2[j]){
                b += 1;
            }
            if((*array)[i+3] == clust_4A1[j]){
                c += 1;
            }
            if((*array)[i+3] == clust_4A2[j]){
                d += 1;
            }
        }        
    }
    if(a == 3 && d == 3){
        if(b == 0 && c == 0){
            return 1;
        }
    }
    if(b == 3 && c == 3){
        if(a == 0 && d == 0){
            return 1;
        }
    }
    else{
        return 0;
    }
}

int overlap_6A_6A_11O(int **array, int *clust_6A1, int *clust_6A2){
    int a = 0;
    int c = 3;
    int e = 6;
    int b,d,i,x;
    for(i = 0; i < 6; ++i){
        b = 0;
        for(x = 0; x < 6; ++x){
            d = 0;
            if (clust_6A1[i] == clust_6A2[x]){
                if(a < 3){
                    (*array)[a] = clust_6A1[i];
                }  
                a += 1;
                b +=1;              
            }
        }
        if(b == 0 && c < 6){
            (*array)[c] = clust_6A1[i];
            c += 1;
        }
    }

    for(i = 0; i < 6; ++i){
        b = 0;
        for(x = 0; x < 6; ++x){
            if (clust_6A2[i] == clust_6A1[x]){
                b +=1;              
            }
        }
        if(b == 0 && c < 9){
            (*array)[c] = clust_6A2[i];
            c += 1;
        }
    }

    if(a == 3){
        return 1;
    }
    if(a == 2){
        return 2;
    }    
    else{
        return 0;
    }
}

int check_unique_11O(int **new_11O_cluster){
    int u;
    for (int old_11O_id = 0; old_11O_id < n11O; ++old_11O_id) {
        u = 0;
        for (int r = 0; r < 11; ++r){
            for (int q = 0; q < 11; ++q){
                if(hc11O[old_11O_id][r] == (*new_11O_cluster)[q]){
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
void add_11O(int **arr_11O) {
    int clusSize = 11;
    //printf("new_particle %i\n", new_particle);
    if (n11O == m11O) {
        hc11O = resize_2D_int(hc11O, m11O, m11O + incrStatic, clusSize, -1);
        m11O = m11O + incrStatic;
    }
    hc11O[n11O][0] = (*arr_11O)[0];
    hc11O[n11O][1] = (*arr_11O)[1];
    hc11O[n11O][2] = (*arr_11O)[2];
    hc11O[n11O][3] = (*arr_11O)[3];
    hc11O[n11O][4] = (*arr_11O)[4];
    hc11O[n11O][5] = (*arr_11O)[5];
    hc11O[n11O][6] = (*arr_11O)[6];
    hc11O[n11O][7] = (*arr_11O)[7];
    hc11O[n11O][8] = (*arr_11O)[8];
    hc11O[n11O][9] = (*arr_11O)[9];
    hc11O[n11O][10] = (*arr_11O)[10];

    for (int i = 0; i < 11; ++i) {
        s11O[hc11O[n11O][i]] = 'B';
    }
    ++n11O;
}

void get_11O(int **arr_11O, int *first_6A_cluster, int *second_6A_cluster, int *first_4A_cluster, int *second_4A_cluster){
    int u, v, w, idx2;
    //put the first 6A into the cluster
    for (int i = 0; i < 6; ++i){
        (*arr_11O)[0] = first_6A_cluster[0];
        (*arr_11O)[1] = first_6A_cluster[1];
        (*arr_11O)[2] = first_6A_cluster[2];
        (*arr_11O)[3] = first_6A_cluster[3];
        (*arr_11O)[4] = first_6A_cluster[4];
        (*arr_11O)[5] = first_6A_cluster[5];
    }
    int idx = 6;
    //find the two non shared ring particles
    for(int k = 0; k < 6; ++k){
        u = 0;
        for(int j = 0; j < 6; ++j){
            if((*arr_11O)[j] == second_6A_cluster[k]){
                u += 1;
            }
        }
        if(u == 0){
            (*arr_11O)[idx] = second_6A_cluster[k];
            idx += 1;
        }
    }
    idx2 = idx;
    for(int a = 0; a < 4; ++a){
        v = 0;
        for(int b = 0; b < idx2; ++b){
            if((*arr_11O)[b] == first_4A_cluster[a]){
                v += 1;

            }
        }
        if(v == 0){
            (*arr_11O)[idx] = first_4A_cluster[a];
            idx +=1;
        }
    }
    idx2 = idx;
    for(int c = 0; c < 4; ++c){
        w = 0;
        for(int d = 0; d < idx2; ++d){
            if((*arr_11O)[d] == second_4A_cluster[c]){
                w += 1;
            }
        }
        if(w == 0){
            (*arr_11O)[idx] = second_4A_cluster[c];
            idx +=1;
        }
    }
    return;
}