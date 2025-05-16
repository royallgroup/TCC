#include <globals.h>
#include <tools.h>
#include "12O.h"
#include "11O.h"
//!  A 12O is a made of two intersecting 6A clusters and 2 4A clusters
// The two 6A clusters share 2 ring particles
// two rings and one spindle of each 6A cluster are in a 4A cluster
// No 4A cluster overlaps with these two 4A clusters


void Clusters_Get12O() {

    int i,j,u,v,w;
    int *array = malloc(12 * sizeof(int));
    //get 6a
    for (int first_6A_id = 0; first_6A_id < nsp4c; ++first_6A_id) {
        int *first_6A_cluster = hcsp4c[first_6A_id];
        for (int second_6A_id = 0; second_6A_id < nsp4c; ++second_6A_id) {
            int *second_6A_cluster = hcsp4c[second_6A_id];
            //6a share two particles
            if(overlap_6A_6A_11O(&array, first_6A_cluster, second_6A_cluster)== 2){
                //first 4a
                for(int first_4A_id = 0; first_4A_id < nsp3b; ++first_4A_id){
                    int *first_4A_cluster = hcsp3b[first_4A_id];
                    if(overlap_4A_12O(first_4A_cluster, first_6A_cluster, second_6A_cluster) == 1){
                        //printf("%i %i\n", array[0], array[1]);
                        //second 4a
                        for(int second_4A_id = 0; second_4A_id < nsp3b; ++second_4A_id){
                        int *second_4A_cluster = hcsp3b[second_4A_id];
                        if(first_4A_id != second_4A_id){
                            
                            if(overlap_4A_4A(first_4A_cluster, second_4A_cluster) == 1){
                                if(overlap_4A_12O(second_4A_cluster, second_6A_cluster, first_6A_cluster)== 1){
                                    int k = 0;
                                //the two 4as do not share any particles with any other 4a
                                for(int third_4A_id = 0; third_4A_id < nsp3b; ++third_4A_id){   
                                    v = 0;
                                    w = 0;
                                    int *third_4A_cluster = hcsp3b[third_4A_id];
                                    u = 0;
                                    //third 4a shares 2 particles with the shared 6a particles
                                    for(i = 0; i < 4; ++i){
                                        if(array[0] == third_4A_cluster[i] || array[1] == third_4A_cluster[i]){
                                            u +=1;
                                        }
                                    }
                                    if(u == 2){
                                        u = 0;
                                        for(i = 0; i < 4; ++i){
                                            for(j = 0; j < 4; ++j){
                                                if(third_4A_cluster[i] == second_4A_cluster[j]){
                                                    v += 1;
                                                }
                                            if(third_4A_cluster[i]== first_4A_cluster[j]){
                                                    w +=1;
                                                }
                                            }                                            
                                        }
                                    }
                                    if(v != 0 && w != 0){
                                        k += 1;
                                    }
                                }
                                
                                if(k == 0){
                                    array[2] = first_4A_cluster[0];
                                    array[3] = first_4A_cluster[1];
                                    array[4] = first_4A_cluster[2];
                                    array[5] = first_4A_cluster[3];

                                    array[6] = second_4A_cluster[0];
                                    array[7] = second_4A_cluster[1];
                                    array[8] = second_4A_cluster[2];
                                    array[9] = second_4A_cluster[3];     
                                    for(i = 0; i < 6; ++i){
                                        u = 0;
                                        v = 0;
                                        for(j = 0; j < 10; ++j){
                                            if(array[j] == first_6A_cluster[i]){
                                                u +=1;
                                            }
                                            if(array[j] == second_6A_cluster[i]){
                                                v +=1;
                                            }
                                        }
                                        if(u == 0){
                                            array[10] = first_6A_cluster[i];
                                        }
                                        if(v == 0){
                                            array[11] = second_6A_cluster[i];
                                        }
                                    }                                       

                                    if(check_unique_12O(&array) == 0){                                
                                        add_12O(&array);
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

int overlap_4A_12O(int *clust_4A, int *clust_6A1, int *clust_6A2){
    int u = 0;
    int v = 0;
    for(int i = 0; i < 4; ++i){
        for(int j = 0; j < 6; ++j){
            if(clust_4A[i] == clust_6A1[j]){
                u +=1;
                }
            if(clust_4A[i] == clust_6A2[j]){
                v +=1;
            }
        }
    }
    if(v == 0 && u == 3){
        return 1;
    }
    else{
        return 0;
    }
}

//put the particles into the array
void get_12O(int **array, int *first_4A_cluster, int *second_4A_cluster, int *first_6A_cluster, int *second_6A_cluster){
    //put the first 6A into the cluster
    int u = 6;
    int i,j,a;
    int v = 6;
    for (i = 0; i < 6; ++i){
        a = 0;
        for (j = 0; j < u; ++j){
            if((*array)[j] == first_6A_cluster[i]){
                a += 1;
            }
        }
        if(a == 0){
            (*array)[v] = first_6A_cluster[i];
            v += 1;
        }
    }
    
    u = v;
    for (i = 0; i < 6; ++i){
        a = 0;
        for (j = 0; j < u; ++j){
            if((*array)[j] == first_6A_cluster[i]){
                a += 1;
            }
        }
        if(a == 0){
            
            (*array)[v] = first_6A_cluster[i];
            v += 1;
        }
    }
    
    u = v;
    for (i = 0; i < 4; ++i){
        a = 0;
        for (j = 0; j < u; ++j){
            if((*array)[j] == first_4A_cluster[i]){
                a += 1;
            }
        }
        if(a == 0){
            (*array)[v] = first_4A_cluster[i];
            v += 1;
        }
    }    
    u = v;
    for (i = 0; i < 4; ++i){
        a = 0;
        for (j = 0; j < u; ++j){
            if((*array)[j] == second_4A_cluster[i]){
                a += 1;
            }
        }
        if(a == 0){
            (*array)[v] = second_4A_cluster[i];
            v += 1;
        }
    }
    return;
}

int overlap_6A_ring(int *first_6A_cluster, int *second_6A_cluster){
    int u = 0;
    int v = 0;
    for (int i = 0; i < 6; ++i){
        for(int j = 0; j < 4; ++j){
            if(first_6A_cluster[i] == second_6A_cluster[j]){
                u += 1;
            }
        }
    }
    for (int y = 4; y < 6; ++y){
        for(int x = 4; x < 6; ++x){
            if(first_6A_cluster[y] == second_6A_cluster[x]){
                v += 1;
            }
        }
    }
    if(u == 2 && v == 0){
        return 1;
    }
    return 0;
}

int overlap_6A_6A_12O(int **array, int *clust_6A1, int *clust_6A2){
    int a = 0;
    for(int i = 0; i < 6; ++i){
        for(int x = 0; x < 6; ++x){
            if (clust_6A1[i] == clust_6A2[x]){
                if(a < 2){
                    (*array)[a] = clust_6A1[i];
                }
                a += 1;
                
            }
        }
    }
    if(a == 2){
        return 1;
    }
    else{
        return 0;
    }
}

int overlap_6A_6A_4A_4A(int **array,int *clust_6A1, int *clust_4A1,int *clust_6A2, int *clust_4A2){
    int a = 0;
    int l = 0;
    int b,c, m,n,p1,p2,p3,p4;
    b = 0;
    c = 0;
    n = 0;
    m = 0;
    for(int i = 0; i < 4; ++i){
        for(int x = 0; x < 2; ++x){
            if ((*array)[x] == clust_4A1[i]){
                a += 1;
                }
            if ((*array)[x] == clust_4A2[i]){
                l += 1;
                }
        }
        for(int j = 0; j < 6; ++j){
            if(clust_6A1[j] == clust_4A1[i] && clust_4A1[i] != (*array)[0] && clust_4A1[i] != (*array)[1]){
                b += 1;
                p1 = clust_6A1[j];
            }
            if(clust_6A2[j] == clust_4A1[i] && clust_4A1[i] != (*array)[0] && clust_4A1[i] != (*array)[1]){
                c += 1;
                p2 = clust_6A2[j];
            }

            if(clust_6A1[j] == clust_4A2[i] && clust_4A2[i] != (*array)[0] && clust_4A2[i] != (*array)[1]){
                m += 1;
                p3 = clust_6A1[j];
            }
            if(clust_6A2[j] == clust_4A2[i] && clust_4A2[i] != (*array)[0] && clust_4A2[i] != (*array)[1]){
                n += 1;
                p4 = clust_6A2[j];
            }
        }
    }

    if(a == 2 && l == 2 && b == 1 && c == 1 && m == 1 && n == 1){
        if(p1 != p3 && p1 != p4 && p2 != p3 && p2 != p4){  
            (*array)[2] = p1;
            (*array)[3] = p2;
            (*array)[4] = p3;
            (*array)[5] = p4;
            return 1;
        }
    }
    else{
        return 0;
    }
}

//6A shares a spindle and two rings with the 4a
int overlap_6A_4A(int *clust_6A, int *clust_4A){
    int a = 0;
    for(int i = 0; i < 6; ++i){
        for(int x = 0; x < 4; ++x){
            if (clust_6A[i] == clust_4A[x]){
                a += 1;
            }
        }
    }
    if(a == 3){
        return 3;
    }
    if(a == 2){
        return 2;
    }
    else{
        return 0;
    }
}

int overlap_4A_4A(int *clust_4A1, int *clust_4A2){
    int a = 0;
    for(int i = 0; i < 4; ++i){
        for(int j = 0; j < 4; ++j){
            if (clust_4A1[i] == clust_4A2[j]){
                a += 1;
            }
        }
    }
    if(a == 0){
        return 1;
    }
    if(a == 2){
        return 2;
    }
    if(a == 4){
        return 4;
    }
    else{
        return 0;
    }
}

int shared_4A(int *third_4A_cluster, int *fourth_4A_cluster, int *first_4A_cluster, int *second_4A_cluster){
    int u, v,w,x;
    u = 0;
    v = 0;
    w = 0;
    x = 0;
    //iterate over third 4a
    for(int k = 0; k < 4; ++k){
        for(int i = 0; i < 4; ++i){
                if(first_4A_cluster[i] == fourth_4A_cluster[k]){
                    u += 1;
                }   
                if(first_4A_cluster[i] == third_4A_cluster[k]){
                    v += 1;
                }     
                if(second_4A_cluster[i] == third_4A_cluster[k]){
                    w += 1;
                }   
                if(second_4A_cluster[i] == fourth_4A_cluster[k]){
                    x += 1;
                }                        
            }
        }
    //printf("%i %i %i %i \n",u,v,w,x);
    if(u == 1  && v == 0){
        if(w == 1  && x == 0){
        return 1;
        }
    }
    if(u == 0  && v == 1){
        if(w == 0  && x == 1){
            return 1;
        }
    }
    return 0;
}

int check_unique_12O(int **new_12O_cluster){
    int u;
    for (int old_12O_id = 0; old_12O_id < n12O; ++old_12O_id) {
        u = 0;
        for (int r = 0; r < 12; ++r){
            for (int q = 0; q < 12; ++q){
                if(hc12O[old_12O_id][q] == (*new_12O_cluster)[r]){
                    u += 1;
                }
            }

        }  
        if(u >12){
            printf("boop\n");
        }  
        if(u == 12){
            return 1;           
        }
    }
    return 0;
}

void add_12O(int **new_12O_cluster) {
    int clusSize = 12;
    //printf("new_particle %i\n", new_particle);
    if (n12O == m12O) {
        hc12O = resize_2D_int(hc12O, m12O, m12O + incrStatic, clusSize, -1);
        m12O = m12O + incrStatic;
    }
    hc12O[n12O][0] = (*new_12O_cluster)[0];
    hc12O[n12O][1] = (*new_12O_cluster)[1];
    hc12O[n12O][2] = (*new_12O_cluster)[2];
    hc12O[n12O][3] = (*new_12O_cluster)[3];
    hc12O[n12O][4] = (*new_12O_cluster)[4];
    hc12O[n12O][5] = (*new_12O_cluster)[5];
    hc12O[n12O][6] = (*new_12O_cluster)[6];
    hc12O[n12O][7] = (*new_12O_cluster)[7];
    hc12O[n12O][8] = (*new_12O_cluster)[8];
    hc12O[n12O][9] = (*new_12O_cluster)[9];
    hc12O[n12O][10] = (*new_12O_cluster)[10];
    hc12O[n12O][11] = (*new_12O_cluster)[11];
    for (int i = 0; i < 12; ++i) {
        s12O[hc12O[n12O][i]] = 'B';
    }
    ++n12O;
}
