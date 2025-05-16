#include <globals.h>
#include <tools.h>
#include "10O.h"
#include "11O.h"
#include "12O.h"
//!  A 10O is a made of two intersecting 6A clusters and 2 4A clusters
// The two 6A clusters share 2 ring particles
// two rings and one spindle of each 6A cluster are in a 4A cluster
// No 4A cluster overlaps with these two 4A clusters


void Clusters_Get10O() {
    int *array = malloc(10 * sizeof(int));
    //get 6a
    for (int first_6A_id = 0; first_6A_id < nsp4c; ++first_6A_id) {
        int *first_6A_cluster = hcsp4c[first_6A_id];
        for (int second_6A_id = 0; second_6A_id < nsp4c; ++second_6A_id) {
            int *second_6A_cluster = hcsp4c[second_6A_id];
            //6a share two particles
            if(overlap_6A_6A_10O(&array, first_6A_cluster, second_6A_cluster)== 1){
                
                if(check_unique_10O(&array) == 0){  
                    //printf("%i %i %i %i %i %i %i %i %i %i\n", array[0], array[1],array[2],array[3],array[4],array[5],array[6],array[7],array[8],array[9]);
                    add_10O(&array);
                }
            }
            }
        }
    }

int overlap_6A_6A_10O(int **array, int *clust_6A1, int *clust_6A2){
    int a = 0;
    int b = 6;
    int c;
    int d;
    (*array)[0] = clust_6A2[0];
    (*array)[1] = clust_6A2[1];
    (*array)[2] = clust_6A2[2];
    (*array)[3] = clust_6A2[3];
    (*array)[4] = clust_6A2[4];
    (*array)[5] = clust_6A2[5];

    for(int i = 0; i < 6; ++i){
        c = 0;
        for(int x = 0; x < 6; ++x){
            if (clust_6A1[i] == clust_6A2[x]){
                //if(a < 2){
                   // (*array)[a] = clust_6A1[i];
               // }
                a += 1;
                c +=1;
            } 
        }
        if(c == 0 && b < 10){ // particle is not shared
                (*array)[b] = clust_6A1[i];
                b += 1;
        }

    }
   
    if(a == 2){
        return 1;
    }
    else{
        return 0;
    }
}

int check_unique_10O(int **new_10O_cluster){
    int u;
    for (int old_10O_id = 0; old_10O_id < n10O; ++old_10O_id) {
        u = 0;
        for (int r = 0; r < 10; ++r){
            for (int q = 0; q < 10; ++q){
                if(hc10O[old_10O_id][q] == (*new_10O_cluster)[r]){
                    u += 1;
                }
            }

        }  
        
        if(u == 10){
            return 1;           
        }
    }
    return 0;
}

void add_10O(int **new_10O_cluster) {
    int clusSize = 10;
    //printf("new_particle %i\n", new_particle);

    if (n10O == m10O) {
        hc10O = resize_2D_int(hc10O, m10O, m10O + incrStatic, clusSize, -1);
        m10O = m10O + incrStatic;
        
    }
    hc10O[n10O][0] = (*new_10O_cluster)[0];
    hc10O[n10O][1] = (*new_10O_cluster)[1];
    hc10O[n10O][2] = (*new_10O_cluster)[2];
    hc10O[n10O][3] = (*new_10O_cluster)[3];
    hc10O[n10O][4] = (*new_10O_cluster)[4];
    hc10O[n10O][5] = (*new_10O_cluster)[5];
    hc10O[n10O][6] = (*new_10O_cluster)[6];
    hc10O[n10O][7] = (*new_10O_cluster)[7];
    hc10O[n10O][8] = (*new_10O_cluster)[8];
    hc10O[n10O][9] = (*new_10O_cluster)[9];
    for (int i = 0; i < 10; ++i) {
        s10O[hc10O[n10O][i]] = 'B';
    }
    ++n10O;
}
