#include <globals.h>
#include <tools.h>
#include "14O.h"
//!  A 14O is a made of two intersecting 6A clusters and 2 4A clusters
// The two 6A clusters share 2 ring particles
// two rings and one spindle of each 6A cluster are in a 4A cluster
// No 4A cluster overlaps with these two 4A clusters


void Clusters_Get14O() {
    int *array = malloc(14 * sizeof(int));
    //get 6a
    for (int first_10O_id = 0; first_10O_id < n10O; ++first_10O_id) {
        int *first_10O_cluster = hcsp4c[first_10O_id];
        for (int second_6A_id = 0; second_6A_id < nsp4c; ++second_6A_id) {
            int *second_6A_cluster = hcsp4c[second_6A_id];
            //6a share two particles
            if(overlap_6A_14O(&array, first_10O_cluster, second_6A_cluster)== 1){
                
                //if(check_unique_14O(&array) == 0){  
                    //printf("%i %i %i %i %i %i %i %i %i %i\n", array[0], array[1],array[2],array[3],array[4],array[5],array[6],array[7],array[8],array[9]);
                    add_14O(&array);
                //}
            }
            }
        }
    }

int overlap_6A_14O(int **array, int *clust_10O, int *clust_6A2){
    int a = 0;
    int b = 10;
    int c;
    int d;

    for(int i = 0; i < 6; ++i){
        c = 0;
        for(int x = 4; x < 10; ++x){
            if (clust_10O[x] == clust_6A2[i]){
                a += 1;
                c +=1;
            } 
        }
        if(c == 0 && b < 14){ // particle is not shared
                (*array)[b] = clust_6A2[i];
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


void add_14O(int **new_14O_cluster) {
    int clusSize = 14;
    //printf("new_particle %i\n", new_particle);

    if (n14O == m14O) {
        hc14O = resize_2D_int(hc14O, m14O, m14O + incrStatic, clusSize, -1);
        m14O = m14O + incrStatic;
        
    }
    hc14O[n14O][0] = (*new_14O_cluster)[0];
    hc14O[n14O][1] = (*new_14O_cluster)[1];
    hc14O[n14O][2] = (*new_14O_cluster)[2];
    hc14O[n14O][3] = (*new_14O_cluster)[3];
    hc14O[n14O][4] = (*new_14O_cluster)[4];
    hc14O[n14O][5] = (*new_14O_cluster)[5];
    hc14O[n14O][6] = (*new_14O_cluster)[6];
    hc14O[n14O][7] = (*new_14O_cluster)[7];
    hc14O[n14O][8] = (*new_14O_cluster)[8];
    hc14O[n14O][9] = (*new_14O_cluster)[9];
    hc14O[n14O][10] = (*new_14O_cluster)[12];
    hc14O[n14O][11] = (*new_14O_cluster)[11];
    hc14O[n14O][12] = (*new_14O_cluster)[12];
    hc14O[n14O][13] = (*new_14O_cluster)[13];
    for (int i = 0; i < 14; ++i) {
        s14O[hc14O[n14O][i]] = 'B';
    }
    ++n14O;
}
