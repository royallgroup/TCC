#include <globals.h>
#include <tools.h>
#include "10S.h"
#include "11S.h"
#include "8PAB.h"
#include "7PAB.h"
#include "9S.h"
#include "7T.h"
//!  A 10S is a 7A with an additional particle
/*!
* A 7A cluster contains two 6Z clusters at the top and bottom of the cluster. 
* A 10S takes these two 6Z clusters and tries to make a 7T from the cluster and an extra particle.
*  Cluster output: BBBBBBBBBB
*  Storage order: original 7A particles x 9, new 5A spindle)
*/

void Clusters_Get10S() {
    int new_part;
    for (int first_7A_id = 0; first_7A_id < nsp5c; ++first_7A_id) {
        int *first_7A_cluster = hcsp5c[first_7A_id];
        for (int spindle_pointer = 0; spindle_pointer < 2; ++spindle_pointer) {
            int primary_spindle = first_7A_cluster[5 + spindle_pointer];
            for (int new_5A_pointer = 0; new_5A_pointer < nmem_sp3c[primary_spindle]; ++new_5A_pointer) {
                int new_5A_id = mem_sp3c[primary_spindle][new_5A_pointer];
                int *new_5A_cluster = hcsp3c[new_5A_id];
                    if (is_particle_spindle_of_5A(primary_spindle, new_5A_cluster) == 1) { // the 8A and 5A share a spindle
                        if(common_ring(first_7A_cluster, new_5A_cluster, 5, 3) == 1){
                                if(new_5A_cluster[3] == primary_spindle){
                                    new_part == new_5A_cluster[4];
                                }
                                if(new_5A_cluster[4] == primary_spindle){
                                    new_part == new_5A_cluster[3];
                                }
                                int *clust_10S;
                                get_10S(&clust_10S, first_7A_cluster, new_5A_cluster);
                                add_10S(&clust_10S);
                                free(clust_10S);
                        }
                    }
                }
            }
        }
    }

void get_10S(int **array, const int *clust1, const int *clust2){
    *array = malloc(10 * sizeof(int));
    if (*array == NULL){
        return;
    }
    for (int i = 0 ; i < 7 ; i++){
        (*array)[i] = clust1[i];
        }
    int count;
    int index = 7;
    for(int j = 0; j < 5; ++j){
        count = 0;
        for(int k = 0; k < 7; ++ k){
            if((*array)[k] == clust2[j]){
                count += 1;
            }
        }
        if(count == 0){
                (*array)[index] = clust2[j];
                index += 1;
            }
    }

}

void add_10S(int **clust) {
    int clusSize = 11;
    //printf("new_particle %i\n", new_particle);
    if (n10S == m10S) {
        hc10S = resize_2D_int(hc10S, m10S, m10S + incrStatic, clusSize, -1);
        m10S = m10S + incrStatic;
    }
    hc10S[n10S][0] = (*clust)[0];
    hc10S[n10S][1] = (*clust)[1];
    hc10S[n10S][2] = (*clust)[2];
    hc10S[n10S][3] = (*clust)[3];
    hc10S[n10S][4] = (*clust)[4];
    hc10S[n10S][5] = (*clust)[5];
    hc10S[n10S][6] = (*clust)[6];
    hc10S[n10S][7] = (*clust)[7];
    hc10S[n10S][8] = (*clust)[8];
    hc10S[n10S][9] = (*clust)[9];

    for (int i = 0; i < 10; ++i) {
        s10S[hc10S[n10S][i]] = 'B';
    }
    ++n10S;
}