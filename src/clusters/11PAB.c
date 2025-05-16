#include <globals.h>
#include <tools.h>
#include "7PAB.h"
#include "7T.h"
#include "8PAA.h"
#include "11PAB.h"
//!  A 11PAB is a 10PAA with an additional particle
/*!
* 11PAB is a 10PAA cluster with an sp3a at one end. 
* Two of the sp3a particles form the bonded and nonbonded pair of one 
* of the 6Zs at either end of the 6Z.
*The other sp3a particle is the additional particle
*  Storage order: bonded, nbonded2, bonded2, nbonded, spindle1, spindle 2,
* 3 extra particles top to bottom, new particle)
*/

    void Clusters_Get11PAB() {
    for (int old_10PAA_id = 0; old_10PAA_id < n10PAA; ++old_10PAA_id) {
        int *old_10PAA_cluster = hc10PAA[old_10PAA_id];
        int clust_10PAA[10][2]; 
        clust_10PAA[0][0] = old_10PAA_cluster[3];
        clust_10PAA[1][0] = old_10PAA_cluster[2];
        clust_10PAA[2][0] = old_10PAA_cluster[5];
        clust_10PAA[3][0] = old_10PAA_cluster[9];
        clust_10PAA[4][0] = old_10PAA_cluster[0];
        clust_10PAA[5][0] = old_10PAA_cluster[4];
        clust_10PAA[6][0] = old_10PAA_cluster[1];
        clust_10PAA[7][0] = old_10PAA_cluster[6];
        clust_10PAA[8][0] = old_10PAA_cluster[7];
        clust_10PAA[9][0] = old_10PAA_cluster[8];
        
        clust_10PAA[0][1] = old_10PAA_cluster[7];
        clust_10PAA[1][1] = old_10PAA_cluster[5];
        clust_10PAA[2][1] = old_10PAA_cluster[2];
        clust_10PAA[3][1] = old_10PAA_cluster[8];
        clust_10PAA[4][1] = old_10PAA_cluster[6];
        clust_10PAA[5][1] = old_10PAA_cluster[1];
        clust_10PAA[6][1] = old_10PAA_cluster[4];
        clust_10PAA[7][1] = old_10PAA_cluster[0];
        clust_10PAA[8][1] = old_10PAA_cluster[3];
        clust_10PAA[9][1] = old_10PAA_cluster[9];

    for(int g = 0; g<2; ++g){
        int spindle_id_2 = clust_10PAA[5][g];
        int spindle_id = clust_10PAA[4][g];
        int nbonded_id = clust_10PAA[3][g];
        int bonded_id = clust_10PAA[0][g]; 
        int nbonded_id_2 = clust_10PAA[1][g];
        int bonded_id_2 = clust_10PAA[2][g]; 
        int final_part1 = clust_10PAA[6][g]; 
        int final_part2 = clust_10PAA[7][g]; 
        int final_part3 = clust_10PAA[8][g]; 
        int final_part4 = clust_10PAA[9][g]; 

        int arr_6Z[6] = {bonded_id, nbonded_id_2, bonded_id_2, nbonded_id, spindle_id, spindle_id_2};
        for (int new_3_id = 0; new_3_id < nsp3a; ++new_3_id ){ //iterate over all sp3a rings
            int *new_3_cluster = hcsp3a[new_3_id];
            if(nring_in_cluster(old_10PAA_cluster, new_3_cluster, 10) == 2){
                if(nring_bonded_nbonded(bonded_id, nbonded_id, new_3_cluster) == 2){ //two particles in ring are in 6Z
                   int bond_counter = check_ring_bonds(new_3_cluster, arr_6Z);
                    if (bond_counter == 9){
                        int new_part = get_new_particle_P2(new_3_cluster, bonded_id, nbonded_id);
                        if (check_unique11PAB(old_10PAA_cluster, new_part, 11) == 0){
                            if(check_5A(arr_6Z, new_part) == 0){
                                add_11PAB(bonded_id, nbonded_id_2, bonded_id_2, nbonded_id,spindle_id, spindle_id_2, final_part1, final_part2, final_part3, final_part4, new_part);
                            }
                            }
                        }
                    }
                }
            }
        }
    } 
} 

void add_11PAB(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int final_part2, int final_part3, int final_part4, int new_particle) {
    int clusSize = 11;
    //printf("new_particle %i\n", new_particle);
    if (n11PAB == m11PAB) {
        hc11PAB = resize_2D_int(hc11PAB, m11PAB, m11PAB + incrStatic, clusSize, -1);
        m11PAB = m11PAB + incrStatic;
    }
    hc11PAB[n11PAB][0] = bonded_id;
    hc11PAB[n11PAB][1] = nbonded_id_2;
    hc11PAB[n11PAB][2] = bonded_id_2;
    hc11PAB[n11PAB][3] = nbonded_id;
    hc11PAB[n11PAB][4] = spindle_id;
    hc11PAB[n11PAB][5] = spindle_id_2;
    hc11PAB[n11PAB][6] = final_part1;
    hc11PAB[n11PAB][7] = final_part2;
    hc11PAB[n11PAB][8] = final_part3;
    hc11PAB[n11PAB][9] = final_part4;
    hc11PAB[n11PAB][10] = new_particle;

    for (int i = 0; i < 11; ++i) {
        s11PAB[hc11PAB[n11PAB][i]] = 'B';
    }
    ++n11PAB;
}

int check_unique11PAB(const int *old_clust, int new_particle, int new_clust_size){
    int u;
    for (int old_11PAB_id = 0; old_11PAB_id < n11PAB; ++old_11PAB_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<10; ++r){
                //printf("%i %i\n", hc11P[old_11P_id][q],old_clust[r]);
                if(hc11PAB[old_11PAB_id][q] == old_clust[r]){
                    u += 1;
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc11PAB[old_11PAB_id][p] == new_particle){
                    u += 1;
                }
        }
    if(u == new_clust_size){
        return 1;            
    }
    }
    for (int old_11PAA_id = 0; old_11PAA_id < n11PAA; ++old_11PAA_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<10; ++r){
                //printf("%i %i\n", hc11P[old_11P_id][q],old_clust[r]);
                if(hc11PAA[old_11PAA_id][q] == old_clust[r]){
                    u += 1;
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc11PAA[old_11PAA_id][p] == new_particle){
                    u += 1;
                }
        }
    if(u == new_clust_size){
        return 1;            
    }
    }
    return 0;
}
