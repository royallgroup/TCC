
#include <globals.h>
#include <tools.h>
#include "10PAA.h"
#include "10PAB.h"
#include "7T.h"
#include "7PAB.h"
//!  A 10PAB is a 9PAA with an additional particle
/*!
* 10PAB is a 9PAA cluster with an sp3a at one end. 
* Two of the sp3a particles form the bonded and nonbonded pair of one 
* of the 6Zs at either end of the 6Z.
*The other sp3a particle is the additional particle
*  Storage order: bonded, nbonded2, bonded2, nbonded, spindle1, spindle 2,
* 3 extra particles top to bottom, new particle)
*/

void Clusters_Get10PAB() {
for (int old_9PAA_id = 0; old_9PAA_id < n9PAA; ++old_9PAA_id) {
    int *old_9PAA_cluster = hc9PAA[old_9PAA_id];
    int clust_7T[9][2]; 
    clust_7T[0][0] = old_9PAA_cluster[3];
    clust_7T[1][0] = old_9PAA_cluster[2];
    clust_7T[2][0] = old_9PAA_cluster[5];
    clust_7T[3][0] = old_9PAA_cluster[8];
    clust_7T[4][0] = old_9PAA_cluster[0];
    clust_7T[5][0] = old_9PAA_cluster[4];
    clust_7T[6][0] = old_9PAA_cluster[1];
    clust_7T[7][0] = old_9PAA_cluster[6];
    clust_7T[8][0] = old_9PAA_cluster[7];
    
    clust_7T[0][1] = old_9PAA_cluster[6];
    clust_7T[1][1] = old_9PAA_cluster[4];
    clust_7T[2][1] = old_9PAA_cluster[5];
    clust_7T[3][1] = old_9PAA_cluster[7];
    clust_7T[4][1] = old_9PAA_cluster[1];
    clust_7T[5][1] = old_9PAA_cluster[2];
    clust_7T[6][1] = old_9PAA_cluster[0];
    clust_7T[7][1] = old_9PAA_cluster[3];
    clust_7T[8][1] = old_9PAA_cluster[8];

    for(int g = 0; g<2; ++g){
        int spindle_id_2 = clust_7T[5][g];
        int spindle_id = clust_7T[4][g];
        int nbonded_id = clust_7T[3][g];
        int bonded_id = clust_7T[0][g]; 
        int nbonded_id_2 = clust_7T[1][g];
        int bonded_id_2 = clust_7T[2][g]; 
        int final_part1 = clust_7T[6][g]; 
        int final_part2 = clust_7T[7][g]; 
        int final_part3 = clust_7T[8][g]; 
        int arr_6Z[6] = {bonded_id, nbonded_id_2, bonded_id_2, nbonded_id, spindle_id, spindle_id_2};
        for (int new_3_id = 0; new_3_id < nsp3a; ++new_3_id ){ //iterate over all sp3a rings
            int *new_3_cluster = hcsp3a[new_3_id];
            if(nring_in_cluster(old_9PAA_cluster, new_3_cluster, 9) == 2){
                if(nring_bonded_nbonded(bonded_id, nbonded_id, new_3_cluster) == 2){ //two particles in ring are in 6Z
                   int bond_counter = check_ring_bonds(new_3_cluster, arr_6Z);
                    if (bond_counter == 9){
                        int new_part = get_new_particle_P2(new_3_cluster, bonded_id, nbonded_id);
                        //printf("6Z %i %i %i %i %i %i \n", bonded_id, nbonded_id_2, bonded_id_2, nbonded_id, spindle_id, spindle_id_2);
                        //printf("new %i\n", new_part);
                        if (check_unique10PAB(old_9PAA_cluster, new_part, 10) == 0){
                            if(check_5A(arr_6Z, new_part) == 0){
                            add_10PAB(bonded_id, nbonded_id_2, bonded_id_2, nbonded_id,spindle_id, spindle_id_2, final_part1, final_part2, final_part3, new_part);
                            }
                            }
                        }
                    }
                }
            }
        }
    }
}

void add_10PAB(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int final_part2, int final_part3, int new_particle) {
    int clusSize = 10;
    if (n10PAB == m10PAB) {
        hc10PAB = resize_2D_int(hc10PAB, m10PAB, m10PAB + incrStatic, clusSize, -1);
        m10PAB = m10PAB + incrStatic;
    }
    hc10PAB[n10PAB][0] = bonded_id;
    hc10PAB[n10PAB][1] = nbonded_id_2;
    hc10PAB[n10PAB][2] = bonded_id_2;
    hc10PAB[n10PAB][3] = nbonded_id;
    hc10PAB[n10PAB][4] = spindle_id;
    hc10PAB[n10PAB][5] = spindle_id_2;
    hc10PAB[n10PAB][6] = final_part1;
    hc10PAB[n10PAB][7] = final_part2;
    hc10PAB[n10PAB][8] = final_part3;
    hc10PAB[n10PAB][9] = new_particle;

    for (int i = 0; i < 10; ++i) {
        s10PAB[hc10PAB[n10PAB][i]] = 'B';
    }
    ++n10PAB;
}

int check_unique10PAB(const int *old_clust, int new_particle, int new_clust_size){
    int u;
    for (int old_10PAB_id = 0; old_10PAB_id < n10PAB; ++old_10PAB_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<9; ++r){
                //printf("%i %i\n", hc11P[old_11P_id][q],old_clust[r]);
                if(hc10PAB[old_10PAB_id][q] == old_clust[r]){
                    u += 1;
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc10PAB[old_10PAB_id][p] == new_particle){
                    u += 1;
                }
        }
    if(u == new_clust_size){
        return 1;            
    }
    }
    for (int old_10PAA_id = 0; old_10PAA_id < n10PAA; ++old_10PAA_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<9; ++r){
                //printf("%i %i\n", hc11P[old_11P_id][q],old_clust[r]);
                if(hc10PAA[old_10PAA_id][q] == old_clust[r]){
                    u += 1;
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc10PAA[old_10PAA_id][p] == new_particle){
                    u += 1;
                }
        }
    if(u == new_clust_size){
        return 1;            
    }
    }
    return 0;
}
