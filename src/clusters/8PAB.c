
#include <globals.h>
#include <tools.h>
#include "8PAB.h"
#include "7T.h"
#include "7PAB.h"
//!  A 8PAB is a 6Z with an additional particle
/*!
* 8PAB is a 6Z cluster with an sp3a at one end. 
* Two of the sp3a particles form the bonded and nonbonded pair of one 
* of the 6Zs at either end of the 6Z.
*The other sp3a particle is the additional particle
*  Storage order: bonded, nbonded2, bonded2, nbonded, spindle1, spindle 2,
* new particle)
*/

void Clusters_Get8PAB() {
for (int old_7T_id = 0; old_7T_id < n7T_a; ++old_7T_id) {
    int *old_7T_cluster = hc7T_a[old_7T_id];
    int clust_7T[7][2]; 
    clust_7T[0][0] = old_7T_cluster[2];
    clust_7T[1][0] = old_7T_cluster[3];
    clust_7T[2][0] = old_7T_cluster[0];
    clust_7T[3][0] = old_7T_cluster[1];
    clust_7T[4][0] = old_7T_cluster[5];
    clust_7T[5][0] = old_7T_cluster[4];
    clust_7T[6][0] = old_7T_cluster[6];
    
    clust_7T[0][1] = old_7T_cluster[3];
    clust_7T[1][1] = old_7T_cluster[2];
    clust_7T[2][1] = old_7T_cluster[5];
    clust_7T[3][1] = old_7T_cluster[6];
    clust_7T[4][1] = old_7T_cluster[0];
    clust_7T[5][1] = old_7T_cluster[4];
    clust_7T[6][1] = old_7T_cluster[1];
    for(int g = 0; g<2; ++g){
        int spindle_id_2 = clust_7T[5][g];
        int spindle_id = clust_7T[4][g];
        int nbonded_id = clust_7T[3][g];
        int bonded_id = clust_7T[0][g]; 
        int nbonded_id_2 = clust_7T[1][g];
        int bonded_id_2 = clust_7T[2][g]; 
        int final_part1 = clust_7T[6][g]; 
        int arr_6Z[6] = {bonded_id, nbonded_id_2, bonded_id_2, nbonded_id, spindle_id, spindle_id_2};
        for (int new_3_id = 0; new_3_id < nsp3a; ++new_3_id ){ //iterate over all sp3a rings
            int *new_3_cluster = hcsp3a[new_3_id];
            if(nring_in_cluster(old_7T_cluster, new_3_cluster, 7) == 2){
                if(nring_bonded_nbonded(bonded_id, nbonded_id, new_3_cluster) == 2){ //two particles in ring are in 6Z
                   int bond_counter = check_ring_bonds(new_3_cluster, arr_6Z);
                    if (bond_counter == 9){
                        int new_part = get_new_particle_P2(new_3_cluster, bonded_id, nbonded_id);
                        if (check_unique8PAB(old_7T_cluster, new_part, 8) == 0){
                            //printf("bodop\n");
                            if(check_5A(arr_6Z, new_part) == 0){
                                
                            add_8PAB(bonded_id, nbonded_id_2, bonded_id_2, nbonded_id,spindle_id, spindle_id_2, final_part1, new_part);
                            }
                            }
                        }
                    }
                }
            }
        }
    }
}

void add_8PAB(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int new_particle) {
    int clusSize = 8;
    if (n8PAB == m8PAB) {
        hc8PAB = resize_2D_int(hc8PAB, m8PAB, m8PAB + incrStatic, clusSize, -1);
        m8PAB = m8PAB + incrStatic;
    }
    hc8PAB[n8PAB][0] = bonded_id;
    hc8PAB[n8PAB][1] = nbonded_id_2;
    hc8PAB[n8PAB][2] = bonded_id_2;
    hc8PAB[n8PAB][3] = nbonded_id;
    hc8PAB[n8PAB][4] = spindle_id;
    hc8PAB[n8PAB][5] = spindle_id_2;
    hc8PAB[n8PAB][6] = final_part1;
    hc8PAB[n8PAB][7] = new_particle;

    for (int i = 0; i < 8; ++i) {
        s8PAB[hc8PAB[n8PAB][i]] = 'B';
    }
    ++n8PAB;
}

int check_unique8PAB(const int *old_clust, int new_particle, int new_clust_size){
    int u;
    for (int old_8PAB_id = 0; old_8PAB_id < n8PAB; ++old_8PAB_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<7; ++r){
                if(hc8PAB[old_8PAB_id][q] == old_clust[r]){
                    u += 1;
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc8PAB[old_8PAB_id][p] == new_particle){
                    u += 1;
                }
        }
    if(u == new_clust_size){
        return 1;            
    }
    }

    for (int old_8PAA_id = 0; old_8PAA_id < n8PAA; ++old_8PAA_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<7; ++r){
                if(hc8PAA[old_8PAA_id][q] == old_clust[r]){
                    u += 1;
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc8PAA[old_8PAA_id][p] == new_particle){
                    u += 1;
                }
        }
    if(u == new_clust_size){
        return 1;            
    }
    }
    return 0;
}
