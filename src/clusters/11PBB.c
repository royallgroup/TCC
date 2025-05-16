#include <globals.h>
#include <tools.h>
#include "7PAB.h"
#include "7T.h"
#include "8PAA.h"
#include "11PBB.h"
//!  A 11PBB is a 10PAB with an additional particle
/*!
* 11PBB is a 10PAB cluster with an sp3a at one end. 
* Two of the sp3a particles form the bonded and nonbonded pair of one 
* of the 6Zs at either end of the 6Z.
*The other sp3a particle is the additional particle
*  Storage order: bonded, nbonded2, bonded2, nbonded, spindle1, spindle 2,
* 3 extra particles top to bottom, new particle)
*/

    void Clusters_Get11PBB() {
    for (int old_10PAB_id = 0; old_10PAB_id < n10PAB; ++old_10PAB_id) {
        int *old_10PAB_cluster = hc10PAB[old_10PAB_id];
        int clust_10PAB[10]; 
        clust_10PAB[0] = old_10PAB_cluster[7];
        clust_10PAB[1] = old_10PAB_cluster[5];
        clust_10PAB[2] = old_10PAB_cluster[2];
        clust_10PAB[3] = old_10PAB_cluster[8];
        clust_10PAB[4] = old_10PAB_cluster[6];
        clust_10PAB[5] = old_10PAB_cluster[1];
        clust_10PAB[6] = old_10PAB_cluster[4];
        clust_10PAB[7] = old_10PAB_cluster[0];
        clust_10PAB[8] = old_10PAB_cluster[3];
        clust_10PAB[9] = old_10PAB_cluster[9];
        int spindle_id_2 = clust_10PAB[5];
        int spindle_id = clust_10PAB[4];
        int nbonded_id = clust_10PAB[3];
        int bonded_id = clust_10PAB[0]; 
        int nbonded_id_2 = clust_10PAB[1];
        int bonded_id_2 = clust_10PAB[2]; 
        int final_part1 = clust_10PAB[6]; 
        int final_part2 = clust_10PAB[7]; 
        int final_part3 = clust_10PAB[8]; 
        int final_part4 = clust_10PAB[9]; 

        int arr_6Z[6] = {bonded_id, nbonded_id_2, bonded_id_2, nbonded_id, spindle_id, spindle_id_2};
        for (int new_3_id = 0; new_3_id < nsp3a; ++new_3_id ){ //iterate over all sp3a rings
            int *new_3_cluster = hcsp3a[new_3_id];
            if(nring_in_cluster(old_10PAB_cluster, new_3_cluster, 11) == 2){
                if(nring_bonded_nbonded(bonded_id, nbonded_id, new_3_cluster) == 2){ //two particles in ring are in 6Z
                   int bond_counter = check_ring_bonds(new_3_cluster, arr_6Z);
                    if (bond_counter == 9){
                        int new_part = get_new_particle_P2(new_3_cluster, bonded_id, nbonded_id);
                        if (check_unique11PBB(old_10PAB_cluster, new_part, 11) == 0){
                            if(check_5A(arr_6Z, new_part) == 0){
                            add_11PBB(bonded_id, nbonded_id_2, bonded_id_2, nbonded_id,spindle_id, spindle_id_2, final_part1, final_part2, final_part3, final_part4, new_part);
                            }
                            }
                        }
                    }
                }
            }
        } 
    } 

void add_11PBB(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int final_part2, int final_part3, int final_part4, int new_particle) {
    int clusSize = 11;
    //printf("new_particle %i\n", new_particle);
    if (n11PBB == m11PBB) {
        hc11PBB = resize_2D_int(hc11PBB, m11PBB, m11PBB + incrStatic, clusSize, -1);
        m11PBB = m11PBB + incrStatic;
    }
    hc11PBB[n11PBB][0] = bonded_id;
    hc11PBB[n11PBB][1] = nbonded_id_2;
    hc11PBB[n11PBB][2] = bonded_id_2;
    hc11PBB[n11PBB][3] = nbonded_id;
    hc11PBB[n11PBB][4] = spindle_id;
    hc11PBB[n11PBB][5] = spindle_id_2;
    hc11PBB[n11PBB][6] = final_part1;
    hc11PBB[n11PBB][7] = final_part2;
    hc11PBB[n11PBB][8] = final_part3;
    hc11PBB[n11PBB][9] = final_part4;
    hc11PBB[n11PBB][10] = new_particle;

    for (int i = 0; i < 11; ++i) {
        s11PBB[hc11PBB[n11PBB][i]] = 'B';
    }
    ++n11PBB;
}

int check_unique11PBB(const int *old_clust, int new_particle, int new_clust_size){
    int u;
    for (int old_11PBB_id = 0; old_11PBB_id < n11PBB; ++old_11PBB_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<10; ++r){
                //printf("%i %i\n", hc11P[old_11P_id][q],old_clust[r]);
                if(hc11PBB[old_11PBB_id][q] == old_clust[r]){
                    u += 1;
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc11PBB[old_11PBB_id][p] == new_particle){
                    u += 1;
                }
        }
    if(u == new_clust_size){
        return 1;            
    }
    }
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
    return 0;
}
