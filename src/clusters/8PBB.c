
#include <globals.h>
#include <tools.h>
#include "6Z.h"
#include "7PAB.h"
#include "8PBB.h"
#include "7T.h"
//!  A 7PAB is a 6Z with an additional particle
/*!
* 7PAB is a 6Z cluster with an sp3a at one end. 
* Two of the sp3a particles form the bonded and nonbonded pair of one 
* of the 6Zs at either end of the 6Z.
*The other sp3a particle is the additional particle
*  Storage order: bonded, nbonded2, bonded2, nbonded, spindle1, spindle 2,
* new particle)
*/

void Clusters_Get8PBB() {
for (int old_7PAB_id = 0; old_7PAB_id < n7PAB; ++old_7PAB_id) {
    int *old_7PAB_cluster = hc7PAB[old_7PAB_id];
    int clust_7PAB[7]; 
    clust_7PAB[0] = old_7PAB_cluster[2];
    clust_7PAB[1] = old_7PAB_cluster[3];
    clust_7PAB[2] = old_7PAB_cluster[0];
    clust_7PAB[3] = old_7PAB_cluster[1];
    clust_7PAB[4] = old_7PAB_cluster[5];
    clust_7PAB[5] = old_7PAB_cluster[4];
    clust_7PAB[6] = old_7PAB_cluster[6];

    int spindle_id_2 = clust_7PAB[5];
    int spindle_id = clust_7PAB[4];
    int nbonded_id = clust_7PAB[3];
    int bonded_id = clust_7PAB[0]; 
    int nbonded_id_2 = clust_7PAB[1];
    int bonded_id_2 = clust_7PAB[2]; 
    int final_part1 = clust_7PAB[6];
    int arr_6Z[6] = {bonded_id, nbonded_id_2, bonded_id_2, nbonded_id, spindle_id, spindle_id_2};
    for (int new_3_id = 0; new_3_id < nsp3a; ++new_3_id ){ //iterate over all sp3a rings
        int *new_3_cluster = hcsp3a[new_3_id];
        if(nring_in_cluster(old_7PAB_cluster, new_3_cluster, 8) == 2){
            if(nring_bonded_nbonded(bonded_id, nbonded_id, new_3_cluster) == 2){ //two particles in ring are in 6Z
                int bond_counter = check_ring_bonds(new_3_cluster, arr_6Z);
                if (bond_counter == 9){
                    int new_part = get_new_particle_P2(new_3_cluster, bonded_id, nbonded_id);
                    if (check_unique8PBB(old_7PAB_cluster, new_part, 8) == 0){
                        if(check_5A(arr_6Z, new_part) == 0){
                        add_8PBB(bonded_id, nbonded_id_2, bonded_id_2, nbonded_id,spindle_id, spindle_id_2, final_part1, new_part);
                        }
                        }
                    }
                }
            }
        }
    }
}

void add_8PBB(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int new_particle) {
    int clusSize = 8;
    if (n8PBB == m8PBB) {
        hc8PBB = resize_2D_int(hc8PBB, m8PBB, m8PBB + incrStatic, clusSize, -1);
        m8PBB = m8PBB + incrStatic;
    }
    hc8PBB[n8PBB][0] = bonded_id;
    hc8PBB[n8PBB][1] = nbonded_id_2;
    hc8PBB[n8PBB][2] = bonded_id_2;
    hc8PBB[n8PBB][3] = nbonded_id;
    hc8PBB[n8PBB][4] = spindle_id;
    hc8PBB[n8PBB][5] = spindle_id_2;
    hc8PBB[n8PBB][7] = new_particle;
    hc8PBB[n8PBB][6] = final_part1;

    for (int i = 0; i < 8; ++i) {
        s8PBB[hc8PBB[n8PBB][i]] = 'B';
    }
    ++n8PBB;
}

int check_unique8PBB(const int *old_clust, int new_particle, int new_clust_size){
    int u;
    for (int old_8PBB_id = 0; old_8PBB_id < n8PBB; ++old_8PBB_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<7; ++r){
                if(hc8PBB[old_8PBB_id][q] == old_clust[r]){
                    u += 1;
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc8PBB[old_8PBB_id][p] == new_particle){
                    u += 1;
                }
        }
    if(u == new_clust_size){
        return 1;            
    }
    }
    return 0;
}
