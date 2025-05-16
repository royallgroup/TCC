
#include <globals.h>
#include <tools.h>
#include "6Z.h"
#include "7PAB.h"
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

void Clusters_Get7PAB() {
for (int old_6Z_id = 0; old_6Z_id < n6Z; ++old_6Z_id) {
    int *old_6Z_cluster = hc6Z[old_6Z_id];
    int clust_6Z[6][2]; 
    clust_6Z[0][0] = old_6Z_cluster[0];
    clust_6Z[1][0] = old_6Z_cluster[1];
    clust_6Z[2][0] = old_6Z_cluster[2];
    clust_6Z[3][0] = old_6Z_cluster[3];
    clust_6Z[4][0] = old_6Z_cluster[4];
    clust_6Z[5][0] = old_6Z_cluster[5];
    
    clust_6Z[0][1] = old_6Z_cluster[2];
    clust_6Z[1][1] = old_6Z_cluster[3];
    clust_6Z[2][1] = old_6Z_cluster[0];
    clust_6Z[3][1] = old_6Z_cluster[1];
    clust_6Z[4][1] = old_6Z_cluster[5];
    clust_6Z[5][1] = old_6Z_cluster[4];
    for(int g = 0; g<2; ++g){
        int spindle_id_2 = clust_6Z[5][g];
        int spindle_id = clust_6Z[4][g];
        int nbonded_id = clust_6Z[3][g];
        int bonded_id = clust_6Z[0][g]; 
        int nbonded_id_2 = clust_6Z[1][g];
        int bonded_id_2 = clust_6Z[2][g]; 
        int arr_6Z[6] = {bonded_id, nbonded_id_2, bonded_id_2, nbonded_id, spindle_id, spindle_id_2};
        for (int new_3_id = 0; new_3_id < nsp3a; ++new_3_id ){ //iterate over all sp3a rings
            int *new_3_cluster = hcsp3a[new_3_id];
            if(nring_in_cluster(old_6Z_cluster, new_3_cluster, 6) == 2){
                if(nring_bonded_nbonded(bonded_id, nbonded_id, new_3_cluster) == 2){ //two particles in ring are in 6Z
                   int bond_counter = check_ring_bonds(new_3_cluster, arr_6Z);
                    if (bond_counter == 9){
                        int new_part = get_new_particle_P2(new_3_cluster, bonded_id, nbonded_id);
                        int c = check_unique7PAB(old_6Z_cluster, new_part, 7);
                        if (check_unique7PAB(old_6Z_cluster, new_part, 7) == 0){
                            if(check_5A(arr_6Z, new_part) == 0){
                            add_7PAB(bonded_id, nbonded_id_2, bonded_id_2, nbonded_id,spindle_id, spindle_id_2, new_part);
                            }
                            }
                        }
                    }
                }
            }
        }
    }
}

int check_5A(const int *arr_6Z, int new_part){
    int i;
    for (int new_5A_pointer = 0; new_5A_pointer < nsp3c; ++new_5A_pointer) {
    i = 0;
  
    int *new_5A_cluster = hcsp3c[new_5A_pointer];
    for(int j = 0; j < 5; ++j){
        for(int k = 0; k < 6; ++k){
            if(new_5A_cluster[j] == arr_6Z[k]){
                i += 1;
                }
            }
        if(new_5A_cluster[j] == new_part){
            i += 2;
        }
        if(i == 6){
            return 1;
        }
    }
    }
    return 0;
}

int nring_bonded_nbonded(int bonded_id, int nbonded_id, const int *clust_2){
    int cross_count =0; 
    if(bonded_id == clust_2[0]){
        cross_count += 1;
    }
    if(bonded_id == clust_2[1]){
        cross_count += 1;
    }
    if(bonded_id == clust_2[2]){
        cross_count += 1;
    }
    
    if(nbonded_id == clust_2[0]){
        cross_count += 1;
    }
    if(nbonded_id == clust_2[1]){
        cross_count += 1;
    }
    if(nbonded_id == clust_2[2]){
        cross_count += 1;
    }
    return cross_count;
}

int nring_in_cluster(const int *clust_1, const int *clust_2, int clust_1_size){
    int s1 = clust_1_size;
    int cross_count =0; 
    for(int n=0; n<s1; ++n){
        for(int j = 0; j<3; ++j){
            if(clust_1[n]==clust_2[j]){
                cross_count += 1;
        }

        }
    }
    return cross_count;
}

void add_7PAB(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int new_particle) {
    int clusSize = 7;
    if (n7PAB == m7PAB) {
        hc7PAB = resize_2D_int(hc7PAB, m7PAB, m7PAB + incrStatic, clusSize, -1);
        m7PAB = m7PAB + incrStatic;
    }
    hc7PAB[n7PAB][0] = bonded_id;
    hc7PAB[n7PAB][1] = nbonded_id_2;
    hc7PAB[n7PAB][2] = bonded_id_2;
    hc7PAB[n7PAB][3] = nbonded_id;
    hc7PAB[n7PAB][4] = spindle_id;
    hc7PAB[n7PAB][5] = spindle_id_2;
    hc7PAB[n7PAB][6] = new_particle;

    for (int i = 0; i < 7; ++i) {
        s7PAB[hc7PAB[n7PAB][i]] = 'B';
    }
    ++n7PAB;
}

int get_new_particle_P2(const int *new_3_ring, int bonded_id, int nbonded_id) {
    for(int ring_id = 0; ring_id <3; ++ring_id){
        int noverlap = 0;
        if(new_3_ring[ring_id] != bonded_id && new_3_ring[ring_id] != nbonded_id){
            return new_3_ring[ring_id];
        }
    }
}

int check_unique7PAB(const int *old_clust, int new_particle, int new_clust_size){
    int u;
    for (int old_7PAB_id = 0; old_7PAB_id < n7PAB; ++old_7PAB_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<6; ++r){
                if(hc7PAB[old_7PAB_id][q] == old_clust[r]){
                    u += 1;
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc7PAB[old_7PAB_id][p] == new_particle){
                    u += 1;
                }
        }
    if(u == new_clust_size){
        return 1;            
    }
    }
    return 0;
}
