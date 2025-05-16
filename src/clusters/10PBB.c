#include <globals.h>
#include <tools.h>
#include "10PBB.h"
#include "9PAB.h"
#include "7PAB.h"
#include "7T.h"
//!  A 10PBB is an 9PAB with an additional particle
/*!
* A 9PAB cluster contains two 6Z clusters at the top and bottom of the cluster. 
* A 10PBB takes these two 6Z clusters and tries to make a 7T from the cluster and an extra particle.
*  Cluster output: BBBBBBBBB
*  Storage order: original 9PAB particles x 8, new 5A spindle)
*/

    void Clusters_Get10PBB() {
    for (int old_9PAB_id = 0; old_9PAB_id < n9PAB; ++old_9PAB_id) {
        int *old_9PAB_cluster = hc9PAB[old_9PAB_id];
        int clust_6Z[8]; 
        
        clust_6Z[0] = old_9PAB_cluster[1];
        clust_6Z[1] = old_9PAB_cluster[0];
        clust_6Z[2] = old_9PAB_cluster[4];
        clust_6Z[3] = old_9PAB_cluster[6];
        clust_6Z[4] = old_9PAB_cluster[2];
        clust_6Z[5] = old_9PAB_cluster[5];
        clust_6Z[6] = old_9PAB_cluster[3];
        clust_6Z[7] = old_9PAB_cluster[7];

        int spindle_id_2 = clust_6Z[5];
        int spindle_id = clust_6Z[4];
        int nbonded_id = clust_6Z[3];
        int bonded_id = clust_6Z[0]; 
        int nbonded_id_2 = clust_6Z[1];
        int bonded_id_2 = clust_6Z[2]; 
        int final_part1 = clust_6Z[6]; 
        int final_part2 = clust_6Z[7]; 
        int final_part3 = clust_6Z[8];
        int arr_6Z[6] = {bonded_id, nbonded_id_2, bonded_id_2, nbonded_id, spindle_id, spindle_id_2};
        for (int new_3_id = 0; new_3_id < nsp3a; ++new_3_id ){ //iterate over all sp3a rings
            int *new_3_cluster = hcsp3a[new_3_id];
            if(nring_in_cluster(old_9PAB_cluster, new_3_cluster, 10) == 2){
                if(nring_bonded_nbonded(bonded_id, nbonded_id, new_3_cluster) == 2){ //two particles in ring are in 6Z
                int bond_counter = check_ring_bonds(new_3_cluster, arr_6Z);
                    if (bond_counter == 9){
                        int new_part = get_new_particle_P2(new_3_cluster, bonded_id, nbonded_id);
                        if (check_unique10PBB(old_9PAB_cluster, new_part, 10) == 0){
                            if(check_5A(arr_6Z, new_part) == 0){
                            add_10PBB(bonded_id, nbonded_id_2, bonded_id_2, nbonded_id,spindle_id, spindle_id_2, final_part1, final_part2, final_part3, new_part);
                            }
                            }
                        }
                    }
                }
            }
        } 
    }


void add_10PBB(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int final_part2, int final_part3, int new_particle) {
    int clusSize = 9;
    //printf("new_particle %i\n", new_particle);
    if (n10PBB == m10PBB) {
        hc10PBB = resize_2D_int(hc10PBB, m10PBB, m10PBB + incrStatic, clusSize, -1);
        m10PBB = m10PBB + incrStatic;
    }
    hc10PBB[n10PBB][0] = bonded_id;
    hc10PBB[n10PBB][1] = nbonded_id_2;
    hc10PBB[n10PBB][2] = bonded_id_2;
    hc10PBB[n10PBB][3] = nbonded_id;
    hc10PBB[n10PBB][4] = spindle_id;
    hc10PBB[n10PBB][5] = spindle_id_2;
    hc10PBB[n10PBB][6] = final_part1;
    hc10PBB[n10PBB][7] = final_part2;
    hc10PBB[n10PBB][8] = final_part2;
    hc10PBB[n10PBB][9] = new_particle;

    for (int i = 0; i < 10; ++i) {
        s10PBB[hc10PBB[n10PBB][i]] = 'B';
    }
    ++n10PBB;
}

int check_unique10PBB(const int *old_clust, int new_particle, int new_clust_size){
    int u;
    for (int old_10PBB_id = 0; old_10PBB_id < n10PBB; ++old_10PBB_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<9; ++r){
                //printf("%i %i\n", hc11P[old_11P_id][q],old_clust[r]);
                if(hc10PBB[old_10PBB_id][q] == old_clust[r]){
                    u += 1;
             
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc10PBB[old_10PBB_id][p] == new_particle){
                    u += 1;
                }
        }

    if(u == new_clust_size){
        return 1;            
    }
    }
    return 0;
}
