#include <globals.h>
#include <tools.h>
#include "9PBB.h"
#include "8PAB.h"
#include "7PAB.h"
#include "7T.h"
//!  A 9PBB is an 8PAB with an additional particle
/*!
* A 8PAB cluster contains two 6Z clusters at the top and bottom of the cluster. 
* A 9PBB takes these two 6Z clusters and tries to make a 7T from the cluster and an extra particle.
*  Cluster output: BBBBBBBBB
*  Storage order: original 8PAB particles x 8, new 5A spindle)
*/

    void Clusters_Get9PBB() {
    for (int old_8PAB_id = 0; old_8PAB_id < n8PAB; ++old_8PAB_id) {
        int *old_8PAB_cluster = hc8PAB[old_8PAB_id];
        int clust_6Z[8]; 
        int new_part, bond_counter;
        clust_6Z[0] = old_8PAB_cluster[1];
        clust_6Z[1] = old_8PAB_cluster[0];
        clust_6Z[2] = old_8PAB_cluster[4];
        clust_6Z[3] = old_8PAB_cluster[6];
        clust_6Z[4] = old_8PAB_cluster[2];
        clust_6Z[5] = old_8PAB_cluster[5];
        clust_6Z[6] = old_8PAB_cluster[3];
        clust_6Z[7] = old_8PAB_cluster[7];

        int spindle_id_2 = clust_6Z[5];
        int spindle_id = clust_6Z[4];
        int nbonded_id = clust_6Z[3];
        int bonded_id = clust_6Z[0]; 
        int nbonded_id_2 = clust_6Z[1];
        int bonded_id_2 = clust_6Z[2]; 
        int final_part1 = clust_6Z[6]; 
        int final_part2 = clust_6Z[7]; 
        int arr_6Z[6] = {bonded_id, nbonded_id_2, bonded_id_2, nbonded_id, spindle_id, spindle_id_2};
        for (int new_3_id = 0; new_3_id < nsp3a; ++new_3_id ){ //iterate over all sp3a rings
            int *new_3_cluster = hcsp3a[new_3_id];
            if(nring_in_cluster(old_8PAB_cluster, new_3_cluster, 9) == 2){
                if(nring_bonded_nbonded(bonded_id, nbonded_id, new_3_cluster) == 2){ //two particles in ring are in 6Z
                bond_counter = check_ring_bonds(new_3_cluster, arr_6Z);
                    if (bond_counter == 9){
                        new_part = get_new_particle_P2(new_3_cluster, bonded_id, nbonded_id);
                        if (check_unique9PBB(old_8PAB_cluster, new_part, 9) == 0){
                            if(check_5A(arr_6Z, new_part) == 0){
                            add_9PBB(bonded_id, nbonded_id_2, bonded_id_2, nbonded_id,spindle_id, spindle_id_2, final_part1, final_part2, new_part);
                            }
                            }
                        }
                    }
                }
            }
        } 
    }


void add_9PBB(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int final_part2, int new_particle) {
    int clusSize = 9;
    //printf("new_particle %i\n", new_particle);
    if (n9PBB == m9PBB) {
        hc9PBB = resize_2D_int(hc9PBB, m9PBB, m9PBB + incrStatic, clusSize, -1);
        m9PBB = m9PBB + incrStatic;
    }
    hc9PBB[n9PBB][0] = bonded_id;
    hc9PBB[n9PBB][1] = nbonded_id_2;
    hc9PBB[n9PBB][2] = bonded_id_2;
    hc9PBB[n9PBB][3] = nbonded_id;
    hc9PBB[n9PBB][4] = spindle_id;
    hc9PBB[n9PBB][5] = spindle_id_2;
    hc9PBB[n9PBB][6] = final_part1;
    hc9PBB[n9PBB][7] = final_part2;
    hc9PBB[n9PBB][8] = new_particle;

    for (int i = 0; i < 9; ++i) {
        s9PBB[hc9PBB[n9PBB][i]] = 'B';
    }
    ++n9PBB;

}

int check_unique9PBB(const int *old_clust, int new_particle, int new_clust_size){
    int u;
    for (int old_9PBB_id = 0; old_9PBB_id < n9PBB; ++old_9PBB_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<8; ++r){
                //printf("%i %i\n", hc11P[old_11P_id][q],old_clust[r]);
                if(hc9PBB[old_9PBB_id][q] == old_clust[r]){
                    u += 1;
             
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc9PBB[old_9PBB_id][p] == new_particle){
                    u += 1;
                }
        }

    if(u == new_clust_size){
        return 1;            
    }
    }
    return 0;
}
