#include <globals.h>
#include <tools.h>
#include "7PAB.h"
#include "7T.h"
#include "8PAA.h"
#include "9PAA.h"
#include "12PBB.h"
//!  A 11P is a 10P with an additional particle
/*!
* A 10P cluster contains two 6Z clusters at the top and bottom of the cluster. 
* A 11P takes these two 6Z clusters and tries to make a 7T from the cluster and an extra particle.
*  Cluster output: BBBBBBBBBB
*  Storage order: original 10P particles x 9, new 5A spindle)
*/

    void Clusters_Get12PBB() {
    for (int old_11PAB_id = 0; old_11PAB_id < n11PAB; ++old_11PAB_id) {
        int *old_11PAB_cluster = hc11PAB[old_11PAB_id];
        int clust_11PAB[11];
        clust_11PAB[0] = old_11PAB_cluster[8];
        clust_11PAB[1] = old_11PAB_cluster[2];
        clust_11PAB[2] = old_11PAB_cluster[1];
        clust_11PAB[3] = old_11PAB_cluster[9];
        clust_11PAB[4] = old_11PAB_cluster[7];
        clust_11PAB[5] = old_11PAB_cluster[6];
        clust_11PAB[6] = old_11PAB_cluster[5];
        clust_11PAB[7] = old_11PAB_cluster[4];
        clust_11PAB[8] = old_11PAB_cluster[0];
        clust_11PAB[9] = old_11PAB_cluster[3];
        clust_11PAB[10] = old_11PAB_cluster[10];


    for(int g = 0; g<2; ++g){
        int spindle_id_2 = clust_11PAB[5];
        int spindle_id = clust_11PAB[4];
        int nbonded_id = clust_11PAB[3];
        int bonded_id = clust_11PAB[0]; 
        int nbonded_id_2 = clust_11PAB[1];
        int bonded_id_2 = clust_11PAB[2]; 
        int final_part1 = clust_11PAB[6]; 
        int final_part2 = clust_11PAB[7]; 
        int final_part3 = clust_11PAB[8]; 
        int final_part4 = clust_11PAB[9]; 
        int final_part5 = clust_11PAB[10]; 

        int arr_6Z[6] = {bonded_id, nbonded_id_2, bonded_id_2, nbonded_id, spindle_id, spindle_id_2};
        for (int new_3_id = 0; new_3_id < nsp3a; ++new_3_id ){ //iterate over all sp3a rings
            int *new_3_cluster = hcsp3a[new_3_id];
            if(nring_in_cluster(old_11PAB_cluster, new_3_cluster, 12) == 2){
                if(nring_bonded_nbonded(bonded_id, nbonded_id, new_3_cluster) == 2){ //two particles in ring are in 6Z
                   int bond_counter = check_ring_bonds(new_3_cluster, arr_6Z);
                    if (bond_counter == 9){
                        int new_part = get_new_particle_P2(new_3_cluster, bonded_id, nbonded_id);
                        if (check_unique12PBB(old_11PAB_cluster, new_part, 12) == 0){
                            if(check_5A(arr_6Z, new_part) == 0){
                            add_12PBB(bonded_id, nbonded_id_2, bonded_id_2, nbonded_id,spindle_id, spindle_id_2, final_part1, final_part2, final_part3, final_part4, final_part5, new_part);
                            }
                            }
                        }
                    }
                }
            }
        }
    }
} 

void add_12PBB(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int final_part2, int final_part3, int final_part4, int final_part5, int new_particle) {
    int clusSize = 12;
    //printf("new_particle %i\n", new_particle);
    if (n12PBB == m12PBB) {
        hc12PBB = resize_2D_int(hc12PBB, m12PBB, m12PBB + incrStatic, clusSize, -1);
        m12PBB = m12PBB + incrStatic;
    }
    hc12PBB[n12PBB][0] = bonded_id;
    hc12PBB[n12PBB][1] = nbonded_id_2;
    hc12PBB[n12PBB][2] = bonded_id_2;
    hc12PBB[n12PBB][3] = nbonded_id;
    hc12PBB[n12PBB][4] = spindle_id;
    hc12PBB[n12PBB][5] = spindle_id_2;
    hc12PBB[n12PBB][6] = final_part1;
    hc12PBB[n12PBB][7] = final_part2;
    hc12PBB[n12PBB][8] = final_part3;
    hc12PBB[n12PBB][9] = final_part4;
    hc12PBB[n12PBB][10] = final_part5;
    hc12PBB[n12PBB][11] = new_particle;

    for (int i = 0; i < 12; ++i) {
        s12PBB[hc12PBB[n12PBB][i]] = 'B';
    }
    ++n12PBB;
}

int check_unique12PBB(const int *old_clust, int new_particle, int new_clust_size){
    int u;
    for (int old_12PBB_id = 0; old_12PBB_id < n12PBB; ++old_12PBB_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<11; ++r){
                //printf("%i %i\n", hc11P[old_11P_id][q],old_clust[r]);
                if(hc12PBB[old_12PBB_id][q] == old_clust[r]){
                    u += 1;
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc12PBB[old_12PBB_id][p] == new_particle){
                    u += 1;
                }
        }
    if(u == new_clust_size){
        return 1;            
    }
    }

    return 0;
}
