#include <globals.h>
#include <tools.h>
#include "9PAA.h"
#include "7T.h"
#include "8PAA.h"
#include "10PAA.h"
//!  A 10PAA is a 9PAA with an additional particle
/*!
* A 9PAA cluster contains two 6Z clusters at the top and bottom of the cluster. 
* A 10PAA takes these two 6Z clusters and tries to make a 7T from the cluster and an extra particle.
*  Cluster output: BBBBBBBBBB
*  Storage order: original 9PAA particles x 9, new 5A spindle)
*/

    void Clusters_Get10PAA() {
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
            for (int new_5A_pointer = 0; new_5A_pointer < nmem_sp3c[spindle_id_2]; ++new_5A_pointer) {
                int new_5A_id = mem_sp3c[spindle_id_2][new_5A_pointer];
                int *new_5A_cluster = hcsp3c[new_5A_id];
                if (is_particle_spindle_of_5A(spindle_id_2, new_5A_cluster) == 1) {
                    if(is_particle_ring_of_5A(bonded_id,new_5A_cluster) == 1){
                        if(is_particle_ring_of_5A(nbonded_id,new_5A_cluster) == 1){
                            if(is_particle_ring_of_5A(spindle_id,new_5A_cluster) == 1){
                                int new_particle_id = get_new_particle(new_5A_cluster, spindle_id_2);
                                    if(crossover_part(old_9PAA_cluster, new_particle_id, 10) == 0){
                                    int arr_6Z[6] = {bonded_id, nbonded_id_2, bonded_id_2, nbonded_id, spindle_id, spindle_id_2};
                                    int bond_counter = check_ring_bonds(new_5A_cluster, arr_6Z);
                                        if (bond_counter == 25 ){
                                            if (check_unique10PAA(old_9PAA_cluster, new_particle_id,10) == 0){
                                            add_10PAA(bonded_id, nbonded_id_2, bonded_id_2, nbonded_id,spindle_id, spindle_id_2, final_part1, final_part2, final_part3, new_particle_id);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    } 
} 

void add_10PAA(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int final_part2, int final_part3, int new_particle) {
    int clusSize = 10;
    //printf("new_particle %i\n", new_particle);
    if (n10PAA == m10PAA) {
        hc10PAA = resize_2D_int(hc10PAA, m10PAA, m10PAA + incrStatic, clusSize, -1);
        m10PAA = m10PAA + incrStatic;
    }
    hc10PAA[n10PAA][0] = bonded_id;
    hc10PAA[n10PAA][1] = nbonded_id_2;
    hc10PAA[n10PAA][2] = bonded_id_2;
    hc10PAA[n10PAA][3] = nbonded_id;
    hc10PAA[n10PAA][4] = spindle_id;
    hc10PAA[n10PAA][5] = spindle_id_2;
    hc10PAA[n10PAA][6] = final_part1;
    hc10PAA[n10PAA][7] = final_part2;
    hc10PAA[n10PAA][8] = final_part3;
    hc10PAA[n10PAA][9] = new_particle;

    for (int i = 0; i < 10; ++i) {
        s10PAA[hc10PAA[n10PAA][i]] = 'B';
    }
    ++n10PAA;
}

int check_unique10PAA(const int *old_clust, int new_particle, int new_clust_size){
    int u;
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