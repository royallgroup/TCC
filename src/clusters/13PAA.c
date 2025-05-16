#include <globals.h>
#include <tools.h>
#include "7T.h"
#include "7PAB.h"
#include "8PAA.h"
#include "12PAA.h"
#include "13PAA.h"
//!  A 13PAA cluster is made of a 12PAA cluster with one extra particle which forms a 3 ring
/*!
*  13PAA is made up of a 12PAA cluster and a sp3a cluster
*  The sp3a cluster shares two particles with the 13PAA cluster
*  One particle of the sp3a is a bonded spindle of the 12PAA, and one is a non- bonded spindle
*
*  Cluster output: BBBBBBBBBBBBB
*  Storage order: original 12PAA particles x 12, new sp3a particle)
*/

    void Clusters_Get13PAA() {
    for (int old_12PAA_id = 0; old_12PAA_id < n12PAA; ++old_12PAA_id) {
        int *old_12PAA_cluster = hc12PAA[old_12PAA_id];
        int clust_12PAA[12][2]; 
        clust_12PAA[0][0] = old_12PAA_cluster[3];
        clust_12PAA[1][0] = old_12PAA_cluster[2];
        clust_12PAA[2][0] = old_12PAA_cluster[5];
        clust_12PAA[3][0] = old_12PAA_cluster[11];
        clust_12PAA[4][0] = old_12PAA_cluster[0];
        clust_12PAA[5][0] = old_12PAA_cluster[4];
        clust_12PAA[6][0] = old_12PAA_cluster[1];
        clust_12PAA[7][0] = old_12PAA_cluster[6];
        clust_12PAA[8][0] = old_12PAA_cluster[7];
        clust_12PAA[9][0] = old_12PAA_cluster[8];
        clust_12PAA[10][0] = old_12PAA_cluster[9];
        clust_12PAA[11][0] = old_12PAA_cluster[10];
        
        clust_12PAA[0][1] = old_12PAA_cluster[9];
        clust_12PAA[1][1] = old_12PAA_cluster[1];
        clust_12PAA[2][1] = old_12PAA_cluster[6];
        clust_12PAA[3][1] = old_12PAA_cluster[10];
        clust_12PAA[4][1] = old_12PAA_cluster[8];
        clust_12PAA[5][1] = old_12PAA_cluster[7];
        clust_12PAA[6][1] = old_12PAA_cluster[2];
        clust_12PAA[7][1] = old_12PAA_cluster[5];
        clust_12PAA[8][1] = old_12PAA_cluster[4];
        clust_12PAA[9][1] = old_12PAA_cluster[0];
        clust_12PAA[10][1] = old_12PAA_cluster[3];
        clust_12PAA[11][1] = old_12PAA_cluster[11];

    for(int g = 0; g<2; ++g){
        int spindle_id_2 = clust_12PAA[5][g];
        int spindle_id = clust_12PAA[4][g];
        int nbonded_id = clust_12PAA[3][g];
        int bonded_id = clust_12PAA[0][g]; 
        int nbonded_id_2 = clust_12PAA[1][g];
        int bonded_id_2 = clust_12PAA[2][g]; 
        int final_part1 = clust_12PAA[6][g]; 
        int final_part2 = clust_12PAA[7][g]; 
        int final_part3 = clust_12PAA[8][g]; 
        int final_part4 = clust_12PAA[9][g]; 
        int final_part5 = clust_12PAA[10][g]; 
        int final_part6 = clust_12PAA[11][g]; 
        //printf("%i %i %i %i %i %i %i %i %i %i %i %i \n", bonded_id, nbonded_id_2, bonded_id_2, nbonded_id,spindle_id, 
        //spindle_id_2, final_part1, final_part2, final_part3, final_part4, final_part5, final_part6);
            for (int new_5A_pointer = 0; new_5A_pointer < nmem_sp3c[spindle_id_2]; ++new_5A_pointer) {
                int new_5A_id = mem_sp3c[spindle_id_2][new_5A_pointer];
                int *new_5A_cluster = hcsp3c[new_5A_id];
                if (is_particle_spindle_of_5A(spindle_id_2, new_5A_cluster) == 1) {
                    if(is_particle_ring_of_5A(bonded_id,new_5A_cluster) == 1){
                        if(is_particle_ring_of_5A(nbonded_id,new_5A_cluster) == 1){
                            if(is_particle_ring_of_5A(spindle_id,new_5A_cluster) == 1){
                                int new_particle_id = get_new_particle(new_5A_cluster, spindle_id_2);
                                    if(crossover_part(old_12PAA_cluster, new_particle_id, 13) == 0){
                                    int arr_6Z[6] = {bonded_id, nbonded_id_2, bonded_id_2, nbonded_id, spindle_id, spindle_id_2};
                                    int bond_counter = check_ring_bonds(new_5A_cluster, arr_6Z);
                                        if (bond_counter == 25 ){
                                            if (check_unique13PAA(old_12PAA_cluster, new_particle_id,13) == 0){
                                            add_13PAA(bonded_id, nbonded_id_2, bonded_id_2, nbonded_id,spindle_id, spindle_id_2, final_part1, final_part2, final_part3, final_part4, final_part5, final_part6, new_particle_id);
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


void add_13PAA(int bonded_id, int nbonded_id_2, int bonded_id_2, int nbonded_id, int spindle_id, int spindle_id_2, int final_part1, int final_part2, int final_part3, int final_part4, int final_part5, int final_part6, int new_particle) {
    int clusSize = 13;
    if (n13PAA == m13PAA) {
        hc13PAA = resize_2D_int(hc13PAA, m13PAA, m13PAA + incrStatic, clusSize, -1);
        m13PAA = m13PAA + incrStatic;
    }
    hc13PAA[n13PAA][0] = bonded_id;
    hc13PAA[n13PAA][1] = nbonded_id_2;
    hc13PAA[n13PAA][2] = bonded_id_2;
    hc13PAA[n13PAA][3] = nbonded_id;
    hc13PAA[n13PAA][4] = spindle_id;
    hc13PAA[n13PAA][5] = spindle_id_2;
    hc13PAA[n13PAA][6] = final_part1;
    hc13PAA[n13PAA][7] = final_part2;
    hc13PAA[n13PAA][8] = final_part3;
    hc13PAA[n13PAA][9] = final_part4;
    hc13PAA[n13PAA][10] = final_part5;
    hc13PAA[n13PAA][11] = final_part6;
    hc13PAA[n13PAA][12] = new_particle;

    for (int i = 0; i < 12; ++i) {
       s13PAA[hc13PAA[n13PAA][i]] = 'B';
    }
    ++n13PAA;
}

int check_unique13PAA(const int *old_clust, int new_particle, int new_clust_size){
    int u;
    for (int old_13PAA_id = 0; old_13PAA_id < n13PAA; ++old_13PAA_id) {
        u = 0;
        for (int q = 0; q<new_clust_size; ++q){
            for (int r = 0; r<12; ++r){
                if(hc13PAA[old_13PAA_id][q] == old_clust[r]){
                    u += 1;
                }
            }

        }
        for (int p = 0; p<new_clust_size; ++p){
            if(hc13PAA[old_13PAA_id][p] == new_particle){
                    u += 1;
                }
        }
    if(u == new_clust_size){
        return 1;            
    }
    }
    return 0;
}