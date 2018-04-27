#include <clusters/simple_cluster_methods.h>
#include "globals.h"
#include "bonds.h"
#include "tools.h"
#include "11C.h"

void Clusters_Get11C() {
    //!  An 11C cluster is the intersection of two 7A clusters with a common spindle
    /*!
   *  Find 11C clusters
   *  An 11C is two 7A clusters where:
   *      - There are two common particles between the two sp5 rings. These are a bonded pair.
   *      - There are two more bonds between two pairs of distinct particles in the sp5 ring.
   *
   *  Cluster output: SOOBBBBBBBB
   *  Storage order: common_spindle x 1, uncommon_spindles x 2, common_ring_particles x 2,
   */

    int common_ring_particles[5],uncommon_spindle[2];
    int k, l, m, common_spindle;
    int break_out;

    int second_7A_pointer;

    common_ring_particles[0]=common_ring_particles[1]=uncommon_spindle[0]=uncommon_spindle[1]=common_spindle=-1;

    for (int first_7A_id = 0; first_7A_id < nsp5c - 1; ++first_7A_id) {
        int *first_7A_cluster = hcsp5c[first_7A_id];
        for(int first_spindle_pointer = 5; first_spindle_pointer < 7; first_spindle_pointer++) {
            int first_spindle_id = first_7A_cluster[first_spindle_pointer];
            for (second_7A_pointer = 0; second_7A_pointer < nmem_sp5c[first_spindle_id]; second_7A_pointer++) {
                int second_7A_id = mem_sp5c[first_spindle_id][second_7A_pointer];
                if(second_7A_id > first_7A_id) {
                    int *second_7A_cluster = hcsp5c[second_7A_id];

                    if (get_11C_spindle_particles(uncommon_spindle, first_7A_id, second_7A_id, &common_spindle) == 0) continue;

                    // need two common particles from SP5 rings
                    if (count_common_ring_particles(first_7A_cluster, second_7A_cluster, 5, 5, common_ring_particles) == 2) {

                        // the common ring particles need to be bonded
                        if (Bonds_BondCheck(common_ring_particles[0], common_ring_particles[1]) == 1) {

                            // two bonds between non-common SP5 ring particles
                            if (count_bonded_ring_particles_11C(common_ring_particles, first_7A_cluster, second_7A_cluster) == 2) {


                                int clusSize = 11;

                                if (n11C == m11C) {
                                    hc11C = resize_2D_int(hc11C, m11C, m11C + incrStatic, clusSize, -1);
                                    m11C = m11C + incrStatic;
                                }

                                // hc11C key: (s_com, s_i, s_j, r_ca, r_cb, d_i, d_i, d_j, d_j, unc_i, unc_j)
                                hc11C[n11C][0] = common_spindle;
                                hc11C[n11C][1] = uncommon_spindle[0];
                                hc11C[n11C][2] = uncommon_spindle[1];
                                hc11C[n11C][3] = common_ring_particles[0];
                                hc11C[n11C][4] = common_ring_particles[1];

                                l = 5;
                                m = 7;
                                break_out = 0;
                                for (k = 0; k < 5; ++k) {
                                    if (Bonds_BondCheck(first_7A_cluster[k], common_ring_particles[0]) &&
                                        first_7A_cluster[k] != common_ring_particles[1]) {
                                        if (l == 7) {
                                            break_out = 1;
                                            break;
                                        }
                                        hc11C[n11C][l] = first_7A_cluster[k];
                                        l++;
                                    }
                                }
                                for (k = 0; k < 5; ++k) {
                                    if (Bonds_BondCheck(first_7A_cluster[k], common_ring_particles[1]) &&
                                        first_7A_cluster[k] != common_ring_particles[0]) {
                                        if (l == 7) {
                                            break_out = 1;
                                            break;
                                        }
                                        hc11C[n11C][l] = first_7A_cluster[k];
                                        l++;
                                    }
                                }
                                for (k = 0; k < 5; ++k) {
                                    if (Bonds_BondCheck(second_7A_cluster[k], common_ring_particles[0]) &&
                                        second_7A_cluster[k] != common_ring_particles[1]) {
                                        if (m == 9) {
                                            break_out = 1;
                                            break;
                                        }
                                        hc11C[n11C][m] = second_7A_cluster[k];
                                        m++;
                                    }
                                }
                                for (k = 0; k < 5; ++k) {
                                    if (Bonds_BondCheck(second_7A_cluster[k], common_ring_particles[1]) &&
                                        second_7A_cluster[k] != common_ring_particles[0]) {
                                        if (m == 9) {
                                            break_out = 1;
                                            break;
                                        }
                                        hc11C[n11C][m] = second_7A_cluster[k];
                                        m++;
                                    }
                                }
                                if (break_out == 1 || l < 7 || m < 9) continue;

                                // Check that the bonded non-common particles are bonded
                                if (Bonds_BondCheck(hc11C[n11C][5], hc11C[n11C][7]) == 0) continue;
                                if (Bonds_BondCheck(hc11C[n11C][6], hc11C[n11C][8]) == 0) continue;

                                // Get the ID's of the non-common particles
                                for (k = 0; k < 5; ++k) {
                                    if (Bonds_BondCheck(first_7A_cluster[k], hc11C[n11C][5]) &&
                                        Bonds_BondCheck(first_7A_cluster[k], hc11C[n11C][6])) {
                                        hc11C[n11C][9] = first_7A_cluster[k];
                                    }
                                    if (Bonds_BondCheck(second_7A_cluster[k], hc11C[n11C][7]) &&
                                        Bonds_BondCheck(second_7A_cluster[k], hc11C[n11C][8])) {
                                        hc11C[n11C][10] = second_7A_cluster[k];
                                    }
                                }
                                quickSort(&hc11C[n11C][1], 2);
                                quickSort(&hc11C[n11C][3], 2);
                                quickSort(&hc11C[n11C][5], 4);
                                quickSort(&hc11C[n11C][9], 2);

                                Cluster_Write_11C();

                                ++n11C;
                            }
                        }
                    }
                }
            }
        }
    }
}

int count_bonded_ring_particles_11C(const int *common_ring_particles, const int *first_7A_cluster, const int *second_7A_cluster) {
    int num_bonded_ring_particles = 0;
    for (int first_ring_pointer = 0; first_ring_pointer < 5; ++first_ring_pointer) {
        int first_particle_id = first_7A_cluster[first_ring_pointer];

        if (is_particle_in_cluster(common_ring_particles, 2, first_particle_id) == 0) {
            for (int second_ring_pointer = 0; second_ring_pointer < 5; ++second_ring_pointer) {
                int second_particle_id = second_7A_cluster[second_ring_pointer];

                if (is_particle_in_cluster(common_ring_particles, 2, second_particle_id) == 0) {
                    if (Bonds_BondCheck(first_particle_id, second_particle_id)) {
                        ++num_bonded_ring_particles;
                    }
                }
            }
        }
    }
    return num_bonded_ring_particles;
}

int get_11C_spindle_particles(int *uncommon_spindle, int id_first_7A, int id_second7A, int *common_spindle) {
    int num_common_spindles = 0;

    if (hcsp5c[id_first_7A][5] == hcsp5c[id_second7A][5]) {
        (*common_spindle) = hcsp5c[id_first_7A][5];
        uncommon_spindle[0] = hcsp5c[id_first_7A][6];
        uncommon_spindle[1] = hcsp5c[id_second7A][6];
        ++num_common_spindles;
    }
    if (hcsp5c[id_first_7A][6] == hcsp5c[id_second7A][6]) {
        (*common_spindle) = hcsp5c[id_first_7A][6];
        uncommon_spindle[0] = hcsp5c[id_first_7A][5];
        uncommon_spindle[1] = hcsp5c[id_second7A][5];
        ++num_common_spindles;
    }
    if (hcsp5c[id_first_7A][5] == hcsp5c[id_second7A][6]) {
        (*common_spindle) = hcsp5c[id_first_7A][5];
        uncommon_spindle[0] = hcsp5c[id_first_7A][6];
        uncommon_spindle[1] = hcsp5c[id_second7A][5];
        ++num_common_spindles;
    }
    if (hcsp5c[id_first_7A][6] == hcsp5c[id_second7A][5]) {
        (*common_spindle) = hcsp5c[id_first_7A][6];
        uncommon_spindle[0] = hcsp5c[id_first_7A][5];
        uncommon_spindle[1] = hcsp5c[id_second7A][6];
        ++num_common_spindles;
    }

    if (num_common_spindles == 1) return 1;
    else return 0;
}

void Cluster_Write_11C() {
    int i;

    s11C[hc11C[n11C][0]] = 'S';
    if(s11C[hc11C[n11C][1]] != 'S') s11C[hc11C[n11C][1]] = 'O';
    if(s11C[hc11C[n11C][2]] != 'S') s11C[hc11C[n11C][2]] = 'O';
    for(i=3; i< 11; i++) {
        if (s11C[hc11C[n11C][i]] == 'C') s11C[hc11C[n11C][i]] = 'B';
    }
}