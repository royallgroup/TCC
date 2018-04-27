#include <globals.h>
#include <tools.h>
#include "9K.h"

void Clusters_Get9K() {

    //!  An 9K cluster is the intersection of two 6A clusters sharing a spindle and two ring particles.
    /*!
   *  Find 9K clusters
   *  An 9K is 2 6A clusters where:
   *      - There is one common spindle.
   *      - The uncommon spindles are not in the opposite cluster.
   *      - There are two common particles between the two 6A rings
   *      - The uncommon sp4 ring particle are bonded
   *
   *  Cluster output: BBBBBBOOO
   *  Storage order: common_ring_particles x 2, uncommon_ring_particles x 4, uncommon_spindles x 2, common_spindle
   */

    int j, k, l, m;
    int cp[2], scom, sother[2];
    int trial[9];


    cp[0]=cp[1]=scom=sother[0]=sother[1]=-1;


    for(int first_6A_id = 0; first_6A_id < nsp4c - 1; ++first_6A_id) {
        int *first_6A_cluster = hcsp4c[first_6A_id];
        for (int spindle_pointer = 4; spindle_pointer < 6; spindle_pointer++) {
            for (j = 0; j < nmem_sp4c[first_6A_cluster[spindle_pointer]]; ++j) {
                int second_6A_id = mem_sp4c[first_6A_cluster[spindle_pointer]][j];
                int *second_6A_cluster = hcsp4c[second_6A_id];
                if (second_6A_id > first_6A_id) {

                    m = 0;        // sp4c_i and sp4c_mem_sp4c[sp4c[i][j2]][j] have exactly one common spindle
                    for (k = 4; k < 6; ++k) {
                        for (l = 4; l < 6; ++l) {
                            if (first_6A_cluster[k] == second_6A_cluster[l]) {
                                m++;
                                scom = first_6A_cluster[k];
                            }
                        }
                    }
                    if (m != 1) continue;

                    if (first_6A_cluster[4] == scom) sother[0] = first_6A_cluster[5];
                    else sother[0] = first_6A_cluster[4];
                    if (second_6A_cluster[4] == scom) sother[1] = second_6A_cluster[5];
                    else sother[1] = second_6A_cluster[4];

                    m = 0;        // check sother[0] is not in cluster sp4c_mem_sp4c[sp4c[i][j2]][j]
                    for (k = 0; k < 6; ++k) {
                        if (sother[0] == second_6A_cluster[k] || sother[1] == first_6A_cluster[k]) {
                            m++;
                        }
                    }
                    if (m != 0) continue;

                    m = 0;        // SP4 ring from sp4c_i and SP4 ring from sp4c_mem_sp4c[sp4c[i][j2]][j] have exactly two common particles
                    for (k = 0; k < 4; ++k) {
                        for (l = 0; l < 4; ++l) {
                            if (first_6A_cluster[k] == second_6A_cluster[l]) {
                                if (m >= 2) {
                                    m = 3;
                                    break;
                                }
                                cp[m] = first_6A_cluster[k];
                                m++;
                            }
                        }
                        if (m > 2) break;
                    }
                    if (m != 2) continue;

                    trial[0] = cp[0];
                    trial[1] = cp[1];
                    trial[6] = sother[0];
                    trial[7] = sother[1];
                    trial[8] = scom;
                    quickSort(&trial[0], 2);
                    quickSort(&trial[6], 2);

                    m = 2;        // find uncommon SP4 ring particles in sp4c[i][k]
                    for (k = 0; k < 4; ++k) {
                        for (l = 0; l < 2; l++) {
                            if (first_6A_cluster[k] == cp[l]) break;
                        }
                        if (l == 2) {
                            if (m >= 4) {
                                m++;
                                break;
                            }
                            trial[m] = first_6A_cluster[k];
                            m++;
                        }
                    }
                    if (m != 4) continue;

                    for (k = 0; k < 4; ++k) { // find uncommon SP4 ring particles in sp4c[mem_sp4c[sp4c[i][j2]][j]][k]
                        for (l = 0; l < 2; l++) {
                            if (second_6A_cluster[k] == cp[l]) break;
                        }
                        if (l == 2) {
                            if (m >= 6) {
                                m++;
                                break;
                            }
                            trial[m] = second_6A_cluster[k];
                            m++;
                        }
                    }
                    if (m != 6) continue;

                    quickSort(&trial[2], 4);

                    Cluster_Write_9k(trial);
                }
            }
        }
    }
}

void Cluster_Write_9k(const int trial[]) {
    int i;
    int clusSize = 9;

    if(n9K == m9K) {
        hc9K= resize_2D_int(hc9K, m9K, m9K + incrStatic, clusSize, -1);
        m9K= m9K + incrStatic;
    }
    for (i=0; i<9; i++) hc9K[n9K][i]=trial[i];

    for(i=0; i<6; i++) {
        if (s9K[hc9K[n9K][i]] == 'C') s9K[hc9K[n9K][i]] = 'B';
    }
    s9K[hc9K[n9K][6]] = 'O';
    s9K[hc9K[n9K][7]] = 'O';
    s9K[hc9K[n9K][8]] = 'O';

    ++n9K;
}