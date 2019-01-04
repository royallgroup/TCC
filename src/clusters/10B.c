#include <globals.h>
#include <bonds.h>
#include <tools.h>
#include <clusters/simple_cluster_methods.h>
#include "10B.h"

//!  An 10B cluster is the intersection of a 9B and a 7A cluster.
/*!
*  Find 10B clusters
*  An 10A is a 9B and 7A cluster where:
*      -  One spindle from 7A is common to the common spindle particle of the 9B cluster
*      -  The other spindle from 7A is bonded to the two distinct spindles of 9B.
*      -  Two sp5 ring particles from 7A are common with the distinct spindles of 9B.
*      -  Two sp5 ring particles from 7A are common with the distinct sp5 particles of 9B.
*      -  The final sp5 ring particle from 7A is distinct from the 9B cluster.
*
*  Cluster output: BBBBBBOOOS
*  Storage order: ordered_shell_particles x 6, spindles x 3, common_spindle
*/
void Clusters_Get10B(int j) {        // Return 1 if 9B is also 10B cluster
    int l,m;
    int flg1, flg2;
    int trial[10];
    int break_out;
    int clusSize=10;

    for (int second_7A = j + 1; second_7A < nsp5c; ++second_7A) {  // loop over all 7A_k
        int *second_7A_cluster = hcsp5c[second_7A];
        int *parent_9B = hc9B[n9B];
        // At least one 7A spindle must be the common spindle of 9B
        if (second_7A_cluster[5] == parent_9B[8]) {
            // The parent 9B has two spindle particles
            int first_9B_spindle = parent_9B[6];
            int second_9B_spindle = parent_9B[7];

            // Check the 9B spindles are bonded to the 7A spindle.
            if (Bonds_BondCheck(second_7A_cluster[6], first_9B_spindle) == 0) continue;
            if (Bonds_BondCheck(second_7A_cluster[6], second_9B_spindle) == 0) continue;

            // Check that both of the 9B spindles are in the 7A ring
            if (is_particle_in_cluster(second_7A_cluster, 5, first_9B_spindle) == 0) continue;
            if (is_particle_in_cluster(second_7A_cluster, 5, second_9B_spindle) == 0) continue;

            // Build the trial cluster
            trial[6] = first_9B_spindle;
            trial[7] = second_9B_spindle;
            trial[8] = second_7A_cluster[6];
            trial[9] = parent_9B[8];

            m = 0;
            break_out = 0;
            for (l = 0; l < 6; l++) {
                if (parent_9B[l] == second_7A_cluster[6]) continue;
                if (m == 5) {
                    m++;
                    break_out = 1;
                    break;
                }
                trial[m] = parent_9B[l];
                m++;
            }
            if (break_out == 1 || m != 5) continue;

            break_out = 0;
            for (l = 0; l < 5; l++) {
                if (second_7A_cluster[l] == first_9B_spindle) continue;
                if (second_7A_cluster[l] == second_9B_spindle) continue;
                for (m = 0; m < 5; m++) {
                    if (second_7A_cluster[l] == trial[m]) break;
                }
                if (m == 5) {
                    trial[5] = second_7A_cluster[l];
                    break_out++;
                }
            }
            if (break_out != 1) continue;

            if (n10B == m10B) {
                hc10B = resize_2D_int(hc10B, m10B, m10B + incrStatic, clusSize, -1);
                m10B = m10B + incrStatic;
            }
            // Now we have found the 10B C3v cluster
            // ###### NOTE #####
            // we have sterically assumed that
            // 1) one member of the SP5 ring of 7A_k is common with one member of SP5 ring of 7A_i
            // (this member of 7A_i was uncommon to the SP5 ring of 7A_j)
            // 2) one member of the SP5 ring of 7A_k is common with one member of SP5 ring of 7A_j
            // (this member of 7A_j was uncommon to the SP5 ring of 7A_i)

            quickSort(&trial[0], 6);
            quickSort(&trial[6], 3);
            for (l = 0; l < 10; l++) hc10B[n10B][l] = trial[l];

            Cluster_Write_10B();
        }

        if (second_7A_cluster[6] == parent_9B[8]) {    // check one spindle of 7A_k is the common spindle of 9B (hc9B[id9B][.] at this point)
            if (Bonds_BondCheck(second_7A_cluster[5], parent_9B[6]) == 0) continue;  // check other spindle of 7A_k is bonded to spindle d1 of 9B
            if (Bonds_BondCheck(second_7A_cluster[5], parent_9B[7]) == 0) continue;  // check other spindle of 7A_k is bonded to spindle d2 of 9B

            flg1 = 0;
            flg2 = 0;
            for (l = 0; l < 5; l++) {
                if (second_7A_cluster[l] == parent_9B[6]) {
                    flg1 = 1;
                    continue;
                }
                if (second_7A_cluster[l] == parent_9B[7]) {
                    flg2 = 1;
                    continue;
                }
            }
            if (flg1 == 0 || flg2 == 0) continue;
            trial[6] = parent_9B[6];
            trial[7] = parent_9B[7];
            trial[8] = second_7A_cluster[5];
            trial[9] = parent_9B[8];

            m = 0;
            break_out = 0;
            for (l = 0; l < 6; l++) {
                if (parent_9B[l] == second_7A_cluster[5]) continue;
                if (m == 5) {
                    m++;
                    break_out = 1;
                    break;
                }
                trial[m] = parent_9B[l];
                m++;
            }
            if (break_out == 1 || m != 5) continue;

            break_out = 0;
            for (l = 0; l < 5; l++) {
                if (second_7A_cluster[l] == parent_9B[6]) continue;
                if (second_7A_cluster[l] == parent_9B[7]) continue;
                for (m = 0; m < 5; m++) {
                    if (second_7A_cluster[l] == trial[m]) break;
                }
                if (m == 5) {
                    trial[5] = second_7A_cluster[l];
                    break_out++;
                }
            }
            if (break_out != 1) continue;

            if (n10B == m10B) {
                hc10B = resize_2D_int(hc10B, m10B, m10B + incrStatic, clusSize, -1);
                m10B = m10B + incrStatic;
            }
            // Now we have found the 10B C3v cluster
            // ###### NOTE #####
            // we have sterically assumed that
            // 1) one member of the SP5 ring of 7A_k is common with one member of SP5 ring of 7A_i
            // (this member of 7A_i was uncommon to the SP5 ring of 7A_j)
            // 2) one member of the SP5 ring of 7A_k is common with one member of SP5 ring of 7A_j
            // (this member of 7A_j was uncommon to the SP5 ring of 7A_i)

            quickSort(&trial[0], 6);
            quickSort(&trial[6], 3);
            for (l = 0; l < 10; l++) hc10B[n10B][l] = trial[l];

            Cluster_Write_10B();
        }
    }
}

void Cluster_Write_10B() {

    for(int i=0; i<6; i++) {
        if (s10B[hc10B[n10B][i]] == 'C') s10B[hc10B[n10B][i]] = 'B';
    }
    for(int i=6; i<9; i++) {
        if (s10B[hc10B[n10B][i]] != 'S') s10B[hc10B[n10B][i]] = 'O';
    }
    s10B[hc10B[n10B][9]] = 'S';
    ++n10B;
}