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
*      -  Both of the distinct spindles of 9B are in the ring of the 7A.
*      -  Two 7A ring particles are common with the distinct sp5 particles of 9B.
*      -  There is one 7A ring particle which is distinct from the 9B cluster.
*
*  Cluster output: BBBBBBOOOS
*  Storage order: ordered_shell_particles x 6, spindles x 3, common_spindle
*/
void Clusters_Get10B(int j) {
    int trial[10];
    int secondary_7A_spindle;

    for (int second_7A = j + 1; second_7A < nsp5c; ++second_7A) {
        int *second_7A_cluster = hcsp5c[second_7A];
        int *parent_9B = hc9B[n9B];

        // The parent 9B has two spindle particles
        int first_9B_spindle = parent_9B[6];
        int second_9B_spindle = parent_9B[7];

        // At least one 7A spindle must be the common spindle of 9B - this is the primary spindle
        if (second_7A_cluster[5] == parent_9B[8]) {
            secondary_7A_spindle = second_7A_cluster[6];
        }
        else if (second_7A_cluster[6] == parent_9B[8]) {
            secondary_7A_spindle = second_7A_cluster[5];
        }
        else {
            // If the 7A has no common spindles then try the next 7A
            continue;
        }

        // Check the other 7A spindle is bonded to both 9B spindles.
        if (Bonds_BondCheck(secondary_7A_spindle, first_9B_spindle) == 0) continue;
        if (Bonds_BondCheck(secondary_7A_spindle, second_9B_spindle) == 0) continue;

        // Check that both of the 9B spindles are in the 7A ring
        if (is_particle_in_cluster(second_7A_cluster, 5, first_9B_spindle) == 0) continue;
        if (is_particle_in_cluster(second_7A_cluster, 5, second_9B_spindle) == 0) continue;

        // Build the trial cluster
        trial[6] = first_9B_spindle;
        trial[7] = second_9B_spindle;
        trial[8] = secondary_7A_spindle;
        trial[9] = parent_9B[8];

        // Get the 7A ring particles from the 9B.
        int particle_count = 0;
        for (int pointer_9B = 0; pointer_9B < 6; pointer_9B++) {
            if (parent_9B[pointer_9B] != secondary_7A_spindle) {
                trial[particle_count] = parent_9B[pointer_9B];
                particle_count++;
            }
        }
        if (particle_count != 5) continue;

        // Find the distinct particle from the 7A
        int particle_found = 0;
        for (int l = 0; l < 5; l++) {
            if (second_7A_cluster[l] == first_9B_spindle) continue;
            if (second_7A_cluster[l] == second_9B_spindle) continue;

            if (is_particle_in_cluster(trial, 5, second_7A_cluster[l]) == 0) {
                trial[5] = second_7A_cluster[l];
                particle_found++;
            }
        }

        // Now we have found the 10B C3v cluster
        // ###### NOTE #####
        // we have sterically assumed that
        // 1) one member of the SP5 ring of 7A_k is common with one member of SP5 ring of 7A_i
        // (this member of 7A_i was uncommon to the SP5 ring of 7A_j)
        // 2) one member of the SP5 ring of 7A_k is common with one member of SP5 ring of 7A_j
        // (this member of 7A_j was uncommon to the SP5 ring of 7A_i)

        if (particle_found == 1) {
            Cluster_Write_10B(trial);
        }
    }
}

void Cluster_Write_10B(int trial[10]) {
    int clusSize=10;

    if (n10B == m10B) {
        hc10B = resize_2D_int(hc10B, m10B, m10B + incrStatic, clusSize, -1);
        m10B = m10B + incrStatic;
    }

    for (int i = 0; i < 10; i++) {
        hc10B[n10B][i] = trial[i];
    }

    for (int i = 0; i < 6; i++) {
        if (s10B[hc10B[n10B][i]] == 'C') s10B[hc10B[n10B][i]] = 'B';
    }
    for (int i = 6; i < 9; i++) {
        if (s10B[hc10B[n10B][i]] != 'S') s10B[hc10B[n10B][i]] = 'O';
    }
    s10B[hc10B[n10B][9]] = 'S';
    ++n10B;
}