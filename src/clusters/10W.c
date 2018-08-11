#include "10W.h"
#include "globals.h"
#include "tools.h"

void Clusters_Get10W() {

    //!  An 10W cluster is the intersection six sp5b which all share one spindle.
    /*!
   *  Find 10W clusters
   *  An 10W is six sp5b clusters where:
   *      - All clusters share the same spindle
   *      - The coordination number of this common spindle is 9.
   *  The reason this construction is chosen over simpler ones is that it is robust to different values of the voronoi parameter from 0.82 to 1.
   *
   *  Cluster output: BBBBBBBBBS
   *  Storage order: central_spindle_particle, shell_particles x 9
   */

    int first_sp5b_id, j, k, l, m;
    int shell_parts[9];
    int clusSize=10;

    for (first_sp5b_id=0; first_sp5b_id < nsp5b; ++first_sp5b_id) {
        int *first_sp5b_cluster = hcsp5b[first_sp5b_id];
        int center_id = first_sp5b_cluster[5];  // The id of the shared spindle at the center of the 10W
        if (num_bonds[center_id] == 9) {   // central particle must have coordination number 9

            int num_shared_sp5b = 0;    // find 5 other sp5b's with spindle in common with sp5b_i
            for (int other_sp5b_pointer = 0; other_sp5b_pointer < nmem_sp5b[center_id]; ++other_sp5b_pointer) {
                int other_sp5_id = mem_sp5b[center_id][other_sp5b_pointer];
                if (other_sp5_id > first_sp5b_id) {
                    num_shared_sp5b++;
                }
            }
            if (num_shared_sp5b == 5) {
                // now found exactly 5 sp5b clusters common to spindle of sp5b_i
                for (j = 0; j < 5; j++) {
                    shell_parts[j] = first_sp5b_cluster[j];
                }

                m = 5;
                for (j = 0; j < 5; j++) {
                    for (k = 0; k < 5; k++) {
                        for (l = 0; l < m; l++) {
                            if (shell_parts[l] == hcsp5b[mem_sp5b[center_id][j]][k]) break;
                        }
                        if (l == m) {
                            if (m >= 9) {
                                m++;
                                break;
                            }
                            shell_parts[m] = hcsp5b[mem_sp5b[center_id][j]][k];
                            m++;
                        }
                    }
                    if (m >= 10) break;
                }
                if (m != 9)
                    continue; // not all coordination shell particles of sp5b[first_sp5b_id][5] are in the SP5 rings of the 5xsp5b clusters we found

                if (n10W == m10W) {
                    hc10W = resize_2D_int(hc10W, m10W, m10W + incrStatic, clusSize, -1);
                    m10W = m10W + incrStatic;
                }
                hc10W[n10W][0] = center_id;
                for (j = 0; j < 9; j++) hc10W[n10W][j + 1] = shell_parts[j];
                quickSort(&hc10W[n10W][1], 9);
                Cluster_Write_10W();
            }
        }
    }
}

void Cluster_Write_10W() {
    int i;

    for(i=1; i<10; i++) {
        if (s10W[hc10W[n10W][i]] == 'C') s10W[hc10W[n10W][i]] = 'B';
    }
    s10W[hc10W[n10W][0]] = 'S';

    ++n10W;
}