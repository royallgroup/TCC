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
   *  Cluster output: SBBBBBBBB
   *  Storage order: central_spindle_particle, shell_particles x 9
   */

    int shell_ids[9], neighbouring_sp5_ids[5];


    for (int first_sp5b_id = 0; first_sp5b_id < nsp5b; ++first_sp5b_id) {
        int *first_sp5b_cluster = hcsp5b[first_sp5b_id];
        int center_id = first_sp5b_cluster[5];  // The id of the shared spindle at the center of the 10W
        if (num_bonds[center_id] == 9) {   // central particle must have coordination number 9

            int num_shared_sp5b = 0;    // find 5 other sp5b's with spindle in common with sp5b_i
            for (int other_sp5b_pointer = 0; other_sp5b_pointer < nmem_sp5b[center_id]; ++other_sp5b_pointer) {
                int other_sp5_id = mem_sp5b[center_id][other_sp5b_pointer];
                if (other_sp5_id > first_sp5b_id) {
                    if (num_shared_sp5b < 9) {
                        neighbouring_sp5_ids[num_shared_sp5b] = other_sp5_id;
                    }
                    num_shared_sp5b++;
                }
            }
            if (num_shared_sp5b == 5) {
                for (int i = 0; i < 5; i++) {
                    shell_ids[i] = first_sp5b_cluster[i];
                }

                int num_shell_particles = get_shell_particle_ids(shell_ids, neighbouring_sp5_ids);
                if (num_shell_particles == 9) {
                    Cluster_Write_10W(center_id, shell_ids);
                }
            }
        }
    }
}

int get_shell_particle_ids(int *shell_ids, const int *neighbouring_sp5_ids) {
    int num_shell_particles = 5;
    for (int neighbouring_sp5_pointer = 0; neighbouring_sp5_pointer < 5; neighbouring_sp5_pointer++) {
        int *neighbouring_cluster = hcsp5b[neighbouring_sp5_ids[neighbouring_sp5_pointer]];
        for (int k = 0; k < 5; k++) {
            int shell_pointer;
            for (shell_pointer = 0; shell_pointer < num_shell_particles; shell_pointer++) {
                if (shell_ids[shell_pointer] == neighbouring_cluster[k]) {
                    break;
                }
            }
            if (shell_pointer == num_shell_particles) {
                if (num_shell_particles >= 9) {
                    num_shell_particles++;
                    break;
                }
                shell_ids[num_shell_particles] = neighbouring_cluster[k];
                num_shell_particles++;
            }
        }
        if (num_shell_particles >= 10) break;
    }
    return num_shell_particles;
}

void Cluster_Write_10W(int center_id, int *shell_ids) {
    int clusSize=10;

    if (n10W == m10W) {
        hc10W = resize_2D_int(hc10W, m10W, m10W + incrStatic, clusSize, -1);
        m10W = m10W + incrStatic;
    }

    hc10W[n10W][0] = center_id;
    for (int i = 0; i < 9; i++) {
        hc10W[n10W][i + 1] = shell_ids[i];
    }

    for (int i = 1; i < 10; i++) {
        if (s10W[hc10W[n10W][i]] == 'C') s10W[hc10W[n10W][i]] = 'B';
    }
    s10W[hc10W[n10W][0]] = 'S';

    ++n10W;
}