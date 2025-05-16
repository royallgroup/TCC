#include <globals.h>
#include <tools.h>
#include "9S.h"
#include "8PAB.h"
#include "7PAB.h"
#include "7T.h"
#include "11S.h"
//!  A 9S is a 8B with an additional particle
/*!
* A 9S cluster contains a a 8B and a 5A
*
* A 9S takes these two 6Z clusters and tries to make a 7T from the cluster and an extra particle.
*  Cluster output: BBBBBBBBB
*  Storage order: original 8PAB particles x 8, new 5A spindle)
*/

    void Clusters_Get9S() {
    for (int new_8B_id = 0; new_8B_id < n8B; ++new_8B_id) {
        int new_part;
        int *new_8B_cluster = hc8B[new_8B_id];
            for (int new_5A_id = 0; new_5A_id< nsp3c; ++new_5A_id) {
                int *new_5A_cluster = hcsp3c[new_5A_id];
                    if (common_spindle(new_8B_cluster, new_5A_cluster, 5, 3) == 1) { // the 8A and 5A share a spindle
                    if (common_ring(new_8B_cluster, new_5A_cluster, 5, 3) == 2) {
                        if (is_particle_ring_of_5A(new_8B_cluster[7], new_5A_cluster) == 1) { // the extra particle of the 8B is a ring of the 5A
                            new_part = get_new_part(new_8B_cluster, new_5A_cluster);
                            add_9S(new_8B_cluster, new_part);
                            }
                        }
                    }
                }
            }
        }

int get_new_part(const int *clust1, const int *clust2){
    if(clust1[5] == clust2[3] || clust1[6] == clust2[3]){
        return clust2[4];
    }
    if(clust1[5] == clust2[4] || clust1[6] == clust2[4]){
        return clust2[3];
    }
}

void add_9S(const int *old_8B, int new_part) {
    int clusSize = 9;
    if (n9S == m9S) {
        hc9S = resize_2D_int(hc9S, m9S, m9S + incrStatic, clusSize, -1);
        m9S = m9S + incrStatic;
    }
    hc9S[n9S][0] = old_8B[0];
    hc9S[n9S][1] =old_8B[1];
    hc9S[n9S][2] = old_8B[2];
    hc9S[n9S][3] = old_8B[3];
    hc9S[n9S][4] = old_8B[4];
    hc9S[n9S][5] = old_8B[5];
    hc9S[n9S][6] = old_8B[6];
    hc9S[n9S][7] = old_8B[7];
    hc9S[n9S][8] = new_part;

    for (int i = 0; i < 9; ++i) {
        s9S[hc9S[n9S][i]] = 'B';
    }
    ++n9S;
}
