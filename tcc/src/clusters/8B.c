#include "simple_cluster_methods.h"
#include "8B.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

void Clusters_Get8B() {

    //!  An 8B cluster is an 7A cluster with an extra particle attached.
    /*!
   *  Find 8B clusters
   *  An 8B is a 7A cluster with an extra particle bonded to two of the ring particles and a spindle of the 7A cluster
   *
   *  Cluster output: BBBBOOB
   *  Storage order: as for 7A x 7, extra_particle)
   */

    int first_7A_id;
    int* first_7A_cluster;
    int primary_spindle;

    int new_particle_pointer;
    int new_particle_id;

    for (first_7A_id = 0; first_7A_id < nsp5c; ++first_7A_id) {
        first_7A_cluster = hcsp5c[first_7A_id];
        for (int spindle_pointer = 0; spindle_pointer < 2; ++spindle_pointer) {
            primary_spindle = first_7A_cluster[5 + spindle_pointer];

            for (new_particle_pointer = 0; new_particle_pointer < num_bonds[primary_spindle]; ++new_particle_pointer) {
                new_particle_id = bNums[primary_spindle][new_particle_pointer];

                if(is_particle_in_cluster(first_7A_cluster, 7, new_particle_id) == 0) {

                    if (count_bonds_to_7A_ring(first_7A_id, new_particle_id) == 2) {

                        Cluster_Write_8B(first_7A_cluster, new_particle_id);
                    }
                }
            }
        }
    }
}

int count_bonds_to_7A_ring(int first_7A_id, int new_particle_id) {
    int num_bonds_to_ring = 0;
    for (int ring_pointer = 0; ring_pointer < 5; ++ring_pointer) {
        if (Bonds_BondCheck(new_particle_id, hcsp5c[first_7A_id][ring_pointer])) {
            ++num_bonds_to_ring;
        }
    }
    return num_bonds_to_ring;
}

void Cluster_Write_8B(int *first_7A_cluster, int new_particle_id) {
    int clusSize = 8;

    // Now we have found the 8B Cs cluster
    if (n8B == m8B) {
        hc8B = resize_2D_int(hc8B, m8B, m8B + incrStatic, clusSize, -1);
        m8B = m8B + incrStatic;
    }

    for (int i = 0; i < 7; i++) {
        hc8B[n8B][i] = first_7A_cluster[i];
    }
    hc8B[n8B][7] = new_particle_id;

    if (s8B[hc8B[n8B][7]] == 'C') s8B[hc8B[n8B][7]] = 'B';
    if (s8B[hc8B[n8B][0]] == 'C') s8B[hc8B[n8B][0]] = 'B';
    if (s8B[hc8B[n8B][1]] == 'C') s8B[hc8B[n8B][1]] = 'B';
    if (s8B[hc8B[n8B][2]] == 'C') s8B[hc8B[n8B][2]] = 'B';
    if (s8B[hc8B[n8B][3]] == 'C') s8B[hc8B[n8B][3]] = 'B';
    if (s8B[hc8B[n8B][4]] == 'C') s8B[hc8B[n8B][4]] = 'B';
    s8B[hc8B[n8B][5]] = 'O';
    s8B[hc8B[n8B][6]] = 'O';
    ++n8B;
}