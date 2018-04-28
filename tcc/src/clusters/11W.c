#include "simple_cluster_methods.h"
#include "11W.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

void Clusters_Get11W() {

    //!  An 11W cluster is a 10B with an extra particle
    /*!
   *  Find 11W clusters
   *  An 11W is constructed from a 10B and an extra particle where:
   *      - The common spindle of the 10B cluster has coordination number 10.
   *      - The additional particle is not bonded to any of the distinct spindles of the 7A clusters
   *        constituting the 10B cluster.
   *
   *  Cluster output: BBBBBBBBBBO
   *  Storage order: as_for_10B x 10, extra_particle
   *
   */

    for(int first_10B_id = 0; first_10B_id < n10B; first_10B_id++) {
        int *first_10B_cluster = hc10B[first_10B_id];
        int first_10B_spindle_id = first_10B_cluster[9];
        if (num_bonds[first_10B_spindle_id] == 10) {

            int extra_particle = get_11W_extra_particle(first_10B_cluster, first_10B_spindle_id);

            // extra particle must not be bonded to three 7A spindles in shell of 10B
            if (is_particle_bonded_to_7As(first_10B_id, extra_particle) == 0) {
                Write_11W(first_10B_id, extra_particle);
            }
        }
    }
}

int is_particle_bonded_to_7As(int id_10B, int extra_particle) {

    for(int i = 6; i < 9; i++) {
        if (Bonds_BondCheck(extra_particle, hc10B[id_10B][i]) == 1){
            return 1;
        }
    }
    return 0;
}

int get_11W_extra_particle(int *parent_10B_cluster, int spindle_10B) {

    for (int i = 0; i < num_bonds[spindle_10B]; ++i) {
        int extra_particle = bNums[spindle_10B][i];
        if (is_particle_in_cluster(parent_10B_cluster, 10, extra_particle) == 0) {
            return extra_particle;
        }
    }
    Error("11W extra particle not found");
    return 0;
}

void Write_11W(int id_10B, int extra_particle) {
    int clusSize=11;

    if (n11W == m11W) {
        hc11W = resize_2D_int(hc11W, m11W, m11W + incrStatic, clusSize, -1);
        m11W = m11W + incrStatic;
    }

    for (int i = 0; i < 10; i++){
        hc11W[n11W][i] = hc10B[id_10B][i];
    }
    hc11W[n11W][10] = extra_particle;

    for(int i = 0; i < 10; i++) {
        if (s11W[hc11W[n11W][i]] == 'C') s11W[hc11W[n11W][i]] = 'B';
    }
    s11W[hc11W[n11W][10]] = 'O';

    n11W++;
}

