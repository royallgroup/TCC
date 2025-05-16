#include <clusters/simple_cluster_methods.h>
#include "11B.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

//!  An 11B cluster is 9B with two additional particles.
/*!
*  Find 11B clusters
*  An 11A is a 9B and two extra particles where:
*      - The common spindle particle from the 9B cluster has coordination number 10.
*      - The two additional particles (a1 and a2) are bonded to each other and to the common spindle particle of 9B.
*      - a1 is bonded to two particles in the shell of the 9B cluster (b1 and b2).
*      - a2 is bonded to two particles in the shell of the 9B cluster (c1 and c2).
*      - b1 and b2 are not bonded. c1 and c2 are not bonded.
*      - One of b1 or b2 is bonded to one of c1 or c2. The other pair are also bonded.
*
*  Cluster output: BBBBBBOOSBB
*  Storage order: as_for_9B x 9, extra_particles x 2
*/
void Clusters_Get11B() {

    for (int parent_9B_id = 0; parent_9B_id < n9B; ++parent_9B_id) {

        int *parent_9B_cluster = hc9B[parent_9B_id];

        // Check that the central particle of the 9B has exactly 10 neighbours
        if(num_bonds[parent_9B_cluster[8]] != 10) {
            continue;
        }

        // Check that there are exactly 2 particles bonded to the 9B center that are not in the 9B.
        int *parent_9B_center_neighbours = bond_list[parent_9B_cluster[8]];
        int extra_particles[18];
        if(count_uncommon_particles(parent_9B_center_neighbours, parent_9B_cluster, 10, 8, extra_particles) != 2) {
            continue;
        }

        // Check the two extra particles are bonded to each other
        if (Bonds_BondCheck(extra_particles[0], extra_particles[1]) == 0) {
            continue;
        }

        // Check that both of the extra particles are bonded to exactly two 9B shell particles.
        int bonded_b[8], bonded_c[8], num_bonded_b, num_bonded_c;
        num_bonded_b = count_cluster_bonds_to_particle(extra_particles[0], parent_9B_cluster, 8, bonded_b);
        num_bonded_c = count_cluster_bonds_to_particle(extra_particles[1], parent_9B_cluster, 8, bonded_c);
        if(num_bonded_b != 2 || num_bonded_c != 2) {
            continue;
        }

        // Particles bonded to the extra particles must be distinct.
        if (are_clusters_distinct(bonded_b, bonded_c, 2, 2) == 0) {
            continue;
        }

        // Pairs of bonded particles must not be bonded.
        if(Bonds_BondCheck(bonded_b[0], bonded_b[1]) || Bonds_BondCheck(bonded_c[0], bonded_c[1])) {
            continue;
        }

        // There must be either bonds between (b1-c1 and b2-c2)
        // or bonds between b1-c2 and b2-c1
        int bond_b1_c1, bond_b1_c2, bond_b2_c1, bond_b2_c2;
        bond_b1_c1 = Bonds_BondCheck(bonded_b[0], bonded_c[0]);
        bond_b1_c2 = Bonds_BondCheck(bonded_b[0], bonded_c[1]);
        bond_b2_c1 = Bonds_BondCheck(bonded_b[1], bonded_c[0]);
        bond_b2_c2 = Bonds_BondCheck(bonded_b[1], bonded_c[1]);

        if(bond_b1_c1 && bond_b2_c2 && !bond_b1_c2 && !bond_b2_c1) {
            Cluster_Write_11B(extra_particles, parent_9B_id);
        }
        else if(bond_b1_c2 && bond_b2_c1 && !bond_b1_c1 && !bond_b2_c2) {
            Cluster_Write_11B(extra_particles, parent_9B_id);
        }
    }
}

void Cluster_Write_11B(const int *extra_particles, int parent_9B_id) {
    int clusSize=11;

    if(n11B == m11B) {
        hc11B=resize_2D_int(hc11B,m11B,m11B+incrStatic,clusSize,-1);
        m11B=m11B+incrStatic;
    }
    for (int i = 0; i < 9; ++i) {
        hc11B[n11B][i] = hc9B[parent_9B_id][i];
    }

    hc11B[n11B][9] = extra_particles[0];
    hc11B[n11B][10] = extra_particles[1];

    for(int i = 0; i < 6; i++) {
        if (s11B[hc11B[n11B][i]] == 'C') s11B[hc11B[n11B][i]] = 'B';
    }
    if(s11B[hc11B[n11B][6]] != 'S') s11B[hc11B[n11B][6]] = 'O';
    if(s11B[hc11B[n11B][7]] != 'S') s11B[hc11B[n11B][7]] = 'O';
    s11B[hc11B[n11B][8]] = 'S';
    if(s11B[hc11B[n11B][9]] == 'C') s11B[hc11B[n11B][9]] = 'B';
    if(s11B[hc11B[n11B][10]] == 'C') s11B[hc11B[n11B][10]] = 'B';
    n11B++;
}

