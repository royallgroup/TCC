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
*      - The two additional particles are bonded to each other and to the common spindle particle of 9B.
*      - Each additional particle is bonded to two more particles in the shell of the 9B cluster,
*        leading to a total of four bonds between the additional particles and 9B.
*      - For each additional particle, the two shell particles to which they are bonded are not themselves bonded.
*      - The four shell particles of the 9B cluster that are bonded to the two additional particles form two pairs that are neighbours.
*
*  Cluster output: BBBBBBOOSBB
*  Storage order: as_for_9B x 9, extra_particles x 2
*/
int Clusters_Get11B() {

    int b1[2], b2[2], nb1, nb2;
    int extra_particles[10]; // The id's of the two extra particles
    int flg11, flg12, flg21, flg22;
    int break_out;

    int *parent_9B_cluster = hc9B[n9B];
    int center_9B_id = parent_9B_cluster[8];

    if(num_bonds[parent_9B_cluster[8]] != 10) {
        return 0; // s_com has 10 bonds in total (all forming the shell)
    }

    int *potential_11B = bond_list[center_9B_id];

    // Check that there are exactly 2 particles bonded to the 9B center that are not in the 9B.
    if(count_uncommon_particles(potential_11B, parent_9B_cluster, 10, 8, extra_particles) != 2) {
        return 0;
    }

    // Check the two extra particles are bonded to each other
    if (Bonds_BondCheck(extra_particles[0], extra_particles[1]) == 0) {
        return 0;
    }

    // Check that both of the extra particles are bonded to exactly two 9B shell particles.
    nb1 = nb2 = 0;
    for(int k=0; k<8; ++k){
        if(Bonds_BondCheck(parent_9B_cluster[k], extra_particles[0])){
            if(nb1 == 2) return 0;  // extra particle 1 bonded to 2 members of 9B shell
            b1[nb1]=parent_9B_cluster[k];
            nb1++;
        }
        if(Bonds_BondCheck(parent_9B_cluster[k], extra_particles[1])){
            if(nb2 == 2) return 0;  // extra particle 2 bonded to 2 members of 9B shell
            b2[nb2]=parent_9B_cluster[k];
            nb2++;
        }
    }
    if(nb1 != 2 || nb2 != 2) return 0;  // extra particles 1 and 2 bonded to 2 members of 9B shell

    // Particles bonded to extra 2 particles b[] must be distinct.
    if (b1[0] == b2[0] || b1[0] == b2[1] || b1[1] == b2[0] || b1[1] == b2[1]) {
        return 0;
    }
    if(Bonds_BondCheck(b1[0], b1[1]) || Bonds_BondCheck(b2[0], b2[1])) {
        return 0;
    }

    flg11 = Bonds_BondCheck(b1[0], b2[0]);
    flg12 = Bonds_BondCheck(b1[0], b2[1]);
    flg21 = Bonds_BondCheck(b1[1], b2[0]);
    flg22 = Bonds_BondCheck(b1[1], b2[1]);
    if(!((flg11 && !flg12) || (!flg11 && flg12))) return 0;
    if(!((flg21 && !flg22) || (!flg21 && flg22))) return 0;

    Cluster_Write_11B(extra_particles);

    return 1;
}

void Cluster_Write_11B(int *extra_particles) {
    int clusSize=11;

    if(n11B == m11B) {
        hc11B=resize_2D_int(hc11B,m11B,m11B+incrStatic,clusSize,-1);
        m11B=m11B+incrStatic;
    }
    for (int i = 0; i < 9; ++i) {
        hc11B[n11B][i] = hc9B[n9B][i];
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
}
