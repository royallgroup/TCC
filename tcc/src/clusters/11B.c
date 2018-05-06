#include "11B.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

int Clusters_Get11B() {

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

    int k, l, m;
    int b1[2], b2[2], nb1, nb2;
    int ep[2]; // The two extra particles
    int flg11, flg12, flg21, flg22;
    int break_out;
    int clusSize=11;

    if(num_bonds[hc9B[n9B][8]]!= 10) return 0; // s_com has 10 bonds in total (all forming the shell)

    m = 0;
    break_out=0;
    for(k=0; k<10; ++k) {
        for(l=0; l<8; ++l) {
            if(bond_list[hc9B[n9B][8]][k] == hc9B[n9B][l]) break;
        }
        if(l==8){
            if(m==2) {
                break_out=1;
                break;
            }
            ep[m++] = bond_list[hc9B[n9B][8]][k];    // two extra particles
        }
    }
    if(break_out==1 || m<2) return 0;

    if(Bonds_BondCheck(ep[0], ep[1])==0) return 0;  // extra particles must be bonded
    nb1 = nb2 = 0;
    for(k=0; k<8; ++k){
        if(Bonds_BondCheck(hc9B[n9B][k], ep[0])){
            if(nb1 == 2) return 0;  // extra particle 1 bonded to 2 members of 9B shell
            b1[nb1]=hc9B[n9B][k];
            nb1++;
        }
        if(Bonds_BondCheck(hc9B[n9B][k], ep[1])){
            if(nb2 == 2) return 0;  // extra particle 2 bonded to 2 members of 9B shell
            b2[nb2]=hc9B[n9B][k];
            nb2++;
        }
    }
    if(nb1 != 2 || nb2 != 2) return 0;  // extra particles 1 and 2 bonded to 2 members of 9B shell

    flg11 = b1[0] == b2[0] || b1[0] == b2[1]; // Particles bonded to extra 2 particles b[]
    flg11 = flg11 || b1[1] == b2[0] || b1[1] == b2[1]; // must be distinct.
    if(flg11) return 0;
    flg11 = Bonds_BondCheck(b1[0], b1[1]); // paritcles b1[] mustn't be bonded
    flg22 = Bonds_BondCheck(b2[0], b2[1]);
    if(flg11 || flg22) return 0;
    flg11 = Bonds_BondCheck(b1[0], b2[0]);
    flg12 = Bonds_BondCheck(b1[0], b2[1]);
    flg21 = Bonds_BondCheck(b1[1], b2[0]);
    flg22 = Bonds_BondCheck(b1[1], b2[1]);
    if(!((flg11 && !flg12) || (!flg11 && flg12))) return 0;
    if(!((flg21 && !flg22) || (!flg21 && flg22))) return 0;

    if(n11B == m11B) {
        hc11B=resize_2D_int(hc11B,m11B,m11B+incrStatic,clusSize,-1);
        m11B=m11B+incrStatic;
    }
    hc11B[n11B][0] = hc9B[n9B][0];
    hc11B[n11B][1] = hc9B[n9B][1];
    hc11B[n11B][2] = hc9B[n9B][2];
    hc11B[n11B][3] = hc9B[n9B][3];
    hc11B[n11B][4] = hc9B[n9B][4];
    hc11B[n11B][5] = hc9B[n9B][5];
    hc11B[n11B][6] = hc9B[n9B][6];
    hc11B[n11B][7] = hc9B[n9B][7];
    hc11B[n11B][8] = hc9B[n9B][8];
    if (b1[0]==1) {
        hc11B[n11B][9] = ep[0];
        hc11B[n11B][10] = ep[1];
    }
    else {
        hc11B[n11B][9] = ep[1];
        hc11B[n11B][10] = ep[0];
    }
    Cluster_Write_11B();

    return 1;
}

void Cluster_Write_11B() {

    int i;
    for(i=0; i<6; i++) {
        if (s11B[hc11B[n11B][i]] == 'C') s11B[hc11B[n11B][i]] = 'B';
    }
    if(s11B[hc11B[n11B][6]] != 'S') s11B[hc11B[n11B][6]] = 'O';
    if(s11B[hc11B[n11B][7]] != 'S') s11B[hc11B[n11B][7]] = 'O';
    s11B[hc11B[n11B][8]] = 'S';
    if(s11B[hc11B[n11B][9]] == 'C') s11B[hc11B[n11B][9]] = 'B';
    if(s11B[hc11B[n11B][10]] == 'C') s11B[hc11B[n11B][10]] = 'B';
}
