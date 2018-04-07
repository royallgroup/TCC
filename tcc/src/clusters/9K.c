#include <globals.h>
#include <tools.h>
#include "9K.h"

void Clusters_Get9K() {
    // 9K are made from 2 6A clusters with a common spindle and two common SP4 ring particles
    int j2, j, k, l, m;
    int cp[2], scom, sother[2];
    int trial[9];
    int id_first_6A, id_second_6A;

    cp[0]=cp[1]=scom=sother[0]=sother[1]=-1;


    for(id_first_6A=0; id_first_6A < nsp4c - 1; ++id_first_6A) {   // loop over all sp4c_i
        for (j2=4; j2<6; j2++) {    // loop over all spindles of sp4c_i
            for (j=0; j < nmem_sp4c[hcsp4c[id_first_6A][j2]]; ++j) {
                id_second_6A = mem_sp4c[hcsp4c[id_first_6A][j2]][j];
                if (id_second_6A<=id_first_6A) continue; // don't find again 9K twice

                m=0;        // sp4c_i and sp4c_mem_sp4c[sp4c[i][j2]][j] have exactly one common spindle
                for(k=4; k<6; ++k) {
                    for(l=4; l<6; ++l) {
                        if(hcsp4c[id_first_6A][k] == hcsp4c[id_second_6A][l]) {
                            m++;
                            scom= hcsp4c[id_first_6A][k];
                        }
                    }
                }
                if(m!=1) continue;

                if (hcsp4c[id_first_6A][4] == scom) sother[0]= hcsp4c[id_first_6A][5];
                else sother[0]= hcsp4c[id_first_6A][4];
                if (hcsp4c[id_second_6A][4] == scom) sother[1]= hcsp4c[id_second_6A][5];
                else sother[1]= hcsp4c[id_second_6A][4];

                m=0;        // check sother[0] is not in cluster sp4c_mem_sp4c[sp4c[i][j2]][j]
                for(k=0; k<6; ++k) {
                    if(sother[0] == hcsp4c[id_second_6A][k] || sother[1] == hcsp4c[id_first_6A][k]) {
                        m++;
                    }
                }
                if(m!=0) continue;

                m=0;        // SP4 ring from sp4c_i and SP4 ring from sp4c_mem_sp4c[sp4c[i][j2]][j] have exactly two common particles
                for(k=0; k<4; ++k) {
                    for(l=0; l<4; ++l) {
                        if(hcsp4c[id_first_6A][k] == hcsp4c[id_second_6A][l]) {
                            if (m>=2) {
                                m=3;
                                break;
                            }
                            cp[m]= hcsp4c[id_first_6A][k];
                            m++;
                        }
                    }
                    if (m>2) break;
                }
                if(m!=2) continue;

                // hc9K key: (common_SP4_1, common_SP4_2, other_SP4*4, other_spindle_1, other_spindle_2, scom)

                trial[0]=cp[0];
                trial[1]=cp[1];
                trial[6]=sother[0];
                trial[7]=sother[1];
                trial[8]=scom;
                quickSort(&trial[0],2);
                quickSort(&trial[6],2);

                m=2;        // find uncommon SP4 ring particles in sp4c[i][k]
                for(k=0; k<4; ++k) {
                    for (l=0; l<2; l++) {
                        if (hcsp4c[id_first_6A][k] == cp[l]) break;
                    }
                    if (l==2) {
                        if (m>=4) {
                            m++;
                            break;
                        }
                        trial[m]= hcsp4c[id_first_6A][k];
                        m++;
                    }
                }
                if(m!=4) continue;

                for(k=0; k<4; ++k) { // find uncommon SP4 ring particles in sp4c[mem_sp4c[sp4c[i][j2]][j]][k]
                    for (l=0; l<2; l++) {
                        if (hcsp4c[id_second_6A][k] == cp[l]) break;
                    }
                    if (l==2) {
                        if (m>=6) {
                            m++;
                            break;
                        }
                        trial[m]= hcsp4c[id_second_6A][k];
                        m++;
                    }
                }
                if(m!=6) continue;

                quickSort(&trial[2],4);

                Cluster_Write_9k(trial);
            }
        }
    }
}

void Cluster_Write_9k(const int trial[]) {
    int i;
    int clusSize = 9;

    if(n9K == m9K) {
        hc9K= resize_2D_int(hc9K, m9K, m9K + incrStatic, clusSize, -1);
        m9K= m9K + incrStatic;
    }
    for (i=0; i<9; i++) hc9K[n9K][i]=trial[i];

    for(i=0; i<6; i++) {
        if (s9K[hc9K[n9K][i]] == 'C') s9K[hc9K[n9K][i]] = 'B';
    }
    s9K[hc9K[n9K][6]] = 'O';
    s9K[hc9K[n9K][7]] = 'O';
    s9K[hc9K[n9K][8]] = 'O';

    ++n9K;
}