#include <clusters/simple_cluster_methods.h>
#include "8A.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"
#include "string.h"

void Clusters_Get8A() {

    //!  An 8A cluster is one of 3 possible topological combinations of sp5b/c clusters.
    /*!
   *  Find 8A clusters
   *  There are 3 methods used for 8A detection
   *  - A pair of sp5b where the spindles are distinct and share 4 ring particles
   *  - A pair of 7A clusters where
   *      - Both 7Ai spindle particles are common with the 7Aj spindles.
   *      - There are four common particles between sp5 rings of 7Ai and 7Aj .
   *  - A sp5b cluster and a 7A cluster where:
   *      - One 7A spindle is common with the sp5b spindle.
   *      - The other 7A spindle is distinct from all the sp5b particles.
   *      - There are four common particles between sp5 rings of sp5b and 7A.
   *
   *  Cluster output: BBBBBBOO
   *  Storage order: spindles x 4, not spindles x 4)
   */

    method_1();
    method_2();
    method_3();
}

void method_1() {

    int uncommon_ring_particle_ids[10];
    int com[4];
    int second_sp5b_pointer, first_sp5b_ring_pointer, k, l;

    int flg;

    int trial[8];
    int *used_sp5b;

    used_sp5b = malloc(nsp5b * sizeof(int));
    if (used_sp5b==NULL) {
        Error("Clusters_Get8A(): used_sp5b[] malloc out of memory\n");
    }

    for (int first_sp5b_id = 0; first_sp5b_id < nsp5b; ++first_sp5b_id) {  // loop over all sp5b_i
        int *first_sp5b_cluster = hcsp5b[first_sp5b_id];
        memset(used_sp5b, 0, nsp5b * sizeof(*used_sp5b));
        used_sp5b[first_sp5b_id] = 1;
        for (first_sp5b_ring_pointer = 0; first_sp5b_ring_pointer < 5; ++first_sp5b_ring_pointer) {
            for (second_sp5b_pointer=0; second_sp5b_pointer < nmem_sp5b[first_sp5b_cluster[first_sp5b_ring_pointer]]; ++second_sp5b_pointer) {  // loop over all sp5b_j
                int second_sp5b_id = mem_sp5b[first_sp5b_cluster[first_sp5b_ring_pointer]][second_sp5b_pointer];
                int *second_sp5b_cluster = hcsp5b[second_sp5b_id];


                if (second_sp5b_id > first_sp5b_id) {
                    if (used_sp5b[second_sp5b_id] == 0) {
                        used_sp5b[second_sp5b_id] = 1;

                        // Check rings share 4 particles
                        if (count_common_ring_particles(first_sp5b_cluster, second_sp5b_cluster, 5, com) == 4) {

                            // Check for distinct spindles
                            if (first_sp5b_cluster[5] != second_sp5b_cluster[5]) {

                                count_uncommon_ring_particles(first_sp5b_cluster, second_sp5b_cluster, 5, 5, uncommon_ring_particle_ids);

                                // build up trial cluster
                                trial[0] = first_sp5b_cluster[5];
                                trial[1] = second_sp5b_cluster[5];
                                trial[2] = uncommon_ring_particle_ids[0];
                                trial[3] = uncommon_ring_particle_ids[1];
                                trial[4] = com[0];
                                trial[5] = com[1];
                                trial[6] = com[2];
                                trial[7] = com[3];

                                int clusSize = 8;


                                quickSort(&trial[0], 8);
                                flg = 0;  // check trial cluster not already found
                                for (k = 0; k < n8A; ++k) {
                                    for (l = 0; l < 8; ++l) {
                                        if (trial[l] != hc8A[k][l]) break;
                                    }
                                    if (l == 8) {
                                        flg = 1;
                                    }
                                }
                                if (flg == 0) {

                                    // Now we have found the 8A D2d cluster
                                    if (n8A == m8A) {
                                        hc8A = resize_2D_int(hc8A, m8A, m8A + incrStatic, clusSize, -1);
                                        m8A = m8A + incrStatic;
                                    }

                                    for (k = 0; k < 8; ++k) {
                                        hc8A[n8A][k] = trial[k];
                                    }

                                    Cluster_Write_8A();

                                }
                            }
                        }
                    }
                }
            }
        }
    }
    free(used_sp5b);
}

void method_2() {

    int unc[2];
    int spindle_ids[2];
    int com[5];
    int first_7A_id, second_7A_pointer, first_7A_spindle_pointer, k, l, m;
    int cnt;
    int flg;
    int break_out;
    int trial[8];
    int clusSize = 8;

    for (first_7A_id = 0; first_7A_id < nsp5c; ++first_7A_id) {
        int *first_7A_cluster = hcsp5c[first_7A_id];
        for (first_7A_spindle_pointer = 5; first_7A_spindle_pointer < 6; ++first_7A_spindle_pointer) {
            for (second_7A_pointer = 0; second_7A_pointer < nmem_sp5c[first_7A_cluster[first_7A_spindle_pointer]]; ++second_7A_pointer) {  // loop over all 7A_j
                int second_7A_id = mem_sp5c[first_7A_cluster[first_7A_spindle_pointer]][second_7A_pointer];
                int *second_7A_cluster = hcsp5c[second_7A_id];

                if (second_7A_id <= first_7A_id) {

                    // exactly four members of the SP5 rings of 7A_i and 7A_j in common
                    if (count_common_ring_particles(first_7A_cluster, second_7A_cluster, 5, com) == 4) {

                        if (count_common_spindle_particles(first_7A_cluster, second_7A_cluster, 7, 7, spindle_ids) == 2) {
                            
                            for (k = 0; k < 5; ++k) {
                                m = 0;
                                for (l = 0; l < 4; ++l) {
                                    if (first_7A_cluster[k] == com[l]) m++;
                                }
                                if (m == 0) unc[0] = first_7A_cluster[k];
                            }
                            for (k = 0; k < 5; ++k) {
                                m = 0;
                                for (l = 0; l < 4; ++l) {
                                    if (second_7A_cluster[k] == com[l]) m++;
                                }
                                if (m == 0) unc[1] = second_7A_cluster[k];
                            }

                            // Now we have found the 8A D2d cluster
                            if (n8A == m8A) {
                                hc8A = resize_2D_int(hc8A, m8A, m8A + incrStatic, clusSize, -1);
                                m8A = m8A + incrStatic;
                            }
                            trial[0] = first_7A_cluster[5];    // build up trial cluster
                            trial[1] = first_7A_cluster[6];
                            trial[4] = unc[0];
                            trial[5] = unc[1];

                            cnt = 2;
                            break_out = 0;
                            for (k = 0; k < 5; ++k) {
                                if (Bonds_BondCheck(first_7A_cluster[k], trial[4]) == 1 &&
                                    first_7A_cluster[k] != trial[4] &&
                                    first_7A_cluster[k] != trial[5]) {
                                    if (cnt == 4) {
                                        break_out = 1;
                                        break;
                                    }
                                    trial[cnt] = first_7A_cluster[k];
                                    cnt++;
                                }
                            }
                            if (break_out == 1 || cnt < 4) continue;

                            for (k = 0; k < 5; ++k) {
                                if (Bonds_BondCheck(first_7A_cluster[k], trial[2]) == 1 &&
                                    first_7A_cluster[k] != trial[2] &&
                                    first_7A_cluster[k] != trial[4] && first_7A_cluster[k] != trial[5]) {
                                    trial[6] = first_7A_cluster[k];
                                }
                                if (Bonds_BondCheck(first_7A_cluster[k], trial[3]) == 1 &&
                                    first_7A_cluster[k] != trial[3] &&
                                    first_7A_cluster[k] != trial[4] && first_7A_cluster[k] != trial[5]) {
                                    trial[7] = first_7A_cluster[k];
                                }
                            }

                            quickSort(&trial[0], 4);
                            quickSort(&trial[4], 4);
                            flg = 0;  // check trial cluster not already found
                            for (k = 0; k < n8A; ++k) {
                                for (l = 0; l < 8; ++l) {
                                    if (trial[l] != hc8A[k][l]) break;
                                }
                                if (l == 8) flg = 1;
                            }
                            if (flg == 0) {
                                for (k = 0; k < 8; ++k) hc8A[n8A][k] = trial[k];

                                Cluster_Write_8A();
                            }
                        }
                    }
                }
            }
        }
    }
}

void method_3() {

    int unc[2];
    int com[4];
    int i, j, j2, k, l, m;
    int cnt;
    int flg;
    int break_out;
    int trial[8];
    int clusSize = 8;

    for (i=0; i < nsp5b; ++i) {    // loop over all sp5b_i
        for (j2=5; j2<6; ++j2) {
            for (j=0; j<nmem_sp5c[hcsp5b[i][j2]]; ++j) {  // loop over all 7A_j
                m = 0;
                for (k=0; k<5; ++k) {
                    for (l=0; l<5; ++l) {
                        if (hcsp5b[i][k] == hcsp5c[mem_sp5c[hcsp5b[i][j2]][j]][l]) {
                            if (m<5) com[m]=hcsp5b[i][k];
                            ++m;
                        }
                    }
                }
                if (m!=4) continue; // exactly four members of the SP5 rings of sp5b_i and 7A_j in common

                flg = hcsp5b[i][5] == hcsp5c[mem_sp5c[hcsp5b[i][j2]][j]][5] || hcsp5b[i][5] == hcsp5c[mem_sp5c[hcsp5b[i][j2]][j]][6];
                if (flg!=1) continue;   // sp5b_i spindle common with one of 7A_j spindles

                for (k=0; k<5; ++k) {
                    m=0;
                    for (l=0; l<4; ++l) {
                        if (hcsp5b[i][k]==com[l]) m++;
                    }
                    if (m==0) unc[0]=hcsp5b[i][k];
                }
                for (k=0; k<5; ++k) {
                    m=0;
                    for (l=0; l<4; ++l) {
                        if (hcsp5c[mem_sp5c[hcsp5b[i][j2]][j]][k]==com[l]) m++;
                    }
                    if (m==0) unc[1]=hcsp5c[mem_sp5c[hcsp5b[i][j2]][j]][k];
                }

                // Now we have found the 8A D2d cluster
                if (n8A==m8A) {
                    hc8A=resize_2D_int(hc8A,m8A,m8A+incrStatic,clusSize,-1);
                    m8A=m8A+incrStatic;
                }
                trial[0]=hcsp5c[mem_sp5c[hcsp5b[i][j2]][j]][5]; // build up trial cluster
                trial[1]=hcsp5c[mem_sp5c[hcsp5b[i][j2]][j]][6];
                trial[4]=unc[0];
                trial[5]=unc[1];

                cnt=2;
                break_out=0;
                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(hcsp5b[i][k],trial[4])==1 && hcsp5b[i][k]!=trial[4] && hcsp5b[i][k]!=trial[5]) {
                        if (cnt==4) {
                            break_out=1;
                            break;
                        }
                        trial[cnt]=hcsp5b[i][k];
                        cnt++;
                    }
                }
                if (break_out==1 || cnt<4) continue;

                for (k=0; k<5; ++k) {
                    if (Bonds_BondCheck(hcsp5b[i][k],trial[2])==1 && hcsp5b[i][k]!=trial[2] && hcsp5b[i][k]!=trial[4] && hcsp5b[i][k]!=trial[5]) {
                        trial[6]=hcsp5b[i][k];
                    }
                    if (Bonds_BondCheck(hcsp5b[i][k],trial[3])==1 && hcsp5b[i][k]!=trial[3] && hcsp5b[i][k]!=trial[4] && hcsp5b[i][k]!=trial[5]) {
                        trial[7]=hcsp5b[i][k];
                    }
                }

                quickSort(&trial[0],4);
                quickSort(&trial[4],4);
                flg=0;  // check trial cluster not already found
                for (k=0; k<n8A; ++k) {
                    for (l=0; l<8; ++l) {
                        if (trial[l]!=hc8A[k][l]) break;
                    }
                    if (l==8) flg=1;
                }
                if (flg==0) {
                    for (k=0; k<8; ++k) hc8A[n8A][k]=trial[k];
                    Cluster_Write_8A();
                }
            }
        }
    }
}

void Cluster_Write_8A() {// hc8A key: (4 of 8A_possible_spindles increasing, 4 of 8A_not_possible_spindles increasing)
    if (s8A[hc8A[n8A][0]] == 'C') s8A[hc8A[n8A][0]] = 'B';
    if (s8A[hc8A[n8A][1]] == 'C') s8A[hc8A[n8A][1]] = 'B';
    if (s8A[hc8A[n8A][2]] == 'C') s8A[hc8A[n8A][2]] = 'B';
    if (s8A[hc8A[n8A][3]] == 'C') s8A[hc8A[n8A][3]] = 'B';
    if (s8A[hc8A[n8A][4]] == 'C') s8A[hc8A[n8A][4]] = 'B';
    if (s8A[hc8A[n8A][5]] == 'C') s8A[hc8A[n8A][5]] = 'B';
    s8A[hc8A[n8A][6]] = 'O';
    s8A[hc8A[n8A][7]] = 'O';

    ++n8A;
}
