#include "10W.h"
#include "globals.h"
#include "tools.h"

void Clusters_Get10W() { // Detect 10W clusters
    int i, j, k, l, m;
    int sp5b_clusts[5], shell_parts[9];
    int clusSize=10;

    sp5b_clusts[0]=sp5b_clusts[1]=sp5b_clusts[2]=sp5b_clusts[3]=sp5b_clusts[4]=-1;


    for (i=0; i < nsp5b; ++i) { // loop over all sp5b
        if (num_bonds[hcsp5b[i][5]] != 9) continue;   // central particle must have coordination number 9

        k=0;    // find 5 other sp5b's with spindle in common with sp5b_i
        for (j=0; j < nmem_sp5b[hcsp5b[i][5]]; ++j) { // note check that spindle of sp5b_i and sp5b_j must be common by later check
            if (mem_sp5b[hcsp5b[i][5]][j] <= i) continue;   // i for sp5b must be lowest of all sp5b indices
            // ERROR !! need to check that spindle of sp5b_j is spindle of sp5b_i
            if (k>=5) {
                k++;
                break;
            }
            sp5b_clusts[k]= mem_sp5b[hcsp5b[i][5]][j];
            k++;
        }
        if (k!=5) continue; // not correct number of sp5b clusters
        // now found exactly 5 sp5b clusters common to spindle of sp5b_i
        for (j=0; j<5; j++) {
            shell_parts[j]= hcsp5b[i][j];
        }

        m=5;
        for (j=0; j<5; j++) {
            for (k=0; k<5; k++) {
                for (l=0; l<m; l++) {
                    if (shell_parts[l] == hcsp5b[mem_sp5b[hcsp5b[i][5]][j]][k]) break;
                }
                if (l==m) {
                    if (m>=9) {
                        m++;
                        break;
                    }
                    shell_parts[m]= hcsp5b[mem_sp5b[hcsp5b[i][5]][j]][k];
                    m++;
                }
            }
            if (m>=10) break;
        }
        if (m!=9) continue; // not all coordination shell particles of sp5b[i][5] are in the SP5 rings of the 5xsp5b clusters we found

        if (n10W == m10W) {
            hc10W= resize_2D_int(hc10W, m10W, m10W + incrStatic, clusSize, -1);
            m10W= m10W + incrStatic;
        }
        // hc10W key: (sp5bs_common_central_spindle_particle, sp5bs_SP5_ring_shell_particles)
        hc10W[n10W][0] = hcsp5b[i][5];
        for (j=0; j<9; j++) hc10W[n10W][j + 1]=shell_parts[j];
        quickSort(&hc10W[n10W][1], 9);
        Cluster_Write_10W();
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