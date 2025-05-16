#include <globals.h>
#include <tools.h>
#include "8MW.h"
#include "simple_cluster_methods.h"



// an 8MW is a cube with the 8 particles at the vertices
// check to see if the four four-membered rings each share two particles
//



void Clusters_Get8MW() {
    for (int fst_4a_id = 0; fst_4a_id < nsp4a; ++fst_4a_id) {
        int *fst_4a_cluster = hcsp4a[fst_4a_id];
            for (int scd_4a_id = 0; scd_4a_id < nsp4a; ++scd_4a_id) {
                int *scd_4a_cluster = hcsp4a[scd_4a_id];
                if(shared_4a(fst_4a_cluster, scd_4a_cluster) == 2){ // then we share two particles between fst and scd
                    for (int thd_4a_id = 0; thd_4a_id < nsp4a; ++thd_4a_id) {
                        int *thd_4a_cluster = hcsp4a[thd_4a_id];

                        if(shared_4a(fst_4a_cluster, thd_4a_cluster) == 0 && shared_4a(scd_4a_cluster, thd_4a_cluster) == 2){
                            for (int fth_4a_id = 0; fth_4a_id < nsp4a; ++fth_4a_id) {
                                int *fth_4a_cluster = hcsp4a[fth_4a_id];

                                if(shared_4a(fst_4a_cluster, fth_4a_cluster) == 2 && shared_4a(scd_4a_cluster, fth_4a_cluster) == 0
                                && shared_4a(thd_4a_cluster, fth_4a_cluster) == 2){
                                    //printf("%i %i %i %i \n", fst_4a_id, scd_4a_id, thd_4a_id, fth_4a_id);
                                    add_8MW(fst_4a_cluster, thd_4a_cluster);

                                }
                            }
                        }
                    }
                }
            }
        }
    }



/*int shared_4a(int *clust4a1, int *clust4a2){
    int u = 0;
    for (int i = 0; i < 4; ++i){
        for (int j = 0; j < 4; ++j){
            if(clust4a1[i] == clust4a2[j]){
                u +=1;
            }
        }
    }
    /*if (u == 2){
        return 1;
    }
    if (u == 0){
        return 0;
    }
    return 2;*/
/*
    return u;
}*/


void add_8MW(const int *clust4a1, int *clust4a2) {
    int clusSize = 8;

    if (n8MW == m8MW) {
        hc8MW = resize_2D_int(hc8MW, m8MW, m8MW + incrStatic, clusSize, -1);
        m8MW = m8MW + incrStatic;
    }

    hc8MW[n8MW][0] = clust4a1[0];
    hc8MW[n8MW][1] = clust4a1[1];
    hc8MW[n8MW][2] = clust4a1[2];
    hc8MW[n8MW][3] = clust4a1[3];
    hc8MW[n8MW][4] = clust4a2[0];
    hc8MW[n8MW][5] = clust4a2[1];
    hc8MW[n8MW][6] = clust4a2[2];
    hc8MW[n8MW][7] = clust4a2[3];


    for (int i = 0; i < 8; ++i) {
        s8MW[hc8MW[n8MW][i]] = 'B';
    }
    ++n8MW;
}
