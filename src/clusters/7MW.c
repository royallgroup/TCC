#include <globals.h>
#include <tools.h>
#include "7MW.h"
#include "simple_cluster_methods.h"



// a 7MW is three sides of a cube each being an sp4a
// find 3 sp4as, each of which shares two particles with both of the other two
// and require one particle to be shared with both of the other two


void Clusters_Get7MW() {

    for (int fst_4a_id = 0; fst_4a_id < nsp4a; ++fst_4a_id) {
        int *fst_4a = hcsp4a[fst_4a_id];
            for (int scd_4a_id = 0; scd_4a_id < nsp4a; ++scd_4a_id) {
                int *scd_4a = hcsp4a[scd_4a_id];
                if(shared_4a(fst_4a, scd_4a) == 2){ // then we share two particles between fst and scd
                    for (int thd_4a_id = 0; thd_4a_id < nsp4a; ++thd_4a_id) {
                        int *thd_4a = hcsp4a[thd_4a_id];
                        if(shared_4a(fst_4a, thd_4a) == 2 && shared_4a(scd_4a, thd_4a) == 2 && shared_4a3(fst_4a, scd_4a, thd_4a) == 1){
                            add_7MW(fst_4a, scd_4a, thd_4a);
                        }
                    }
                }
            }
        }
    }



void add_7MW(const int *clust4a1, const int *clust4a2, const int *clust4a3){
    int clusSize = 7;
    int in2in1[4];     // identify the particles in the first and second sp4a
    int in3nin2nin1[4];


    if (n7MW == m7MW) {
        hc7MW = resize_2D_int(hc7MW, m7MW, m7MW + incrStatic, clusSize, -1);
        m7MW = m7MW + incrStatic;
    }

    hc7MW[n7MW][0] = clust4a1[0];
    hc7MW[n7MW][1] = clust4a1[1];
    hc7MW[n7MW][2] = clust4a1[2];
    hc7MW[n7MW][3] = clust4a1[3];

    for (int i = 0; i < 4; ++i){ // identify particles in the second sp4a which are also in the first
        in2in1[i] = 0;
        for (int j = 0; j < 4; ++j){
            if(clust4a2[i] == clust4a2[j]){
              in2in1[i] = 1;
            }
        }
    }
    int j = 0;
    for (int i = 0; i < 4; ++i){
        if(in2in1[i] == 0) // these are particles which are not in the first sp4a but are not in the second 4a
          hc7MW[n7MW][4+j] = clust4a2[i];
          ++j;
    }
    //printf("\ncounter j should be 2, it is %d\n", j);

    // get the particle in the third sp4a which is in neither the first sp4a nor the second sp4a

    for (int i = 0; i < 4; ++i)
        in3nin2nin1[0] = 0;

    for (int i = 0; i < 4; ++i){
        for (int j = 0; j < 4; ++j)
            if(clust4a3[i] == clust4a1[j])
              in3nin2nin1[i] = 1;

        for (int j = 0; j < 4; ++j)
            if(clust4a3[i] == clust4a2[j])
              in3nin2nin1[i] = 1;
    }

    // the remaining 0 of in3nin2nin1 is the particle in third sp4a which is in neither of the other two sp4a
    for (int i = 0; i < 4; ++i){
        if(in3nin2nin1[i] == 0) hc7MW[n7MW][6] = clust4a3[i];
    }


    for (int i = 0; i < 7; ++i) {
        s7MW[hc7MW[n7MW][i]] = 'B';
    }
    ++n7MW;
}
