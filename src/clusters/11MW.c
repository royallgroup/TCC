#include <globals.h>
#include <tools.h>
#include "11MW.h"
#include "simple_cluster_methods.h"



// a 11MW is an 10MW with one particle attached
// (in the same way as 9MW to 8MW)
// so there are 4 sp5as and 3 sp4as
// the first two sp5as share 3 particles (these correspond to the sp5as for 9MW)
// and these first two each share 3 particles with the 2 "face" sp5as 3 and 4 which do not share any particles with each other

// then there are 3 sp4as
// the first shares 2 particles with sp5a 1,3,4
// the second shares 2 particles with sp5a 2,3,4
// the third shares 2 particles with sp5a 3,4 and sp4a 1 and 2



void Clusters_Get11MW() {


    for (int fst_5a_id = 0; fst_5a_id < nsp5a; ++fst_5a_id){
        int *fst_5a = hcsp5a[fst_5a_id];
        for (int scd_5a_id = 0; scd_5a_id < nsp5a; ++scd_5a_id){
            int *scd_5a = hcsp5a[scd_5a_id];
            if(shared_5a(fst_5a, scd_5a) == 3){ // then we share 3 particles between fst and scd sp5a
                for (int thd_5a_id = 0; thd_5a_id < nsp5a; ++thd_5a_id){
                    int *thd_5a = hcsp5a[thd_5a_id];
                    if(shared_5a(fst_5a, thd_5a) == 2 && shared_5a(scd_5a, thd_5a) == 2){
                        // and 3 particles betweem each of the fst scd and thd
                        for (int for_5a_id = 0; for_5a_id < nsp5a; ++for_5a_id){
                            int *for_5a = hcsp5a[for_5a_id];
                            if(shared_5a(fst_5a, for_5a) == 2 && shared_5a(scd_5a, for_5a) == 2
                            && shared_5a(thd_5a, for_5a) == 0){
                                // and 3 particles betweem each of the fst scd and for but none between the thd and for
                                //printf("\nwe are good for sp5as\n");
                                // time for the sp4as
                                for (int fst_4a_id = 0; fst_4a_id < nsp4a; ++fst_4a_id){
                                    int *fst_4a = hcsp4a[fst_4a_id];
                                    if(shared_4a5a(fst_4a, fst_5a) == 2 && shared_4a5a(fst_4a, thd_5a) == 2
                                    && shared_4a5a(fst_4a, for_5a)==2){ // first sp4a is good
                                        for (int scd_4a_id = 0; scd_4a_id < nsp4a; ++scd_4a_id){
                                            int *scd_4a = hcsp4a[scd_4a_id];
                                            if(shared_4a5a(scd_4a, scd_5a) == 2 && shared_4a5a(scd_4a, thd_5a) == 2
                                            && shared_4a5a(scd_4a, for_5a)==2){ // second sp4a is good
                                                for (int thd_4a_id = 0; thd_4a_id < nsp4a; ++thd_4a_id){
                                                    int *thd_4a = hcsp4a[thd_4a_id];
                                                    if(shared_4a5a(thd_4a, thd_5a) == 2 && shared_5a(thd_4a, for_5a) == 2
                                                    && shared_4a(thd_4a, fst_4a)   == 2 && shared_4a(thd_4a, scd_4a) == 2){ // third sp4a is good
                                                        add_11MW(fst_5a, thd_5a, for_5a);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


}


void add_11MW(const int *clust5a1, const int *clust5a3, const int *clust5a4){
// sp5a 3 and 4 are the "sides"
// so all their particles are counted
// the remaining particle is that which is in the first two sp5as 1 and 2 which is not
// in either of the side sp5as


    int clusSize = 11;

    if (n11MW == m11MW) {
        hc11MW = resize_2D_int(hc11MW, m11MW, m11MW + incrStatic, clusSize, -1);
        m11MW  = m11MW + incrStatic;
    }

    // add the third sp5a3 (one of the face sp5as)
    hc11MW[n11MW][0] = clust5a3[0];
    hc11MW[n11MW][1] = clust5a3[1];
    hc11MW[n11MW][2] = clust5a3[2];
    hc11MW[n11MW][3] = clust5a3[3];
    hc11MW[n11MW][4] = clust5a3[4];


    // add the fourth sp5a4 (the other of the face sp5as)
    hc11MW[n11MW][5] = clust5a4[0];
    hc11MW[n11MW][6] = clust5a4[1];
    hc11MW[n11MW][7] = clust5a4[2];
    hc11MW[n11MW][8] = clust5a4[3];
    hc11MW[n11MW][9] = clust5a4[4];


    // and now we need to find the particle in the first sp5a which is in neither of the other two
    for(int i = 0; i < 5; ++i){
        if(is_particle_in_cluster(clust5a3, 5, clust5a1[i])==0 && is_particle_in_cluster(clust5a4, 5, clust5a1[i])==0) {
            hc11MW[n11MW][10] = clust5a1[i];
        }
    }


    for (int i = 0; i < 11; ++i) {
        s11MW[hc11MW[n11MW][i]] = 'B';
    }
    ++n11MW;
}
