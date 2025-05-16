#include <globals.h>
#include <tools.h>
#include "10MW.h"
#include "simple_cluster_methods.h"
//#include "8MW.h"


// find two five-membered rings which don't share any particles
// but which have two particles each commmon to five four-membered rings



void Clusters_Get10MW() {
    for (int fst_5a_id = 0; fst_5a_id < nsp5a; ++fst_5a_id) {
        int *fst_5a_cluster = hcsp5a[fst_5a_id];
        //printf("\n sp5a "); for (int i = 0; i < 5; ++i) printf("%i", fst_5a_cluster[i]);
        for (int scd_5a_id = 0; scd_5a_id < nsp5a; ++scd_5a_id) {
            int *scd_5a_cluster = hcsp5a[scd_5a_id];
            if(shared_5a(fst_5a_cluster, scd_5a_cluster) == 0){ // the five-membered rings share no particles
                for (int fst_4a_id = 0; fst_4a_id < nsp4a; ++fst_4a_id) {
                    int *fst_4a_cluster = hcsp4a[fst_4a_id];
                    //printf("\n fst_sp4a "); for (int i = 0; i < 4; ++i) printf("%i", fst_4a_cluster[i]);
                    if(shared_4a5a(fst_4a_cluster, fst_5a_cluster) == 2 && shared_4a5a(fst_4a_cluster, scd_5a_cluster) == 2){
                        //printf("\nhave 2 sp5a sharing 1 sp4a\n");
                        for (int scd_4a_id = 0; scd_4a_id < nsp4a; ++scd_4a_id) {
                            int *scd_4a_cluster = hcsp4a[scd_4a_id];
                            //printf("\n scd_sp4a "); for (int i = 0; i < 4; ++i) printf("%i", scd_4a_cluster[i]);
                            if(shared_4a5a(scd_4a_cluster, fst_5a_cluster) == 2 && shared_4a5a(scd_4a_cluster, scd_5a_cluster) == 2 && shared_4a(scd_4a_cluster, fst_4a_cluster) == 2){
                                //printf("\nhave 2 sp5a sharing 2 sp4a\n");
                                for (int thd_4a_id = 0; thd_4a_id < nsp4a; ++thd_4a_id) {
                                    int *thd_4a_cluster = hcsp4a[thd_4a_id];
                                    //printf("\n thd_sp4a "); for (int i = 0; i < 4; ++i) printf("%i", thd_4a_cluster[i]);
                                    if(shared_4a5a(thd_4a_cluster, fst_5a_cluster) == 2 && shared_4a5a(thd_4a_cluster, scd_5a_cluster) == 2 && shared_4a(thd_4a_cluster, scd_4a_cluster) == 2){
                                        //printf("\nhave 2 sp5a sharing 3 sp4a\n");
                                        for (int for_4a_id = 0; for_4a_id < nsp4a; ++for_4a_id) {
                                            int *for_4a_cluster = hcsp4a[for_4a_id];
                                            if(shared_4a5a(for_4a_cluster, fst_5a_cluster) == 2 && shared_4a5a(for_4a_cluster, scd_5a_cluster) == 2 && shared_4a(for_4a_cluster, thd_4a_cluster) == 2){
                                                for (int fif_4a_id = 0; fif_4a_id < nsp4a; ++fif_4a_id) {
                                                    int *fif_4a_cluster = hcsp4a[fif_4a_id];
                                                    if(shared_4a5a(fif_4a_cluster, fst_5a_cluster) == 2 && shared_4a5a(fif_4a_cluster, scd_5a_cluster) == 2 && shared_4a(fif_4a_cluster, for_4a_cluster) == 2 && shared_4a(fif_4a_cluster, fst_4a_cluster) == 2){
                                                      //printf("found a 10MW %i %i \n", fst_5a_id, scd_5a_id);
                                                      add_10MW(fst_5a_cluster, scd_5a_cluster);
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





void add_10MW(const int *clust5a1, int *clust5a2) {
    int clusSize = 10;

    if (n10MW == m10MW) {
        hc10MW = resize_2D_int(hc10MW, m10MW, m10MW + incrStatic, clusSize, -1);
        m10MW = m10MW + incrStatic;
    }

    hc10MW[n10MW][0] = clust5a1[0];
    hc10MW[n10MW][1] = clust5a1[1];
    hc10MW[n10MW][2] = clust5a1[2];
    hc10MW[n10MW][3] = clust5a1[3];
    hc10MW[n10MW][4] = clust5a1[4];
    hc10MW[n10MW][5] = clust5a2[0];
    hc10MW[n10MW][6] = clust5a2[1];
    hc10MW[n10MW][7] = clust5a2[2];
    hc10MW[n10MW][8] = clust5a2[3];
    hc10MW[n10MW][9] = clust5a2[4];


    for (int i = 0; i < 10; ++i) {
        s10MW[hc10MW[n10MW][i]] = 'B';
    }
    ++n10MW;
}
