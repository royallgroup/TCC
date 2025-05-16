#include "simple_cluster_methods.h"
#include "6MW.h"
#include "10MW.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"


// find two five-membered rings, each of whuch share four particles
// and one four-membered ring which shares three particles with each five-membered ring


void Clusters_Get6MW() {
    for (int fst_5a_id = 0; fst_5a_id < nsp5a; ++fst_5a_id) {
        int *fst_5a_cluster = hcsp5a[fst_5a_id];
        //printf("\n sp5a "); for (int i = 0; i < 5; ++i) printf("%i", fst_5a_cluster[i]);
        for (int scd_5a_id = 0; scd_5a_id < nsp5a; ++scd_5a_id) {
            int *scd_5a_cluster = hcsp5a[scd_5a_id];
            if(shared_5a(fst_5a_cluster, scd_5a_cluster) == 4){
                for (int fst_4a_id = 0; fst_4a_id < nsp4a; ++fst_4a_id) {
                    int *fst_4a_cluster = hcsp4a[fst_4a_id];
                    //printf("\n sp4a "); for (int i = 0; i < 4; ++i) printf("%i", fst_4a_cluster[i]);
                    if(shared_4a5a(fst_4a_cluster, fst_5a_cluster) == 3 && shared_4a5a(fst_4a_cluster, scd_5a_cluster) == 3){
                        //printf("found a 6MW fst_5a %i scd_5a %i fst_4a %i \n", fst_5a_id, scd_5a_id, fst_4a_id);
                        add_6MW(fst_5a_cluster, scd_5a_cluster);
                    }
                }
            }
        }
    }
}

/*int shared_5a(int *clust5a1, int *clust5a2){
    int u = 0;
    for (int i = 0; i < 5; ++i){
        for (int j = 0; j < 5; ++j){
            if(clust5a1[i] == clust5a2[j]){
                u +=1;
            }
        }
    }
    //if (u == 2){
    //    return 0;
    //}
    return u;
}*/



/*int shared_4a5a(int *clust4a, int *clust5a){
    int u = 0;
    for (int i = 0; i < 4; ++i){
        for (int j = 0; j < 5; ++j){
            if(clust4a[i] == clust5a[j]){
                ++u; // +=1;
            }
        }
    }
    /*if (u == 2){
        printf("shared_4a5a \n 4a ");
        for (int i = 0; i < 4; ++i) printf("%i ", clust4a[i]);
        printf("\n 5a ");
        for (int i = 0; i < 5; ++i) printf("%i ", clust5a[i]);
        return 1;
    }
    if (u == 0){
        return 0;
    }*/
    /*return u;
}*/



void add_6MW(const int *clust5a1, int *clust5a2) {
    int clusSize = 6;
    int in2in1[5];

    if (n6MW == m6MW) {
        hc6MW = resize_2D_int(hc6MW, m6MW, m6MW + incrStatic, clusSize, -1);
        m6MW = m6MW + incrStatic;
    }

    hc6MW[n6MW][0] = clust5a1[0];
    hc6MW[n6MW][1] = clust5a1[1];
    hc6MW[n6MW][2] = clust5a1[2];
    hc6MW[n6MW][3] = clust5a1[3];
    hc6MW[n6MW][4] = clust5a1[4];

    for (int i = 0; i < 5; ++i){
        in2in1[i] = 0;
        for (int j = 0; j < 5; ++j){
            if(clust5a2[i] == clust5a1[j]){
              in2in1[i] = 1;
            }
        }
    }
    for (int i = 0; i < 5; ++i){
        if(in2in1[i] == 0) // this is the particle which is not in the first 5 membered ring but is in the second
          hc6MW[n6MW][5] = clust5a2[i];
    }


    for (int i = 0; i < 6; ++i) {
        s6MW[hc6MW[n6MW][i]] = 'B';
    }
    ++n6MW;
}
