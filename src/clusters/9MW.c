#include <globals.h>
#include <tools.h>
#include "9MW.h"
#include "simple_cluster_methods.h"



// a 9MW is an 8MW with one particle attached
// so there are 4 sp4as
// if 1 and 3 are the sides, 1 and 3 share no particles
// but 1,2; 2,3; 2,4; 3,4 and 1,4 all share two particles
// is_particle_in_cluster each of which share two particles
// and two sp5a each of which have 2 particles with sp4a and 3 particles shared between the two sp5a
// each sp5a will share two particles with three out of the four sp4as


void Clusters_Get9MW() {

    for (int fst_4a_id = 0; fst_4a_id < nsp4a; ++fst_4a_id){
        int *fst_4a = hcsp4a[fst_4a_id];
        for (int scd_4a_id = 0; scd_4a_id < nsp4a; ++scd_4a_id){
            int *scd_4a = hcsp4a[scd_4a_id];
            if(shared_4a(fst_4a, scd_4a) == 2){ // then we share two particles between fst and scd sp4a
                for (int thd_4a_id = 0; thd_4a_id < nsp4a; ++thd_4a_id){
                    int *thd_4a = hcsp4a[thd_4a_id];
                    if(shared_4a(scd_4a, thd_4a) == 2){
                        for (int for_4a_id = 0; for_4a_id < nsp4a; ++for_4a_id){
                            int *for_4a = hcsp4a[for_4a_id];
                            if(shared_4a(for_4a, fst_4a) == 2 && shared_4a(for_4a, scd_4a) == 2 && shared_4a(for_4a, thd_4a) == 2){
                                //printf("\nGet9MW: all four sp4as seem to be OK\n");
                                for (int fst_5a_id = 0; fst_5a_id < nsp5a; ++fst_5a_id){
                                    int *fst_5a = hcsp5a[fst_5a_id];
                                    if(shared_4a5a(fst_4a, fst_5a) + shared_4a5a(scd_4a, fst_5a) + shared_4a5a(thd_4a, fst_5a)  + shared_4a5a(for_4a, fst_5a) == 6){
                                        // we will suppose this to be good enough
                                            for (int scd_5a_id = 0; scd_5a_id < nsp5a; ++scd_5a_id){
                                                int *scd_5a = hcsp5a[scd_5a_id];
                                                if(shared_4a5a(fst_4a, scd_5a) + shared_4a5a(scd_4a, scd_5a) + shared_4a5a(thd_4a, scd_5a)  + shared_4a5a(for_4a, scd_5a) == 6
                                                  && shared_5a(fst_5a, scd_5a)==3){
                                                    add_9MW(fst_5a, scd_5a, fst_4a, scd_4a);
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

void add_9MW(const int *clust5a1, const int *clust5a2, const int *clust4a1, const int *clust4a2){
// two sp5a share 3, so we add those exlcuding the repeated 3
// we then need from the sp4a the two particles that are shared with neiter of the sp5a

    int clusSize = 9;
    int nfound = 0;
    int in2in1[5];     // identify the particles in the first and second sp4a

    if (n9MW == m9MW) {
        hc9MW = resize_2D_int(hc9MW, m9MW, m9MW + incrStatic, clusSize, -1);
        m9MW = m9MW + incrStatic;
    }

    // add the first sp5a1
    hc9MW[n7MW][0] = clust5a1[0];
    hc9MW[n7MW][1] = clust5a1[1];
    hc9MW[n7MW][2] = clust5a1[2];
    hc9MW[n7MW][3] = clust5a1[3];
    hc9MW[n7MW][4] = clust5a1[4];
    nfound = 5;


    for (int i = 0; i < 5; ++i){  // identify particles in the second sp5a which are also in the first
        in2in1[i] = 0;            // which we don't want...but we need to know we don't want!
        for (int j = 0; j < 5; ++j){
            if(clust5a2[i] == clust5a1[j]){
              in2in1[i] = 1;
            }
        }
    }

    for (int i = 0; i < 5; ++i){
        if(in2in1[i] == 0){ // then we have a particle which is in the second sp5a which is not in the first
            hc9MW[n7MW][nfound] = clust5a2[i];
            ++nfound;
        }
    }



    // finally want the two sp4a particles which are not in either of the sp5as
    for(int i = 0; i < 4; ++i){
      if(is_particle_in_cluster(clust5a1, 5, clust4a1[i])==0 && is_particle_in_cluster(clust5a2, 5, clust4a1[i])==0) {
        hc9MW[n9MW][nfound] = clust4a1[i];
        ++nfound;
      }
    }

    if(nfound == 8){ // then we only found one particle in the sp4a so try with another 4a
      for(int i = 0; i < 4; ++i){
          if(is_particle_in_cluster(clust5a1, 5, clust4a2[i])==0 && is_particle_in_cluster(clust5a2, 5, clust4a2[i])==0
          && hc9MW[n9MW][7] != clust4a2[i]) {
            hc9MW[n9MW][nfound] = clust4a2[i];
            //++notinsp5a;
            ++nfound;
          }
      }
    }

  
    for (int i = 0; i < 9; ++i) {
        s9MW[hc9MW[n9MW][i]] = 'B';
    }
    ++n9MW;
}
