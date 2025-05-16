#include <globals.h>
#include <tools.h>
#include "13MW.h"
#include "simple_cluster_methods.h"



// a 12MW is an irregular octahedron whose faces are comprised of four sp4as and five sp5as
// sp51 shares 2 with sp43, sp52, sp41, sp53, sp42
// sp52 shares 2 with sp45, sp43, sp41, sp54, [sp51]
// sp53 shares 2 with sp41, sp54, sp44, sp42, [sp51]
// sp54 shares 2 with sp41, sp45, sp44, [sp52, sp53]
// sp41 shares 2 with  [sp51, sp52, sp53, sp54]
// sp42 shares 2 with  sp44, sp43, [sp51, sp53]
// sp43 shares 2 with sp45, [sp42, sp51, sp52]
// sp44 shares 2 with sp45 [sp42, sp53, sp54]
// sp45 shares 2 with [sp43, sp44, sp52, sp54]

// here the square brackets denote that the edges shared have already been specified above



void Clusters_Get13MW() {

    for (int sp5a1_id = 0; sp5a1_id < nsp5a; ++sp5a1_id){
      int *sp5a1 = hcsp5a[sp5a1_id]; // we work backwards through the list above starting with sp5a4
      for (int sp4a1_id = 0; sp4a1_id < nsp4a; ++sp4a1_id){
        int *sp4a1 = hcsp4a[sp4a1_id];
        if(shared_4a5a(sp4a1, sp5a1) == 2){
          for (int sp4a2_id = 0; sp4a2_id < nsp4a; ++sp4a2_id){
            int *sp4a2 = hcsp4a[sp4a2_id];
            if(shared_4a5a(sp4a2, sp5a1) == 2){
              for (int sp4a3_id = 0; sp4a3_id < nsp4a; ++sp4a3_id){
                int *sp4a3 = hcsp4a[sp4a3_id];
                if(shared_4a5a(sp4a3, sp5a1) == 2){
                  for (int sp5a2_id = 0; sp5a2_id < nsp5a; ++sp5a2_id){
                    int *sp5a2 = hcsp5a[sp5a2_id];
                    if(shared_5a(sp5a2, sp5a1) == 2){
                      for (int sp5a3_id = 0; sp5a3_id < nsp5a; ++sp5a3_id){
                        int *sp5a3 = hcsp5a[sp5a3_id];
                        if(shared_5a(sp5a3, sp5a1) == 2){ // then sp5a1 is OK
                          //printf("\n Clusters_Get13MW sp5a1 is OK \n");
                          //sp52 shares 2 with sp45, sp43, sp41, sp54, [sp51]
                          if(shared_4a5a(sp4a3, sp5a2) == 2 && shared_4a5a(sp4a1, sp5a2) == 2){
                            for (int sp4a5_id = 0; sp4a5_id < nsp4a; ++sp4a5_id){
                              int *sp4a5 = hcsp4a[sp4a5_id];
                              if(shared_4a5a(sp4a5, sp5a2) == 2){
                                for (int sp5a4_id = 0; sp5a4_id < nsp5a; ++sp5a4_id){
                                  int *sp5a4 = hcsp5a[sp5a4_id];
                                  if(shared_5a(sp5a4, sp5a2) == 2){  // then sp5a2 is OK
                                    //printf("\n Clusters_Get13MW sp5a2 is OK \n");
                                    //sp53 shares 2 with sp41, sp54, sp44, sp42, [sp51]
                                    if(shared_4a5a(sp4a1, sp5a3) == 2 && shared_5a(sp5a4, sp5a3) == 2
                                    && shared_4a5a(sp4a2, sp5a3) == 2){
                                      for(int sp4a4_id = 0; sp4a4_id < nsp4a; ++sp4a4_id){
                                        int *sp4a4 = hcsp4a[sp4a4_id];
                                        if(shared_4a5a(sp4a4, sp5a3) == 2){  // then sp5a3 is OK
                                          //printf("\n Clusters_Get13MW sp5a3 is OK \n");
                                          //sp54 shares 2 with sp41, sp45, sp44, [sp52, sp53]
                                          if(shared_4a5a(sp4a1, sp5a4) == 2 && shared_5a(sp4a5, sp5a4) == 2
                                          && shared_4a5a(sp4a4, sp5a4) == 2){  // then sp5a4 is OK sp4a1 is dealt with anyway
                                            //printf("\n Clusters_Get13MW sp5a4 is OK \n");
                                            //sp42 shares 2 with  sp44, sp43, [sp51, sp53]
                                            if(shared_4a(sp4a2, sp4a4) == 2 && shared_4a(sp4a2, sp4a3) == 2){ // then sp4a2 is OK
                                              //printf("\n Clusters_Get13MW sp4a2 is OK \n");
                                              if(shared_4a(sp4a3, sp4a5) == 2 && shared_4a(sp4a4, sp4a5) == 2) // then sp4a3 and sp4a4 are OK
                                                add_13MW(sp5a1, sp5a4, sp4a2, sp4a3);
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
                  }
                }
              }
            }
          }
        }
      }



}


void add_13MW(const int *sp5a1, const int *sp5a4, const int *sp4a2, const int *sp4a3){

// sp5a 1 and 4 give us 10
// so all their particles are counted
// and sp4a2 and sp4a3 habe three particles, one each not in sp5a1 and one shared between them not in sp5a1



    int clusSize = 13;

    if (n13MW == m13MW) {
        hc13MW = resize_2D_int(hc13MW, m13MW, m13MW + incrStatic, clusSize, -1);
        m13MW  = m13MW + incrStatic;
    }

    // add the third sp5a3 (one of the face sp5as)
    // nov 23 -- n12MW --> n13MW seems more likely to be right!
    hc13MW[n13MW][0] = sp5a1[0];
    hc13MW[n13MW][1] = sp5a1[1];
    hc13MW[n13MW][2] = sp5a1[2];
    hc13MW[n13MW][3] = sp5a1[3];
    hc13MW[n13MW][4] = sp5a1[4];


    // add the fourth sp5a4 (the other of the face sp5as)
    // nov 23 -- n11MW --> n13MW seems more likely to be right!
    hc13MW[n13MW][5] = sp5a4[0];
    hc13MW[n13MW][6] = sp5a4[1];
    hc13MW[n13MW][7] = sp5a4[2];
    hc13MW[n13MW][8] = sp5a4[3];
    hc13MW[n13MW][9] = sp5a4[4];


    // and now we need to find the particle in the first sp5a which is in neither of the other two
    int sp4aextra=10;

    // get the particle in sp4a2 not in sp4a3 or s@5a1
    for(int i = 0; i < 4; ++i){
        if(is_particle_in_cluster(sp5a1, 5, sp4a2[i]) == 0 && is_particle_in_cluster(sp4a3, 4, sp4a2[i]) == 0) {
            // nov 23 -- n11MW --> n13MW seems more likely to be right!
            hc13MW[n13MW][sp4aextra] = sp4a2[i];
            ++sp4aextra;
        }
    }

    // get the particle in sp4a3 not in sp4a2 or s@5a1
    for(int i = 0; i < 4; ++i){
        if(is_particle_in_cluster(sp5a1, 5, sp4a3[i]) == 0 && is_particle_in_cluster(sp4a2, 4, sp4a3[i]) == 0) {
            hc13MW[n13MW][sp4aextra] = sp4a3[i];
            ++sp4aextra;
        }
    }

    // get the particle in both sp4a3 and in sp4a2 but not in s@5a1
    for(int i = 0; i < 4; ++i){
        if(is_particle_in_cluster(sp5a1, 5, sp4a3[i]) == 0 && is_particle_in_cluster(sp4a2, 4, sp4a3[i]) == 1) {
            hc13MW[n13MW][sp4aextra] = sp4a3[i];
            ++sp4aextra;
        }
    }


    //printf("\n add_13MW: sp4aextra should be 13, it is %d\n", sp4aextra);

    for (int i = 0; i < 13; ++i) {
        s13MW[hc13MW[n13MW][i]] = 'B';
    }
    ++n13MW;
}
