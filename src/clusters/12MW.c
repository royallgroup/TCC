#include <globals.h>
#include <tools.h>
#include "12MW.h"
#include "simple_cluster_methods.h"



// a 12MW is an irregular octahedron whose faces are comprised of four sp4as and five sp5as
// sp51 shares 2 with sp41, sp42, sp43, sp52, sp54
// sp41 shares 2 with sp43, sp52, sp53, [sp51]
// sp42 shares 2 with sp54, sp44, sp52, [sp51]
// sp52 shares 2 with sp44, [sp51, sp42, sp41, sp53]
// sp43 shares 2 with sp53, sp54, [sp41, sp51]
// sp44 shares 2 with sp54, sp53, [sp42, sp52]
// sp53 shares 2 with sp54, [sp44, sp52, sp41, sp43]
// sp54 shares 2 with [sp53, sp43, sp51, sp42, sp44]
// here the square brackets denote that the edges shared have already been specified above



void Clusters_Get12MW() {

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
                      for (int sp5a4_id = 0; sp5a4_id < nsp5a; ++sp5a4_id){
                        int *sp5a4 = hcsp5a[sp5a4_id];
                        if(shared_5a(sp5a4, sp5a1) == 2){ // then sp5a1 is OK
                          if(shared_4a(sp4a1, sp4a3) == 2 && shared_4a5a(sp4a1, sp5a2) == 2){
                            for (int sp5a3_id = 0; sp5a3_id < nsp5a; ++sp5a3_id){
                              int *sp5a3 = hcsp5a[sp5a3_id];
                              if(shared_4a5a(sp4a1, sp5a3) == 2){ // then sp4a1 is OK
                                if(shared_4a5a(sp4a2, sp5a4) == 2 && shared_4a5a(sp4a2, sp5a2) == 2){
                                  for (int sp4a4_id = 0; sp4a4_id < nsp4a; ++sp4a4_id){
                                    int *sp4a4 = hcsp4a[sp4a4_id];
                                    if(shared_4a(sp4a2, sp4a4) == 2){  // then sp4a2 is OK
                                      if(shared_4a5a(sp4a4, sp5a2) == 2){  // then sp5a2 is OK
                                        if(shared_4a5a(sp4a3, sp5a3) == 2 && shared_4a5a(sp4a3, sp5a4) == 2){  // then sp4a3 is OK
                                          if(shared_4a5a(sp4a4, sp5a3) == 2 && shared_4a5a(sp4a4, sp5a4) == 2){  // then sp4a4 is OK
                                            if(shared_5a(sp5a3, sp5a4) == 2){  // then sp5a3 is OK and sp5a4 is already dealt with
                                              add_12MW(sp5a2, sp5a4, sp4a1);
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


void add_12MW(const int *sp5a2, const int *sp5a4, const int *sp4a1){

// sp5a 2 and 4 give us 10
// so all their particles are counted
// and sp4a1 has two particlrs that are not bonded to sp5a2



    int clusSize = 12;

    if (n12MW == m12MW) {
        hc12MW = resize_2D_int(hc12MW, m12MW, m12MW + incrStatic, clusSize, -1);
        m12MW  = m12MW + incrStatic;
    }

    // add the third sp5a3 (one of the face sp5as)
    hc12MW[n12MW][0] = sp5a2[0];
    hc12MW[n12MW][1] = sp5a2[1];
    hc12MW[n12MW][2] = sp5a2[2];
    hc12MW[n12MW][3] = sp5a2[3];
    hc12MW[n12MW][4] = sp5a2[4];


    // add the fourth sp5a4 (the other of the face sp5as)
    hc12MW[n12MW][5] = sp5a4[0];
    hc12MW[n12MW][6] = sp5a4[1];
    hc12MW[n12MW][7] = sp5a4[2];
    hc12MW[n12MW][8] = sp5a4[3];
    hc12MW[n12MW][9] = sp5a4[4];


    // and now we need to find the particle in the first sp5a which is in neither of the other two
    int sp4a1extra=10;
    for(int i = 0; i < 4; ++i){
        if(is_particle_in_cluster(sp5a2, 5, sp4a1[i]) == 0) {
            hc12MW[n12MW][sp4a1extra] = sp4a1[i];
            ++sp4a1extra;
        }
    }

    //printf("\n add_12MW: sp4a1extra should be 12, it is %d\n", sp4a1extra);

    for (int i = 0; i < 12; ++i) {
        s12MW[hc12MW[n12MW][i]] = 'B';
    }
    ++n12MW;
}
