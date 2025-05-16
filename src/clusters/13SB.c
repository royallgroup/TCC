#include <globals.h>
#include <tools.h>
#include "13SB.h"


//!  A 11SB is a made of a 11SB cluster and two 5A clusters
/*!
* One spindle of each 11SB is a ring of the other 11SB
* There is one common ring particle
*/

void Clusters_Get13SB() {
    int spindle, spindle2, v,u,v2,u2;
    for (int old_11SB_id = 0; old_11SB_id < n11SB; ++old_11SB_id) {
        int *old_11SB_cluster = hc11SB[old_11SB_id];
        for (int first_5A_id = 0; first_5A_id < nsp3c; ++first_5A_id) {
            int *first_5A_cluster = hcsp3c[first_5A_id];
            if (first_5A_cluster[3] == old_11SB_cluster[6] ||first_5A_cluster[3] == old_11SB_cluster[7] ||
                first_5A_cluster[4] == old_11SB_cluster[6] ||first_5A_cluster[4] == old_11SB_cluster[7]){
                    
                if (overlap_13SB(old_11SB_cluster, first_5A_cluster) == 1 || overlap_13SB(old_11SB_cluster, first_5A_cluster) == 2){ // the ring particles of 5a have the right overlap
                    if (first_5A_cluster[3] == old_11SB_cluster[6] ||first_5A_cluster[3] == old_11SB_cluster[7]){
                        spindle = first_5A_cluster[4];
                        }
                    else if (first_5A_cluster[4] == old_11SB_cluster[6] ||first_5A_cluster[4] == old_11SB_cluster[7]){
                        spindle = first_5A_cluster[3];
                    }    
                                
                    for (int first_sp3a_id = 0; first_sp3a_id < nsp3a; ++first_sp3a_id) { //loop over 3 rings
                        int *first_sp3a_cluster = hcsp3a[first_sp3a_id];
                        
                        if(overlap_13SB_sp3a(old_11SB_cluster, first_sp3a_cluster,spindle) == 2){ //sp3a shares 2 particles with 11sb and spindle
                            v = 0;
                            u = 0;
                            
                            for(int i = 0; i < 3; ++i){
                                u2 = u;
                                v2 = u;
                                if(first_sp3a_cluster[i] == spindle){
                                    v += 1;
                                    v2 +=1;
                                }
                                if(first_sp3a_cluster[i] == old_11SB_cluster[9] ||first_sp3a_cluster[i] == old_11SB_cluster[10]){
                                    u += 1;
                                    u2 +=1;
                                } 
                                if(u-u2 == 0 && v-v2 ==0){
                                    spindle2 = first_sp3a_cluster[i];
                                }                           
                            }
                            
                            printf("%i %i\n",v,u);
                            if(v == 1 && u == 1){ 
                                
                                //if(check_5A_13SB(spindle, spindle2, old_11SB_cluster) == 1){ // check there is no 5a made up of spindle2 and the 11sb
                                    
                                    if(check_unique_13SB(old_11SB_cluster, spindle, spindle2) == 0){
                                        add_13SB(old_11SB_cluster, spindle, spindle2);
                                    }
                                //}
                            }
                        }
                    }
                }
            }
        }
    }
}

int check_5A_13SB(int spindle, int spindle2, int *clust_11SB){
    for (int first_5A_id = 0; first_5A_id < nsp3c; ++first_5A_id) {
        int u = 0;
        int v = 0;
        int *first_5A_cluster = hcsp3c[first_5A_id];
        for(int j = 0; j < 5; ++j){
            for(int i = 0; i < 11; ++i){
                if(first_5A_cluster[j] == clust_11SB[i]){
                    u +=1;
                }
            }            
        }
        for(int j = 0; j < 5; ++j){
            if(first_5A_cluster[j] == spindle){
                u += 1;
            }
        }
        for(int k = 3; k < 5; ++k){
            if(spindle2 == first_5A_cluster[k]){
                v +=1;
            }
        }
        if (v == 1 && u ==4){
            return 0;// have found a 5a cluster where the extra particle is a spindle and all the other particles are in the 12 cluster
        }
    }   
    return 1;
}

int overlap_13SB(int *clust_11SB, int *clust_5A){ // the ring particles in the 5a overlap in the right way with the 11sb
    int v = 0;
    int u = 0;
    int u2 = 0;
    int v2 = 0;
    for (int i = 0; i < 3; ++i){ //loop over ring
        if(clust_5A[i] == clust_11SB[0]){
            u +=1;
        }
        if(clust_5A[i] == clust_11SB[1]){
            u +=1;
        }
        if(clust_5A[i] == clust_11SB[9]){
            u2 +=1;
        }

        if(clust_5A[i] == clust_11SB[3]){
            v +=1;
        }
        if(clust_5A[i] == clust_11SB[2]){
            v +=1;
        }
        if(clust_5A[i] == clust_11SB[10]){
            v2 +=1;
        }
    }
    
    if(v == 2 && (v2 == 1 || u2 ==1)){
        return 2;
    }
    if(u == 2 && (v2 == 1 || u2 ==1)){
        return 1;
    }
    return 0;
}

int overlap_13SB_sp3a(int *clust_11SB, int *clust_sp3a,int spindle){
    int v = 0;
    int u = 0;
    for (int i = 0; i < 11; ++i){ //loop over ring
        for (int j = 0; j < 3; ++j){ //loop over ring
            if(clust_11SB[i] == clust_sp3a[j]){
                u +=1;
            }
        }
    }
    for (int k = 0; k < 3; ++k){ //loop over ring
        if(spindle == clust_sp3a[k]){
            v +=1;
           
        }
    }
    printf("%i %i\n",u,v);
    if(u == 1 && v == 1){
        return 2;
    }
    return 0;
}

int check_unique_13SB(int *clust11SB, int s1, int s2){
    int u;
    for (int old_13SB_id = 0; old_13SB_id < n13SB; ++old_13SB_id) {
        u = 0;
        for (int q = 0; q < 11; ++q){
            for (int r = 0; r < 13; ++r){
                if(hc13SB[old_13SB_id][r] == clust11SB[q]){
                    u += 1;
                }
            }

        }
            for (int r = 0; r < 13; ++r){
                if(hc13SB[old_13SB_id][r] == s1 || hc13SB[old_13SB_id][r] == s2){
                    u += 1;
                }
            }
    if(u == 13){
        return 1;            
    }
    }
    return 0;
}

void add_13SB(int *clust11SB, int spindle1, int spindle2) {
    int clusSize = 13;
    //printf("new_particle %i\n", new_particle);
    if (n13SB == m13SB) {
        hc13SB = resize_2D_int(hc13SB, m13SB, m13SB + incrStatic, clusSize, -1);
        m13SB = m13SB + incrStatic;
    }
    hc13SB[n13SB][0] = clust11SB[0];
    hc11SB[n13SB][1] = clust11SB[1];
    hc11SB[n13SB][2] = clust11SB[2];
    hc11SB[n13SB][3] = clust11SB[3];
    hc11SB[n13SB][4] = clust11SB[4];
    hc11SB[n13SB][5] = clust11SB[5];
    hc11SB[n13SB][6] = clust11SB[6];
    hc11SB[n13SB][7] = clust11SB[7];
    hc11SB[n13SB][8] = clust11SB[8];
    hc11SB[n13SB][9] = clust11SB[8];
    hc11SB[n13SB][10] = clust11SB[8];
    hc11SB[n13SB][11] = spindle1;
    hc11SB[n13SB][12] = spindle2;

    for (int i = 0; i < 13; ++i) {
        s13SB[hc13SB[n13SB][i]] = 'B';
    }
    ++n13SB;
}