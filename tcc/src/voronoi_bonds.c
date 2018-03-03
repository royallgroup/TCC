#include "globals.h"
#include "tools.h"
#include "bonds.h"
#include "math.h"

void Get_Bonds_With_Voronoi() {
    char error_message[200];
    int particle_1, particle_2, k, l, m;
    const int nBs = 4 * nB;
    int cnbs, cnbs2;
    int *S, *S2;
    double *Sr, *Sr2;
    double x1, x2, squared_distance;
    double rijx, rijy, rijz, rikx, riky, rikz, rjkx, rjky, rjkz;
    double *store_dr2;
    int *Sb;

    S = malloc(nBs*sizeof(int));
    S2 = malloc(nBs*sizeof(int));
    Sb = malloc(nBs*sizeof(int));
    Sr = malloc(nBs*sizeof(double));
    Sr2 = malloc(nBs*sizeof(double));
    store_dr2 = malloc(current_frame_particle_number*sizeof(double));

    printf("Vor: N%d rcut2 %.15lg\n",current_frame_particle_number,rcutAA2);

    for (particle_1=0; particle_1<current_frame_particle_number; ++particle_1) {
        cnb[particle_1] = 0;
    }

    for (particle_1=0; particle_1<current_frame_particle_number; ++particle_1) {
        cnbs = 0;
        for (particle_2=0; particle_2<current_frame_particle_number; ++particle_2) {
            if (particle_1 != particle_2) {
                squared_distance = Get_Interparticle_Distance(particle_1, particle_2);
                if (squared_distance < rcutAA2) {
                    if (cnbs < nBs) {
                        k = cnbs;
                        S[k] = particle_2;
                        Sb[k] = 1;
                        Sr[k] = squared_distance;
                        store_dr2[particle_2] = squared_distance;
                        cnbs++;
                    }
                    else {
                        Too_Many_Bonds(particle_1, particle_2, __func__);
                    }
                }
            }
        }
        cnbs2 = 0;
        for (particle_2=0; particle_2<cnbs; ++particle_2){
            for(k=0; k<cnbs2; ++k){ // find spot to insert S[particle_2]
                if (Sr[particle_2] < Sr2[k]){
                    for (l=cnbs2; l>k; --l) {
                        S2[l] = S2[l-1];
                        Sr2[l] = Sr2[l-1];
                    }
                    S2[k] = S[particle_2];
                    Sr2[k] = Sr[particle_2];
                    break;
                }
            }
            if (k==cnbs2){
                S2[cnbs2] = S[particle_2];
                Sr2[cnbs2] = Sr[particle_2];
            }
            ++cnbs2;
        } // Now sorted the list in order of distance from particle_1

        if (cnbs!=cnbs2) {
            sprintf(error_message, "Bonds_GetBondsV(): part %d - cnbs %d does not equal cnbs2 %d \n",particle_1,cnbs,cnbs2);
            Error(error_message);
        }

        for (particle_2=0; particle_2<cnbs2; ++particle_2) Sb[particle_2] = 1;

        for (l=0; l<cnbs2-1; ++l){
            k = S2[l];
            for (m=l+1; m<cnbs2; ++m) {
                particle_2 = S2[m];
                rijx = x[particle_1] - x[particle_2];
                rijy = y[particle_1] - y[particle_2];
                rijz = z[particle_1] - z[particle_2];
                rikx = x[particle_1] - x[k];
                riky = y[particle_1] - y[k];
                rikz = z[particle_1] - z[k];
                rjkx = x[particle_2] - x[k];
                rjky = y[particle_2] - y[k];
                rjkz = z[particle_2] - z[k];

                if (PBCs==1)  {
                    enforce_PBCs(&rijx, &rijy, &rijz);
                    enforce_PBCs(&rikx, &riky, &rikz);
                    enforce_PBCs(&rjkx, &rjky, &rjkz);
                }

                x1 = rijx * rikx + rijy * riky + rijz * rikz;
                x1 -= rijx * rjkx + rijy * rjky + rijz * rjkz;
                x2 = rikx * rikx + riky * riky + rikz * rikz;
                x2 += rjkx * rjkx + rjky * rjky + rjkz * rjkz;
                x1 = x1 / x2;
                if (x1-fc > EPS) { // Eliminate particle_2 from S
                    Sb[m] = 0;
                }
            }
        }

        for (l=0; l<cnbs2; ++l){
            particle_2 = S2[l];
            if (particle_type[particle_1]==2 && particle_type[particle_2]==2) {
                if (Sr2[l]>rcutBB2) {
                    Sb[l]=0;
                }
            }
            else if (particle_type[particle_1]==2 || particle_type[particle_2]==2) {
                if (Sr2[l]>rcutAB2) {
                    Sb[l]=0;
                }
            }
        }

        for (l=0; l<cnbs2; ++l) {
            if (Sb[l]) {
                particle_2 = S2[l];
                if (cnb[particle_1] < nB && cnb[particle_2] < nB) {  // max number of bonds, do ith particle
                    k = cnb[particle_1]++;
                    bNums[particle_1][k] = particle_2;
                    bondlengths[particle_1][k]=sqrt(store_dr2[particle_2]);
                }
                else {    // list is now full
                    Too_Many_Bonds(particle_1, particle_2, __func__);
                }
            }
        }
    }

    free(store_dr2);
    free(S);
    free(S2);
    free(Sb);
    free(Sr);
    free(Sr2);
}