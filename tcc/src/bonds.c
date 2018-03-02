#include "math.h"
#include "bonds.h"
#include "globals.h"
#include "tools.h"
#include "cell_list.h"
#include "string.h"

double Get_Interparticle_Distance(int i, int j) {
    // Returns the squared interparticle distance between i and j
    double dx, dy, dz;
    double total_distance;

    dx = x[i] - x[j];
    dy = y[i] - y[j];
    dz = z[i] - z[j];
    total_distance = dx * dx + dy * dy + dz * dz;
    return total_distance;
}

double Get_Interparticle_Distance_With_PBCs(int i, int j) {
    // Returns the PBC wrapped squared interparticle distance between i and j
    double dx, dy, dz;

    dx = x[i] - x[j];
    dy = y[i] - y[j];
    dz = z[i] - z[j];

    if(box_type==2)
    {
        if (dx<-half_sidex) dx+=sidex;
        else if (dx>half_sidex) dx-=sidex;
        if (dy<-half_sidey) dy+=sidey;
        else if (dy>half_sidey) dy-=sidey;
        if (dz<-half_sidez) dz+=sidez;
        else if (dz>half_sidez) dz-=sidez;
        return dx * dx + dy * dy + dz * dz;
    }

    else {
        // if it is a triclinic periodic box...
        if (dz > sidez*0.5) {
                dz -= sidez;
                dy -= tiltyz;
                dx -= tiltxz;
        }
        if (dz < -sidez*0.5) {
            dz += sidez;
            dy += tiltyz;
            dx += tiltxz;
        }
            //deal with y, which affects x
        if (dy > sidey*0.5) {
                dx-=tiltxy;
                dy -= sidey;
        }
        if (dy < -sidey*0.5) {
                dx+=tiltxy ;
                dy += sidey;
        }
            //deal with x
        if (dx > sidex*0.5) {
            dx -= sidex;
        }
        if (dx < -sidex*0.5) {
            dx+= sidex;
        }
        return dx * dx + dy * dy + dz * dz;
    }

}

void Are_All_Bonds_Symmetric() {
    int i, j, k;

    for (i=0; i<current_frame_particle_number; ++i) {
        for (j=0; j<cnb[i]; ++j) {
            for (k=0; k<cnb[bNums[i][j]]; k++) {
                if (i==bNums[bNums[i][j]][k]) break;
            }
            if (k==cnb[bNums[i][j]]) {
                bNums[bNums[i][j]][k]=i;
                cnb[bNums[i][j]]++;
                bondlengths[bNums[i][j]][k]=bondlengths[i][j];
                correctedBonds++;
            }
        }
    }
}

void Get_Bonds() {

    if (Vor==1) {
        if (USELIST==0) Bonds_GetBondsV();
        else Bonds_GetBondsV_CellList();
        Are_All_Bonds_Symmetric();
    }
    else {
        Get_Simple_Bonds();
    }

    printf("\n");
    printf("Got Bonds\n");
}

void Get_Simple_Bonds() {
    // Get bonds using simple lengths
    int particle_1, particle_2, k;
    double squared_distance;

    printf("Simple: N%d rcut2_AA %.15lg rcutAB2 %.15lg rcutBB2 %.15lg\n",current_frame_particle_number,rcutAA2,rcutAB2,rcutBB2);

    memset(cnb, 0, current_frame_particle_number* sizeof(int));

    for (particle_1=0; particle_1<current_frame_particle_number; ++particle_1) {
        for(particle_2=particle_1+1; particle_2<current_frame_particle_number; ++particle_2) {
            if (PBCs == 1) {
                squared_distance = Get_Interparticle_Distance_With_PBCs(particle_1, particle_2);
            }
            else {
                squared_distance = Get_Interparticle_Distance(particle_1, particle_2);
            }
            if (particle_type[particle_1]==1 && particle_type[particle_2]==1 && squared_distance < rcutAA2){
                if (cnb[particle_1] < nB && cnb[particle_2] < nB){  // max number of bonds, do ith particle
                    k = cnb[particle_1]++;
                    bNums[particle_1][k] = particle_2;
                    bondlengths[particle_1][k]=sqrt(squared_distance);
                    k = cnb[particle_2]++;
                    bNums[particle_2][k] = particle_1;
                    bondlengths[particle_2][k]=sqrt(squared_distance);
                }
                else{    // list is now full
                    printf("Get_Bonds(): nB %d number of bonds per particle is not big enough: particle particle_1 %d or particle_2% d has too many bonds\nThis is probably because rcutAA is too large\n",nB,particle_1,particle_2);
                    exit(1);
                }
            }
            else if (particle_type[particle_1]==2 && particle_type[particle_2]==2 && squared_distance < rcutBB2){
                if (cnb[particle_1] < nB && cnb[particle_2] < nB){  // max number of bonds, do ith particle
                    k = cnb[particle_1]++;
                    bNums[particle_1][k] = particle_2;
                    bondlengths[particle_1][k]=sqrt(squared_distance);
                    k = cnb[particle_2]++;
                    bNums[particle_2][k] = particle_1;
                    bondlengths[particle_2][k]=sqrt(squared_distance);
                }
                else{    // list is now full
                    printf("Get_Bonds(): nB %d number of bonds per particle is not big enough: particle particle_1 %d or particle_2% d has too many bonds\nThis is probably because rcutAA is too large\n",nB,particle_1,particle_2);
                    exit(1);
                }
            }
            else if (squared_distance < rcutAB2) {
                if ((particle_type[particle_1]==1 && particle_type[particle_2]==2) || (particle_type[particle_1]==2 && particle_type[particle_2]==1)) {
                    if (cnb[particle_1] < nB && cnb[particle_2] < nB){  // max number of bonds, do ith particle
                        k = cnb[particle_1]++;
                        bNums[particle_1][k] = particle_2;
                        bondlengths[particle_1][k]=sqrt(squared_distance);
                        k = cnb[particle_2]++;
                        bNums[particle_2][k] = particle_1;
                        bondlengths[particle_2][k]=sqrt(squared_distance);
                    }
                    else{    // list is now full
                        printf("Get_Bonds(): nB %d number of bonds per particle is not big enough: particle particle_1 %d or particle_2% d has too many bonds\nThis is probably because rcutAA is too large\n",nB,particle_1,particle_2);
                        exit(1);
                    }
                }
            }
        }
        if (PRINTINFO==1) if (!((particle_1+1)%1000)) printf("Get_Bonds(): particle %d of %d done\n",particle_1+1,current_frame_particle_number);
    }
}

void Bonds_GetBondsV()  {  // Get bonds using Voronoi
    int i, j, k, l, m;
    const int nBs = 4 * nB;
    int cnbs, cnbs2;
    int *S, *S2;
    double *Sr, *Sr2;
    double x1, x2, dr2;
    double rijx, rijy, rijz, rikx, riky, rikz, rjkx, rjky, rjkz;
    double *store_dr2;
    int *Sb;
    char errMsg[1000];

    S = malloc(nBs*sizeof(int));
    S2 = malloc(nBs*sizeof(int));
    Sb = malloc(nBs*sizeof(int));
    Sr = malloc(nBs*sizeof(double));
    Sr2 = malloc(nBs*sizeof(double));
    store_dr2 = malloc(current_frame_particle_number*sizeof(double));   if (store_dr2==NULL) { sprintf(errMsg,"Bonds_GetBondsV(): store_dr2[] malloc out of memory\n"); Error(errMsg); }

    printf("Vor: N%d rcut2 %.15lg\n",current_frame_particle_number,rcutAA2);
   
    if (PRINTINFO==1) { 
        printf("Voronoi fc %lg rcutAA %lg\n",fc,rcutAA);
        if (PBCs==0) printf("No bonds through edge of box\n\n");
        else  printf("Periodic Boundary Conditions - PBC bonds\n\n");
    }
    for (i=0; i<current_frame_particle_number; ++i) {
        cnb[i] = 0;
        store_dr2[i]=-1.0;
    }
    
    for (i=0; i<current_frame_particle_number; ++i) {
        cnbs = 0;
        for (j=0; j<current_frame_particle_number; ++j) {
            store_dr2[j]=-1.0;
        }
        for (j=0; j<current_frame_particle_number; ++j) {
            if (i==j) continue;
            if (PBCs == 1) dr2 = Get_Interparticle_Distance_With_PBCs(i, j);
            else dr2 = Get_Interparticle_Distance(i, j);
            if (dr2 < rcutAA2) {
                if (cnbs < nBs) {  // max number of bonds, do ith particle
                    k = cnbs++;
                    S[k] = j;
                    Sb[k] = 1;
                    Sr[k] = dr2;
                    
                    store_dr2[j]=dr2;
                }
                else {    // list is now full
                    printf("Bonds_GetBondsV(): nBs %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",nBs,i,j);
                    exit(1); 
                }
            }
        } // We've now filled up the initial S
        cnbs2 = 0;
        for (j=0; j<cnbs; ++j){
            for(k=0; k<cnbs2; ++k){ // find spot to insert S[j]
                if (Sr[j] < Sr2[k]){
                    for (l=cnbs2; l>k; --l) {
                        S2[l] = S2[l-1];
                        Sr2[l] = Sr2[l-1];
                    }
                    S2[k] = S[j];
                    Sr2[k] = Sr[j];
                    break;
                }
            }
            if (k==cnbs2){
                S2[cnbs2] = S[j];
                Sr2[cnbs2] = Sr[j];
            }
            ++cnbs2;  
        } // Now sorted the list in order of distance from i
        
        if (cnbs!=cnbs2) {
            printf("Bonds_GetBondsV(): part %d - cnbs %d does not equal cnbs2 %d \n",i,cnbs,cnbs2);
            exit(1); 
        }
            
        for (j=0; j<cnbs2; ++j) Sb[j] = 1;
        
        for (l=0; l<cnbs2-1; ++l){       
            k = S2[l];
            for (m=l+1; m<cnbs2; ++m) {
                j = S2[m];
                rijx = x[i] - x[j];
                rijy = y[i] - y[j];
                rijz = z[i] - z[j];
                rikx = x[i] - x[k];
                riky = y[i] - y[k];
                rikz = z[i] - z[k];
                rjkx = x[j] - x[k];
                rjky = y[j] - y[k];
                rjkz = z[j] - z[k];
                if (box_type!=3){
                    if (PBCs==1) { // if PBCs are being used
                        if (rijx>half_sidex) rijx-=sidex;
                        else if (rijx<-half_sidex) rijx+=sidex;
                        if (rijy>half_sidey) rijy-=sidey;
                        else if (rijy<-half_sidey) rijy+=sidey;
                        if (rijz>half_sidez) rijz-=sidez;
                        else if (rijz<-half_sidez) rijz+=sidez;

                        if (rikx>half_sidex) rikx-=sidex;
                        else if (rikx<-half_sidex) rikx+=sidex;
                        if (riky>half_sidey) riky-=sidey;
                        else if (riky<-half_sidey) riky+=sidey;
                        if (rikz>half_sidez) rikz-=sidez;
                        else if (rikz<-half_sidez) rikz+=sidez;

                        if (rjkx>half_sidex) rjkx-=sidex;
                        else if (rjkx<-half_sidex) rjkx+=sidex;
                        if (rjky>half_sidey) rjky-=sidey;
                        else if (rjky<-half_sidey) rjky+=sidey;
                        if (rjkz>half_sidez) rjkz-=sidez;
                        else if (rjkz<-half_sidez) rjkz+=sidez;
                    }
                }
                else {//if triclinc PBC are used
                    // printf("tilt  %g\n", tilt);
                    if (rijz<-half_sidez) {
                        rijz+=sidez;
                        rijy += tiltyz;
                        rijx += tiltxz;
                        }
                    else if (rijz>half_sidez) {
                        rijz-=sidez;
                        rijy -= tiltyz;
                        rijx -= tiltxz;
                    }
                    if (rijy<-half_sidey){
                            rijx+=tiltxy;
                            rijy+=sidey;}

                    else if (rijy>half_sidey) {
                        rijx-=tiltxy;
                        rijy-=sidey;
                        }      
                    if (rijx<-half_sidex) rijx+=sidex;
                    else if (rijx>half_sidex) rijx-=sidex;

                    if (rikz<-half_sidez) {
                        rikz+=sidez;
                        riky += tiltyz;
                        rikx += tiltxz;

                    }
                    else if (rikz>half_sidez) rikz-=sidez;
                    if (riky<-half_sidey){
                            rikx+=tiltxy;
                            riky+=sidey;}
                    else if (riky>half_sidey) {
                        rikx-=tiltxy;
                        riky-=sidey;
                        }      
                    if (rikx<-half_sidex) rikx+=sidex;
                    else if (rikx>half_sidex) rikx-=sidex;

                    if (rjkz<-half_sidez) {
                        rjkz+=sidez;
                        rjky += tiltyz;
                        rjkx += tiltxz;
                    }
                    else if (rjkz>half_sidez) {
                        rjkz-=sidez;
                        rjky -= tiltyz;
                        rjkx -= tiltxz;
                    }
                    if (rjky<-half_sidey){
                            rjkx+=tiltxy;
                            rjky+=sidey;}
                    else if (rjky>half_sidey) {
                        rjkx-=tiltxy;
                        rjky-=sidey;
                        }      
                    if (rjkx<-half_sidex) rjkx+=sidex;
                    else if (rjkx>half_sidex) rjkx-=sidex;

                }

                x1 = rijx * rikx + rijy * riky + rijz * rikz;
                x1 -= rijx * rjkx + rijy * rjky + rijz * rjkz;  
                x2 = rikx * rikx + riky * riky + rikz * rikz;
                x2 += rjkx * rjkx + rjky * rjky + rjkz * rjkz;
                x1 = x1 / x2;
                if (x1-fc > EPS) { // Eliminate j from S
                    Sb[m] = 0;
                }
            }
        }
        
        for (l=0; l<cnbs2; ++l){ 
            j = S2[l];
            if (particle_type[i]==2 && particle_type[j]==2) {
                if (Sr2[l]>rcutBB2) {
                    Sb[l]=0;
                }
            }
            else if (particle_type[i]==2 || particle_type[j]==2) {
                if (Sr2[l]>rcutAB2) {
                    Sb[l]=0;
                }
            }
        }
        
        for (l=0; l<cnbs2; ++l) {
            if (Sb[l]) {
                j = S2[l];
                if (cnb[i] < nB && cnb[j] < nB) {  // max number of bonds, do ith particle
                    k = cnb[i]++;
                    bNums[i][k] = j;
                    bondlengths[i][k]=sqrt(store_dr2[j]);
                }
                else {    // list is now full
                    printf("Bonds_GetBondsV(): nB %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",nB,i,j);
                    exit(1); 
                }
            }
        }
        if (PRINTINFO==1) if (!((i+1)%1000)) printf("Bonds_GetBondsV(): particle %d of %d done\n",i+1,current_frame_particle_number);
    }
    
    free(store_dr2);
    free(S);
    free(S2);
    free(Sb);
    free(Sr);
    free(Sr2);
}



int Bonds_BondCheck(int i, int j) { // Returns 1 if i & j are bonded; 0 otherwise
    int k;

    for (k=0; k<cnb[i]; ++k) {
        if (bNums[i][k] == j) return 1;
    } 
    return 0;
}
