#include "math.h"
#include "bonds.h"
#include "globals.h"
#include "tools.h"

double Bonds_GetR2(int i, int j) {  // get separation between particles i and j
    double dx, dy, dz;
    double total_distance;

    dx = x[i] - x[j];
    dy = y[i] - y[j];
    dz = z[i] - z[j];
    total_distance = dx * dx + dy * dy + dz * dz;
    return total_distance;
}

double Bonds_GetR2_PBCs(int i, int j) { // get PBC wrapped separation between particles i and j
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

    else // if it is a triclinic periodic box...
    {
    if (dz > sidez*0.5) 
        {
            dz -= sidez;
            dy -= tiltyz;
            dx -= tiltxz;
        }  
    if (dz < -sidez*0.5)  {
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
    if (dx < -sidex*0.5)  {
        dx+= sidex;
    }
    
    return dx * dx + dy * dy + dz * dz;
        
    }

}

void Bonds_CheckSymmetric() {
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

void Bonds_GetBonds() {    // Get bonds using simple lengths
    int i, j, k;
    double dr2;

    if (Vor==1) {
        if (USELIST==0) Bonds_GetBondsV();
        else Bonds_GetBondsV_CellList();
        Bonds_CheckSymmetric();
        printf("Got Bonds\n");
        return;
    }
    
    printf("Simple: N%d rcut2_AA %.15lg rcutAB2 %.15lg rcutBB2 %.15lg\n",current_frame_particle_number,rcutAA2,rcutAB2,rcutBB2);
    
    if (PRINTINFO==1) { 
        printf("Simple Bond Length rcutAA %lg rcutAB %lg rcutBB %lg\n",rcutAA,rcutAB,rcutBB);
        if (PBCs==0) printf("No bonds through edge of box\n\n");
        else  printf("Periodic Boundary Conditions - PBC bonds\n\n");
    }
    for (i=0; i<current_frame_particle_number; ++i) cnb[i] = 0;
    // POSSIBLE IMPROVEMENT: add cell list here
    for (i=0; i<current_frame_particle_number; ++i) {
        for(j=i+1; j<current_frame_particle_number; ++j) {
            if (PBCs == 1) dr2 = Bonds_GetR2_PBCs(i,j);
            else dr2 = Bonds_GetR2(i,j);
            if (particle_type[i]==1 && particle_type[j]==1 && dr2 < rcutAA2){
                if (cnb[i] < nB && cnb[j] < nB){  // max number of bonds, do ith particle
                    k = cnb[i]++;
                    bNums[i][k] = j;
                    bondlengths[i][k]=sqrt(dr2);
                    k = cnb[j]++;
                    bNums[j][k] = i;
                    bondlengths[j][k]=sqrt(dr2);
                }
                else{    // list is now full
                    printf("Bonds_GetBonds(): nB %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",nB,i,j);
                    exit(1); 
                }
            }
            else if (particle_type[i]==2 && particle_type[j]==2 && dr2 < rcutBB2){
                if (cnb[i] < nB && cnb[j] < nB){  // max number of bonds, do ith particle
                    k = cnb[i]++;
                    bNums[i][k] = j;
                    bondlengths[i][k]=sqrt(dr2);
                    k = cnb[j]++;
                    bNums[j][k] = i;
                    bondlengths[j][k]=sqrt(dr2);
                }
                else{    // list is now full
                    printf("Bonds_GetBonds(): nB %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",nB,i,j);
                    exit(1); 
                }
            }
            else if (dr2 < rcutAB2) {
                if ((particle_type[i]==1 && particle_type[j]==2) || (particle_type[i]==2 && particle_type[j]==1)) {
                    if (cnb[i] < nB && cnb[j] < nB){  // max number of bonds, do ith particle
                        k = cnb[i]++;
                        bNums[i][k] = j;
                        bondlengths[i][k]=sqrt(dr2);
                        k = cnb[j]++;
                        bNums[j][k] = i;
                        bondlengths[j][k]=sqrt(dr2);
                    }
                    else{    // list is now full
                        printf("Bonds_GetBonds(): nB %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",nB,i,j);
                        exit(1); 
                    }
                }
            }
        }
        if (PRINTINFO==1) if (!((i+1)%1000)) printf("Bonds_GetBonds(): particle %d of %d done\n",i+1,current_frame_particle_number);
    }
    printf("\n");
    printf("Got Bonds\n");
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
            if (PBCs == 1) dr2 = Bonds_GetR2_PBCs(i,j);
            else dr2 = Bonds_GetR2(i,j);
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

void Bonds_GetBondsV_CellList() {  // Get bonds using Voronoi
    int i, j, k, l, m;
    int ic, jcell0, jcell,nabor;    // various counters
    const int nBs = 4 * nB;
    int cnbs, cnbs2;
    int *S, *S2, *Sb;
    double *Sr, *Sr2;
    double x1, x2, dr2;
    double rijx, rijy, rijz, rikx, riky, rikz, rjkx, rjky, rjkz;
    double *store_dr2;
    int *temp_cnb, **temp_bNums;
    char errMsg[1000];
    

    S = malloc(nBs*sizeof(int));
    S2 = malloc(nBs*sizeof(int));
    Sb = malloc(nBs*sizeof(int));
    Sr = malloc(nBs*sizeof(double));
    Sr2 = malloc(nBs*sizeof(double));

    store_dr2 = malloc(current_frame_particle_number*sizeof(double));   if (store_dr2==NULL) { sprintf(errMsg,"Bonds_GetBondsV_CellList(): store_dr2[] malloc out of memory\n");    Error(errMsg); }
    temp_cnb = malloc(current_frame_particle_number*sizeof(int));   if (temp_cnb==NULL) { sprintf(errMsg,"Bonds_GetBondsV_CellList(): temp_cnb[] malloc out of memory\n");  Error(errMsg); }
    temp_bNums = malloc(current_frame_particle_number*sizeof(int *));   if (temp_bNums==NULL) { sprintf(errMsg,"Bonds_GetBondsV_CellList(): temp_bNums[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<current_frame_particle_number; ++j) { temp_bNums[j] = malloc(nBs*sizeof(int));  if (temp_bNums[j]==NULL) { sprintf(errMsg,"Bonds_GetBondsV_CellList(): temp_bNums[][] malloc out of memory\n"); Error_no_free(errMsg); } }
            
    printf("Vor: N%d rcut2 %.15lg\n",current_frame_particle_number,rcutAA2);
   
    if (PRINTINFO==1) { 
        printf("Voronoi fc %lg rcutAA %lg\n",fc,rcutAA);
        if (PBCs==0) printf("No bonds through edge of box\n\n");
        else  printf("Periodic Boundary Conditions - PBC bonds\n\n");
    }
    llist[0]=-1;
    for (i=0; i<current_frame_particle_number; ++i) {
        llist[i+1]=-1;
        cnb[i] = 0;
        temp_cnb[i]=0;
        store_dr2[i]=-1.0;
        for (j=0; j<nBs; j++) temp_bNums[i][j]=0;
    }
    for (i=0; i<(n_cells_total+1); ++i) head[i]=-1;
    links();
    for (ic=1;ic<=n_cells_total;ic++) {        // loop over all cells
        i=head[ic];     // head of list particle for cell ic    
        while (i>0) {   // loop over all particles in ic
            
            j=llist[i]; // next particle in current cell ic
            while (j>0) {   // loop over all particles in cell ic
                if (PBCs == 1) dr2 = Bonds_GetR2_PBCs(i-1,j-1);
                else dr2 = Bonds_GetR2(i-1,j-1);
                if (dr2 < rcutAA2) {
                    if (temp_cnb[i-1] < nBs && temp_cnb[j-1] < nBs) {  // max number of bonds, do ith particle
                        temp_bNums[i-1][temp_cnb[i-1]]=j-1;
                        temp_bNums[j-1][temp_cnb[j-1]]=i-1;
                        temp_cnb[i-1]++;
                        temp_cnb[j-1]++;
                    }
                    else {    // list is now full
                        printf("Bonds_GetBondsV_CellList(): nBs %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",nBs,i-1,j-1);
                        exit(1); 
                    }
                }
                j=llist[j]; // loop over next particle in cell ic
            }
            jcell0=13*(ic-1);       // now loop over adjacent cells to cell ic
            for (nabor=1;nabor<=13;nabor++) {
                jcell=map[jcell0+nabor];    
                j=head[jcell];  // head of cell for jcell
                while (j>0) {   // loop over head of cell and all other particles in jcell
                    if (PBCs == 1) dr2 = Bonds_GetR2_PBCs(i-1,j-1);
                    else dr2 = Bonds_GetR2(i-1,j-1);
                    if (dr2 < rcutAA2) {
                        if (temp_cnb[i-1] < nBs && temp_cnb[j-1] < nBs) {  // max number of bonds, do ith particle
                            temp_bNums[i-1][temp_cnb[i-1]]=j-1;
                            temp_bNums[j-1][temp_cnb[j-1]]=i-1;
                            temp_cnb[i-1]++;
                            temp_cnb[j-1]++;
                        }
                        else {    // list is now full
                            printf("Bonds_GetBondsV_CellList(): nBs %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",nBs,i-1,j-1);
                            exit(1); 
                        }
                    }
                    j=llist[j]; // next particle in jcell
                }
            }
            i=llist[i]; // next particle in ic cell
        }
    }
    
    for (i=0; i<current_frame_particle_number; ++i) {
        cnbs = 0;
        for (j=0; j<current_frame_particle_number; ++j) {
            store_dr2[j]=-1.0;
        }
        for (j=0; j<temp_cnb[i]; ++j) {
            if (PBCs == 1) dr2 = Bonds_GetR2_PBCs(i,temp_bNums[i][j]);
            else dr2 = Bonds_GetR2(i,temp_bNums[i][j]);
            k = cnbs++;
            S[k] = temp_bNums[i][j];
            Sb[k] = 1;
            Sr[k] = dr2;
                    
            store_dr2[temp_bNums[i][j]]=dr2;
            
        } // We've now filled up the initial S
        cnbs2 = 0;
        for (j=0; j<cnbs; ++j) {
            for(k=0; k<cnbs2; ++k) { // find spot to insert S[j]
                if (Sr[j] < Sr2[k]) {
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
            printf("Bonds_GetBondsV_CellList(): part %d - cnbs %d does not equal cnbs2 %d \n",i,cnbs,cnbs2);
            exit(1); 
        }
        cnb[i]=0;
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

                if(box_type!=3){
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
                    if (rijz<-half_sidez) {
                        rijz +=sidez;
                        rijy +=tiltyz;
                        rijx +=tiltxz;
                    }
                    else if (rijz>half_sidez) {
                        rijz-=sidez;
                        rijy -=tiltyz;
                        rijx -=tiltxz;  
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

                    //k

                    if (rikz<-half_sidez) {
                        rikz+=sidez;
                        rikz +=tiltyz;
                        rikz +=tiltxz;

                    }
                    else if (rikz>half_sidez) {
                        rikz-=sidez;
                        rikz -=tiltyz;
                        rikz -=tiltxz;

                    }
                    if (riky<-half_sidey){
                            rikx+=tiltxy;
                            riky+=sidey;}
                    else if (riky>half_sidey) {
                        rikx-=tiltxy;
                        riky-=sidey;
                        }      
                    if (rikx<-half_sidex) rikx+=sidex;
                    else if (rikx>half_sidex) rikx-=sidex;

                    if (rjkz<-half_sidez) rjkz+=sidez;
                    else if (rjkz>half_sidez) rjkz-=sidez;
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
                    printf("Bonds_GetBondsV_CellList(): nB %d number of bonds per particle is not big enough: particle i %d cnb[i] %d or j %d cnb[j] %d has too many bonds\nThis is probably because rcutAA is too large\n",nB,i,cnb[i],j,cnb[j]);
                    exit(1); 
                }
            }
        }
        if (PRINTINFO==1) if (!((i+1)%10000)) printf("Bonds_GetBondsV_CellList(): particle %d of %d done\n",i+1,current_frame_particle_number);
    } // End i loop
    
    for (i=0; i<current_frame_particle_number; i++) free(temp_bNums[i]);
    free(temp_bNums);
    free(temp_cnb);
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
