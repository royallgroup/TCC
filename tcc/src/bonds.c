#include "bonds.h"
//// START: Bonds routines
double Bonds_GetR2(int i, int j) {  // get separation between particles i and j
    double dx, dy, dz;
    
    dx = x[i] - x[j];
    dy = y[i] - y[j];
    dz = z[i] - z[j];
    return dx * dx + dy * dy + dz * dz;
}

double Bonds_GetR2_PBCs(int i, int j) { // get PBC wrapped separation between particles i and j
    double dx, dy, dz;

    dx = x[i] - x[j];
    dy = y[i] - y[j];
    dz = z[i] - z[j];

    
    if(ISNOTCUBIC!=3)
    {
        if (dx<-halfSidex) dx+=sidex;
        else if (dx>halfSidex) dx-=sidex;
        if (dy<-halfSidey) dy+=sidey;
        else if (dy>halfSidey) dy-=sidey;
        if (dz<-halfSidez) dz+=sidez;
        else if (dz>halfSidez) dz-=sidez;
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




void Bonds_WriteBondNetwork(int f) {
    int i,j,k;
    char input[1000];
    FILE *writeout;
    
    sprintf(input,"bondnetwork.f%d.matrix_bonds",f);
    writeout=fopen(input, "w");
    for (i=0; i<N; ++i) {
        for (j=0; j<N; ++j) {
            if (Bonds_BondCheck(i,j)==1) k=1;
            else k=0;
            fprintf(writeout,"%d",k);
        }
        fprintf(writeout,"\n");
    }
    /*for (i=0; i<N; ++i) {
        for (j=0; j<N; ++j) {
            if (bondnetwork[i][j]!=bondnetwork[j][i]) {
                printf("Bonds_WriteBondNetwork(): i%d rtype %d j%d rtype %d bondnetwork matrix not symmetric! sep2 is %lg i-j%d j-i%d\n",i,rtype[i],j,rtype[j],Bonds_GetR2(i,j),bondnetwork[i][j],bondnetwork[j][i]);
                exit(1);
            }
        }
    }*/
    fclose(writeout);
    printf("d%d Bonds_WriteBondNetwork(): Written bondnetwork matrix to %s\n\n",rank,input);
}

void Bonds_WriteBonds(int f) {
    int i, j, sum;
    char errMsg[1000];
    
    sum=0;
    for (i=0; i<N; ++i) {
        sum+=cnb[i];
    }
    if (sum%2!=0) {
        sprintf(errMsg,"Bonds_WriteBonds(): total number of bonds is not even %d\n",sum);
        exit(1);
    }
    
    fprintf(bondsout,"frame %d  total bonds %d\n",f,sum/2);
    for (i=0; i<N; ++i) {
        fprintf(bondsout,"%d    %d",i,cnb[i]);
        for (j=0; j<cnb[i]; ++j) {
            fprintf(bondsout,"  %d  %.5lg",bNums[i][j],bondlengths[i][j]);
        }
        fprintf(bondsout,"\n");
    }
}

void Bonds_TickBLDistro(double length, int *histo, int *count) {
    int k;
    char errMsg[1000];
    
    k = (int)(length/binWidth);
    if (k>=BLDistroNoBins) {
        sprintf(errMsg,"Bonds_TickBLDistro():  interaction bondlength %lg binwidth %lg k %d nobins %d",length,binWidth,k,BLDistroNoBins);
        Error(errMsg);
    }
    histo[k]=histo[k]+1;
    (*count)=(*count)+1;
}

void Bonds_WriteBLDistro(char *filename, int *histo, int *count, double *meanRtn) {
    int i;
    double lengthtemp, mean;
    char errMsg[1000];
    FILE *writeout;
    
    writeout=fopen(filename,"w");
    if (writeout==NULL)  {
        sprintf(errMsg,"Bonds_WriteBLDistro(): Error opening file %s",filename);    // Always test file open
        Error(errMsg);
    }
    
    if ((*count)%2!=0) {
        sprintf(errMsg,"Bonds_WriteBLDistro(): total no samples %d is not even",*count);    // Always test file open
        Error(errMsg);
    }
    (*count)=(*count)/2;
    
    fprintf(writeout,"%s\n",filename);
    fprintf(writeout,"bond length   frequency   normalized frequency (%d samples)\n",(*count));
    
    lengthtemp=mean=0.0;
    for (i=0; i<BLDistroNoBins; i++) {
        if ((histo[i])%2!=0) {
            sprintf(errMsg,"Bonds_WriteBLDistro(): no samples %d is not even",histo[i]);    // Always test file open
            Error(errMsg);
        }
        histo[i]=histo[i]/2;
        
        fprintf(writeout,"%.15lg    %d  %.15lg\n",lengthtemp,histo[i],(double)(histo[i])/(*count));
        mean+=(lengthtemp+binWidth/2.0)*histo[i];
        lengthtemp+=binWidth;
    }
    (*meanRtn)=mean/(double)(*count);
    fprintf(writeout,"mean nearest neighbour bond length    %.15lg\n",meanBL);
    fclose(writeout);
    printf("d%d Written %s\n",rank,filename);
}

void Bonds_WriteBLDistroClust(char *filename, int *histo, int *count,double *meanRtn) {
    int i;
    double lengthtemp, mean;
    char errMsg[1000];
    FILE *writeout;
    
    writeout=fopen(filename,"w");
    if (writeout==NULL)  {
        sprintf(errMsg,"Bonds_WriteBLDistroClust(): Error opening file %s",filename);   // Always test file open
        Error(errMsg);
    }

    fprintf(writeout,"%s\n",filename);
    fprintf(writeout,"bond length   frequency   normalized frequency (%d samples)\n",(*count));
    
    lengthtemp=mean=0.0;
    for (i=0; i<BLDistroNoBins; i++) {
        fprintf(writeout,"%.15lg    %d  %.15lg\n",lengthtemp,histo[i],(double)(histo[i])/(*count));
        mean+=(lengthtemp+binWidth/2.0)*histo[i];
        lengthtemp+=binWidth;
    }
    (*meanRtn)=mean/(double)(*count);
    fprintf(writeout,"mean nearest neighbour bond length    %.15lg\n",(*meanRtn));
    fclose(writeout);
    printf("d%d Written %s\n",rank,filename);
}

void Bonds_TallynbDistro() {
    int i, j, noA;
    char errMsg[1000];
    
    for (i=0; i<N; i++) {
        if (cnb[i]>=0 && cnb[i]<=nB) {
            nbDistro[cnb[i]]++;
            nbDistroNoSamples++;
        }
        else {
            sprintf(errMsg,"Bonds_TallynbDistro(): cnb[%d] %d is less than 0 or greater than nB %d\n",i,cnb[i],nB);
            Error(errMsg);
        }
        
        if (doBinary==1) {
            noA=0;
            if (rtype[i]==1) {
                for (j=0; j<cnb[i]; j++) {
                    if (rtype[bNums[i][j]]==1) noA++;
                }
                if (noA>=0 && noA<=cnb[i]) {
                    nbDistroAA[noA]++;
                    nbDistroNoSamplesAA++;
                    nbDistroAB[cnb[i]-noA]++;
                    nbDistroNoSamplesAB++;
                }
                else {
                    sprintf(errMsg,"Bonds_TallynbDistro(): noA %d is less than 0 or greater than cnb[%d] %d \n",noA,i,cnb[i]);
                    Error(errMsg);
                }
            }
            else {
                for (j=0; j<cnb[i]; j++) {
                    if (rtype[bNums[i][j]]==1) noA++;
                }
                if (noA>=0 && noA<=cnb[i]) {
                    nbDistroBA[noA]++;
                    nbDistroNoSamplesBA++;
                    nbDistroBB[cnb[i]-noA]++;
                    nbDistroNoSamplesBB++;
                }
                else {
                    sprintf(errMsg,"Bonds_TallynbDistro(): noA %d is less than 0 or greater than cnb[%d] %d \n",noA,i,cnb[i]);  // Always test file open
                    Error(errMsg);
                }
            }
        }
    }
}

void Bonds_WritenbDistro(char *filename, int *histo, int *count, double *meannbRtn) {
    int i;
    int mean;
    char errMsg[1000];
    FILE *writeout;
    
    writeout=fopen(filename,"w");
    if (writeout==NULL)  {
        sprintf(errMsg,"Bonds_WritenbDistro(): Error opening file %s",filename);    // Always test file open
        Error(errMsg);
    }
    
    fprintf(writeout,"%s\n",filename);
    fprintf(writeout,"no Bonds  frequency   normalized frequency (%d samples)\n",(*count));
    
    for (i=0; i<=nB; i++) fprintf(writeout,"%d  %d  %.15lg\n",i,histo[i],(double)(histo[i])/(double)(*count));
    
    mean=0;
    for (i=0; i<=nB; i++) {
        mean+=i*histo[i];
    }
    (*meannbRtn)=(double)(mean)/(double)(*count);
    fprintf(writeout,"mean nB   %.15lg\n",*meannbRtn);
    fclose(writeout);
    printf("d%d Written %s\n",rank,filename);
}

void Bonds_CheckSymmetric() {
    int i, j, k;
    //char errMsg[1000];
    
    for (i=0; i<N; ++i) {
        for (j=0; j<cnb[i]; ++j) {
            for (k=0; k<cnb[bNums[i][j]]; k++) {
                if (i==bNums[bNums[i][j]][k]) break;
            }
            if (k==cnb[bNums[i][j]]) {
                //sprintf(errMsg,"Bonds_CheckSymmetric(): unsymmetric bond network from Voronoi code - i %d bonded to bNums[i][j] %d but bNums[i][j] %d not boned to i %d",i,bNums[i][j],bNums[i][j],i);
                //Error(errMsg);
                bNums[bNums[i][j]][k]=i;
                cnb[bNums[i][j]]++;
                bondlengths[bNums[i][j]][k]=bondlengths[i][j];
                correctedBonds++;
            }
        }
    }
}
void Bonds_GetBonds(int f) {    // Get bonds using simple lengths
    int i, j, k;
    double dr2;

    if (Vor==1) {
        if (USELIST==0) Bonds_GetBondsV(f);
        else Bonds_GetBondsV_CellList(f);
        Bonds_CheckSymmetric();
        //if (f==0) Bonds_WriteBondNetwork(f);
        //Bonds_ReadVoro(f);
        if (doBLDistros==1) {
            for (i=0; i<N; ++i) {
                for (j=0; j<cnb[i]; ++j) {
                    Bonds_TickBLDistro(bondlengths[i][j],BLDistro,&BLDistroNoSamples);
                    if (doBinary==1) {
                        if (rtype[i]==1 && rtype[bNums[i][j]]==1) Bonds_TickBLDistro(bondlengths[i][j],BLDistroAA,&BLDistroNoSamplesAA);
                        else if (rtype[i]==2 && rtype[bNums[i][j]]==2) Bonds_TickBLDistro(bondlengths[i][j],BLDistroBB,&BLDistroNoSamplesBB);
                        else Bonds_TickBLDistro(bondlengths[i][j],BLDistroAB,&BLDistroNoSamplesAB);
                    }
                }
            }
        }
        if (doWriteBonds==1) Bonds_WriteBonds(f);
        if (donbDistros==1) Bonds_TallynbDistro();
        printf("d%d Got Bonds\n",rank);
        return;
    }
    
    printf("d%d Simple: N%d NA%d rcut2_AA %.15lg rcutAB2 %.15lg rcutBB2 %.15lg\n",rank,N,NA,rcutAA2,rcutAB2,rcutBB2);
    
    if (PRINTINFO==1) { 
        printf("d%d Simple Bond Length rcutAA %lg rcutAB %lg rcutBB %lg\n",rank,rcutAA,rcutAB,rcutBB);
        if (PBCs==0) printf("d%d No bonds through edge of box\n\n",rank);
        else  printf("d%d Periodic Boundary Conditions - PBC bonds\n\n",rank);
    }
    for (i=0; i<N; ++i) cnb[i] = 0;
    // POSSIBLE IMPROVEMENT: add cell list here
    for (i=0; i<N; ++i) {
        for(j=i+1; j<N; ++j) {
            if (PBCs == 1) dr2 = Bonds_GetR2_PBCs(i,j);
            else dr2 = Bonds_GetR2(i,j);
            if (rtype[i]==1 && rtype[j]==1 && dr2 < rcutAA2){
                if (cnb[i] < nB && cnb[j] < nB){  // max number of bonds, do ith particle
                    k = cnb[i]++;
                    bNums[i][k] = j;
                    bondlengths[i][k]=sqrt(dr2);
                    k = cnb[j]++;
                    bNums[j][k] = i;
                    bondlengths[j][k]=sqrt(dr2);
                }
                else{    // list is now full
                    printf("d%d Bonds_GetBonds(): nB %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",rank,nB,i,j);
                    exit(1); 
                }
            }
            else if (rtype[i]==2 && rtype[j]==2 && dr2 < rcutBB2){
                if (cnb[i] < nB && cnb[j] < nB){  // max number of bonds, do ith particle
                    k = cnb[i]++;
                    bNums[i][k] = j;
                    bondlengths[i][k]=sqrt(dr2);
                    k = cnb[j]++;
                    bNums[j][k] = i;
                    bondlengths[j][k]=sqrt(dr2);
                }
                else{    // list is now full
                    printf("d%d Bonds_GetBonds(): nB %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",rank,nB,i,j);
                    exit(1); 
                }
            }
            else if (dr2 < rcutAB2) {
                if ((rtype[i]==1 && rtype[j]==2) || (rtype[i]==2 && rtype[j]==1)) {
                    if (cnb[i] < nB && cnb[j] < nB){  // max number of bonds, do ith particle
                        k = cnb[i]++;
                        bNums[i][k] = j;
                        bondlengths[i][k]=sqrt(dr2);
                        k = cnb[j]++;
                        bNums[j][k] = i;
                        bondlengths[j][k]=sqrt(dr2);
                    }
                    else{    // list is now full
                        printf("d%d Bonds_GetBonds(): nB %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",rank,nB,i,j);
                        exit(1); 
                    }
                }
            }
        }
        if (PRINTINFO==1) if (!((i+1)%1000)) printf("d%d Bonds_GetBonds(): particle %d of %d done\n",rank,i+1,N);
    }
    printf("\n");
    
    //if (f==0) Bonds_WriteBondNetwork(f);
    if (doBLDistros==1) {
        for (i=0; i<N; ++i) {
            for (j=0; j<cnb[i]; ++j) {
                Bonds_TickBLDistro(bondlengths[i][j],BLDistro,&BLDistroNoSamples);
                if (doBinary==1) {
                    if (rtype[i]==1 && rtype[bNums[i][j]]==1) Bonds_TickBLDistro(bondlengths[i][j],BLDistroAA,&BLDistroNoSamplesAA);
                    else if (rtype[i]==2 && rtype[bNums[i][j]]==2) Bonds_TickBLDistro(bondlengths[i][j],BLDistroBB,&BLDistroNoSamplesBB);
                    else Bonds_TickBLDistro(bondlengths[i][j],BLDistroAB,&BLDistroNoSamplesAB);
                }
            }
        }
    }
    if (doWriteBonds==1) Bonds_WriteBonds(f);
    if (donbDistros==1) Bonds_TallynbDistro();
    printf("d%d Got Bonds\n",rank);
}

void Bonds_GetBondsV(int f)  {  // Get bonds using Voronoi
    int i, j, k, l, m;
    const int nBs = 4 * nB;
    int cnbs, cnbs2;
    int S[nBs], S2[nBs];
    double Sr[nBs], Sr2[nBs];
    double x1, x2, dr2;
    double rijx, rijy, rijz, rikx, riky, rikz, rjkx, rjky, rjkz;
    double *store_dr2;
    int Sb[nBs];
    char errMsg[1000];
    
    store_dr2 = malloc(N*sizeof(double));   if (store_dr2==NULL) { sprintf(errMsg,"Bonds_GetBondsV(): store_dr2[] malloc out of memory\n"); Error(errMsg); }

    printf("d%d Vor: N%d NA%d rcut2 %.15lg\n",rank,N,NA,rcutAA2);
   
    if (PRINTINFO==1) { 
        printf("d%d Voronoi fc %lg rcutAA %lg\n",rank,fc,rcutAA);
        if (PBCs==0) printf("d%d No bonds through edge of box\n\n",rank);
        else  printf("d%d Periodic Boundary Conditions - PBC bonds\n\n",rank);
    }
    for (i=0; i<N; ++i) {
        cnb[i] = 0;
        store_dr2[i]=-1.0;
    }
    
    for (i=0; i<N; ++i) {
        cnbs = 0;
        for (j=0; j<N; ++j) {
            store_dr2[j]=-1.0;
        }
        for (j=0; j<N; ++j) {
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
                    printf("d%d Bonds_GetBondsV(): nBs %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",rank,nBs,i,j);
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
            printf("d%d Bonds_GetBondsV(): part %d - cnbs %d does not equal cnbs2 %d \n",rank,i,cnbs,cnbs2);
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
                if (ISNOTCUBIC!=3){
                    if (PBCs==1) { // if PBCs are being used
                        if (rijx>halfSidex) rijx-=sidex;
                        else if (rijx<-halfSidex) rijx+=sidex;
                        if (rijy>halfSidey) rijy-=sidey;
                        else if (rijy<-halfSidey) rijy+=sidey;
                        if (rijz>halfSidez) rijz-=sidez;
                        else if (rijz<-halfSidez) rijz+=sidez;

                        if (rikx>halfSidex) rikx-=sidex;
                        else if (rikx<-halfSidex) rikx+=sidex;
                        if (riky>halfSidey) riky-=sidey;
                        else if (riky<-halfSidey) riky+=sidey;
                        if (rikz>halfSidez) rikz-=sidez;
                        else if (rikz<-halfSidez) rikz+=sidez;

                        if (rjkx>halfSidex) rjkx-=sidex;
                        else if (rjkx<-halfSidex) rjkx+=sidex;
                        if (rjky>halfSidey) rjky-=sidey;
                        else if (rjky<-halfSidey) rjky+=sidey;
                        if (rjkz>halfSidez) rjkz-=sidez;
                        else if (rjkz<-halfSidez) rjkz+=sidez;
                    }
                }
                else {//if triclinc PBC are used
                    // printf("tilt  %g\n", tilt);
                    if (rijz<-halfSidez) {
                        rijz+=sidez;
                        rijy += tiltyz;
                        rijx += tiltxz;
                        }
                    else if (rijz>halfSidez) {
                        rijz-=sidez;
                        rijy -= tiltyz;
                        rijx -= tiltxz;
                    }
                    if (rijy<-halfSidey){   
                            rijx+=tiltxy;
                            rijy+=sidey;}

                    else if (rijy>halfSidey) {
                        rijx-=tiltxy;
                        rijy-=sidey;
                        }      
                    if (rijx<-halfSidex) rijx+=sidex;
                    else if (rijx>halfSidex) rijx-=sidex;

                    if (rikz<-halfSidez) {
                        rikz+=sidez;
                        riky += tiltyz;
                        rikx += tiltxz;

                    }
                    else if (rikz>halfSidez) rikz-=sidez;
                    if (riky<-halfSidey){   
                            rikx+=tiltxy;
                            riky+=sidey;}
                    else if (riky>halfSidey) {
                        rikx-=tiltxy;
                        riky-=sidey;
                        }      
                    if (rikx<-halfSidex) rikx+=sidex;
                    else if (rikx>halfSidex) rikx-=sidex;

                    if (rjkz<-halfSidez) {
                        rjkz+=sidez;
                        rjky += tiltyz;
                        rjkx += tiltxz;
                    }
                    else if (rjkz>halfSidez) {
                        rjkz-=sidez;
                        rjky -= tiltyz;
                        rjkx -= tiltxz;
                    }
                    if (rjky<-halfSidey){   
                            rjkx+=tiltxy;
                            rjky+=sidey;}
                    else if (rjky>halfSidey) {
                        rjkx-=tiltxy;
                        rjky-=sidey;
                        }      
                    if (rjkx<-halfSidex) rjkx+=sidex;
                    else if (rjkx>halfSidex) rjkx-=sidex;

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
            if (rtype[i]==2 && rtype[j]==2) {
                if (Sr2[l]>rcutBB2) {
                    Sb[l]=0;
                }
            }
            else if (rtype[i]==2 || rtype[j]==2) {
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
                    printf("d%d Bonds_GetBondsV(): nB %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",rank,nB,i,j);
                    exit(1); 
                }
            }
        }
        if (PRINTINFO==1) if (!((i+1)%1000)) printf("d%d Bonds_GetBondsV(): particle %d of %d done\n",rank,i+1,N);
    } // End i loop
    
    free(store_dr2);
}

void Bonds_GetBondsV_CellList(int f) {  // Get bonds using Voronoi
    int i, j, k, l, m;
    int ic, jcell0, jcell,nabor;    // various counters
    const int nBs = 4 * nB;
    int cnbs, cnbs2;
    int S[nBs], S2[nBs];
    double Sr[nBs], Sr2[nBs];
    double x1, x2, dr2;
    double rijx, rijy, rijz, rikx, riky, rikz, rjkx, rjky, rjkz;
    double *store_dr2;
    int *temp_cnb, **temp_bNums;
    int Sb[nBs];
    char errMsg[1000];
    
    cnbs=0;
    
    store_dr2 = malloc(N*sizeof(double));   if (store_dr2==NULL) { sprintf(errMsg,"Bonds_GetBondsV_CellList(): store_dr2[] malloc out of memory\n");    Error(errMsg); }
    temp_cnb = malloc(N*sizeof(int));   if (temp_cnb==NULL) { sprintf(errMsg,"Bonds_GetBondsV_CellList(): temp_cnb[] malloc out of memory\n");  Error(errMsg); }
    temp_bNums = malloc(N*sizeof(int *));   if (temp_bNums==NULL) { sprintf(errMsg,"Bonds_GetBondsV_CellList(): temp_bNums[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<N; ++j) { temp_bNums[j] = malloc(nBs*sizeof(int));  if (temp_bNums[j]==NULL) { sprintf(errMsg,"Bonds_GetBondsV_CellList(): temp_bNums[][] malloc out of memory\n"); Error_no_free(errMsg); } }
            
    printf("d%d Vor: N%d NA%d rcut2 %.15lg\n",rank,N,NA,rcutAA2);
   
    if (PRINTINFO==1) { 
        printf("d%d Voronoi fc %lg rcutAA %lg\n",rank,fc,rcutAA);
        if (PBCs==0) printf("d%d No bonds through edge of box\n\n",rank);
        else  printf("d%d Periodic Boundary Conditions - PBC bonds\n\n",rank);
    }
    llist[0]=-1;
    for (i=0; i<N; ++i) {
        llist[i+1]=-1;
        cnb[i] = 0;
        temp_cnb[i]=0;
        store_dr2[i]=-1.0;
        for (j=0; j<nBs; j++) temp_bNums[i][j]=0;
    }
    for (i=0; i<(ncells+1); ++i) head[i]=-1; 
    links();
    for (ic=1;ic<=ncells;ic++) {        // loop over all cells
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
                        printf("d%d Bonds_GetBondsV_CellList(): nBs %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",rank,nBs,i-1,j-1);
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
                            printf("d%d Bonds_GetBondsV_CellList(): nBs %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",rank,nBs,i-1,j-1);
                            exit(1); 
                        }
                    }
                    j=llist[j]; // next particle in jcell
                }
            }
            i=llist[i]; // next particle in ic cell
        }
    }
    
    for (i=0; i<N; ++i) {
        cnbs = 0;
        for (j=0; j<N; ++j) {
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
            printf("d%d Bonds_GetBondsV_CellList(): part %d - cnbs %d does not equal cnbs2 %d \n",rank,i,cnbs,cnbs2);
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

                if(ISNOTCUBIC!=3){
                if (PBCs==1) { // if PBCs are being used
                    if (rijx>halfSidex) rijx-=sidex;
                    else if (rijx<-halfSidex) rijx+=sidex;
                    if (rijy>halfSidey) rijy-=sidey;
                    else if (rijy<-halfSidey) rijy+=sidey;
                    if (rijz>halfSidez) rijz-=sidez;
                    else if (rijz<-halfSidez) rijz+=sidez;
                    if (rikx>halfSidex) rikx-=sidex;
                    else if (rikx<-halfSidex) rikx+=sidex;
                    if (riky>halfSidey) riky-=sidey;
                    else if (riky<-halfSidey) riky+=sidey;
                    if (rikz>halfSidez) rikz-=sidez;
                    else if (rikz<-halfSidez) rikz+=sidez;
                    if (rjkx>halfSidex) rjkx-=sidex;
                    else if (rjkx<-halfSidex) rjkx+=sidex;
                    if (rjky>halfSidey) rjky-=sidey;
                    else if (rjky<-halfSidey) rjky+=sidey;
                    if (rjkz>halfSidez) rjkz-=sidez;
                    else if (rjkz<-halfSidez) rjkz+=sidez;
                }
            }
                else {//if triclinc PBC are used
                    if (rijz<-halfSidez) {
                        rijz +=sidez;
                        rijy +=tiltyz;
                        rijx +=tiltxz;
                    }
                    else if (rijz>halfSidez) {
                        rijz-=sidez;
                        rijy -=tiltyz;
                        rijx -=tiltxz;  
                    }
                    if (rijy<-halfSidey){   
                            rijx+=tiltxy;
                            rijy+=sidey;}
                    else if (rijy>halfSidey) {
                        rijx-=tiltxy;
                        rijy-=sidey;
                        }

                    if (rijx<-halfSidex) rijx+=sidex;
                    else if (rijx>halfSidex) rijx-=sidex;

                    //k

                    if (rikz<-halfSidez) {
                        rikz+=sidez;
                        rikz +=tiltyz;
                        rikz +=tiltxz;

                    }
                    else if (rikz>halfSidez) {
                        rikz-=sidez;
                        rikz -=tiltyz;
                        rikz -=tiltxz;

                    }
                    if (riky<-halfSidey){   
                            rikx+=tiltxy;
                            riky+=sidey;}
                    else if (riky>halfSidey) {
                        rikx-=tiltxy;
                        riky-=sidey;
                        }      
                    if (rikx<-halfSidex) rikx+=sidex;
                    else if (rikx>halfSidex) rikx-=sidex;

                    if (rjkz<-halfSidez) rjkz+=sidez;
                    else if (rjkz>halfSidez) rjkz-=sidez;
                    if (rjky<-halfSidey){   
                            rjkx+=tiltxy;
                            rjky+=sidey;}
                    else if (rjky>halfSidey) {
                        rjkx-=tiltxy;
                        rjky-=sidey;
                        }      
                    if (rjkx<-halfSidex) rjkx+=sidex;
                    else if (rjkx>halfSidex) rjkx-=sidex;

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
            if (rtype[i]==2 && rtype[j]==2) {
                if (Sr2[l]>rcutBB2) {
                    Sb[l]=0;
                }
            }
            else if (rtype[i]==2 || rtype[j]==2) {
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
                    printf("d%d Bonds_GetBondsV_CellList(): nB %d number of bonds per particle is not big enough: particle i %d cnb[i] %d or j %d cnb[j] %d has too many bonds\nThis is probably because rcutAA is too large\n",rank,nB,i,cnb[i],j,cnb[j]);
                    exit(1); 
                }
            }
        }
        if (PRINTINFO==1) if (!((i+1)%10000)) printf("d%d Bonds_GetBondsV_CellList(): particle %d of %d done\n",rank,i+1,N);
    } // End i loop
    
    for (i=0; i<N; i++) free(temp_bNums[i]);
    free(temp_bNums);
    free(temp_cnb);
    free(store_dr2);
}

int Bonds_BondCheck(int i, int j) { // Returns 1 if i & j are bonded; 0 otherwise
    int k;

    for (k=0; k<cnb[i]; ++k) {
        if (bNums[i][k] == j) return 1;
    } 
    return 0;
}

int Bonds_cnb_j(int i, int j) { // Returns number k of if bNums[i][k]=j; -1 if i and j unbonded
    int k;

    for (k=0; k<cnb[i]; ++k) {
        if (bNums[i][k] == j) return k;
    }
    
    return -1;
}
