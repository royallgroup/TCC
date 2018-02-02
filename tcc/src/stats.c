#include "stats.h"

void Stats_Init() {
    int i;
    char errMsg[1000];

    a5=malloc(N*sizeof(int));	if (a5==NULL) { sprintf(errMsg,"Stats_Init(): a5[] malloc out of memory\n");	Error(errMsg); }
    a6=malloc(N*sizeof(int));	if (a6==NULL) { sprintf(errMsg,"Stats_Init(): a6[] malloc out of memory\n");	Error(errMsg); }
    a7=malloc(N*sizeof(int));	if (a7==NULL) { sprintf(errMsg,"Stats_Init(): a7[] malloc out of memory\n");	Error(errMsg); }
    a8=malloc(N*sizeof(int));	if (a8==NULL) { sprintf(errMsg,"Stats_Init(): a8[] malloc out of memory\n");	Error(errMsg); }
    a9=malloc(N*sizeof(int));	if (a9==NULL) { sprintf(errMsg,"Stats_Init(): a9[] malloc out of memory\n");	Error(errMsg); }
    a10=malloc(N*sizeof(int));	if (a10==NULL) { sprintf(errMsg,"Stats_Init(): a10[] malloc out of memory\n");	Error(errMsg); }
    a11=malloc(N*sizeof(int));	if (a11==NULL) { sprintf(errMsg,"Stats_Init(): a11[] malloc out of memory\n");	Error(errMsg); }
    a12=malloc(N*sizeof(int));	if (a12==NULL) { sprintf(errMsg,"Stats_Init(): a12[] malloc out of memory\n");	Error(errMsg); }
    a13=malloc(N*sizeof(int));	if (a13==NULL) { sprintf(errMsg,"Stats_Init(): a13[] malloc out of memory\n");	Error(errMsg); }
    a15=malloc(N*sizeof(int));	if (a15==NULL) { sprintf(errMsg,"Stats_Init(): a15[] malloc out of memory\n");	Error(errMsg); }

    for(i=0;i<N;i++) {
        a5[i]=a6[i]=a7[i]=a8[i]=a9[i]=a10[i]=a11[i]=a12[i]=a13[i]=a15[i]=0;
    }

    ngsp3=ngsp3a=ngsp3b=ng5A=0;
    ngsp4=ngsp4a=ngsp4b=ng6A=0;
    ngsp5=ngsp5a=ngsp5b=ng7A=0;
    ng6Z=ng7K=0;
    ng8A=ng8B=ng8K=0;
    ng9A=ng9B=ng9K=0;
    ng10A=ng10B=ng10K=ng10W=0;
    ng11A=ng11B=ng11C=ng11E=ng11F=ng11W=0;	// number of particles in clusers gross
    ng12A=ng12B=ng12D=ng12E=ng12K=0;
    ng13A=ng13B=ng13K=0;
    ngFCC=ngHCP=ngBCC_9=ngBCC_15=0;

    nn5A=0;
    nn6A=0;
    nn7A=0;
    nn6Z=nn7K=0;
    nn8A=nn8B=nn8K=0;
    nn9A=nn9B=nn9K=0;
    nn10A=nn10B=nn10K=nn10W=0;
    nn11A=nn11B=nn11C=nn11E=nn11F=nn11W=0;	// number of particles in clusers gross
    nn12A=nn12B=nn12D=nn12E=nn12K=0;
    nn13A=nn13B=nn13K=0;
    nnFCC=nnHCP=nnBCC_9=nnBCC_15=0;

    ncsp3=ncsp3a=ncsp3b=nc5A=0;
    ncsp4=ncsp4a=ncsp4b=nc6A=0;
    ncsp5=ncsp5a=ncsp5b=nc7A=0;
    nc6Z=nc7K=0;
    nc8A=nc8B=nc8K=0;
    nc9A=nc9B=nc9K=0;
    nc10A=nc10B=nc10K=nc10W=0;
    nc11A=nc11B=nc11C=nc11E=nc11F=nc11W=0;	// number of particles in clusers gross
    nc12A=nc12B=nc12D=nc12E=nc12K=0;
    nc13A=nc13B=nc13K=0;
    ncFCC=ncHCP=ncBCC_9=ncBCC_15=0;

    ncsp3_excess_spindles=ncsp4_excess_spindles=ncsp5_excess_spindles=0;	// total number of _excess_spindlesed basic clusters
    ncsp3c_spindlebonds=ncsp4c_spindlebonds=nc6A_spindlebonds=ncsp5c_spindlebonds=0;
}

void Stats_Reset() {
    int i;

    for(i=0;i<N;i++) {
        a5[i]=a6[i]=a7[i]=a8[i]=a9[i]=a10[i]=a11[i]=a12[i]=a13[i]=a15[i]=0;
    }
}

void Stats_FreeMem() {
    free(a5);
    free(a6);
    free(a7);
    free(a8);
    free(a9);
    free(a10);
    free(a11);
    free(a12);
    free(a13);
    free(a15);
}

void Stats_SetA() { // Set arrays to true if the ith particle is a member of any clusters with this or a larger number of particles
    int flg1, flg2;
    int i;

    for (i=0;i<N;i++) {
        flg1 = sFCC[i] != 'C' || sHCP[i] != 'C' || sBCC_15[i] != 'C' || s13A[i] != 'C' || s13B[i] != 'C' || s13K[i] != 'C';
        if(flg1==1) a13[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = s12A[i] != 'C' || s12B[i] != 'C' || s12D[i] != 'C' || s12E[i] != 'C' || s12K[i] != 'C';
        flg2 = a13[i];
        if(flg1==1 || flg2==1) a12[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = s11A[i] != 'C' || s11B[i] != 'C' || s11C[i] != 'C' || s11E[i] != 'C' || s11F[i] != 'C' || s11W[i] != 'C';
        flg2 = a12[i];
        if(flg1==1 || flg2==1) a11[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = s10A[i] != 'C' || s10B[i] != 'C' || s10K[i] != 'C' || s10W[i] != 'C';
        flg2 = a11[i];
        if(flg1==1 || flg2==1) a10[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = s9A[i] != 'C' || s9B[i] != 'C' || s9K[i] != 'C' || sBCC_9[i] != 'C';
        flg2 = a10[i];
        if(flg1==1 || flg2==1) a9[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = s8A[i] != 'C' || s8B[i] != 'C' || s8K[i] != 'C';
        flg2 = a9[i];
        if(flg1==1 || flg2==1) a8[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = s7A[i] != 'C' || s7K[i] != 'C';
        flg2 = a8[i];
        if(flg1==1 || flg2==1) a7[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = s6A[i] != 'C' || s6Z[i] != 'C';
        flg2 = a7[i];
        if(flg1==1 || flg2==1) a6[i] = 1;
    }
}

void Stats_Analyse() {
    int i, Nclus;
    int flg;

    for(i=0; i<N; ++i){
        if(ssp3[i] != 'C') ++ngsp3;
        if(ssp3a[i] != 'C') ++ngsp3a;
        if(ssp3b[i] != 'C') ++ngsp3b;
        if(s5A[i] != 'C') ++ng5A;
        if(ssp4[i] != 'C') ++ngsp4;
        if(ssp4a[i] != 'C') ++ngsp4a;
        if(ssp4b[i] != 'C') ++ngsp4b;
        if(s6A[i] != 'C') ++ng6A;
        if(s6Z[i] != 'C') ++ng6Z;
        if(s7K[i] != 'C') ++ng7K;
        if(ssp5[i] != 'C') ++ngsp5;
        if(ssp5a[i] != 'C') ++ngsp5a;
        if(ssp5b[i] != 'C') ++ngsp5b;
        if(s7A[i] != 'C') ++ng7A;
        if(s8A[i] != 'C') ++ng8A;
        if(s8B[i] != 'C') ++ng8B;
        if(s8K[i] != 'C') ++ng8K;
        if(s9A[i] != 'C') ++ng9A;
        if(s9B[i] != 'C') ++ng9B;
        if(s9K[i] != 'C') ++ng9K;
        if(s10A[i] != 'C') ++ng10A;
        if(s10B[i] != 'C') ++ng10B;
        if(s10K[i] != 'C') ++ng10K;
        if(s10W[i] != 'C') ++ng10W;
        if(s11A[i] != 'C') ++ng11A;
        if(s11B[i] != 'C') ++ng11B;
        if(s11C[i] != 'C') ++ng11C;
        if(s11E[i] != 'C') ++ng11E;
        if(s11F[i] != 'C') ++ng11F;
        if(s11W[i] != 'C') ++ng11W;
        if(s12A[i] != 'C') ++ng12A;
        if(s12B[i] != 'C') ++ng12B;
        if(s12D[i] != 'C') ++ng12D;
        if(s12E[i] != 'C') ++ng12E;
        if(s12K[i] != 'C') ++ng12K;
        if(s13A[i] != 'C') ++ng13A;
        if(s13B[i] != 'C') ++ng13B;
        if(s13K[i] != 'C') ++ng13K;
        if(sFCC[i] != 'C') ++ngFCC;
        if(sHCP[i] != 'C') ++ngHCP;
        if(sBCC_9[i] != 'C') ++ngBCC_9;
        if(sBCC_15[i] != 'C') ++ngBCC_15;
    }
    Stats_SetA();
    for(i=0; i<N; ++i){
        if(s5A[i] != 'C' && !a6[i]) ++nn5A;
        if(s6A[i] != 'C' && !a7[i]) ++nn6A;
        if(s6Z[i] != 'C' && !a7[i]) ++nn6Z;
        if(s7K[i] != 'C' && !a7[i]) ++nn7K;
        if(s7A[i] != 'C' && !a8[i]) ++nn7A;
        if(s8A[i] != 'C' && !a9[i]) ++nn8A;
        if(s8B[i] != 'C' && !a9[i]) ++nn8B;
        if(s8K[i] != 'C' && !a9[i]) ++nn8K;
        if(s9A[i] != 'C' && !a10[i]) ++nn9A;
        if(s9B[i] != 'C' && !a10[i]) ++nn9B;
        if(s9K[i] != 'C' && !a10[i]) ++nn9K;
        if(s10A[i] != 'C' && !a11[i]) ++nn10A;
        if(s10B[i] != 'C' && !a11[i]) ++nn10B;
        if(s10K[i] != 'C' && !a11[i]) ++nn10K;
        if(s10W[i] != 'C' && !a11[i]) ++nn10W;
        if(s11A[i] != 'C' && !a12[i]) ++nn11A;
        if(s11B[i] != 'C' && !a12[i]) ++nn11B;
        if(s11C[i] != 'C' && !a12[i]) ++nn11C;
        if(s11E[i] != 'C' && !a12[i]) ++nn11E;
        if(s11F[i] != 'C' && !a12[i]) ++nn11F;
        if(s11W[i] != 'C' && !a12[i]) ++nn11W;
        if(s12A[i] != 'C' && !a13[i]) ++nn12A;
        if(s12B[i] != 'C' && !a13[i]) ++nn12B;
        if(s12D[i] != 'C' && !a13[i]) ++nn12D;
        if(s12E[i] != 'C' && !a13[i]) ++nn12E;
        if(s12K[i] != 'C' && !a13[i]) ++nn12K;
        if(s13A[i] != 'C') ++nn13A;
        if(s13B[i] != 'C') ++nn13B;
        if(s13K[i] != 'C') ++nn13K;
        if(sFCC[i] != 'C') ++nnFCC;
        if(sHCP[i] != 'C') ++nnHCP;
        if(sBCC_9[i] != 'C') ++nnBCC_9;
        if(sBCC_15[i] != 'C') ++nnBCC_15;
    }
    Nclus = 0;
    for(i=0; i<N; ++i){
        flg = 0;
        if(s5A[i] != 'C') flg=1;
        if(s6A[i] != 'C') flg=1;
        if(s6Z[i] != 'C') flg=1;
        if(s7K[i] != 'C') flg=1;
        if(s7A[i] != 'C') flg=1;
        if(s8A[i] != 'C') flg=1;
        if(s8B[i] != 'C') flg=1;
        if(s8K[i] != 'C') flg=1;
        if(s9A[i] != 'C') flg=1;
        if(s9B[i] != 'C') flg=1;
        if(s9K[i] != 'C') flg=1;
        if(s10A[i] != 'C') flg=1;
        if(s10B[i] != 'C') flg=1;
        if(s10K[i] != 'C') flg=1;
        if(s10W[i] != 'C') flg=1;
        if(s11A[i] != 'C') flg=1;
        if(s11B[i] != 'C') flg=1;
        if(s11C[i] != 'C') flg=1;
        if(s11E[i] != 'C') flg=1;
        if(s11F[i] != 'C') flg=1;
        if(s11W[i] != 'C') flg=1;
        if(s12A[i] != 'C') flg=1;
        if(s12B[i] != 'C') flg=1;
        if(s12D[i] != 'C') flg=1;
        if(s12E[i] != 'C') flg=1;
        if(s12K[i] != 'C') flg=1;
        if(s13A[i] != 'C') flg=1;
        if(s13B[i] != 'C') flg=1;
        if(s13K[i] != 'C') flg=1;
        if(sFCC[i] != 'C') flg=1;
        if(sHCP[i] != 'C') flg=1;
        if(sBCC_9[i] != 'C') flg=1;
        if(sBCC_15[i] != 'C') flg=1;
        if(flg) ++Nclus;
    }
    totNclus+=Nclus;
    if (PRINTINFO==1) printf("Nclus = %d\n",Nclus);
}

void Stats_Report(char *filename) {
    int f;
    int clusTot=0;
    double temp;
    char errMsg[1000];
    FILE *writeout;

    for (f=0; f<FRAMES; f++) {
        ncsp3+=nsp3[f];
        ncsp3a+=nsp3a[f];
        ncsp3b+=nsp3b[f];
        nc5A+=nsp3c[f];
        ncsp3_excess_spindles+=nsp3_excess_spindles[f];
        ncsp3c_spindlebonds+=nsp3c_spindlebonds[f];
        ncsp4+=nsp4[f];
        ncsp4a+=nsp4a[f];
        ncsp4b+=nsp4b[f];
        ncsp4c+=nsp4c[f];
        ncsp4_excess_spindles+=nsp4_excess_spindles[f];
        ncsp4c_spindlebonds+=nsp4c_spindlebonds[f];
        nc6A+=n6A[f];
        nc6A_spindlebonds+=n6A_spindlebonds[f];
        nc6Z+=n6Z[f];
        nc7K+=n7K[f];
        ncsp5+=nsp5[f];
        ncsp5a+=nsp5a[f];
        ncsp5b+=nsp5b[f];
        nc7A+=nsp5c[f];
        ncsp5_excess_spindles+=nsp5_excess_spindles[f];
        ncsp5c_spindlebonds+=nsp5c_spindlebonds[f];
        nc8A+=n8A[f];
        nc8B+=n8B[f];
        nc8K+=n8K[f];
        nc9A+=n9A[f];
        nc9B+=n9B[f];
        nc9K+=n9K[f];
        nc10A+=n10A[f];
        nc10B+=n10B[f];
        nc10K+=n10K[f];
        nc10W+=n10W[f];
        nc11A+=n11A[f];
        nc11B+=n11B[f];
        nc11C+=n11C[f];
        nc11E+=n11E[f];
        nc11F+=n11F[f];
        nc11W+=n11W[f];
        nc12A+=n12A[f];
        nc12B+=n12B[f];
        nc12D+=n12D[f];
        nc12E+=n12E[f];
        nc12K+=n12K[f];
        nc13A+=n13A[f];
        nc13B+=n13B[f];
        nc13K+=n13K[f];
        ncFCC+=nFCC[f];
        ncHCP+=nHCP[f];
        ncBCC_9+=nBCC_9[f];
        ncBCC_15+=nBCC_15[f];
        clusTot+=nsp3a[f]+nsp3b[f]+nsp3c[f]+nsp4a[f]+nsp4b[f]+n6A[f]+nsp5a[f]+nsp5b[f]+nsp5c[f];
        clusTot+=n6Z[f]+n7K[f]+n8A[f]+n8B[f]+n8K[f]+n9A[f]+n9B[f]+n9K[f]+n10A[f]+n10B[f]+n10K[f]+n10W[f];
        clusTot+=n11A[f]+n11B[f]+n11C[f]+n11E[f]+n11F[f]+n11W[f];
        clusTot+=n12A[f]+n12B[f]+n12D[f]+n12E[f]+n12K[f];
        clusTot+=n13A[f]+n13B[f]+n13K[f]+nFCC[f]+nHCP[f]+nBCC_9[f]+nBCC_15[f];
    }

    printf("Clust	Number	Gross	Net	Mean Pop Per Frame\n");
    printf("sp3	%d	%d	%d	%.5lg\n",ncsp3,ngsp3,0,mean_pop_per_frame_sp3);
    printf("sp3a	%d	%d	%d	%.5lg\n",ncsp3a,ngsp3a,0,mean_pop_per_frame_sp3a);
    printf("sp3b	%d	%d	%d	%.5lg\n",ncsp3b,ngsp3b,0,mean_pop_per_frame_sp3b);
    printf("5A_D3h	%d	%d	%d	%.5lg\n",nc5A,ng5A,nn5A,mean_pop_per_frame_sp3c);
    printf("sp4	%d	%d	%d	%.5lg\n",ncsp4,ngsp4,0,mean_pop_per_frame_sp4);
    printf("sp4a	%d	%d	%d	%.5lg\n",ncsp4a,ngsp4a,0,mean_pop_per_frame_sp4a);
    printf("sp4b	%d	%d	%d	%.5lg\n",ncsp4b,ngsp4b,0,mean_pop_per_frame_sp4b);
    printf("sp4c	%d	%d	%d	%.5lg\n",ncsp4c,ng6A,0,mean_pop_per_frame_6A);
    printf("6A_Oh	%d	%d	%d	%.5lg\n",nc6A,ng6A,nn6A,mean_pop_per_frame_6A);
    printf("6Z_C2v	%d	%d	%d	%.5lg\n",nc6Z,ng6Z,nn6Z,mean_pop_per_frame_6Z);
    printf("sp5	%d	%d	%d	%.5lg\n",ncsp5,ngsp5,0,mean_pop_per_frame_sp5);
    printf("sp5a	%d	%d	%d	%.5lg\n",ncsp5a,ngsp5a,0,mean_pop_per_frame_sp5a);
    printf("sp5b	%d	%d	%d	%.5lg\n",ncsp5b,ngsp5b,0,mean_pop_per_frame_sp5b);
    printf("7A_D5h	%d	%d	%d	%.5lg\n",nc7A,ng7A,nn7A,mean_pop_per_frame_sp5c);
    printf("7K	%d	%d	%d	%.5lg\n",nc7K,ng7K,nn7K,mean_pop_per_frame_7K);
    printf("8A_D2d	%d	%d	%d	%.5lg\n",nc8A,ng8A,nn8A,mean_pop_per_frame_8A);
    printf("8B_Cs	%d	%d	%d	%.5lg\n",nc8B,ng8B,nn8B,mean_pop_per_frame_8B);
    printf("8K	%d	%d	%d	%.5lg\n",nc8K,ng8K,nn8K,mean_pop_per_frame_8K);
    printf("9A_D3h	%d	%d	%d	%.5lg\n",nc9A,ng9A,nn9A,mean_pop_per_frame_9A);
    printf("9B_C2v	%d	%d	%d	%.5lg\n",nc9B,ng9B,nn9B,mean_pop_per_frame_9B);
    printf("9K	%d	%d	%d	%.5lg\n",nc9K,ng9K,nn9K,mean_pop_per_frame_9K);
    printf("10A_D4d	%d	%d	%d	%.5lg\n",nc10A,ng10A,nn10A,mean_pop_per_frame_10A);
    printf("10B_C3v	%d	%d	%d	%.5lg\n",nc10B,ng10B,nn10B,mean_pop_per_frame_10B);
    printf("10K	%d	%d	%d	%.5lg\n",nc10K,ng10K,nn10K,mean_pop_per_frame_10K);
    printf("10W	%d	%d	%d	%.5lg\n",nc10W,ng10W,nn10W,mean_pop_per_frame_10W);
    printf("11A_D4d	%d	%d	%d	%.5lg\n",nc11A,ng11A,nn11A,mean_pop_per_frame_11A);
    printf("11B_C2v	%d	%d	%d	%.5lg\n",nc11B,ng11B,nn11B,mean_pop_per_frame_11B);
    printf("11CD	%d	%d	%d	%.5lg\n",nc11C,ng11C,nn11C,mean_pop_per_frame_11C);
    printf("11E_C2	%d	%d	%d	%.5lg\n",nc11E,ng11E,nn11E,mean_pop_per_frame_11E);
    printf("11F_C2v	%d	%d	%d	%.5lg\n",nc11F,ng11F,nn11F,mean_pop_per_frame_11F);
    printf("11W_Cs	%d	%d	%d	%.5lg\n",nc11W,ng11W,nn11W,mean_pop_per_frame_11W);
    printf("12A_C2v	%d	%d	%d	%.5lg\n",nc12A,ng12A,nn12A,mean_pop_per_frame_12A);
    printf("12BC	%d	%d	%d	%.5lg\n",nc12B,ng12B,nn12B,mean_pop_per_frame_12B);
    printf("12D_D2d	%d	%d	%d	%.5lg\n",nc12D,ng12D,nn12D,mean_pop_per_frame_12D);
    printf("12E_D3h	%d	%d	%d	%.5lg\n",nc12E,ng12E,nn12E,mean_pop_per_frame_12E);
    printf("12K	%d	%d	%d	%.5lg\n",nc12K,ng12K,nn12K,mean_pop_per_frame_12K);
    printf("13A_Ih	%d	%d	%d	%.5lg\n",nc13A,ng13A,nn13A,mean_pop_per_frame_13A);
    printf("13B_D5h	%d	%d	%d	%.5lg\n",nc13B,ng13B,nn13B,mean_pop_per_frame_13B);
    printf("13K	%d	%d	%d	%.5lg\n",nc13K,ng13K,nn13K,mean_pop_per_frame_13K);
    printf("FCC_m13	%d	%d	%d	%.5lg\n",ncFCC,ngFCC,nnFCC,mean_pop_per_frame_FCC);
    printf("HCP_m13	%d	%d	%d	%.5lg\n",ncHCP,ngHCP,nnHCP,mean_pop_per_frame_HCP);
    printf("BCC_m9	%d	%d	%d	%.5lg\n",ncBCC_9,ngBCC_9,nnBCC_9,mean_pop_per_frame_BCC_9);
    printf("BCC_m15	%d	%d	%d	%.5lg\n",ncBCC_15,ngBCC_15,nnBCC_15,mean_pop_per_frame_BCC_15);
    printf("totals	%d	%d	%d\n",clusTot,totNclus,N*FRAMES);
    printf("maxnB	%d\n",maxnb);
    printf("max to sp3	%d\n",maxto3);
    printf("max to sp4	%d\n",maxto4);
    printf("max to sp5	%d\n",maxto5);
    printf("excess spindles sp3	%d\n",ncsp3_excess_spindles);
    printf("excess spindles sp4	%d\n",ncsp4_excess_spindles);
    printf("excess spindles sp5	%d\n",ncsp5_excess_spindles);

    temp=(double)ncsp3c_spindlebonds/(double)nc5A;
    printf("5A spindle bonds	%d	%.15lg\n",ncsp3c_spindlebonds,temp);
    temp=(double)ncsp4c_spindlebonds/(double)ncsp4c;
    printf("sp4c spindle bonds	%d	%.15lg\n",ncsp4c_spindlebonds,temp);
    temp=(double)nc6A_spindlebonds/(double)nc6A;
    printf("6A spindle bonds	%d	%.15lg\n",nc6A_spindlebonds,temp);
    temp=(double)ncsp5c_spindlebonds/(double)nc7A;
    printf("7A spindle bonds	%d	%.15lg\n",ncsp5c_spindlebonds,temp);
    printf("correctedBonds %d per frame %.5lg per part per frame %.5lg\n",correctedBonds,(double)correctedBonds/FRAMES,(double)correctedBonds/(FRAMES*N));


    writeout=fopen(filename,"w");
    if (writeout==NULL)  {
        sprintf(errMsg,"Stats_Report(): Error opening file %s",filename);	// Always test file open
        Error(errMsg);
    }
    fprintf(writeout,"%s\n",filename);
    fprintf(writeout,"Clust	Number	Gross	Net	Mean Pop Per Frame\n");
    fprintf(writeout,"sp3	%d	%d	%d	%.5lg\n",ncsp3,ngsp3,0,mean_pop_per_frame_sp3);
    fprintf(writeout,"sp3a	%d	%d	%d	%.5lg\n",ncsp3a,ngsp3a,0,mean_pop_per_frame_sp3a);
    fprintf(writeout,"sp3b	%d	%d	%d	%.5lg\n",ncsp3b,ngsp3b,0,mean_pop_per_frame_sp3b);
    fprintf(writeout,"5A_D3h	%d	%d	%d	%.5lg\n",nc5A,ng5A,nn5A,mean_pop_per_frame_sp3c);
    fprintf(writeout,"sp4	%d	%d	%d	%.5lg\n",ncsp4,ngsp4,0,mean_pop_per_frame_sp4);
    fprintf(writeout,"sp4a	%d	%d	%d	%.5lg\n",ncsp4a,ngsp4a,0,mean_pop_per_frame_sp4a);
    fprintf(writeout,"sp4b	%d	%d	%d	%.5lg\n",ncsp4b,ngsp4b,0,mean_pop_per_frame_sp4b);
    fprintf(writeout,"sp4c	%d	%d	%d	%.5lg\n",ncsp4c,ng6A,0,mean_pop_per_frame_6A);
    fprintf(writeout,"6A_Oh	%d	%d	%d	%.5lg\n",nc6A,ng6A,nn6A,mean_pop_per_frame_6A);
    fprintf(writeout,"6Z_C2v	%d	%d	%d	%.5lg\n",nc6Z,ng6Z,nn6Z,mean_pop_per_frame_6Z);
    fprintf(writeout,"sp5	%d	%d	%d	%.5lg\n",ncsp5,ngsp5,0,mean_pop_per_frame_sp5);
    fprintf(writeout,"sp5a	%d	%d	%d	%.5lg\n",ncsp5a,ngsp5a,0,mean_pop_per_frame_sp5a);
    fprintf(writeout,"sp5b	%d	%d	%d	%.5lg\n",ncsp5b,ngsp5b,0,mean_pop_per_frame_sp5b);
    fprintf(writeout,"7A_D5h	%d	%d	%d	%.5lg\n",nc7A,ng7A,nn7A,mean_pop_per_frame_sp5c);
    fprintf(writeout,"7K	%d	%d	%d	%.5lg\n",nc7K,ng7K,nn7K,mean_pop_per_frame_7K);
    fprintf(writeout,"8A_D2d	%d	%d	%d	%.5lg\n",nc8A,ng8A,nn8A,mean_pop_per_frame_8A);
    fprintf(writeout,"8B_Cs	%d	%d	%d	%.5lg\n",nc8B,ng8B,nn8B,mean_pop_per_frame_8B);
    fprintf(writeout,"8K	%d	%d	%d	%.5lg\n",nc8K,ng8K,nn8K,mean_pop_per_frame_8K);
    fprintf(writeout,"9A_D3h	%d	%d	%d	%.5lg\n",nc9A,ng9A,nn9A,mean_pop_per_frame_9A);
    fprintf(writeout,"9B_C2v	%d	%d	%d	%.5lg\n",nc9B,ng9B,nn9B,mean_pop_per_frame_9B);
    fprintf(writeout,"9K	%d	%d	%d	%.5lg\n",nc9K,ng9K,nn9K,mean_pop_per_frame_9K);
    fprintf(writeout,"10A_D4d	%d	%d	%d	%.5lg\n",nc10A,ng10A,nn10A,mean_pop_per_frame_10A);
    fprintf(writeout,"10B_C3v	%d	%d	%d	%.5lg\n",nc10B,ng10B,nn10B,mean_pop_per_frame_10B);
    fprintf(writeout,"10K	%d	%d	%d	%.5lg\n",nc10K,ng10K,nn10K,mean_pop_per_frame_10K);
    fprintf(writeout,"10W	%d	%d	%d	%.5lg\n",nc10W,ng10W,nn10W,mean_pop_per_frame_10W);
    fprintf(writeout,"11A_D4d	%d	%d	%d	%.5lg\n",nc11A,ng11A,nn11A,mean_pop_per_frame_11A);
    fprintf(writeout,"11B_C2v	%d	%d	%d	%.5lg\n",nc11B,ng11B,nn11B,mean_pop_per_frame_11B);
    fprintf(writeout,"11CD	%d	%d	%d	%.5lg\n",nc11C,ng11C,nn11C,mean_pop_per_frame_11C);
    fprintf(writeout,"11E_C2	%d	%d	%d	%.5lg\n",nc11E,ng11E,nn11E,mean_pop_per_frame_11E);
    fprintf(writeout,"11F_C2v	%d	%d	%d	%.5lg\n",nc11F,ng11F,nn11F,mean_pop_per_frame_11F);
    fprintf(writeout,"11W_Cs	%d	%d	%d	%.5lg\n",nc11W,ng11W,nn11W,mean_pop_per_frame_11W);
    fprintf(writeout,"12A_C2v	%d	%d	%d	%.5lg\n",nc12A,ng12A,nn12A,mean_pop_per_frame_12A);
    fprintf(writeout,"12BC	%d	%d	%d	%.5lg\n",nc12B,ng12B,nn12B,mean_pop_per_frame_12B);
    fprintf(writeout,"12D_D2d	%d	%d	%d	%.5lg\n",nc12D,ng12D,nn12D,mean_pop_per_frame_12D);
    fprintf(writeout,"12E_D3h	%d	%d	%d	%.5lg\n",nc12E,ng12E,nn12E,mean_pop_per_frame_12E);
    fprintf(writeout,"12K	%d	%d	%d	%.5lg\n",nc12K,ng12K,nn12K,mean_pop_per_frame_12K);
    fprintf(writeout,"13A_Ih	%d	%d	%d	%.5lg\n",nc13A,ng13A,nn13A,mean_pop_per_frame_13A);
    fprintf(writeout,"13B_D5h	%d	%d	%d	%.5lg\n",nc13B,ng13B,nn13B,mean_pop_per_frame_13B);
    fprintf(writeout,"13K	%d	%d	%d	%.5lg\n",nc13K,ng13K,nn13K,mean_pop_per_frame_13K);
    fprintf(writeout,"FCC_m13	%d	%d	%d	%.5lg\n",ncFCC,ngFCC,nnFCC,mean_pop_per_frame_FCC);
    fprintf(writeout,"HCP_m13	%d	%d	%d	%.5lg\n",ncHCP,ngHCP,nnHCP,mean_pop_per_frame_HCP);
    fprintf(writeout,"BCC_m9	%d	%d	%d	%.5lg\n",ncBCC_9,ngBCC_9,nnBCC_9,mean_pop_per_frame_BCC_9);
    fprintf(writeout,"BCC_m15	%d	%d	%d	%.5lg\n",ncBCC_15,ngBCC_15,nnBCC_15,mean_pop_per_frame_BCC_15);
    fprintf(writeout,"totals	%d	%d	%d\n",clusTot,totNclus,N*FRAMES);
    fprintf(writeout,"maxnB	%d\n",maxnb);
    fprintf(writeout,"max to sp3	%d\n",maxto3);
    fprintf(writeout,"max to sp4	%d\n",maxto4);
    fprintf(writeout,"max to sp5	%d\n",maxto5);
    fprintf(writeout,"excess spindles sp3	%d\n",ncsp3_excess_spindles);
    fprintf(writeout,"excess spindles sp4	%d\n",ncsp4_excess_spindles);
    fprintf(writeout,"excess spindles sp5	%d\n",ncsp5_excess_spindles);

    temp=(double)ncsp3c_spindlebonds/(double)nc5A;
    fprintf(writeout,"5A spindle bonds	%d	%.15lg\n",ncsp3c_spindlebonds,temp);
    temp=(double)ncsp4c_spindlebonds/(double)ncsp4c;
    fprintf(writeout,"sp4c spindle bonds	%d	%.15lg\n",ncsp4c_spindlebonds,temp);
    temp=(double)nc6A_spindlebonds/(double)nc6A;
    fprintf(writeout,"6A spindle bonds	%d	%.15lg\n",nc6A_spindlebonds,temp);
    temp=(double)ncsp5c_spindlebonds/(double)nc7A;
    fprintf(writeout,"7A spindle bonds	%d	%.15lg\n",ncsp5c_spindlebonds,temp);
    fprintf(writeout,"correctedBonds	%d	per frame	%.15lg	per part per frame	%.15lg\n",correctedBonds,(double)correctedBonds/FRAMES,(double)correctedBonds/(FRAMES*N));

    fclose(writeout);
}

void Pop_Per_Frame(int f) {
    int i;

    for(i=0; i<N; ++i){
        if(ssp3[i] != 'C') pop_per_frame_sp3[f]+=1.0;
        if(ssp3a[i] != 'C') pop_per_frame_sp3a[f]+=1.0;
        if(ssp3b[i] != 'C') pop_per_frame_sp3b[f]+=1.0;
        if(s5A[i] != 'C') pop_per_frame_sp3c[f]+=1.0;
        if(ssp4[i] != 'C') pop_per_frame_sp4[f]+=1.0;
        if(ssp4a[i] != 'C') pop_per_frame_sp4a[f]+=1.0;
        if(ssp4b[i] != 'C') pop_per_frame_sp4b[f]+=1.0;
        if(s6A[i] != 'C') pop_per_frame_6A[f]+=1.0;
        if(s6Z[i] != 'C') pop_per_frame_6Z[f]+=1.0;
        if(ssp5[i] != 'C') pop_per_frame_sp5[f]+=1.0;
        if(ssp5a[i] != 'C') pop_per_frame_sp5a[f]+=1.0;
        if(ssp5b[i] != 'C') pop_per_frame_sp5b[f]+=1.0;
        if(s7A[i] != 'C') pop_per_frame_sp5c[f]+=1.0;
        if(s7K[i] != 'C') pop_per_frame_7K[f]+=1.0;
        if(s8A[i] != 'C') pop_per_frame_8A[f]+=1.0;
        if(s8B[i] != 'C') pop_per_frame_8B[f]+=1.0;
        if(s8K[i] != 'C') pop_per_frame_8K[f]+=1.0;
        if(s9A[i] != 'C') pop_per_frame_9A[f]+=1.0;
        if(s9B[i] != 'C') pop_per_frame_9B[f]+=1.0;
        if(s9K[i] != 'C') pop_per_frame_9K[f]+=1.0;
        if(s10A[i] != 'C') pop_per_frame_10A[f]+=1.0;
        if(s10B[i] != 'C') pop_per_frame_10B[f]+=1.0;
        if(s10K[i] != 'C') pop_per_frame_10K[f]+=1.0;
        if(s10W[i] != 'C') pop_per_frame_10W[f]+=1.0;
        if(s11A[i] != 'C') pop_per_frame_11A[f]+=1.0;
        if(s11B[i] != 'C') pop_per_frame_11B[f]+=1.0;
        if(s11C[i] != 'C') pop_per_frame_11C[f]+=1.0;
        if(s11E[i] != 'C') pop_per_frame_11E[f]+=1.0;
        if(s11F[i] != 'C') pop_per_frame_11F[f]+=1.0;
        if(s11W[i] != 'C') pop_per_frame_11W[f]+=1.0;
        if(s12A[i] != 'C') pop_per_frame_12A[f]+=1.0;
        if(s12B[i] != 'C') pop_per_frame_12B[f]+=1.0;
        if(s12D[i] != 'C') pop_per_frame_12D[f]+=1.0;
        if(s12E[i] != 'C') pop_per_frame_12E[f]+=1.0;
        if(s12K[i] != 'C') pop_per_frame_12K[f]+=1.0;
        if(s13A[i] != 'C') pop_per_frame_13A[f]+=1.0;
        if(s13B[i] != 'C') pop_per_frame_13B[f]+=1.0;
        if(s13K[i] != 'C') pop_per_frame_13K[f]+=1.0;
        if(sFCC[i] != 'C') pop_per_frame_FCC[f]+=1.0;
        if(sHCP[i] != 'C') pop_per_frame_HCP[f]+=1.0;
        if(sBCC_9[i] != 'C') pop_per_frame_BCC_9[f]+=1.0;
        if(sBCC_15[i] != 'C') pop_per_frame_BCC_15[f]+=1.0;
    }

    pop_per_frame_sp3[f]/=N;
    pop_per_frame_sp3a[f]/=N;
    pop_per_frame_sp3b[f]/=N;
    pop_per_frame_sp3c[f]/=N;
    pop_per_frame_sp4[f]/=N;
    pop_per_frame_sp4a[f]/=N;
    pop_per_frame_sp4b[f]/=N;
    pop_per_frame_6A[f]/=N;
    pop_per_frame_6Z[f]/=N;
    pop_per_frame_sp5[f]/=N;
    pop_per_frame_sp5a[f]/=N;
    pop_per_frame_sp5b[f]/=N;
    pop_per_frame_sp5c[f]/=N;
    pop_per_frame_7K[f]/=N;
    pop_per_frame_8A[f]/=N;
    pop_per_frame_8B[f]/=N;
    pop_per_frame_8K[f]/=N;
    pop_per_frame_9A[f]/=N;
    pop_per_frame_9B[f]/=N;
    pop_per_frame_9K[f]/=N;
    pop_per_frame_10A[f]/=N;
    pop_per_frame_10B[f]/=N;
    pop_per_frame_10K[f]/=N;
    pop_per_frame_10W[f]/=N;
    pop_per_frame_11A[f]/=N;
    pop_per_frame_11B[f]/=N;
    pop_per_frame_11C[f]/=N;
    pop_per_frame_11E[f]/=N;
    pop_per_frame_11F[f]/=N;
    pop_per_frame_11W[f]/=N;
    pop_per_frame_12A[f]/=N;
    pop_per_frame_12B[f]/=N;
    pop_per_frame_12D[f]/=N;
    pop_per_frame_12E[f]/=N;
    pop_per_frame_12K[f]/=N;
    pop_per_frame_13A[f]/=N;
    pop_per_frame_13B[f]/=N;
    pop_per_frame_13K[f]/=N;
    pop_per_frame_FCC[f]/=N;
    pop_per_frame_HCP[f]/=N;
    pop_per_frame_BCC_9[f]/=N;
    pop_per_frame_BCC_15[f]/=N;

    mean_pop_per_frame_sp3+=pop_per_frame_sp3[f];
    mean_pop_per_frame_sp3a+=pop_per_frame_sp3a[f];
    mean_pop_per_frame_sp3b+=pop_per_frame_sp3b[f];
    mean_pop_per_frame_sp3c+=pop_per_frame_sp3c[f];
    mean_pop_per_frame_sp4+=pop_per_frame_sp4[f];
    mean_pop_per_frame_sp4a+=pop_per_frame_sp4a[f];
    mean_pop_per_frame_sp4b+=pop_per_frame_sp4b[f];
    mean_pop_per_frame_6A+=pop_per_frame_6A[f];
    mean_pop_per_frame_6Z+=pop_per_frame_6Z[f];
    mean_pop_per_frame_sp5+=pop_per_frame_sp5[f];
    mean_pop_per_frame_sp5a+=pop_per_frame_sp5a[f];
    mean_pop_per_frame_sp5b+=pop_per_frame_sp5b[f];
    mean_pop_per_frame_sp5c+=pop_per_frame_sp5c[f];
    mean_pop_per_frame_7K+=pop_per_frame_7K[f];
    mean_pop_per_frame_8A+=pop_per_frame_8A[f];
    mean_pop_per_frame_8B+=pop_per_frame_8B[f];
    mean_pop_per_frame_8K+=pop_per_frame_8K[f];
    mean_pop_per_frame_9A+=pop_per_frame_9A[f];
    mean_pop_per_frame_9B+=pop_per_frame_9B[f];
    mean_pop_per_frame_9K+=pop_per_frame_9K[f];
    mean_pop_per_frame_10A+=pop_per_frame_10A[f];
    mean_pop_per_frame_10B+=pop_per_frame_10B[f];
    mean_pop_per_frame_10K+=pop_per_frame_10K[f];
    mean_pop_per_frame_10W+=pop_per_frame_10W[f];
    mean_pop_per_frame_11A+=pop_per_frame_11A[f];
    mean_pop_per_frame_11B+=pop_per_frame_11B[f];
    mean_pop_per_frame_11C+=pop_per_frame_11C[f];
    mean_pop_per_frame_11E+=pop_per_frame_11E[f];
    mean_pop_per_frame_11F+=pop_per_frame_11F[f];
    mean_pop_per_frame_11W+=pop_per_frame_11W[f];
    mean_pop_per_frame_12A+=pop_per_frame_12A[f];
    mean_pop_per_frame_12B+=pop_per_frame_12B[f];
    mean_pop_per_frame_12D+=pop_per_frame_12D[f];
    mean_pop_per_frame_12E+=pop_per_frame_12E[f];
    mean_pop_per_frame_12K+=pop_per_frame_12K[f];
    mean_pop_per_frame_13A+=pop_per_frame_13A[f];
    mean_pop_per_frame_13B+=pop_per_frame_13B[f];
    mean_pop_per_frame_13K+=pop_per_frame_13K[f];
    mean_pop_per_frame_FCC+=pop_per_frame_FCC[f];
    mean_pop_per_frame_HCP+=pop_per_frame_HCP[f];
    mean_pop_per_frame_BCC_9+=pop_per_frame_BCC_9[f];
    mean_pop_per_frame_BCC_15+=pop_per_frame_BCC_15[f];
}

void Normalise_Populations() {
    mean_pop_per_frame_sp3/=FRAMES;
    mean_pop_per_frame_sp3a/=FRAMES;
    mean_pop_per_frame_sp3b/=FRAMES;
    mean_pop_per_frame_sp3c/=FRAMES;
    mean_pop_per_frame_sp4/=FRAMES;
    mean_pop_per_frame_sp4a/=FRAMES;
    mean_pop_per_frame_sp4b/=FRAMES;
    mean_pop_per_frame_6A/=FRAMES;
    mean_pop_per_frame_6Z/=FRAMES;
    mean_pop_per_frame_sp5/=FRAMES;
    mean_pop_per_frame_sp5a/=FRAMES;
    mean_pop_per_frame_sp5b/=FRAMES;
    mean_pop_per_frame_sp5c/=FRAMES;
    mean_pop_per_frame_7K/=FRAMES;
    mean_pop_per_frame_8A/=FRAMES;
    mean_pop_per_frame_8B/=FRAMES;
    mean_pop_per_frame_8K/=FRAMES;
    mean_pop_per_frame_9A/=FRAMES;
    mean_pop_per_frame_9B/=FRAMES;
    mean_pop_per_frame_9K/=FRAMES;
    mean_pop_per_frame_10A/=FRAMES;
    mean_pop_per_frame_10B/=FRAMES;
    mean_pop_per_frame_10K/=FRAMES;
    mean_pop_per_frame_10W/=FRAMES;
    mean_pop_per_frame_11A/=FRAMES;
    mean_pop_per_frame_11B/=FRAMES;
    mean_pop_per_frame_11C/=FRAMES;
    mean_pop_per_frame_11E/=FRAMES;
    mean_pop_per_frame_11F/=FRAMES;
    mean_pop_per_frame_11W/=FRAMES;
    mean_pop_per_frame_12A/=FRAMES;
    mean_pop_per_frame_12B/=FRAMES;
    mean_pop_per_frame_12D/=FRAMES;
    mean_pop_per_frame_12E/=FRAMES;
    mean_pop_per_frame_12K/=FRAMES;
    mean_pop_per_frame_13A/=FRAMES;
    mean_pop_per_frame_13B/=FRAMES;
    mean_pop_per_frame_13K/=FRAMES;
    mean_pop_per_frame_FCC/=FRAMES;
    mean_pop_per_frame_HCP/=FRAMES;
    mean_pop_per_frame_BCC_9/=FRAMES;
    mean_pop_per_frame_BCC_15/=FRAMES;
}