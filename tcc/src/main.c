/* Alex Malins - alex.malins@gmail.com */
/* TCC: A topological cluster classification code with temporal tracking of clusters. */
/* Not for general consumption */	


#include <string.h>
#include "globals.h"
#include "setup.h"
#include "tools.h"
#include "bonds.h"
#include "rings.h"
#include "clusters.h"

int icell(int tix, int tiy, int tiz) { 	// returns cell number (from 1 to ncells) for given (tix,tiy,tiz) coordinate
	return 1 + (tix-1+M)%M + M*((tiy-1+M)%M) + M*M*((tiz-1+M)%M); 
}

void Rings_gSP3(int f, int n0) {	// get SP3/4/5 rings including particle n0
	int i,j;
	int n1, n2;
	
	for (i=0; i<cnb[n0]-1; i++){
		n1=bNums[n0][i];
		if (n1 < n0) continue;
		for (j=i+1; j<cnb[n0]; ++j){
			n2=bNums[n0][j];
			if (n2<n0) continue;
			if (Bonds_BondCheck(n1,n2)) { // is n1 bonded to n2
				if (n1<n2) Rings_aSP3(f,n0,n1,n2); // SP3 found, check type and store
				else Rings_aSP3(f,n0,n2,n1); // SP3 found, check type and store
			} 
			else { // not SP3, search for SP4 & SP5
				if (dosp4==1) {
					if (n1<n2) Rings_gSP4(f,n0,n1,n2);
					else Rings_gSP4(f,n0,n2,n1);
				}
			}
		} 		
	}
}


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


void Write_Raw_Init() {
	char errMsg[1000];
	char output[1000];
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp3",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_sp3=fopen(output,"w");
	if (file_raw_sp3==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_sp3,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp3a",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_sp3a=fopen(output,"w");
	if (file_raw_sp3a==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_sp3a,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp3b",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_sp3b=fopen(output,"w");
	if (file_raw_sp3b==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_sp3b,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_5A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_5A=fopen(output,"w");
	if (file_raw_5A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_5A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp4",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_sp4=fopen(output,"w");
	if (file_raw_sp4==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_sp4,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp4a",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_sp4a=fopen(output,"w");
	if (file_raw_sp4a==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_sp4a,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp4b",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_sp4b=fopen(output,"w");
	if (file_raw_sp4b==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_sp4b,"%s\n",output);

	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_6A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_6A=fopen(output,"w");
	if (file_raw_6A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_6A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_6Z",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_6Z=fopen(output,"w");
	if (file_raw_6Z==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_6Z,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp5",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_sp5=fopen(output,"w");
	if (file_raw_sp5==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_sp5,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp5a",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_sp5a=fopen(output,"w");
	if (file_raw_sp5a==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_sp5a,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp5b",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_sp5b=fopen(output,"w");
	if (file_raw_sp5b==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_sp5b,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_7A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_7A=fopen(output,"w");
	if (file_raw_7A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_7A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_7K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_7K=fopen(output,"w");
	if (file_raw_7K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_7K,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_8A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_8A=fopen(output,"w");
	if (file_raw_8A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_8A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_8B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_8B=fopen(output,"w");
	if (file_raw_8B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_8B,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_8K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_8K=fopen(output,"w");
	if (file_raw_8K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_8K,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_9A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_9A=fopen(output,"w");
	if (file_raw_9A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_9A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_9B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_9B=fopen(output,"w");
	if (file_raw_9B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_9B,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_9K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_9K=fopen(output,"w");
	if (file_raw_9K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_9K,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_10A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_10A=fopen(output,"w");
	if (file_raw_10A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_10A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_10B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_10B=fopen(output,"w");
	if (file_raw_10B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_10B,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_10K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_10K=fopen(output,"w");
	if (file_raw_10K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_10K,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_10W",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_10W=fopen(output,"w");
	if (file_raw_10W==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_10W,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_11A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_11A=fopen(output,"w");
	if (file_raw_11A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_11A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_11B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_11B=fopen(output,"w");
	if (file_raw_11B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_11B,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_11C",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_11C=fopen(output,"w");
	if (file_raw_11C==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_11C,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_11E",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_11E=fopen(output,"w");
	if (file_raw_11E==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_11E,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_11F",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_11F=fopen(output,"w");
	if (file_raw_11F==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_11F,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_11W",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_11W=fopen(output,"w");
	if (file_raw_11W==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_11W,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_12A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_12A=fopen(output,"w");
	if (file_raw_12A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_12A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_12B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_12B=fopen(output,"w");
	if (file_raw_12B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_12B,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_12D",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_12D=fopen(output,"w");
	if (file_raw_12D==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_12D,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_12E",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_12E=fopen(output,"w");
	if (file_raw_12E==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_12E,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_12K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_12K=fopen(output,"w");
	if (file_raw_12K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_12K,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_13A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_13A=fopen(output,"w");
	if (file_raw_13A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_13A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_13B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_13B=fopen(output,"w");
	if (file_raw_13B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_13B,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_13K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_13K=fopen(output,"w");
	if (file_raw_13K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_13K,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_FCC",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_FCC=fopen(output,"w");
	if (file_raw_FCC==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_FCC,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_HCP",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_HCP=fopen(output,"w");
	if (file_raw_HCP==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_HCP,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_BCC_9",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_BCC_9=fopen(output,"w");
	if (file_raw_BCC_9==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_BCC_9,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_BCC_15",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	file_raw_BCC_15=fopen(output,"w");
	if (file_raw_BCC_15==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(file_raw_BCC_15,"%s\n",output);
}

void Write_Raw_Xmol(int f, FILE *thefile, char *sarr) {
	int i;
	
	fprintf(thefile,"%d\nframe %d of %d\n",N,f,TOTALFRAMES);
	for(i=0; i<N; i++) {
		if (sarr[i]!='C') {
			if (rtype[i]==1) fprintf(thefile,"C\n");
			else fprintf(thefile,"D\n");
		}
		else if (sarr[i]=='C') {
			if (rtype[i]==1) fprintf(thefile,"A\n");
			else fprintf(thefile,"B\n");
		}
	}
}

void Write_11A_cen_xmol(int f) {
	int i;
	
	int n11A_cen=0;
	
	for(i=0; i<N; i++) if(s11A_cen[i]=='O') ++n11A_cen;         // get total number of 13A centres 
		
	fprintf(file_11A_cen_xmol,"%d\nframe %d of %d\n",n11A_cen,f,TOTALFRAMES);
	for(i=0; i<N; i++) {
		if (s11A_cen[i]=='O') fprintf(file_11A_cen_xmol ,"O\t%.5lg\t%.5lg\t%.5lg\n", x[i], y[i], z[i]);
	}
}


void Write_13A_cen_xmol(int f) {
	int i;
	
	int n13A_cen=0;
	
	for(i=0; i<N; i++) if(s13A_cen[i]=='O') ++n13A_cen;         // get total number of 13A centres 
		
	fprintf(file_13A_cen_xmol,"%d\nframe %d of %d\n",n13A_cen,f,FRAMES);
	for(i=0; i<N; i++) {
		if (s13A_cen[i]=='O') fprintf(file_13A_cen_xmol ,"O\t%.5lg\t%.5lg\t%.5lg\n", x[i], y[i], z[i]);
	}
}

        

void Write_Raw(int f) {
	Write_Raw_Xmol(f,file_raw_sp3,&ssp3[0]);
	Write_Raw_Xmol(f,file_raw_sp3a,&ssp3a[0]);
	Write_Raw_Xmol(f,file_raw_sp3b,&ssp3b[0]);
	Write_Raw_Xmol(f,file_raw_5A,&s5A[0]);
	Write_Raw_Xmol(f,file_raw_sp4,&ssp4[0]);
	Write_Raw_Xmol(f,file_raw_sp4a,&ssp4a[0]);
	Write_Raw_Xmol(f,file_raw_sp4b,&ssp4b[0]);
	Write_Raw_Xmol(f,file_raw_6A,&s6A[0]);
	Write_Raw_Xmol(f,file_raw_6Z,&s6Z[0]);
	Write_Raw_Xmol(f,file_raw_7K,&s7K[0]);
	Write_Raw_Xmol(f,file_raw_sp5,&ssp5[0]);
	Write_Raw_Xmol(f,file_raw_sp5a,&ssp5a[0]);
	Write_Raw_Xmol(f,file_raw_sp5b,&ssp5b[0]);
	Write_Raw_Xmol(f,file_raw_7A,&s7A[0]);
	Write_Raw_Xmol(f,file_raw_8A,&s8A[0]);
	Write_Raw_Xmol(f,file_raw_8B,&s8B[0]);
	Write_Raw_Xmol(f,file_raw_8K,&s8K[0]);
	Write_Raw_Xmol(f,file_raw_9A,&s9A[0]);
	Write_Raw_Xmol(f,file_raw_9B,&s9B[0]);
	Write_Raw_Xmol(f,file_raw_9K,&s9K[0]);
	Write_Raw_Xmol(f,file_raw_10A,&s10A[0]);
	Write_Raw_Xmol(f,file_raw_10B,&s10B[0]);
	Write_Raw_Xmol(f,file_raw_10K,&s10K[0]);
	Write_Raw_Xmol(f,file_raw_10W,&s10W[0]);
	Write_Raw_Xmol(f,file_raw_11A,&s11A[0]);
	Write_Raw_Xmol(f,file_raw_11B,&s11B[0]);
	Write_Raw_Xmol(f,file_raw_11C,&s11C[0]);
	Write_Raw_Xmol(f,file_raw_11E,&s11E[0]);
	Write_Raw_Xmol(f,file_raw_11F,&s11F[0]);
	Write_Raw_Xmol(f,file_raw_11W,&s11W[0]);
	Write_Raw_Xmol(f,file_raw_12A,&s12A[0]);
	Write_Raw_Xmol(f,file_raw_12B,&s12B[0]);
	Write_Raw_Xmol(f,file_raw_12D,&s12D[0]);
	Write_Raw_Xmol(f,file_raw_12E,&s12E[0]);
	Write_Raw_Xmol(f,file_raw_12K,&s12K[0]);
	Write_Raw_Xmol(f,file_raw_13A,&s13A[0]);
	Write_Raw_Xmol(f,file_raw_13B,&s13B[0]);
	Write_Raw_Xmol(f,file_raw_13K,&s13K[0]);
	Write_Raw_Xmol(f,file_raw_FCC,&sFCC[0]);
	Write_Raw_Xmol(f,file_raw_HCP,&sHCP[0]);
	Write_Raw_Xmol(f,file_raw_BCC_9,&sBCC_9[0]);
	Write_Raw_Xmol(f,file_raw_BCC_15,&sBCC_15[0]);
}

void Write_Raw_Close() {
	fclose(file_raw_sp3);
	fclose(file_raw_sp3a);
	fclose(file_raw_sp3b);
	fclose(file_raw_5A);
	fclose(file_raw_sp4);
	fclose(file_raw_sp4a);
	fclose(file_raw_sp4b);
	fclose(file_raw_6A);
	fclose(file_raw_6Z);
	fclose(file_raw_7K);
	fclose(file_raw_sp5);
	fclose(file_raw_sp5a);
	fclose(file_raw_sp5b);
	fclose(file_raw_7A);
	fclose(file_raw_8A);
	fclose(file_raw_8B);
	fclose(file_raw_8K);
	fclose(file_raw_9A);
	fclose(file_raw_9B);
	fclose(file_raw_9K);
	fclose(file_raw_10A);
	fclose(file_raw_10B);
	fclose(file_raw_10K);
	fclose(file_raw_10W);
	fclose(file_raw_11A);
	fclose(file_raw_11B);
	fclose(file_raw_11C);
	fclose(file_raw_11E);
	fclose(file_raw_11F);
	fclose(file_raw_11W);
	fclose(file_raw_12A);
	fclose(file_raw_12B);
	fclose(file_raw_12D);
	fclose(file_raw_12E);
	fclose(file_raw_12K);
	fclose(file_raw_13A);
	fclose(file_raw_13B);
	fclose(file_raw_13K);
	fclose(file_raw_FCC);
	fclose(file_raw_HCP);
	fclose(file_raw_BCC_9);
	fclose(file_raw_BCC_15);
}


void Write_Cluster_Init() {
	char errMsg[1000];
	char output[1000];
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp3",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	wsp3=fopen(output,"w");
	if (wsp3==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(wsp3,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp3a",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	wsp3a=fopen(output,"w");
	if (wsp3a==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(wsp3a,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp3b",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	wsp3b=fopen(output,"w");
	if (wsp3b==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(wsp3b,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_5A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w5A=fopen(output,"w");
	if (w5A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w5A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp4",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	wsp4=fopen(output,"w");
	if (wsp4==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(wsp4,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp4a",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	wsp4a=fopen(output,"w");
	if (wsp4a==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(wsp4a,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp4b",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	wsp4b=fopen(output,"w");
	if (wsp4b==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(wsp4b,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp4c",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	wsp4c=fopen(output,"w");
	if (wsp4c==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(wsp4c,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_6A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w6A=fopen(output,"w");
	if (w6A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w6A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_6Z",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w6Z=fopen(output,"w");
	if (w6Z==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w6Z,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp5",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	wsp5=fopen(output,"w");
	if (wsp5==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(wsp5,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp5a",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	wsp5a=fopen(output,"w");
	if (wsp5a==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(wsp5a,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp5b",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	wsp5b=fopen(output,"w");
	if (wsp5b==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(wsp5b,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_7A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w7A=fopen(output,"w");
	if (w7A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w7A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_7K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w7K=fopen(output,"w");
	if (w7K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w7K,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_8A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w8A=fopen(output,"w");
	if (w8A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w8A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_8B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w8B=fopen(output,"w");
	if (w8B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w8B,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_8K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w8K=fopen(output,"w");
	if (w8K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w8K,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_9A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w9A=fopen(output,"w");
	if (w9A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w9A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_9B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w9B=fopen(output,"w");
	if (w9B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w9B,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_9K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w9K=fopen(output,"w");
	if (w9K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w9K,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_10A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w10A=fopen(output,"w");
	if (w10A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w10A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_10B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w10B=fopen(output,"w");
	if (w10B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w10B,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_10K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w10K=fopen(output,"w");
	if (w10K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w10K,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_10W",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w10W=fopen(output,"w");
	if (w10W==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w10W,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_11A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w11A=fopen(output,"w");
	if (w11A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w11A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_11B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w11B=fopen(output,"w");
	if (w11B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w11B,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_11C",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w11C=fopen(output,"w");
	if (w11C==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w11C,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_11E",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w11E=fopen(output,"w");
	if (w11E==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w11E,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_11F",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w11F=fopen(output,"w");
	if (w11F==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w11F,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_11W",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w11W=fopen(output,"w");
	if (w11W==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w11W,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_12A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w12A=fopen(output,"w");
	if (w12A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w12A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_12B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w12B=fopen(output,"w");
	if (w12B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w12B,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_12D",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w12D=fopen(output,"w");
	if (w12D==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w12D,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_12E",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w12E=fopen(output,"w");
	if (w12E==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w12E,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_12K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w12K=fopen(output,"w");
	if (w12K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w12K,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_13A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w13A=fopen(output,"w");
	if (w13A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w13A,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_13B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w13B=fopen(output,"w");
	if (w13B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w13B,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_13K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	w13K=fopen(output,"w");
	if (w13K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(w13K,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_FCC",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	wFCC=fopen(output,"w");
	if (wFCC==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(wFCC,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_HCP",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	wHCP=fopen(output,"w");
	if (wHCP==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(wHCP,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_BCC_9",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	wBCC_9=fopen(output,"w");
	if (wBCC_9==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(wBCC_9,"%s\n",output);
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_BCC_15",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	wBCC_15=fopen(output,"w");
	if (wBCC_15==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
	fprintf(wBCC_15,"%s\n",output);
}

void Write_Cluster_Close() {
	fclose(wsp3);
	fclose(wsp3a);
	fclose(wsp3b);
	fclose(w5A);
	fclose(wsp4);
	fclose(wsp4a);
	fclose(wsp4b);
	fclose(wsp4c);
	fclose(w6A);
	fclose(w6Z);
	fclose(w7K);
	fclose(wsp5);
	fclose(wsp5a);
	fclose(wsp5b);
	fclose(w7A);
	fclose(w8A);
	fclose(w8B);
	fclose(w8K);
	fclose(w9A);
	fclose(w9B);
	fclose(w9K);
	fclose(w10A);
	fclose(w10B);
	fclose(w10K);
	fclose(w10W);
	fclose(w11A);
	fclose(w11B);
	fclose(w11C);
	fclose(w11E);
	fclose(w11F);
	fclose(w11W);
	fclose(w12A);
	fclose(w12B);
	fclose(w12D);
	fclose(w12E);
	fclose(w12K);
	fclose(w13A);
	fclose(w13B);
	fclose(w13K);
	fclose(wFCC);
	fclose(wHCP);
	fclose(wBCC_9);
	fclose(wBCC_15);
}

void Write_Cluster(int f, FILE *writeout, int *n, int **hc, int clusSize) {
	int i,j;
	
	fprintf(writeout,"%d\n",n[f]);
	for (i=0;i<n[f];i++) {
		fprintf(writeout,"%d",hc[i][0]);
		for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",hc[i][j]);
		fprintf(writeout,"	%d\n",hc[i][clusSize-1]);
	}
}


void Write_Cluster_sp3(int f, FILE *writeout) {
	int i,j;
	int clusSize=3;
	
	fprintf(writeout,"%d\n",nsp3[f]);
	for (i=0;i<nsp3a[f];i++) {
		fprintf(writeout,"%d",sp3a[i][0]);
		for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp3a[i][j]);
		fprintf(writeout,"	%d\n",sp3a[i][clusSize-1]);
	}
	for (i=0;i<nsp3b[f];i++) {
		fprintf(writeout,"%d",sp3b[i][0]);
		for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp3b[i][j]);
		fprintf(writeout,"	%d\n",sp3b[i][clusSize-1]);
	}
	for (i=0;i<nsp3c[f];i++) {
		fprintf(writeout,"%d",sp3c[i][0]);
		for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp3c[i][j]);
		fprintf(writeout,"	%d\n",sp3c[i][clusSize-1]);
	}
}

void Write_Cluster_sp4(int f, FILE *writeout) {
	int i,j;
	int clusSize=4;
	
	fprintf(writeout,"%d\n",nsp4[f]);
	for (i=0;i<nsp4a[f];i++) {
		fprintf(writeout,"%d",sp4a[i][0]);
		for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp4a[i][j]);
		fprintf(writeout,"	%d\n",sp4a[i][clusSize-1]);
	}
	for (i=0;i<nsp4b[f];i++) {
		fprintf(writeout,"%d",sp4b[i][0]);
		for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp4b[i][j]);
		fprintf(writeout,"	%d\n",sp4b[i][clusSize-1]);
	}
	for (i=0;i<nsp4c[f];i++) {
		fprintf(writeout,"%d",sp4c[i][0]);
		for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp4c[i][j]);
		fprintf(writeout,"	%d\n",sp4c[i][clusSize-1]);
	}
}

void Write_Cluster_sp5(int f, FILE *writeout) {
	int i,j;
	int clusSize=5;
	
	fprintf(writeout,"%d\n",nsp5[f]);
	for (i=0;i<nsp5a[f];i++) {
		fprintf(writeout,"%d",sp5a[i][0]);
		for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp5a[i][j]);
		fprintf(writeout,"	%d\n",sp5a[i][clusSize-1]);
	}
	for (i=0;i<nsp5b[f];i++) {
		fprintf(writeout,"%d",sp5b[i][0]);
		for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp5b[i][j]);
		fprintf(writeout,"	%d\n",sp5b[i][clusSize-1]);
	}
	for (i=0;i<nsp5c[f];i++) {
		fprintf(writeout,"%d",sp5c[i][0]);
		for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp5c[i][j]);
		fprintf(writeout,"	%d\n",sp5c[i][clusSize-1]);
	}
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

void Norm_Write_bonded_to_cen_distro(char *output, int *the_array, int clusSize, double normFactor) {
	FILE *fOut;
	int i;
	char errMsg[1000];
	double tempratio;
	
	fOut=fopen(output,"w");
	if (fOut==NULL)  {
		sprintf(errMsg,"Norm_Write_bonded_to_cen_distro(): Error opening file %s",output);	// Always test file open
		Error(errMsg);
	}
	fprintf(fOut,"%s\n",output);
	
	fprintf(fOut,"N_A in clus	N_B in clus	Freq	Norm\n");
	for (i=0; i<nB+1; i++) {
		tempratio=(double)the_array[i]/(normFactor);
		fprintf(fOut,"%d	%d %d	%.15lg\n",i,clusSize-i,the_array[i],tempratio);
	}
	
	fclose(fOut);
	printf("Written %s\n",output);
}

//// START: main() routine
int main(int argc, char **argv) {
	int e, f, i;
	int write, remainder;
	char errMsg[1000], output[1000], other[1000];
	int ix, iy, iz;
	int imap;
	FILE *rXmol;
	FILE *rSizes;

	sprintf(fInputParamsName,"inputparameters.ini");
	Setup_ReadIniFile(fInputParamsName);	// read input params
	printf("box size file: %s\n",fBoxSizeName);

	if (ISNOTCUBIC!=0){
		if (USELIST==1) {
		sprintf(errMsg,"main() : Error! Need switch cell list off for non-cubic/NPT system");	// Always test file open
		Error(errMsg);
		}
	}
	//read in box data if noncubic/NPT
	if  (ISNOTCUBIC!=0){
		printf("reading box size data from %s\n",fBoxSizeName);
		rSizes=fopen(fBoxSizeName,"r");
		if(rSizes==NULL)  {
			sprintf(errMsg,"main() : Error opening boxfile %s",fBoxSizeName);
			Error_no_free(errMsg);
		}
		fgets(other,1000,rSizes); //reads first line
		if  (ISNOTCUBIC==1) {
			Setup_ReadBox(rSizes);
			printf("sidex: %f, sidey: %f, sidez: %f\n", sidex,sidey,sidez);
		}
		if  (ISNOTCUBIC==3) {
            printf("======> Triclinic Box \n");
        }
    }

	printf("reading coordinate frames from %s\n\n",fXmolName);
	rXmol=fopen(fXmolName,"r");	// open xmol trajecotry

	if (rXmol==NULL)  {
		sprintf(errMsg,"main() : Error opening file %s",fXmolName);	// Always test file open
		Error_no_free(errMsg);
	}
	
	if (doWriteBonds==1) {
		sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.bonds",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
		bondsout=fopen(output, "w");
		if (bondsout==NULL)  {
			sprintf(errMsg,"main() : Error opening file %s",output);	// Always test file open
			Error_no_free(errMsg);
		}
	}
	
	if (USELIST==1) {
		M = (int)(side/rcutAA);	// number of cells along box side
		if (M<3) Error_no_free("main(): M<3, too few cells");
		ncells = M*M*M;	// total number of cells
	}
	
	Setup_InitStaticVars();
	
	if (USELIST==1) {
		cellSide = side/M;	// length of cells
		invcellSide = 1.0/cellSide;	// invcellSide
		printf("m %d ncells %d cellside %.15lg\n", M, ncells, cellSide);
		// routine to create the thirteen nearest neighbours array map[] of each cell 
		for (iz=1; iz<=M; iz++) {
			for (iy=1; iy<=M; iy++) {
				for (ix=1; ix<=M; ix++) {
					imap = (icell(ix,iy,iz)-1)*13;
					map[imap+1 ]=icell(ix+1,iy	,iz	);
					map[imap+2 ]=icell(ix+1,iy+1,iz	);
					map[imap+3 ]=icell(ix	 ,iy+1,iz	);
					map[imap+4 ]=icell(ix-1 ,iy+1,iz	);
					map[imap+5 ]=icell(ix+1,iy	,iz-1	);
					map[imap+6 ]=icell(ix+1,iy+1,iz-1	);
					map[imap+7 ]=icell(ix	 ,iy+1,iz-1	);
					map[imap+8 ]=icell(ix-1 ,iy+1,iz-1	);
					map[imap+9 ]=icell(ix+1,iy	,iz+1	);
					map[imap+10]=icell(ix+1,iy+1,iz+1	);
					map[imap+11]=icell(ix	 ,iy+1,iz+1	);
					map[imap+12]=icell(ix-1 ,iy+1,iz+1);
					map[imap+13]=icell(ix	 ,iy	,iz+1	);
				}
			}
		}
	}
	
	printf("initializing static variables...");
	Stats_Init();
	printf("completed\n");
	
	if (doWriteClus==1) {
		printf("\ninitializing cluster files...");
		Write_Cluster_Init();
		printf("completed\n");
	}
	
	if (doWriteRaw==1) {
		printf("\ninitializing raw cluster xmol files...");
		Write_Raw_Init();
		printf("completed\n");
	}
	
	if (do11AcenXmol==1) {
		printf("\ninitializing 11A centre xmol files...");
		sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_11A_cen.xmol",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
		file_11A_cen_xmol=fopen(output, "w");
		printf("completed\n");
	}
	
	if (do13AcenXmol==1) {
		printf("\ninitializing 13A centre xmol files...");
		sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_13A_cen.xmol",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
		file_13A_cen_xmol=fopen(output, "w");
		printf("completed\n");
	}

	printf("begin main loop\n");

	f=0;
	for (e=0;e<TOTALFRAMES;e++) {
		remainder=e%SAMPLEFREQ;
		if (remainder==0 && f<FRAMES && e>=STARTFROM) {
			write=1;
		}
		else write=0;
		if (write==1) Setup_ResetStaticVars(f);

		if (ISNOTCUBIC>=2) {
			Setup_ReadBox(rSizes);
		}
		Setup_Readxyz(e,write,f,rXmol);
		
		if (write==1) {
			Bonds_GetBonds(f);

			for(i=0; i<N; i++) {
				if (cnb[i]>maxnb) maxnb=cnb[i];
				if (dosp3==1) Rings_gSP3(f,i);
			}
			if (dosp3==1) Rings_setSP3c(f);			
			if (dosp4==1) Rings_setSP4c(f);
			if (dosp5==1) Rings_setSP5c(f);
			if (do6Z==1) Clusters_Get6Z_C2v(f);
			if (do7K==1) Clusters_Get7K(f);
			if (do8A==1) Clusters_Get8A_D2d(f);
			if (do8B==1) Clusters_Get8B_Cs(f);
			if (do8K==1) Clusters_Get8K(f);
			if (do9A==1) Clusters_Get9A_D3h(f);
			if (do9B==1) Clusters_Get9B_10B_11B_11E_12D(f);
			if (do9K==1) Clusters_Get9K_10K(f);
			if (do10A==1) Clusters_Get10A_C3v(f);
			if (do10W==1) Clusters_Get10W(f);
			if (do11A==1) Clusters_Get11A_12K(f);
			if (do11C==1) Clusters_Get11C_12A(f);
			if (do11F==1) Clusters_Get11F_12E_13K(f);
			if (do12B==1) Clusters_Get12B_13A(f);
			if (do13B==1) Clusters_Get13B_D5h(f);
			if (doFCC==1) Clusters_GetFCC(f);
			if (doHCP==1) Clusters_GetHCP(f);
			if (doBCC9==1) Clusters_GetBCC_9(f);
			if (doBCC15==1) Clusters_GetBCC_15(f);

			if (doWriteClus==1) {
				Write_Cluster_sp3(f,wsp3);
				Write_Cluster(f,wsp3a,nsp3a,sp3a,3);
				Write_Cluster(f,wsp3b,nsp3b,sp3b,4);
				Write_Cluster(f,w5A,nsp3c,sp3c,5);
				Write_Cluster_sp4(f,wsp4);
				Write_Cluster(f,wsp4a,nsp4a,sp4a,4);
				Write_Cluster(f,wsp4b,nsp4b,sp4b,5);
				Write_Cluster(f,wsp4c,nsp4c,sp4c,6);
				Write_Cluster(f,w6A,n6A,hc6A,6);
				Write_Cluster_sp5(f,wsp5);
				Write_Cluster(f,wsp5a,nsp5a,sp5a,5);
				Write_Cluster(f,wsp5b,nsp5b,sp5b,6);
				Write_Cluster(f,w7A,nsp5c,sp5c,7);
				Write_Cluster(f,w6Z,n6Z,hc6Z,6);
				Write_Cluster(f,w7K,n7K,hc7K,7);
				Write_Cluster(f,w8A,n8A,hc8A,8);
				Write_Cluster(f,w8B,n8B,hc8B,8);
				Write_Cluster(f,w8K,n8K,hc8K,8);
				Write_Cluster(f,w9A,n9A,hc9A,9);
				Write_Cluster(f,w9B,n9B,hc9B,9);
				Write_Cluster(f,w9K,n9K,hc9K,9);
				Write_Cluster(f,w10A,n10A,hc10A,10);
				Write_Cluster(f,w10B,n10B,hc10B,10);
				Write_Cluster(f,w10K,n10K,hc10K,10);
				Write_Cluster(f,w10W,n10W,hc10W,10);
				Write_Cluster(f,w11A,n11A,hc11A,11);
				Write_Cluster(f,w11B,n11B,hc11B,11);
				Write_Cluster(f,w11C,n11C,hc11C,11);
				Write_Cluster(f,w11E,n11E,hc11E,11);
				Write_Cluster(f,w11F,n11F,hc11F,11);
				Write_Cluster(f,w11W,n11W,hc11W,11);
				Write_Cluster(f,w12A,n12A,hc12A,12);
				Write_Cluster(f,w12B,n12B,hc12B,12);
				Write_Cluster(f,w12D,n12D,hc12D,12);
				Write_Cluster(f,w12E,n12E,hc12E,12);
				Write_Cluster(f,w12K,n12K,hc12K,12);
				Write_Cluster(f,w13A,n13A,hc13A,13);
				Write_Cluster(f,w13B,n13B,hc13B,13);
				Write_Cluster(f,w13K,n13K,hc13K,13);
				Write_Cluster(f,wFCC,nFCC,hcFCC,13);
				Write_Cluster(f,wHCP,nHCP,hcHCP,13);
				Write_Cluster(f,wBCC_9,nBCC_9,hcBCC_9,9);
				Write_Cluster(f,wBCC_15,nBCC_15,hcBCC_15,15);
			}
			if (PRINTINFO==1) printf("\n");

			if (doWriteRaw==1) Write_Raw(f);
			
			if (do11AcenXmol==1) Write_11A_cen_xmol(f);
			if (do13AcenXmol==1) Write_13A_cen_xmol(f);
		    
			
			Stats_Reset();
			Stats_Analyse();
			Pop_Per_Frame(f);

			printf("f%d complete\n",f);
			f++;
		}
		if (f==FRAMES) break;
	}


	fclose(rXmol);
	if (doWriteBonds==1) fclose(bondsout);
	
	if (f!=FRAMES) {
		printf("\n\n\n!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!\n\n\n");
		printf("Analysed frames %d less than expected number of FRAMES %d from %s\n\n",f,FRAMES,fInputParamsName);
	}
	
	if (doWriteClus==1) Write_Cluster_Close();
	
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
	
	if (doWritePopPerFrame==1) {
		printf("Writing pop_per_frame %s ....",output);
		sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.pop_per_frame",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
		fPopPerFrame=fopen(output, "w");
		if (fPopPerFrame==NULL)  {
			sprintf(errMsg,"main() : Error opening file %s",output);	// Always test file open
			Error(errMsg);
		}
		fprintf(fPopPerFrame,"%s\n",output);
		
		fprintf(fPopPerFrame,"frame	time	time_norm_t_a	sp3	sp3a	sp3b	5A	sp4	sp4a	sp4b	6A	6Z	sp5	sp5a	sp5b	7A	7K	8A	8B	8K	9A	9B	9K	10A	10B	10K	10W");
		fprintf(fPopPerFrame,"	11A	11B	11C	11E	11F	11W	12A	12B	12D	12E	12K	13A	13B	13K	FCC	HCP	BCC_9	BCC_15\n");
		
		
		fprintf(fPopPerFrame,"mean	-	-	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg",mean_pop_per_frame_sp3,mean_pop_per_frame_sp3a,mean_pop_per_frame_sp3b,mean_pop_per_frame_sp3c,mean_pop_per_frame_sp4,mean_pop_per_frame_sp4a,mean_pop_per_frame_sp4b,mean_pop_per_frame_6A,mean_pop_per_frame_6Z,mean_pop_per_frame_sp5,mean_pop_per_frame_sp5a,mean_pop_per_frame_sp5b,mean_pop_per_frame_sp5c,mean_pop_per_frame_7K,mean_pop_per_frame_8A,mean_pop_per_frame_8B,mean_pop_per_frame_8K,mean_pop_per_frame_9A,mean_pop_per_frame_9B,mean_pop_per_frame_9K,mean_pop_per_frame_10A,mean_pop_per_frame_10B,mean_pop_per_frame_10K,mean_pop_per_frame_10W);
		fprintf(fPopPerFrame,"	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg\n",mean_pop_per_frame_11A,mean_pop_per_frame_11B,mean_pop_per_frame_11C,mean_pop_per_frame_11E,mean_pop_per_frame_11F,mean_pop_per_frame_11W,mean_pop_per_frame_12A,mean_pop_per_frame_12B,mean_pop_per_frame_12D,mean_pop_per_frame_12E,mean_pop_per_frame_12K,mean_pop_per_frame_13A,mean_pop_per_frame_13B,mean_pop_per_frame_13K,mean_pop_per_frame_FCC,mean_pop_per_frame_HCP,mean_pop_per_frame_BCC_9,mean_pop_per_frame_BCC_15);
		for (f=0;f<FRAMES;f++) {
			fprintf(fPopPerFrame,"%d	%.15lg	%.15lg",f,(double)f*FRAMETSTEP*SAMPLEFREQ,(double)f*FRAMETSTEP*SAMPLEFREQ/talpha);
			fprintf(fPopPerFrame,"	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg",pop_per_frame_sp3[f],pop_per_frame_sp3a[f],pop_per_frame_sp3b[f],pop_per_frame_sp3c[f],pop_per_frame_sp4[f],pop_per_frame_sp4a[f],pop_per_frame_sp4b[f],pop_per_frame_6A[f],pop_per_frame_6Z[f],pop_per_frame_sp5[f],pop_per_frame_sp5a[f],pop_per_frame_sp5b[f],pop_per_frame_sp5c[f],pop_per_frame_7K[f],pop_per_frame_8A[f],pop_per_frame_8B[f],pop_per_frame_8K[f],pop_per_frame_9A[f],pop_per_frame_9B[f],pop_per_frame_9K[f],pop_per_frame_10A[f],pop_per_frame_10B[f],pop_per_frame_10K[f],pop_per_frame_10W[f]);
			fprintf(fPopPerFrame,"	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg\n",pop_per_frame_11A[f],pop_per_frame_11B[f],pop_per_frame_11C[f],pop_per_frame_11E[f],pop_per_frame_11F[f],pop_per_frame_11W[f],pop_per_frame_12A[f],pop_per_frame_12B[f],pop_per_frame_12D[f],pop_per_frame_12E[f],pop_per_frame_12K[f],pop_per_frame_13A[f],pop_per_frame_13B[f],pop_per_frame_13K[f],pop_per_frame_FCC[f],pop_per_frame_HCP[f],pop_per_frame_BCC_9[f],pop_per_frame_BCC_15[f]);
		}
		fclose(fPopPerFrame);
		printf("Closed file %s\n\n",output);
	}	
	
	if (doWriteRaw==1) {
		printf("Closing raw cluster xmol files....");
		Write_Raw_Close();
		printf("closed!\n\n");
	}
	
	if (do11AcenXmol==1) {
		printf("Closing 11A centre xmol files....");
		fclose(file_11A_cen_xmol);
		printf("closed!\n\n");
	}
	
	if (do13AcenXmol==1) {
		printf("=Closing 13A centre xmol files....");
		fclose(file_13A_cen_xmol);
		printf("closed!\n\n");
	}
	
	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.static_clust",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	printf("\n");
	Stats_Report(output);
	printf("\nWritten %s\n\n",output);

	
	Setup_FreeStaticVars();
	Stats_FreeMem();
	if (ISNOTCUBIC > 0) {
        fclose(rSizes);
    }
	printf("\n\nFIN \n\n");
	return 0;
}
