/* Alex Malins - alex.malins@gmail.com */
/* TCC: A topological cluster classification code with temporal tracking of clusters. */
/* Not for general consumption */

#include "globals.h"
#include "bonds.h"
#include "rings.h"
#include "clusters.h"
#include "output.h"
#include "stats.h"

int icell(int tix, int tiy, int tiz) { 	// returns cell number (from 1 to ncells) for given (tix,tiy,tiz) coordinate
	return 1 + (tix-1+M)%M + M*((tiy-1+M)%M) + M*M*((tiz-1+M)%M); 
}

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

            // Write output files
			if (doWriteClus==1) Write_Cluster(f);
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
	
    Normalise_Populations();
	
	if (doWritePopPerFrame==1) {
        Write_Pop_Per_Frame(f);
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
		printf("Closing 13A centre xmol files....");
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
