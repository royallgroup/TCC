/* Alex Malins - alex.malins@gmail.com */
/* TCC: A topological cluster classification code with temporal tracking of clusters. */
/* Not for general consumption */

#include "globals.h"
#include "bonds.h"
#include "rings.h"
#include "clusters.h"
#include "output.h"
#include "stats.h"
#include "tools.h"

char* cluster_names[] = {"sp3", "sp3a", "sp3b", "sp3c", "sp4", "sp4a", "sp4b", "sp4c", "sp5", "sp5a",
                         "sp5b", "sp5c", "6Z", "7K", "8A", "8B", "8K", "9A", "9B", "9K", "10A", "10B",
                         "10K", "10W", "11A", "11B", "11C", "11E", "11F", "11W", "12A", "12B", "12D",
                         "12E", "12K", "13A", "13B", "13K", "FCC", "HCP", "BCC_9", "BCC_15"};
int num_cluster_types = 42;
int* do_cluster_list[] = {&dosp3, &dosp3a, &dosp3b, &dosp3c, &dosp4, &dosp4a, &dosp4b, &dosp4c, &dosp5, &dosp5a,
                        &dosp5b, &dosp5c, &do6Z, &do7K, &do8A, &do8B, &do8K, &do9A, &do9B, &do9K, &do10A, &do10B,
                        &do10K, &do10W, &do11A, &do11B, &do11C, &do11E, &do11F, &do11W, &do12A, &do12B, &do12D,
                        &do12E, &do12K, &do13A, &do13B, &do13K, &doFCC, &doHCP, &doBCC9, &doBCC15};

int main(int argc, char **argv) {
	int e, f, i;
	int write, remainder;
	char errMsg[1000], output[1000], other[1000];
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
	
	Setup_InitStaticVars();
	
    Setup_Cell_List();
	
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

    // Open the files to print out cluster centers
    Setup_Centers_Files();

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
			if (do9K==1) Clusters_Get9K(f);
			if (do10A==1) Clusters_Get10A_C3v(f);
            if (do10K==1) Clusters_Get10K(f);
			if (do10W==1) Clusters_Get10W(f);
			if (do11A==1) Clusters_Get11A(f);
			if (do11C==1) Clusters_Get11C(f);
			if (do11F==1) Clusters_Get11F_12E_13K(f);
			if (do11W==1) Clusters_Get11W(f);
            if (do12A==1) Clusters_Get12A(f);
			if (do12B==1) Clusters_Get12B_13A(f);
			if (do12K==1) Clusters_Get12K(f);
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
	
    Close_Centers_Files();

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
