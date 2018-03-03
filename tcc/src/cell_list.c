#include "cell_list.h"
#include "globals.h"
#include "tools.h"
#include "bonds.h"
#include "math.h"

void Get_Bonds_With_Voronoi_And_Cell_List() {  // Get bonds using Voronoi
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

    store_dr2 = malloc(current_frame_particle_number*sizeof(double));   if (store_dr2==NULL) { sprintf(errMsg,"Get_Bonds_With_Voronoi_And_Cell_List(): store_dr2[] malloc out of memory\n");    Error(errMsg); }
    temp_cnb = malloc(current_frame_particle_number*sizeof(int));   if (temp_cnb==NULL) { sprintf(errMsg,"Get_Bonds_With_Voronoi_And_Cell_List(): temp_cnb[] malloc out of memory\n");  Error(errMsg); }
    temp_bNums = malloc(current_frame_particle_number*sizeof(int *));   if (temp_bNums==NULL) { sprintf(errMsg,"Get_Bonds_With_Voronoi_And_Cell_List(): temp_bNums[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<current_frame_particle_number; ++j) { temp_bNums[j] = malloc(nBs*sizeof(int));  if (temp_bNums[j]==NULL) { sprintf(errMsg,"Get_Bonds_With_Voronoi_And_Cell_List(): temp_bNums[][] malloc out of memory\n"); Error_no_free(errMsg); } }

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
                dr2 = Get_Interparticle_Distance(i - 1, j - 1);
                if (dr2 < rcutAA2) {
                    if (temp_cnb[i-1] < nBs && temp_cnb[j-1] < nBs) {  // max number of bonds, do ith particle
                        temp_bNums[i-1][temp_cnb[i-1]]=j-1;
                        temp_bNums[j-1][temp_cnb[j-1]]=i-1;
                        temp_cnb[i-1]++;
                        temp_cnb[j-1]++;
                    }
                    else {    // list is now full
                        printf("Get_Bonds_With_Voronoi_And_Cell_List(): nBs %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",nBs,i-1,j-1);
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
                    dr2 = Get_Interparticle_Distance(i - 1, j - 1);
                    if (dr2 < rcutAA2) {
                        if (temp_cnb[i-1] < nBs && temp_cnb[j-1] < nBs) {  // max number of bonds, do ith particle
                            temp_bNums[i-1][temp_cnb[i-1]]=j-1;
                            temp_bNums[j-1][temp_cnb[j-1]]=i-1;
                            temp_cnb[i-1]++;
                            temp_cnb[j-1]++;
                        }
                        else {    // list is now full
                            printf("Get_Bonds_With_Voronoi_And_Cell_List(): nBs %d number of bonds per particle is not big enough: particle i %d or j% d has too many bonds\nThis is probably because rcutAA is too large\n",nBs,i-1,j-1);
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
            dr2 = Get_Interparticle_Distance(i, temp_bNums[i][j]);
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
            printf("Get_Bonds_With_Voronoi_And_Cell_List(): part %d - cnbs %d does not equal cnbs2 %d \n",i,cnbs,cnbs2);
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
                    printf("Get_Bonds_With_Voronoi_And_Cell_List(): nB %d number of bonds per particle is not big enough: particle i %d cnb[i] %d or j %d cnb[j] %d has too many bonds\nThis is probably because rcutAA is too large\n",nB,i,cnb[i],j,cnb[j]);
                    exit(1);
                }
            }
        }
        if (PRINTINFO==1) if (!((i+1)%10000)) printf("Get_Bonds_With_Voronoi_And_Cell_List(): particle %d of %d done\n",i+1,current_frame_particle_number);
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

void links() {  // sorts all the particles into cells, result given by head-of-chain and linked list arrays
    int i, ic;
    for (ic=1;ic<=n_cells_total;ic++) head[ic]=0;
    for (i=1;i<=current_frame_particle_number;i++) {
        ic = 1 + (int)((x[i-1]+ half_sidex)*inv_cell_len_x)
             + n_cells_y*((int)((y[i-1]+half_sidex)*inv_cell_len_y))
             + n_cells_z*n_cells_z*((int)((z[i-1]+half_sidex)*inv_cell_len_z));
        if (ic > n_cells_total || ic <= 0) {
            printf("i %d r_x %lg r_y %lg r_z %lg side %lg halfSide %lg ic %d ncells %d\n",i-1,x[i-1],y[i-1],z[i-1],sidex,half_sidex,ic,n_cells_total);
            Error("links(): ic > ncells, i.e. particle coord no longer in simulation box!!\n");
        }
        llist[i]=head[ic];
        head[ic]=i;
    }
}

int icell(int tix, int tiy, int tiz) { 	// returns cell number (from 1 to ncells) for given (tix,tiy,tiz) coordinate
    return 1 + (tix - 1 + n_cells_x) % n_cells_x
           + n_cells_y * ((tiy - 1 + n_cells_y) % n_cells_y)
           + n_cells_z * n_cells_z * ((tiz - 1 + n_cells_z) % n_cells_z);
}

void Setup_Cell_List() {

    int i;
    int ix, iy, iz;
    int imap;
    char errMsg[1000];

    n_cells_x = (int)(sidex/rcutAA);
    n_cells_y = (int)(sidey/rcutAA);
    n_cells_z = (int)(sidez/rcutAA);
    if (n_cells_x < 3 || n_cells_y < 3 || n_cells_z < 3) Error_no_free("main(): M<3, too few cells");
    n_cells_total = n_cells_x*n_cells_y*n_cells_z;
    inv_cell_len_x = n_cells_x/sidex;
    inv_cell_len_y = n_cells_y/sidey;
    inv_cell_len_z = n_cells_z/sidez;

    printf("x_cells %d y_cells %d z_cells %d total_cells %d\n", n_cells_x, n_cells_y, n_cells_z, n_cells_total);

    head=malloc((n_cells_total+1)*sizeof(int));    if (head==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): head[] malloc out of memory\n");  Error_no_free(errMsg); }
    map=malloc((13*n_cells_total+1)*sizeof(int));  if (map==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): map[] malloc out of memory\n");    Error_no_free(errMsg); }
    llist=malloc((current_frame_particle_number+1)*sizeof(int)); if (llist==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): llist[] malloc out of memory\n");    Error_no_free(errMsg); }

    for (i=0; i<n_cells_total+1; i++) head[i]=0;
    for (i=0; i<13*n_cells_total+1; i++) map[i]=0;
    for (i=0; i<current_frame_particle_number+1; i++) llist[i]=0;


    // routine to create the thirteen nearest neighbours array map[] of each cell
    for (iz=1; iz<=n_cells_z; iz++) {
        for (iy=1; iy<=n_cells_y; iy++) {
            for (ix=1; ix<=n_cells_x; ix++) {
                imap = (icell(ix,iy,iz)-1)*13;
                map[imap+1 ]=icell(ix+1, iy,   iz);
                map[imap+2 ]=icell(ix+1, iy+1, iz);
                map[imap+3 ]=icell(ix,   iy+1, iz);
                map[imap+4 ]=icell(ix-1, iy+1, iz);
                map[imap+5 ]=icell(ix+1, iy,   iz-1);
                map[imap+6 ]=icell(ix+1, iy+1, iz-1);
                map[imap+7 ]=icell(ix,   iy+1, iz-1);
                map[imap+8 ]=icell(ix-1, iy+1, iz-1);
                map[imap+9 ]=icell(ix+1, iy,   iz+1);
                map[imap+10]=icell(ix+1, iy+1, iz+1);
                map[imap+11]=icell(ix,   iy+1, iz+1);
                map[imap+12]=icell(ix-1, iy+1, iz+1);
                map[imap+13]=icell(ix,   iy,   iz+1);
            }
        }
    }
}