#include "cell_list.h"
#include "voronoi_bonds.h"
#include "globals.h"
#include "tools.h"
#include "bonds.h"

void get_particle_bonds_in_cell_list(const int max_allowed_bonds, int *temp_cnb, int **temp_bNums, int i, int j);

int cell_list_get_particle_1_neighbours(int i, int num_particle_1_neighbours, int *particle_1_neighbours,
                                        int *particle_1_bonds, double *particle_1_bond_lengths, double *store_dr2,
                                        const int *temp_cnb, int *const *temp_bNums) {
    int j;
    double squared_distance;

    num_particle_1_neighbours = 0;
    for (j = 0; j < temp_cnb[i]; ++j) {
        squared_distance = Get_Interparticle_Distance(i, temp_bNums[i][j]);
        particle_1_neighbours[num_particle_1_neighbours] = temp_bNums[i][j];
        particle_1_bonds[num_particle_1_neighbours] = 1;
        particle_1_bond_lengths[num_particle_1_neighbours] = squared_distance;
        store_dr2[temp_bNums[i][j]] = squared_distance;
        num_particle_1_neighbours++;
    }
    return num_particle_1_neighbours;
}

void cell_list_get_neigbours(const int max_allowed_bonds, int *temp_cnb, int **temp_bNums) {
    int i, j;
    int ic, jcell0, jcell, nabor;    // various counters


    links();

    for (ic= 1; ic <= n_cells_total; ic++) {        // loop over all cells
        i = head[ic];     // head of list particle for cell ic
        while (i > 0) {   // loop over all particles in ic
            j = llist[i]; // next particle in current cell ic
            get_particle_bonds_in_cell_list(max_allowed_bonds, temp_cnb, temp_bNums, i, j);
            jcell0 = 13 * (ic - 1);       // now loop over adjacent cells to cell ic
            for (nabor = 1; nabor <= 13; nabor++) {
                jcell = map[jcell0 + nabor];
                j = head[jcell];  // head of cell for jcell
                get_particle_bonds_in_cell_list(max_allowed_bonds, temp_cnb, temp_bNums, i, j);
            }
            i=llist[i]; // next particle in ic cell
        }
    }
}

void get_particle_bonds_in_cell_list(const int max_allowed_bonds, int *temp_cnb, int **temp_bNums, int i, int j) {
    double squared_distance;

    while (j > 0) {   // loop over head of cell and all other particles in jcell
        squared_distance = Get_Interparticle_Distance(i - 1, j - 1);
        if (squared_distance < rcutAA2) {
            if (temp_cnb[i - 1] < max_allowed_bonds &&
                temp_cnb[j - 1] < max_allowed_bonds) {  // max number of bonds, do ith particle
                temp_bNums[i - 1][temp_cnb[i - 1]] = j - 1;
                temp_bNums[j - 1][temp_cnb[j - 1]] = i - 1;
                temp_cnb[i - 1]++;
                temp_cnb[j - 1]++;
            } else {
                Too_Many_Bonds(i - 1, j - 1, __func__);
            }
        }
        j = llist[j]; // next particle in jcell
    }
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
