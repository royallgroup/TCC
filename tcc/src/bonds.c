#include "math.h"
#include "bonds.h"
#include "globals.h"
#include "tools.h"
#include "cell_list.h"
#include "string.h"
#include "voronoi_bonds.h"

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
        if (USELIST==0) {
            Get_Bonds_With_Voronoi();
        }
        else {
            Get_Bonds_With_Voronoi_And_Cell_List();
        }
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
    int particle_1, particle_2;
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

            Check_For_Valid_Bond(particle_1, particle_2, squared_distance);
        }
    }
}

void Check_For_Valid_Bond(int particle_1, int particle_2, double squared_distance) {
    if (squared_distance < rcutAA2 && particle_type[particle_1] == 1 && particle_type[particle_2] == 1){
        Check_Num_Bonds(particle_1, particle_2, squared_distance);
    }
    else if (squared_distance < rcutBB2 && particle_type[particle_1]==2 && particle_type[particle_2]==2){
        Check_Num_Bonds(particle_1, particle_2, squared_distance);
    }
    else if (squared_distance < rcutAB2) {
        Check_Num_Bonds(particle_1, particle_2, squared_distance);
    }
}

void Check_Num_Bonds(int particle_1, int particle_2, double squared_distance) {
    if (cnb[particle_1] < nB && cnb[particle_2] < nB){
        Add_New_Bond(particle_1, particle_2, squared_distance);
    }
    else {
        Too_Many_Bonds(particle_1, particle_2);
    }
}

void Too_Many_Bonds(int particle_1, int particle_2) {
    char error_message[200];

    sprintf(error_message, "Get_Simple_Bonds(): Too many bonds to particle %d or particle_2 %d.\n"
            "This is probably because rcutAA or rcutBB is too large\n", particle_1,particle_2);
    Error(error_message);
}

void Add_New_Bond(int particle_1, int particle_2, double squared_distance) {
    int k;

    k = cnb[particle_1]++;
    bNums[particle_1][k] = particle_2;
    bondlengths[particle_1][k]=sqrt(squared_distance);
    k = cnb[particle_2]++;
    bNums[particle_2][k] = particle_1;
    bondlengths[particle_2][k]=sqrt(squared_distance);
}

int Bonds_BondCheck(int i, int j) { // Returns 1 if i & j are bonded; 0 otherwise
    int k;

    for (k=0; k<cnb[i]; ++k) {
        if (bNums[i][k] == j) return 1;
    } 
    return 0;
}
