#include "potentials.h"

void BLJ() {    // binary Lennard-Jones potential
    double dx,dy,dz,sep2,invsep2,invsep6;   // sep for separations
    double en;
    int i,j;
    
    for (i=0;i<N-1;i++) {
        for (j=i+1;j<N;j++) {
            dx=x[i]-x[j];
            dy=y[i]-y[j];
            dz=z[i]-z[j];
            if(ISNOTCUBIC!=3){
                if (PBCs==1 ) {
                    if (dx<-halfSidex) { dx+=sidex; }
                    else if (dx>halfSidex)   { dx-=sidex; }
                    if (dy<-halfSidey) { dy+=sidey; }
                    else if (dy>halfSidey)   { dy-=sidey; }
                    if (dz<-halfSidez) { dz+=sidez; }
                    else if (dz>halfSidez)   { dz-=sidez; }
                }
            }
            else {
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

            }
            sep2=dx*dx+dy*dy+dz*dz;
            if (rtype[i]!=rtype[j]) {   // AB interaction
                if (sep2 < rcut_AB2) {  // only calculate force if separation is less than rcut
                    invsep2=1.0/sep2;
                    invsep6=invsep2*invsep2*invsep2;
                    en=epsilon_AB*sigma_AB6*invsep6*(sigma_AB6*invsep6-1)-quarteruTail_AB;
                    en=4*en;
                    potential+=en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
            }
            else if (rtype[i]==1) {     // AA interaction
                if (sep2 < rcut2) { // only calculate force if separation is less than rcut
                    invsep2=1.0/sep2;
                    invsep6=invsep2*invsep2*invsep2;
                    en=invsep6*(invsep6-1)-quarteruTail_AA;
                    en=4*en;
                    potential+=en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
            }
            else {              // BB interaction
                if (sep2 < rcut_BB2) {  // only calculate force if separation is less than rcut
                    invsep2=1.0/sep2;
                    invsep6=invsep2*invsep2*invsep2;
                    en=epsilon_BB*sigma_BB6*invsep6*(sigma_BB6*invsep6-1)-quarteruTail_BB;
                    en=4*en;
                    potential+=en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
            }
        }
    }
}

void listBLJ() {    // binary Lennard-Jones potential using cell list
    double dx,dy,dz,sep2,invsep2,invsep6;   // sep for separations
    double en=0.0;
    int ic, i,j, jcell0, jcell,nabor;   // various counters

    for (ic=1;ic<=ncells_pot;ic++) {        // loop over all cells
        i=head_pot[ic];     // head of list particle for cell ic    
        while (i>0) {   // loop over all particles in ic
            j=llist_pot[i]; // next particle in current cell ic
            while (j>0) {   // loop over all particles in cell ic
                dx=x[i-1]-x[j-1];   // note use of i-1 to refer to particles counted from 0 in pos/velo/force arrays
                dy=y[i-1]-y[j-1];
                dz=z[i-1]-z[j-1];
                if (dx<-halfSidex) { dx+=sidex; }
                else if (dx>halfSidex)   { dx-=sidex; }
                if (dy<-halfSidey) { dy+=sidey; }
                else if (dy>halfSidey)   { dy-=sidey; }
                if (dz<-halfSidez) { dz+=sidez; }
                else if (dz>halfSidez)   { dz-=sidez; }
                sep2=dx*dx+dy*dy+dz*dz;
                if (rtype[i-1]!=rtype[j-1]) {   // AB interaction
                    if (sep2 < rcut_AB2) {  // only calculate force if separation is less than rcut
                        invsep2=1.0/sep2;
                        invsep6=invsep2*invsep2*invsep2;
                        en=epsilon_AB*sigma_AB6*invsep6*(sigma_AB6*invsep6-1)-quarteruTail_AB;
                        en=4*en;
                        potential+=en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                }
                else if (rtype[i-1]==1) {       // AA interaction
                    if (sep2 < rcut2) { // only calculate force if separation is less than rcut
                        invsep2=1.0/sep2;
                        invsep6=invsep2*invsep2*invsep2;
                        en=invsep6*(invsep6-1)-quarteruTail_AA;
                        en=4*en;
                        potential+=en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                }
                else {              // BB interaction
                    if (sep2 < rcut_BB2) {  // only calculate force if separation is less than rcut
                        invsep2=1.0/sep2;
                        invsep6=invsep2*invsep2*invsep2;
                        en=epsilon_BB*sigma_BB6*invsep6*(sigma_BB6*invsep6-1)-quarteruTail_BB;
                        en=4*en;
                        potential+=en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                }
                j=llist_pot[j]; // loop over next particle in cell ic
            }
            jcell0=13*(ic-1);       // now loop over adjacent cells to cell ic
            for (nabor=1;nabor<=13;nabor++) {
                jcell=map_pot[jcell0+nabor];    
                j=head_pot[jcell];  // head of cell for jcell
                while (j>0) {   // loop over head of cell and all other particles in jcell
                    dx=x[i-1]-x[j-1];   // note use of i-1 to refer to particles counted from 0 in pos/velo/force arrays
                    dy=y[i-1]-y[j-1];
                    dz=z[i-1]-z[j-1];
                    if (dx<-halfSidex) { dx+=sidex; }
                    else if (dx>halfSidex)   { dx-=sidex; }
                    if (dy<-halfSidey) { dy+=sidey; }
                    else if (dy>halfSidey)   { dy-=sidey; }
                    if (dz<-halfSidez) { dz+=sidez; }
                    else if (dz>halfSidez)   { dz-=sidez; }
                    sep2=dx*dx+dy*dy+dz*dz;
                    if (rtype[i-1]!=rtype[j-1]) {   // AB interaction
                        if (sep2 < rcut_AB2) {  // only calculate force if separation is less than rcut
                            invsep2=1.0/sep2;
                            invsep6=invsep2*invsep2*invsep2;
                            en=epsilon_AB*sigma_AB6*invsep6*(sigma_AB6*invsep6-1)-quarteruTail_AB;
                            en=4*en;
                            potential+=en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                    }
                    else if (rtype[i-1]==1) {       // AA interaction
                        if (sep2 < rcut2) { // only calculate force if separation is less than rcut
                            invsep2=1.0/sep2;
                            invsep6=invsep2*invsep2*invsep2;
                            en=invsep6*(invsep6-1)-quarteruTail_AA;
                            en=4*en;
                            potential+=en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                    }
                    else {              // BB interaction
                        if (sep2 < rcut_BB2) {  // only calculate force if separation is less than rcut
                            invsep2=1.0/sep2;
                            invsep6=invsep2*invsep2*invsep2;
                            en=epsilon_BB*sigma_BB6*invsep6*(sigma_BB6*invsep6-1)-quarteruTail_BB;
                            en=4*en;
                            potential+=en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                    }
                    j=llist_pot[j]; // next particle in jcell
                }
            }
            i=llist_pot[i]; // next particle in ic cell
        }
    }
}

void BLJSF() {  // binary Lennard-Jones potential with Stoddard-Ford truncation
    double dx,dy,dz,sep2,invsep2,invsep6;   // sep for separations
    double en;
    int i,j;
    
    for (i=0;i<N-1;i++) {
        for (j=i+1;j<N;j++) {
            dx=x[i]-x[j];
            dy=y[i]-y[j];
            dz=z[i]-z[j];
            if (PBCs==1) {
                if (dx<-halfSidex) { dx+=sidex; }
                else if (dx>halfSidex)   { dx-=sidex; }
                if (dy<-halfSidey) { dy+=sidey; }
                else if (dy>halfSidey)   { dy-=sidey; }
                if (dz<-halfSidez) { dz+=sidez; }
                else if (dz>halfSidez)   { dz-=sidez; }
            }
            if(ISNOTCUBIC==3){
                if (dz<-halfSidez) dz+=sidez;
                else if (dz>halfSidez) dz-=sidez;

                if (dy<-halfSidey){
                        dx-=tiltxy;
                        dy+=sidey;
                    }
                else if (dy>halfSidey) {
                    dx+=tiltxy;
                    dy-=sidey;
                    }

                if (dx<-halfSidex) dx+=sidex;
                else if (dx>halfSidex) dx-=sidex;
            }
            sep2=dx*dx+dy*dy+dz*dz;
            if (rtype[i]!=rtype[j]) {   // AB interaction
                if (sep2 < rcut_AB2) {  // only calculate force if separation is less than rcut
                    invsep2=1.0/sep2;
                    invsep6=invsep2*invsep2*invsep2;
                    en=epsilon_AB*sigma_AB6*invsep6*(sigma_AB6*invsep6-1)-quarteruTail_AB+(sep2-rcut_AB2)*stoddardford_AB;
                    en=4*en;
                    potential+=en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
            }
            else if (rtype[i]==1) {     // AA interaction
                if (sep2 < rcut2) { // only calculate force if separation is less than rcut
                    invsep2=1.0/sep2;
                    invsep6=invsep2*invsep2*invsep2;
                    en=invsep6*(invsep6-1)-quarteruTail_AA+(sep2-rcut2)*stoddardford_AA;
                    en=4*en;
                    potential+=en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
            }
            else {              // BB interaction
                if (sep2 < rcut_BB2) {  // only calculate force if separation is less than rcut
                    invsep2=1.0/sep2;
                    invsep6=invsep2*invsep2*invsep2;
                    en=epsilon_BB*sigma_BB6*invsep6*(sigma_BB6*invsep6-1)-quarteruTail_BB+(sep2-rcut_BB2)*stoddardford_BB;
                    en=4*en;
                    potential+=en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
            }
        }
    }
}

void listBLJSF() {  // binary Lennard-Jones potential with Stoddard-Ford truncation using cell list
    double dx,dy,dz,sep2,invsep2,invsep6;   // sep for separations
    double en=0.0;
    int ic, i,j, jcell0, jcell,nabor;   // various counters

    for (ic=1;ic<=ncells_pot;ic++) {        // loop over all cells
        i=head_pot[ic];     // head of list particle for cell ic    
        while (i>0) {   // loop over all particles in ic
            j=llist_pot[i]; // next particle in current cell ic
            while (j>0) {   // loop over all particles in cell ic
                dx=x[i-1]-x[j-1];   // note use of i-1 to refer to particles counted from 0 in pos/velo/force arrays
                dy=y[i-1]-y[j-1];
                dz=z[i-1]-z[j-1];
                if (dx<-halfSidex) { dx+=sidex; }
                else if (dx>halfSidex)   { dx-=sidex; }
                if (dy<-halfSidey) { dy+=sidey; }
                else if (dy>halfSidey)   { dy-=sidey; }
                if (dz<-halfSidez) { dz+=sidez; }
                else if (dz>halfSidez)   { dz-=sidez; }
                sep2=dx*dx+dy*dy+dz*dz;
                if (rtype[i-1]!=rtype[j-1]) {   // AB interaction
                    if (sep2 < rcut_AB2) {  // only calculate force if separation is less than rcut
                        invsep2=1.0/sep2;
                        invsep6=invsep2*invsep2*invsep2;
                        en=epsilon_AB*sigma_AB6*invsep6*(sigma_AB6*invsep6-1)-quarteruTail_AB+(sep2-rcut_AB2)*stoddardford_AB;
                        en=4*en;
                        potential+=en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                }
                else if (rtype[i-1]==1) {       // AA interaction
                    if (sep2 < rcut2) { // only calculate force if separation is less than rcut
                        invsep2=1.0/sep2;
                        invsep6=invsep2*invsep2*invsep2;
                        en=invsep6*(invsep6-1)-quarteruTail_AA+(sep2-rcut2)*stoddardford_AA;
                        en=4*en;
                        potential+=en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                }
                else {              // BB interaction
                    if (sep2 < rcut_BB2) {  // only calculate force if separation is less than rcut
                        invsep2=1.0/sep2;
                        invsep6=invsep2*invsep2*invsep2;
                        en=epsilon_BB*sigma_BB6*invsep6*(sigma_BB6*invsep6-1)-quarteruTail_BB+(sep2-rcut_BB2)*stoddardford_BB;
                        en=4*en;
                        potential+=en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                }
                j=llist_pot[j]; // loop over next particle in cell ic
            }
            jcell0=13*(ic-1);       // now loop over adjacent cells to cell ic
            for (nabor=1;nabor<=13;nabor++) {
                jcell=map_pot[jcell0+nabor];    
                j=head_pot[jcell];  // head of cell for jcell
                while (j>0) {   // loop over head of cell and all other particles in jcell
                    dx=x[i-1]-x[j-1];   // note use of i-1 to refer to particles counted from 0 in pos/velo/force arrays
                    dy=y[i-1]-y[j-1];
                    dz=z[i-1]-z[j-1];
                    if (dx<-halfSidex) { dx+=sidex; }
                    else if (dx>halfSidex)   { dx-=sidex; }
                    if (dy<-halfSidey) { dy+=sidey; }
                    else if (dy>halfSidey)   { dy-=sidey; }
                    if (dz<-halfSidez) { dz+=sidez; }
                    else if (dz>halfSidez)   { dz-=sidez; }
                    sep2=dx*dx+dy*dy+dz*dz;
                    if (rtype[i-1]!=rtype[j-1]) {   // AB interaction
                        if (sep2 < rcut_AB2) {  // only calculate force if separation is less than rcut
                            invsep2=1.0/sep2;
                            invsep6=invsep2*invsep2*invsep2;
                            en=epsilon_AB*sigma_AB6*invsep6*(sigma_AB6*invsep6-1)-quarteruTail_AB+(sep2-rcut_AB2)*stoddardford_AB;
                            en=4*en;
                            potential+=en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                    }
                    else if (rtype[i-1]==1) {       // AA interaction
                        if (sep2 < rcut2) { // only calculate force if separation is less than rcut
                            invsep2=1.0/sep2;
                            invsep6=invsep2*invsep2*invsep2;
                            en=invsep6*(invsep6-1)-quarteruTail_AA+(sep2-rcut2)*stoddardford_AA;
                            en=4*en;
                            potential+=en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                    }
                    else {              // BB interaction
                        if (sep2 < rcut_BB2) {  // only calculate force if separation is less than rcut
                            invsep2=1.0/sep2;
                            invsep6=invsep2*invsep2*invsep2;
                            en=epsilon_BB*sigma_BB6*invsep6*(sigma_BB6*invsep6-1)-quarteruTail_BB+(sep2-rcut_BB2)*stoddardford_BB;
                            en=4*en;
                            potential+=en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                    }
                    j=llist_pot[j]; // next particle in jcell
                }
            }
            i=llist_pot[i]; // next particle in ic cell
        }
    }
}

void MorYuk() { // Morse + Yukawa potential
    double dx,dy,dz,sep2,scaledsep,scaledsep2,invscaledsep,expon,yukexpon,en;   // sep for separations
    double invpsigmaij;
    int i,j,calcdsep;   // various counters
    
    en=scaledsep=scaledsep2=invscaledsep=0.0;   // reset energies
    calcdsep=0;
    
    for (i=0;i<N-1;i++) {       // loop over all cells
        for (j=i+1;j<N;j++) {       // loop over all cells
            dx=x[i]-x[j];   // note use of i-1 to refer to particles counted from 0 in pos/velo/force arrays
            dy=y[i]-y[j];
            dz=z[i]-z[j];
            if (PBCs==1) {
                if (dx<-halfSidex) { dx+=sidex; }
                else if (dx>halfSidex)   { dx-=sidex; }
                if (dy<-halfSidey) { dy+=sidey; }
                else if (dy>halfSidey)   { dy-=sidey; }
                if (dz<-halfSidez) { dz+=sidez; }
                else if (dz>halfSidez)   { dz-=sidez; }
            }
            sep2=dx*dx+dy*dy+dz*dz;
            invpsigmaij = 2.0 / (psigma[i] + psigma[j]);
            scaledsep2 = sep2*invpsigmaij*invpsigmaij;
            if (scaledsep2 < yukcut2 && yukcut2!=0.0) {
                scaledsep=pow(scaledsep2,0.5);
                invscaledsep=1.0/scaledsep;
                dx=dx*invpsigmaij;
                dy=dy*invpsigmaij;
                dz=dz*invpsigmaij;
                calcdsep=1;
                yukexpon=exp(-KAPPA*(scaledsep-1.0))*invscaledsep;
                en=yukepsilon*yukexpon - uYukTail;
                potential+=en;
                part_pot[i]+=en;
                part_pot[j]+=en;
            }
            if (scaledsep2 < mcut2) {   // only calculate force if separation is less than MCUT
                if (calcdsep==0) {
                    scaledsep=pow(scaledsep2,0.5);
                    invscaledsep=1.0/scaledsep;
                    dx=dx*invpsigmaij;
                    dy=dy*invpsigmaij;
                    dz=dz*invpsigmaij;
                }
                expon=exp(RHO0*(1.0-scaledsep));
                en=mepsilon*(expon*(expon-2.0)) - uMorseTail;
                potential+=en;
                part_pot[i]+=en;
                part_pot[j]+=en;
            }
            calcdsep=0;
        }
    }
}

void listMorYuk() { // Morse + Yukawa potential cell list
    double dx,dy,dz,sep2,scaledsep,scaledsep2,invscaledsep,expon,yukexpon,en;   // sep for separations
    double invpsigmaij;
    int ic, i,j, jcell0, jcell,nabor,calcdsep;  // various counters
    
    en=scaledsep=scaledsep2=invscaledsep=0.0;   // reset energies
    calcdsep=0;
    
    for (ic=1;ic<=ncells_pot;ic++) {        // loop over all cells
        i=head_pot[ic];     // head of list particle for cell ic    
        while (i>0) {   // loop over all particles in ic
            j=llist_pot[i]; // next particle in current cell ic
            while (j>0) {   // loop over all particles in cell ic
                dx=x[i-1]-x[j-1];   // note use of i-1 to refer to particles counted from 0 in pos/velo/force arrays
                dy=y[i-1]-y[j-1];
                dz=z[i-1]-z[j-1];
                if (dx<-halfSidex) { dx+=sidex; }
                else if (dx>halfSidex)   { dx-=sidex; }
                if (dy<-halfSidey) { dy+=sidey; }
                else if (dy>halfSidey)   { dy-=sidey; }
                if (dz<-halfSidez) { dz+=sidez; }
                else if (dz>halfSidez)   { dz-=sidez; }
                sep2=dx*dx+dy*dy+dz*dz;
                invpsigmaij = 2 / (psigma[i-1] + psigma[j-1]);
                scaledsep2 = sep2*invpsigmaij*invpsigmaij;
                if (scaledsep2 < yukcut2 && yukcut2!=0.0) {
                    scaledsep=pow(scaledsep2,0.5);
                    invscaledsep=1.0/scaledsep;
                    dx=dx*invpsigmaij;
                    dy=dy*invpsigmaij;
                    dz=dz*invpsigmaij;
                    calcdsep=1;
                    yukexpon=exp(-KAPPA*(scaledsep-1.0))*invscaledsep;
                    en=yukepsilon*yukexpon - uYukTail;
                    potential+=en;
                    part_pot[i-1]+=en;
                    part_pot[j-1]+=en;
                }
                if (scaledsep2 < mcut2) {   // only calculate force if separation is less than MCUT
                    if (calcdsep==0) {
                        scaledsep=pow(scaledsep2,0.5);
                        invscaledsep=1.0/scaledsep;
                        dx=dx*invpsigmaij;
                        dy=dy*invpsigmaij;
                        dz=dz*invpsigmaij;
                    }
                    expon=exp(RHO0*(1.0-scaledsep));
                    en=mepsilon*(expon*(expon-2.0)) - uMorseTail;
                    potential+=en;
                    part_pot[i-1]+=en;
                    part_pot[j-1]+=en;
                }
                calcdsep=0;
                j=llist_pot[j]; // loop over next particle in cell ic
            }
            jcell0=13*(ic-1);       // now loop over adjacent cells to cell ic
            for (nabor=1;nabor<=13;nabor++) {
                jcell=map_pot[jcell0+nabor];    
                j=head_pot[jcell];  // head of cell for jcell
                while (j>0) {   // loop over head of cell and all other particles in jcell
                    dx=x[i-1]-x[j-1];   // note use of i-1 to refer to particles counted from 0 in pos/velo/force arrays
                    dy=y[i-1]-y[j-1];
                    dz=z[i-1]-z[j-1];
                    if (dx<-halfSidex) { dx+=sidex; }
                    else if (dx>halfSidex)   { dx-=sidex; }
                    if (dy<-halfSidey) { dy+=sidey; }
                    else if (dy>halfSidey)   { dy-=sidey; }
                    if (dz<-halfSidez) { dz+=sidez; }
                    else if (dz>halfSidez)   { dz-=sidez; }
                    sep2=dx*dx+dy*dy+dz*dz;
                    invpsigmaij = 2 / (psigma[i-1] + psigma[j-1]);
                    scaledsep2 = sep2*invpsigmaij*invpsigmaij;
                    if (scaledsep2 < yukcut2 && yukcut2!=0.0) {
                        //printf("i%d j%d\n",i,j);
                        scaledsep=pow(scaledsep2,0.5);
                        invscaledsep=1.0/scaledsep;
                        dx=dx*invpsigmaij;
                        dy=dy*invpsigmaij;
                        dz=dz*invpsigmaij;
                        calcdsep=1;
                        yukexpon=exp(-KAPPA*(scaledsep-1.0))*invscaledsep;
                        en=yukepsilon*yukexpon - uYukTail;
                        potential+=en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                    if (scaledsep2 < mcut2) {   // only calculate force if separation is less than MCUT
                        if (calcdsep==0) {
                            scaledsep=pow(scaledsep2,0.5);
                            invscaledsep=1.0/scaledsep;
                            dx=dx*invpsigmaij;
                            dy=dy*invpsigmaij;
                            dz=dz*invpsigmaij;
                        }
                        expon=exp(RHO0*(1.0-scaledsep));
                        en=mepsilon*(expon*(expon-2.0)) - uMorseTail;
                        potential+=en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                    calcdsep=0;
                    j=llist_pot[j]; // next particle in jcell
                }
            }
            i=llist_pot[i]; // next particle in ic cell
        }
    }
    potential=en;   // the potential
}

void BIPL() {   // binary IPL potential
    double dx,dy,dz,sep2,invsep2,ff;    // sep for separations
    double en=0.0;
    int i,j;
    
    for (i=0;i<N-1;i++) {
        for (j=i+1;j<N;j++) {
            dx=x[i]-x[j];
            dy=y[i]-y[j];
            dz=z[i]-z[j];
            if (PBCs==1) {
                if (dx<-halfSidex) { dx+=sidex; }
                else if (dx>halfSidex)   { dx-=sidex; }
                if (dy<-halfSidey) { dy+=sidey; }
                else if (dy>halfSidey)   { dy-=sidey; }
                if (dz<-halfSidez) { dz+=sidez; }
                else if (dz>halfSidez)   { dz-=sidez; }
            }
            if(ISNOTCUBIC==3){
                    if (dz<-halfSidez) dz+=sidez;
                    else if (dz>halfSidez) dz-=sidez;

                if (dy<-halfSidey){
                        dx-=tiltxy;
                        dy+=sidey;
                    }
                else if (dy>halfSidey) {
                    dx+=tiltxy;
                    dy-=sidey;
                }   

                if (dx<-halfSidex) dx+=sidex;
                else if (dx>halfSidex) dx-=sidex;
            }
            sep2=dx*dx+dy*dy+dz*dz;
            if (rtype[i]!=rtype[j]) {   // AB interaction
                if (sep2 < rcut_AB2) {  // only calculate force if separation is less than rcut
                    invsep2=1.0/sep2;
                    ff=epsilon_AB*pow((sigma_AB2*invsep2),half_ipl_exp);
                    ff=invsep2*ff;
                    en=ff-uTail_AB;
                    en=ipl_pre*en;
                    potential+=en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
            }
            else if (rtype[i]==1) {     // AA interaction
                if (sep2 < rcut2) { // only calculate force if separation is less than rcut
                    invsep2=1.0/sep2;
                    ff=pow((invsep2),half_ipl_exp);
                    ff=invsep2*ff;
                    en=ff-uTail_AA;
                    en=ipl_pre*en;
                    potential+=en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
            }
            else {              // BB interaction
                if (sep2 < rcut_BB2) {  // only calculate force if separation is less than rcut
                    invsep2=1.0/sep2;
                    ff=epsilon_BB*pow((sigma_BB2*invsep2),half_ipl_exp);
                    ff=invsep2*ff;
                    en=ff-uTail_BB;
                    en=ipl_pre*en;
                    potential+=en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
            }
        }
    }
}

void listBIPL() {   // binary IPL potential with cell list
    double dx,dy,dz,sep2,invsep2, ff;   // sep for separations
    double en=0.0;
    int ic, i,j, jcell0, jcell,nabor;   // various counters
    
    for (ic=1;ic<=ncells_pot;ic++) {        // loop over all cells
        i=head_pot[ic];     // head of list particle for cell ic    
        while (i>0) {   // loop over all particles in ic
            j=llist_pot[i]; // next particle in current cell ic
            while (j>0) {   // loop over all particles in cell ic
                dx=x[i-1]-x[j-1];   // note use of i-1 to refer to particles counted from 0 in pos/velo/force arrays
                dy=y[i-1]-y[j-1];
                dz=z[i-1]-z[j-1];
                if (dx<-halfSidex) { dx+=sidex; }
                else if (dx>halfSidex)   { dx-=sidex; }
                if (dy<-halfSidey) { dy+=sidey; }
                else if (dy>halfSidey)   { dy-=sidey; }
                if (dz<-halfSidez) { dz+=sidez; }
                else if (dz>halfSidez)   { dz-=sidez; }
                sep2=dx*dx+dy*dy+dz*dz;
                if (rtype[i-1]!=rtype[j-1]) {   // AB interaction
                    if (sep2 < rcut_AB2) {  // only calculate force if separation is less than rcut
                        invsep2=1.0/sep2;
                        ff=epsilon_AB*pow((sigma_AB2*invsep2),half_ipl_exp);
                        ff=invsep2*ff;
                        en=ff-uTail_AB;
                        en=ipl_pre*en;
                        potential+=en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                }
                else if (rtype[i-1]==1) {       // AA interaction
                    if (sep2 < rcut2) { // only calculate force if separation is less than rcut
                        invsep2=1.0/sep2;
                        ff=pow((invsep2),half_ipl_exp);
                        ff=invsep2*ff;
                        en=ff-uTail_AA;
                        en=ipl_pre*en;
                        potential+=en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                }
                else {              // BB interaction
                    if (sep2 < rcut_BB2) {  // only calculate force if separation is less than rcut
                        invsep2=1.0/sep2;
                        ff=epsilon_BB*pow((sigma_BB2*invsep2),half_ipl_exp);
                        ff=invsep2*ff;
                        en=ff-uTail_BB;
                        en=ipl_pre*en;
                        potential+=en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                }
                j=llist_pot[j]; // loop over next particle in cell ic
            }
            jcell0=13*(ic-1);       // now loop over adjacent cells to cell ic
            for (nabor=1;nabor<=13;nabor++) {
                jcell=map_pot[jcell0+nabor];    
                j=head_pot[jcell];  // head of cell for jcell
                while (j>0) {   // loop over head of cell and all other particles in jcell
                    dx=x[i-1]-x[j-1];   // note use of i-1 to refer to particles counted from 0 in pos/velo/force arrays
                    dy=y[i-1]-y[j-1];
                    dz=z[i-1]-z[j-1];
                    if (dx<-halfSidex) { dx+=sidex; }
                    else if (dx>halfSidex)   { dx-=sidex; }
                    if (dy<-halfSidey) { dy+=sidey; }
                    else if (dy>halfSidey)   { dy-=sidey; }
                    if (dz<-halfSidez) { dz+=sidez; }
                    else if (dz>halfSidez)   { dz-=sidez; }
                    sep2=dx*dx+dy*dy+dz*dz;
                    if (rtype[i-1]!=rtype[j-1]) {   // AB interaction
                        if (sep2 < rcut_AB2) {  // only calculate force if separation is less than rcut
                            invsep2=1.0/sep2;
                            ff=epsilon_AB*pow((sigma_AB2*invsep2),half_ipl_exp);
                            ff=invsep2*ff;
                            en=ff-uTail_AB;
                            en=ipl_pre*en;
                            potential+=en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                    }
                    else if (rtype[i-1]==1) {       // AA interaction
                        if (sep2 < rcut2) { // only calculate force if separation is less than rcut
                            invsep2=1.0/sep2;
                            ff=pow((invsep2),half_ipl_exp);
                            ff=invsep2*ff;
                            en=ff-uTail_AA;
                            en=ipl_pre*en;
                            potential+=en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                    }
                    else {              // BB interaction
                        if (sep2 < rcut_BB2) {  // only calculate force if separation is less than rcut
                            invsep2=1.0/sep2;
                            ff=epsilon_BB*pow((sigma_BB2*invsep2),half_ipl_exp);
                            ff=invsep2*ff;
                            en=ff-uTail_BB;
                            en=ipl_pre*en;
                            potential+=en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                    }
                    j=llist_pot[j]; // next particle in jcell
                }
            }
            i=llist_pot[i]; // next particle in ic cell
        }
    }
}

void SFBIPL() { // Stoddard-Ford binary IPL potential
    double dx,dy,dz,sep2,invsep2,ff;    // sep for separations
    double en=0.0;
    int i,j;
    
    for (i=0;i<N-1;i++) {
        for (j=i+1;j<N;j++) {
            dx=x[i]-x[j];
            dy=y[i]-y[j];
            dz=z[i]-z[j];
            if (dx<-halfSidex) { dx+=sidex; }
            else if (dx>halfSidex)   { dx-=sidex; }
            if (dy<-halfSidey) { dy+=sidey; }
            else if (dy>halfSidey)   { dy-=sidey; }
            if (dz<-halfSidez) { dz+=sidez; }
            else if (dz>halfSidez)   { dz-=sidez; }
            sep2=dx*dx+dy*dy+dz*dz;
            if (rtype[i]!=rtype[j]) {   // AB interaction
                if (sep2 < rcut_AB2) {  // only calculate force if separation is less than rcut
                    invsep2=1.0/sep2;
                    ff=epsilon_AB*pow((sigma_AB2*invsep2),half_ipl_exp);
                    en=en+ff-uTail_AB+(sep2-rcut_AB2)*stoddardford_AB;
                    en=ipl_pre*en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
            }
            else if (rtype[i]==1) {     // AA interaction
                if (sep2 < rcut2) { // only calculate force if separation is less than rcut
                    invsep2=1.0/sep2;
                    ff=pow((invsep2),half_ipl_exp);
                    en=en+ff-uTail_AA+(sep2-rcut2)*stoddardford_AA;
                    en=ipl_pre*en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
            }
            else {              // BB interaction
                if (sep2 < rcut_BB2) {  // only calculate force if separation is less than rcut
                    invsep2=1.0/sep2;
                    ff=epsilon_BB*pow((sigma_BB2*invsep2),half_ipl_exp);
                    en=en+ff-uTail_BB+(sep2-rcut_BB2)*stoddardford_BB;
                    en=ipl_pre*en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
            }
        }
    }
}

void listSFBIPL() { // Stoddard-Ford binary IPL potential with cell list
    double dx,dy,dz,sep2,invsep2,ff;    // sep for separations
    double en=0.0;
    int ic, i,j, jcell0, jcell,nabor;   // various counters
    
    for (ic=1;ic<=ncells;ic++) {        // loop over all cells
        i=head[ic];     // head of list particle for cell ic    
        while (i>0) {   // loop over all particles in ic
            j=llist[i]; // next particle in current cell ic
            while (j>0) {   // loop over all particles in cell ic
                dx=x[i-1]-x[j-1];   // note use of i-1 to refer to particles counted from 0 in pos/velo/force arrays
                dy=y[i-1]-y[j-1];
                dz=z[i-1]-z[j-1];
                if (dx<-halfSidex) { dx+=sidex; }
                else if (dx>halfSidex)   { dx-=sidex; }
                if (dy<-halfSidey) { dy+=sidey; }
                else if (dy>halfSidey)   { dy-=sidey; }
                if (dz<-halfSidez) { dz+=sidez; }
                else if (dz>halfSidez)   { dz-=sidez; }
                sep2=dx*dx+dy*dy+dz*dz;
                if (rtype[i-1]!=rtype[j-1]) {   // AB interaction
                    if (sep2 < rcut_AB2) {  // only calculate force if separation is less than rcut
                        invsep2=1.0/sep2;
                        ff=epsilon_AB*pow((sigma_AB2*invsep2),half_ipl_exp);
                        en=en+ff-uTail_AB+(sep2-rcut_AB2)*stoddardford_AB;
                        en=ipl_pre*en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                }
                else if (rtype[i-1]==1) {       // AA interaction
                    if (sep2 < rcut2) { // only calculate force if separation is less than rcut
                        invsep2=1.0/sep2;
                        ff=pow((invsep2),half_ipl_exp);
                        en=en+ff-uTail_AA+(sep2-rcut2)*stoddardford_AA;
                        en=ipl_pre*en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                }
                else {              // BB interaction
                    if (sep2 < rcut_BB2) {  // only calculate force if separation is less than rcut
                        invsep2=1.0/sep2;
                        ff=epsilon_BB*pow((sigma_BB2*invsep2),half_ipl_exp);
                        en=en+ff-uTail_BB+(sep2-rcut_BB2)*stoddardford_BB;
                        en=ipl_pre*en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                }
                j=llist[j]; // loop over next particle in cell ic
            }
            jcell0=13*(ic-1);       // now loop over adjacent cells to cell ic
            for (nabor=1;nabor<=13;nabor++) {
                jcell=map[jcell0+nabor];    
                j=head[jcell];  // head of cell for jcell
                while (j>0) {   // loop over head of cell and all other particles in jcell
                    dx=x[i-1]-x[j-1];   // note use of i-1 to refer to particles counted from 0 in pos/velo/force arrays
                    dy=y[i-1]-y[j-1];
                    dz=z[i-1]-z[j-1];
                    if (dx<-halfSidex) { dx+=sidex; }
                    else if (dx>halfSidex)   { dx-=sidex; }
                    if (dy<-halfSidey) { dy+=sidey; }
                    else if (dy>halfSidey)   { dy-=sidey; }
                    if (dz<-halfSidez) { dz+=sidez; }
                    else if (dz>halfSidez)   { dz-=sidez; }
                    sep2=dx*dx+dy*dy+dz*dz;
                    if (rtype[i-1]!=rtype[j-1]) {   // AB interaction
                        if (sep2 < rcut_AB2) {  // only calculate force if separation is less than rcut
                            invsep2=1.0/sep2;
                            ff=epsilon_AB*pow((sigma_AB2*invsep2),half_ipl_exp);
                            en=en+ff-uTail_AB+(sep2-rcut_AB2)*stoddardford_AB;
                            en=ipl_pre*en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                    }
                    else if (rtype[i-1]==1) {       // AA interaction
                        if (sep2 < rcut2) { // only calculate force if separation is less than rcut
                            invsep2=1.0/sep2;
                            ff=pow((invsep2),half_ipl_exp);
                            en=en+ff-uTail_AA+(sep2-rcut2)*stoddardford_AA;
                            en=ipl_pre*en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                    }
                    else {              // BB interaction
                        if (sep2 < rcut_BB2) {  // only calculate force if separation is less than rcut
                            invsep2=1.0/sep2;
                            ff=epsilon_BB*pow((sigma_BB2*invsep2),half_ipl_exp);
                            en=en+ff-uTail_BB+(sep2-rcut_BB2)*stoddardford_BB;
                            en=ipl_pre*en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                    }
                    j=llist[j]; // next particle in jcell
                }
            }
            i=llist[i]; // next particle in ic cell
        }
    }
}

void BLJ_WCA_s() {  // cubic smoothed WCA potential from Coslovich 2011
    double dx,dy,dz,sep,sep2,sqno,invsep2,invsep6;  // sep for separations
    double en=0.0;
    int i,j;
    
    for (i=0;i<N-1;i++) {
        for (j=i+1;j<N;j++) {
            dx=x[i]-x[j];
            dy=y[i]-y[j];
            dz=z[i]-z[j];
            if (PBCs==1) {
                if (dx<-halfSidex) { dx+=sidex; }
                else if (dx>halfSidex)   { dx-=sidex; }
                if (dy<-halfSidey) { dy+=sidey; }
                else if (dy>halfSidey)   { dy-=sidey; }
                if (dz<-halfSidez) { dz+=sidez; }
                else if (dz>halfSidez)   { dz-=sidez; }
            }
            if(ISNOTCUBIC==3){
                if (dz<-halfSidez) dz+=sidez;
                else if (dz>halfSidez) dz-=sidez;

                if (dy<-halfSidey){
                        dx-=tiltxy;
                        dy+=sidey;
                    }
                else if (dy>halfSidey) {
                    dx+=tiltxy;
                    dy-=sidey;
                    }

                if (dx<-halfSidex) dx+=sidex;
                else if (dx>halfSidex) dx-=sidex;
            }
            sep2=dx*dx+dy*dy+dz*dz;
            if (rtype[i]!=rtype[j]) {   // AB interaction
                if (sep2 < cubic_a_AB2) {   // only calculate force if separation is less than rcut
                    invsep2=1.0/sep2;
                    invsep6=invsep2*invsep2*invsep2;
                    en=en+4.0*epsilon_AB*sigma_AB6*invsep6*(sigma_AB6*invsep6-1)+cubicA_AB;
                    potential+=en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
                else if (sep2 < rcut_AB2) { // only calculate force if separation is less than rcut
                    sep=pow(sep2,0.5);
                    sqno=(rcut_AB-sep)*(rcut_AB-sep);
                    en=cubicB_AB*sqno*(rcut_AB-sep);
                    potential+=en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
            }
            else if (rtype[i]==1) {     // AA interaction
                if (sep2 < cubic_a_AA2) {   // only calculate force if separation is less than rcut
                    invsep2=1.0/sep2;
                    invsep6=invsep2*invsep2*invsep2;
                    en=en+4.0*invsep6*(invsep6-1)+cubicA_AA;
                    potential+=en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
                else if (sep2 < rcut2) {    // only calculate force if separation is less than rcut
                    sep=pow(sep2,0.5);
                    sqno=(rcut-sep)*(rcut-sep);
                    en=cubicB_AA*sqno*(rcut-sep);
                    potential+=en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
            }
            else {              // BB interaction
                if (sep2 < cubic_a_BB2) {   // only calculate force if separation is less than rcut
                    invsep2=1.0/sep2;
                    invsep6=invsep2*invsep2*invsep2;
                    en=en+4.0*epsilon_BB*sigma_BB6*invsep6*(sigma_BB6*invsep6-1)+cubicA_BB;
                    potential+=en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
                else if (sep2 < rcut_BB2) { // only calculate force if separation is less than rcut
                    sep=pow(sep2,0.5);
                    sqno=(rcut_BB-sep)*(rcut_BB-sep);
                    en=cubicB_BB*sqno*(rcut_BB-sep);
                    potential+=en;
                    part_pot[i]+=en;
                    part_pot[j]+=en;
                }
            }
        }
    }
}

void listBLJ_WCA_s() {  // cubic smoothed WCA potential with cell list from Coslovich 2011  
    double dx,dy,dz,sep,sep2,sqno,invsep2,invsep6;  // sep for separations
    double en=0.0;
    int ic, i,j, jcell0, jcell,nabor;   // various counters
    
    for (ic=1;ic<=ncells_pot;ic++) {        // loop over all cells
        i=head_pot[ic];     // head of list particle for cell ic    
        while (i>0) {   // loop over all particles in ic
            j=llist_pot[i]; // next particle in current cell ic
            while (j>0) {   // loop over all particles in cell ic
                dx=x[i-1]-x[j-1];   // note use of i-1 to refer to particles counted from 0 in pos/velo/force arrays
                dy=y[i-1]-y[j-1];
                dz=z[i-1]-z[j-1];
                if (dx<-halfSidex) { dx+=sidex; }
                else if (dx>halfSidex)   { dx-=sidex; }
                if (dy<-halfSidey) { dy+=sidey; }
                else if (dy>halfSidey)   { dy-=sidey; }
                if (dz<-halfSidez) { dz+=sidez; }
                else if (dz>halfSidez)   { dz-=sidez; }
                sep2=dx*dx+dy*dy+dz*dz;
                if (rtype[i-1]!=rtype[j-1]) {   // AB interaction
                    if (sep2 < cubic_a_AB2) {   // only calculate force if separation is less than rcut
                        invsep2=1.0/sep2;
                        invsep6=invsep2*invsep2*invsep2;
                        en=en+4.0*epsilon_AB*sigma_AB6*invsep6*(sigma_AB6*invsep6-1)+cubicA_AB;
                        potential+=en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                    else if (sep2 < rcut_AB2) { // only calculate force if separation is less than rcut
                        sep=pow(sep2,0.5);
                        sqno=(rcut_AB-sep)*(rcut_AB-sep);
                        en=cubicB_AB*sqno*(rcut_AB-sep);
                        potential+=en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                }
                else if (rtype[i-1]==1) {       // AA interaction
                    if (sep2 < cubic_a_AA2) {   // only calculate force if separation is less than rcut
                        invsep2=1.0/sep2;
                        invsep6=invsep2*invsep2*invsep2;
                        en=en+4.0*invsep6*(invsep6-1)+cubicA_AA;
                        potential+=en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                    else if (sep2 < rcut2) {    // only calculate force if separation is less than rcut
                        sep=pow(sep2,0.5);
                        sqno=(rcut-sep)*(rcut-sep);
                        en=cubicB_AA*sqno*(rcut-sep);
                        potential+=en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                }
                else {              // BB interaction
                    if (sep2 < cubic_a_BB2) {   // only calculate force if separation is less than rcut
                        invsep2=1.0/sep2;
                        invsep6=invsep2*invsep2*invsep2;
                        en=en+4.0*epsilon_BB*sigma_BB6*invsep6*(sigma_BB6*invsep6-1)+cubicA_BB;
                        potential+=en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                    else if (sep2 < rcut_BB2) { // only calculate force if separation is less than rcut
                        sep=pow(sep2,0.5);
                        sqno=(rcut_BB-sep)*(rcut_BB-sep);
                        en=cubicB_BB*sqno*(rcut_BB-sep);
                        potential+=en;
                        part_pot[i-1]+=en;
                        part_pot[j-1]+=en;
                    }
                }
                j=llist_pot[j]; // loop over next particle in cell ic
            }
            jcell0=13*(ic-1);       // now loop over adjacent cells to cell ic
            for (nabor=1;nabor<=13;nabor++) {
                jcell=map_pot[jcell0+nabor];    
                j=head_pot[jcell];  // head of cell for jcell
                while (j>0) {   // loop over head of cell and all other particles in jcell
                    dx=x[i-1]-x[j-1];   // note use of i-1 to refer to particles counted from 0 in pos/velo/force arrays
                    dy=y[i-1]-y[j-1];
                    dz=z[i-1]-z[j-1];
                    if (dx<-halfSidex) { dx+=sidex; }
                    else if (dx>halfSidex)   { dx-=sidex; }
                    if (dy<-halfSidey) { dy+=sidey; }
                    else if (dy>halfSidey)   { dy-=sidey; }
                    if (dz<-halfSidez) { dz+=sidez; }
                    else if (dz>halfSidez)   { dz-=sidez; }
                    sep2=dx*dx+dy*dy+dz*dz;
                    if (rtype[i-1]!=rtype[j-1]) {   // AB interaction
                        if (sep2 < cubic_a_AB2) {   // only calculate force if separation is less than rcut
                            invsep2=1.0/sep2;
                            invsep6=invsep2*invsep2*invsep2;
                            en=en+4.0*epsilon_AB*sigma_AB6*invsep6*(sigma_AB6*invsep6-1)+cubicA_AB;
                            potential+=en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                        else if (sep2 < rcut_AB2) { // only calculate force if separation is less than rcut
                            sep=pow(sep2,0.5);
                            sqno=(rcut_AB-sep)*(rcut_AB-sep);
                            en=cubicB_AB*sqno*(rcut_AB-sep);
                            potential+=en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                    }
                    else if (rtype[i-1]==1) {       // AA interaction
                        if (sep2 < cubic_a_AA2) {   // only calculate force if separation is less than rcut
                            invsep2=1.0/sep2;
                            invsep6=invsep2*invsep2*invsep2;
                            en=en+4.0*invsep6*(invsep6-1)+cubicA_AA;
                            potential+=en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                        else if (sep2 < rcut2) {    // only calculate force if separation is less than rcut
                            sep=pow(sep2,0.5);
                            sqno=(rcut-sep)*(rcut-sep);
                            en=cubicB_AA*sqno*(rcut-sep);
                            potential+=en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                    }
                    else {              // BB interaction
                        if (sep2 < cubic_a_BB2) {   // only calculate force if separation is less than rcut
                            invsep2=1.0/sep2;
                            invsep6=invsep2*invsep2*invsep2;
                            en=en+4.0*epsilon_BB*sigma_BB6*invsep6*(sigma_BB6*invsep6-1)+cubicA_BB;
                            potential+=en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                        else if (sep2 < rcut_BB2) { // only calculate force if separation is less than rcut
                            sep=pow(sep2,0.5);
                            sqno=(rcut_BB-sep)*(rcut_BB-sep);
                            en=cubicB_BB*sqno*(rcut_BB-sep);
                            potential+=en;
                            part_pot[i-1]+=en;
                            part_pot[j-1]+=en;
                        }
                    }
                    j=llist_pot[j]; // next particle in jcell
                }
            }
            i=llist_pot[i]; // next particle in ic cell
        }
    }
}

void CRVT() { // CRVT potential scaled so minimum at r=1 and well-depth \varepsilon=1
    double dx,dy,dz,sep,sep2;   // sep for separations
    double en=0.0;
    int i,j;
    
    for (i=0;i<N-1;i++) {
        for (j=i+1;j<N;j++) {
            dx=x[i]-x[j];
            dy=y[i]-y[j];
            dz=z[i]-z[j];
            if (PBCs==1) {
                if (dx<-halfSidex) { dx+=sidex; }
                else if (dx>halfSidex)   { dx-=sidex; }
                if (dy<-halfSidey) { dy+=sidey; }
                else if (dy>halfSidey)   { dy-=sidey; }
                if (dz<-halfSidez) { dz+=sidez; }
                else if (dz>halfSidez)   { dz-=sidez; }
            }
            sep2=dx*dx+dy*dy+dz*dz;
            if (sep2 < rcut2) {
                sep=sqrt(sep2);
                en=(437.96*exp(-2.2322*4.3706*sep)-0.18382*exp(-0.214*(sep*4.3706-3.5344)*(sep*4.3706-3.5344)))/0.132896398;
                potential+=en;
                part_pot[i]+=en;
                part_pot[j]+=en;
            }
        }
    }
}