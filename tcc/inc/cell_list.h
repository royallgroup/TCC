#ifndef TCC_CELL_LIST_H
#define TCC_CELL_LIST_H

void Get_Bonds_With_Voronoi_And_Cell_List();

void links();  // sorts all the particles into cells, result given by head-of-chain and linked list arrays

int icell(int tix, int tiy, int tiz);

void Setup_Cell_List();

#endif

