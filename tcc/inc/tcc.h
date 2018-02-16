#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int icell(int tix, int tiy, int tiz);

void Bonds_GetBonds(int f);	// Get bonds using simple lengths
int Bonds_BondCheck(int i, int j);	// Returns 1 if i & j are bonded; 0 otherwise

