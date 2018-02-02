/* Alex Malins - alex.malins@gmail.com */
/* TCC: A topological cluster classification code with temporal tracking of clusters. */
/* Not for general consumption */	

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tools.h"
#include "setup.h"

int icell(int tix, int tiy, int tiz);

void Bonds_GetBonds(int f);	// Get bonds using simple lengths
int Bonds_BondCheck(int i, int j);	// Returns 1 if i & j are bonded; 0 otherwise

void Stats_Init();	// initialize Stats routine
void Stats_Reset();	// reset Cluster routine variables
void Stats_FreeMem();	// free memory from stats variables
void Stats_SetA();	// Set arrays to true if the ith particle is a member of any clusters with this or a larger number of particles
void Stats_Analyse();	// output Cluster statistics to file
void Stats_Report(char *filename);
void Pop_Per_Frame(int f);