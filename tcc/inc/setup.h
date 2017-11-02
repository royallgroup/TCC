#ifndef SETUP_H
#define SETUP_H
#include <math.h>
#include "globals.h"
#include <string.h>

//// START: Setup routines
int Setup_GetFirstIntFromLine(FILE *);
double Setup_GetFirstDoubleFromLine(FILE *);
void Setup_ReadIniFile(char *);
void Setup_Readxyz(int e, int write, int f, FILE *);  // output single float valued arrays in gopenmol xyz format

void Setup_InitStaticVars();    // Initialize bond detection
void Setup_ResetStaticVars();   // reset one frame of the variables (memory for which allocated in Bonds_Init())
void Setup_FreeStaticVars();    // Reset bond detection variables
void Setup_InitPotentialVars(char *);   // Initialize potential energy variables
void Setup_print_U_r(); // print the potential file
int icell_pot(int tix, int tiy, int tiz);
void links_pot();
void BLJ();
void listBLJ();
void BLJSF();
void listBLJSF();
void MorYuk();
void listMorYuk();
void BIPL();
void listBIPL();
void SFBIPL();
void listSFBIPL();
void BLJ_WCA_s();
void listBLJ_WCA_s();
void CRVT();
void Setup_InitgsblVars(char *);    // Initialize ground state bond length deviation variables
void Update_bl_mom();
void Setup_InitDynamicVars(char *); // Initialize bond detection
void Setup_FreeDynamicVars();   // Reset bond detection variables
void Setup_FreePotentialVars();
void Setup_ReadBox(FILE *);  

void Setup_ClusComp() ; // zero arrays for cluster compostion analysis

void Setup_InitgsblVars(char *) ;
void Setup_InitgsblVars(char *) ; // Initialize ground state bond length deviation distribution arrays
void Setup_FreeDynamicVars() ;  // Free bond detection variables

#endif