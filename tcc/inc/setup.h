#ifndef SETUP_H
#define SETUP_H

void Setup_Output_Files();
void Setup_InitStaticVars();    // Initialize bond detection
void Setup_ResetStaticVars();   // reset one frame of the variables (memory for which allocated in Bonds_Init())
void Setup_FreeStaticVars();    // Reset bond detection variables
int icell(int tix, int tiy, int tiz);
void Setup_Cell_List();

#endif