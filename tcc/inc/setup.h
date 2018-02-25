#ifndef SETUP_H
#define SETUP_H

void Setup_Output_Files();
void Initialise_Global_Variables();    // Initialize bond detection
void Reset_Frame_Variables();   // reset one frame of the variables (memory for which allocated in Bonds_Init())
void Free_All_Variables();    // Reset bond detection variables
int icell(int tix, int tiy, int tiz);
void Setup_Cell_List();

#endif