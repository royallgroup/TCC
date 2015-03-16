#ifndef RINGS_H
#define RINGS_H 

void Rings_gSP4(int, int , int , int ) ;   // {n0,n1,n2} is not an SP3 ring, is it an SP4 or SP5 ring?

void Rings_gSP5(int, int , int , int , int ) ;   // {n0,n1,n2,n3} is not an SP4 ring, is it an SP5 ring?

void Rings_aSP3(int, int , int , int ) ;   // Take {n0,n1,n2}, check SP3 ring and if so detect SP3a/b/c cluster
void Rings_aSP4(int, int , int , int , int ) ;   // Take {n0,n1,n2,n3}, check SP4 ring and if so detect SP4a/b/c cluster
void Rings_aSP5(int, int , int , int , int , int ) ;   // Take {n0,n1,n2,n3,n4}, check SP5 ring and if so detect SP5a/b/c cluster

void Rings_setSP3c(int) ;// store cluster 5A D3h from Bonds_aSP3

void Rings_setSP4c(int) ;// store cluster 6A Oh from Bonds_aSP4()

void Rings_setSP5c(int) ;// store cluster 7A D5h from Bonds_aSP5()

#endif
