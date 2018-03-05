#ifndef RINGS_H
#define RINGS_H 

void get_basic_clusters();	// get SP3/4/5 rings including particle n0
void get_basic_sp4_rings(int, int, int);   // {n0,n1,n2} is not an SP3 ring, is it an SP4 or SP5 ring?
void get_basic_sp5_rings(int, int, int, int);   // {n0,n1,n2,n3} is not an SP4 ring, is it an SP5 ring?

void get_sp3_clusters(int, int, int);   // Take {n0,n1,n2}, check SP3 ring and if so detect SP3a/b/c cluster
void get_sp4_clusters(int, int, int, int);   // Take {n0,n1,n2,n3}, check SP4 ring and if so detect SP4a/b/c cluster
void get_sp5_clusters(int, int, int, int, int);   // Take {n0,n1,n2,n3,n4}, check SP5 ring and if so detect SP5a/b/c cluster

void add_mem_sp3b(int particle_ID);
void add_mem_sp3c(int particle_ID);
void add_mem_sp4b(int particle_ID);
void add_mem_sp4c(int particle_ID);
void add_mem_sp5b(int particle_ID);
void add_mem_sp5c(int particle_ID);


#endif
