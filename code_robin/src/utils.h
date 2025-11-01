#ifndef UTILS_H
#define UTILS_H

// Déclaration (extern) des variables globales
extern int Nx;
extern int Ny;
extern double Lx;
extern double Ly;
extern double D;
extern int Np;
extern double dx;
extern double dy;
extern int iterations;
extern double Tfinal;
extern double dt;
extern double alphaRobin;
extern double betaRobin;
extern double coeffRobin;
extern int cas;
extern double tolSchwarz;
extern int iterationSchwarzMax;
extern int recouvrement;

// Prototype de la fonction charge
void charge(int rang, int Ny, int Np, int *ibeg, int *iend);

#endif // UTILS_H