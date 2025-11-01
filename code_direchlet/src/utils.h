// utils.h
#ifndef UTILS_H
#define UTILS_H

#include <vector>

// Déclarations des variables globales
extern int Nx;
extern int Ny;
extern int Np; // Nombre de processus
extern int me; // Rang du processus actuel
extern double Lx;
extern double Ly;
extern double D;
extern double dx;
extern double dy;
extern int iterations;
extern double Tfinal;
extern double dt;
extern int cas;
extern double tolSchwarz;
extern int iterationSchwarzMax;
extern int recouvrement;

// Prototype de la fonction charge
void charge(int rang, int Ny, int Np, int *ibeg, int *iend);
#endif // UTILS_H
