// utils.cpp
#include "utils.h"
#include <cmath>
#include <stdexcept>

// Définition des variables globales
int Nx = 50;
int Ny = 50;
int Np = 0; // Nombre de processus
int me = 0; // Rang du processus 
double Lx = 1.0;
double Ly = 1.0;
double D = 1.0;
double dx = Lx / (Nx + 1);
double dy = Ly / (Ny + 1);
int iterations = 2; //pour le cas stationnaire prendre iterations=2
double Tfinal = 1;
double dt = Tfinal / iterations;

// Choix du cas test :
int cas = 1;

// Données pour le bi-gradient conjugué
double tolSchwarz = 1e-6;
int iterationSchwarzMax = 10000;

// Longueur du recouvrement
int recouvrement = 5; 

// Fonction pour répartir les lignes entre les processus
void charge(int rang, int Ny, int Np, int *ibeg, int *iend)
{
    if (Np <= 0) {
        throw std::invalid_argument("Le nombre de processus doit être supérieur à 0.");
    }

    int i = Ny % Np;
    int j = Ny / Np;

    if (rang < i)
    {
        *ibeg = rang * (j + 1) + 1;
        *iend = (rang + 1) * (j + 1);
    }
    else
    {
        *ibeg = i + rang * j + 1;
        *iend = *ibeg + j - 1;
    }

    // Vérification des bornes
    if (*iend > Ny) {
        *iend = Ny;
    }
}