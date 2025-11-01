#include "utils.h"
#include <stdexcept>

// Définition des variables globales
int Nx = 50;
int Ny = 50;
double Lx = 1.0;
double Ly = 1.0;
double D  = 1.0;
int Np    = 0;
double dx = Lx / (Nx + 1);
double dy = Ly / (Ny + 1);
// Paramètres temporels
int iterations = 100; //pour le cas stationnaire prendre iterations=2
double Tfinal = 1.0;
double dt = Tfinal / iterations;
// Paramètres pour la condition de Robin
double alphaRobin = 1.0;
double betaRobin  = 0.0;
double coeffRobin = 2.0 * dt * betaRobin * D / (alphaRobin * dy);
// Choix du cas test
int cas = 4;
// Paramètres pour l'itération de Schwarz
double tolSchwarz = 1e-6;
int iterationSchwarzMax = 10000;
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