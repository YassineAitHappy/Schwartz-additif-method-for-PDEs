#include "matrix_operations.h"
#include "utils.h"
#include <cmath>
#include <stdexcept>

// Fonction pour calculer le produit matriciel avec les conditions de Robin
std::vector<double> produit_mat_vect(int recouvrement, double alpha, double beta, double gamma, const std::vector<double>& U, int me)
{
    int ibeg, iend;
    charge(me, Ny, Np, &ibeg, &iend);

    std::vector<double> resultat;

    // Domaine bas
    if (me == 0)
    {
        int N0 = (iend - ibeg + 1) + recouvrement;
        int size0 = N0 * Nx;
        resultat.resize(size0, 0.0);

        for (int i = 0; i < size0; i++)
        {
            if (i >= 0 && i < Nx) // Première ligne
            {
                if (i == 0)
                {
                    resultat[i] = alpha * U[i] + beta * U[i + 1] + gamma * U[i + Nx];
                }
                else if (i == Nx - 1)
                {
                    resultat[i] = beta * U[i - 1] + alpha * U[i] + gamma * U[i + Nx];
                }
                else
                {
                    resultat[i] = alpha * U[i] + beta * (U[i + 1] + U[i - 1]) + gamma * U[i + Nx];
                }
            }
            else if (i >= Nx && i < size0 - Nx) // Corps du domaine
            {
                if (i % Nx == 0)
                {
                    resultat[i] = alpha * U[i] + beta * U[i + 1] + gamma * (U[i + Nx] + U[i - Nx]);
                }
                else if (i % Nx == Nx - 1)
                {
                    resultat[i] = beta * U[i - 1] + alpha * U[i] + gamma * (U[i + Nx] + U[i - Nx]);
                }
                else
                {
                    resultat[i] = alpha * U[i] + beta * (U[i + 1] + U[i - 1]) + gamma * (U[i + Nx] + U[i - Nx]);
                }
            }
            else // Bordure supérieure avec recouvrement
            {
                resultat[i] = (1.0 - gamma) * U[i];
            }
        }

        return resultat;
    }

    // Domaine haut
    else if (me == Np - 1)
    {
        int index = iend - ibeg + 1;
        int sizeNp = index * Nx;
        resultat.resize(sizeNp, 0.0);

        for (int i = 0; i < sizeNp; i++)
        {
            if (i >= 0 && i < Nx) // Première ligne
            {
                resultat[i] = (1.0 - gamma) * U[i];
            }
            else if (i >= Nx && i < sizeNp - Nx) // Corps du domaine
            {
                if (i % Nx == 0)
                {
                    resultat[i] = alpha * U[i] + beta * U[i + 1] + gamma * (U[i + Nx] + U[i - Nx]);
                }
                else if (i % Nx == Nx - 1)
                {
                    resultat[i] = beta * U[i - 1] + alpha * U[i] + gamma * (U[i + Nx] + U[i - Nx]);
                }
                else
                {
                    resultat[i] = alpha * U[i] + beta * (U[i + 1] + U[i - 1]) + gamma * (U[i + Nx] + U[i - Nx]);
                }
            }
            else // Dernière ligne avec recouvrement
            {
                if (i == sizeNp - Nx)
                {
                    resultat[i] = alpha * U[i] + beta * U[i + 1] + gamma * U[i - Nx];
                }
                else if (i == sizeNp - 1)
                {
                    resultat[i] = beta * U[i - 1] + alpha * U[i] + gamma * U[i - Nx];
                }
                else
                {
                    resultat[i] = alpha * U[i] + beta * (U[i + 1] + U[i - 1]) + gamma * U[i - Nx];
                }
            }
        }

        return resultat;
    }

    // Domaine milieu
    else
    {
        int index = (iend - ibeg + 1) + recouvrement;
        int sizeme = index * Nx;
        resultat.resize(sizeme, 0.0);

        for (int i = 0; i < sizeme; i++)
        {
            if (i >= 0 && i < Nx) // Première ligne
            {
                resultat[i] = (1.0 - gamma) * U[i];
            }
            else if (i >= Nx && i < sizeme - Nx) // Corps du domaine
            {
                if (i % Nx == 0)
                {
                    resultat[i] = alpha * U[i] + beta * U[i + 1] + gamma * (U[i + Nx] + U[i - Nx]);
                }
                else if (i % Nx == Nx - 1)
                {
                    resultat[i] = beta * U[i - 1] + alpha * U[i] + gamma * (U[i + Nx] + U[i - Nx]);
                }
                else
                {
                    resultat[i] = alpha * U[i] + beta * (U[i + 1] + U[i - 1]) + gamma * (U[i + Nx] + U[i - Nx]);
                }
            }
            else // Bordure supérieure avec recouvrement
            {
                resultat[i] = (1.0 - gamma) * U[i];
            }
        }

        return resultat;
    }
}
