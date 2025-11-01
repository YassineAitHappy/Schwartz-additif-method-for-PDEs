#include "matrix_operations.h"
#include "utils.h"
#include <vector>

std::vector<double> produit_mat_vect(int recouvrement, double alpha, double beta, 
                                double gamma, double coeff_Robin,
                                const std::vector<double>& U, int me)
{
    int ibeg, iend;
    charge(me, Ny, Np, &ibeg, &iend);

    // domaine bas
    if (me == 0)
    {
        int index = (iend - ibeg + 1) + recouvrement;
        int size0 = index * Nx;

        std::vector<double> resultat(size0, 0.0);

        for (int i = 0; i < size0; i++)
        {
            // Première ligne
            if (i < Nx)
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
            // Corps du domaine
            else if (i >= Nx && i < size0 - Nx)
            {
                if (i % Nx == 0)
                {
                    resultat[i] = alpha * U[i] + beta * U[i + 1] 
                              + gamma * (U[i + Nx] + U[i - Nx]);
                }
                else if (i % Nx == Nx - 1)
                {
                    resultat[i] = beta * U[i - 1] + alpha * U[i] 
                              + gamma * (U[i + Nx] + U[i - Nx]);
                }
                else
                {
                    resultat[i] = alpha * U[i] 
                              + beta * (U[i + 1] + U[i - 1]) 
                              + gamma * (U[i + Nx] + U[i - Nx]);
                }
            }
            // Bordure supérieure du domaine bas
            else
            {
                // La dernière ligne du domaine bas reçoit la condition de Robin
                if (i == size0 - Nx)
                {
                    resultat[i] = (alpha + coeff_Robin) * U[i] 
                              + beta * U[i + 1] 
                              + 2.0 * gamma * U[i - Nx];
                }
                else if (i == size0 - 1)
                {
                    resultat[i] = beta * U[i - 1] 
                              + (alpha + coeff_Robin) * U[i] 
                              + 2.0 * gamma * U[i - Nx];
                }
                else
                {
                    resultat[i] = (alpha + coeff_Robin) * U[i] 
                              + beta * (U[i + 1] + U[i - 1]) 
                              + 2.0 * gamma * U[i - Nx];
                }
            }
        }
        return resultat;
    }

    // domaine haut
    else if (me == Np - 1)
    {
        int index = (iend - ibeg + 1);
        int sizeNp = index * Nx;

        std::vector<double> resultat(sizeNp, 0.0);

        for (int i = 0; i < sizeNp; i++)
        {
            // Première ligne (bord inférieur du domaine haut) : Robin
            if (i < Nx)
            {
                if (i == 0)
                {
                    resultat[i] = (alpha + coeff_Robin) * U[i] 
                              + beta * U[i + 1] 
                              + 2.0 * gamma * U[i + Nx];
                }
                else if (i == Nx - 1)
                {
                    resultat[i] = beta * U[i - 1] 
                              + (alpha + coeff_Robin) * U[i] 
                              + 2.0 * gamma * U[i + Nx];
                }
                else
                {
                    resultat[i] = (alpha + coeff_Robin) * U[i] 
                              + beta * (U[i + 1] + U[i - 1]) 
                              + 2.0 * gamma * U[i + Nx];
                }
            }
            // Corps du domaine
            else if (i >= Nx && i < sizeNp - Nx)
            {
                if (i % Nx == 0)
                {
                    resultat[i] = alpha * U[i] + beta * U[i + 1] 
                              + gamma * (U[i + Nx] + U[i - Nx]);
                }
                else if (i % Nx == Nx - 1)
                {
                    resultat[i] = beta * U[i - 1] + alpha * U[i] 
                              + gamma * (U[i + Nx] + U[i - Nx]);
                }
                else
                {
                    resultat[i] = alpha * U[i] 
                              + beta * (U[i + 1] + U[i - 1]) 
                              + gamma * (U[i + Nx] + U[i - Nx]);
                }
            }
            // Dernière ligne (bord supérieur du domaine haut) : Dirichlet / Neumann « classique »
            else
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
                    resultat[i] = alpha * U[i] 
                              + beta * (U[i + 1] + U[i - 1]) 
                              + gamma * U[i - Nx];
                }
            }
        }
        return resultat;
    }

    // domaine milieu
    else
    {
        int index = (iend - ibeg + 1) + recouvrement;
        int sizeMe = index * Nx;

        std::vector<double> resultat(sizeMe, 0.0);

        for (int i = 0; i < sizeMe; i++)
        {
            // Première ligne du sous-domaine (bord inférieur) : Robin
            if (i < Nx)
            {
                if (i == 0)
                {
                    resultat[i] = (alpha + coeff_Robin) * U[i] 
                              + beta * U[i + 1] 
                              + 2.0 * gamma * U[i + Nx];
                }
                else if (i == Nx - 1)
                {
                    resultat[i] = beta * U[i - 1] 
                              + (alpha + coeff_Robin) * U[i] 
                              + 2.0 * gamma * U[i + Nx];
                }
                else
                {
                    resultat[i] = (alpha + coeff_Robin) * U[i] 
                              + beta * (U[i + 1] + U[i - 1]) 
                              + 2.0 * gamma * U[i + Nx];
                }
            }
            // Corps du domaine
            else if (i >= Nx && i < sizeMe - Nx)
            {
                if (i % Nx == 0)
                {
                    resultat[i] = alpha * U[i] + beta * U[i + 1] 
                              + gamma * (U[i + Nx] + U[i - Nx]);
                }
                else if (i % Nx == Nx - 1)
                {
                    resultat[i] = beta * U[i - 1] + alpha * U[i] 
                              + gamma * (U[i + Nx] + U[i - Nx]);
                }
                else
                {
                    resultat[i] = alpha * U[i] 
                              + beta * (U[i + 1] + U[i - 1]) 
                              + gamma * (U[i + Nx] + U[i - Nx]);
                }
            }
            // Dernière ligne (bord supérieur) : Robin
            else
            {
                if (i == sizeMe - Nx)
                {
                    resultat[i] = (alpha + coeff_Robin) * U[i] 
                              + beta * U[i + 1] 
                              + 2.0 * gamma * U[i - Nx];
                }
                else if (i == sizeMe - 1)
                {
                    resultat[i] = beta * U[i - 1] 
                              + (alpha + coeff_Robin) * U[i] 
                              + 2.0 * gamma * U[i - Nx];
                }
                else
                {
                    resultat[i] = (alpha + coeff_Robin) * U[i] 
                              + beta * (U[i + 1] + U[i - 1]) 
                              + 2.0 * gamma * U[i - Nx];
                }
            }
        }
        return resultat;
    }
}
