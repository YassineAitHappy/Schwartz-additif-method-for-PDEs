#include "functions.h"
#include "utils.h"
#include "vector_operations.h"
#include "matrix_operations.h"
#include <cmath>
#include <vector>
#include <stdexcept> 

// Fonction f
double f(double x, double y, double t, int cas)
{
    double resultat = 0.0; // Initialisé par défaut
    if (cas == 1)
    {
        resultat = 2 * (y - std::pow(y, 2) + x - std::pow(x, 2));
    }
    else if (cas == 2)
    {
        resultat = std::sin(x) + std::cos(y);
    }
    else if (cas == 3)
    {
        resultat = std::exp(-std::pow(x - Lx / 2.0, 2)) * std::exp(-std::pow(y - Ly / 2.0, 2)) * std::cos(M_PI * t / 2.0);
    }
   else if (cas == 4)
{
    resultat = 0.5 * M_PI * std::cos(0.5 * M_PI * t) * x * (1-x) * y * (1-y) + 2 * std::sin(0.5 * M_PI * t) * (x * (1-x) + y * (1-y));
    return resultat;
}

    else
    {
        throw std::invalid_argument("Cas non supporté dans la fonction f.");
    }
    return resultat;
}

// Fonction g
double g(double x, double y, int cas)
{
    double resultat = 0.0; // Initialisé par défaut
    if (cas == 1)
    {
        resultat = 0.0;
    }
    else if (cas == 2)
    {
        resultat = std::sin(x) + std::cos(y);
    }
    else if (cas == 3)
    {
        resultat = 0.0;
    }
    else if (cas == 4)
    {
        resultat = 0.0;
    }
    else
    {
        throw std::invalid_argument("Cas non supporté dans la fonction g.");
    }
    return resultat;
}

// Fonction h
double h(double x, double y, int cas)
{
    double resultat = 0.0; // Initialisé par défaut
    if (cas == 1)
    {
        resultat = 0.0;
    }
    else if (cas == 2)
    {
        resultat = std::sin(x) + std::cos(y);
    }
    else if (cas == 3)
    {
        resultat = 1.0;
    }
    else if (cas == 4)
    {
        resultat = 0.0;
    }
    else
    {
        throw std::invalid_argument("Cas non supporté dans la fonction h.");
    }
    return resultat;
}

// Fonction solution exacte
double solution_exacte(double x, double y, int cas, double t)
{
    if (cas == 1)
    {
        return x * (1.0 - x) * y * (1.0 - y);
    }
    else if (cas == 2)
    {
        return std::sin(x) + std::cos(y);
    }
    else if (cas == 3)
    {
        return 0;
    }
    else if (cas == 4)
    {
        return std::sin(M_PI*0.5*t)*x*(1-x)*(1-y)*y;
    }
    else
    {
        return 0.0;
    }
}

double condition_initiale(double x, double y)
{
    return x*y*(1-x)*(1-y);
}

// Fonction terme_source
std::vector<double> terme_source(double dt, double beta, double gamma, double t, 
                                  const std::vector<double>& stencil_recouv_haut, 
                                  const std::vector<double>& stencil_recouv_bas, 
                                  int me, const std::vector<double>& U, 
                                  int recouvrement, int cas)
{
    int jbeg, jend;
    charge(me, Ny, Np, &jbeg, &jend);

    int taille = U.size();
    std::vector<double> resultat(taille, 0.0); // Initialisé à 0

    if (me == 0)
    {
        for (int k = 0; k < taille; k++)
        {
            int i = k % Nx + 1;
            int j = k / Nx + 1;
            resultat[k] = U[k] + dt * f(i * dx, (jbeg + j - 1) * dy, t, cas);

            if (j == 1)
            {
                if (i >= 1 && i < Nx + 1)
                {
                    resultat[(j - 1) * Nx + (i - 1)] += -gamma * g(i * dx, 0.0, cas);

                    if (i == 1)
                    {
                        resultat[(j - 1) * Nx + (i - 1)] += -beta * h(0.0, (jbeg + j - 1) * dy, cas);
                    }
                    else if (i == Nx)
                    {
                        resultat[(j - 1) * Nx + (i - 1)] += -beta * h(Lx, (jbeg + j - 1) * dy, cas);
                    }
                }
            }
            else if (j >= 2 && j < (jend - jbeg + 1) + recouvrement)
            {
                if (i >= 1 && i < Nx + 1)
                {
                    if (i == 1)
                    {
                        resultat[(j - 1) * Nx + (i - 1)] += -beta * h(0.0, (jbeg + j - 1) * dy, cas);
                    }
                    else if (i == Nx)
                    {
                        resultat[(j - 1) * Nx + (i - 1)] += -beta * h(Lx, (jbeg + j - 1) * dy, cas);
                    }
                }
            }
            else
            {
                if (i >= 1 && i < Nx + 1)
                {
                    resultat[(j - 1) * Nx + (i - 1)] = stencil_recouv_haut[i - 1];
                }
            }
        }

        return resultat;
    }

    else if (me == Np - 1)
    {
        for (int k = 0; k < taille; k++)
        {
            int i = k % Nx + 1;
            int j = k / Nx + 1;
            resultat[k] = U[k] + dt * f(i * dx, (jbeg + j - 1) * dy, t, cas);

            if (j == 1)
            {
                if (i >= 1 && i < Nx + 1)
                {
                    resultat[(j - 1) * Nx + (i - 1)] = stencil_recouv_bas[i - 1];
                }
            }
            else if (j >= 2 && j < (jend - jbeg + 1))
            {
                if (i >= 1 && i < Nx + 1)
                {
                    if (i == 1)
                    {
                        resultat[(j - 1) * Nx + (i - 1)] += -beta * h(0.0, (jbeg + j - 1) * dy, cas);
                    }
                    else if (i == Nx)
                    {
                        resultat[(j - 1) * Nx + (i - 1)] += -beta * h(Lx, (jbeg + j - 1) * dy, cas);
                    }
                }
            }
            else
            {
                if (i >= 1 && i < Nx + 1)
                {
                    resultat[(j - 1) * Nx + (i - 1)] += -gamma * g(i * dx, Ly, cas);

                    if (i == 1)
                    {
                        resultat[(j - 1) * Nx + (i - 1)] += -beta * h(0.0, (jbeg + j - 1) * dy, cas);
                    }
                    else if (i == Nx)
                    {
                        resultat[(j - 1) * Nx + (i - 1)] += -beta * h(Lx, (jbeg + j - 1) * dy, cas);
                    }
                }
            }
        }

        return resultat;
    }

    else
    {
        for (int k = 0; k < taille; k++)
        {
            int i = k % Nx + 1;
            int j = k / Nx + 1;
            resultat[k] = U[k] + dt * f(i * dx, (jbeg + j - 1) * dy, t, cas);

            if (j == 1)
            {
                if (i >= 1 && i < Nx + 1)
                {
                    resultat[(j - 1) * Nx + (i - 1)] = stencil_recouv_bas[i - 1];
                }
            }
            else if (j >= 2 && j < (jend - jbeg + 1) + recouvrement)
            {
                if (i >= 1 && i < Nx + 1)
                {
                    if (i == 1)
                    {
                        resultat[(j - 1) * Nx + (i - 1)] += -beta * h(0.0, (jbeg + j - 1) * dy, cas);
                    }
                    else if (i == Nx)
                    {
                        resultat[(j - 1) * Nx + (i - 1)] += -beta * h(Lx, (jbeg + j - 1) * dy, cas);
                    }
                }
            }
            else
            {
                if (i >= 1 && i < Nx + 1)
                {
                    resultat[(j - 1) * Nx + (i - 1)] = stencil_recouv_haut[i - 1];
                }
            }
        }

        return resultat;
    }
}

// Fonction BiCGStab
std::vector<double> BiCGStab(int overlap, double betaCoeff, double alphaCoeff, double gammaCoeff, 
                             const std::vector<double>& rightHandSide, int me)
{
    int vectorSize = rightHandSide.size();
    std::vector<double> residual(vectorSize, 0.0), residualHat(vectorSize, 0.0), 
                        solution(vectorSize, 0.0), A_v(vectorSize, 0.0), 
                        searchDirection(vectorSize, 0.0), sVector(vectorSize, 0.0), 
                        A_s(vectorSize, 0.0);
    
    int currentIteration = 0;
    const int maxIterations = 10000;
    const double tolerance = 1e-8;
    
    double rhoCurrent = 1.0, rhoPrevious = 1.0;
    double alpha = 1.0, beta = 0.0, omega = 1.0;

    double residualNorm = 0.0, rightHandSideNorm = 0.0;

    // Initialisation des vecteurs
    std::vector<double> matrixProduct = produit_mat_vect(overlap, alphaCoeff, betaCoeff, gammaCoeff, solution, me);
    residual = soustraction(rightHandSide, matrixProduct);
    residualHat = residual;

    residualNorm = std::sqrt(produit_scalaire(residual, residual));
    rightHandSideNorm = std::sqrt(produit_scalaire(rightHandSide, rightHandSide));

    while (residualNorm > tolerance * rightHandSideNorm && currentIteration < maxIterations)
    {
        rhoPrevious = rhoCurrent;
        rhoCurrent = produit_scalaire(residualHat, residual);
        if (rhoCurrent == 0)
        {
            throw std::runtime_error("Breakdown: rho = 0");
        }
        beta = (rhoCurrent / rhoPrevious) * (alpha / omega);

        std::vector<double> omegaV = multiplication_scalaire(omega, A_v);
        std::vector<double> directionUpdate = soustraction(searchDirection, omegaV);
        std::vector<double> scaledDirectionUpdate = multiplication_scalaire(beta, directionUpdate);
        searchDirection = somme(residual, scaledDirectionUpdate);

        std::vector<double> matrixP = produit_mat_vect(overlap, alphaCoeff, betaCoeff, gammaCoeff, searchDirection, me);
        A_v = matrixP;
        double denominatorAlpha = produit_scalaire(residualHat, A_v);
        if (denominatorAlpha == 0)
        {
            throw std::runtime_error("Breakdown: denom_alpha = 0");
        }
        alpha = rhoCurrent / denominatorAlpha;

        std::vector<double> scaledA_v = multiplication_scalaire(alpha, A_v);
        sVector = soustraction(residual, scaledA_v);

        std::vector<double> matrixS = produit_mat_vect(overlap, alphaCoeff, betaCoeff, gammaCoeff, sVector, me);
        A_s = matrixS;
        double denominatorOmega = produit_scalaire(A_s, A_s);
        if (denominatorOmega == 0)
        {
            throw std::runtime_error("Breakdown: denom_omega = 0");
        }
        omega = produit_scalaire(A_s, sVector) / denominatorOmega;

        std::vector<double> alphaP = multiplication_scalaire(alpha, searchDirection);
        std::vector<double> omegaS = multiplication_scalaire(omega, sVector);
        std::vector<double> updateSolution = somme(alphaP, omegaS);
        solution = somme(solution, updateSolution);

        std::vector<double> omegaA_s = multiplication_scalaire(omega, A_s);
        residual = soustraction(sVector, omegaA_s);

        residualNorm = std::sqrt(produit_scalaire(residual, residual));
        rightHandSideNorm = std::sqrt(produit_scalaire(rightHandSide, rightHandSide));

        currentIteration++;
    }

    return solution;
}
