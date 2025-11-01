#include "functions.h"
#include "utils.h"
#include "vector_operations.h"
#include "matrix_operations.h"
#include <cmath>
#include <stdexcept>

double f(double x, double y, double t, int cas)
{
    double resultat = 0.0;
    if (cas == 1)
    {
        resultat = 2.0 * (y - std::pow(y, 2) + x - std::pow(x, 2));
    }
    else if (cas == 2)
    {
        resultat = std::sin(x) + std::cos(y);
    }
    else if (cas == 3)
    {
        resultat = std::exp(-std::pow(x - Lx / 2.0, 2))
              * std::exp(-std::pow(y - Ly / 2.0, 2))
              * std::cos(M_PI * t / 2.0);
    }
    else if (cas == 4)
    {
        resultat = 0.5 * M_PI * std::cos(0.5 * M_PI * t) * x * (1-x) * y * (1-y) + 2 * std::sin(0.5 * M_PI * t) * (x * (1-x) + y * (1-y));
    }
    return resultat;
}

double g(double x, double y, int cas)
{
    double resultat = 0.0;
    if (cas == 1)
    {
        resultat = 0.0;
    }
    else if (cas == 2)
    {
        resultat = std::sin(x) + std::cos(y);
    }
    else if (cas == 3 || cas == 4)
    {
        resultat = 0.0;
    }
    return resultat;
}

double h(double x, double y, int cas)
{
    double resultat = 0.0;
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


std::vector<double> terme_source(double dt, double beta, double gamma, double t,
                                 const std::vector<double>& stencil_recouv_haut,
                                 const std::vector<double>& stencil_recouv_bas,
                                 int me, const std::vector<double>& U,
                                 int recouvrement, int cas)
{
    int ibeg, iend;
    charge(me, Ny, Np, &ibeg, &iend);

    int size = U.size();
    std::vector<double> resultat(size, 0.0);

    if (me == 0)
    {
        for (int k = 0; k < size; k++)
        {
            int i = (k % Nx) + 1;
            int j = (k / Nx) + 1;

            resultat[k] = U[k] + dt * f(i * dx, (ibeg + j - 1) * dy, t, cas);

            if (j == 1)
            {
                resultat[k] += -gamma * g(i * dx, 0.0, cas);
                if (i == 1) {
                    resultat[k] += -beta * h(0.0, (ibeg + j - 1) * dy, cas);
                }
                else if (i == Nx) {
                    resultat[k] += -beta * h(Lx, (ibeg + j - 1) * dy, cas);
                }
            }
            else if (j >= 2 && j < (iend - ibeg + 1) + recouvrement)
            {
                // Condition sur x=0 ou x=Lx
                if (i == 1) {
                    resultat[k] += -beta * h(0.0, (ibeg + j - 1) * dy, cas);
                }
                else if (i == Nx) {
                    resultat[k] += -beta * h(Lx, (ibeg + j - 1) * dy, cas);
                }
            }
            else
            {
                resultat[k] += stencil_recouv_haut[i - 1];
                if (i == 1) {
                    resultat[k] += -beta * h(0.0, (ibeg + j - 1) * dy, cas);
                }
                else if (i == Nx) {
                    resultat[k] += -beta * h(Lx, (ibeg + j - 1) * dy, cas);
                }
            }
        }
        return resultat;
    }
    else if (me == Np - 1)
    {
        for (int k = 0; k < size; k++)
        {
            int i = (k % Nx) + 1;
            int j = (k / Nx) + 1;

            resultat[k] = U[k] + dt * f(i * dx, (ibeg + j - 1) * dy, t, cas);

            if (j == 1)
            {
                resultat[k] += stencil_recouv_bas[i - 1];
                if (i == 1) {
                    resultat[k] += -beta * h(0.0, (ibeg + j - 1) * dy, cas);
                }
                else if (i == Nx) {
                    resultat[k] += -beta * h(Lx, (ibeg + j - 1) * dy, cas);
                }
            }
            else if (j >= 2 && j < (iend - ibeg + 1))
            {
                if (i == 1) {
                    resultat[k] += -beta * h(0.0, (ibeg + j - 1) * dy, cas);
                }
                else if (i == Nx) {
                    resultat[k] += -beta * h(Lx, (ibeg + j - 1) * dy, cas);
                }
            }
            else
            {
                resultat[k] += -gamma * g(i * dx, Ly, cas);
                if (i == 1) {
                    resultat[k] += -beta * h(0.0, (ibeg + j - 1) * dy, cas);
                }
                else if (i == Nx) {
                    resultat[k] += -beta * h(Lx, (ibeg + j - 1) * dy, cas);
                }
            }
        }
        return resultat;
    }
    else
    {
        for (int k = 0; k < size; k++)
        {
            int i = (k % Nx) + 1;
            int j = (k / Nx) + 1;

            resultat[k] = U[k] + dt * f(i * dx, (ibeg + j - 1) * dy, t, cas);

            if (j == 1)
            {
                resultat[k] += stencil_recouv_bas[i - 1];
                if (i == 1) {
                    resultat[k] += -beta * h(0.0, (ibeg + j - 1) * dy, cas);
                }
                else if (i == Nx) {
                    resultat[k] += -beta * h(Lx, (ibeg + j - 1) * dy, cas);
                }
            }
            else if (j >= 2 && j < (iend - ibeg + 1) + recouvrement)
            {
                if (i == 1) {
                    resultat[k] += -beta * h(0.0, (ibeg + j - 1) * dy, cas);
                }
                else if (i == Nx) {
                    resultat[k] += -beta * h(Lx, (ibeg + j - 1) * dy, cas);
                }
            }
            else
            {
                resultat[k] += stencil_recouv_haut[i - 1];
                if (i == 1) {
                    resultat[k] += -beta * h(0.0, (ibeg + j - 1) * dy, cas);
                }
                else if (i == Nx) {
                    resultat[k] += -beta * h(Lx, (ibeg + j - 1) * dy, cas);
                }
            }
        }
        return resultat;
    }
}

std::vector<double> BiCGStab(int overlap, double betaCoeff, double alphaCoeff, double gammaCoeff, 
                             double robinCoefficient, const std::vector<double>& rightHandSide, int me)
{
    int vectorSize = rightHandSide.size();
    std::vector<double> residual(vectorSize), residualHat(vectorSize), solution(vectorSize), 
                        A_p(vectorSize), searchDirection(vectorSize), sVector(vectorSize), A_s(vectorSize);
    
    int currentIteration = 0;
    int maxIterations = 10000;
    double tolerance = 1e-8;
    
    double rhoCurrent = 1.0, rhoPrevious = 1.0;
    double alpha = 1.0, beta = 0.0, omega = 1.0;

    double residualNorm = 0.0, rightHandSideNorm = 0.0;

    // Initialisation des vecteurs
    for (int i = 0; i < vectorSize; i++)
    {
        solution[i] = 0.0;
        A_p[i] = 0.0;
        searchDirection[i] = 0.0;
    }

    residual = soustraction(rightHandSide, produit_mat_vect(overlap, alphaCoeff, betaCoeff, gammaCoeff, robinCoefficient, solution, me));
    residualHat = residual;

    residualNorm = std::sqrt(produit_scalaire(residual, residual));
    rightHandSideNorm = std::sqrt(produit_scalaire(rightHandSide, rightHandSide));

    while (residualNorm > tolerance * rightHandSideNorm && currentIteration < maxIterations)
    {
        rhoPrevious = rhoCurrent;
        rhoCurrent = produit_scalaire(residualHat, residual);

        beta = (rhoCurrent / rhoPrevious) * (alpha / omega);

        std::vector<double> omegaV = multiplication_scalaire(omega, A_p);
        std::vector<double> directionUpdate = soustraction(searchDirection, omegaV);
        std::vector<double> scaledDirectionUpdate = multiplication_scalaire(beta, directionUpdate);
        searchDirection = somme(residual, scaledDirectionUpdate);

        A_p = produit_mat_vect(overlap, alphaCoeff, betaCoeff, gammaCoeff, robinCoefficient, searchDirection, me);
        alpha = rhoCurrent / produit_scalaire(residualHat, A_p);

        std::vector<double> scaledA_p = multiplication_scalaire(alpha, A_p);
        sVector = soustraction(residual, scaledA_p);

        A_s = produit_mat_vect(overlap, alphaCoeff, betaCoeff, gammaCoeff, robinCoefficient, sVector, me);
        double denominatorA_s = produit_scalaire(A_s, A_s);
        if (denominatorA_s < 1e-30) { 
            //  éviter la division par zéro
            break;
        }
        omega = produit_scalaire(A_s, sVector) / denominatorA_s;

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
