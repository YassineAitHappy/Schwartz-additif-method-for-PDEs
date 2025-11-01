#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector> // Ajouté pour utiliser std::vector

// Déclaration des fonctions mathématiques et des solutions exactes
double f(double x, double y, double t, int cas);
double g(double x, double y, int cas);
double h(double x, double y, int cas);
double solution_exacte(double x, double y, int cas, double t);
double condition_initiale(double x, double y);

// Déclaration des fonctions de calcul
std::vector<double> terme_source(double dt, double beta, double gamma, double t, 
                                  const std::vector<double>& stencil_recouv_haut, 
                                  const std::vector<double>& stencil_recouv_bas, 
                                  int me, const std::vector<double>& U, 
                                  int recouvrement, int cas);

std::vector<double> BiCGStab(int overlap, double betaCoeff, double alphaCoeff, double gammaCoeff, 
                             const std::vector<double>& rightHandSide, int me);

#endif // FUNCTIONS_H
