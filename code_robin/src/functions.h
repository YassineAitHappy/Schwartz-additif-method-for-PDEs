#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>

// Fonctions physiques / conditions aux bords
double f(double x, double y, double t, int cas);
double g(double x, double y, int cas);
double h(double x, double y, int cas);

// Solution exacte
double solution_exacte(double x, double y, int cas, double t);

// Construction du second membre (terme source + bords)
std::vector<double> terme_source(double dt, double beta, double gamma, double t,
                                 const std::vector<double>& stencil_recouv_haut,
                                 const std::vector<double>& stencil_recouv_bas,
                                 int me, const std::vector<double>& U,
                                 int recouvrement, int cas);

// Solveur BiCGStab
std::vector<double> BiCGStab(int overlap, double betaCoeff, double alphaCoeff, double gammaCoeff, 
                             double robinCoefficient, const std::vector<double>& rightHandSide, int me);

#endif // FUNCTIONS_H
