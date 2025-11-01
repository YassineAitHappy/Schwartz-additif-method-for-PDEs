#ifndef VECTOR_OPERATIONS_H
#define VECTOR_OPERATIONS_H

#include <vector>

// Déclarations des opérations sur les vecteurs
std::vector<double> soustraction(const std::vector<double>& u, const std::vector<double>& v);
std::vector<double> somme(const std::vector<double>& u, const std::vector<double>& v);
std::vector<double> multiplication_scalaire(double lambda, const std::vector<double>& u);
double produit_scalaire(const std::vector<double>& u, const std::vector<double>& v);
double norme(const std::vector<double>& u, const std::vector<double>& v);

#endif // VECTOR_OPERATIONS_H
