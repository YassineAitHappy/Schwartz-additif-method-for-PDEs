#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include <vector>

// Déclaration de la fonction pour le produit matriciel avec les conditions de Robin
std::vector<double> produit_mat_vect(int recouvrement, double alpha, double beta, double gamma, const std::vector<double>& U, int me)
;

#endif // MATRIX_OPERATIONS_H
