#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include <vector>

// Déclaration de la fonction pour le produit matriciel (conditions de Robin)
std::vector<double> produit_mat_vect(int recouvrement, double alpha, double beta, 
                                double gamma, double coeff_Robin,
                                const std::vector<double>& U, int me);
                            

#endif // MATRIX_OPERATIONS_H
