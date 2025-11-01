#include "vector_operations.h"
#include <stdexcept>
#include <math.h>

// Fonction pour soustraire deux vecteurs
std::vector<double> soustraction(const std::vector<double>& u, const std::vector<double>& v)
{
    if (u.size() != v.size())
    {
        throw std::invalid_argument("Les vecteurs doivent avoir la même taille pour la soustraction.");
    }

    int n = u.size();
    std::vector<double> s(n);
    for (int i = 0; i < n; i++)
    {
        s[i] = u[i] - v[i];
    }
    return s;
}

// Fonction pour additionner deux vecteurs
std::vector<double> somme(const std::vector<double>& u, const std::vector<double>& v)
{
    if (u.size() != v.size())
    {
        throw std::invalid_argument("Les vecteurs doivent avoir la même taille pour l'addition.");
    }

    int n = u.size();
    std::vector<double> resultat(n);
    for (int i = 0; i < n; i++)
    {
        resultat[i] = u[i] + v[i];
    }
    return resultat;
}

// Fonction pour multiplier un vecteur par un scalaire
std::vector<double> multiplication_scalaire(double lambda, const std::vector<double>& u)
{
    int n = u.size();
    std::vector<double> resultat(n);
    for (int i = 0; i < n; i++)
    {
        resultat[i] = lambda * u[i];
    }
    return resultat;
}

// Fonction pour calculer le produit scalaire de deux vecteurs
double produit_scalaire(const std::vector<double>& u, const std::vector<double>& v)
{
    if (u.size() != v.size())
    {
        throw std::invalid_argument("Les vecteurs doivent avoir la même taille pour le produit scalaire.");
    }

    double resultat = 0.0;
    for (size_t i = 0; i < u.size(); i++)
    {
        resultat += u[i] * v[i];
    }
    return resultat;
}

// Fonction pour calculer la norme L2 de la différence entre deux vecteurs
double norme(const std::vector<double>& u, const std::vector<double>& v)
{
    if (u.size() != v.size())
    {
        throw std::invalid_argument("Les vecteurs doivent avoir la même taille pour calculer la norme L2.");
    }
    int n = u.size();
    std::vector<double> E(n);
    E=soustraction(u,v);
    return std::sqrt(produit_scalaire(E,E));
}
