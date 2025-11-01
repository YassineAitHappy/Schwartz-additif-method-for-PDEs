#!/bin/bash

# Définir un tableau avec le nombre de processus souhaités
np_values=(2 3 4)

# Boucle sur chaque valeur de processus
for np in "${np_values[@]}"; do
    echo "Lancement de la simulation avec $np processus..."
    mpiexec -np $np simulation
    echo "Simulation avec $np processus terminée."
done

echo "Toutes les simulations ont été exécutées."

