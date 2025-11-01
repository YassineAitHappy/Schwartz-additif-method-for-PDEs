#include "utils.h"
#include "vector_operations.h"
#include "matrix_operations.h"
#include "functions.h"
#include <mpi.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <stdexcept>

int main(int argc, char **argv)
{   
    

    double alpha = 1.0 + D * dt * (2.0 / (dx * dx) + 2.0 / (dy * dy));
    double beta = -D * dt / (dx * dx);
    double gamma = -D * dt / (dy * dy);
   
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &Np);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    int jbeg, jend;
    charge(me, Ny, Np, &jbeg, &jend);

    MPI_Status Status;

    double debut = MPI_Wtime();

    int size0 = 0;
    std::vector<double> source0, U0;

    int sizeNp = 0;
    std::vector<double> sourceNp, UNp;

    int sizeme = 0;
    std::vector<double> sourceme, Ume;

    if (me == Np - 1)
    {
        sizeNp = jend - jbeg + 1;
        UNp.resize(sizeNp * Nx, 0.0);
        sourceNp.resize(sizeNp * Nx, 0.0);
    }
    else
    {
        size0 = jend - jbeg + 1 + recouvrement;
        U0.resize(size0 * Nx, 0.0);
        source0.resize(size0 * Nx, 0.0);

        sizeme = jend - jbeg + 1 + recouvrement;
        Ume.resize(sizeme * Nx, 0.0);
        sourceme.resize(sizeme * Nx, 0.0);
    }

    std::vector<double> stencil_envoi_haut(Nx, 0.0), stencil_envoi_bas(Nx, 0.0);
    std::vector<double> stencil_reception_haut(Nx, 0.0), stencil_reception_bas(Nx, 0.0);

    double t ;
    double erreurSchwarz;
    int iterationSchwarz=0;
    for (int k = 0; k < iterations; k++)
    {
        t = (k + 1) * dt;
        int iter = 0;
        double erreur_schwarz = 1.0;
        while (iter < iterationSchwarzMax && erreur_schwarz > tolSchwarz)
        {
            if (me == 0)
            {
                std::vector<double> vect_bas(Nx, 0.0);
                for (int i = 0; i < Nx; i++)
                {
                    vect_bas[i] = U0[i + size0 * Nx - Nx];
                }

                source0 = terme_source(dt, beta, gamma, t , stencil_reception_haut, stencil_reception_bas, me, U0, recouvrement, cas);
                U0 = BiCGStab(recouvrement, beta, alpha, gamma, source0, me);

                for (int i = 0; i < Nx; i++)
                {
                    stencil_envoi_haut[i] = (1.0 - gamma) * U0[(size0 - recouvrement) * Nx + i];
                }

                double erreurSchwarzMax0 = 0.0;
                for (int i = 0; i < Nx; i++)
                {
                    erreurSchwarz = std::fabs(vect_bas[i] - U0[i + size0 * Nx - Nx]);
                    erreurSchwarzMax0 = std::fmax(erreurSchwarzMax0, erreurSchwarz);
                }
                erreurSchwarz = erreurSchwarzMax0;

                MPI_Send(&stencil_envoi_haut[0], Nx, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
                MPI_Recv(&stencil_reception_haut[0], Nx, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &Status);
            }
            else if (me == Np - 1)
            {
                std::vector<double> hautnp(Nx, 0.0);
                for (int i = 0; i < Nx; i++)
                {
                    hautnp[i] = UNp[i];
                }

                sourceNp = terme_source(dt, beta, gamma, t , stencil_reception_haut, stencil_reception_bas, me, UNp, recouvrement, cas);
                UNp = BiCGStab(recouvrement, beta, alpha, gamma, sourceNp, me);

                for (int i = 0; i < Nx; i++)
                {
                    stencil_envoi_bas[i] = (1.0 - gamma) * UNp[i + (recouvrement - 1) * Nx];
                }

                double erreurSchwarzMaxNp = 0.0;
                for (int i = 0; i < Nx; i++)
                {
                    erreurSchwarz = std::fabs(hautnp[i] - UNp[i]);
                    erreurSchwarzMaxNp = std::fmax(erreurSchwarzMaxNp, erreurSchwarz);
                }
                erreurSchwarz = erreurSchwarzMaxNp;

                MPI_Send(&stencil_envoi_bas[0], Nx, MPI_DOUBLE, Np - 2, 0, MPI_COMM_WORLD);
                MPI_Recv(&stencil_reception_bas[0], Nx, MPI_DOUBLE, Np - 2, 0, MPI_COMM_WORLD, &Status);
            }
            else
            {
                double erreur1 = 0.0, erreur2 = 0.0;
                std::vector<double> basMe(Nx, 0.0);
                std::vector<double> hautMe(Nx, 0.0);

                for (int i = 0; i < Nx; i++)
                {
                    basMe[i] = Ume[i + sizeme * Nx - Nx];
                    hautMe[i] = Ume[i];
                }

                sourceme = terme_source(dt, beta, gamma, t , stencil_reception_haut, stencil_reception_bas, me, Ume, recouvrement, cas);
                Ume = BiCGStab(recouvrement, beta, alpha, gamma, sourceme, me);

                for (int i = 0; i < Nx; i++)
                {
                    stencil_envoi_haut[i] = (1.0 - gamma) * Ume[(sizeme - recouvrement) * Nx + i];
                    stencil_envoi_bas[i] = (1.0 - gamma) * Ume[i + (recouvrement - 1) * Nx];
                }

                double erreurSchwarzMaxMe = 0.0;
                for (int i = 0; i < Nx; i++)
                {
                    erreur1 = std::fabs(Ume[i + sizeme * Nx - Nx] - basMe[i]);
                    erreur2 = std::fabs(Ume[i] - hautMe[i]);
                    erreurSchwarz = std::fmax(erreur1, erreur2);
                    erreurSchwarzMaxMe = std::fmax(erreurSchwarzMaxMe, erreurSchwarz);
                }
                erreurSchwarz = erreurSchwarzMaxMe;

                MPI_Send(&stencil_envoi_bas[0], Nx, MPI_DOUBLE, me - 1, 0, MPI_COMM_WORLD);
                MPI_Send(&stencil_envoi_haut[0], Nx, MPI_DOUBLE, me + 1, 0, MPI_COMM_WORLD);

                MPI_Recv(&stencil_reception_haut[0], Nx, MPI_DOUBLE, me + 1, 0, MPI_COMM_WORLD, &Status);
                MPI_Recv(&stencil_reception_bas[0], Nx, MPI_DOUBLE, me - 1, 0, MPI_COMM_WORLD, &Status);
            }
            MPI_Allreduce(&erreurSchwarz, &erreur_schwarz, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            iter++;
        }
        if (me==0 and k==1){
            printf("Iterations schwarz me 0 = %d\n", iter);
            iterationSchwarz = iter;
        }
        // Écriture des résultats locaux
        std::string filename = "res_" + std::to_string(cas) + "_" + std::to_string(me) + "_" + std::to_string(k) + ".dat";
        std::ofstream outfile(filename);
        if (!outfile)
        {
            std::cerr << "Erreur lors de l'ouverture du fichier " << filename << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (me == 0)
        {
            for (int i = 1; i <= Nx; ++i)
            {
                for (int j = 1; j <= size0; j++)
                {
                    outfile << i * dx << " " << (jbeg + j - 1) * dy << " "
                            << U0[(j - 1) * Nx + i - 1] << " "
                            << solution_exacte(i * dx, (jbeg + j - 1) * dy, cas, t ) << std::endl;
                }
            }
        }
        else if (me == Np - 1)
        {
            for (int i = 1; i <= Nx; ++i)
            {
                for (int j = 1; j <= sizeNp; j++)
                {
                    outfile << i * dx << " " << (jbeg + j - 1) * dy << " "
                            << UNp[(j - 1) * Nx + i - 1] << " "
                            << solution_exacte(i * dx, (jbeg + j - 1) * dy, cas, t ) << std::endl;
                }
            }
        }
        else
        {
            for (int i = 1; i <= Nx; ++i)
            {
                for (int j = 1; j <= sizeme; j++)
                {
                    outfile << i * dx << " " << (jbeg + j - 1) * dy << " "
                            << Ume[(j - 1) * Nx + i - 1] << " "
                            << solution_exacte(i * dx, (jbeg + j - 1) * dy, cas, t ) << std::endl;
                }
            }
        }
        outfile.close();
    }
    
    // Calcul de l'erreur pour ce domaine
    MPI_Barrier(MPI_COMM_WORLD);

    double erreurLocal = 0.0;
    if (me == 0)
    {
        std::vector<double> U0_ex(size0 * Nx, 0.0);
        for (int i = 1; i <= Nx; ++i)
        {
            for (int j = 1; j <= size0; j++)
            {
                int k = (j - 1) * Nx + i - 1;
                U0_ex[k] = solution_exacte(i * dx, (jbeg + j - 1) * dy, cas, (iterations)*dt);
            }
        }
        erreurLocal = norme(U0, U0_ex)*sqrt(dx*dy);
    }

    if (me == Np - 1)
    {
        std::vector<double> UNp_ex(sizeNp * Nx, 0.0);
        for (int i = 1; i <= Nx; ++i)
        {
            for (int j = 1; j <= sizeNp; j++)
            {
                int k = (j - 1) * Nx + i - 1;
                UNp_ex[k] = solution_exacte(i * dx, (jbeg + j - 1) * dy, cas, (iterations)*dt);
            }
        }
        erreurLocal = norme(UNp, UNp_ex)*sqrt(dx*dy);
    }

    if (me >= 1 && me <= Np - 2)
    {
        std::vector<double> Ume_ex(sizeme * Nx, 0.0);
        for (int i = 1; i <= Nx; ++i)
        {
            for (int j = 1; j <= sizeme; j++)
            {
                int k = (j - 1) * Nx + i - 1;
                Ume_ex[k] = solution_exacte(i * dx, (jbeg + j - 1) * dy, cas, (iterations)*dt);
            }
        }
        erreurLocal = norme(Ume, Ume_ex)*sqrt(dx*dy);
    }

    // On récupère l'erreur maximale globale (par exemple) pour avoir une idée du pire cas
    double erreurGlobal;
    MPI_Reduce(&erreurLocal, &erreurGlobal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    double fin = MPI_Wtime();
    double tempsLocal= fin - debut;

    // On récupère le temps maximal sur tous les processus
    double tempsGlobal;
    MPI_Reduce(&tempsLocal, &tempsGlobal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (me == 0)
    {
        printf("Temps global pour Np=%d: %f s\n", Np, tempsGlobal);
        printf("Erreur globale (max) = %.15f\n", erreurGlobal);

        // On enregistre dans un fichier CSV pour le post-traitement
        // Le fichier contient : Np, recouvrement, Nx, Ny, temps, erreur
        std::ofstream csvfile("results.csv", std::ios::app);
        if(!csvfile.good())
        {
            std::cerr << "Impossible d'ouvrir le fichier results.csv pour ecriture.\n";
        }
        else
        {
            // Si le fichier est vide, on écrit l'entête
            static bool header_written = false;
            if (csvfile.tellp() == 0) {
                csvfile << "cas,Np,recouvrement,Nx,Ny,temps,erreur,iterationSchwarz\n";
            }
            csvfile << cas << "," << Np << "," << recouvrement << "," << Nx << "," << Ny << "," << tempsGlobal << "," << erreurGlobal << "," << iterationSchwarz << "\n";
        }
    }

    MPI_Finalize();
    return 0;
}
