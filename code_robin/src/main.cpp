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
   
    // Paramètres 
    double alpha = 1.0 + D * dt * (2.0 / (dx * dx) + 2.0 / (dy * dy));
    double beta  = -D * dt / (dx * dx);
    double gamma = -D * dt / (dy * dy);
  
    // Initialisation MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &Np);

    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    int jbeg, jend;
    charge(me, Ny, Np, &jbeg, &jend);
    MPI_Status Status;

    double debut = MPI_Wtime();

    // Allocation des vecteurs
    std::vector<double> source0, U0;
    std::vector<double> sourcenp, Unp;
    std::vector<double> sourceme, Ume;

    int size0 = 0;
    int sizenp = 0;
    int sizeme = 0;

    if (me == Np - 1)  
    {
        sizenp = (jend - jbeg + 1);
        Unp.resize(sizenp * Nx, 0.0);
        sourcenp.resize(sizenp * Nx, 0.0);
    }
    else  
    {
        size0 = (jend - jbeg + 1) + recouvrement;
        U0.resize(size0 * Nx, 0.0);
        source0.resize(size0 * Nx, 0.0);

        sizeme = (jend - jbeg + 1) + recouvrement;
        Ume.resize(sizeme * Nx, 0.0);
        sourceme.resize(sizeme * Nx, 0.0);
    }

    std::vector<double> stencil_envoi_haut(Nx, 0.0),   stencil_envoi_bas(Nx, 0.0);
    std::vector<double> stencil_reception_haut(Nx, 0.0), stencil_reception_bas(Nx, 0.0);

    int iterationSchwarz = 0;
    double t = 0;
    // Boucle en temps
    for(int k = 0; k < iterations; k++)
    {
        t = (k + 1) * dt;

        int count = 0;
        double erreurSchwarz = 1.0;

        // Itération de Schwarz
        while (count < iterationSchwarzMax && erreurSchwarz > tolSchwarz)
        {
            double erreurSchwarzLocal = 0.0;

            if (me == 0)
            {
                std::vector<double> vect_bas(Nx);
                for (int i = 0; i < Nx; i++)
                {
                    vect_bas[i] = U0[i + size0 * Nx - Nx];
                }

                source0 = terme_source(dt, beta, gamma, t, stencil_reception_haut, stencil_reception_bas,
                                    me, U0, recouvrement, cas);

                U0 = BiCGStab(recouvrement, beta, alpha, gamma, coeffRobin, source0, me);

                for (int i = 0; i < Nx; i++)
                {
                    stencil_envoi_haut[i] = -gamma *
                        (-U0[(size0 - (recouvrement + 1)) * Nx + 2*Nx + i] 
                          + U0[(size0 - (recouvrement + 1)) * Nx + i])
                        + coeffRobin 
                        * U0[(size0 - (recouvrement + 1)) * Nx + Nx + i];
                }

                // Calcul erreur Schwarz
                double erreurSchwarzMax0 = 0.0;
                for (int i = 0; i < Nx; i++)
                {
                    double error = std::fabs(vect_bas[i] - U0[i + size0 * Nx - Nx]);
                    erreurSchwarzMax0 = std::fmax(erreurSchwarzMax0, error);
                }
                erreurSchwarzLocal = erreurSchwarzMax0;

                // Échange MPI
                MPI_Send(&stencil_envoi_haut[0], Nx, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
                MPI_Recv(&stencil_reception_haut[0], Nx, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &Status);
            }
            else if (me == Np - 1)
            {
                std::vector<double> hautnp(Nx);
                for (int i = 0; i < Nx; i++)
                {
                    hautnp[i] = Unp[i];
                }

                sourcenp = terme_source(dt, beta, gamma, t, stencil_reception_haut, stencil_reception_bas,
                                     me, Unp, recouvrement, cas);

                Unp = BiCGStab(recouvrement, beta, alpha, gamma, coeffRobin, sourcenp, me);

                for (int i = 0; i < Nx; i++)
                {
                    stencil_envoi_bas[i] = -gamma * 
                        (-Unp[i + (recouvrement - 2) * Nx] 
                          + Unp[i + (recouvrement - 2) * Nx + 2*Nx])
                        + coeffRobin 
                        * Unp[i + (recouvrement - 2) * Nx + Nx];
                }

                double erreurSchwarzMaxNp = 0.0;
                for (int i = 0; i < Nx; i++)
                {
                    double error = std::fabs(hautnp[i] - Unp[i]);
                    erreurSchwarzMaxNp = std::fmax(erreurSchwarzMaxNp, error);
                }
                erreurSchwarzLocal = erreurSchwarzMaxNp;

                // Échange MPI
                MPI_Send(&stencil_envoi_bas[0], Nx, MPI_DOUBLE, Np - 2, 0, MPI_COMM_WORLD);
                MPI_Recv(&stencil_reception_bas[0], Nx, MPI_DOUBLE, Np - 2, 0, MPI_COMM_WORLD, &Status);
            }
            else
            {
                std::vector<double> basMe(Nx), hautMe(Nx);
                for (int i = 0; i < Nx; i++)
                {
                    basMe[i] = Ume[i + sizeme * Nx - Nx];
                    hautMe[i]   = Ume[i];
                }

                sourceme = terme_source(dt, beta, gamma, t, stencil_reception_haut, stencil_reception_bas,
                                     me, Ume, recouvrement, cas);

                Ume = BiCGStab(recouvrement, beta, alpha, gamma, coeffRobin, sourceme, me);

                for (int i = 0; i < Nx; i++)
                {
                    stencil_envoi_haut[i] = -gamma *
                        (-Ume[(sizeme - (recouvrement + 1)) * Nx + 2*Nx + i]
                          + Ume[(sizeme - (recouvrement + 1)) * Nx + i])
                        + coeffRobin
                        * Ume[(sizeme - (recouvrement + 1)) * Nx + Nx + i];

                    stencil_envoi_bas[i] = -gamma *
                        (-Ume[i + (recouvrement - 2) * Nx]
                          + Ume[i + (recouvrement - 2) * Nx + 2*Nx])
                        + coeffRobin
                        * Ume[i + (recouvrement - 2) * Nx + Nx];
                }

                double erreurSchwarzMaxme = 0.0;
                for (int i = 0; i < Nx; i++)
                {
                    double erreur1 = std::fabs(Ume[i + sizeme * Nx - Nx] - basMe[i]);
                    double erreur2 = std::fabs(Ume[i] - hautMe[i]);
                    double error = std::fmax(erreur1, erreur2);
                    erreurSchwarzMaxme = std::fmax(erreurSchwarzMaxme, error);
                }
                erreurSchwarzLocal = erreurSchwarzMaxme;

                // Échanges MPI
                MPI_Send(&stencil_envoi_bas[0], Nx, MPI_DOUBLE, me - 1, 0, MPI_COMM_WORLD);
                MPI_Send(&stencil_envoi_haut[0], Nx, MPI_DOUBLE, me + 1, 0, MPI_COMM_WORLD);

                MPI_Recv(&stencil_reception_haut[0], Nx, MPI_DOUBLE, me + 1, 0, MPI_COMM_WORLD, &Status);
                MPI_Recv(&stencil_reception_bas[0], Nx, MPI_DOUBLE, me - 1, 0, MPI_COMM_WORLD, &Status);
            }

            MPI_Allreduce(&erreurSchwarzLocal, &erreurSchwarz, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            count++;
        }

        if (me == 0 && k == 0)
        {
            iterationSchwarz = count;
            std::cout << "Nombre d'iterations Schwarz (pas de temps k=0): " << iterationSchwarz << std::endl;
        }

        // Écriture des solutions locales
        if (me == 0)
        {
            std::string filename = "res_" + std::to_string(cas) + "_" 
                                   + std::to_string(me) + "_" + std::to_string(k) + ".dat";
            std::ofstream outfile(filename);
            for (int i = 1; i <= Nx; ++i)
            {
                for (int j = 1; j <= size0; j++)
                {
                    outfile << i * dx << " " << (jbeg + j - 1) * dy << " "
                            << U0[(j - 1) * Nx + (i - 1)] << " "
                            << solution_exacte(i * dx, (jbeg + j - 1) * dy, cas, t) << std::endl;
                }
            }
            outfile.close();
        }
        else if (me == Np - 1)
        {
            std::string filename = "res_" + std::to_string(cas) + "_" 
                                   + std::to_string(me) + "_" + std::to_string(k) + ".dat";
            std::ofstream outfile(filename);
            for (int i = 1; i <= Nx; ++i)
            {
                for (int j = 1; j <= (jend - jbeg + 1); j++)
                {
                    outfile << i * dx << " " << (jbeg + j - 1) * dy << " "
                            << Unp[(j - 1) * Nx + (i - 1)] << " "
                            << solution_exacte(i * dx, (jbeg + j - 1) * dy, cas, t) << std::endl;
                }
            }
            outfile.close();
        }
        else
        {
            std::string filename = "res_" + std::to_string(cas) + "_" 
                                   + std::to_string(me) + "_" + std::to_string(k) + ".dat";
            std::ofstream outfile(filename);
            for (int i = 1; i <= Nx; ++i)
            {
                for (int j = 1; j <= sizeme; j++)
                {
                    outfile << i * dx << " " << (jbeg + j - 1) * dy << " "
                            << Ume[(j - 1) * Nx + (i - 1)] << " "
                            << solution_exacte(i * dx, (jbeg + j - 1) * dy, cas, t) << std::endl;
                }
            }
            outfile.close();
        }
    }

    double erreurLocal = 0.0;
    if (me == 0)
    {
        std::vector<double> U0ex(size0 * Nx, 0.0);
        for (int i = 1; i <= Nx; ++i)
        {
            for (int j = 1; j <= size0; j++)
            {
                int k = (j - 1) * Nx + (i - 1);
                U0ex[k] = solution_exacte(i * dx, (jbeg + j - 1) * dy, cas, t);
            }
        }
        double err0 = norme(U0, U0ex)*sqrt(dx*dy);
        erreurLocal = err0;
        std::cout << "me " << me << ", L'erreur est " << err0 << std::endl;
    }

    if (me == Np - 1)
    {
        std::vector<double> Unpex((jend - jbeg + 1) * Nx, 0.0);
        for (int i = 1; i <= Nx; ++i)
        {
            for (int j = 1; j <= (jend - jbeg + 1); j++)
            {
                int k = (j - 1) * Nx + (i - 1);
                Unpex[k] = solution_exacte(i * dx, (jbeg + j - 1) * dy, cas, t);
            }
        }
        double erreurNp = norme(Unp, Unpex)*sqrt(dx*dy);
        erreurLocal = erreurNp;
        std::cout << "me " << me << ", L'erreur est " << erreurNp << std::endl;
    }

    if (me >= 1 && me <= Np - 2)
    {
        std::vector<double> UmeEx((jend - jbeg + 1 + recouvrement)*Nx, 0.0);
        for (int i = 1; i <= Nx; ++i)
        {
            for (int j = 1; j <= (jend - jbeg + 1 + recouvrement); j++)
            {
                int k = (j - 1) * Nx + (i - 1);
                UmeEx[k] = solution_exacte(i * dx, (jbeg + j - 1) * dy, cas, t);
            }
        }
        double erreurMe = norme(Ume, UmeEx)*sqrt(dx*dy);
        erreurLocal = erreurMe;
        std::cout << "me " << me << ", L'erreur est " << erreurMe << std::endl;
    }

    double fin = MPI_Wtime();
    double tempsLocal = fin - debut;
    std::cout << "me " << me << ", le temps final est " << tempsLocal << " s" << std::endl;
    double erreurGlobal = 0.0;
    MPI_Reduce(&erreurLocal, &erreurGlobal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // On récupère le temps maximal global 
    double tempsGlobal = 0.0;
    MPI_Reduce(&tempsLocal, &tempsGlobal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (me == 0)
    {
        std::cout << "Temps global (max) pour np=" << Np << " : " << tempsGlobal << " s" << std::endl;
        std::cout << "Erreur globale (max) = " << erreurGlobal << std::endl;

        // Append dans un fichier CSV
        std::ofstream csvfile("results.csv", std::ios::app);
        if (!csvfile.good())
        {
            std::cerr << "Impossible d'ouvrir le fichier results.csv en écriture.\n";
        }
        else
        {
            // S'il est vide, on écrit l'entête
            static bool header_written = false;
            if (csvfile.tellp() == 0)
            {
                csvfile << "cas,Np,recouvrement,Nx,Ny,temps,erreur,iterationSchwarz,alpha,beta\n";
            }
            csvfile << cas << ","
                    << Np << ","
                    << recouvrement << ","
                    << Nx << ","
                    << Ny << ","
                    << tempsGlobal << ","
                    << erreurGlobal << ","
                    << iterationSchwarz << ","
                    << alphaRobin << ","
                    << betaRobin << "\n";
            csvfile.close();
        }
    }

    MPI_Finalize();
    return 0;
}
