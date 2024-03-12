Statistics overall execution : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)   avg time         (min,max)                avg count (min,max)                  avg qty (min,max)              
                                       Temps total (01)       1.72  100.0%                                  0                                  0                        
                                   Preparer calcul (01)       1.72  100.0%                                  0                                  0                        
                          Resoudre (timestep loop) (01)       1.38   80.5%                                  1                                  0                        
                      SolveurSys::resoudre_systeme (01)       0.11    6.6%                                 20                                633                        
                    Solveurpetsc::resoudre_systeme (01)       0.00    0.0%                                  0                                  0                        
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00    0.0%                                  0                                  0                        
                  Operateur::calculer_pas_de_temps (01)       0.00    0.0%                                  0                                  0                        
                  ModeleTurbulence*::mettre_a_jour (01)       0.00    0.0%                                  0                                  0                        
                  Operateur_Conv::ajouter/calculer (01)       0.00    0.0%                                  0                                  0                        
                  Operateur_Diff::ajouter/calculer (01)       0.00    0.0%                                  0                                  0                        
                  Operateur_Decr::ajouter/calculer (01)       0.00    0.0%                                  0                                  0                        
                  Operateur_Grad::ajouter/calculer (01)       0.00    0.0%                                  0                                  0                        
                   Operateur_Div::ajouter/calculer (01)       0.00    0.0%                                  0                                  0                        
                          Source::ajouter/calculer (01)       0.00    0.0%                                  0                                  0                        
                               Pb_base::postraiter (01)       0.06    3.6%                                  2                                  0                        
                             Probleme_base::sauver (01)       0.00    0.0%                                  0                                  0                        
                                         temporary (01)       0.00    0.0%                                  0                                  0                        
                             Assembleur::assembler (01)       0.00    0.0%                                  0                                  0                        
                                            Divers (01)      -1.41  -82.1%                                  0                                  0                        
                                                m1 (01)       0.00    0.0%                                  0                                  0                        
                                                m2 (01)       0.00    0.0%                                  0                                  0                        
                                                m3 (01)       0.00    0.0%                                  0                                  0                        
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.00    0.1%                              2e+03                                  0                        
                        MPI_sendrecv:MPI_send_recv (02)       0.00    0.0%                                  0                                  0                        
                             MPI_sendrecv:MPI_send (02)       0.00    0.0%                                  0                                  0                        
                             MPI_sendrecv:MPI_recv (02)       0.00    0.0%                                  0                                  0                        
                        MPI_sendrecv:MPI_broadcast (02)       0.00    0.0%                                  0                                  0                        
                         MPI_sendrecv:MPI_alltoall (02)       0.00    0.0%                                  0                                  0                        
                      MPI_allreduce:MPI_partialsum (02)       0.00    0.0%                                  0                                  0                        
                       MPI_allreduce:MPI_sumdouble (02)       0.00    0.0%                                  0                                  0                        
                       MPI_allreduce:MPI_mindouble (02)       0.00    0.0%                                  0                                  0                        
                       MPI_allreduce:MPI_maxdouble (02)       0.00    0.0%                                  0                                  0                        
                          MPI_allreduce:MPI_sumint (02)       0.00    0.0%                                  0                                  0                        
                          MPI_allreduce:MPI_minint (02)       0.00    0.0%                                  0                                  0                        
                          MPI_allreduce:MPI_maxint (02)       0.00    0.0%                                  0                                  0                        
                         MPI_allreduce:MPI_barrier (02)       0.00    0.0%                                  0                                  0                        
                               io:MPI_send_recv_io (02)       0.00    0.0%                                  0                                  0                        
                                           Scatter (02)       0.00    0.0%                                  0                                  0                        
                                      ijk_exchange (01)       0.02    0.9%                                 80                                  0                        
                                      ijk_exchange (01)       0.03    1.5%                                498                                  0                        
                                  calcul dtstab QC (01)       0.00    0.1%                                  1                                  0                        
                               update statistiques (01)       0.04    2.3%                                  1                                  0                        
                         calcul terme acceleration (01)       0.00    0.0%                                  1                                  0                        
                                     checkpointing (01)       0.02    1.3%                                  1                                  0                        
                                         TF update (01)       0.00    0.0%                                  1                                  0                        
                                 qc convection rho (01)       0.02    1.3%                                  3                                  0                        
                     qc calcul pthermo,t,rho_v,etc (01)       0.01    0.8%                                  3                                  0                        
                              qc diffusion vitesse (01)       0.02    1.4%                                  3                                  0                        
                   qc diffusion turbulente vitesse (01)       0.00    0.0%                                  0                                  0                        
            qc diffusion vitesse modele structurel (01)       0.00    0.0%                                  0                                  0                        
                             qc convection vitesse (01)       0.02    1.4%                                  6                                  0                        
                                     qc rk3 update (01)       0.01    0.7%                                  3                                  0                        
                          qc diffusion temperature (01)       0.01    0.5%                                  3                                  0                        
                          qc turbulent diffusivity (01)       0.00    0.0%                                  0                                  0                        
               qc turbulent diffusivity structurel (01)       0.00    0.0%                                  0                                  0                        
                                    qc prepare rhs (01)       0.01    0.3%                                  3                                  0                        
                               qc resoudre systeme (01)       1.17   67.9%                                  3                                  0                        
                                 qc ajouter grad p (01)       0.01    0.3%                                  3                                  0                        
                                  jacobi_fine_mesh (01)       0.20   11.9%                                 52                                  0                        
                                        jacobi_all (01)       0.90   52.3%                                166                                  0                        
                                     coarse_solver (01)       0.11    6.7%                                 20                                  0                        
                                  jacobi_fine_mesh (01)       0.01    0.4%                                  3                                  0                        
                                        jacobi_all (01)       0.01    0.4%                                  3                                  0                        
                                      MPI_sendrecv (00)       0.00    0.0%                                  0                                  0                        
                                     MPI_allreduce (00)       0.00    0.0%                                  0                                  0                        
                                                io (00)       0.00    0.0%                                  0                                  0                        

