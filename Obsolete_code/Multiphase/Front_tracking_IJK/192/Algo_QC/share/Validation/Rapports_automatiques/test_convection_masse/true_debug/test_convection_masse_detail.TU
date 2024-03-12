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

Statistics overall execution : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)   avg time         (min,max)                avg count (min,max)                  avg qty (min,max)              
                                       Temps total (01)       0.67  100.0% (0.67,0.67)                      0                                  0                        
                                   Preparer calcul (01)       0.67  100.0% (0.67,0.67)                      0                                  0                        
                          Resoudre (timestep loop) (01)       0.33   49.3% (0.33,0.34)                      1                                  0                        
                      SolveurSys::resoudre_systeme (01)       0.06    9.2% (0.06,0.08)                     20                           1.04e+03                        
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
                               Pb_base::postraiter (01)       0.02    3.6% (0.02,0.04)                      2                                  0                        
                             Probleme_base::sauver (01)       0.00    0.0%                                  0                                  0                        
                                         temporary (01)       0.00    0.0%                                  0                                  0                        
                             Assembleur::assembler (01)       0.00    0.0%                                  0                                  0                        
                                            Divers (01)      -0.35  -51.2% (-0.36,-0.33)                    0                                  0                        
                                                m1 (01)       0.00    0.0%                                  0                                  0                        
                                                m2 (01)       0.00    0.0%                                  0                                  0                        
                                                m3 (01)       0.00    0.0%                                  0                                  0                        
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.04    6.4% (0.04,0.05)               3.22e+03                                  0                        
                        MPI_sendrecv:MPI_send_recv (02)       0.10   14.3% (0.07,0.13)               3.56e+04                           6.48e+07 (6.48e+07,6.48e+07)    
                             MPI_sendrecv:MPI_send (02)       0.03    5.0% (0.00,0.05)                   62.3 (34,64)                      6e+05 (442,6.49e+05)         
                             MPI_sendrecv:MPI_recv (02)       0.00    0.2% (0.00,0.02)                   62.3 (2,1.09e+03)                 6e+05 (26,1.08e+07)          
                        MPI_sendrecv:MPI_broadcast (02)       0.02    2.4% (0.01,0.02)                    294                           3.76e+03                        
                         MPI_sendrecv:MPI_alltoall (02)       0.00    0.1% (0.00,0.00)                      4                                 16                        
                      MPI_allreduce:MPI_partialsum (02)       0.00    0.0%                                  0                                  0                        
                       MPI_allreduce:MPI_sumdouble (02)       0.04    6.5% (0.03,0.05)               2.42e+03                                  0                        
                       MPI_allreduce:MPI_mindouble (02)       0.00    0.0% (0.00,0.00)                      1                                  0                        
                       MPI_allreduce:MPI_maxdouble (02)       0.00    0.0%                                  0                                  0                        
                          MPI_allreduce:MPI_sumint (02)       0.00    0.0% (0.00,0.00)                     10                                  0                        
                          MPI_allreduce:MPI_minint (02)       0.00    0.0%                                  0                                  0                        
                          MPI_allreduce:MPI_maxint (02)       0.00    0.0%                                  0                                  0                        
                         MPI_allreduce:MPI_barrier (02)       0.06    8.7% (0.00,0.08)                     12                                  0                        
                               io:MPI_send_recv_io (02)       0.00    0.0%                                  0                                  0                        
                                           Scatter (02)       0.00    0.0%                                  0                                  0                        
                                      ijk_exchange (01)       0.03    4.7% (0.02,0.05)                     80                                  0                        
                                      ijk_exchange (01)       0.04    6.1% (0.03,0.06)                    498                                  0                        
                                  calcul dtstab QC (01)       0.00    0.0% (0.00,0.00)                      1                                  0                        
                               update statistiques (01)       0.01    0.9% (0.01,0.01)                      1                                  0                        
                         calcul terme acceleration (01)       0.00    0.0% (0.00,0.00)                      1                                  0                        
                                     checkpointing (01)       0.04    6.1% (0.04,0.04)                      1                                  0                        
                                         TF update (01)       0.00    0.0% (0.00,0.00)                      1                                  0                        
                                 qc convection rho (01)       0.00    0.3% (0.00,0.00)                      3                                  0                        
                     qc calcul pthermo,t,rho_v,etc (01)       0.00    0.4% (0.00,0.00)                      3                                  0                        
                              qc diffusion vitesse (01)       0.00    0.5% (0.00,0.00)                      3                                  0                        
                   qc diffusion turbulente vitesse (01)       0.00    0.0%                                  0                                  0                        
            qc diffusion vitesse modele structurel (01)       0.00    0.0%                                  0                                  0                        
                             qc convection vitesse (01)       0.00    0.3% (0.00,0.00)                      6                                  0                        
                                     qc rk3 update (01)       0.00    0.2% (0.00,0.00)                      3                                  0                        
                          qc diffusion temperature (01)       0.00    0.1% (0.00,0.00)                      3                                  0                        
                          qc turbulent diffusivity (01)       0.00    0.0%                                  0                                  0                        
               qc turbulent diffusivity structurel (01)       0.00    0.0%                                  0                                  0                        
                                    qc prepare rhs (01)       0.00    0.1% (0.00,0.00)                      3                                  0                        
                               qc resoudre systeme (01)       0.24   35.3% (0.24,0.24)                      3                                  0                        
                                 qc ajouter grad p (01)       0.00    0.2% (0.00,0.00)                      3                                  0                        
                                  jacobi_fine_mesh (01)       0.02    2.9% (0.02,0.02)                     52                                  0                        
                                        jacobi_all (01)       0.11   17.0% (0.11,0.12)                    166                                  0                        
                                     coarse_solver (01)       0.08   11.9% (0.08,0.08)                     20                                  0                        
                                  jacobi_fine_mesh (01)       0.00    0.3% (0.00,0.00)                      3                                  0                        
                                        jacobi_all (01)       0.00    0.4% (0.00,0.00)                      3                                  0                        
                                      MPI_sendrecv (00)       0.15   21.9% (0.11,0.19)               3.61e+04 (3.6e+04,3.71e+04)         6.6e+07 (6.54e+07,7.56e+07)    
                                     MPI_allreduce (00)       0.10   15.3% (0.04,0.13)               2.44e+03                                  0                        
                                                io (00)       0.00    0.0%                                  0                                  0                        

Statistics execution [PE 0] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.34  50.1%          1          0
                      SolveurSys::resoudre_systeme (01)       0.08  11.9%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.04   5.3%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.36 -53.5%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.04   6.0%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.08  12.2%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.00   0.0%         34        442
                             MPI_sendrecv:MPI_recv (02)       0.02   3.2%   1.09e+03   1.08e+07
                        MPI_sendrecv:MPI_broadcast (02)       0.01   1.2%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.04   6.4%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.00   0.0%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.02   3.7%         80          0
                                      ijk_exchange (01)       0.04   5.3%        498          0
                                  calcul dtstab QC (01)       0.00   0.1%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.3%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.5%          3          0
                              qc diffusion vitesse (01)       0.00   0.5%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.3%          6          0
                                     qc rk3 update (01)       0.00   0.2%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.3%          3          0
                                 qc ajouter grad p (01)       0.00   0.3%          3          0
                                  jacobi_fine_mesh (01)       0.02   3.0%         52          0
                                        jacobi_all (01)       0.11  16.8%        166          0
                                     coarse_solver (01)       0.08  12.0%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.4%          3          0
                                        jacobi_all (01)       0.00   0.4%          3          0
                                      MPI_sendrecv (00)       0.11  16.7%   3.71e+04   7.56e+07
                                     MPI_allreduce (00)       0.04   6.5%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 1] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.33  49.2%          1          0
                      SolveurSys::resoudre_systeme (01)       0.06   9.1%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.02   3.4%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.35 -51.6%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.04   6.5%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.11  16.4%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.02   3.1%         64   6.27e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.02   2.3%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.1%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.05   7.6%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.08  11.3%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.03   5.2%         80          0
                                      ijk_exchange (01)       0.05   7.5%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.3%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.3%          3          0
                              qc diffusion vitesse (01)       0.00   0.5%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.2%          6          0
                                     qc rk3 update (01)       0.00   0.1%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.2%          3          0
                                 qc ajouter grad p (01)       0.00   0.1%          3          0
                                  jacobi_fine_mesh (01)       0.02   2.5%         52          0
                                        jacobi_all (01)       0.11  16.9%        166          0
                                     coarse_solver (01)       0.08  11.9%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.1%          3          0
                                        jacobi_all (01)       0.00   0.2%          3          0
                                      MPI_sendrecv (00)       0.15  21.8%    3.6e+04   6.54e+07
                                     MPI_allreduce (00)       0.13  19.0%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 2] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.33  49.2%          1          0
                      SolveurSys::resoudre_systeme (01)       0.06   9.1%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.02   3.4%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.34 -50.5%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.04   6.1%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.09  13.1%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.02   3.5%         64   6.37e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.02   2.5%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.1%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.04   6.5%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.07  11.0%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.03   4.4%         80          0
                                      ijk_exchange (01)       0.04   5.5%        498          0
                                  calcul dtstab QC (01)       0.00   0.1%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.3%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.5%          3          0
                              qc diffusion vitesse (01)       0.00   0.5%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.3%          6          0
                                     qc rk3 update (01)       0.00   0.2%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.3%          3          0
                                 qc ajouter grad p (01)       0.00   0.3%          3          0
                                  jacobi_fine_mesh (01)       0.02   3.0%         52          0
                                        jacobi_all (01)       0.11  16.9%        166          0
                                     coarse_solver (01)       0.08  11.9%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.4%          3          0
                                        jacobi_all (01)       0.00   0.4%          3          0
                                      MPI_sendrecv (00)       0.13  19.2%    3.6e+04   6.54e+07
                                     MPI_allreduce (00)       0.12  17.5%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 3] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.33  49.2%          1          0
                      SolveurSys::resoudre_systeme (01)       0.06   9.1%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.02   3.5%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.34 -50.4%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.04   6.1%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.09  12.8%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.02   3.6%         64   6.27e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.02   2.3%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.04   6.5%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.07  10.8%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.03   4.3%         80          0
                                      ijk_exchange (01)       0.04   5.4%        498          0
                                  calcul dtstab QC (01)       0.00   0.1%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.3%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.5%          3          0
                              qc diffusion vitesse (01)       0.00   0.5%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.3%          6          0
                                     qc rk3 update (01)       0.00   0.2%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.3%          3          0
                                 qc ajouter grad p (01)       0.00   0.3%          3          0
                                  jacobi_fine_mesh (01)       0.02   3.0%         52          0
                                        jacobi_all (01)       0.11  16.9%        166          0
                                     coarse_solver (01)       0.08  11.9%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.4%          3          0
                                        jacobi_all (01)       0.00   0.4%          3          0
                                      MPI_sendrecv (00)       0.13  18.8%    3.6e+04   6.54e+07
                                     MPI_allreduce (00)       0.12  17.3%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 4] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.33  49.2%          1          0
                      SolveurSys::resoudre_systeme (01)       0.06   9.1%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.02   3.5%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.35 -52.1%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.04   6.7%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.12  17.2%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.03   4.1%         64   6.27e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.02   2.4%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.1%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.05   7.6%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.07  10.4%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.04   5.4%         80          0
                                      ijk_exchange (01)       0.05   7.7%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.2%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.3%          3          0
                              qc diffusion vitesse (01)       0.00   0.5%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.2%          6          0
                                     qc rk3 update (01)       0.00   0.1%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.2%          3          0
                                 qc ajouter grad p (01)       0.00   0.1%          3          0
                                  jacobi_fine_mesh (01)       0.02   2.5%         52          0
                                        jacobi_all (01)       0.11  16.9%        166          0
                                     coarse_solver (01)       0.08  11.9%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.1%          3          0
                                        jacobi_all (01)       0.00   0.2%          3          0
                                      MPI_sendrecv (00)       0.16  23.9%    3.6e+04   6.54e+07
                                     MPI_allreduce (00)       0.12  18.1%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 5] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.33  49.2%          1          0
                      SolveurSys::resoudre_systeme (01)       0.06   9.1%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.02   3.4%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.34 -50.4%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.04   6.0%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.09  13.0%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.03   4.3%         64   6.37e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.02   2.4%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.1%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.04   6.6%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.07  10.2%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.03   4.4%         80          0
                                      ijk_exchange (01)       0.04   5.5%        498          0
                                  calcul dtstab QC (01)       0.00   0.1%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.3%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.5%          3          0
                              qc diffusion vitesse (01)       0.00   0.5%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.3%          6          0
                                     qc rk3 update (01)       0.00   0.2%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.3%          3          0
                                 qc ajouter grad p (01)       0.00   0.3%          3          0
                                  jacobi_fine_mesh (01)       0.02   3.0%         52          0
                                        jacobi_all (01)       0.11  16.9%        166          0
                                     coarse_solver (01)       0.08  11.9%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.4%          3          0
                                        jacobi_all (01)       0.00   0.4%          3          0
                                      MPI_sendrecv (00)       0.13  19.7%    3.6e+04   6.54e+07
                                     MPI_allreduce (00)       0.11  16.8%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 6] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.33  49.2%          1          0
                      SolveurSys::resoudre_systeme (01)       0.06   9.1%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.02   3.5%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.34 -50.3%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.04   6.0%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.09  12.8%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.03   4.4%         64   6.37e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.02   2.4%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.04   6.5%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.07  10.0%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.03   4.1%         80          0
                                      ijk_exchange (01)       0.04   5.5%        498          0
                                  calcul dtstab QC (01)       0.00   0.1%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.3%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.5%          3          0
                              qc diffusion vitesse (01)       0.00   0.5%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.3%          6          0
                                     qc rk3 update (01)       0.00   0.2%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.3%          3          0
                                 qc ajouter grad p (01)       0.00   0.3%          3          0
                                  jacobi_fine_mesh (01)       0.02   3.0%         52          0
                                        jacobi_all (01)       0.11  17.0%        166          0
                                     coarse_solver (01)       0.08  11.9%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.4%          3          0
                                        jacobi_all (01)       0.00   0.4%          3          0
                                      MPI_sendrecv (00)       0.13  19.6%    3.6e+04   6.54e+07
                                     MPI_allreduce (00)       0.11  16.6%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 7] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.33  49.2%          1          0
                      SolveurSys::resoudre_systeme (01)       0.06   9.1%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.02   3.5%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.34 -50.6%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.04   6.1%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.09  13.2%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.03   4.6%         64   6.37e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.02   2.6%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.1%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.04   6.6%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.07   9.8%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.03   4.7%         80          0
                                      ijk_exchange (01)       0.04   5.4%        498          0
                                  calcul dtstab QC (01)       0.00   0.1%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.3%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.5%          3          0
                              qc diffusion vitesse (01)       0.00   0.5%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.3%          6          0
                                     qc rk3 update (01)       0.00   0.2%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.3%          3          0
                                 qc ajouter grad p (01)       0.00   0.3%          3          0
                                  jacobi_fine_mesh (01)       0.02   3.0%         52          0
                                        jacobi_all (01)       0.11  16.9%        166          0
                                     coarse_solver (01)       0.08  11.9%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.4%          3          0
                                        jacobi_all (01)       0.00   0.4%          3          0
                                      MPI_sendrecv (00)       0.14  20.5%    3.6e+04   6.54e+07
                                     MPI_allreduce (00)       0.11  16.4%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 8] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.33  49.2%          1          0
                      SolveurSys::resoudre_systeme (01)       0.06   9.1%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.02   3.5%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.34 -50.4%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.04   6.5%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.09  13.3%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.03   5.0%         64   6.47e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.02   2.5%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.04   6.1%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.06   9.6%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.03   4.4%         80          0
                                      ijk_exchange (01)       0.04   5.3%        498          0
                                  calcul dtstab QC (01)       0.00   0.1%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.3%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.5%          3          0
                              qc diffusion vitesse (01)       0.00   0.5%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.3%          6          0
                                     qc rk3 update (01)       0.00   0.2%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.3%          3          0
                                 qc ajouter grad p (01)       0.00   0.3%          3          0
                                  jacobi_fine_mesh (01)       0.02   3.0%         52          0
                                        jacobi_all (01)       0.11  16.9%        166          0
                                     coarse_solver (01)       0.08  11.9%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.4%          3          0
                                        jacobi_all (01)       0.00   0.4%          3          0
                                      MPI_sendrecv (00)       0.14  20.8%    3.6e+04   6.55e+07
                                     MPI_allreduce (00)       0.11  15.7%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 9] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.33  49.2%          1          0
                      SolveurSys::resoudre_systeme (01)       0.06   9.1%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.02   3.5%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.34 -51.1%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.04   6.4%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.09  13.8%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.03   5.0%         64   6.29e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.02   2.6%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.04   6.1%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.06   9.5%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.03   4.8%         80          0
                                      ijk_exchange (01)       0.04   5.6%        498          0
                                  calcul dtstab QC (01)       0.00   0.1%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.3%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.5%          3          0
                              qc diffusion vitesse (01)       0.00   0.5%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.3%          6          0
                                     qc rk3 update (01)       0.00   0.2%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.3%          3          0
                                 qc ajouter grad p (01)       0.00   0.3%          3          0
                                  jacobi_fine_mesh (01)       0.02   3.0%         52          0
                                        jacobi_all (01)       0.11  17.0%        166          0
                                     coarse_solver (01)       0.08  11.9%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.4%          3          0
                                        jacobi_all (01)       0.00   0.4%          3          0
                                      MPI_sendrecv (00)       0.14  21.4%    3.6e+04   6.55e+07
                                     MPI_allreduce (00)       0.11  15.6%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 10] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.33  49.2%          1          0
                      SolveurSys::resoudre_systeme (01)       0.06   9.1%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.02   3.5%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.34 -50.7%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.04   6.6%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.09  13.8%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.04   5.5%         64   6.29e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.02   2.6%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.04   5.9%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.06   8.9%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.03   4.4%         80          0
                                      ijk_exchange (01)       0.04   5.8%        498          0
                                  calcul dtstab QC (01)       0.00   0.1%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.3%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.5%          3          0
                              qc diffusion vitesse (01)       0.00   0.5%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.3%          6          0
                                     qc rk3 update (01)       0.00   0.2%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.3%          3          0
                                 qc ajouter grad p (01)       0.00   0.3%          3          0
                                  jacobi_fine_mesh (01)       0.02   2.9%         52          0
                                        jacobi_all (01)       0.11  16.9%        166          0
                                     coarse_solver (01)       0.08  11.9%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.4%          3          0
                                        jacobi_all (01)       0.00   0.4%          3          0
                                      MPI_sendrecv (00)       0.15  22.0%    3.6e+04   6.55e+07
                                     MPI_allreduce (00)       0.10  14.9%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 11] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.33  49.2%          1          0
                      SolveurSys::resoudre_systeme (01)       0.06   9.1%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.02   3.5%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.36 -53.5%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.04   6.3%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.13  18.6%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.04   5.9%         64   6.39e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.02   2.6%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.1%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.05   8.1%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.06   8.7%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.05   6.8%         80          0
                                      ijk_exchange (01)       0.05   8.1%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.2%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.2%          3          0
                              qc diffusion vitesse (01)       0.00   0.4%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.2%          6          0
                                     qc rk3 update (01)       0.00   0.1%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.2%          3          0
                                 qc ajouter grad p (01)       0.00   0.1%          3          0
                                  jacobi_fine_mesh (01)       0.02   2.5%         52          0
                                        jacobi_all (01)       0.11  16.9%        166          0
                                     coarse_solver (01)       0.08  11.9%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.1%          3          0
                                        jacobi_all (01)       0.00   0.2%          3          0
                                      MPI_sendrecv (00)       0.18  27.3%    3.6e+04   6.55e+07
                                     MPI_allreduce (00)       0.11  16.9%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 12] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.33  49.2%          1          0
                      SolveurSys::resoudre_systeme (01)       0.06   9.0%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.02   3.5%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.33 -49.3%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.04   6.4%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.07  10.7%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.04   6.0%         64   6.29e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.02   2.5%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.03   4.8%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.06   8.5%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.02   3.0%         80          0
                                      ijk_exchange (01)       0.03   4.4%        498          0
                                  calcul dtstab QC (01)       0.00   0.1%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.3%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.5%          3          0
                              qc diffusion vitesse (01)       0.00   0.6%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.3%          6          0
                                     qc rk3 update (01)       0.00   0.2%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.3%          3          0
                                 qc ajouter grad p (01)       0.00   0.3%          3          0
                                  jacobi_fine_mesh (01)       0.02   3.4%         52          0
                                        jacobi_all (01)       0.12  17.4%        166          0
                                     coarse_solver (01)       0.08  11.9%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.4%          3          0
                                        jacobi_all (01)       0.00   0.5%          3          0
                                      MPI_sendrecv (00)       0.13  19.3%    3.6e+04   6.55e+07
                                     MPI_allreduce (00)       0.09  13.3%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 13] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.33  49.2%          1          0
                      SolveurSys::resoudre_systeme (01)       0.06   9.0%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.02   3.5%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.34 -49.7%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.04   6.2%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.08  11.4%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.04   6.4%         64   6.29e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.02   2.7%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.04   5.2%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.05   8.0%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.02   3.3%         80          0
                                      ijk_exchange (01)       0.03   5.1%        498          0
                                  calcul dtstab QC (01)       0.00   0.1%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.3%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.5%          3          0
                              qc diffusion vitesse (01)       0.00   0.6%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.3%          6          0
                                     qc rk3 update (01)       0.00   0.2%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.2%          3          0
                                 qc ajouter grad p (01)       0.00   0.3%          3          0
                                  jacobi_fine_mesh (01)       0.02   3.1%         52          0
                                        jacobi_all (01)       0.12  17.2%        166          0
                                     coarse_solver (01)       0.08  11.9%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.4%          3          0
                                        jacobi_all (01)       0.00   0.4%          3          0
                                      MPI_sendrecv (00)       0.14  20.6%    3.6e+04   6.55e+07
                                     MPI_allreduce (00)       0.09  13.3%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 14] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.33  49.2%          1          0
                      SolveurSys::resoudre_systeme (01)       0.06   9.0%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.02   3.5%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.36 -52.7%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.04   6.6%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.12  18.2%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.04   6.6%         64   6.39e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.02   2.7%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.1%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.05   7.9%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.05   7.8%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.04   6.0%         80          0
                                      ijk_exchange (01)       0.05   8.1%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.2%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.2%          3          0
                              qc diffusion vitesse (01)       0.00   0.4%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.2%          6          0
                                     qc rk3 update (01)       0.00   0.1%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.2%          3          0
                                 qc ajouter grad p (01)       0.00   0.1%          3          0
                                  jacobi_fine_mesh (01)       0.02   2.5%         52          0
                                        jacobi_all (01)       0.11  16.9%        166          0
                                     coarse_solver (01)       0.08  11.9%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.1%          3          0
                                        jacobi_all (01)       0.00   0.2%          3          0
                                      MPI_sendrecv (00)       0.19  27.6%    3.6e+04   6.55e+07
                                     MPI_allreduce (00)       0.11  15.7%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 15] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.33  49.2%          1          0
                      SolveurSys::resoudre_systeme (01)       0.06   9.0%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.02   3.5%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.36 -53.3%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.04   6.7%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.13  18.7%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.05   6.9%         64   6.39e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.02   2.6%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.1%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.05   7.7%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.05   7.7%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.04   6.4%         80          0
                                      ijk_exchange (01)       0.06   8.2%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.2%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.2%          3          0
                              qc diffusion vitesse (01)       0.00   0.4%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.2%          6          0
                                     qc rk3 update (01)       0.00   0.1%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.2%          3          0
                                 qc ajouter grad p (01)       0.00   0.1%          3          0
                                  jacobi_fine_mesh (01)       0.02   2.5%         52          0
                                        jacobi_all (01)       0.11  17.0%        166          0
                                     coarse_solver (01)       0.08  11.9%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.1%          3          0
                                        jacobi_all (01)       0.00   0.2%          3          0
                                      MPI_sendrecv (00)       0.19  28.3%    3.6e+04   6.55e+07
                                     MPI_allreduce (00)       0.10  15.4%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 16] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.33  49.2%          1          0
                      SolveurSys::resoudre_systeme (01)       0.06   9.0%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.02   3.5%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.34 -50.7%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.05   7.2%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.10  14.2%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.05   7.0%         64   6.39e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.02   2.6%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.03   5.0%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.05   7.4%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.03   4.6%         80          0
                                      ijk_exchange (01)       0.04   5.4%        498          0
                                  calcul dtstab QC (01)       0.00   0.1%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.3%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.5%          3          0
                              qc diffusion vitesse (01)       0.00   0.5%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.3%          6          0
                                     qc rk3 update (01)       0.00   0.2%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.3%          3          0
                                 qc ajouter grad p (01)       0.00   0.3%          3          0
                                  jacobi_fine_mesh (01)       0.02   3.0%         52          0
                                        jacobi_all (01)       0.11  17.0%        166          0
                                     coarse_solver (01)       0.08  11.9%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.4%          3          0
                                        jacobi_all (01)       0.00   0.5%          3          0
                                      MPI_sendrecv (00)       0.16  23.9%    3.6e+04   6.55e+07
                                     MPI_allreduce (00)       0.08  12.4%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 17] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       0.67 100.0%          0          0
                                   Preparer calcul (01)       0.67 100.0%          0          0
                          Resoudre (timestep loop) (01)       0.33  49.2%          1          0
                      SolveurSys::resoudre_systeme (01)       0.06   9.1%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.02   3.5%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -0.34 -50.3%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.05   6.8%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.09  13.6%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.05   7.3%         64   6.49e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.02   2.4%        294   3.76e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          4         16
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.04   5.9%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%          1          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         10          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          0          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          0          0
                         MPI_allreduce:MPI_barrier (02)       0.05   7.2%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.03   4.2%         80          0
                                      ijk_exchange (01)       0.04   5.6%        498          0
                                  calcul dtstab QC (01)       0.00   0.1%          1          0
                               update statistiques (01)       0.01   0.9%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.04   6.1%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.00   0.3%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.00   0.5%          3          0
                              qc diffusion vitesse (01)       0.00   0.5%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.00   0.3%          6          0
                                     qc rk3 update (01)       0.00   0.2%          3          0
                          qc diffusion temperature (01)       0.00   0.1%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.00   0.1%          3          0
                               qc resoudre systeme (01)       0.24  35.2%          3          0
                                 qc ajouter grad p (01)       0.00   0.2%          3          0
                                  jacobi_fine_mesh (01)       0.02   3.0%         52          0
                                        jacobi_all (01)       0.11  16.9%        166          0
                                     coarse_solver (01)       0.08  11.9%         20          0
                                  jacobi_fine_mesh (01)       0.00   0.3%          3          0
                                        jacobi_all (01)       0.00   0.4%          3          0
                                      MPI_sendrecv (00)       0.16  23.4%    3.6e+04   6.55e+07
                                     MPI_allreduce (00)       0.09  13.2%   2.44e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics overall execution : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)   avg time         (min,max)                avg count (min,max)                  avg qty (min,max)              
                                       Temps total (01)       5.75  100.0% (5.75,5.75)                      0                                  0                        
                                   Preparer calcul (01)       5.75  100.0% (5.75,5.75)                      0                                  0                        
                          Resoudre (timestep loop) (01)       5.51   95.8% (5.51,5.51)                      1                                  0                        
                      SolveurSys::resoudre_systeme (01)       0.15    2.6% (0.15,0.15)                     20                           1.04e+03                        
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
                               Pb_base::postraiter (01)       0.06    1.0% (0.05,0.07)                      2                                  0                        
                             Probleme_base::sauver (01)       0.00    0.0%                                  0                                  0                        
                                         temporary (01)       0.00    0.0%                                  0                                  0                        
                             Assembleur::assembler (01)       0.00    0.0%                                  0                                  0                        
                                            Divers (01)      -5.84  -101.5% (-5.85,-5.82)                    0                                  0                        
                                                m1 (01)       0.00    0.0%                                  0                                  0                        
                                                m2 (01)       0.00    0.0%                                  0                                  0                        
                                                m3 (01)       0.00    0.0%                                  0                                  0                        
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.05    0.9% (0.05,0.06)               3.22e+03                                  0                        
                        MPI_sendrecv:MPI_send_recv (02)       0.20    3.4% (0.17,0.21)               3.56e+04                           6.48e+07 (6.48e+07,6.48e+07)    
                             MPI_sendrecv:MPI_send (02)       0.07    1.3% (0.00,0.11)                    145 (34,152)                  6.01e+05 (442,6.51e+05)         
                             MPI_sendrecv:MPI_recv (02)       0.00    0.0% (0.00,0.04)                    145 (2,2.58e+03)              6.01e+05 (26,1.08e+07)          
                        MPI_sendrecv:MPI_broadcast (02)       0.04    0.7% (0.02,0.06)                    618                           5.05e+03                        
                         MPI_sendrecv:MPI_alltoall (02)       0.00    0.0% (0.00,0.00)                      7                                 28                        
                      MPI_allreduce:MPI_partialsum (02)       0.00    0.0%                                  0                                  0                        
                       MPI_allreduce:MPI_sumdouble (02)       0.12    2.0% (0.09,0.13)               2.42e+03                                  0                        
                       MPI_allreduce:MPI_mindouble (02)       0.00    0.0% (0.00,0.00)                     10                                  0                        
                       MPI_allreduce:MPI_maxdouble (02)       0.00    0.0% (0.00,0.00)                      9                                  0                        
                          MPI_allreduce:MPI_sumint (02)       0.00    0.0% (0.00,0.00)                     13                                  0                        
                          MPI_allreduce:MPI_minint (02)       0.00    0.1% (0.00,0.01)                      7                                  0                        
                          MPI_allreduce:MPI_maxint (02)       0.00    0.0% (0.00,0.00)                      7                                  0                        
                         MPI_allreduce:MPI_barrier (02)       0.05    0.8% (0.00,0.08)                     12                                  0                        
                               io:MPI_send_recv_io (02)       0.00    0.0%                                  0                                  0                        
                                           Scatter (02)       0.00    0.0%                                  0                                  0                        
                                      ijk_exchange (01)       0.05    0.9% (0.04,0.06)                     80                                  0                        
                                      ijk_exchange (01)       0.29    5.1% (0.27,0.30)                    498                                  0                        
                                  calcul dtstab QC (01)       0.00    0.0% (0.00,0.00)                      1                                  0                        
                               update statistiques (01)       0.03    0.5% (0.03,0.03)                      1                                  0                        
                         calcul terme acceleration (01)       0.00    0.0% (0.00,0.00)                      1                                  0                        
                                     checkpointing (01)       0.05    0.9% (0.05,0.05)                      1                                  0                        
                                         TF update (01)       0.00    0.0% (0.00,0.00)                      1                                  0                        
                                 qc convection rho (01)       0.09    1.6% (0.09,0.09)                      3                                  0                        
                     qc calcul pthermo,t,rho_v,etc (01)       0.01    0.2% (0.01,0.01)                      3                                  0                        
                              qc diffusion vitesse (01)       0.13    2.2% (0.13,0.13)                      3                                  0                        
                   qc diffusion turbulente vitesse (01)       0.00    0.0%                                  0                                  0                        
            qc diffusion vitesse modele structurel (01)       0.00    0.0%                                  0                                  0                        
                             qc convection vitesse (01)       0.15    2.5% (0.14,0.15)                      6                                  0                        
                                     qc rk3 update (01)       0.02    0.4% (0.02,0.02)                      3                                  0                        
                          qc diffusion temperature (01)       0.05    0.9% (0.05,0.05)                      3                                  0                        
                          qc turbulent diffusivity (01)       0.00    0.0%                                  0                                  0                        
               qc turbulent diffusivity structurel (01)       0.00    0.0%                                  0                                  0                        
                                    qc prepare rhs (01)       0.01    0.3% (0.01,0.02)                      3                                  0                        
                               qc resoudre systeme (01)       4.86   84.5% (4.86,4.86)                      3                                  0                        
                                 qc ajouter grad p (01)       0.02    0.4% (0.02,0.02)                      3                                  0                        
                                  jacobi_fine_mesh (01)       0.92   16.1% (0.92,0.93)                     52                                  0                        
                                        jacobi_all (01)       4.25   73.9% (4.23,4.27)                    166                                  0                        
                                     coarse_solver (01)       0.15    2.7% (0.15,0.15)                     20                                  0                        
                                  jacobi_fine_mesh (01)       0.02    0.3% (0.02,0.02)                      3                                  0                        
                                        jacobi_all (01)       0.02    0.3% (0.02,0.02)                      3                                  0                        
                                      MPI_sendrecv (00)       0.32    5.5% (0.24,0.35)               3.66e+04 (3.64e+04,3.89e+04)        6.6e+07 (6.54e+07,7.56e+07)    
                                     MPI_allreduce (00)       0.17    2.9% (0.12,0.21)               2.48e+03                                  0                        
                                                io (00)       0.00    0.0%                                  0                                  0                        

Statistics execution [PE 0] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.9%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.7%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.07   1.2%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.84 -101.5%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.05   0.9%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.18   3.2%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.00   0.0%         34        442
                             MPI_sendrecv:MPI_recv (02)       0.04   0.7%   2.58e+03   1.08e+07
                        MPI_sendrecv:MPI_broadcast (02)       0.02   0.3%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.11   1.9%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.00   0.0%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.05   0.8%         80          0
                                      ijk_exchange (01)       0.28   4.9%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.2%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.14   2.5%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.01   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.4%          3          0
                                  jacobi_fine_mesh (01)       0.93  16.1%         52          0
                                        jacobi_all (01)       4.26  74.0%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.24   4.2%   3.89e+04   7.56e+07
                                     MPI_allreduce (00)       0.12   2.0%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 1] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.8%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.6%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.06   1.0%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.84 -101.6%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.05   0.9%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.19   3.4%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.04   0.6%        152   6.28e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.06   1.1%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.12   2.0%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.08   1.3%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.05   0.9%         80          0
                                      ijk_exchange (01)       0.29   5.1%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.2%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.14   2.5%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.01   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.4%          3          0
                                  jacobi_fine_mesh (01)       0.92  16.1%         52          0
                                        jacobi_all (01)       4.25  74.0%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.30   5.2%   3.64e+04   6.54e+07
                                     MPI_allreduce (00)       0.20   3.4%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 2] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.8%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.6%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.06   1.0%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.84 -101.5%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.05   0.9%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.21   3.6%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.04   0.6%        152   6.38e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.05   0.9%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.13   2.2%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.01   0.1%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.07   1.3%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.06   1.1%         80          0
                                      ijk_exchange (01)       0.30   5.2%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.2%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.14   2.5%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.01   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.3%          3          0
                                  jacobi_fine_mesh (01)       0.92  16.0%         52          0
                                        jacobi_all (01)       4.24  73.8%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.30   5.1%   3.64e+04   6.54e+07
                                     MPI_allreduce (00)       0.21   3.6%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 3] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.8%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.6%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.06   1.0%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.82 -101.3%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.05   0.9%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.17   2.9%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.06   1.0%        152   6.28e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.05   0.9%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.11   1.9%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.07   1.2%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.05   0.9%         80          0
                                      ijk_exchange (01)       0.27   4.7%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.2%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.15   2.5%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.01   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.4%          3          0
                                  jacobi_fine_mesh (01)       0.92  16.1%         52          0
                                        jacobi_all (01)       4.26  74.1%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.28   4.8%   3.64e+04   6.54e+07
                                     MPI_allreduce (00)       0.18   3.1%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 4] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.8%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.6%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.06   1.0%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.82 -101.2%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.06   1.0%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.20   3.6%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.06   1.1%        152   6.28e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.05   0.8%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.13   2.3%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.07   1.2%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.05   0.9%         80          0
                                      ijk_exchange (01)       0.30   5.1%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.2%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.14   2.5%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.02   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.4%          3          0
                                  jacobi_fine_mesh (01)       0.92  16.0%         52          0
                                        jacobi_all (01)       4.23  73.6%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.31   5.5%   3.64e+04   6.54e+07
                                     MPI_allreduce (00)       0.20   3.5%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 5] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.8%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.6%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.06   1.0%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.83 -101.5%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.05   0.9%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.20   3.5%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.07   1.1%        152   6.38e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.05   0.8%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.13   2.2%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.06   1.1%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.05   0.9%         80          0
                                      ijk_exchange (01)       0.29   5.1%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.2%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.14   2.5%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.01   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.3%          3          0
                                  jacobi_fine_mesh (01)       0.92  16.1%         52          0
                                        jacobi_all (01)       4.24  73.8%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.32   5.5%   3.64e+04   6.54e+07
                                     MPI_allreduce (00)       0.19   3.3%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 6] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.8%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.6%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.06   1.0%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.85 -101.7%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.05   0.9%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.21   3.7%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.07   1.2%        152   6.38e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.05   0.9%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.12   2.1%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.01   0.1%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.06   1.1%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.06   1.1%         80          0
                                      ijk_exchange (01)       0.30   5.2%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.2%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.14   2.5%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.01   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.3%          3          0
                                  jacobi_fine_mesh (01)       0.93  16.1%         52          0
                                        jacobi_all (01)       4.25  73.9%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.33   5.8%   3.64e+04   6.54e+07
                                     MPI_allreduce (00)       0.19   3.2%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 7] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.8%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.6%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.06   1.0%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.84 -101.5%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.05   0.9%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.21   3.6%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.06   1.1%        152   6.38e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.03   0.6%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.12   2.1%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.1%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.06   1.0%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.06   1.0%         80          0
                                      ijk_exchange (01)       0.30   5.2%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.2%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.14   2.5%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.01   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.4%          3          0
                                  jacobi_fine_mesh (01)       0.92  15.9%         52          0
                                        jacobi_all (01)       4.25  73.9%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.30   5.3%   3.64e+04   6.54e+07
                                     MPI_allreduce (00)       0.19   3.2%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 8] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.8%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.6%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.06   1.0%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.83 -101.5%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.05   0.9%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.20   3.5%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.08   1.5%        152   6.48e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.04   0.7%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.12   2.1%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.1%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.05   0.9%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.06   1.0%         80          0
                                      ijk_exchange (01)       0.29   5.0%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.2%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.14   2.5%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.01   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.3%          3          0
                                  jacobi_fine_mesh (01)       0.93  16.1%         52          0
                                        jacobi_all (01)       4.25  73.9%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.32   5.7%   3.64e+04   6.55e+07
                                     MPI_allreduce (00)       0.18   3.1%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 9] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.8%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.6%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.05   1.0%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.85 -101.7%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.05   0.9%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.21   3.7%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.07   1.2%        152    6.3e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.05   0.9%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.12   2.0%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.1%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.05   0.9%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.06   1.0%         80          0
                                      ijk_exchange (01)       0.30   5.2%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.2%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.14   2.5%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.01   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.4%          3          0
                                  jacobi_fine_mesh (01)       0.92  16.0%         52          0
                                        jacobi_all (01)       4.25  73.9%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.34   5.9%   3.64e+04   6.55e+07
                                     MPI_allreduce (00)       0.17   3.0%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 10] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.8%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.6%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.05   1.0%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.84 -101.5%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.05   1.0%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.20   3.4%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.09   1.5%        152    6.3e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.03   0.6%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.11   2.0%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.05   0.8%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.05   0.9%         80          0
                                      ijk_exchange (01)       0.29   5.1%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.3%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.15   2.5%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.01   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.4%          3          0
                                  jacobi_fine_mesh (01)       0.92  16.1%         52          0
                                        jacobi_all (01)       4.25  73.9%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.32   5.5%   3.64e+04   6.55e+07
                                     MPI_allreduce (00)       0.16   2.9%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 11] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.8%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.6%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.05   1.0%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.85 -101.7%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.05   0.9%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.20   3.5%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.09   1.5%        152    6.4e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.03   0.6%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.11   1.9%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.05   0.8%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.05   0.9%         80          0
                                      ijk_exchange (01)       0.30   5.2%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.2%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.15   2.5%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.01   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.4%          3          0
                                  jacobi_fine_mesh (01)       0.92  16.0%         52          0
                                        jacobi_all (01)       4.26  74.0%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.32   5.6%   3.64e+04   6.55e+07
                                     MPI_allreduce (00)       0.16   2.7%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 12] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.8%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.6%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.05   1.0%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.84 -101.6%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.06   1.0%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.20   3.5%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.10   1.7%        152    6.3e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.04   0.7%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.11   1.9%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.04   0.7%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.05   0.9%         80          0
                                      ijk_exchange (01)       0.30   5.1%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.2%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.15   2.5%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.01   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.4%          3          0
                                  jacobi_fine_mesh (01)       0.92  16.0%         52          0
                                        jacobi_all (01)       4.26  74.0%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.34   5.9%   3.64e+04   6.55e+07
                                     MPI_allreduce (00)       0.15   2.7%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 13] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.8%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.6%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.05   1.0%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.84 -101.5%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.05   0.9%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.20   3.5%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.10   1.7%        152    6.3e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.04   0.7%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.11   2.0%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.04   0.7%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.05   0.9%         80          0
                                      ijk_exchange (01)       0.30   5.2%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.2%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.15   2.5%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.01   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.4%          3          0
                                  jacobi_fine_mesh (01)       0.92  16.0%         52          0
                                        jacobi_all (01)       4.25  73.9%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.34   5.9%   3.64e+04   6.55e+07
                                     MPI_allreduce (00)       0.16   2.7%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 14] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.8%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.6%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.06   1.0%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.83 -101.4%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.05   0.9%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.20   3.5%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.10   1.8%        152    6.4e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.04   0.7%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.12   2.2%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.04   0.6%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.05   0.9%         80          0
                                      ijk_exchange (01)       0.30   5.1%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.2%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.15   2.5%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.01   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.4%          3          0
                                  jacobi_fine_mesh (01)       0.92  16.0%         52          0
                                        jacobi_all (01)       4.25  73.8%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.34   6.0%   3.64e+04   6.55e+07
                                     MPI_allreduce (00)       0.16   2.8%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 15] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.8%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.6%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.05   1.0%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.83 -101.4%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.05   1.0%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.20   3.5%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.10   1.7%        152    6.4e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.05   0.9%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.12   2.1%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.03   0.6%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.05   0.9%         80          0
                                      ijk_exchange (01)       0.30   5.1%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.2%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.15   2.5%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.01   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.4%          3          0
                                  jacobi_fine_mesh (01)       0.92  16.0%         52          0
                                        jacobi_all (01)       4.25  73.8%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.35   6.1%   3.64e+04   6.55e+07
                                     MPI_allreduce (00)       0.16   2.8%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 16] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.8%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.6%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.06   1.0%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.85 -101.8%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.06   1.0%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.17   3.0%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.11   1.8%        152    6.4e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.03   0.6%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.09   1.5%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.03   0.5%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.04   0.7%         80          0
                                      ijk_exchange (01)       0.28   4.8%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.3%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.15   2.6%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.01   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.4%          3          0
                                  jacobi_fine_mesh (01)       0.93  16.2%         52          0
                                        jacobi_all (01)       4.27  74.3%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.31   5.5%   3.64e+04   6.55e+07
                                     MPI_allreduce (00)       0.12   2.1%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

Statistics execution [PE 17] : Statistiques de post resolution
 the time was measured by the following method :
 gettimeofday (elapsed CPU time)
                                        Counter (level)    time(s)             count   quantity
                                       Temps total (01)       5.75 100.0%          0          0
                                   Preparer calcul (01)       5.75 100.0%          0          0
                          Resoudre (timestep loop) (01)       5.51  95.8%          1          0
                      SolveurSys::resoudre_systeme (01)       0.15   2.6%         20   1.04e+03
                    Solveurpetsc::resoudre_systeme (01)       0.00   0.0%          0          0
        Equation_base::Gradient_conjugue_diff_impl (01)       0.00   0.0%          0          0
                  Operateur::calculer_pas_de_temps (01)       0.00   0.0%          0          0
                  ModeleTurbulence*::mettre_a_jour (01)       0.00   0.0%          0          0
                  Operateur_Conv::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Diff::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Decr::ajouter/calculer (01)       0.00   0.0%          0          0
                  Operateur_Grad::ajouter/calculer (01)       0.00   0.0%          0          0
                   Operateur_Div::ajouter/calculer (01)       0.00   0.0%          0          0
                          Source::ajouter/calculer (01)       0.00   0.0%          0          0
                               Pb_base::postraiter (01)       0.05   1.0%          2          0
                             Probleme_base::sauver (01)       0.00   0.0%          0          0
                                         temporary (01)       0.00   0.0%          0          0
                             Assembleur::assembler (01)       0.00   0.0%          0          0
                                            Divers (01)      -5.84 -101.6%          0          0
                                                m1 (01)       0.00   0.0%          0          0
                                                m2 (01)       0.00   0.0%          0          0
                                                m3 (01)       0.00   0.0%          0          0
        DoubleVect/IntVect::echange_espace_virtuel (02)       0.06   1.0%   3.22e+03          0
                        MPI_sendrecv:MPI_send_recv (02)       0.20   3.4%   3.56e+04   6.48e+07
                             MPI_sendrecv:MPI_send (02)       0.10   1.7%        152   6.51e+05
                             MPI_sendrecv:MPI_recv (02)       0.00   0.0%          2         26
                        MPI_sendrecv:MPI_broadcast (02)       0.04   0.6%        618   5.05e+03
                         MPI_sendrecv:MPI_alltoall (02)       0.00   0.0%          7         28
                      MPI_allreduce:MPI_partialsum (02)       0.00   0.0%          0          0
                       MPI_allreduce:MPI_sumdouble (02)       0.10   1.8%   2.42e+03          0
                       MPI_allreduce:MPI_mindouble (02)       0.00   0.0%         10          0
                       MPI_allreduce:MPI_maxdouble (02)       0.00   0.0%          9          0
                          MPI_allreduce:MPI_sumint (02)       0.00   0.0%         13          0
                          MPI_allreduce:MPI_minint (02)       0.00   0.0%          7          0
                          MPI_allreduce:MPI_maxint (02)       0.00   0.0%          7          0
                         MPI_allreduce:MPI_barrier (02)       0.03   0.5%         12          0
                               io:MPI_send_recv_io (02)       0.00   0.0%          0          0
                                           Scatter (02)       0.00   0.0%          0          0
                                      ijk_exchange (01)       0.05   0.9%         80          0
                                      ijk_exchange (01)       0.29   5.0%        498          0
                                  calcul dtstab QC (01)       0.00   0.0%          1          0
                               update statistiques (01)       0.03   0.5%          1          0
                         calcul terme acceleration (01)       0.00   0.0%          1          0
                                     checkpointing (01)       0.05   0.9%          1          0
                                         TF update (01)       0.00   0.0%          1          0
                                 qc convection rho (01)       0.09   1.6%          3          0
                     qc calcul pthermo,t,rho_v,etc (01)       0.01   0.2%          3          0
                              qc diffusion vitesse (01)       0.13   2.3%          3          0
                   qc diffusion turbulente vitesse (01)       0.00   0.0%          0          0
            qc diffusion vitesse modele structurel (01)       0.00   0.0%          0          0
                             qc convection vitesse (01)       0.15   2.5%          6          0
                                     qc rk3 update (01)       0.02   0.4%          3          0
                          qc diffusion temperature (01)       0.05   0.9%          3          0
                          qc turbulent diffusivity (01)       0.00   0.0%          0          0
               qc turbulent diffusivity structurel (01)       0.00   0.0%          0          0
                                    qc prepare rhs (01)       0.01   0.3%          3          0
                               qc resoudre systeme (01)       4.86  84.5%          3          0
                                 qc ajouter grad p (01)       0.02   0.4%          3          0
                                  jacobi_fine_mesh (01)       0.93  16.1%         52          0
                                        jacobi_all (01)       4.26  74.0%        166          0
                                     coarse_solver (01)       0.15   2.7%         20          0
                                  jacobi_fine_mesh (01)       0.02   0.3%          3          0
                                        jacobi_all (01)       0.02   0.3%          3          0
                                      MPI_sendrecv (00)       0.33   5.8%   3.64e+04   6.55e+07
                                     MPI_allreduce (00)       0.13   2.3%   2.48e+03          0
                                                io (00)       0.00   0.0%          0          0

