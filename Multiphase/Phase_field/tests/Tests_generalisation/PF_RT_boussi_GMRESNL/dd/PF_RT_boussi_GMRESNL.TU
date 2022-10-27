Statistiques d'initialisation du calcul

Temps total                       0.013386

Statistiques de resolution du probleme

Temps total                       0.132794


Timesteps                         10
Secondes / pas de temps           0.0132704
Dont solveurs Ax=B                0.0070309 52% (3 appels/pas de temps)
Dont mettre_a_jour                0.0001706 1% (1 appel/pas de temps)
Dont operateurs convection        0.0014137 10% (6 appels/pas de temps)
Dont operateurs diffusion         0.0003448 2% (6 appels/pas de temps)
Dont operateurs gradient          0.0001967 1% (11 appels/pas de temps)
Dont operateurs divergence        5.17e-05 0% (6 appels/pas de temps)
Dont operateurs source            0.0001827 1% (9 appels/pas de temps)
Dont operations postraitement     0.0001959 1% (1 appel/pas de temps)
Dont calcul dt                    0.0001348 1% (12 appels/pas de temps)
Dont calcul divers                0.0035486 26% (0 appels/pas de temps)
Nb echange_espace_virtuel / pas de temps 771.2
Nb solveur / pas de temps         3
Secondes / solveur                0.00234363
Iterations / solveur              65.9667
I/O:

Timesteps = number of time steps
Nb solveur = number of linear system resolutions
Nb assemblage implicite = number of matrix assemblies for the implicit scheme
Iterations = average number of iterations of the solver
Communications = fraction of the time spent
                 in communications between processors (excluding io files)
Network latency = time of one mpsum measured by an internal bench over 0.1s
Network bandwidth = maximum on all processors
                    of the average bandwidth of send_recv operations
Waiting time = estimation of the waiting time of the different processors

Max_waiting_time big    => probably due to a bad partitioning
Communications > 30%    => too many processors or network too slow

Statistiques de post resolution

Temps total                       0.015807

