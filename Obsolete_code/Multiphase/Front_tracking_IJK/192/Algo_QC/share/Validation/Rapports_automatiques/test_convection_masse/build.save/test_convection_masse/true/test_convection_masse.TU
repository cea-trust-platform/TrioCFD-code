Statistiques de post resolution

Temps total                       1.71659

Timesteps                         1
Secondes / pas de temps           1.38135
Dont solveurs Ax=B                0.114127 8% (20 appels/pas de temps)
Dont operations postraitement     0.061226 4% (2 appels/pas de temps)
Dont calcul divers                -1.40995 -102%
Nb echange_espace_virtuel / pas de temps 1999
Nb solveur / pas de temps         20
Secondes / solveur                0.00570635
Iterations / solveur              31.65

Timesteps = nombre de pas de temps
Nb solveur = nombre de resolutions de systeme lineaire
Nb assemblage implicite = nombre d'assemblage de matrice pour le schema implicite
Iterations = nombre moyen d'iterations du solveur
Communications = fraction du temps passe dans les
                 communications entre processeurs (hors io fichiers)
Network latency = temps d'un mpsum mesure par un bench interne sur 0.1s
Network bandwidth = maximum sur l'ensemble des processeurs
                    de la bande passante moyenne des operations send_recv
Waiting time = estimation du temps d'attente des differents processeurs

Max_waiting_time grand  => probablement mauvais decoupage
Communications > 30%    => trop de processeurs ou reseau trop lent

