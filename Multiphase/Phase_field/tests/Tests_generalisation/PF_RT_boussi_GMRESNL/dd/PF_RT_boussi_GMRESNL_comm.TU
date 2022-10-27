*************************************************************************************************
Statistiques d'initialisation du calcul
Block divided into 2 parts:
Part 1: for each function, gives the percentage consumed by communications
Part 2: for each type of communication, gives the percentage of time spent in each function
        Note that some counters are nested within each other, 
        that's why the sum of the percentages does not necessarily add up to 100% in this part.
*************************************************************************************************


-----------------------------------------------------------------
DoubleVect/IntVect::echange_espace_virtuel
-----------------------------------------------------------------
Temps total : 4.2e-05s
	dont Preparer calcul : 4.2e-05s (100%) 
	dont SolveurSys::resoudre_systeme : 3e-06s (7%) 
	dont Operateur_Grad::ajouter/calculer : 2e-06s (5%) 
	dont Pb_base::postraiter : 5e-06s (12%) 


-----------------------------------------------------------------
                    Famille MPI_allreduce                         
-----------------------------------------------------------------
Temps total : 0s

-----------------------------------------------------------------
                    Famille MPI_sendrecv                           
-----------------------------------------------------------------
Temps total : 0s


*************************************************************************************************
Statistiques de resolution du probleme
Block divided into 2 parts:
Part 1: for each function, gives the percentage consumed by communications
Part 2: for each type of communication, gives the percentage of time spent in each function
        Note that some counters are nested within each other, 
        that's why the sum of the percentages does not necessarily add up to 100% in this part.
*************************************************************************************************


-----------------------------------------------------------------
DoubleVect/IntVect::echange_espace_virtuel
-----------------------------------------------------------------
Temps total : 0.00149s
	dont Resoudre (timestep loop) : 0.00149s (100%) 
	dont SolveurSys::resoudre_systeme : 0.001082s (73%) 
	dont Operateur_Grad::ajouter/calculer : 5.3e-05s (4%) 
	dont Source::ajouter/calculer : 1.6e-05s (1%) 
	dont ::mettre_a_jour : 9e-06s (1%) 


-----------------------------------------------------------------
                    Famille MPI_allreduce                         
-----------------------------------------------------------------
Temps total : 0s

-----------------------------------------------------------------
                    Famille MPI_sendrecv                           
-----------------------------------------------------------------
Temps total : 0s


*************************************************************************************************
Statistiques de post resolution
Block divided into 2 parts:
Part 1: for each function, gives the percentage consumed by communications
Part 2: for each type of communication, gives the percentage of time spent in each function
        Note that some counters are nested within each other, 
        that's why the sum of the percentages does not necessarily add up to 100% in this part.
*************************************************************************************************



-----------------------------------------------------------------
                    Famille MPI_allreduce                         
-----------------------------------------------------------------
Temps total : 0s

-----------------------------------------------------------------
                    Famille MPI_sendrecv                           
-----------------------------------------------------------------
Temps total : 0s


