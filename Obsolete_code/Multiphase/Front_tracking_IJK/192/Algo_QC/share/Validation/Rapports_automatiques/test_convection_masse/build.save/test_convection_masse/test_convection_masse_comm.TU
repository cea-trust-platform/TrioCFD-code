*************************************************************************************************
Statistiques de post resolution
Block divided into 2 parts:
Part 1: for each function, gives the percentage consumed by communications
Part 2: for each type of communication, gives the percentage of time spent in each function
        Note that some counters are nested within each other, 
        that's why the sum of the percentages does not necessarily add up to 100% in this part.
*************************************************************************************************


-----------------------------------------------------------------
DoubleVect/IntVect::echange_espace_virtuel
-----------------------------------------------------------------
Temps total : 0.039105s
	dont Preparer calcul : 0.039105s (100%) 
	dont Resoudre (timestep loop) : 0.035821s (92%) 
	dont SolveurSys::resoudre_systeme : 0.000422s (1%) 
	dont qc diffusion vitesse : 0.001212s (3%) 
	dont qc resoudre systeme : 0.026672s (68%) 
	dont projection: multigrille : 0.022514s (58%) 
	dont multigrille : residu jacobi : 0.02032s (52%) 
	dont multigrille : residu jacobi : 0.000672s (2%) 


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
DoubleVect/IntVect::echange_espace_virtuel
-----------------------------------------------------------------
Temps total : 0.038357s
	dont Preparer calcul : 0.038357s (100%) 
	dont Resoudre (timestep loop) : 0.035419s (92%) 
	dont SolveurSys::resoudre_systeme : 0.000489s (1%) 
	dont qc diffusion vitesse : 0.001143s (3%) 
	dont qc resoudre systeme : 0.02666s (70%) 
	dont projection: multigrille : 0.022535s (59%) 
	dont multigrille : residu jacobi : 0.020301s (53%) 
	dont multigrille : residu jacobi : 0.000655s (2%) 


-----------------------------------------------------------------
                    Famille MPI_allreduce                         
-----------------------------------------------------------------
Temps total : 0s

-----------------------------------------------------------------
                    Famille MPI_sendrecv                           
-----------------------------------------------------------------
Temps total : 0s


