#########################
### Author : ES243900 ###
###      Mai 2020     ###
#########################

This README file is to present the latest options/possibilities integrated in the
 Discontinuous Ftont Tracking (FTD) Baltik of TrioCFD, V1.8.1.

==================================
| Introduction and configuration |
==================================

Steps to follow in order to use this baltik :

	1- Source env_TRUST.sh
	2- baltik_build_configure -execute
	3- make optim debug
	4- source env_Front_tracking_discontinu.sh 

Once done, the user can test the existing test case by running
	
	1- make check_optim (for test cases found in tests/Reference)
	2- make check_deps_optim (for all test cases with dependency to FTD)

================================
| NEW options that can be used |
================================

1. Wall Hydraulic standard law 

	It is possible to run a simulation by prescribing a wall boundary law.
	To use, just put "turbulence_paroi loi_standard_hydr" in the block devoted to
	"modele_turbulence" in "Navier_Stokes_FT_Disc" reading zone. 
	
	This law can be used for both VDF and VEF discretizations. If the medium defined
	in the data file is "FLUIDE_DIPHASIQUE", the law takes into account the 
	viscosity of each phase, thanks to the indicator function. Otherwise, the medium
	if "FLUIDE_INCOMPRESSIBLE" and thus the standard law is applied.

	Recall that for FTD works with LES and not RANS. Moreover, it is advised to use the
	WALE SGS model (sous_maille_wale).

2- Turbulent Convection/Diffusion Concentration problem
	
	It is possible to define a laminar or a turbulent Convection/Diffusion 
	Concentration equation. For a laminar problem, the user should define simply the
	type "Convection_Diffusion_Concentration_FT_Disc".
	
	For a turbulent resolution, "Convection_Diffusion_Concentration_Turbulent_FT_Disc"
	is to be used. The constant Schmidt SGS model is to be used.
	 ATTENTION : in the last case, defining the turbulence model is a
	MUST and not a CHOICE. Otherwise, the code stops !

3- Different reaction rate formulations when dealing with instantaneous reactive flows

	DFT offers 4 different reaction rate models. To identify what model to be used,
	the user should define the "modele_cinetique" keyword in the zone where the
	Convection/Diffusion Concentration problem is read. The keyword "modele_cinetique" 
	should be followed by an unsigned integer 1 to 4, depending on which model is 
	desired.

	> modele_cinetique 1 : means that the reaction model to be used follows a simple-rate 
	formulation (omega = alpha C_A*C_B, for two reactive species A and B). This model can
	be used, both for a laminar and a turbulent Convection/Diffusion Concentration problem.

	> modele_cinetique 2 : this is an LES diffusive model, known by eddy disspertion 
	concept (EDC). Refer to  Bertrand et al. 2016 Cemical Engineering Journal for further 
	information. This model requires the SGS diffusion and so to be used ONLY WITH A 
	TURBULENT Convection/Diffusion Concentration problem. Otherwise, DFT stops !
	
	> modele_cinetique 3 : this is an LES variance model dependent on the SGS viscocity. 
	It can be used with both laminar and turbulent Convection/Diffusion Concentration 
	problem because "Navier_Stokes_FT_Disc" is turbulent by nature (it is generic) !

	> modele_cinetique 4 : mix between model 2 and 3. As for 3, this model is dependent
	on the SGS viscocity and can be used both for laminar and turbulent Convection/Diffusion
	Concentration problems.



