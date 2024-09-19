/****************************************************************************
* Copyright (c) 2024, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
//////////////////////////////////////////////////////////////////////////////
//
// File:        Navier_Stokes_FT_Disc_interne.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/22
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Navier_Stokes_FT_Disc_interne_included
#define Navier_Stokes_FT_Disc_interne_included

class Navier_Stokes_FT_Disc_interne
{
public:
  Navier_Stokes_FT_Disc_interne() :
    correction_courbure_ordre_(0), // Par defaut, pas de prise en compte de la courbure pour corriger le champ etendu delta_vitesse
    mpoint_inactif(0),   // Par defaut, mpoint cree un saut de vitesse
    mpointv_inactif(0),   // Par defaut, mpointv  cree un saut de vitesse
    matrice_pression_invariante(0),   // Par defaut, recalculer la matrice pression
    clipping_courbure_interface(1e40),   // Par defaut, pas de clipping
    terme_gravite_(GRAVITE_GRAD_I),   // Par defaut terme gravite ft sans courants parasites
    is_explicite(1),                  // Par defaut, calcul explicite de vpoint etape predicition
    is_boussinesq_(0),                // Par defaut, l'hypothese de Boussinesq n'est pas utilisee pour la flottabilite dans les phases.
    new_mass_source_(0),              // Par defaut, on utilise la methode historique pour imposer le saut de vitesse du changement de phase.
    type_interpol_indic_pour_dI_dt_(INTERP_STANDARD), // Default is the historical interpolation
    OutletCorrection_pour_dI_dt_(NO_CORRECTION),   // Default is the historical
    is_penalized(0),                  // Par defaut, pas de penalisation L2 du forcage
    eta(1.0),                         // Par defaut, coefficient de penalisation L2 = 1.
    p_ref_pena(-1.e40),               // Par defaut, pas de penalisation L2 de la pression sinon valeur reference
    is_pfl_flottant(0),               // Traitement local Dirichlet pression si les CL pression sont toutes en Neumann homogene
    x_pfl_imp(-1.e40),                // Par defaut, x, y, z du point de modification de la pression fluide
    y_pfl_imp(-1.e40), z_pfl_imp(-1.e40)
  { }

  int correction_courbure_ordre_;
  int mpoint_inactif;
  int mpointv_inactif;

  Champ_Fonc second_membre_projection;
  Champ_Fonc second_membre_projection_jump_;
  Champ_Fonc derivee_u_etoile;
  Champ_Fonc gradient_pression;
  Champ_Fonc terme_diffusion;
  Champ_Fonc terme_convection;
  Champ_Fonc terme_source;
  Champ_Fonc terme_source_interfaces;
  Champ_Fonc indicatrice_p1b;
  Champ_Fonc gradient_indicatrice;
  Champ_Fonc potentiel_faces;
  Champ_Fonc potentiel_elements;
  // delta_u_interface = la partie "saut de vitesse" du champ de vitesse a l'interface
  OWN_PTR(Champ_Inc_base) delta_u_interface;
  Champ_Fonc laplacien_d;
  Champ_Fonc mpoint;
  Champ_Fonc mpoint_vap;
  // Variation temporelle indicatrice de phase
  Champ_Fonc derivee_temporelle_indicatrice;
  Champ_Fonc ai; // Eulerian interfacial area.
  OWN_PTR(Champ_Inc_base) vitesse_jump0_; // Extended Velocity of phase 0.

  LIST(REF(Champ_base)) liste_champs_compris;

  // Si matrice_pression_invariante != 0,
  //   on ne recalcule pas la matrice de pression a chaque pas de temps.
  int matrice_pression_invariante;
  // Si on veut ajouter une interface a vitesse imposee :
  //  reference a l'equation de transport correspondante :
  VECT(REF(Transport_Interfaces_FT_Disc)) ref_eq_interf_vitesse_imposee;
  // Si le fluide est diphasique, c'est l'indicatrice de l'equation suivante
  // qui est utilisee pour determiner les proprietes du fluide:
  // (masse volumique, viscosite, tension superficielle, ...)
  REF(Transport_Interfaces_FT_Disc) ref_eq_interf_proprietes_fluide;
  // Si le fluide est diphasique, la reference au fluide:
  REF(Fluide_Diphasique) ref_fluide_diphasique;

  REF(Convection_Diffusion_Temperature_FT_Disc) ref_equation_mpoint_;
  REF(Convection_Diffusion_Temperature_FT_Disc) ref_equation_mpoint_vap_;

  // Valeur maximale de courbure autorisee pour calculer le
  // terme source de tension de surface (clipping si valeur superieur)
  double clipping_courbure_interface;

  enum Terme_Gravite
  {
    GRAVITE_RHO_G, GRAVITE_GRAD_I
  };
  Terme_Gravite terme_gravite_;
  Noms equations_concentration_source_fluide_;
  // Si is_explicite != 0,
  //   on calcul vpoint de facon explicite dans l etape de prediction des vitesses.
  int is_explicite;
  // Si is_boussinesq_ != 0, on calcul une force par le modele de boussinesq
  int is_boussinesq_;
  // Flag pour la method du saut de vitesse :
  int new_mass_source_;

  enum Type_interpol_indic_pour_dI_dt
  {
    INTERP_STANDARD, INTERP_MODIFIEE, INTERP_AI_BASED, INTERP_STANDARD_UVEXT, INTERP_MODIFIEE_UVEXT, INTERP_AI_BASED_UVEXT, INTERP_STANDARD_UIEXT, INTERP_MODIFIEE_UIEXT, INTERP_AI_BASED_UIEXT
  };
  Type_interpol_indic_pour_dI_dt type_interpol_indic_pour_dI_dt_;

  enum OutletCorrection_pour_dI_dt
  {
    NO_CORRECTION, CORRECTION_GHOST_INDIC, ZERO_NET_FLUX_ON_MIXED_CELLS, ZERO_OUT_FLUX_ON_MIXED_CELLS
  };
  OutletCorrection_pour_dI_dt OutletCorrection_pour_dI_dt_;

  // Si is_penalized != 0,
  //   on penalise L2 le terme de forcage.
  int is_penalized;
  // Valeur de l'inverse du coefficient de penalisation L2 du terme de forcage.
  double eta;
  // Valeur pour la penalisation L2 de la pression.
  double p_ref_pena;
  // Point de penalisation L2 de la pression du fluide
  int is_pfl_flottant; // Traitement local Dirichlet pression si les CL pression sont toutes en Neumann homogene
  double x_pfl_imp;
  double y_pfl_imp;
  double z_pfl_imp;
};

#endif /* Navier_Stokes_FT_Disc_interne_included */
