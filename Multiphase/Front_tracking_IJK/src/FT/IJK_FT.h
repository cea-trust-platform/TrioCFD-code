/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : IJK_FT.h
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_FT_included
#define IJK_FT_included

#include <FixedVector.h>
#include <IJK_Field.h>
#include <IJK_Splitting.h>
#include <OpDiffTurbIJK.h>
#include <OpConvIJKQuickSharp.h>
#include <OpCentre4IJK.h>
#include <OpConvIJKAmont.h>
#include <Multigrille_Adrien.h>
#include <Interprete.h>
#include <Linear_algebra_tools.h>
#include <Boundary_Conditions.h>
#include <IJK_Interfaces.h>
#include <Ref_Probleme_base.h>
#include <Redistribute_Field.h>
#include <Parser.h>
#include <IJK_FT_Post.h>
#include <IJK_Thermique.h>
#include <init_forcage_THI.h>
#include <corrections_qdm.h>
#include <Force_sp.h>
#include <IJK_Energie.h>

Declare_liste(IJK_Thermique);
Declare_liste(IJK_Energie);

// #define SMOOTHING_RHO

/*! @brief : class IJK_FT
 *
 *  <Description of class IJK_FT>
 *
 *
 *  La classe IJK_FT_double herite de la classe Interprete.
 *
 */
class IJK_FT_double : public Interprete
{
  // We take too much advantage of it ...:
  friend class IJK_Thermique;
  friend class Statistiques_dns_ijk_FT;
  Declare_instanciable(IJK_FT_double) ;

public :
  enum TimeScheme { EULER_EXPLICITE, RK3_FT };
  Entree& interpreter(Entree&) override;

  // Methodes d'acces :
  TimeScheme get_time_scheme() const;
  const IJK_Field_double& get_IJK_field(const Nom& nom) const;
  const FixedVector<IJK_Field_double, 3>& get_velocity() const
  {
    return velocity_;
  };
  const FixedVector<IJK_Field_double, 3>& get_velocity_ft() const
  {
    return velocity_ft_;
  };
  const double& get_current_time() const
  {
    return current_time_ ;
  }
  const double& get_timestep() const
  {
    return timestep_ ;
  }
  const int& get_splitting_extension() const
  {
    return ijk_splitting_ft_extension_ ;
  }
  const IJK_Splitting& get_splitting_ft() const
  {
    return splitting_ft_ ;
  }
  const IJK_Splitting& get_splitting_ns() const
  {
    return splitting_ ;
  }
  const Probleme_base& probleme(const IJK_Splitting& split) const
  {
    if (split == splitting_ft_)
      return refprobleme_ft_disc_.valeur();
    if (split == splitting_)
      //  return refprobleme_ns_.valeur();
      Cerr << "Unrecognized splitting provided" << finl;
    Process::exit();
    return refprobleme_ns_.valeur();
  }
  double t_debut_statistiques() const
  {
    return post_.t_debut_statistiques();
  }
  int get_reprise() const
  {
    return reprise_;
  }
  double get_rho_l() const
  {
    return rho_liquide_;
  }
  double get_rho_v() const
  {
    return rho_vapeur_;
  }
//  const int& nb_thermal_equations() const
//  {
//      return thermique_.
//    static int nb=1;
//    return nb;
//  }

  void run();
  void euler_time_step(ArrOfDouble& var_volume_par_bulle);
  // MODIF Aymeric : c'est plus pratique de mettre cette methode publique,
  // c'est une methode generique cont.
  void euler_explicit_update(const IJK_Field_double& dv, IJK_Field_double& v,
                             const int k_layer) const;

  //MODIF GAB 26 nov 2020
  void write_check_etapes_et_termes(const int rk_step);
  //GAB
  double calculer_true_moyenne_de_phase_vap(const IJK_Field_double& vx);
  double calculer_moyenne_de_phase_vap(const IJK_Field_double& vx);
  double calculer_true_moyenne_de_phase_liq(const IJK_Field_double& vx);
  double calculer_moyenne_de_phase_liq(const IJK_Field_double& vx);
  //void compute_and_add_source_qdm_gr(const double threshold_0, const double threshold_1, const double threshold_2, const double threshold_3);
  void set_time_for_corrections();
  void compute_and_add_qdm_corrections();
  void compute_and_add_qdm_corrections_monophasic();
  IJK_Field_double scalar_product(const FixedVector<IJK_Field_double, 3>& V1, const FixedVector<IJK_Field_double, 3>& V2);
  FixedVector<IJK_Field_double, 3> scalar_times_vector(const IJK_Field_double& Sca, const FixedVector<IJK_Field_double, 3>& Vec);
  IJK_Field_double scalar_fields_product(const IJK_Field_double& S1, const IJK_Field_double& S2, int dir);
  // SURCHARGE DES OPERATEURS : dans FixedVector, on ajoute le produit_scalaire
  //  mais il faut que les operateur * et += soient definis pour IJK_FT_double
  //  /!\ operator* est une fonction, pas une methode de la classe
  //      operator*= est une methode de la classe

  // IJK_FT_double EST UNE COLLECTION DE CHAMPS : CA N'A AUCUN SENS DE SOMMER CETTE COLLECTION

  //IJK_FT_double operator+=(const IJK_FT_double& Champ1, const IJK_FT_double& Champ2)
  //{
  //  int ni(Champ1.ni()),nj(Champ1.nj()),nk(Champ1.nk());
  //  int i,j,k;
  //  for (i=0; i<ni; i++)
  //    {
  //      for (j=0; j<nj; j++)
  //        {
  //          for (k=0; k<nk; k++)
  //            {
  //              resu(i,j,k)=Champ1(i,j,k)+Champ2(i,j,k)
  //  IJK_FT_double resu;
  //  resu =
  //}
  //IJK_FT_double operator*=(const IJK_FT_double& )
  // const IJK_Field_double& get_indicatrice_ft() const {return interfaces_.indicatrice_ft(); }
  // const IJK_Field_double& get_indicatrice_ft_next() const {return interfaces_.indicatrice_ft_next(); }
  const IJK_Interfaces& itfce() const {return interfaces_;}
  int disable_diphasique() const {return disable_diphasique_;}
  // Redistribute_Field& redistrib_to_ft_elem() const {return redistribute_to_splitting_ft_elem_;}
  // FixedVector<Redistribute_Field, 3>& redistrib_to_ft_face() const {return redistribute_to_splitting_ft_faces_;}

  // TODO: aym attention, cette methode n'est pas constante
  Redistribute_Field& redistrib_from_ft_elem() {return redistribute_from_splitting_ft_elem_;}
  // FixedVector<Redistribute_Field, 3>& redistrib_from_ft_face() const {return redistribute_from_splitting_ft_faces_;}
  void get_redistribute_from_splitting_ft_faces(
    const FixedVector<IJK_Field_double, 3>& faces_ft,
    FixedVector<IJK_Field_double, 3>& faces_ns
  )
  {
    for (int dir = 0; dir < 3; dir++)
      {
        redistribute_from_splitting_ft_faces_[dir].redistribute(
          faces_ft[dir],
          faces_ns[dir]);
      }
  }

  const IJK_Grid_Geometry& get_geometry() const
  {
    return splitting_.get_grid_geometry();
  }
  const IJK_FT_Post& get_post() const {return post_;}

protected :
  // Interdit constructeur par copie et operateur copie
  IJK_FT_double(const IJK_FT_double& x);

  const IJK_FT_double& operator=(const IJK_FT_double&);
  void update_rho_v();
  int initialise();
  void terme_source_gravite(IJK_Field_double& dv, int k_index, int dir) const;
  void rk3_sub_step(const int rk_step, const double total_timestep,
                    const double fractionnal_timestep, const double time);
  void calculer_dv(const double timestep, const double time, const int rk_step);

  //ab-sauv/repr-deb
  void ecrire_donnees(const FixedVector<IJK_Field_double, 3>& f3compo, SFichier& le_fichier, const int compo, bool binary) const;
  void dumpxyz_vector(const FixedVector<IJK_Field_double, 3>& f3compo, const char * filename, bool binary) const;
  void sauvegarder_probleme(const char *fichier_sauvegarde); //  const;
  void reprendre_probleme(const char *fichier_reprise);
  //ab-sauv/repr-fin

  // These methods are static in order to make it clear that all data used and modified by the method are explicitely passed as arguments:
  static void force_entry_velocity(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz, double v_imposed);
  static void force_upstream_velocity(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                      double v_imposed,const IJK_Interfaces& interfaces, double nb_diam );
  double find_timestep(const double max_timestep, const double cfl, const double fo, const double oh);
  void parcourir_maillage();
  void calculer_rho_mu_indicatrice(const bool parcourir = true);
  void maj_indicatrice_rho_mu(const bool parcourir = true);
  void ajout_dTrustine();
  void deplacer_interfaces(const double timestep, const int rk_step, ArrOfDouble& var_volume_par_bulle);
  void deplacer_interfaces_rk3(const double timestep, const int rk_step, ArrOfDouble& var_volume_par_bulle);
  void calculer_gradient_indicatrice_et_repul_ns(const IJK_Field_double& indic);

  //ab-forcage-control-ecoulement-deb
  void calculer_terme_source_acceleration(IJK_Field_double& vx, const double time, const double timestep,
                                          const int rk_step);
  void compute_correction_for_momentum_balance(const int rk_step);
  void transfer_ft_to_ns();
  void compute_add_external_forces(const int dir);
  // GAB
  void compute_add_THI_force(const FixedVector<IJK_Field_double, 3>& vitesse,
                             const int time_iteration,
                             const double dt,
                             const double current_time,
                             const IJK_Splitting& my_splitting
                            );
  void compute_add_THI_force_sur_d_velocity(const FixedVector<IJK_Field_double, 3>& vitesse,
                                            const int time_iteration,
                                            const double dt,
                                            const double current_time,
                                            const IJK_Splitting& my_splitting,
                                            const int facteur
                                           );

  void fill_variable_source_and_potential_phi(const double time);

  // GAB, rotation
  int get_direction(const ArrOfDouble& vecteur);
  //
  // GAB, qdm
  Vecteur3 calculer_inv_rho_grad_p_moyen(const IJK_Field_double& inv_rho, const IJK_Field_double& pression);
  Vecteur3 calculer_grad_p_moyen(const IJK_Field_double& pression);
  Vecteur3 calculer_grad_p_over_rho_moyen(const IJK_Field_double& pression);
  FixedVector<IJK_Field_double, 3> terme_convection_mass_solver_;
  FixedVector<IJK_Field_double, 3> terme_diffusion_mass_solver_;
  FixedVector<IJK_Field_double, 3> rho_u_euler_av_prediction_champ;
  FixedVector<IJK_Field_double, 3> rho_du_euler_ap_prediction_champ;
  FixedVector<IJK_Field_double, 3> rho_u_euler_ap_projection_champ;
  FixedVector<IJK_Field_double, 3> rho_du_euler_ap_projection_champ;
  FixedVector<IJK_Field_double, 3> rho_u_euler_av_rho_mu_ind_champ;
  FixedVector<IJK_Field_double, 3> rho_u_euler_ap_rho_mu_ind_champ;
  FixedVector<IJK_Field_double, 3> terme_diffusion_local;
  FixedVector<IJK_Field_double, 3> terme_pression_local;
  FixedVector<IJK_Field_double, 3> terme_pression_in_ustar_local;
  FixedVector<IJK_Field_double, 3> d_v_diff_et_conv;
  Vecteur3 rho_u_euler_av_prediction;
  Vecteur3 rho_du_euler_ap_prediction;
  Vecteur3 rho_u_euler_ap_projection;
  Vecteur3 rho_du_euler_ap_projection;
  Vecteur3 rho_u_euler_av_rho_mu_ind;
  Vecteur3 rho_u_euler_ap_rho_mu_ind;
  Vecteur3 u_euler_ap_rho_mu_ind;
  Vecteur3 terme_diffusion;
  Vecteur3 terme_convection;
  Vecteur3 terme_pression;
  Vecteur3 terme_pression_bis;
  Vecteur3 terme_pression_ter;
  Vecteur3 terme_interfaces;
  Vecteur3 terme_pression_in_ustar;
  Vecteur3 terme_interfaces_bf_mass_solver;
  Vecteur3 terme_interfaces_bf_mass_solver_bis;
  Vecteur3 terme_interfaces_af_mass_solver;
  Vecteur3 temre_intf_conv_diff_mass_solver;
  Vecteur3 terme_moyen_convection_mass_solver_;
  Vecteur3 terme_moyen_diffusion_mass_solver_;
  double pression_ap_proj;
  //
  // GAB qdm patch a posteriori
  //TODO :  enum corrections_qdm::type_dict_ { GB, GR };
  //int patch_qdm_gr_;  // flag
  //ArrOfDouble qdm_patch_correction_ = 0.; // correction value
  //ArrOfDouble qdm_patch_correction_deux_;
  // GAB source qdm a posteriori
  /*
  int source_qdm_gr_;
  ArrOfDouble liste_instants_;
  */
  //ArrOfDouble liste_vap_dl_;  // liste des v_v * dt
  //ArrOfDouble liste_liq_dl_;  // liste des v_l * dt
  /*
  int size_listes_source_;
  bool slide_ = false;
  double time_interval_;
  double v_terminale_;
  double source_qdm_update_dt_ = 3e-2;
  double last_source_qdm_update_time_ = 0.;
  */
  double vap_velocity_tmoy_ = 0.;
  double reprise_vap_velocity_tmoy_ = 0.;
  double liq_velocity_tmoy_ = 0.;
  double reprise_liq_velocity_tmoy_ = 0.;
  /*
  double qdm_source_ = 0.;
  double v_target_ = 0.;
  int list_index_ = 0.;
  int offset_list_index_ = 0.;
  double reprise_v_target_ = 0.;
  double reprise_qdm_source_ = 0.;
  */
  //
  Nom expression_derivee_acceleration_;
  Parser parser_derivee_acceleration_;
  Noms expression_variable_source_; // on attend trois expressions
  Nom expression_potential_phi_; // source variable formulee en gradient
  Vecteur3 store_rhov_moy_;
  Vecteur3 integrated_residu_;
  // terme source qdm pour pousser le fluide dans le canal (en m/s/s)
  double terme_source_acceleration_;

  // Vecteurs de taille 3 a lire dans le jeu de donnees :
  ArrOfDouble terme_source_correction_; // Valeur de la force de correction moyenne a appliquer
  ArrOfInt correction_force_; // 3 flags d'activation de la correction ou non

  // Storage for rhov for the evaluation of the acceleration source term :
  FixedVector<IJK_Field_double, 3> rho_v_;
  FixedVector<IJK_Field_double, 3> psi_velocity_; // for storage.
  //ab-forcage-control-ecoulement-fin

  FixedVector<IJK_Field_double, 3> variable_source_;
  Nom expression_derivee_facteur_variable_source_;
  Parser parser_derivee_facteur_variable_source_;
  double facteur_variable_source_; // ArrOfDouble? vecteur de taille 3 a lire dans le jeu de donnees

  IJK_Field_double potential_phi_;


  int check_divergence_;
  int rk_step_;

  Nom check_stop_file_; // Nom du fichier stop

  //ab-sauv/repr-deb

  Nom fichier_post_; // Nom du fichier post
  int dt_sauvegarde_;
  int sauvegarder_xyz_; // drapeau 0 ou 1
  Nom nom_sauvegarde_;
  Nom nom_reprise_;
  int reprise_;// flag pour indiquer si on fait une reprise
  // Le jeu de donnees doit fournir soit des fichiers de reprise: ..
  Nom fichier_reprise_vitesse_;
  int timestep_reprise_vitesse_;
  // ... soit des expressions f(x,y,z)
  Noms expression_vitesse_initiale_; // on attend trois expressions
  Nom expression_pression_initiale_; // useless, unless post-pro OR pressure_increment.
  //ab-sauv/repr-fin

  IJK_Splitting splitting_;
  // Le probleme ft disc qui porte le maillage vdf pour les algorithmes front-tracking
  REF(Probleme_base) refprobleme_ft_disc_;
  // Creation d'un probleme sur le domaine d'origine pour les sondes et pour faire leur VDF...
  REF(Probleme_base) refprobleme_ns_;

  ArrOfDouble_with_ghost delta_z_local_;

  Boundary_Conditions boundary_conditions_;

  // Inconnues du probleme (a sauvegarder et a reprendre)
  // Velocity field:
  FixedVector<IJK_Field_double, 3> velocity_;

  // Masse volumique:
  IJK_Field_double rho_field_;
  IJK_Field_double inv_rho_field_;

  // Pour les cas a bulles fixes
  double coef_immobilisation_;
  double coef_ammortissement_;
  double coef_mean_force_;
  double coef_force_time_n_;
  double coef_rayon_force_rappel_;
  double p_seuil_max_;
  double p_seuil_min_;
  FixedVector<IJK_Field_double, 3> force_rappel_;
  FixedVector<IJK_Field_double, 3> force_rappel_ft_;

  double vol_bulle_monodisperse_; // Pour imposer le volume des bulles
  ArrOfDouble vol_bulles_;   // Le volume impose individuellement a chaque bulle.

  // Field only needed for the option type_velocity_convection_form_== Nom("non_conservative_rhou")
  IJK_Field_double div_rhou_;

  // Temporary storage for the derivative
  FixedVector<IJK_Field_double, 3> d_velocity_;
  // Temporary storage for the fractional timestep in rk3 :
  FixedVector<IJK_Field_double, 3> RK3_F_velocity_;

  // To work in pressure increment and potentially help the solver for high density ratios
  IJK_Field_double d_pressure_;
  // Only if pressure increment formulation AND RK3 time_scheme.
  IJK_Field_double RK3_F_pressure_;

  // Celui la est discretise sur le maillage etendu:
  FixedVector<IJK_Field_double, 3> terme_source_interfaces_ft_;
  FixedVector<IJK_Field_double, 3> terme_repulsion_interfaces_ft_;
  FixedVector<IJK_Field_double, 3> terme_abs_repulsion_interfaces_ft_;

  // Celui la est discretise sur le maillage ns:
  FixedVector<IJK_Field_double, 3> terme_source_interfaces_ns_; // [ dt (rho u )/dt = N/m^3 ]
  FixedVector<IJK_Field_double, 3> terme_repulsion_interfaces_ns_;
  FixedVector<IJK_Field_double, 3> terme_abs_repulsion_interfaces_ns_;

  // Pour l'option alternative du calcul de la diffusion :
  IJK_Field_double unit_;
  FixedVector<IJK_Field_double, 3> laplacien_velocity_;

#ifdef BIDOUILLE
  // Test pour corriger le bilan de qdm :
  IJK_Field_double rho_old_;
  FixedVector<IJK_Field_double, 3> velocity_old_;
  FixedVector<IJK_Field_double, 3> psi_velocity_;
#endif

  // Booleen pour savoir si on a mis a jour l'indicatrice avec indicatrice next
  // bool indicatrice_is_indicatrice_next_;

  // Pressure field
  IJK_Field_double pressure_;
  // Molecular diffusivity (see diffusion operator)
  IJK_Field_double molecular_mu_;
  // right hand side for pressure solver
  IJK_Field_double pressure_rhs_;
  // Operators and pressure solver
  OpDiffIJK_double velocity_diffusion_op_simple_;
  OpDiffStdWithLaminarTransposeIJK_double velocity_diffusion_op_full_;
  Nom type_velocity_diffusion_form_; /* simple_arithmetic : div(mu grad(u))
  	  	  	  	  	  	  	  	  	 * full_arithmetic    : Tenseur des contraintes complet : div[mu (grad(u)+grad^T(u))]
  	  	  	  	  	  	  	  	  	 *                      mu : moyenne arithmetique
  	  	  	  	  	  	  	  	  	 * full_adaptative    : Tenseur des contraintes complet : div[mu (grad(u)+grad^T(u))]
  	  	  	  	  	  	  	  	  	 *     mu : switch from arithmetic to geometric mean depending on the direction (Not available yet)
  	  	  	  	  	  	  	  	  	 */
  OpConvIJKQuickSharp_double velocity_convection_op_sharp_;
  OpConvCentre4IJK_double velocity_convection_op_centre_;
  OpConvAmontIJK_double velocity_convection_op_amont_;

  Nom type_velocity_convection_form_; /* non_conservative_simple : rho div(u u)
  	  	  	  	  	  	  	  	  	 * non_conservative_rhou     : div(rho u u) - u div(rho u)
  	  	  	  	  	  	  	  	  	 * conservative              : div(rho u u)
  	  	  	  	  	  	  	  	  	 */
  //OpConvAmontIJK_double velocity_convection_op_;
  int type_velocity_convection_op_; // 0 : Quick  / 1 : Centre / 2 : Amont

  Multigrille_Adrien poisson_solver_;
  // Simulation parameters
  int nb_timesteps_;
  double timestep_;
  double timestep_facsec_;
  double cfl_, fo_, oh_;

  double vitesse_entree_;
  double vitesse_upstream_;
  double nb_diam_upstream_;
  double rho_liquide_;
  double rho_vapeur_;
  // GAB, pour THI (23.08.21)
  double rho_moyen_;
  //
  double mu_liquide_;
  double mu_vapeur_;
  double sigma_;
  ArrOfDouble gravite_; // vecteur de taille 3 a lire dans le jeu de donnees
  int direction_gravite_;
#ifdef SMOOTHING_RHO
  int smooth_density_;
  double ratio_density_max_;
  // Redistribute_Field redistribute_to_splitting_ft_elem_; // pour le smoothing de rho..
  IJK_Field_double rho_field_ft_;
#endif
  Redistribute_Field redistribute_to_splitting_ft_elem_;
  int disable_solveur_poisson_;
  int disable_diffusion_qdm_;
  int disable_convection_qdm_;
  int disable_source_interf_;
  int frozen_velocity_;
  int velocity_reset_;

  // Pour la premiere projection, on initialise la pression au champ de pression a l'equilibre diphasique :
  int improved_initial_pressure_guess_;
  // travail en increment de pression pour aider le solveur :
  int include_pressure_gradient_in_ustar_;
  // Discretisation du champ (1/rho) et utilisation dans le calcul de rho*v*v, dans le mass_solver
  //    et dans pressure_projection_with_inv_rho :
  int use_inv_rho_for_mass_solver_and_calculer_rho_v_;
  int use_inv_rho_in_poisson_solver_;
  int use_inv_rho_;

  int correction_bilan_qdm_;

  int refuse_patch_conservation_QdM_RK3_source_interf_;
  // GAB, qdm
  int test_etapes_et_bilan;
  //
  // GAB, champ de reprise + champ initial
  int add_initial_field_;
  //
  int diffusion_alternative_;
  int suppression_rejetons_;
  // Supprime l'appel a quelques fonctions qui n'ont pas de sens en monophasique :
  // comme par exemple : deplacer_interfaces,
  //    calculer_rho_mu_indicatrice, ecrire_statistiques_bulles
  int disable_diphasique_;
  int time_scheme_;
  double store_RK3_source_acc_;
  double store_RK3_fac_sv_;
  double current_time_;
  double current_time_at_rk3_step_;
  int tstep_; // The iteration number

  // GAB
  init_forcage_THI forcage_;
  corrections_qdm qdm_corrections_;
  // Le maillage des interfaces:
  IJK_Interfaces interfaces_;

  // Maillage etendu pour les interfaces
  // Nombre de mailles etendues dans les directions periodiques
  int ijk_splitting_ft_extension_;
  IJK_Splitting splitting_ft_;
  // Classe outil pour passer entre splitting_ et splitting_ft_
  // Une instance par direction des faces:
  FixedVector<Redistribute_Field, 3> redistribute_to_splitting_ft_faces_;
  FixedVector<Redistribute_Field, 3> redistribute_from_splitting_ft_faces_;
  Redistribute_Field redistribute_from_splitting_ft_elem_;

  // Champs volumiques sur le maillage FT (pour convection des interfaces, etc)
  FixedVector<IJK_Field_double, 3> velocity_ft_;

  // Dealing with thermal aspects:
  LIST(IJK_Thermique) thermique_;
  LIST(IJK_Energie) energie_;

  //
  // Post-processing helper class:
  //
  friend class IJK_FT_Post;
  IJK_FT_Post post_;
};

#endif /* IJK_FT_included */
