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
// File      : IJK_FT_Post.h
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////


#ifndef IJK_FT_Post_included
#define IJK_FT_Post_included

#include <IJK_Field.h>
#include <Motcle.h>
#include <TRUST_Vector.h>
#include <Statistiques_dns_ijk_FT.h>
#include <Statistiques_dns_ijk_monophasique.h>
#include <Sondes_IJK.h>
#include <IJK_Interfaces.h>
#include <Multigrille_Adrien.h>


class IJK_FT;
class IJK_Splitting;
class IJK_Thermique;
class IJK_Energie;
class List_IJK_Thermique;
class list_curseurIJK_Thermique;
class const_list_curseurIJK_Thermique;
class List_IJK_Energie;
class list_curseurIJK_Energie;
class const_list_curseurIJK_Energie;
/**
 * All the post-processing stuff of IJK_FT delegated into this helper class:
 */
class IJK_FT_Post
{

  friend class Statistiques_dns_ijk_FT;

public:
  IJK_FT_Post(IJK_FT_double& ijk_ft);
  void complete_interpreter(Param& param, Entree& e);
  int initialise(int reprise);
  void complete(int reprise);
  int initialise_stats(IJK_Splitting& splitting, ArrOfDouble& vol_bulles, const double vol_bulle_monodisperse);
  void init_indicatrice_non_perturbe();

  void posttraiter_champs_instantanes(const char * lata_name, double time, int time_iteration);
  void posttraiter_statistiques_plans(double time);
  void ecrire_statistiques_bulles(int reset, const Nom& nom_cas, const ArrOfDouble& gravite, const double current_time) const;
  void update_stat_ft(const double dt);
  void get_update_lambda2();
  void get_update_lambda2_and_rot_and_curl();

  IJK_Field_double& rebuilt_indic()
  {
    return rebuilt_indic_;
  }
  IJK_Field_double& potentiel()
  {
    return potentiel_;
  }
  FixedVector<IJK_Field_double, 3>& coords()
  {
    return coords_;
  }
  IJK_Field_double& integrated_timescale()
  {
    return integrated_timescale_;
  }
  bool postraiter_sous_pas_de_temps() const
  {
    return postraiter_sous_pas_de_temps_;
  }
  int dt_post() const
  {
    return dt_post_;
  }
  int post_par_paires() const
  {
    return post_par_paires_;
  }
  double t_debut_statistiques() const
  {
    return t_debut_statistiques_;
  }

  inline int sondes_demande()
  {
    return sondes_demande_;
  }
  const IJK_Field_double& get_IJK_field(const Nom& nom) const;
  const int& get_IJK_flag(const Nom& nom) const;

  const FixedVector<IJK_Field_double, 3>& get_IJK_vector_field(const Nom& nom) const;

  inline FixedVector<IJK_Field_double, 3>& get_grad_I_ns()
  {
    return grad_I_ns_;
  }

  void sauvegarder_post(const Nom& lata_name);
  void sauvegarder_post_maitre(const Nom& lata_name, SFichier& fichier) const;
  void reprendre_post(Param& param);

  void fill_op_conv();
  void fill_surface_force(FixedVector<IJK_Field_double, 3>& the_field_you_know);//const Nom lata_name, double instant, int iteration);
  void fill_surface_force_bis(const char * lata_name, double time, int time_iteration);
  FixedVector<IJK_Field_double, 3> get_rho_Ssigma();

  void calculer_gradient_indicatrice_et_pression(const IJK_Field_double& indic);

  // Part of the run() method in IJK_FT:
  int alloc_fields();
  int alloc_velocity_and_co(bool flag_variable_source);
  void completer_sondes();
  void postraiter_sondes();
  void improved_initial_pressure_guess(bool imp);
  void postraiter_ci(const Nom& lata_name, const double current_time);
  void postraiter_fin(bool stop, int tstep, double current_time, double timestep, const Nom& lata_name,
                      const ArrOfDouble& gravite, const Nom& nom_cas);
  //void ijk_interpolate_implementation_bis(const IJK_Field_double& field, const DoubleTab& coordinates, ArrOfDouble& result,
  //                                        int skip_unknown_points, double value_for_bad_points,const IJK_Field_double& indic);
  //void  ijk_interpolate_skip_unknown_points_bis(const IJK_Field_double& field, const DoubleTab& coordinates, ArrOfDouble& result,
  //                                              const double value_for_bad_points,const IJK_Field_double& indic);
  void compute_extended_pressures(const Maillage_FT_IJK& mesh);
  //IJK_Field_double& extended_p);
  void posttraiter_tous_champs_thermique(Motcles& liste,  const int idx) const;
  void posttraiter_tous_champs_energie(Motcles& liste,  const int idx) const;
  int posttraiter_champs_instantanes_thermique(const Motcles& liste_post_instantanes,
                                               const char *lata_name,
                                               const int lata_step, const double current_time,
                                               IJK_Thermique& ,  const int idx);
  int posttraiter_champs_instantanes_energie(const Motcles& liste_post_instantanes,
                                             const char *lata_name,
                                             const int lata_step, const double current_time,
                                             IJK_Energie& ,  const int idx);
  int posttraiter_champs_instantanes_thermique_interfaciaux(const Motcles& liste_post_instantanes,
                                                            const char *lata_name,
                                                            const int lata_step, const double current_time,
                                                            IJK_Thermique& ,  const int idx);
  int posttraiter_champs_instantanes_energie_interfaciaux(const Motcles& liste_post_instantanes,
                                                          const char *lata_name,
                                                          const int lata_step, const double current_time,
                                                          IJK_Energie& ,  const int idx);
//  void calculer_gradient_temperature(const IJK_Field_double& temperature, FixedVector<IJK_Field_double, 3>& grad_T);

  Motcles get_liste_post_instantanes() const
  {
    return liste_post_instantanes_;
  }
protected:
  void compute_phase_pressures_based_on_poisson(const int phase);
  Statistiques_dns_ijk_FT statistiques_FT_;
  int dt_post_;
  int dt_post_stats_plans_; // intervalle de posttraitement des donnees par plan (pour les statistiques de canal)
  int dt_post_stats_bulles_; // intervalle de posttraitement des donnees par bulles
  Motcles liste_post_instantanes_; // liste des champs instantanes a postraiter
  // Pour numeroter les fichiers .lata il faut compter combien on en a ecrit:
  int compteur_post_instantanes_;
  int postraiter_sous_pas_de_temps_; // drapeau 0 ou 1
  // Pour reconstruire au post-traitement la grandeur du/dt, on peut choisir de relever u^{dt_post} et u^{dt_post+1} :
  int post_par_paires_; // drapeau 0 ou 1

  // Pour des fiches de validation, on post-traite le champ analytique attendu dans le lata pour calcul de l'erreur:
  Noms expression_vitesse_analytique_; // on attend trois expressions
  Nom expression_pression_analytique_; // on attend une expression
  Noms expression_dvitesse_analytique_; // on attend trois expressions

  // Pour check_stats (and_grads)
  int check_stats_;
  Noms expression_gradP_analytique_; // on attend trois expressions
  Noms expression_gradU_analytique_; // on attend trois expressions
  Noms expression_gradV_analytique_; // on attend trois expressions
  Noms expression_gradW_analytique_; // on attend trois expressions
  // Second gradient (laplacian_P or Hessienne
  Noms expression_grad2P_analytique_; // on attend 6 expressions (car tens sym)
  // And for the 3 components of velocity :
  Noms expression_grad2U_analytique_; // on attend 6 expressions (car tens sym)
  Noms expression_grad2V_analytique_; // on attend 6 expressions (car tens sym)
  Noms expression_grad2W_analytique_; // on attend 6 expressions (car tens sym)

  // -------------------------------------------------
  // Statistiques temporelles
  // Pour les groupes de bulles : on cree un vecteur d'objet stat de taille dimensionnee a 0 puis nb_groups():
  VECT(Statistiques_dns_ijk_FT) groups_statistiques_FT_;
  //Statistiques_dns_ijk_FT statistiques_FT_;
  //
  // 2020.03.12. CHOIX : Meme en disable_diphasique, on fait appel a la classe fille stats FT.
  // La classe de stats monophasique n'est plus maintenue. Suppression du membre.
  // Statistiques_dns_ijk_monophasique statistiques_;
  double t_debut_statistiques_;
  // -------------------------------------------------

  // Pour les cas a bulles fixes
  FixedVector<IJK_Field_double, 3> integrated_velocity_;
  IJK_Field_double integrated_pressure_;
  IJK_Field_double indicatrice_non_perturbe_;
  IJK_Field_double integrated_timescale_;
  // Pour la reprise bulles fixes, parametres de lecture de champ de condition initiale pour variables de post-tt:
  Nom fichier_reprise_integrated_velocity_;
  Nom fichier_reprise_integrated_pressure_;
  Nom fichier_reprise_indicatrice_non_perturbe_;
  Nom fichier_reprise_integrated_timescale_;

  // Temporary storage for the coords (for postprocessing) :
  FixedVector<IJK_Field_double, 3> coords_;
  FixedVector<IJK_Field_double, 3> velocity_ana_;
  FixedVector<IJK_Field_double, 3> ecart_ana_;
  FixedVector<IJK_Field_double, 3> op_conv_;
  FixedVector<IJK_Field_double, 3> cell_op_conv_;
  FixedVector<IJK_Field_double, 3> rho_Ssigma_;
  FixedVector<IJK_Field_double, 3> cell_rho_Ssigma_;

  FixedVector<IJK_Field_double, 3> d_velocity_ana_;
  IJK_Field_double pressure_ana_,ecart_p_ana_;

  // Celui la est discretise sur le maillage etendu:
  FixedVector<IJK_Field_double, 3> grad_I_ft_;

  // Pour postraitement :
  IJK_Field_double rebuilt_indic_;
  IJK_Field_double potentiel_;
  IJK_Field_double ai_ft_;
  int extended_pressure_computed_;
  IJK_Field_double pressure_ft_;
  IJK_Field_double extended_pl_ft_;
  IJK_Field_double extended_pv_ft_;
  IJK_Field_double extended_pl_;
  IJK_Field_double extended_pv_;
  //For the liquid pressure gradient
  // FixedVector<IJK_Field_double 3> dP_ft_;
  // FixedVector<IJK_Field_double 3> dP_;
  // Pour le calcul des stats  :
  IJK_Field_double kappa_ai_ft_;
  FixedVector<IJK_Field_double, 3> normale_cell_ft_;
  IJK_Field_double ai_ns_;
  IJK_Field_double kappa_ai_ns_;
  FixedVector<IJK_Field_double, 3> normale_cell_ns_;
  // The following three fields are needed too for the gradient extension
// /IJK_Field_double dudy_;
  //IJK_Field_double dvdx_;//
  //IJK_Field_double dwdy_;
  // For lambda and curl
  IJK_Field_double dudx_;
  IJK_Field_double dvdy_;
  IJK_Field_double dwdx_;
  IJK_Field_double dudz_;
  IJK_Field_double dvdz_;
  IJK_Field_double dwdz_;
  IJK_Field_double critere_Q_;
  FixedVector<IJK_Field_double, 3> rot_;
  FixedVector<IJK_Field_double, 3> grad_I_ns_;
  FixedVector<IJK_Field_double, 3> grad_P_;
  //  FixedVector<IJK_Field_double, 3> grad_U_ns_;
  //  FixedVector<IJK_Field_double, 3> grad_V_ns_;
  //  FixedVector<IJK_Field_double, 3> grad_W_ns_;
  IJK_Field_double num_compo_ft_;

  // Pour la verification des stats :
  // Le gradient de pression aux faces :
  //  FixedVector<IJK_Field_double, 3> gradP_;
  // Les gradients des compo de vitesses aux elems : (sont finalement stockes dans statistiques_FT_ si besoin)
  //FixedVector<IJK_Field_double, 3> dUd_;
  //FixedVector<IJK_Field_double, 3> dVd_;
  //FixedVector<IJK_Field_double, 3> dWd_;
  // Et leurs solutions analytiques :
  FixedVector<IJK_Field_double, 3> ana_gradP_;
  FixedVector<IJK_Field_double, 3> ana_dUd_;
  FixedVector<IJK_Field_double, 3> ana_dVd_;
  FixedVector<IJK_Field_double, 3> ana_dWd_;
  // Pour les deriv secondes :
  FixedVector<IJK_Field_double, 3> ana_grad2Pi_; // Partie diagonale de la jacobienne
  FixedVector<IJK_Field_double, 3> ana_grad2Pc_; // contient les deriv croisees
  FixedVector<IJK_Field_double, 3> ana_grad2Ui_; // Partie diagonale de la jacobienne
  FixedVector<IJK_Field_double, 3> ana_grad2Uc_; // contient les deriv croisees
  FixedVector<IJK_Field_double, 3> ana_grad2Vi_; // Partie diagonale de la jacobienne
  FixedVector<IJK_Field_double, 3> ana_grad2Vc_; // contient les deriv croisees
  FixedVector<IJK_Field_double, 3> ana_grad2Wi_; // Partie diagonale de la jacobienne
  FixedVector<IJK_Field_double, 3> ana_grad2Wc_; // contient les deriv croisees

  // GAB
  IJK_Field_double IFd_source_spectraleX_;
  IJK_Field_double AOD_source_spectraleX_;
  IJK_Field_double source_spectraleY_;
  IJK_Field_double source_spectraleZ_;
  // Pour post-traitement :
  IJK_Field_double lambda2_, dudy_, dvdx_, dwdy_;
  FixedVector<IJK_Field_double, 3> cell_velocity_;
  FixedVector<IJK_Field_double, 3> cell_source_spectrale_;
  FixedVector<IJK_Field_double, 3> cell_bk_tsi_ns_;
  //  FixedVector<IJK_Field_double, 3> cell_source_interface_totale_;   // non-const because some echange_espace_virtuel()
  FixedVector<IJK_Field_double, 3> cell_grad_p_;
  FixedVector<IJK_Field_double, 3> cell_source_interface_; // toujours en _ns_
  FixedVector<IJK_Field_double, 3> cell_backup_source_interface_; // toujours en _ns_
  FixedVector<IJK_Field_double, 3> cell_repulsion_interface_; // toujours en _ns_


  int sondes_demande_;
  Sondes_IJK les_sondes_;           // Sondes a traiter

  //
  // References to various members of IJK_FT_double that are heavily used in the post:
  //
  IJK_FT_double& ref_ijk_ft_;

  const int& disable_diphasique_;    // yes a ref, not a const value.
  const IJK_Interfaces& interfaces_;
  IJK_Field_double& kappa_ft_;
  IJK_Field_double& pressure_;                   // non-const because some echange_espace_virtuel()
  FixedVector<IJK_Field_double, 3>& velocity_;   // non-const because some echange_espace_virtuel()
  FixedVector<IJK_Field_double, 3>& source_spectrale_;   // non-const because some echange_espace_virtuel()
  FixedVector<IJK_Field_double, 3>& bk_tsi_ns_;
  FixedVector<IJK_Field_double, 3> source_interface_ft_;   // non-const because some echange_espace_virtuel()
  FixedVector<IJK_Field_double, 3> source_interface_ns_;   // non-const because some echange_espace_virtuel()
  FixedVector<IJK_Field_double, 3> repulsion_interface_ns_;   // non-const because some echange_espace_virtuel()
  const FixedVector<IJK_Field_double, 3>& d_velocity_;

  IJK_Splitting& splitting_;
  IJK_Splitting& splitting_ft_;
  LIST(IJK_Thermique)& thermique_;
  LIST(IJK_Energie)& energie_;
  /* IJK_Field_double temperature_ana_, ecart_t_ana_;
    Nom expression_T_ana_;

    IJK_Field_double source_temperature_ana_, ecart_source_t_ana_; */

// FixedVector<IJK_Field_double, 3> grad_T_;


  Multigrille_Adrien poisson_solver_post_;
};


#endif /* IJK_FT_Post_included */
