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
// File      : IJK_Interfaces.h
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Interfaces_included
#define IJK_Interfaces_included

#include <Objet_U.h>
#include <Maillage_FT_IJK.h>
#include <FixedVector.h>
#include <Parcours_interface.h>
#include <Connectivite_frontieres.h>
#include <Remaillage_FT_IJK.h>
#include <Ref_Zone_dis.h>
#include <SFichier.h>
#include <Vecteur3.h>
#include <Linear_algebra_tools_impl.h>
#include <Ref_IJK_FT_double.h>

// #define SMOOTHING_RHO
static const int max_authorized_nb_of_groups_ = 3;

/////////////////////////////////////////////////////////////////////////////
// .NAME        : IJK_Interfaces
// .DESCRIPTION : class IJK_Interfaces
//
// Cette classe rassemble tous les algorithmes de gestion des interfaces pour le ijk
// (le maillage, les algo de remaillage, sauvegarde, reprise, etc)
//
/////////////////////////////////////////////////////////////////////////////

class IJK_Interfaces : public Objet_U
{

  Declare_instanciable(IJK_Interfaces) ;

public :
  void initialize(const IJK_Splitting& splitting_FT,
                  const IJK_Splitting& splitting_NS,
                  const Zone_dis& zone_dis);
  void associer(const IJK_FT_double& ijk_ft);
  void posttraiter_tous_champs(Motcles& liste) const;
  int posttraiter_champs_instantanes(const Motcles& liste_post_instantanes,
                                     const char *lata_name,
                                     const int lata_step) const;
  void sauvegarder_interfaces(const char *lata_name); // const;
  void calculer_color(ArrOfInt& color) const ;
  void postraiter_colors(Sortie& os, const double current_time) const;

  void transporter_maillage(const double dt_tot,
                            ArrOfDouble& dvol,
                            const int rk_step,
                            const double temps);
  void calculer_bounding_box_bulles(DoubleTab& bounding_box) const;
  void preparer_duplicata_bulles(const DoubleTab& bounding_box_of_bubbles,
                                 const DoubleTab& authorized_bounding_box,
                                 ArrOfInt& masque_duplicata_pour_compo);
  void dupliquer_bulle_perio(ArrOfInt& masque_duplicata_pour_compo);
  void creer_duplicata_bulles();
  void supprimer_duplicata_bulles();
//void supprimer_bulle(ArrOfInt & liste_bulle_pour_suppression);
  void supprimer_certaines_bulles_reelles();

  void deplacer_bulle_perio(const ArrOfInt& masque_deplacement_par_compo);
  void transferer_bulle_perio(); // Les bulles trop proche du bord sont deplacees a l'oppose,
//                                vers l'autre bord perio.

  void compute_vinterp();
  // methode pour bulles fixes
  void compute_external_forces_(FixedVector<IJK_Field_double, 3>& rappel_ft,
                                FixedVector<IJK_Field_double, 3>& rappel,
                                const FixedVector<IJK_Field_double, 3>& vitesse,
                                const IJK_Field_double& indic_ns,
                                const IJK_Field_double& indic_ft,
                                const double coef_immo,
                                const int tstep,
                                const double current_time,
                                const double coef_ammortissement,
                                const double coef_rayon_force_rappel,
                                double compteur,
                                double coef_mean_force,
                                double coef_force_time_n);
  void compute_external_forces_parser(FixedVector<IJK_Field_double, 3>& rappel,
                                      const IJK_Field_double& indic_ns,
                                      const DoubleTab& individual_forces,
                                      const ArrOfDouble& volume_reel,
                                      const DoubleTab& position,
                                      const double coef_rayon_force_rappel);
  void compute_external_forces_color_function(FixedVector<IJK_Field_double, 3>& rappel,
                                              const IJK_Field_double& indic_ns,
                                              const IJK_Field_double& indic_ft,
                                              DoubleTab& individual_forces,
                                              const ArrOfDouble& volume_reel,
                                              const DoubleTab& position);

  void compute_indicatrice_non_perturbe(IJK_Field_double& indic_np,
                                        IJK_Field_double& indic,
                                        const ArrOfDouble& volume_reel,
                                        const DoubleTab& position) const;

  // fin de methode pour bulles fixes
  int get_nb_bulles_reelles() const
  {
    return nb_bulles_reelles_;
  };

  int get_flag_positions_reference() const
  {
    return flag_positions_reference_;
  };

  int is_frozen() const
  {
    return frozen_;
  };
  void freeze()
  {
    frozen_=1;
  };

  int get_nb_bulles_ghost(const int print=0) const
  {
    if (print)
      Cerr << "IJK_Interfaces::get_nb_bulles_ghost : On a " << nb_bulles_ghost_ << " bulles ghost." << finl;
    return nb_bulles_ghost_;
  };

  int get_forcing_method() const
  {
    return parser_;
  };
  int get_recompute_indicator() const
  {
    return recompute_indicator_;
  };

  void set_recompute_indicator(int i)
  {
    recompute_indicator_ = i;
  };

  int follow_colors() const
  {
    return follow_colors_;
  };

  const ArrOfInt& get_colors() const
  {
    return through_yminus_;
  };
  int ghost_compo_converter(const int i) const
  {
    return ghost_compo_converter_[i];
  };

// Methodes d'acces aux parametres de la repulsion a la paroi :
  double portee_wall_repulsion() const
  {
    return portee_wall_repulsion_;
  };
  double delta_p_wall_max_repulsion() const
  {
    return delta_p_wall_max_repulsion_;
  };

//
  void set_reprise(const int i)
  {
    reprise_ = i;
    return ;
  };

  void set_fichier_sauvegarde(const char *lataname)
  {
    fichier_sauvegarde_interface_ = lataname;
  };

  void set_fichier_reprise(const char *lataname)
  {
    fichier_reprise_interface_ = lataname;
  };
  Nom get_fichier_reprise()
  {
    return fichier_reprise_interface_;
  };
// methode non-const. Pour appel depuis l'exterieur (comme IJK_switch)
  void lire_maillage_ft_dans_lata()
  {
    maillage_ft_ijk_.lire_maillage_ft_dans_lata(fichier_reprise_interface_,
                                                timestep_reprise_interface_,
                                                lata_interfaces_meshname_);
    return;
  };

// methodes non-const.
  void parcourir_maillage()
  {
    maillage_ft_ijk_.parcourir_maillage();
    return;
  };

  void RK3_G_store_vi_resize(int n, int n2)
  {
    RK3_G_store_vi_.resize(n, n2);
    return;
  };

  void RK3_G_store_vi_echange_esp_vect()
  {
    maillage_ft_ijk_.desc_sommets().echange_espace_virtuel(RK3_G_store_vi_);
    return;
  };

// methode non-const. Remplit le tableau des positions de reference
// vers lesquelles la force de rappel attire les bulles):
// (a partir de la position actuelle des bulles)
// pour calculs a bulles fixes
  void set_positions_reference()
  {
    ArrOfDouble vol;
    calculer_volume_bulles(vol, positions_reference_ );
    flag_positions_reference_ = 1;
    return;
  };

  const Maillage_FT_IJK& maillage_ft_ijk() const
  {
    return maillage_ft_ijk_;
  }
  const DoubleTab& RK3_G_store_vi() const
  {
    return  RK3_G_store_vi_;
  }

  const IntVect& get_num_compo() const
  {
    return num_compo_;
  }
  void calculer_indicatrice(IJK_Field_double& indic);
  void calculer_indicatrice_optim(IJK_Field_double& indic);
  void calculer_indicatrices(FixedVector<IJK_Field_double, 3>&   indic);
  void calculer_indicatrices_optim(FixedVector<IJK_Field_double, 3>&   indic);

  // Methode qui parcourt tous les elements de indic et met a jour uniquement ceux
  // qui etaient traverses par l'interface a l'iteration precedente
  // et qui ne le sont plus a l'iteration courante
  // ces elements ont ete marques au prealable
  int update_indicatrice(IJK_Field_double& indic);

  // Methode qui parcourt les facettes contenues dans la cellule num_elem
  // pour trouver la phase de sa cellule voisine
  // Donnees d'entree:
  //   - num_elem   : indice de l'element traverse par l'interface
  //   - direction  : precise dans quelle direction se trouve le voisin dont on cherche la phase
  //   - face_plus  : precise dans quel sens se trouve le voisin dont on cherche la phase
  //                  (+1 pour le voisin d'indice plus eleve, -1 pour l'autre )
  int compute_cell_phase_with_interface_normal(int num_elem, int direction, int face_plus);

  void calculer_normales_et_aires_interfaciales(IJK_Field_double& ai,
                                                IJK_Field_double& kappa_ai,
                                                FixedVector<IJK_Field_double, 3>& normale_cell,
                                                const int igroup) const;
  void calculer_normale_et_bary_element_pour_compo(const int icompo,
                                                   const int elem,
                                                   const  Maillage_FT_IJK& maillage,
                                                   Vecteur3& normale,
                                                   Vecteur3& bary_facettes_dans_elem) const;
  int compute_list_compo_connex_in_element(const Maillage_FT_IJK& mesh,
                                           const int elem,
                                           ArrOfInt& liste_composantes_connexes_dans_element) const;
  const int& nb_groups() const
  {
    return nb_groups_ ;
  };
  const ArrOfDoubleFT& get_distance_autres_interfaces() const
  {
    return distance_autres_interfaces_;
  };
  int get_ghost_number_from_compo(const int compo) const;

  void calculer_surface_bulles(ArrOfDouble& surfaces) const;
  void compute_surface_average_per_bubble(const ArrOfDouble& surfaces, const ArrOfDouble& in, ArrOfDouble& out) const;
  void calculer_volume_bulles(ArrOfDouble& volumes, DoubleTab& centre_gravite) const;
  void calculer_poussee_bulles(const ArrOfDouble& gravite,  DoubleTab& poussee) const;
  void calculer_aire_interfaciale(IJK_Field_double& ai) const;
  void calculer_normale_et_aire_interfaciale(IJK_Field_double& ai,
                                             IJK_Field_double& kappa_ai,
                                             FixedVector<IJK_Field_double, 3>& normale_cell) const;

  void compute_drapeaux_vapeur_v2(const IntVect& vecteur_composantes,
                                  ArrOfInt& drapeau_liquide) const ;
  void compute_drapeaux_vapeur_v3(const Maillage_FT_IJK& mesh,
                                  const IJK_Splitting& split,
                                  const IntVect& vecteur_composantes,
                                  ArrOfInt& drapeau_vapeur) const;
  void compute_drapeaux_vapeur_v4(const IntVect& vecteur_composantes,
                                  ArrOfInt& drapeau_vapeur) const;

  void convert_to_IntVect(const ArrOfInt& in, IntVect& out) const;
#ifdef SMOOTHING_RHO
  void ajouter_terme_source_interfaces(FixedVector<IJK_Field_double, 3>& vpoint, double sigma,
                                       const ArrOfDouble& gravite,
                                       const double delta_rho,
                                       /*Pour post-traitement : */
                                       IJK_Field_double& field_indicatrice,
                                       IJK_Field_double& field_potentiel,
                                       /* Pour le smoothing : */
                                       IJK_Field_double& rho_field,
                                       const double rho_v,
                                       const int smooth_density,
                                       const double time,
                                       const int num_iter);
#else
  void ajouter_terme_source_interfaces(FixedVector<IJK_Field_double, 3>& vpoint,
                                       FixedVector<IJK_Field_double, 3>& vrepul,
                                       FixedVector<IJK_Field_double, 3>& vabsrepul,
                                       double sigma,
                                       const ArrOfDouble& gravite,
                                       const double delta_rho,
                                       /*Pour post-traitement : */
                                       IJK_Field_double& field_indicatrice,
                                       IJK_Field_double& field_potentiel,
                                       const double time,
                                       const int num_iter);
#endif


  void remailler_interface(const double temps,
                           Maillage_FT_IJK& maillage,
                           ArrOfDouble& var_volume,
                           Remaillage_FT_IJK& algo_remaillage_local);

  int is_terme_gravite_rhog() const;
  void detecter_et_supprimer_rejeton(bool duplicatas_etaient_presents) ;
  void update_surface_normale() const;

protected:
  void calculer_vmoy_composantes_connexes(const Maillage_FT_IJK& maillage,
                                          const ArrOfDouble& surface_facette,
                                          const ArrOfDouble& surface_par_bulle,
                                          const ArrOfInt& compo_connexes_facettes,
                                          const int nbulles_reelles,
                                          const int nbulles_ghost,
                                          const DoubleTab& vitesse_sommets,
                                          DoubleTab& vitesses) const;
  void recursive_calcul_distance_chez_voisin(DoubleTab& vinterp_tmp,
                                             int dir,
                                             const Maillage_FT_IJK& mesh,
                                             DoubleTab& coord_sommets,
                                             ArrOfInt& compo_sommet,
                                             ArrOfDouble& distance,
                                             DoubleTab& v_closer,
                                             double distmax);
  void calculer_distance_autres_compo_connexe2(ArrOfDouble& distance,
                                               DoubleTab& v_closer);
  void calculer_distance_autres_compo_connexe(const DoubleTab& sommets_a_tester,
                                              const ArrOfInt& compo_connexe_sommets,
                                              const DoubleTab& vinterp_tmp,
                                              const Maillage_FT_IJK& mesh,
                                              ArrOfDouble& distance,
                                              DoubleTab& v_closer,
                                              const double distmax);


// reference vers le splitting_ft_ pour les interfaces :
  REF(IJK_Splitting) ref_splitting_;
  REF(Zone_dis) refzone_dis_;
  REF(IJK_FT_double) ref_ijk_ft_;
// Interdit le constructeur par copie (car constructeurs par copie interdits pour parcours_ et autres
  IJK_Interfaces(const IJK_Interfaces& x) : Objet_U(x)
  {
    Cerr << "Erreur IJK_Interfaces(const IJK_Interfaces&)" << finl;
    Process::exit();
  }
  const IJK_Interfaces& operator=(const IJK_Interfaces&)
  {
    Cerr << "Erreur IJK_Interfaces& operator=" << finl;
    Process::exit();
    return *this;
  }
  // Ou lire le maillage initial, dans la methode initialize():
  // On peut reprendre un fichier lata ou sauv.lata :

  IntVect num_compo_;
  int nb_compo_in_num_compo_;
  // variales pour calcul a bulles fixes
  DoubleTab mean_force_;
  DoubleTab force_time_n_;
  DoubleTab positions_reference_;
  int flag_positions_reference_;

  Nom fichier_reprise_interface_;
  int timestep_reprise_interface_;
  Nom lata_interfaces_meshname_;

  // Pour ecrire dans le fichier sauv :
  Nom fichier_sauvegarde_interface_;
  int timestep_sauvegarde_interface_;

  // Activation du suivi des couleurs des bulles
  int follow_colors_;

  // Activer la repulsion aux parois :
  int active_repulsion_paroi_;

  // Modification de l'evaluation du potentiel :
  int correction_gradient_potentiel_;

  int compute_distance_autres_interfaces_;
  // Nombres de bulles reeles :
  int nb_bulles_reelles_;
  int nb_bulles_ghost_;
  int nb_bulles_ghost_before_;
  int recompute_indicator_;
  int parser_;
  // Stockage du maillage:
  Maillage_FT_IJK maillage_ft_ijk_;

  // Tableau intermediaire pour le deplacement des marqueurs en RK3 :
  DoubleTab RK3_G_store_vi_;
  DoubleTab vinterp_;

  // Algorithmes de parcours de l'interface (intersections Eulerien/Lagrangien)
  Parcours_interface parcours_;
  // Le parcours a besoin d'une connectivite des faces de bords du maillage eulerien
  Connectivite_frontieres  connectivite_frontieres_;
  // Algorithme de remaillage local
  Remaillage_FT_IJK remaillage_ft_ijk_;

  // Donnees sur la maillage NS : sa bounding_box et la periodicite :
  DoubleTab bounding_box_NS_domain_;
  bool perio_NS_[3];

  // Domaine autorise pour les bulles :
  // c'est le geom du splitting_FT reduit de ncells_forbidden_
  // dans toutes les direction ou le domaine NS est perio.
  int ncells_forbidden_;
  int ncells_deleted_;
  int frozen_; // flag to disable the interfaces motion.
  DoubleTab bounding_box_forbidden_criteria_;
  DoubleTab bounding_box_delete_criteria_;

  // Domaine a l'interieur duquel la duplication d'une bulle est inutile
  DoubleTab bounding_box_duplicate_criteria_;

  // Distance max en metres a laquelle agit la force de repulsion entre les bulles
  double portee_force_repulsion_;
  // delta de pression maxi cree par la force de repulsion
  // (pour l'instant lineaire, valeur max quand la distance est nulle)
  double delta_p_max_repulsion_;

  // Si souhaite, une valeur differente pour les parois :
  double portee_wall_repulsion_;
  double delta_p_wall_max_repulsion_;

  ArrOfDoubleFT distance_autres_interfaces_;

  // Pour chaque bulle, a-t-elle traverse la frontiere perio y?
  ArrOfInt through_yminus_; // 0 : On ne sait pas
  // 1 : entree par yminus
  // -1 : entree par yplus

  ArrOfInt ghost_compo_converter_;

  int reprise_; // Flag indiquant si on fait une reprise

  enum Terme_Gravite { GRAVITE_RHO_G, GRAVITE_GRAD_I };
  //Terme_Gravite terme_gravite_;
  int terme_gravite_;

  int nb_groups_; // Nombre de groupes/classes de bulles.
  ArrOfInt compo_to_group_; // Tableau de conversion: numero_de_groupe = compo_to_group_[icompo]
};

#endif /* IJK_Interfaces_included */
