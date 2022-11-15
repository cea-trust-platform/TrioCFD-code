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

#include <Connectivite_frontieres.h>
#include <FixedVector.h>
#include <IJK_Field.h> // est-ce que j'en ai vraiment besoin ?
#include <Linear_algebra_tools_impl.h>
#include <Maillage_FT_IJK.h>
#include <Objet_U.h>
#include <Parcours_interface.h>
#include <Ref_IJK_FT_double.h>
#include <Ref_Zone_dis.h>
#include <Remaillage_FT_IJK.h>
#include <SFichier.h>
#include <Vecteur3.h>
#include <Linear_algebra_tools_impl.h>
#include <SurfaceVapeurIJKComputation.h>
#include <ComputeValParCompoInCell.h>

#define VERIF_INDIC 0

/*! @brief : class IJK_Interfaces
 *
 *  Cette classe rassemble tous les algorithmes de gestion des interfaces pour le ijk
 *  (le maillage, les algo de remaillage, sauvegarde, reprise, etc)
 *
 *
 *
 */
class IJK_Interfaces : public Objet_U
{

  Declare_instanciable_sans_constructeur(IJK_Interfaces);

public :
  IJK_Interfaces();
  void initialize(const IJK_Splitting& splitting_FT,
                  const IJK_Splitting& splitting_NS,
                  const Zone_dis& zone_dis,
                  const bool compute_vint=true);
  void associer(const IJK_FT_double& ijk_ft);
  void posttraiter_tous_champs(Motcles& liste) const;
  int posttraiter_champs_instantanes(const Motcles& liste_post_instantanes,
                                     const char *lata_name,
                                     const int lata_step) const;
  void sauvegarder_interfaces(const char *lata_name); // const;
  void calculer_color(ArrOfInt& color) const;
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
  void supprimer_certaines_bulles_reelles();

  void deplacer_bulle_perio(const ArrOfInt& masque_deplacement_par_compo);
  void transferer_bulle_perio(); // Les bulles trop proche du bord sont
  // deplacees a l'oppose, vers l'autre bord perio.

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
                                        const IJK_Field_double& indic,
                                        const ArrOfDouble& volume_reel,
                                        const DoubleTab& position) const;

  int lire_motcle_non_standard(const Motcle& un_mot, Entree& is) override;
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
    frozen_ = 1;
  };

  int get_nb_bulles_ghost(const int print = 0) const
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

  void set_reprise(const int i)
  {
    reprise_ = i;
    return;
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
    calculer_volume_bulles(vol, positions_reference_);
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

  // Methode qui parcourt les facettes contenues dans la cellule num_elem
  // pour trouver la phase de sa cellule voisine
  // Donnees d'entree:
  //   - num_elem   : indice de l'element traverse par l'interface
  //   - direction  : precise dans quelle direction se trouve le voisin dont on
  //   cherche la phase
  //   - face_plus  : precise dans quel sens se trouve le voisin dont on cherche
  //   la phase
  //                  (+1 pour le voisin d'indice plus eleve, -1 pour l'autre )
  int compute_cell_phase_with_interface_normal(int num_elem, int direction, int face_plus);

  void calculer_normales_et_aires_interfaciales(IJK_Field_double& ai,
                                                IJK_Field_double& kappa_ai,
                                                FixedVector<IJK_Field_double, 3>& normale_cell,
                                                const int igroup) const;

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
  void calculer_poussee_bulles(const ArrOfDouble& gravite, DoubleTab& poussee) const;
  void calculer_aire_interfaciale(IJK_Field_double& ai) const;
  void calculer_normale_et_aire_interfaciale(IJK_Field_double& ai,
                                             IJK_Field_double& kappa_ai,
                                             FixedVector<IJK_Field_double, 3>& normale_cell) const;

  void compute_drapeaux_vapeur_v2(const IntVect& vecteur_composantes,
                                  ArrOfInt& drapeau_liquide) const;
  void compute_drapeaux_vapeur_v3(const Maillage_FT_IJK& mesh,
                                  const IJK_Splitting& split,
                                  const IntVect& vecteur_composantes,
                                  ArrOfInt& drapeau_vapeur) const;
  void compute_drapeaux_vapeur_v4(const IntVect& vecteur_composantes,
                                  ArrOfInt& drapeau_vapeur) const;

  void convert_to_IntVect(const ArrOfInt& in, IntVect& out) const;

  void ajouter_terme_source_interfaces(
    FixedVector<IJK_Field_double, 3>& vpoint,
    FixedVector<IJK_Field_double, 3>& vrepul,
    FixedVector<IJK_Field_double, 3>& vabsrepul
  ) const;

  void remailler_interface(const double temps,
                           Maillage_FT_IJK& maillage,
                           ArrOfDouble& var_volume,
                           Remaillage_FT_IJK& algo_remaillage_local);

  int is_terme_gravite_rhog() const;
  void detecter_et_supprimer_rejeton(bool duplicatas_etaient_presents);
  void update_surface_normale() const;

  // Permet de recuperer un mcu qui est une vue de maillage_bulles_ft_ijk. Il
  // n'y a pas de copie memoire, seulement des passages de case memoire.
  // splitting.get_grid_geometry().get_origin, get_constant_delta(DIRECTION_X),
  // ...
  static void get_maillage_MED_from_IJK_FT(MEDCouplingUMesh *maillage_bulles_mcu,
                                           const Maillage_FT_IJK& maillage_bulles_ft_ijk);

  // Getter des surfaces par face
  const FixedVector<IJK_Field_double, 3>& get_surface_vapeur_par_face_ft() const
  {
    return surface_vapeur_par_face_[old()];
  }
  const FixedVector<IJK_Field_double, 3>& get_surface_vapeur_par_face() const
  {
    return surface_vapeur_par_face_ns_[old()];
  }
  // Getter des surfaces par face
  // void get_surface_vapeur_par_face_ns(FixedVector<IJK_Field_double, 3> &surfs) const ;
  // Getter des barycentres par face
  const FixedVector<FixedVector<IJK_Field_double, 3>, 3>& get_barycentre_vapeur_par_face_ft() const
  {
    return barycentre_vapeur_par_face_[old()];
  }
  const FixedVector<FixedVector<IJK_Field_double, 3>, 3>& get_barycentre_vapeur_par_face() const
  {
    return barycentre_vapeur_par_face_ns_[old()];
  }

  int get_nb_face_mouillees() const { return n_faces_mouilles_[old()]; }

  // TODO: utiliser le allocate de allocate_velociter dans
  // IJK_Navier_Stockes_tools.cpp utiliser le pslitting de NS, pas le FT

  //////////////////////////////////////////////////////////////////////
  // API indicatrice et variables du maillage par cellule et par face //
  //////////////////////////////////////////////////////////////////////

  static double mean_over_compo(
    const FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& field_for_compo,
    const IJK_Field_int& nb_compo_traversante,
    const int i,
    const int j,
    const int k
  )
  {
    double res = 0.;
    for (int compo = 0; compo <= nb_compo_traversante(i, j, k); compo++)
      {
        res += field_for_compo[compo](i, j, k);
      }
    if (nb_compo_traversante(i,j,k) > 0)
      return res / nb_compo_traversante(i, j, k);
    else
      return 0.;
  }

  static void mean_over_compo(
    const FixedVector<IJK_Field_double, 3 * max_authorized_nb_of_components_>& field_for_compo,
    const IJK_Field_int& nb_compo_traversante,
    FixedVector<IJK_Field_double, 3>& mean_par_compo_field
  )
  {
    const int ni = nb_compo_traversante.ni();
    const int nj = nb_compo_traversante.nj();
    const int nk = nb_compo_traversante.nk();
    for (int dir=0; dir < 3; dir ++)
      for (int i=0; i < ni; i++)
        for (int j=0; j < nj; j++)
          for (int k=0; k < nk; k++)
            {
              double res = 0.;
              for (int compo = 0; compo <= nb_compo_traversante(i, j, k); compo++)
                {
                  int idx = 3 * compo + dir;
                  const double last_val = field_for_compo[idx](i, j, k);
                  res += last_val;
                }
              if (nb_compo_traversante(i,j,k) > 0)
                res = res / nb_compo_traversante(i, j, k);
              else
                res = 0.;
              mean_par_compo_field[dir](i,j,k) = res;
            }
  }

  static double mean_over_compo(
    const FixedVector<IJK_Field_double, 3 * max_authorized_nb_of_components_>& field_for_compo,
    const IJK_Field_int& nb_compo_traversante,
    const int dir,
    const int i,
    const int j,
    const int k
  )
  {
    double res = 0.;
    for (int compo = 0; compo <= nb_compo_traversante(i, j, k); compo++)
      {
        int idx = 3 * compo + dir;
        const double last_val = field_for_compo[idx](i, j, k);
        res += last_val;
      }
    if (nb_compo_traversante(i,j,k) > 0)
      return res / nb_compo_traversante(i, j, k);
    else
      return 0.;
  }

  int old() const { return 1 - old_en_premier_; }
  int next() const { return old_en_premier_; }

  // TODO: retirer l'acces a FT
  const FixedVector<IJK_Field_double, max_authorized_nb_of_components_ * 3>& get_bary_par_compo_itfc_in_cell_ft() const
  {
    return bary_par_compo_[old()];
  }
  // TODO: retirer l'acces a FT
  const FixedVector<IJK_Field_double, max_authorized_nb_of_components_ * 3>& get_norm_par_compo_itfc_in_cell_ft() const
  {
    return normale_par_compo_[old()];
  }

  // TODO: retirer l'acces a FT
  const IJK_Field_double& I_ft() const { return indicatrice_ft_[old()]; }
  const double& I_ft(const int i, const int j, const int k) const { return indicatrice_ft_[old()](i, j, k); }
  const IJK_Field_double& In_ft() const { return indicatrice_ft_[next()]; }

  const IJK_Field_double& I() const { return indicatrice_ns_[old()]; }
  const IJK_Field_double& In() const { return indicatrice_ns_[next()]; }
  const double& I(const int i, const int j, const int k) const { return indicatrice_ns_[old()](i, j, k); }
  const double& In(const int i, const int j, const int k) const { return indicatrice_ns_[next()](i, j, k); }

  const double& SI(const int compo, const int i, const int j, const int k) const
  {
    return surface_par_compo_[old()][compo](i, j, k);
  }
  double SI(const int i, const int j, const int k) const
  {
    double res = mean_over_compo(surface_par_compo_[old()], nb_compo_traversante_[old()], i, j, k);
    return res;
  }

  const FixedVector<IJK_Field_double, max_authorized_nb_of_groups_>& groups_indicatrice_ft() const
  {
    return groups_indicatrice_ft_[old()];
  }
  const FixedVector<IJK_Field_double, max_authorized_nb_of_groups_>& groups_indicatrice_ns() const
  {
    return groups_indicatrice_ns_[old()];
  }
  const FixedVector<IJK_Field_double, max_authorized_nb_of_groups_>& groups_indicatrice_n_ft() const
  {
    return groups_indicatrice_ft_[next()];
  }
  const FixedVector<IJK_Field_double, max_authorized_nb_of_groups_>& groups_indicatrice_n_ns() const
  {
    return groups_indicatrice_ns_[next()];
  }
  // const double& SIn(const int compo, const int i, const int j, const int k)
  // const {
  //   return surface_par_compo_next_[compo](i,j,k);
  // }
  // double SIn(const int i, const int j, const int k) const {
  //   double res = mean_over_compo(surface_par_compo_next_, i,j,k);
  //   return res;
  // }

  // Les composantes (dir) du vecteur normal moyen pour chaque compo comnnexe
  // (c) de bulle dans le maille (compo = max_nb_of_compo_ * c + dir)
  const double& nI(const int compo, const int i, const int j, const int k) const
  {
    return normal_of_interf_ns_[old()][compo](i, j, k);
  }
  Vecteur3 nI(const int i, const int j, const int k) const
  {
    const Vecteur3 res(normal_of_interf_ns_[old()][0](i, j, k),
                       normal_of_interf_ns_[old()][1](i, j, k),
                       normal_of_interf_ns_[old()][2](i, j, k));
    return res;
  }
  const double& nIn(const int compo, const int i, const int j, const int k) const
  {
    return normal_of_interf_ns_[next()][compo](i, j, k);
  }
  Vecteur3 nIn(const int i, const int j, const int k) const
  {
    const Vecteur3 res(normal_of_interf_ns_[next()][0](i, j, k),
                       normal_of_interf_ns_[next()][1](i, j, k),
                       normal_of_interf_ns_[next()][2](i, j, k));
    return res;
  }

  // Les composantes (dir) du barycentre de l'interface pour chaque compo
  // comnnexe (c) de bulle dans le maille (compo = max_nb_of_compo_ * c + dir)
  const double& xI(const int compo, const int i, const int j, const int k) const
  {
    return bary_of_interf_ns_[old()][compo](i, j, k);
  }
  Vecteur3 xIn(const int i, const int j, const int k) const
  {
    const Vecteur3 res(bary_of_interf_ns_[next()][0](i, j, k),
                       bary_of_interf_ns_[next()][1](i, j, k),
                       bary_of_interf_ns_[next()][2](i, j, k));
    return res;
  }
  Vecteur3 xI(const int i, const int j, const int k) const
  {
    const Vecteur3 res(bary_of_interf_ns_[old()][0](i, j, k),
                       bary_of_interf_ns_[old()][1](i, j, k),
                       bary_of_interf_ns_[old()][2](i, j, k));
    return res;
  }

  // Surface de la vapeur sur la face dans la direction compo de la cellule ijk
  const double& Sf(const int compo, const int i, const int j, const int k) const
  {
    return surface_vapeur_par_face_ns_[old()][compo](i, j, k);
  }
  // Surface de la vapeur sur la face dans la direction compo de la cellule ijk
  // au temps n+1
  const double& Sfn(const int compo, const int i, const int j, const int k) const
  {
    return surface_vapeur_par_face_ns_[next()][compo](i, j, k);
  }
  // Surface de la vapeur sur la face dans la direction compo de la cellule ijk
  // moyennee entre n et n+1.
  // TODO: faire un tableau surface_vapeur_par_face_moy et le calculer en
  // faisant la moyenne sur plusieurs petits pas.
  double Sfm(const int compo, const int i, const int j, const int k) const
  {
    return (surface_vapeur_par_face_ns_[old()][compo](i, j, k) + surface_vapeur_par_face_ns_[next()][compo](i, j, k)) * 0.5;
  }

  void switch_indicatrice_next_old();
  void calculer_indicatrice_next(
    IJK_Field_double& field_repulsion,
    const ArrOfDouble& gravite,
    const double delta_rho,
    const double sigma,
    const double time,
    const int itstep,
    const bool parcourir = true
  );
  void set_compute_surfaces_mouillees() { surface_vapeur_par_face_computation_.set_compute_surfaces_mouillees(); }

  const int& nb_compo_traversantes(const int i, const int j, const int k) const
  {
    return nb_compo_traversante_[next()](i,j,k);
  }

protected:
  // Met à jour les valeurs de surface_vapeur_par_face_ et barycentre_vapeur_par_face_
  SurfaceVapeurIJKComputation surface_vapeur_par_face_computation_;

  ComputeValParCompoInCell val_par_compo_in_cell_computation_;

  void verif_indic() ;

  void calculer_phi_repuls_sommet(
    ArrOfDouble& potentiels_sommets,
    ArrOfDouble& repulsions_sommets,
    const ArrOfDouble& gravite,
    const double delta_rho,
    const double sigma,
    const double time,
    const int itstep
  );

  void calculer_phi_repuls_par_compo(
    FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& phi_par_compo,
    FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& repuls_par_compo,
    IJK_Field_double& field_repulsion,
    const ArrOfDouble& gravite,
    const double delta_rho,
    const double sigma,
    const double time,
    const int itstep
  );

  void calculer_indicatrice(IJK_Field_double& indic);
  void calculer_indicatrice_optim(IJK_Field_double& indic);
  void calculer_indicatrices(FixedVector<IJK_Field_double, 3>& indic);
  void calculer_indicatrices_optim(FixedVector<IJK_Field_double, 3>& indic);

  // Methode qui parcourt tous les elements de indic et met a jour uniquement
  // ceux qui etaient traverses par l'interface a l'iteration precedente et qui
  // ne le sont plus a l'iteration courante ces elements ont ete marques au
  // prealable
  int update_indicatrice(IJK_Field_double& indic);


  // Cette methode appelle la methode statique get_maillage_MED_from_IJK_FT sur
  // ses propres membres. Elle met donc a jour le maillage maillage_bulles_med_.

  // TODO: utiliser le allocate de allocate_velocity dans IJK_Navier_Stokes_tools.cpp utiliser le pslitting de NS, pas le FT

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
  // Interdit le constructeur par copie (car constructeurs par copie interdits
  // pour parcours_ et autres
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
  // Le parcours a besoin d'une connectivite des faces de bords du maillage
  // eulerien
  Connectivite_frontieres connectivite_frontieres_;
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

  // Distance max en metres a laquelle agit la force de repulsion entre les
  // bulles
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
  // Terme_Gravite terme_gravite_;
  int terme_gravite_;

  int nb_groups_;           // Nombre de groupes/classes de bulles.
  ArrOfInt compo_to_group_; // Tableau de conversion: numero_de_groupe =
  // compo_to_group_[icompo]

  /////////////////////////////////////////////////////////////////
  // Tableau des valeurs aux faces mouillees (calcules avec med) //
  /////////////////////////////////////////////////////////////////

  bool compute_surf_mouillees_; // active seulement dans le cas
  // ou il y a des champs thermique ou d energie.
  // attention, ca desactive seulement le calcul, pas l'allocation.
  FixedVector<int,2> n_faces_mouilles_;

  bool is_diphasique_;

  // Surfaces vapeur des faces du maillage IJK
  FixedVector<FixedVector<IJK_Field_double, 3>, 2> surface_vapeur_par_face_;
  FixedVector<FixedVector<IJK_Field_double, 3>, 2> surface_vapeur_par_face_ns_;

  // Normale de l'interface par maille ijk sur domaine NS
  FixedVector<FixedVector<IJK_Field_double, 3>, 2> normal_of_interf_;
  FixedVector<FixedVector<IJK_Field_double, 3>, 2> normal_of_interf_ns_;
  // Barycentre de l'interface par maille ijk sur domaine NS
  FixedVector<FixedVector<IJK_Field_double, 3>, 2> bary_of_interf_;
  FixedVector<FixedVector<IJK_Field_double, 3>, 2> bary_of_interf_ns_;

  // Barycentres vapeur des faces du maillage IJK
  FixedVector<FixedVector<FixedVector<IJK_Field_double, 3>, 3>, 2> barycentre_vapeur_par_face_;
  FixedVector<FixedVector<FixedVector<IJK_Field_double, 3>, 3>, 2> barycentre_vapeur_par_face_ns_;

  /////////////////////////////////////
  // indicatrice et var moy par cell //
  /////////////////////////////////////

  int n_cell_diph_;
  bool old_en_premier_;

  FixedVector<IJK_Field_double, 2> indicatrice_ns_;
  FixedVector<IJK_Field_double, 2> indicatrice_ft_;

  // On prevoie un tableau assez grand pour contenir tous les groupes.
  FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_groups_>, 2> groups_indicatrice_ft_;
  FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_groups_>, 2> groups_indicatrice_ns_;

#if VERIF_INDIC
  // pour verifier le calcul optimise de l'indicatrice
  IJK_Field_double indicatrice_ft_test_;
  FixedVector<IJK_Field_double, max_authorized_nb_of_groups_> groups_indicatrice_ft_test_;
#endif

  // Vecteur des composantes normale dans chaque cellule par composante connexe
  FixedVector<IJK_Field_int, 2> nb_compo_traversante_;
  FixedVector<FixedVector<IJK_Field_int, max_authorized_nb_of_components_>, 2> compos_traversantes_;
  FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_components_>, 2> indicatrice_par_compo_;
  FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_components_>, 2> surface_par_compo_;
  FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_components_>, 2> courbure_par_compo_;
  FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_components_>, 2> phi_par_compo_;
  FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_components_>, 2> repuls_par_compo_;
  FixedVector<FixedVector<IJK_Field_double, 3 * max_authorized_nb_of_components_>, 2> normale_par_compo_;
  FixedVector<FixedVector<IJK_Field_double, 3 * max_authorized_nb_of_components_>, 2> bary_par_compo_;
};

#endif /* IJK_Interfaces_included */
