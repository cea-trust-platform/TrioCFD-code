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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Transport_Interfaces_FT_Disc.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/38
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Transport_Interfaces_FT_Disc_included
#define Transport_Interfaces_FT_Disc_included

#include <Equation_base.h>
#include <Transport_Interfaces_base.h>
#include <Ref_Probleme_base.h>
#include <Postraitement_base.h>
#include <Champ_Inc.h>
#include <Remaillage_FT.h>
#include <Parcours_interface.h>
#include <Marching_Cubes.h>
#include <Connectivite_frontieres.h>
#include <Topologie_Maillage_FT.h>
#include <Algorithmes_Transport_FT_Disc.h>
#include <Champ_Fonc.h>
#include <Navier_Stokes_FT_Disc.h>
#include <Ref_Milieu_base.h>
#include <Ref_Loi_horaire.h>
#include <Proprietes_part_vol.h>
#include <Ref_Navier_Stokes_FT_Disc.h>
#include <TRUSTTabFT_forward.h>

class Transport_Interfaces_FT_Disc_interne;
template <typename titi> class TRUSTTab;
using FloatTab = TRUSTTab<float>;

class Transport_Interfaces_FT_Disc : public Transport_Interfaces_base
{
  Declare_instanciable_sans_constructeur(Transport_Interfaces_FT_Disc);
public:

  Transport_Interfaces_FT_Disc();
  //
  void set_param(Param& titi) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  // Methodes virtuelles pures de Equation_base
  //
  int            nombre_d_operateurs(void) const override; // Zero, y'en a pas.
  const Operateur& operateur(int i) const override;    // Erreur
  Operateur&        operateur(int i) override;         // Erreur
  const Champ_Inc& inconnue(void) const override;         // C'est l'indicatrice
  Champ_Inc&        inconnue(void) override;
  //
  // Methodes surchargees de Equation_base
  //
  void                associer_milieu_base(const Milieu_base& milieu) override;

  void                 associer_equation_ns(const Navier_Stokes_FT_Disc& ns);

  Milieu_base&        milieu() override;       // Erreur
  const Milieu_base& milieu() const override;  // Erreur
  void    associer_pb_base(const Probleme_base& probleme) override;
  void    discretiser(void) override;
  Entree& lire_cond_init(Entree& is) override;
  int  verif_Cl() const override;
  double  calculer_pas_de_temps(void) const override;
  DoubleTab& derivee_en_temps_inco(DoubleTab& derivee) override;
  void assembler( Matrice_Morse& mat_morse, const DoubleTab& present, DoubleTab& secmem) override ;

  void    mettre_a_jour(double temps) override;
  int  sauvegarder(Sortie& ) const override;
  int  reprendre(Entree&) override;
  int impr(Sortie& os) const override;
  void update_critere_statio();

  //
  // Nouvelles methodes
  //
  virtual void                            lire_maillage_ft_cao(Entree& is);
  int                          preparer_calcul() override;
  virtual void                            preparer_pas_de_temps();
  const Maillage_FT_Disc&                 maillage_interface() const;
  const Champ_base&               get_update_indicatrice() override;
  virtual const Champ_base&               get_indicatrice_faces();
  virtual const Champ_base&               get_compute_indicatrice_faces();
  virtual const Parcours_interface&       parcours_interface() const;
  virtual const Marching_Cubes&           marching_cubes() const;
  virtual const Algorithmes_Transport_FT_Disc& algorithmes_transport() const;
  virtual const Connectivite_frontieres& connectivite_frontieres() const;
  Remaillage_FT&                          remaillage_interface();
  const Remaillage_FT&                    remaillage_interface() const;
  const Topologie_Maillage_FT&            topologie_interface() const;
  virtual double calculer_integrale_indicatrice(const DoubleVect& indicatrice) const;

  const Proprietes_part_vol&           proprietes_particules() const;
  const Maillage_FT_Disc&              maillage_inject() const;
  const Proprietes_part_vol&           proprietes_inject() const;

  void nettoyer_proprietes_particules(const ArrOfInt& som_utilises);

  virtual void calculer_vitesse_transport_interpolee(const Champ_base& champ_vitesse,
                                                     const Maillage_FT_Disc&,
                                                     DoubleTab& vitesse_noeuds,
                                                     int nv_calc) const;
  void calculer_scalaire_interpole(const Champ_base& ch_scal,
                                   const Maillage_FT_Disc&,
                                   DoubleTab& ch_scal_noeuds,
                                   int nv_calc) const;

  virtual void remailler_interface();

  //methodes utilisees pour le post-traitement
  virtual int get_champ_post_FT(const Motcle& champ, Postraitement_base::Localisation loc, FloatTab *ftab = 0) const;
  virtual int get_champ_post_FT(const Motcle& champ, Postraitement_base::Localisation loc, IntTab   *itab = 0) const;
  virtual const Maillage_FT_Disc& maillage_interface_pour_post() const;
  int get_mesh_tag() const override
  {
    return maillage_interface_pour_post().get_mesh_tag();
  };

  //Methode d acces au probleme
  const Probleme_base& get_probleme_base() const;

  //Modifie vpoint (pour N_S) pour imposer au fluide la vitesse de l interface
  void modifier_vpoint_pour_imposer_vit(const DoubleTab& inco_val,DoubleTab& vpoint0,
                                        DoubleTab& vpoint,const DoubleTab& rho_faces,
                                        DoubleTab& terme_source,const double temps, const double dt,
                                        const int is_explicite,const double eta) override;

  //Methode outil utilisee par modifier_vpoint_pour_imposer_vit(...)
  void calcul_source(const DoubleTab& inco_val,
                     const DoubleTab& vpoint,
                     const DoubleTab& rho_faces,
                     DoubleTab& source_val,
                     const DoubleTab& vit_imposee,
                     const DoubleTab& indicatrice_faces,
                     const int is_QC,
                     const double dt,
                     const int is_explicite,
                     const double eta);
  void modifie_source(DoubleTab& so_modif,const DoubleTab& so_val,const DoubleTab& rho_faces,
                      const int n,const int m, const int is_QC,
                      const DoubleVect& vol_entrelaces,const Solveur_Masse& solv_masse);

  void calcul_effort_fluide_interface(const DoubleTab& vpoint,const DoubleTab& rho_faces,
                                      DoubleTab& source_val,const int is_explicite,const double eta);

  void impr_effort_fluide_interface( DoubleTab& source_val, DoubleTab& pressure_part, DoubleTab& friction_part ) ;

  //Calcul la vitesse imposee a l interface a partir de expression_vitesse_imposee
  virtual void calcul_vitesse(DoubleTab& vitesse_imp, const DoubleTab& champ_vitesse,
                              const DoubleTab& vpoint, const double temps, const double dt);
  virtual void get_expression_vitesse_imposee(DoubleTab& vitesse_imp);
  //Effectue l integration d un ensemble de points (sans notion de facettes)
  void integrer_ensemble_lagrange(const double temps) override;

  virtual void interpoler_vitesse_face(const DoubleTab& distance_interface,
                                       const int phase, const int stencil_width,
                                       DoubleTab& champ, DoubleTab& gradient,
                                       const double t, const double dt ) ;

  virtual void calcul_nb_traverse(  const DoubleTab& xe, const double dx,
                                    const int dim, const int ori,
                                    Maillage_FT_Disc& maillage, int elem,
                                    int& traverse ) ;
  virtual void calcul_OutElemFa7( Maillage_FT_Disc& maillage,
                                  const DoubleTab& indicatrice,
                                  const int nb_elem,
                                  int& nb_fa7_accepted,
                                  IntList& OutElem, IntList& OutFa7 ) ;
  virtual void PPP_face_interface( Maillage_FT_Disc& maillage, const DoubleTab& indicatrice,
                                   const DoubleTab& indicatrice_face, DoubleTab& Vertex ) ;

  virtual void PPP_face_interface_voisin( const DoubleTab& indicatrice, const DoubleTab& indicatrice_face,
                                          DoubleTab& Vertex, DoubleTab& PPP ) ;
  virtual void PPP_face_voisin( const DoubleTab& indicatrice, const DoubleTab& indicatrice_face, DoubleTab& PPP ) ;

  virtual void calcul_maxfa7( Maillage_FT_Disc& maillage, const DoubleTab& indicatrice,
                              const int nb_elem, int& max_fa7, const int exec_planfa7existan ) ;
  virtual void RenumFa7( DoubleTab& Barycentre, DoubleTab& Tab110,DoubleTab& Tab111,
                         DoubleTab& Tab112, IntTab& Tab12, IntTab& CptFacette,
                         const int nb_facettes, const int nb_facettes_dim ) ;
  virtual void StockageFa7( Maillage_FT_Disc& maillage, IntTab& CptFacette, DoubleTab& Tab100,
                            DoubleTab& Tab101,DoubleTab& Tab102, DoubleTab& Tab103,
                            DoubleTab& Tab110,DoubleTab& Tab111,DoubleTab& Tab112,
                            IntTab& Tab12, DoubleTab& Barycentre, const DoubleTab& indicatrice,
                            IntList& OutElem, ArrOfBit& fa7, const int exec_planfa7existant) ;
  virtual void StockageFa7( Maillage_FT_Disc& maillage, DoubleTab& Tab100, DoubleTab& Tab101,
                            DoubleTab& Tab102, DoubleTab& Tab103, DoubleTab& Tab110,
                            DoubleTab& Tab111, DoubleTab& Tab112, IntTab& Tab12,
                            DoubleTab& Barycentre, IntList& OutElem, IntTab& TabOutFa7, ArrOfBit& fa7 ) ;
  virtual void BaryFa7( Maillage_FT_Disc& maillage, const int i_facette, DoubleTab& Barycentre ) ;
  virtual void plan_facette_existant( Maillage_FT_Disc& maillage,
                                      DoubleList A, DoubleList B, DoubleList C,
                                      DoubleList D, const int i_facette,
                                      int& test_liste ) ;
  virtual void calcul_eq_plan_facette(Maillage_FT_Disc& maillage, const int i_facette,
                                      double& a, double& b, double& c, double& d);
  virtual void calcul_eq_plan_facette(const int i_facette,
                                      double& a, double& b, double& c, double& d);

  virtual void calcul_tolerance_projete_monophasique( const int i_face, const int ori, const int voisin0,
                                                      const int voisin1, const DoubleTab& indicatrice_face,
                                                      const DoubleTab& indicatrice, double& tol ) ;

  virtual void calcul_tolerance_projete_diphasique( const int i_face, const int ori, const int voisin0,
                                                    const int voisin1, const DoubleTab& indicatrice, double& tol ) ;

  void verifprojete(const int monophasique,const double Lref, double d, const DoubleTab& x,
                    const DoubleTab& V, DoubleTab& coord_projete, int& cpt ) ;


  virtual void uzawa(const double d,const DoubleTab& matrice, const DoubleTab& x,
                     const DoubleTab& secmem, DoubleTab& solution) const ;

  virtual void projete_point_face_fluide( int& nb_proj_modif, const int dim_fa7,
                                          const DoubleTab& indicatrice_face, const DoubleTab& indicatrice,
                                          const DoubleTab& dist_face, const double t, const double dt,
                                          DoubleTab& Tab100, DoubleTab& Tab101,DoubleTab& Tab102,
                                          DoubleTab& Tab103, IntTab& Tab12, IntTab& CptFacette,
                                          DoubleTab& v_imp, DoubleTab& Vertex,
                                          Parser& parser_x, Parser& parser_y,Parser& parser_z );
  virtual void projete_point_face_interface(   int& nb_proj_modif,const int dim_fa7,
                                               const DoubleTab& indicatrice_face,
                                               const DoubleTab& indicatrice,
                                               const DoubleTab& dist_face, const double t,
                                               const double dt, DoubleTab& Tab100,
                                               DoubleTab& Tab101,DoubleTab& Tab102,
                                               DoubleTab& Tab103, IntTab& Tab12,
                                               IntTab& CptFacette, DoubleTab& v_imp, DoubleTab& Vertex,
                                               Parser& parser_x, Parser& parser_y,Parser& parser_z) ;

  virtual void transporter_sans_changement_topologie(DoubleTab& vitesse,
                                                     const double coeff,const double temps);

  virtual int calculer_composantes_connexes_pour_suppression(IntVect& num_compo);
  virtual double suppression_interfaces(const IntVect& num_compo, const ArrOfInt& flags_compo_a_supprimer);

  const int& get_vimp_regul() const;

  virtual const Champ_base& get_update_distance_interface() const;
  virtual const Champ_base& get_update_distance_interface_faces() const;
  virtual const Champ_base& get_update_normale_interface() const;
  // renvoie DoubleTab parce qu'il n'existe pas de champ aux sommets en VDF ! Zut...
  virtual const DoubleTab&   get_update_distance_interface_sommets() const;
  void ramasse_miettes(const Maillage_FT_Disc& maillage,
                       DoubleVect& flux,
                       DoubleVect& valeurs);

protected:

  virtual void calculer_vmoy_composantes_connexes(const Maillage_FT_Disc& maillage,
                                                  const ArrOfInt& compo_connexes_facettes,
                                                  const int nb_compo_tot,
                                                  const DoubleTab& vitesse_sommets,
                                                  DoubleTab& vitesses,
                                                  DoubleTab& positions) const;

  virtual void deplacer_maillage_ft_v_fluide(const double temps);

  virtual void calculer_distance_interface(const Maillage_FT_Disc& maillage,
                                           DoubleTab& distance_elements,
                                           DoubleTab& normale_elements,
                                           const int n_iter) const;

  virtual void calculer_distance_interface_sommets(const DoubleTab& dist_elem,
                                                   const DoubleTab& normale_elem,
                                                   DoubleTab&        dist_som) const;


  virtual void calculer_vitesse_repere_local(const Maillage_FT_Disc& maillage,
                                             DoubleTab& deplacement,
                                             DoubleTab& Positions,
                                             DoubleTab& Vitesses) const;
  virtual void test_suppression_interfaces_sous_zone();


  virtual void calculer_distance_interface_faces(const DoubleTab& dist_elem,
                                                 const DoubleTab& normale_elem,
                                                 DoubleTab&        dist_faces) const;
  // Nouvelle methodes
  // Methode outil utilisee par modifier_vpoint_pour_imposer_vit(...)
  // Calcul l'indicatrice sur chaque face
  void calcul_indicatrice_faces(const DoubleTab& indicatrice,
                                const IntTab& face_voisins);

  REF(Probleme_base) probleme_base_;
  REF(Navier_Stokes_FT_Disc) equation_ns_;
  // L'inconnue du probleme
  Champ_Inc indicatrice_;
  Champ_Inc indicatrice_faces_;
  // Utiliser ces accesseurs :
  Maillage_FT_Disc& maillage_interface();
  // Utiliser ces accesseurs :
  Marching_Cubes& marching_cubes();
  // Utiliser ces accesseurs :
  Topologie_Maillage_FT& topologie_interface();

  Proprietes_part_vol&         proprietes_particules();
  Maillage_FT_Disc&            maillage_inject();
  Proprietes_part_vol&         proprietes_inject();

  DoubleTab& tableaux_positions();
  IntVect& vecteur_elements();
  DoubleTab&    deplacement_som();

  // On utilise des DERIV() pour ne pas avoir a inclure la definition
  // de ces classes (pour reduire les dependances).
  static void transfert_conservatif_eulerien_vers_lagrangien_sommets(const Maillage_FT_Disc& maillage,
                                                                     const DoubleVect& valeurs_euler,
                                                                     ArrOfDouble& valeurs_lagrange);

  Nom suppression_interfaces_sous_zone_;

  Champ_Fonc vitesse_imp_interp_;


private:
  // Variables internes a la methode de transport
  Transport_Interfaces_FT_Disc_interne *variables_internes_;


  double temps_debut_;

  REF(Milieu_base) ref_milieu_;

  int interpolation_repere_local_;
  ArrOfDouble force_;
  ArrOfDouble moment_;
};

class Transport_Interfaces_FT_Disc_interne : public Objet_U
{
  Declare_instanciable_sans_constructeur_ni_destructeur(Transport_Interfaces_FT_Disc_interne);
public:
  Transport_Interfaces_FT_Disc_interne() :
    indicatrice_cache_tag(-1),
    iterations_correction_volume(0),
    volume_impose_phase_1(-1.),
    n_iterations_distance(3),
    n_iterations_interpolation_ibc(5),
    distance_normale_cache_tag(-1),
    distance_sommets_cache_tag(-1),
    distance_faces_cache_tag(-1),
    methode_transport(INDEFINI),
    methode_interpolation_v(VALEUR_A_ELEM),
    injection_interfaces_last_time_(0.),
    interpolation_champ_face(BASE),
    vf_explicite(0),
    is_extra_diphasique_solide(0),
    is_extra_solide(0),
    is_distance_projete_face(0),
    nomb_fa7_accepted(3),
    seuil_uzawa(1.e-8),
    nb_iter_uzawa(30),
    vimp_regul(1),
    indic_faces_modif(0),
    modified_indic_faces_position(0.),
    modified_indic_faces_thickness(1.),
    type_vitesse_imposee(UNIFORME),
    type_distance_calculee(DIST_INITIALE),
    type_projete_calcule(PROJETE_INITIAL),
    expression_vitesse_imposee(Objet_U::dimension)


  {};
  ~Transport_Interfaces_FT_Disc_interne() override
  {};
private:
  Transport_Interfaces_FT_Disc_interne(const Transport_Interfaces_FT_Disc_interne& a): Objet_U(a)
  {}; // Interdit
  const Transport_Interfaces_FT_Disc_interne& operator=(const Transport_Interfaces_FT_Disc_interne& a)
  {
    return *this;
  }; // Interdit

public:
  int sauvegarder(Sortie& os) const override;
  int reprendre(Entree& is) override;

  // Les membres suivantes sont sauvegardes et repris:
  Champ_Inc        indicatrice_cache;     // L'indicatrice calculee par get_update_indicatrice
  int           indicatrice_cache_tag; // Le tag du maillage correspondant
  Maillage_FT_Disc maillage_interface;          // Objet qui peut se reduire a un ensemble de sommets
  // quand il represente les positions de particules
  Remaillage_FT    remaillage_interface_;
  // Fin des membres sauvegardes / repris

  Proprietes_part_vol proprietes_particules_; //Proprietes physiques de particules

  Maillage_FT_Disc maillage_inject_;              //Ensemble de particules a injecter periodiquement
  Proprietes_part_vol proprietes_inject_;     //Proprietes physiques des particules injectees

  // Si iterations_correction_volume == 0, le maillage est transporte par le fluide
  // a l'aide du champ de vitesse L2 interpole (pas de conservation du volume).
  // Si iterations_correction_volume > 0, on calcule une correction de volume aux
  // sommets de l'interface et on l'etale par autant d'iterations d'un lisseur
  // Voir Transport_Interfaces_FT_Disc::mettre_a_jour
  int iterations_correction_volume;
  // Rustine introduite pour corriger les pertes de masse:
  // Si cette valeur est positive, on deplace toute l'interface d'une certaine distance
  // pour que le volume de la phase 1 reste toujours egal a la valeur prescrite
  double volume_impose_phase_1;
  // Si non nul, le calcul de l'integrale pour le volume impose porte sur cette sous-zone
  Nom    nom_zone_volume_impose_;

  Champ_Inc vitesse_filtree;
  DoubleTab doubletab_pos;
  DoubleTab doubletab_vitesses;
  IntVect   intvect_elements;

  Champ_Fonc distance_interface; // Distance a l'interface (aux elements)
  Champ_Fonc normale_interface;  // Une normale etalee
  Champ_Fonc surface_interface;  // GB : La surface d'interface dans chaque element
  Champ_Fonc tmp_flux;           // Tableau temporaire pour le ramasse-miettes
  Champ_Fonc distance_interface_faces; // CF : Distance a l'interface (aux faces)
  DoubleTab  distance_interface_sommets; // Distance a l'interface (aux sommets)
  Champ_Fonc distance_interface_faces_corrigee; // CI : Distance a l'interface corrigee (aux faces)
  Champ_Fonc distance_interface_faces_difference; // CI : Distance a l'interface corrigee - Distance a l'interface calculee (aux faces)
  Champ_Fonc index_element; // CI : indexation des elements
  Champ_Fonc nelem_par_direction; // CI : nombre d'elements par direction
  // Note de B.M. : zut, y'a pas de champ aux sommets en VDF... donc DoubleTab et
  // du coup on ne peut pas postraiter facilement. C'est trop con...
  int     n_iterations_distance;
  int     n_iterations_interpolation_ibc;
  int     distance_normale_cache_tag;
  int     distance_sommets_cache_tag;
  int     distance_faces_cache_tag;
  // Le maillage postraite n'est pas forcement le maillage
  Maillage_FT_Disc maillage_pour_post;
  DoubleTabFT   deplacement_sommets;

  enum Methode_transport { INDEFINI, VITESSE_IMPOSEE, LOI_HORAIRE, VITESSE_INTERPOLEE };
  Methode_transport      methode_transport;
  REF(Navier_Stokes_std) refequation_vitesse_transport;

  enum Methode_interpolation_v { VALEUR_A_ELEM, VDF_LINEAIRE };
  Methode_interpolation_v methode_interpolation_v;


  // Injecteur d'interfaces au cours du temps
  ArrOfDouble injection_interfaces_temps_;
  ArrOfInt    injection_interfaces_phase_;
  Noms        injection_interfaces_expressions_;
  // temps physique reel auquel a eu lieu la derniere injection (ce n'est pas le temps demande)
  double      injection_interfaces_last_time_;

  enum Interpolation_champ_face { BASE, LINEAIRE };
  Interpolation_champ_face interpolation_champ_face;

  int vf_explicite, is_extra_diphasique_solide, is_extra_solide, is_distance_projete_face;
  int nomb_fa7_accepted ;
  double seuil_uzawa ;
  int nb_iter_uzawa ;
  int vimp_regul ;
  int indic_faces_modif ;
  double modified_indic_faces_position ;
  double modified_indic_faces_thickness ;

  enum Type_vitesse_imposee { UNIFORME, ANALYTIQUE };
  Type_vitesse_imposee type_vitesse_imposee;

  enum Type_distance_calculee { DIST_INITIALE, DIST_MODIFIEE };
  Type_distance_calculee type_distance_calculee;

  enum Type_projete_calcule { PROJETE_INITIAL, PROJETE_MODIFIE, PROJETE_SIMPLIFIE };
  Type_projete_calcule type_projete_calcule;

  Noms                      expression_vitesse_imposee;
  // Reference a une loi horaire eventuelle
  REF(Loi_horaire)         loi_horaire_;

  // Integration de la vitesse a partir du point de depart (x,y,z)
  //  pendant un temps dt.
  void integrer_vitesse_imposee(double temps, double dt, double& x, double& y, double& z) const;

  // Les objets-algorithmes :
  Parcours_interface      parcours_interface_;
  Marching_Cubes          marching_cubes_;
  Connectivite_frontieres connectivite_frontieres_;
  Topologie_Maillage_FT   topologie_interface_;
  // Cet objet est type en fonction de la discretisation:
  DERIV(Algorithmes_Transport_FT_Disc) algorithmes_transport_;
};
#endif
