/****************************************************************************
 * Copyright (c) 2019, CEA
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 *modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *this list of conditions and the following disclaimer in the documentation
 *and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *may be used to endorse or promote products derived from this software without
 *specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *****************************************************************************/
/////////////////////////////////////////////////////////////////////////////
//
// File      : Corrige_flux_FT.h
// Directory : $IJK_ROOT/src/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Corrige_flux_FT_included
#define Corrige_flux_FT_included

#include <MonofluidVar.h>
#include <Boundary_Conditions_Thermique.h>
#include <IJK_Splitting.h>
#include <Objet_U.h>
#include <Parser.h>
#include <IJK_Interfaces.h>
#include <IJK_Lata_writer.h>
#include <Intersection_Interface_ijk.h>
#include <Ouvrir_fichier.h>
#include <Ref_IJK_FT_double.h>

/*! @brief : class Corrige_flux_FT
 *
 *  <Description of class Corrige_flux_FT>
 *
 *
 *
 */
class ParcoursIJKDir
{
public:
  ParcoursIJKDir() {next_elem_.resize(3); indices_to_keep_.resize(2);};
  ~ParcoursIJKDir() {};
  // void initialize();
  FixedVector<int, 3> elem( const int q ) const
  {
    FixedVector<int, 3> res;
    res[0] = i_ + q * next_elem_[0];
    res[1] = j_ + q * next_elem_[1];
    res[2] = k_ + q * next_elem_[2];
    return res;
  }
  const int& i() const {return i_;}
  const int& j() const {return j_;}
  const int& k() const {return k_;}
  const int& face() const {return dir_;}
  void set_dir(const int dir)
  {
    dir_ = dir;
    set_next_elem();
    set_indices_to_keep();
  }
  void set_elem(const int i1, const int j1, const int k1)
  {
    i_ = i1;
    j_ = j1;
    k_ = k1;
  }
  const ArrOfInt& get_indices_to_keep() const { return indices_to_keep_;}
  const ArrOfInt& get_normale_vec() const {return next_elem_;}

  double calculer_surface_face(const IJK_Grid_Geometry& geom) const
  {
    double surface = 1.;
    for (int c = 0; c < 2; c++)
      surface *= geom.get_constant_delta(indices_to_keep_[c]);
    return surface;
  }

protected:
  void set_next_elem();
  void set_indices_to_keep();
  int dir_;
  int i_;
  int j_;
  int k_;
  ArrOfInt next_elem_;
  ArrOfInt indices_to_keep_;
};

// Declare_ref(Corrige_flux_FT); TODO: a mettre en place éventuellement

// API pour modifier un champ de flux à partir de donnees à l'interface. Cette
// classe est abstraite, elle est vouée à être héritée pour être adaptée aux
// différents cas / conditions aux limites. Par ex sur la température pour
// imposer la continuité de lda grad T à l'interface. Mais ça pourrait être
// pareil pour la vitesse pour imposer la continuité tangentielle et normale
// dans le cas sans changement de phase, ou pour imposer la bonne différence de
// vitesse normale dans le cas du changement de phase.
class Corrige_flux_FT // : public Objet_U
{
  // Declare_instanciable( Corrige_flux_FT ) ;
public:
  Corrige_flux_FT() {};
  virtual ~Corrige_flux_FT() {};
  void initialize(const IJK_Splitting& splitting,
                  const IJK_Field_double& field,
                  const IJK_Interfaces& interfaces, const IJK_FT_double& ijk_ft);
  // On va calculer sur la grille IJ de la layer k tous les flux a proximite de
  // l'interface. On remplace les flux donnes en entree par ces flux la.
  virtual void corrige_flux_faceIJ(IJK_Field_local_double *const flux,
                                   const int k_layer, const int dir) = 0;
  virtual void update() = 0;

protected:
  const IJK_Interfaces *interfaces_;
  const IJK_Field_double *field_;
  const IJK_Splitting *splitting_;
  REF(IJK_FT_double) ref_ijk_ft_;

  // TODO: mettre ces méthodes dans une petite classe pour parcourir
  // les trois directions.

  bool test_if_stencil_inclut_bout_interface_liquide() const ;
  bool test_if_stencil_inclut_bout_interface_vapeur() const ;
  ParcoursIJKDir parcours_;
};


// Surclasse Corrige Flux FT dans le cas de la température sans changement de
// phase. Cad avec continuite normale de lda grad T. Ici le flux corrige est le
// flux convectif.
class Corrige_flux_FT_temperature_conv : public Corrige_flux_FT
{
public:
  Corrige_flux_FT_temperature_conv() {};
  virtual ~Corrige_flux_FT_temperature_conv() {};

  void initialize(const IJK_Splitting& splitting, const IJK_Field_double& field,
                  const IJK_Interfaces& interfaces, const IJK_FT_double& ijk_ft,
                  const double rhocp_l, const double rhocp_v, const double lda_l,
                  const double lda_v);

  void corrige_flux_faceIJ(IJK_Field_local_double *const flux,
                           const int k_layer, const int dir) override;

  void update() override;

  // Cette méthode calcule la temperature et le flux a l'interface aux points
  // donnes en entree. Dans le cas classique elle est utilisee aux points qui
  // sont les projections des centres de gravites des faces mouillees coupees
  // par l'interface. Elle est aussi utilisee sur tous les centres de facette
  // pour les post traitement de ces valeurs.
  static void calcul_temperature_flux_interface(
    const IJK_Field_double& temperature, const double ldal, const double ldav,
    const double dist, const DoubleTab& positions, const DoubleTab& normale,
    ArrOfDouble& temperature_interp, ArrOfDouble& flux_normal_interp,
    ArrOfDouble& temp_liqu, ArrOfDouble& temp_vap, DoubleTab& coo_liqu,
    DoubleTab& coo_vap);

protected:
  // Calcule effectivement sur les points a l'interface les valeurs de temp et
  // flux. Met a jour les tableaux temp_interface_ et q_interface_.
  void calcul_temp_flux_interf_pour_bary_face(ArrOfDouble& temp_vap, ArrOfDouble& temp_liqu);

  // Calcule effectivement sur les points a l'interface les valeurs de temp et
  // flux. Met a jour les tableaux temp_interface_ et q_interface_.
  void calcul_temp_flux_interf_pour_bary_cell(ArrOfDouble& temp_vap, ArrOfDouble& temp_liqu);

  // Ici on suppose que les methodes d interpolation de la temperature ont deja
  // ete executees a ce pas de temps et on interpole dans l autre sens pour
  // recuperer la valeur aux barycentres des faces. Dans cette methode pas
  // besoin de faire comme pour la premiere partie qui sert aussi au post
  // traitement. Ici on ne s'en sert strictement que dans l'operateur.
  void interp_back_to_bary_faces(const ArrOfDouble& temp_vap, const ArrOfDouble& temp_liqu);

  // Méthode qui met à jour temperature_ghost_.
  void update_temperature_ghost(const ArrOfDouble& temp_vap, const ArrOfDouble& temp_liqu);

  // Comme ca on ne peut plus utiliser par erreur l'init de la classe abstraite
  // en cas de polymorphisme.
  using Corrige_flux_FT::initialize;

  // On fait référence au champ thermique pour accéder surtout aux valeurs des
  // propriétés thermiques.
  // TODO: Il faudrait créer un REF(IJK_Energie)
  double rhocp_l_, rhocp_v_;
  double lda_l_, lda_v_;
  IJK_MonofluidVar rho_cp_;
  IJK_MonofluidVar lda_;

  // C'est un peu bizarre de le mettre ici, mais c'est probablement le plus
  // simple. Cet objet calcul / stocke les coo de projection entre les barys des
  // faces mouillées et l'interface.
  Intersection_Interface_ijk_face intersection_ijk_face_;
  // Les tableaux pour stocker les valeurs aux points de part et d'autre de
  // l'interface.
  // ArrOfDouble temp_vap1_, temp_liqu1_;
  // On prepare des tableaux pour les coordonnees des points d'interpolation et
  // des valeurs. Ces valeurs sont stockées pour éventuel post-traitement.
  ArrOfDouble temp_interface_face_;
  ArrOfDouble q_interface_face_;

  ArrOfDouble temp_interface_cell_;
  ArrOfDouble q_interface_cell_;
  // Ce tableau est de taille (n_diph, 2)
  // La premiere compo donne la température de la phase liquide et la deuxième
  // celle de la phase vapeur.
  DoubleTab temperature_barys_;

  // Cet objet calcul / stocke les coo de projection des centres des
  // cellules diphasique sur l'interface.
  Intersection_Interface_ijk_cell intersection_ijk_cell_;
  // Dans ce tableau on stocke les température liquid et vapeur ghost des cell
  // diph.
  DoubleTab temperature_ghost_;

  void multiplie_par_rho_cp_de_la_face_monophasique(
    const double frac_liquide,
    IJK_Field_local_double * const flux
  ) const;

  void remplace_flux_par_quick_ghost_amont_1(
    const double frac_liquide,
    const double s_face,
    IJK_Field_local_double * const flux
  ) const ;
  void remplace_flux_par_somme_rhocpf_Tf_v_Sf(
    const double frac_liquide,
    const double s_face,
    IJK_Field_local_double  * const flux) const ;
  double quick(
    const double Tim1,
    const double Ti,
    const double Tip1,
    const double Tip2,
    const double velocity) const ;
  const double& get_ghost_temp_if_cell_is_diph(
    const FixedVector<int, 3>& elem,
    const bool from_liqu_phase) const ;
  double extrapolation_amont_1_depuis_l_interface(
    const double frac_liquide, const double decal) const ;
  double interpolation_quick_avec_1_ghost(
    const double frac_liquide,
    const double decal) const ;
  bool is_flux_upwind_from_interface(const double decal) const ;

};

#endif /* Corrige_flux_FT_included */
