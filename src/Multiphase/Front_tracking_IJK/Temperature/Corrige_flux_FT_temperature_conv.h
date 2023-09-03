/****************************************************************************
* Copyright (c) 2023, CEA
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
// File      : Corrige_flux_FT_temperature_conv.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Corrige_flux_FT_temperature_conv_included
#define Corrige_flux_FT_temperature_conv_included

#include <Corrige_flux_FT_base.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Corrige_flux_FT_temperature_conv
//
// <Description of class Corrige_flux_FT_temperature_conv>
//
/////////////////////////////////////////////////////////////////////////////

class Corrige_flux_FT_temperature_conv : public Corrige_flux_FT_base
{

  Declare_instanciable( Corrige_flux_FT_temperature_conv ) ;

public:

  void initialize(const IJK_Splitting& splitting,
                  const IJK_Field_double& field,
                  const IJK_Interfaces& interfaces,
                  const IJK_FT_double& ijk_ft,
                  Intersection_Interface_ijk_face& intersection_ijk_face,
                  Intersection_Interface_ijk_cell& intersection_ijk_cell) override;

  void corrige_flux_faceIJ(IJK_Field_local_double *const flux,
                           const int k_layer, const int dir) override;

  void update() override;

  /* Cette méthode calcule la temperature et le flux a l'interface aux points
   * donnes en entree. Dans le cas classique elle est utilisee aux points qui
   * sont les projections des centres de gravites des faces mouillees coupees
   * par l'interface. Elle est aussi utilisee sur tous les centres de facette
   * pour les post traitement de ces valeurs.
   */
  void calcul_temperature_flux_interface(
    const IJK_Field_double& temperature, const double ldal, const double ldav,
    const double dist, const DoubleTab& positions, const DoubleTab& normale,
    ArrOfDouble& temperature_interp, ArrOfDouble& flux_normal_interp,
    ArrOfDouble& temp_liqu, ArrOfDouble& temp_vap, DoubleTab& coo_liqu,
    DoubleTab& coo_vap) const override;

protected:
  /* Calcule effectivement sur les points a l'interface les valeurs de temp et
   * flux. Met a jour les tableaux temp_interface_ et q_interface_.
   */
  void calcul_temp_flux_interf_pour_bary_face(ArrOfDouble& temp_vap, ArrOfDouble& temp_liqu);

  /* Calcule effectivement sur les points a l'interface les valeurs de temp et
   * flux. Met a jour les tableaux temp_interface_ et q_interface_.
   */
  void calcul_temp_flux_interf_pour_bary_cell(ArrOfDouble& temp_vap, ArrOfDouble& temp_liqu);

  /* Ici on suppose que les methodes d interpolation de la temperature ont deja
   * ete executees a ce pas de temps et on interpole dans l autre sens pour
   * recuperer la valeur aux barycentres des faces. Dans cette methode pas
   * besoin de faire comme pour la premiere partie qui sert aussi au post
   * traitement. Ici on ne s'en sert strictement que dans l'operateur.
   * void interp_back_to_bary_faces(const ArrOfDouble& temp_vap, const ArrOfDouble& temp_liqu);
   */
  void interp_back_to_bary_faces(const ArrOfDouble& temp_vap, const ArrOfDouble& temp_liqu);

  // Méthode qui met à jour temperature_ghost_.
  void update_temperature_ghost(const ArrOfDouble& temp_vap, const ArrOfDouble& temp_liqu);

  /* Comme ca on ne peut plus utiliser par erreur l'init de la classe abstraite
   * en cas de polymorphisme.
   */
//  using Corrige_flux_FT::initialize;

  /* On fait référence au champ thermique pour accéder surtout aux valeurs des
   * propriétés thermiques.
   * TODO: Il faudrait créer un REF(IJK_Energie)
   */
  IJK_MonofluidVar rho_cp_;
  IJK_MonofluidVar lda_;


  // FIXME: Move Intersection_Interface_ijk_face to IJK_Interfaces
  /* C'est un peu bizarre de le mettre ici, mais c'est probablement le plus
   * simple. Cet objet calcul / stocke les coo de projection entre les barys des
   * faces mouillées et l'interface.
   */
  // Intersection_Interface_ijk_face intersection_ijk_face_;

  /* Les tableaux pour stocker les valeurs aux points de part et d'autre de
   * l'interface.
   * ArrOfDouble temp_vap1_, temp_liqu1_;
   * On prepare des tableaux pour les coordonnees des points d'interpolation et
   * des valeurs. Ces valeurs sont stockées pour éventuel post-traitement.
   */
  ArrOfDouble temp_interface_face_;
  ArrOfDouble q_interface_face_;

  ArrOfDouble temp_interface_cell_;
  ArrOfDouble q_interface_cell_;
  /* Ce tableau est de taille (n_diph, 2)
   * La premiere compo donne la température de la phase liquide et la deuxième
   * celle de la phase vapeur.
   */
  DoubleTab temperature_barys_;

  // FIXME: Move Intersection_Interface_ijk_cell to IJK_Interfaces
  /* Cet objet calcul / stocke les coo de projection des centres des
   * cellules diphasique sur l'interface.
   */
  // Intersection_Interface_ijk_cell intersection_ijk_cell_;
  /* Dans ce tableau on stocke les température liquid et vapeur ghost des cell
   * diph.
   */
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
