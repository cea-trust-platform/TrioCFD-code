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
// File      : ComputeValParCompoInCell.h
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#ifndef ComputeValParCompoInCell_included
#define ComputeValParCompoInCell_included

#include <Connectivite_frontieres.h>
#include <FixedVector.h>
#include <IJK_Field.h> // est-ce que j'en ai vraiment besoin ?
#include <Linear_algebra_tools_impl.h>
#include <Maillage_FT_IJK.h>
#include <Parcours_interface.h>
#include <Vecteur3.h>
#include <Linear_algebra_tools_impl.h>
#include <SurfaceVapeurIJKComputation.h>

class ComputeValParCompoInCell
{
public:
  ComputeValParCompoInCell() { mesh_ = nullptr;}
  ~ComputeValParCompoInCell() {}
  void initialize(const IJK_Splitting& splitting, const Maillage_FT_IJK& maillage_ft_ijk)
  {
    ref_splitting_ = splitting;
    mesh_ = &maillage_ft_ijk;
  }

  void calculer_valeur_par_compo(
#ifdef SMOOTHING_RHO
    const double delta_rho,
#endif
    const double time,
    const int itstep,
    IJK_Field_int& nb_compo_trav,
    FixedVector<IJK_Field_int, max_authorized_nb_of_components_>& compos_trav,
    FixedVector<IJK_Field_double, 3*max_authorized_nb_of_components_>& normale_par_compo,
    FixedVector<IJK_Field_double, 3*max_authorized_nb_of_components_>& bary_par_compo,
    FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& indicatrice_par_compo,
    FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& surface_par_compo,
    FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& courbure_par_compo
  );

  void calculer_moy_field_sommet_par_compo(
    const ArrOfDouble& val_on_sommet,
    FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& field_par_compo
  ) const;

  void calculer_moy_field_fa7_par_compo(
    const ArrOfDouble& val_on_fa7,
    FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& field_par_compo
  ) const;

protected:
  // Calcul des tableau de surface, normal et bary par compo
  void calculer_moyennes_interface_element_pour_compo(
    const int num_compo,
    const int elem,
    double& surface,
    Vecteur3& normale,
    Vecteur3& bary
  ) const;

  int calculer_indic_elem_pour_compo(const int icompo, const int elem, double& indic) const;

  void calculer_moy_par_compo(
#ifdef SMOOTHING_RHO
    const double delta_rho,
#endif
    IJK_Field_int& nb_compo_traversante,
    FixedVector<IJK_Field_int, max_authorized_nb_of_components_>& compos_traversantes,
    FixedVector<IJK_Field_double, 3 * max_authorized_nb_of_components_>& normale_par_compo,
    FixedVector<IJK_Field_double, 3 * max_authorized_nb_of_components_>& bary_par_compo,
    FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& indic_par_compo,
    FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& surface_par_compo
  ) const;

  int compute_list_compo_connex_in_element(
    const int elem,
    ArrOfInt& liste_composantes_connexes_dans_element) const;

  OBS_PTR(IJK_Splitting) ref_splitting_;
  const Maillage_FT_IJK * mesh_;
};

#endif
