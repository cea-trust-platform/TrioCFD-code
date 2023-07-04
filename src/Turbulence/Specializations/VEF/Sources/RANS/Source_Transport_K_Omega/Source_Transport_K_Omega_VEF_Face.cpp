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
// File:        Source_Transport_K_Omega_VEF_Face.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_K_Omega_VEF_Face.h>
#include <Operateur_Grad.h>
#include <TRUSTTabs_forward.h>
#include <Source_Transport_K_Omega_VEF_Face.h>
#include <Modele_turbulence_hyd_K_Omega.h>
#include <Transport_K_Omega.h>
#include <Navier_Stokes_Turbulent.h>
#include <Milieu_base.h>
#include <Domaine_VEF.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Omega_VEF_Face,
                                          "Source_Transport_K_Omega_VEF_P1NC",
                                          Source_Transport_K_Omega_VEF_Face_base);

Sortie& Source_Transport_K_Omega_VEF_Face::printOn(Sortie& s) const { return s << que_suis_je() ; }

Entree& Source_Transport_K_Omega_VEF_Face::readOn(Entree& is)
{
  Source_Transport_K_Omega_VEF_Face_base::verifier_pb_komega(mon_equation->probleme(),que_suis_je());
  return Source_Transport_K_Omega_VEF_Face_base::readOn(is);
}

void Source_Transport_K_Omega_VEF_Face::associer_pb(const Probleme_base& pb)
{
  Source_Transport_K_Omega_VEF_Face_base::associer_pb(pb);
  eqn_K_Omega = ref_cast(Transport_K_Omega, equation());
}
const DoubleTab& Source_Transport_K_Omega_VEF_Face::get_visc_turb() const
{
  return eqn_K_Omega->modele_turbulence().viscosite_turbulente().valeurs();
}

const DoubleTab& Source_Transport_K_Omega_VEF_Face::get_cisaillement_paroi() const
{
  const Modele_turbulence_hyd_K_Omega& mod = ref_cast(Modele_turbulence_hyd_K_Omega,
                                                      eqn_K_Omega->modele_turbulence());
  return mod.loi_paroi().valeur().Cisaillement_paroi();
}

const DoubleTab& Source_Transport_K_Omega_VEF_Face::get_K_pour_production() const
{
  return eqn_K_Omega->inconnue().valeurs();
}

const Nom Source_Transport_K_Omega_VEF_Face::get_type_paroi() const
{
  const Modele_turbulence_hyd_K_Omega& mod = ref_cast(Modele_turbulence_hyd_K_Omega,
                                                      eqn_K_Omega->modele_turbulence());
  return mod.loi_paroi().valeur().que_suis_je();
}

void Source_Transport_K_Omega_VEF_Face::compute_cross_diffusion(DoubleTab gradKgradOmega) const
{
  // Same structure than Source_WC_Chaleur_VEF.cpp in TRUST
  // A mettre dans la base ?

  // = First, store K and Omega in two separated Tab for the gradient computation, not ideal
  const DoubleTab& K_Omega = eqn_K_Omega->inconnue().valeurs();
  DoubleTab enerK; // field on faces
  DoubleTab omega; // field on faces
  enerK.resize(K_Omega.dimension_tot(0));
  omega.resize(K_Omega.dimension_tot(0));
  for (int face = 0; face < le_dom_VEF->nb_faces(); ++face)
    {
      enerK(face) = K_Omega(face, 0);
      omega(face) = K_Omega(face, 1);
    }

  // = Definition of the gradients
  DoubleTab gradK; // field on elements
  DoubleTab gradOmega; // field on elements

  // Get number of componants, for gradient resize
  // We use the velocity field to resize the gradient as the velocity is on faces
  const Navier_Stokes_Turbulent& eqHyd = ref_cast(Navier_Stokes_Turbulent,
                                                  mon_equation->probleme().equation(0));
  const DoubleTab& la_vitesse = eqHyd.vitesse().valeurs(); // Velocity on faces
  const int nb_compo = la_vitesse.line_size();
  const DoubleTab &pressure = eqHyd.pression().valeurs();
  const int nb_tot = pressure.size_totale(); // find a better name than nb_tot
  gradK.resize(nb_tot, nb_compo);
  gradOmega.resize(nb_tot, nb_compo);
  // resize_gradient_tab(gradK);
  // resize_gradient_tab(gradOmega);

  // Compute the two gradients
  const Operateur_Grad Op_Grad_komega = eqn_K_Omega.gradient_operateur_komega();
  Op_Grad_komega.calculer(enerK, gradK);
  Op_Grad_komega.calculer(omega, gradOmega);

  // Correction on the boundaries? Elie put the pressure at zero.

  // Face to elem
  const Domaine_dis_base& domaine_dis = mon_equation->inconnue().domaine_dis_base();
  const Domaine_VF& domaine = ref_cast(Domaine_VF, domaine_dis);
  DoubleTab gradK_face(la_vitesse); // gradK on faces is similar to the velocity
  elem_to_face(domaine, gradK, gradK_face);
  DoubleTab gradOmega_face(la_vitesse);
  elem_to_face(domaine, gradOmega, gradOmega_face);

  // Product gradKgradOmega
  // gradKgradOmega.resize(nb_tot, nb_compo);
  // Si on a gradKgradOmega en argument, il faut s'assurer qu'il est bien dimensionné à l'initialisation
  for (int face = 0; face < le_dom_VEF->nb_faces(); ++face)
    for (int ncompo = 0; ncompo < nb_compo; ++ncompo)
      gradKgradOmega(face, ncompo) += gradK(face, ncompo)*gradOmega_face(face, ncompo);

}

// cAlan, 2023-06-23: salement copié de Source_Chaleur_WC_VEF. À mutualiser.
void Source_Transport_K_Omega_VEF_Face::elem_to_face(const Domaine_VF& domaine,
                                                     const DoubleTab& grad_elems,
                                                     DoubleTab& grad_faces) const
{
  const DoubleVect& vol = domaine.volumes();
  const IntTab& elem_faces = domaine.elem_faces();
  const int nb_face_elem = elem_faces.line_size();
  const int nb_elem_tot = domaine.nb_elem_tot();
  const int nb_comp = grad_faces.line_size();

  assert (grad_elems.dimension_tot(0) == nb_elem_tot);
  assert (grad_faces.dimension_tot(0) == domaine.nb_faces_tot());
  assert (grad_elems.line_size() == nb_comp);

  grad_faces = 0.;
  for (int elem = 0; elem < nb_elem_tot; ++elem)
    for (int s = 0; s < nb_face_elem; ++s)
      {
      const int face = elem_faces(elem, s);
      for (int comp = 0; comp < nb_comp; ++comp)
        grad_faces(face, comp) += grad_elems(elem, comp) * vol(elem);
      }

  const DoubleVect& volumes_entrelaces = le_dom_VEF->volumes_entrelaces();
  for (int f = 0; f < domaine.nb_faces_tot(); ++f)
    for (int comp = 0; comp < nb_comp; ++comp)
      grad_faces(f, comp) /= volumes_entrelaces(f)*nb_face_elem;
}

// cAlan: Tried to get a dedicated function to resize a tab. Make a function for this?
void Source_Transport_K_Omega_VEF_Face::resize_gradient_tab(DoubleTab &grad) const
{
  const Navier_Stokes_Turbulent& eqHyd = ref_cast(Navier_Stokes_Turbulent,
                                                  mon_equation->probleme().equation(0));
  const DoubleTab& la_vitesse = eqHyd.vitesse().valeurs();
  const int nb_compo = la_vitesse.line_size();
  const DoubleTab &pressure = eqHyd.pression().valeurs();
  const int nb_tot = pressure.size_totale(); // find a better name than nb_tot
  grad.resize(nb_tot, nb_compo);
}

void Source_Transport_K_Omega_VEF_Face::fill_resu(const DoubleVect &volumes_entrelaces,
                                                  const DoubleTrav &ProdK, DoubleTab &resu) const {
  const DoubleTab& K_Omega = eqn_K_Omega->inconnue().valeurs();
  const double LeK_MIN = eqn_K_Omega->modele_turbulence().get_K_MIN();
  // cAlan : mettre le calcul du gradKgradOmega ici ?
  for (int face = 0; face < le_dom_VEF->nb_faces(); face++)
    {
      resu(face, 0) += (ProdK(face) - BETA_K*K_Omega(face, 0)*K_Omega(face, 1))*volumes_entrelaces(face);
      if (K_Omega(face, 0) >= LeK_MIN)
        resu(face, 1) += (ALPHA_OMEGA*ProdK(face)*K_Omega(face, 1)/K_Omega(face, 0)
                          - BETA_OMEGA*K_Omega(face, 1)*K_Omega(face, 1))*volumes_entrelaces(face);
    }
}

DoubleTab& Source_Transport_K_Omega_VEF_Face::ajouter(DoubleTab& resu) const
{
  return Source_Transport_K_Omega_VEF_Face_base::ajouter_komega(resu);
}

void Source_Transport_K_Omega_VEF_Face::contribuer_a_avec(const DoubleTab& a,
                                                          Matrice_Morse& matrice) const
{
  const DoubleTab& K_Omega = equation().inconnue().valeurs();
  const double LeK_MIN = eqn_K_Omega->modele_turbulence().get_K_MIN();
  const DoubleVect& porosite_face = eqn_K_Omega->milieu().porosite_face();
  const DoubleVect& volumes_entrelaces = le_dom_VEF->volumes_entrelaces();

  // on implicite le -eps et le -eps^2/k
  // cAlan: to be adapted for k_omega
  for (int face = 0; face < K_Omega.dimension(0); face++)
    if (K_Omega(face, 0) >= LeK_MIN) // -eps*vol  donne +vol dans la bonne case
      {
        const double volporo = porosite_face(face) * volumes_entrelaces(face);

        double coef_k = K_Omega(face, 1)/K_Omega(face, 0)*volporo;
        matrice(face*2, face*2) += coef_k;

        double coef_omega = ALPHA_OMEGA*coef_k;
        matrice(face*2 + 1, face*2 + 1) += coef_omega;
      }
}
