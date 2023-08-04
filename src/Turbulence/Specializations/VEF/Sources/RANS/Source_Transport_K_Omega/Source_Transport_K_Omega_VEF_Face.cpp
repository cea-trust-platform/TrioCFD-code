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

#include <Operateur_Grad.h>
#include <TRUSTTabs_forward.h>
#include <Source_Transport_K_Omega_VEF_Face.h>
#include <Modele_turbulence_hyd_K_Omega.h>
#include <Transport_K_Omega.h>
#include <Navier_Stokes_Turbulent.h>
#include <Pb_Hydraulique_Turbulent.h>
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
  turbulence_model = ref_cast(Modele_turbulence_hyd_K_Omega, eqn_K_Omega->modele_turbulence());
}
const DoubleTab& Source_Transport_K_Omega_VEF_Face::get_visc_turb() const
{
  return eqn_K_Omega->modele_turbulence().viscosite_turbulente().valeurs();
}

const DoubleTab& Source_Transport_K_Omega_VEF_Face::get_cisaillement_paroi() const
{
  return turbulence_model->loi_paroi().valeur().Cisaillement_paroi();
}

const DoubleTab& Source_Transport_K_Omega_VEF_Face::get_K_pour_production() const
{
  return eqn_K_Omega->inconnue().valeurs();
}

const Nom Source_Transport_K_Omega_VEF_Face::get_type_paroi() const
{
  return turbulence_model->loi_paroi().valeur().que_suis_je();
}

void Source_Transport_K_Omega_VEF_Face::compute_blending_F1(DoubleTab& gradKgradOmega)
{

  const DoubleTab& K_Omega = eqn_K_Omega->inconnue().valeurs();
  const DoubleTab& kinematic_viscosity = get_visc_turb();
  const DoubleTab& distmin = le_dom_VEF->y_faces(); // Minimum distance to the edge
  // DoubleTab& tmpF1 = turbulence_model->get_blenderF1(); // The blending field F1
  DoubleTab& tmpF1 = turbulence_model->get_blenderF1();

  // Loop on faces
  for (int face = 0; face < le_dom_VEF->nb_faces(); face++)
    {
      double const dmin = distmin(face);
      double const enerK = K_Omega(face, 0);
      double const omega = K_Omega(face, 1);

      double const tmp1 = sqrt(enerK)/(BETA_K*omega*dmin);
      double const tmp2 = 500.0*kinematic_viscosity(face)/(omega*dmin*dmin);
      double const tmp3 = 4.0*SIGMA_OMEGA2*enerK/(gradKgradOmega(face)*dmin);
      double const arg1 = std::min(std::max(tmp1, tmp2), tmp3); // Common name of the variable
      tmpF1(face) = std::tanh(arg1*arg1*arg1*arg1);
    }
}

double Source_Transport_K_Omega_VEF_Face::blender(double const val1, double const val2,
                                                  int const face) const
{
  const DoubleTab& F1 = turbulence_model->get_blenderF1();
  return F1(face)*val1 + (1 - F1(face))*val2;
}


void Source_Transport_K_Omega_VEF_Face::compute_cross_diffusion(DoubleTab& gradKgradOmega) const
{
  // Same structure than Source_WC_Chaleur_VEF.cpp in TRUST
  // A mettre dans la base ?

  // = First, store K and Omega in two separated Tab for the gradient computation, not ideal
  const DoubleTab& K_Omega = eqn_K_Omega->inconnue().valeurs();
  DoubleTab enerK; // field on faces
  DoubleTab omega; // field on faces
  enerK.resize(K_Omega.dimension_tot(0));
  omega.resize(K_Omega.dimension_tot(0));
  for (int num_face = 0; num_face < le_dom_VEF->nb_faces(); ++num_face)
    {
      enerK(num_face) = K_Omega(num_face, 0);
      omega(num_face) = K_Omega(num_face, 1);
      gradKgradOmega(num_face) = 0; // mise à zéro pour être certain
    }

  // = Definition of the gradients
  DoubleTab gradK_elem; // field on elements
  DoubleTab gradOmega_elem; // field on elements

  // Get number of componants, for gradient resize
  // We use the velocity field to resize the gradient as the velocity is on faces
  // const auto& eqHyd = mon_equation->probleme().equation(0);
  // const Probleme_base& mon_pb = mon_equation->probleme();
  // const auto& eqHyd = mon_pb.equation();
  const Navier_Stokes_Turbulent& eqHyd = ref_cast(Navier_Stokes_Turbulent,
                                                  ref_cast(Pb_Hydraulique_Turbulent,
                                                           mon_equation->probleme()).equation(0));
  const DoubleTab& velocity_field_face = eqHyd.vitesse().valeurs(); // Velocity on faces
  const int nbr_velocity_components = velocity_field_face.dimension(1);
  const DoubleTab& pressure = eqHyd.pression().valeurs();
  const int total_number_of_faces = pressure.dimension(0); // find a better name than nb_tot
  gradK_elem.resize(total_number_of_faces, nbr_velocity_components);
  gradOmega_elem.resize(total_number_of_faces, nbr_velocity_components);
  // resize_gradient_tab(gradK);
  // resize_gradient_tab(gradOmega);

  // Compute the two gradients
  const Operateur_Grad& Op_Grad_komega = eqn_K_Omega->gradient_operator_komega();
  Op_Grad_komega.calculer(enerK, gradK_elem);
  Op_Grad_komega.calculer(omega, gradOmega_elem);

  // Correction on the boundaries? Elie put the pressure at zero.

  // Interpolate from elem to face
  const Domaine_dis_base& domaine_dis = mon_equation->inconnue().domaine_dis_base();
  const Domaine_VF& domaine = ref_cast(Domaine_VF, domaine_dis);
  DoubleTab gradK_face(velocity_field_face); // gradK on faces is similar to the velocity
  elem_to_face(domaine, gradK_elem, gradK_face);
  DoubleTab gradOmega_face(velocity_field_face);
  elem_to_face(domaine, gradOmega_elem, gradOmega_face);

  // Dot Product gradKgradOmega
  for (int num_face = 0; num_face < le_dom_VEF->nb_faces(); ++num_face)
    for (int ncompo = 0; ncompo < nbr_velocity_components; ++ncompo)
      gradKgradOmega(num_face) +=
        gradK_face(num_face, ncompo) * gradOmega_face(num_face, ncompo);
}

// cAlan: Tried to get a dedicated function to resize a tab. Make a function for this?
// void Source_Transport_K_Omega_VEF_Face::resize_gradient_tab(DoubleTab& grad) const
// {
//   const Navier_Stokes_Turbulent& eqHyd = ref_cast(Navier_Stokes_Turbulent,
//                                                   mon_equation->probleme().equation(0));
//   const DoubleTab& la_vitesse = eqHyd.vitesse().valeurs();
//   const int nb_compo = la_vitesse.line_size();
//   const DoubleTab& pressure = eqHyd.pression().valeurs();
//   const int nb_tot = pressure.size_totale(); // find a better name than nb_tot
//   grad.resize(nb_tot, nb_compo);
// }

void Source_Transport_K_Omega_VEF_Face::fill_resu(const DoubleVect& volumes_entrelaces,
                                                  const DoubleTrav& ProdK,
                                                  const DoubleTab& gradKgradOmega,
                                                  DoubleTab& resu) const
{
  const DoubleTab& K_Omega = eqn_K_Omega->inconnue().valeurs();
  const double LeK_MIN = eqn_K_Omega->modele_turbulence().get_K_MIN();

  for (int face = 0; face < le_dom_VEF->nb_faces(); face++)
    {
      resu(face, 0) += (ProdK(face) - BETA_K*K_Omega(face, 0)*K_Omega(face, 1))*volumes_entrelaces(face);
      if (K_Omega(face, 0) >= LeK_MIN)
        {
          double cALPHA = ALPHA_OMEGA;
          double cBETA = BETA_OMEGA;
          double cSIGMA = (gradKgradOmega(face) > 0) ? 1/8 : 0;
          if (turbulence_model->get_model_variant() == "SST")
            {
              cALPHA = blender(GAMMA1, GAMMA2, face);
              cBETA = blender(BETA1, BETA2, face);
              cSIGMA = 2*(1 - turbulence_model->get_blenderF1()(face)*SIGMA_OMEGA2);
            }

          resu(face, 1) += cALPHA*ProdK(face)*K_Omega(face, 1)/K_Omega(face, 0); // production
          resu(face, 1) += - cBETA*K_Omega(face, 1)*K_Omega(face, 1); // dissipation
          resu(face, 1) += cSIGMA/K_Omega(face, 1)*gradKgradOmega(face); // cross diffusion
          resu(face, 1) *= volumes_entrelaces(face);
        }
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

  // cAlan: to be adapted for k_omega? Copy of the k-eps method
  for (int face = 0; face < K_Omega.dimension(0); face++)
    if (K_Omega(face, 0) >= LeK_MIN)
      {
        const double volporo = porosite_face(face) * volumes_entrelaces(face);

        double coef_k = K_Omega(face, 1)/K_Omega(face, 0)*volporo;
        matrice(face*2, face*2) += coef_k;

        double coef_omega = ALPHA_OMEGA*coef_k;
        matrice(face*2 + 1, face*2 + 1) += coef_omega;
      }
}
