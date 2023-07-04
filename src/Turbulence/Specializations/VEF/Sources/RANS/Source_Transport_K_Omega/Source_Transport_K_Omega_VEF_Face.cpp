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

#include <TRUSTTabs_forward.h>
#include <Source_Transport_K_Omega_VEF_Face.h>
#include <Modele_turbulence_hyd_K_Omega.h>
#include <Transport_K_Omega.h>
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

void Source_Transport_K_Omega_VEF_Face::fill_resu(const DoubleVect& volumes_entrelaces,
                                                  const DoubleTrav& P, DoubleTab& resu) const
{
  const DoubleTab& K_Omega = mon_eq_transport_K_Omega->inconnue().valeurs();
  const double LeK_MIN = mon_eq_transport_K_Omega->modele_turbulence().get_K_MIN();
  for (int fac = 0; fac < le_dom_VEF->nb_faces(); fac++)
    {
      resu(fac, 0) += (P(fac) - BETA_K*K_Omega(fac, 1))*volumes_entrelaces(fac);
      if (K_Omega(fac, 0) >= LeK_MIN)
        resu(fac, 1) += (ALPHA_OMEGA*P(fac)*K_Omega(fac, 1)/K_Omega(fac, 0)
                         - BETA_OMEGA*K_Omega(fac, 1)*K_Omega(fac, 1))*volumes_entrelaces(fac);
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
