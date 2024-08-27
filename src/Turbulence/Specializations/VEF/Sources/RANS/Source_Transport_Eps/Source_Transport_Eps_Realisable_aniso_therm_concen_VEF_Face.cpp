/****************************************************************************
* Copyright (c) 2022, CEA
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
// File      : Source_Transport_Eps_Realisable_aniso_therm_concen_VEF_Face.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VEF/Sources/RANS
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_Eps_Realisable_aniso_therm_concen_VEF_Face.h>
#include <Modele_turbulence_hyd_K_Eps_Realisable_Bicephale.h>
#include <Domaine_VEF.h>

Implemente_instanciable_sans_constructeur(Source_Transport_Eps_Realisable_aniso_therm_concen_VEF_Face,"Source_Transport_Eps_Realisable_aniso_therm_concen_VEF_P1NC",Source_Transport_Eps_Realisable_VEF_Face);

Sortie& Source_Transport_Eps_Realisable_aniso_therm_concen_VEF_Face::printOn(Sortie& s) const { return s << que_suis_je() ; }
Entree& Source_Transport_Eps_Realisable_aniso_therm_concen_VEF_Face::readOn(Entree& is) { return Source_Transport_Eps_Realisable_VEF_Face::readOn_anisotherme_concen_real(is,que_suis_je()); }

void Source_Transport_Eps_Realisable_aniso_therm_concen_VEF_Face::associer_pb(const Probleme_base& pb)
{
  Source_Transport_Eps_Realisable_VEF_Face::verifier_milieu_anisotherme_concen(pb,que_suis_je());
  Source_Transport_Eps_Realisable_VEF_Face::associer_pb(pb);
  Source_Transport_Eps_Realisable_VEF_Face::associer_pb_anisotherme_concen(pb);
}

DoubleTab& Source_Transport_Eps_Realisable_aniso_therm_concen_VEF_Face::ajouter(DoubleTab& resu) const
{
  Source_Transport_Eps_Realisable_VEF_Face::ajouter(resu);
  return Source_Transport_Eps_Realisable_VEF_Face::ajouter_anisotherme_concen(resu);
}

void Source_Transport_Eps_Realisable_aniso_therm_concen_VEF_Face::fill_resu_anisotherme_concen(const DoubleTrav& G_t, const DoubleTrav& G_c, const DoubleVect& volumes_entrelaces, DoubleTab& resu) const
{
  const DoubleTab& K = eqn_k_Rea->inconnue()->valeurs(), &eps = eqn_eps_Rea->inconnue()->valeurs();
  // C1 value is not a constant in Realizable K-Epsilon model but here, we take the default value of C1 used in standard K-Epsilon, as proposed by litterature
  double C3_loc, G_sum, C1_loc = C1__, LeK_MIN = eqn_k_Rea->modele_turbulence().get_K_MIN();
  for (int face = 0; face < le_dom_VEF->nb_faces(); face++)
    {
      G_sum = G_t(face) + G_c(face);
      if (K(face) >= LeK_MIN)
        {
          C3_loc = G_sum > 0. ? 0. : C3;
          resu(face) += C1_loc * (1 - C3_loc) * G_sum * volumes_entrelaces(face) * eps(face) / K(face);
        }
    }
}
