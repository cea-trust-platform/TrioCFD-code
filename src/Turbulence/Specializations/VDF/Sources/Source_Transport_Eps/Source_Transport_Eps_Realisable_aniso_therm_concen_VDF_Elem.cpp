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
// File      : Source_Transport_Eps_Realisable_aniso_therm_concen_VDF_Elem.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VDF/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_Eps_Realisable_aniso_therm_concen_VDF_Elem.h>
#include <Modele_turbulence_hyd_K_Eps_Realisable_Bicephale.h>

Implemente_instanciable_sans_constructeur(Source_Transport_Eps_Realisable_aniso_therm_concen_VDF_Elem,"Source_Transport_Eps_Realisable_aniso_therm_concen_VDF_P0_VDF",Source_Transport_Eps_Realisable_VDF_Elem);

Sortie& Source_Transport_Eps_Realisable_aniso_therm_concen_VDF_Elem::printOn(Sortie& s) const { return s << que_suis_je(); }
Entree& Source_Transport_Eps_Realisable_aniso_therm_concen_VDF_Elem::readOn(Entree& is) { return Source_Transport_Eps_Realisable_VDF_Elem::readOn_anisotherme_concen_real(is,que_suis_je()); }

void Source_Transport_Eps_Realisable_aniso_therm_concen_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  Source_Transport_Eps_Realisable_VDF_Elem::verifier_milieu_anisotherme_concen(pb,que_suis_je());
  Source_Transport_Eps_Realisable_VDF_Elem::associer_pb(pb);
  Source_Transport_Eps_Realisable_VDF_Elem::associer_pb_anisotherme_concen(pb);
}

void Source_Transport_Eps_Realisable_aniso_therm_concen_VDF_Elem::fill_resu_anisotherme_concen(const DoubleVect& G_t, const DoubleVect& G_c, const DoubleVect& volumes, const DoubleVect& porosite_vol, DoubleTab& resu) const
{
  const DoubleTab& K_Rea = eqn_k_Rea->inconnue().valeurs(), &eps_Rea = eqn_eps_Rea->inconnue().valeurs();
  // C1 value is not a constant in Realizable K-Epsilon model but here, we take the default value of C1 used in standard K-Epsilon, as proposed by litterature
  double C3_loc, G_sum, C1_loc = C1__, LeK_MIN = eqn_k_Rea->modele_turbulence().get_K_MIN();

  for (int elem = 0; elem < le_dom_VDF->nb_elem(); elem++)
    {
      G_sum = G_t(elem) + G_c(elem) ;
      if (K_Rea(elem) >= LeK_MIN)
        {
          C3_loc = G_sum > 0. ? 0. : C3;
          resu(elem) += C1_loc*(1.-C3_loc)*G_sum *volumes(elem)*porosite_vol(elem)*eps_Rea(elem)/K_Rea(elem);
        }
    }
}

void Source_Transport_Eps_Realisable_aniso_therm_concen_VDF_Elem::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  Source_Transport_Eps_Realisable_VDF_Elem::ajouter_blocs(matrices, secmem, semi_impl);
  Source_Transport_Eps_Realisable_VDF_Elem::ajouter_anisotherme_concen(secmem);
}
