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
// File      : Source_Transport_K_Eps_aniso_therm_concen_VDF_Elem.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VDF/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_K_Eps_aniso_therm_concen_VDF_Elem.h>
#include <Transport_K_Eps.h>
#include <Zone_VDF.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Eps_aniso_therm_concen_VDF_Elem,"Source_Transport_K_Eps_aniso_therm_concen_VDF_P0_VDF",Source_Transport_K_Eps_VDF_Elem);

Sortie& Source_Transport_K_Eps_aniso_therm_concen_VDF_Elem::printOn(Sortie& s) const { return s << que_suis_je(); }
Entree& Source_Transport_K_Eps_aniso_therm_concen_VDF_Elem::readOn(Entree& is)
{
  Source_Transport_K_Eps_VDF_Elem::verifier_pb_keps_anisotherme_concen(mon_equation->probleme(),que_suis_je());
  return Source_Transport_K_Eps_VDF_Elem::readOn_anisotherme_concen(is);
}

void Source_Transport_K_Eps_aniso_therm_concen_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  Source_Transport_K_Eps_VDF_Elem::verifier_milieu_anisotherme_concen(pb,que_suis_je());
  Source_Transport_K_Eps_VDF_Elem::associer_pb(pb);
  Source_Transport_K_Eps_VDF_Elem::associer_pb_anisotherme_concen(pb);
}

void Source_Transport_K_Eps_aniso_therm_concen_VDF_Elem::fill_resu_anisotherme_concen(const DoubleVect& G_t, const DoubleVect& G_c, const DoubleVect& volumes, const DoubleVect& porosite_vol, DoubleTab& resu) const
{
  const DoubleTab& K_eps = mon_eq_transport_K_Eps->inconnue().valeurs();
  double C3_loc, G_sum, LeK_MIN = mon_eq_transport_K_Eps->modele_turbulence().get_LeK_MIN();

  for (int elem = 0; elem < la_zone_VDF->nb_elem(); elem++)
    {
      G_sum = G_t(elem) + G_c(elem) ;
      resu(elem,0) += G_sum *volumes(elem)*porosite_vol(elem);
      if (K_eps(elem,0) >= LeK_MIN)
        {
          C3_loc = G_sum > 0. ? 0. : C3;
          resu(elem,1) += (1.-C3_loc)*G_sum *volumes(elem)*porosite_vol(elem)*K_eps(elem,1)/K_eps(elem,0);
        }
    }
}

DoubleTab& Source_Transport_K_Eps_aniso_therm_concen_VDF_Elem::ajouter(DoubleTab& resu) const
{
  Source_Transport_K_Eps_VDF_Elem::ajouter(resu);
  return Source_Transport_K_Eps_VDF_Elem::ajouter_anisotherme_concen(resu);
}
