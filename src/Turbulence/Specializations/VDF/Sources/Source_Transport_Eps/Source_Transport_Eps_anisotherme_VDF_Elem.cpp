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
// File      : Source_Transport_Eps_anisotherme_VDF_Elem.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VDF/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_Eps_anisotherme_VDF_Elem.h>
#include <Modele_turbulence_hyd_K_Eps_Bicephale.h>

Implemente_instanciable_sans_constructeur(Source_Transport_Eps_anisotherme_VDF_Elem,"Source_Transport_Eps_anisotherme_VDF_P0_VDF",Source_Transport_Eps_VDF_Elem);

Sortie& Source_Transport_Eps_anisotherme_VDF_Elem::printOn(Sortie& s) const { return s << que_suis_je() ; }
Entree& Source_Transport_Eps_anisotherme_VDF_Elem::readOn(Entree& is)
{
  Source_Transport_Eps_VDF_Elem::verifier_pb_keps_anisotherme(mon_equation->probleme(),que_suis_je());
  return Source_Transport_VDF_Elem_base::readOn_anisotherme(is,que_suis_je());
}

void Source_Transport_Eps_anisotherme_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  Source_Transport_Eps_VDF_Elem::verifier_milieu_anisotherme(pb,que_suis_je());
  Source_Transport_Eps_VDF_Elem::associer_pb(pb);
  Source_Transport_Eps_VDF_Elem::associer_pb_anisotherme(pb);
}

void Source_Transport_Eps_anisotherme_VDF_Elem::fill_resu_anisotherme(const DoubleVect& G, const DoubleVect& volumes, const DoubleVect& porosite_vol, DoubleTab& resu) const
{
  const DoubleTab& K = mon_eq_transport_K->inconnue().valeurs(), &Eps = mon_eq_transport_Eps->inconnue().valeurs();
  double C3_loc, LeK_MIN = mon_eq_transport_Eps->modele_turbulence().get_K_MIN() ;
  for (int elem = 0; elem < le_dom_VDF->nb_elem(); elem++)
    {
      if (K(elem) >= LeK_MIN)
        {
          C3_loc = G(elem) > 0. ? 0. : C3;
          resu(elem) += (1-C3_loc)*G(elem)*volumes(elem)*porosite_vol(elem)*Eps(elem)/K(elem);
        }
    }
}

void Source_Transport_Eps_anisotherme_VDF_Elem::ajouter_blocs(matrices_t matrices, DoubleTab& resu, const tabs_t& semi_impl) const
{
  Source_Transport_Eps_VDF_Elem::ajouter_blocs(matrices, resu, semi_impl); // VB : plutot que de calculer P on appelle ajouter de la classe mere
  Source_Transport_Eps_VDF_Elem::ajouter_anisotherme(resu);
}
