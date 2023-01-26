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
// File      : Source_Transport_K_Omega_VDF_Elem_base.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VDF/Sources/new
//
/////////////////////////////////////////////////////////////////////////////

#include <Pb_Hydraulique_Concentration_Turbulent.h>
#include <Pb_Thermohydraulique_Turbulent.h>
#include <Source_Transport_K_Omega_VDF_Elem_base.h>
#include <Modele_turbulence_hyd_K_Omega.h>
#include <Champ_Uniforme.h>
#include <Fluide_base.h>
#include <Zone_Cl_VDF.h>
#include <Constituant.h>
#include <Champ_Face_VDF.h>
#include <Debog.h>

Implemente_base_sans_constructeur(Source_Transport_K_Omega_VDF_Elem_base,
                                  "Source_Transport_K_Omega_VDF_Elem_base",
                                  Source_base);

Sortie& Source_Transport_K_Omega_VDF_Elem_base::printOn( Sortie& os ) const { return os << que_suis_je(); }
Entree& Source_Transport_K_Omega_VDF_Elem_base::readOn( Entree& is ) { return Source_Transport_proto::readOn_proto(is, que_suis_je()); }

// cAlan : mutualisable
void Source_Transport_K_Omega_VDF_Elem_base::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis&  zone_Cl_dis)
{
  la_zone_VDF = ref_cast(Zone_VDF, zone_dis.valeur());
  la_zone_Cl_VDF = ref_cast(Zone_Cl_VDF,zone_Cl_dis.valeur());
}

void Source_Transport_K_Omega_VDF_Elem_base::associer_pb(const Probleme_base& pb) { Source_Transport_proto::associer_pb_proto(pb); }

DoubleTab& Source_Transport_K_Omega_VDF_Elem_base::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}

DoubleTab& Source_Transport_K_Omega_VDF_Elem_base::ajouter_komega(DoubleTab& resu) const
{
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const DoubleTab& visco_turb = get_visc_turb(); // voir les classes filles
  const DoubleTab& vit = eq_hydraulique->inconnue().valeurs();
  const Champ_Face_VDF& ch_vit = ref_cast(Champ_Face_VDF, eq_hydraulique->inconnue().valeur());

  DoubleVect P; // Ajout d'un espace virtuel au tableau P
  zone_VDF.zone().creer_tableau_elements(P);
  calculer_terme_production(ch_vit, visco_turb, vit, P); // voir les classes filles
  fill_resu(P, resu);
  resu.echange_espace_virtuel();
  return resu;
}
