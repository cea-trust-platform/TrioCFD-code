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
// File      : Source_Transport_VEF_K_Omega_Face_base.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VEF/Sources/RANS/Base
//
/////////////////////////////////////////////////////////////////////////////

#include <Pb_Thermohydraulique_Concentration_Turbulent.h>
#include <Pb_Thermohydraulique_Turbulent.h>
#include <Source_Transport_VEF_K_Omega_Face_base.h>
#include <Modele_turbulence_hyd_K_Omega.h>
#include <Champ_Uniforme.h>
#include <Constituant.h>
#include <Fluide_base.h>
#include <Domaine_Cl_VEF.h>
#include <Champ_P1NC.h>
#include <Domaine_VEF.h>

Implemente_base_sans_constructeur(Source_Transport_VEF_K_Omega_Face_base,
				  "Source_Transport_VEF_K_Omega_Face_base",
				  Source_base);

Sortie& Source_Transport_VEF_K_Omega_Face_base::printOn(Sortie& os) const { return os << que_suis_je(); }

Entree& Source_Transport_VEF_K_Omega_Face_base::readOn(Entree& is) { return Source_Transport_proto::readOn_proto(is, que_suis_je()); }

void Source_Transport_VEF_K_Omega_Face_base::associer_domaines(const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis)
{
  le_dom_VEF = ref_cast(Domaine_VEF, domaine_dis.valeur());
  le_dom_Cl_VEF = ref_cast(Domaine_Cl_VEF, domaine_Cl_dis.valeur());
}

void Source_Transport_VEF_K_Omega_Face_base::associer_pb(const Probleme_base& pb) { Source_Transport_proto::associer_pb_proto(pb); }

DoubleTab& Source_Transport_VEF_K_Omega_Face_base::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}

DoubleTab& Source_Transport_VEF_K_Omega_Face_base::ajouter_komega(DoubleTab& resu) const
{
  const Domaine_Cl_VEF& domaine_Cl_VEF = ref_cast(Domaine_Cl_VEF,
						  eq_hydraulique->domaine_Cl_dis().valeur());
  const DoubleTab& visco_turb = get_visc_turb(); // voir les classes filles
  const DoubleTab& vit = eq_hydraulique->inconnue().valeurs();
  const DoubleVect& volumes_entrelaces = le_dom_VEF->volumes_entrelaces();
  const DoubleTab& tab = get_cisaillement_paroi(); // voir les classes filles
  const int nb_faces_ = le_dom_VEF->nb_faces();
  DoubleTrav P{nb_faces_};
  const DoubleTab& K = get_K_pour_production(); // voir les classes filles
  calculer_terme_production_K(le_dom_VEF.valeur(), domaine_Cl_VEF, P, K,
			      vit, visco_turb, _interpolation_viscosite_turbulente);

  fill_resu(volumes_entrelaces, P, resu); // voir les classes filles

  return resu;
}


