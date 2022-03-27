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
// File      : Source_Transport_Realisable_VDF_Elem_base.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VDF/Sources/new
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_Realisable_VDF_Elem_base.h>
#include <Modele_Fonc_Realisable_base.h>
#include <Champ_Uniforme.h>
#include <Fluide_base.h>
#include <TRUSTTrav.h>
#include <Champ_Face.h>
#include <Param.h>

Implemente_base_sans_constructeur(Source_Transport_Realisable_VDF_Elem_base, "Source_Transport_Realisable_VDF_Elem_base",Source_Transport_VDF_Elem_base) ;

Sortie& Source_Transport_Realisable_VDF_Elem_base::printOn(Sortie& os) const { return os << que_suis_je(); }

Entree& Source_Transport_Realisable_VDF_Elem_base::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("C2_eps", &C2);
  param.lire_avec_accolades(is);
  Cerr << "C2_eps = " << C2 << finl;
  return is;
}

Entree& Source_Transport_Realisable_VDF_Elem_base::readOn_anisotherme_real(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("C2_eps", &C2);
  param.ajouter("C3_eps", &C3);
  param.lire_avec_accolades(is);
  Cerr << "C2_eps = " << C2 << finl;
  Cerr << "C3_eps = " << C3 << finl;
  return is;
}

Entree& Source_Transport_Realisable_VDF_Elem_base::readOn_concen_real(Entree& is) { return readOn_anisotherme_real(is); }

Entree& Source_Transport_Realisable_VDF_Elem_base::readOn_anisotherme_concen_real(Entree& is) { return readOn_anisotherme_real(is); }

DoubleTab& Source_Transport_Realisable_VDF_Elem_base::ajouter_keps_real(DoubleTab& resu) const
{
  const DoubleTab& visco_turb = get_visc_turb(); // voir les classes filles
  const Modele_Fonc_Realisable_base& mon_modele_fonc = get_modele_fonc(); // voir les classes filles
  const Fluide_base& fluide = ref_cast(Fluide_base,eq_hydraulique->milieu());
  const DoubleTab& vit = eq_hydraulique->inconnue().valeurs();
  Champ_Face& vitesse = ref_cast_non_const(Champ_Face,eq_hydraulique->inconnue().valeur());
  const Champ_Don& ch_visco_cin = fluide.viscosite_cinematique();
  const DoubleTab& tab_visco = ch_visco_cin->valeurs();
  const int is_visco_const = sub_type(Champ_Uniforme,ch_visco_cin.valeur());

  double visco = -1.;
  if (is_visco_const) visco = std::max(tab_visco(0,0),DMINFLOAT);

  const int nb_elem = la_zone_VDF->nb_elem();
  DoubleTrav P(visco_turb), CC1(nb_elem), S(nb_elem);
  CC1 = mon_modele_fonc.get_C1();
  S = mon_modele_fonc.get_S();

  calculer_terme_production_real(vitesse,visco_turb,vit,P); // voir les classes filles
  P.echange_espace_virtuel();
  fill_resu_real(is_visco_const,tab_visco,P,CC1,S,visco,resu); // voir les classes filles
  return resu;
}

void Source_Transport_Realisable_VDF_Elem_base::contribuer_a_avec(const DoubleTab& a,  Matrice_Morse& matrice) const
{
  const Fluide_base& fluide = ref_cast(Fluide_base,eq_hydraulique->milieu());
  const Champ_Don& ch_visco_cin = fluide.viscosite_cinematique();
  const DoubleTab& tab_visco = ch_visco_cin->valeurs();
  const DoubleVect& porosite = la_zone_VDF->porosite_elem(), &volumes=la_zone_VDF->volumes();
  const int is_visco_const = sub_type(Champ_Uniforme,ch_visco_cin.valeur());
  double visco = -1.;
  if (is_visco_const) visco = std::max(tab_visco(0,0),DMINFLOAT);
  fill_coeff_matrice(is_visco_const,tab_visco,volumes,porosite,visco,matrice); // voir les classes filles
}
