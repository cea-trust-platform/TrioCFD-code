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
// File:        Y_plus_Champ_Face.cpp
// Directory:   $TRUST_ROOT/src/VDF/Champs
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Y_plus_Champ_Face.h>
#include <Zone_Cl_VDF.h>
#include <Champ_Face.h>
#include <Equation_base.h>
#include <Milieu_base.h>
#include <Mod_turb_hyd_base.h>
Implemente_instanciable(Y_plus_Champ_Face,"Y_plus_Champ_Face",Champ_Fonc_P0_VDF);


//     printOn()
/////

Sortie& Y_plus_Champ_Face::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

//// readOn
//

Entree& Y_plus_Champ_Face::readOn(Entree& s)
{
  return s ;
}

void Y_plus_Champ_Face::associer_champ(const Champ_Face& un_champ)
{
  mon_champ_= un_champ;
}

void Y_plus_Champ_Face::me_calculer(double tps)
{
  const Milieu_base& mil = mon_champ ().equation().milieu(); // returns Fluide_Diphasique or Fluide_Incompressible
  const RefObjU& modele_turbulence = mon_champ ().equation().get_modele(TURBULENCE);
  const Mod_turb_hyd_base& mod_turb = ref_cast(Mod_turb_hyd_base,modele_turbulence.valeur());
  const Turbulence_paroi_base& loipar = mod_turb.loi_paroi();
  const Nom& nom_loipar = loipar.que_suis_je();
  const Nom& nom_eq = mon_champ ().equation().que_suis_je();
  const Nom& nom_mil = mil.que_suis_je();

  if ( nom_eq == "Navier_Stokes_FT_Disc" && nom_mil == "Fluide_Diphasique" && nom_loipar =="loi_standard_hydr_diphasique_VDF")
    mon_champ_->calcul_y_plus_diphasique(valeurs(), la_zone_Cl_VDF.valeur());
  else
    mon_champ_->calcul_y_plus(valeurs(), la_zone_Cl_VDF.valeur());
}

const Zone_Cl_dis_base& Y_plus_Champ_Face::zone_Cl_dis_base() const
{
  return la_zone_Cl_VDF.valeur();
}


