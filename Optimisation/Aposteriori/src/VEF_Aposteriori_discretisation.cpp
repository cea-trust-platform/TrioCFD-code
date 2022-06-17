/****************************************************************************
* Copyright (c) 2021, CEA
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

#include <VEF_Aposteriori_discretisation.h>
#include <estimateur_aposteriori_P0_VEF.h>
#include <Champ_P1_isoP1Bulle.h>
#include <Zone_Cl_dis.h>
#include <Zone_Cl_VEF.h>
#include <Zone_VEF.h>
#include <Champ_P0_VEF.h>
#include <Champ_Fonc.h>
#include <Champ_P1NC.h>

Implemente_instanciable(VEF_Aposteriori_discretisation,"VEF_Aposteriori_discretisation",VEF_discretisation);

Entree& VEF_Aposteriori_discretisation::readOn(Entree& s) { return s ; }

Sortie& VEF_Aposteriori_discretisation::printOn(Sortie& s) const { return s ; }

void VEF_Aposteriori_discretisation::estimateur_aposteriori(const Zone_dis& z, const Zone_Cl_dis& zcl, const Champ_Inc& ch_vitesse, const Champ_Inc& ch_pression, const Champ_Don& viscosite_cinematique, Champ_Fonc& champ) const
{
  const Zone_VEF_PreP1b& zone_vef=ref_cast(Zone_VEF_PreP1b, z.valeur());
  const Zone_Cl_VEF& zone_cl_vef=ref_cast(Zone_Cl_VEF, zcl.valeur());
  champ.typer("estimateur_aposteriori_P0_VEF");
  estimateur_aposteriori_P0_VEF& ch=ref_cast(estimateur_aposteriori_P0_VEF,champ.valeur());
  ch.associer_zone_dis_base(zone_vef);
  const Champ_P1NC& vit = ref_cast(Champ_P1NC, ch_vitesse.valeur());
  const Champ_P1_isoP1Bulle& pres = ref_cast(Champ_P1_isoP1Bulle, ch_pression.valeur());
  ch.associer_champ(vit, pres, viscosite_cinematique, zone_cl_vef);
  ch.nommer("estimateur_aposteriori");
  ch.fixer_nb_comp(1);
  ch.fixer_nb_valeurs_nodales(zone_vef.nb_elem());
  ch.fixer_unite("sans");
  ch.changer_temps(ch_vitesse.temps());
}
