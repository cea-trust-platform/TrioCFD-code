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
// File:        LoiParoiHybride_VEF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Lois_Paroi/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <LoiParoiHybride_VEF.h>
#include <Domaine_Cl_VEF.h>
#include <Domaine_VEF.h>
#include <Domaine_dis.h>
#include <Domaine_Cl_dis.h>
#include <LoiParoiHybride.h>

Implemente_instanciable(LoiParoiHybride_VEF,"Loi_Paroi_Hybride_VEF",Paroi_hyd_base_VEF);

Sortie& LoiParoiHybride_VEF::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

Entree& LoiParoiHybride_VEF::readOn(Entree& s)
{
  const int nbord = le_dom_VEF.valeur().nb_front_Cl();
  Noms noms(nbord);
  for (int ibord=0; ibord<nbord; ibord++)
    {
      const Cond_lim& la_cl = le_dom_Cl_VEF->les_conditions_limites(ibord);
      noms[ibord] = la_cl.valeur().frontiere_dis().le_nom();
    }

  LoiParoiHybride::lire(s,noms, mon_modele_turb_hyd.valeur());
  return s ;
}

void LoiParoiHybride_VEF::associer(const Domaine_dis& zd, const Domaine_Cl_dis& zcl)
{
  LoiParoiHybride::associer(zd,zcl);
  le_dom_VEF = ref_cast(Domaine_VEF,zd.valeur());
  le_dom_Cl_VEF = ref_cast(Domaine_Cl_VEF,zcl.valeur());
}


int LoiParoiHybride_VEF::init_lois_paroi()
{
  LoiParoiHybride::init_lois_paroi();

  Cisaillement_paroi_.resize(le_dom_VEF->nb_faces_bord(),dimension);
  return 1;
}

int LoiParoiHybride_VEF::calculer_hyd(DoubleTab& tab1)
{
  LoiParoiHybride::calculer_hyd(tab1, le_dom_VEF.valeur(), le_dom_Cl_VEF.valeur(), Cisaillement_paroi_);

  return 1;
}

int LoiParoiHybride_VEF::calculer_hyd(DoubleTab& tab1, DoubleTab& tab2)
{
  LoiParoiHybride::calculer_hyd(tab1, tab2, le_dom_VEF.valeur(), le_dom_Cl_VEF.valeur(), Cisaillement_paroi_);

  return 1;
}

int LoiParoiHybride_VEF::calculer_hyd_BiK(DoubleTab& tab_k,DoubleTab& tab_eps)
{
  LoiParoiHybride::calculer_hyd_BiK(tab_k, tab_eps, le_dom_VEF.valeur(), le_dom_Cl_VEF.valeur(), Cisaillement_paroi_);

  return 1;
}


