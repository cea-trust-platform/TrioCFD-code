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
// File:        LoiParoiHybride_VDF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Lois_Paroi/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <LoiParoiHybride_VDF.h>
#include <Dirichlet_paroi_fixe.h>
#include <Dirichlet_paroi_defilante.h>
#include <Neumann_paroi.h>
#include <Domaine_Cl_VDF.h>

Implemente_instanciable(LoiParoiHybride_VDF,"Loi_Paroi_Hybride_VDF",Paroi_hyd_base_VDF);

Sortie& LoiParoiHybride_VDF::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

Entree& LoiParoiHybride_VDF::readOn(Entree& s)
{
  const int nbord = le_dom_VDF->nb_front_Cl();
  Noms noms(nbord);
  for (int ibord=0; ibord<nbord; ibord++)
    {
      const Cond_lim& la_cl = le_dom_Cl_VDF->les_conditions_limites(ibord);
      noms[ibord] = la_cl->frontiere_dis().le_nom();
    }

  LoiParoiHybride::lire(s,noms, mon_modele_turb_hyd.valeur());
  return s ;
}

void LoiParoiHybride_VDF::associer(const Domaine_dis_base& zd, const Domaine_Cl_dis_base& zcl)
{
  LoiParoiHybride::associer(zd,zcl);
  Paroi_hyd_base_VDF::associer(zd,zcl);
}


int LoiParoiHybride_VDF::init_lois_paroi()
{
  LoiParoiHybride::init_lois_paroi();
  if (!Cisaillement_paroi_.get_md_vector().non_nul())
    {
      Cisaillement_paroi_.resize(0, dimension);
      le_dom_VDF->creer_tableau_faces_bord(Cisaillement_paroi_);
    }
  return 1;
}

int LoiParoiHybride_VDF::calculer_hyd(DoubleTab& tab1)
{
  LoiParoiHybride::calculer_hyd(tab1, le_dom_VDF.valeur(), le_dom_Cl_VDF.valeur(), Cisaillement_paroi_);

  return 1;
}

int LoiParoiHybride_VDF::calculer_hyd(DoubleTab& tab1, DoubleTab& tab2)
{
  LoiParoiHybride::calculer_hyd(tab1, tab2, le_dom_VDF.valeur(), le_dom_Cl_VDF.valeur(), Cisaillement_paroi_);

  return 1;
}

int LoiParoiHybride_VDF::calculer_hyd_BiK(DoubleTab& tab_k, DoubleTab& tab_eps)
{
  LoiParoiHybride::calculer_hyd_BiK(tab_k, tab_eps, le_dom_VDF.valeur(), le_dom_Cl_VDF.valeur(), Cisaillement_paroi_);

  return 1;
}
