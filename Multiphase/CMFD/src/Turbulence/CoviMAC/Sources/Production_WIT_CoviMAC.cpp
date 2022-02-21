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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Production_WIT_CoviMAC.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/CoviMAC
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Production_WIT_CoviMAC.h>

#include <Zone_CoviMAC.h>
#include <Champ_P0_CoviMAC.h>
#include <Matrix_tools.h>
#include <Probleme_base.h>
#include <grad_Champ_Face_CoviMAC.h>
#include <Champ_Uniforme.h>
#include <Flux_interfacial_base.h>
#include <Milieu_composite.h>
#include <Operateur_Diff.h>
#include <Op_Diff_Turbulent_CoviMAC_Face.h>
#include <Navier_Stokes_std.h>
#include <Viscosite_turbulente_base.h>
#include <Viscosite_turbulente_multiple.h>
#include <ConstDoubleTab_parts.h>


Implemente_instanciable(Production_WIT_CoviMAC,"Production_WIT_P0_CoviMAC", Source_base);

Sortie& Production_WIT_CoviMAC::printOn(Sortie& os) const
{
  return os;
}

Entree& Production_WIT_CoviMAC::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("Reb_c", &Reb_c_);
  param.ajouter("g", &g_);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Production_WIT_CoviMAC::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
// empty : no derivative to add in the blocks
}

void Production_WIT_CoviMAC::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Zone_CoviMAC&                      zone = ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
  const DoubleTab&                      tab_rho = equation().probleme().get_champ("masse_volumique").passe();
  const DoubleTab&                      tab_alp = equation().probleme().get_champ("alpha").passe();
  const DoubleTab&                          vit = equation().probleme().get_champ("vitesse").passe();
  const DoubleTab&                         diam = equation().probleme().get_champ("diametre_bulles").passe();
  const DoubleTab&                           nu = equation().probleme().get_champ("viscosite_cinematique").passe();

  const DoubleVect& pe = zone.porosite_elem(), &ve = zone.volumes();

  int N = equation().inconnue().valeurs().dimension(1), ne = zone.nb_elem(), D = dimension, i_part=-1 ;
  if (N!=1) Process::exit("WIT is only in the liquid phase");
  if (D!=3) Process::exit("WIT is only coded for 3 dimensions");

  ConstDoubleTab_parts p_u(vit); //en CoviMAC, vit contient (nf.u) aux faces, puis (u_i) aux elements
  for (int i = 0; i < p_u.size(); i++) if (p_u[i].get_md_vector() == nu.get_md_vector()) i_part = i; //on cherche une partie ayant le meme support
  if (i_part < 0) Process::exit("Viscosite_turbulente_WIF : inconsistency between velocity and viscosity!");
  const DoubleTab& u = p_u[i_part]; //le bon tableau
  DoubleTrav u_r(u.dimension(0), 1);
  for (int i = 0; i < u_r.dimension(0); i++)
    {
      for (int d = 0; d < D; d++) u_r(i, 0) += (u(i, d, 1) - u(i, d, 0))*(u(i, d, 1) - u(i, d, 0)); // relative speed = gas speed - liquid speed
      u_r(i, 0) = std::sqrt(u_r(i, 0));
    }

  // On calcule le second membre aux elements (implicite uniquement pour le moment)
  for(int e = 0 ; e < ne ; e++)
    {
      double Reb = diam(e,0)*u_r(e, 0)/nu(e,0);
      secmem(e, 0) += ve(e) * pe(e) * tab_alp(e, 0)* tab_rho(e, 0) * tab_alp(e, 1) *(tab_rho(e, 0)-tab_rho(e, 1))/tab_rho(e, 0)*g_*u_r(e,0) *(0.9 - exp(-Reb/Reb_c_));
    }
}
