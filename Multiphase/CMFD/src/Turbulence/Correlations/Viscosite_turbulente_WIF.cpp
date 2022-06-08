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
// File:        Viscosite_turbulente_WIF.cpp
// Directory:   $TRUST_ROOT/src/Turbulence/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Viscosite_turbulente_WIF.h>
#include <Param.h>
#include <Probleme_base.h>
#include <Champ_base.h>
#include <TRUSTTab_parts.h>
#include <TRUSTTrav.h>

Implemente_instanciable(Viscosite_turbulente_WIF, "Viscosite_turbulente_WIF", Viscosite_turbulente_base);

Sortie& Viscosite_turbulente_WIF::printOn(Sortie& os) const
{
  return os;
}

Entree& Viscosite_turbulente_WIF::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("limiter|limiteur", &limiter_);
  param.ajouter("rapport_aspect_moyen|average_aspect_ratio", &gamma_);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Viscosite_turbulente_WIF::eddy_viscosity(DoubleTab& nu_t) const
{
  nu_t = 0; // no eddy viscosity here ; use ggdh for diffusion
}

void Viscosite_turbulente_WIF::reynolds_stress(DoubleTab& R_ij) const // Renvoie <u_i'u_j'>
{
  const DoubleTab& tab_alpha = pb_->get_champ("alpha").passe(), &tab_u = pb_->get_champ("vitesse").passe();
  ConstDoubleTab_parts p_u(tab_u); //en PolyMAC_P0, tab_u contient (nf.u) aux faces, puis (u_i) aux elements

  int i, d, db, D = dimension, i_part = -1, N = tab_alpha.dimension(1);
  if (N!=2) Process::exit("Visocisty_turbulente_WIF is only coded for 2 phases");
  if (D!=3) Process::exit("Visocisty_turbulente_WIF is only coded for 3 dimensions");

  for (i = 0; i < p_u.size(); i++)
    if (p_u[i].get_md_vector() == R_ij.get_md_vector()) i_part = i; //on cherche une partie ayant le meme support
  if (i_part < 0) Process::exit("Viscosite_turbulente_WIF : inconsistency between velocity and Rij!");
  const DoubleTab& u = p_u[i_part]; //le bon tableau
  DoubleTrav u_r(u.dimension(0), D);
  DoubleTrav u_r_carre(u.dimension(0), 1);
  for (i = 0; i < R_ij.dimension(0); i++)
    for (d = 0; d < D; d++)
      {
        u_r(i, d) = u(i, d, 1) - u(i, d, 0); // relative speed = gas speed - liquid speed
        u_r_carre(i,0) += u_r(i, d)*u_r(i, d);
      }

  for (i = 0; i < R_ij.dimension(0); i++)
    for (d = 0; d < D; d++)
      for (db = 0; db < D; db++)
        {
          R_ij(i, 1, d, db) = 0 ; // No WIF for gas phase
          R_ij(i, 0, d, db) = tab_alpha(i, 0) * ( 3/20 * u_r_carre(i,0) * (db==d?1:0) + (1/20 + 1.5 * gamma_*gamma_*gamma_*0.25) * u_r(i, d)* u_r(i, db) ) ;
        }
}

void Viscosite_turbulente_WIF::k_over_eps(DoubleTab& k_sur_eps) const
{
  k_sur_eps =  0;
}

void Viscosite_turbulente_WIF::eps(DoubleTab& eps) const
{
  eps =  0;
}
