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
// File:        Viscosite_turbulente_WIT.cpp
// Directory:   $TRUST_ROOT/src/Turbulence/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Viscosite_turbulente_WIT.h>
#include <Param.h>
#include <Probleme_base.h>
#include <Champ_base.h>
#include <TRUSTTab_parts.h>

Implemente_instanciable(Viscosite_turbulente_WIT, "Viscosite_turbulente_WIT", Viscosite_turbulente_base);

Sortie& Viscosite_turbulente_WIT::printOn(Sortie& os) const
{
  return os;
}

Entree& Viscosite_turbulente_WIT::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("limiter|limiteur", &limiter_);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Viscosite_turbulente_WIT::eddy_viscosity(DoubleTab& nu_t) const
{
  nu_t = 0; // no eddy viscosity here ; use ggdh for diffusion
}

void Viscosite_turbulente_WIT::reynolds_stress(DoubleTab& R_ij) const // Renvoie <u_i'u_j'>
{
  const DoubleTab& tab_k_WIT = pb_->get_champ("k_WIT").passe();

  int i, d, db, D = dimension, N = R_ij.dimension(1);
  if (N!=2) Process::exit("WIT only works when there are 2 phases");

  for (i = 0; i < R_ij.dimension(0); i++)
    for (d = 0; d < D; d++)
      for (db = 0; db < D; db++) //on ne remplit que les phases concernees par k
        {
          R_ij(i, 1, d, db) = 0 ; // No WIT for gas phase
          R_ij(i, 0, d, db) = 2./3. * (db==d?1:0) * tab_k_WIT(i, 0) ;
        }
}

void Viscosite_turbulente_WIT::k_over_eps(DoubleTab& k_sur_eps) const
{
  k_sur_eps =  0;
}

void Viscosite_turbulente_WIT::eps(DoubleTab& eps) const
{
  eps =  0;
}
