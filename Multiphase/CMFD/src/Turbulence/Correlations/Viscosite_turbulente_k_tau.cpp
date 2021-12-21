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
// File:        Viscosite_turbulente_k_tau.cpp
// Directory:   $TRUST_ROOT/src/Turbulence/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Viscosite_turbulente_k_tau.h>
#include <Param.h>
#include <Probleme_base.h>
#include <Champ_base.h>

Implemente_instanciable(Viscosite_turbulente_k_tau, "Viscosite_turbulente_k_tau", Viscosite_turbulente_base);

Sortie& Viscosite_turbulente_k_tau::printOn(Sortie& os) const
{
  return os;
}

Entree& Viscosite_turbulente_k_tau::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("limiter|limiteur", &limiter_);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Viscosite_turbulente_k_tau::eddy_viscosity(DoubleTab& nu_t) const
{
  //tout en explicite
  const DoubleTab& k = pb_->get_champ("k").passe(), &tau = pb_->get_champ("tau").passe(),
                   &nu = pb_->get_champ("viscosite_cinematique").passe();
  //il faut que nu_t et k aient la meme localisation et que nu_t ait au moins autant de composantes que k
  assert(nu_t.dimension_tot(0) == k.dimension_tot(0) && k.dimension(1) <= nu_t.dimension(1));
  //on met 0 pour les composantes au-dela de k.dimension(1) (ex. : vapeur dans Pb_Multiphase)
  for (int i = 0; i < nu_t.dimension_tot(0); i++) for (int n = 0; n < nu_t.dimension(1); n++)
      nu_t(i, n) = n < k.dimension(1) ? max(k(i, n) * tau(i, n), limiter_ * nu(i, n)) : 0;
}

void Viscosite_turbulente_k_tau::reynolds_stress(DoubleTab& R_ij) const
{
  abort(); //TBD
}