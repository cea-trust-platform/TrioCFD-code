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
#include <TRUSTTab_parts.h>

Implemente_instanciable(Viscosite_turbulente_k_tau, "Viscosite_turbulente_k_tau", Viscosite_turbulente_base);
// XD type_diffusion_turbulente_multiphase_k_tau type_diffusion_turbulente_multiphase_deriv k_tau 1 not_set
Sortie& Viscosite_turbulente_k_tau::printOn(Sortie& os) const
{
  return os;
}

Entree& Viscosite_turbulente_k_tau::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("limiter|limiteur", &limiter_); // XD_ADD_P chaine not_set
  param.ajouter("sigma", &sigma_); // XD_ADD_P floattant not_set
  param.ajouter("beta_k", &beta_k_); // XD_ADD_P floattant not_set
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Viscosite_turbulente_k_tau::eddy_viscosity(DoubleTab& nu_t) const
{
  //tout en explicite
  const DoubleTab& k = pb_->get_champ("k").passe(), &tau = pb_->get_champ("tau").passe(),
                   &nu = pb_->get_champ("viscosite_cinematique").passe();
  const int cnu = nu.dimension(0) == 1;
  //il faut que nu_t et k aient la meme localisation et que nu_t ait au moins autant de composantes que k
  assert((k.dimension(1) <= nu_t.dimension(1)));
  //on met 0 pour les composantes au-dela de k.dimension(1) (ex. : vapeur dans Pb_Multiphase)
  for (int i = 0; i < nu_t.dimension(0); i++)
    for (int n = 0; n < nu_t.dimension(1); n++)
      nu_t(i, n) = (n < k.dimension(1)) ? sigma_ * std::max(k(i, n) * tau(i, n), limiter_ * nu(!cnu * i, n)) : 0;
}

void Viscosite_turbulente_k_tau::reynolds_stress(DoubleTab& R_ij) const // Renvoie <u_i'u_j'>
{
  const DoubleTab& k = pb_->get_champ("k").passe(), &tau = pb_->get_champ("tau").passe(),
                   &nu = pb_->get_champ("viscosite_cinematique").passe(), &grad_u = pb_->get_champ("gradient_vitesse").passe();
  ConstDoubleTab_parts p_gu(grad_u); //en PolyMAC_P0, grad_u contient (nf.grad)u_i aux faces, puis (d_j u_i) aux elements
  int i, d, db, D = dimension, i_part = -1, n, N = nu.dimension(1), Nk = k.dimension(1);
  for (i = 0; i < p_gu.size(); i++)
    if (p_gu[i].get_md_vector() == R_ij.get_md_vector()) i_part = i; //on cherche une partie ayant le meme support que k
  if (i_part < 0) Process::exit("Viscosite_turbulente_k_tau : inconsistency between velocity gradient and k!");
  const DoubleTab& gu = p_gu[i_part]; //le bon tableau
  for (i = 0; i < R_ij.dimension(0); i++)
    for (n = 0; n < N; n++)
      {
        double sum_diag = 0.;
        double nut_loc = n < Nk ? std::max(k(i, n) * tau(i, n), limiter_ * nu(i, n)) : 0 ;
        for (d = 0; d < D; d++) sum_diag += gu(i, d, D * n + d) ;
        for (d = 0; d < D; d++)
          for (db = 0; db < D; db++) //on ne remplit que les phases concernees par k
            R_ij(i, n, d, db) = n < Nk ? sigma_ * (2. / 3. * (k(i, n) + nut_loc * sum_diag) * (d ==db) - nut_loc * (gu(i, d, D * n + db) + gu(i, db, D * n + d))) : 0;
      }
}

void Viscosite_turbulente_k_tau::k_over_eps(DoubleTab& k_sur_eps) const
{
  const DoubleTab& tau = pb_->get_champ("tau").passe();
  int i, nl = k_sur_eps.dimension(0), n, N = k_sur_eps.dimension(1), Nt = tau.dimension(1);
  assert(nl == tau.dimension(0) && Nt <= N);
  /* comme tau = 1 / omega et omega = epsilon / (k*beta_k), k / epsilon = tau/beta_k ! */
  for (i = 0; i < nl; i++)
    for (n = 0; n < N; n++) k_sur_eps(i, n) = n < Nt ? tau(i, n)/beta_k_  : 0;
}

void Viscosite_turbulente_k_tau::eps(DoubleTab& eps_) const
{
  const DoubleTab& tau = pb_->get_champ("tau").passe(),
                   &k = pb_->get_champ("k").passe(),
                    &nu = pb_->get_champ("viscosite_cinematique").passe();
  int i, nl = eps_.dimension(0), n, N = eps_.dimension(1), Nt = tau.dimension(1);
  assert(nl == tau.dimension(0) && Nt <= N);
  /* comme tau = 1 / omega et omega = epsilon / (k*beta_k), epsilon = beta_k_ * k / tau ! */
  for (i = 0; i < nl; i++)
    for (n = 0; n < N; n++) eps_(i, n) = beta_k_ * ( ((n < Nt) && (k(i, n)>1.e-8) ) ? k(i, n)*k(i, n)/ std::max(k(i, n) * tau(i, n), limiter_ * nu(i, n)) : 0 );
}
