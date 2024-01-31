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
// File:        Viscosite_turbulente_k_omega.cpp
// Directory:   $TRUST_ROOT/src/Turbulence/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Viscosite_turbulente_k_omega.h>
#include <Pb_Multiphase.h>
#include <Masse_ajoutee_base.h>
#include <TRUSTTab_parts.h>
#include <Probleme_base.h>
#include <Champ_base.h>
#include <Param.h>

Implemente_instanciable(Viscosite_turbulente_k_omega, "Viscosite_turbulente_k_omega", Viscosite_turbulente_base);
// XD type_diffusion_turbulente_multiphase_k_omega type_diffusion_turbulente_multiphase_deriv k_omega 1 not_set

Sortie& Viscosite_turbulente_k_omega::printOn(Sortie& os) const
{
  return os;
}

Entree& Viscosite_turbulente_k_omega::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("limiter|limiteur", &limiter_); // XD_ADD_P chaine not_set
  param.ajouter("sigma", &sigma_); // XD_ADD_P floattant not_set
  param.ajouter("beta_k", &beta_k_); // XD_ADD_P floattant not_set
  param.ajouter("gas_turb", &gas_turb_); // XD_ADD_P flag not_set
  param.lire_avec_accolades_depuis(is);

  return is;
}

void Viscosite_turbulente_k_omega::completer()
{
  const Pb_Multiphase* pbm = sub_type(Pb_Multiphase, pb_.valeur()) ? &ref_cast(Pb_Multiphase, pb_.valeur()) : nullptr ;
  if ( (gas_turb_) && !(pbm) ) Process::exit(que_suis_je() + " : there must be multiphase problem if you want gas phase turbulence !");
  if ( (gas_turb_) && (!pbm->has_correlation("masse_ajoutee")) ) Process::exit(que_suis_je() + " : there must be an added mass correlation if you want gas phase turbulence !");
  if ( (gas_turb_) && (!int(pbm->has_champ("alpha"))) ) Process::exit(que_suis_je() + " : there must be void fraction if you want gas phase turbulence !");
  if (pbm->has_correlation("masse_ajoutee")) correlation_ = pbm->get_correlation("masse_ajoutee");
}

void Viscosite_turbulente_k_omega::eddy_viscosity(DoubleTab& nu_t) const
{
  //tout en explicite
  const DoubleTab& k = pb_->get_champ("k").passe(), &omega = pb_->get_champ("omega").passe(),
                   &nu = pb_->get_champ("viscosite_cinematique").passe(),
                    &rho = pb_->get_champ("masse_volumique").passe(),
                     *alpha = pb_->has_champ("alpha") ? &(pb_->get_champ("alpha").passe()) : nullptr ;
  //il faut que nu_t et k aient la meme localisation et que nu_t ait au moins autant de composantes que k
  assert(k.dimension(1) <= nu_t.dimension(1));
  //on met 0 pour les composantes au-dela de k.dimension(1) (ex. : vapeur dans Pb_Multiphase)
  for (int i = 0; i < nu_t.dimension(0); i++)
    for (int n = 0; n < nu_t.dimension(1); n++)
      //      nu_t(i, n) = sigma_ * ( (omega(i,n*(k.dimension(1)-1)) > 0.) ?
      //                   std::max(k(i, n*(k.dimension(1)-1)) / omega(i, n*(k.dimension(1)-1)), limiter_ * nu(i, n)):
      //                   limiter_ * nu(i, n) )  ;
      nu_t(i, n) = (n<k.dimension(1)) ? sigma_ * ( (omega(i,n) > 0.) ? std::max(k(i, n) / omega(i, n), limiter_ * nu(i, n)): limiter_ * nu(i, n) ) : 0. ;

  if (gas_turb_)
    {
      DoubleTab coeff = DoubleTab(nu_t.dimension(1));
      for (int i = 0; i < nu_t.dimension(0); i++)
        {
          const Masse_ajoutee_base& corr_ma_ = ref_cast(Masse_ajoutee_base, correlation_.valeur());
          corr_ma_.coefficient( & (*alpha)(i,0), &rho(i,0), coeff);
          for (int n = k.dimension(1) ; n < nu_t.dimension(1) ; n++)
            nu_t(i, n) =  nu_t(i, 0.)  * (1. + coeff(n)* rho(i,0)/rho(i,n)) * std::min((*alpha)(i,n)*10, 1.) ;
        }
    }

}

void Viscosite_turbulente_k_omega::reynolds_stress(DoubleTab& R_ij) const // Renvoie <u_i'u_j'>
{
  const DoubleTab& k = pb_->get_champ("k").passe(), &omega = pb_->get_champ("omega").passe(),
                   &nu = pb_->get_champ("viscosite_cinematique").passe(), &grad_u = pb_->get_champ("gradient_vitesse").passe();
  ConstDoubleTab_parts p_gu(grad_u); //en PolyMAC_P0, grad_u contient (nf.grad)u_i aux faces, puis (d_j u_i) aux elements
  int i, d, db, D = dimension, i_part = -1, n, N = nu.dimension(1), Nk = k.dimension(1);
  for (i = 0; i < p_gu.size(); i++)
    if (p_gu[i].get_md_vector() == R_ij.get_md_vector()) i_part = i; //on cherche une partie ayant le meme support que k
  if (i_part < 0) Process::exit("Viscosite_turbulente_k_omega : inconsistency between velocity gradient and k!");
  const DoubleTab& gu = p_gu[i_part]; //le bon tableau
  for (i = 0; i < R_ij.dimension(0); i++)
    for (n = 0; n < N; n++)
      {
        double sum_diag = 0.;
        double nut_loc = n < Nk ? (omega(i,n) > 0.) ? std::max(k(i, n) / omega(i, n), limiter_ * nu(i, n)): limiter_ * nu(i, n) : 0 ;
        for (d = 0; d < D; d++) sum_diag += gu(i, d, D * n + d) ;
        for (d = 0; d < D; d++)
          for (db = 0; db < D; db++) //on ne remplit que les phases concernees par k
            R_ij(i, n, d, db) = (n < Nk)    ?    sigma_* (  2. / 3. * (k(i, n) + nut_loc * sum_diag) * (d ==db) - nut_loc * (gu(i, d, D * n + db) + gu(i, db, D * n + d))    ) : 0;
      }
}

void Viscosite_turbulente_k_omega::k_over_eps(DoubleTab& k_sur_eps) const
{
  const DoubleTab& omega = pb_->get_champ("omega").passe();
  int i, nl = k_sur_eps.dimension(0), n, N = k_sur_eps.dimension(1), Nt = omega.dimension(1);
  assert(nl == omega.dimension(0) && Nt <= N);
  for (i = 0; i < nl; i++)
    for (n = 0; n < N; n++) k_sur_eps(i, n) = (n < Nt) ? ( (omega(i,n) > 0.) ? 1/(omega(i, n)*beta_k_) : 0) : 0;
}

void Viscosite_turbulente_k_omega::eps(DoubleTab& eps_) const
{
  const DoubleTab& omega = pb_->get_champ("omega").passe(),
                   &k = pb_->get_champ("k").passe();
  int i, nl = eps_.dimension(0), n, N = eps_.dimension(1), Nt = omega.dimension(1);
  assert(nl == omega.dimension(0) && Nt <= N);
  for (i = 0; i < nl; i++)
    for (n = 0; n < N; n++) eps_(i, n) = beta_k_ * ((n < Nt) ? k(i,n) * omega(i, n) : 0);
}
