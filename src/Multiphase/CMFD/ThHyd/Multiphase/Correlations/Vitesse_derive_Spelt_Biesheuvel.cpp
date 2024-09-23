/****************************************************************************
* Copyright (c) 2022, CEA
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
// File:        Vitesse_derive_Spelt_Biesheuvel.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Frottement_interfacial_base.h>
#include <Dispersion_bulles_base.h>
#include <Vitesse_derive_Spelt_Biesheuvel.h>
#include <Pb_Multiphase.h>
#include <cmath>

Implemente_instanciable(Vitesse_derive_Spelt_Biesheuvel, "Vitesse_relative_derive_Spelt_Biesheuvel", Vitesse_derive_base);

Sortie& Vitesse_derive_Spelt_Biesheuvel::printOn(Sortie& os) const { return os; }

Entree& Vitesse_derive_Spelt_Biesheuvel::readOn(Entree& is)
{
  Vitesse_derive_base::readOn(is);
  return is;
}

void Vitesse_derive_Spelt_Biesheuvel::completer()
{
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, pb_.valeur());
  if (!pbm.has_correlation("frottement_interfacial")) Process::exit(que_suis_je() + " : there must be an interfacial friction correlation in the problem !");
  if (pbm.has_correlation("dispersion_bulles")) needs_grad_alpha_ = 1;
}

void Vitesse_derive_Spelt_Biesheuvel::evaluate_C0_vg0(const input_t& in) const
{
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, pb_.valeur());
  const int D = dimension, N = in.alpha.dimension(0);
  const double norm_g = sqrt(local_carre_norme_vect(in.g));

  // Newton method to determine dv along gravity
  double dv0 = 0.2, epsilon = 1.e-4; // Initialize dv at random
  int step = 1, iter_max = 20;
  DoubleTab p, T, dv(N, N), coeff(N, N, 2), alpha_l(N);
  const Frottement_interfacial_base& correlation_fi = ref_cast(Frottement_interfacial_base, pbm.get_correlation("frottement_interfacial"));
  for (int n=0; n<N ; n++) alpha_l(n)= std::max(in.alpha(n), 1.e-3);

  do
    {
      dv(n_g,n_l) = dv0;
      dv(n_l,n_g) = dv0;
      correlation_fi.coefficient(alpha_l, p, T, in.rho, in.mu, in.sigma, in.dh, dv, in.d_bulles, coeff);
      dv0 = dv0 - (coeff(n_l, n_g, 0)*dv0 - norm_g*alpha_l(n_g)*(in.rho[n_l]*in.alpha[n_l]+in.rho[n_g]*in.alpha[n_g] - in.rho[n_g])) / (coeff(n_l, n_g, 1)*dv0 + coeff(n_l, n_g, 0));
      step = step+1;
      if(step > iter_max) Process::exit(que_suis_je() + " : Newton algorithm not converging to find relative velocity !");
    }
  while(std::abs(coeff(n_l, n_g, 0)*dv0 - norm_g*alpha_l(n_g)*(in.rho[n_l]*in.alpha[n_l]+in.rho[n_g]*in.alpha[n_g]- in.rho[n_g])) > epsilon);

  /* distribution parameter */
  C0 = 1;

  /* drift velocity along gravity */
  for (int d = 0; d < D; d++) vg0(d) = - dv0 * in.g(d) / norm_g;

  /* Turbulent diffusion ; possible stability issues in the future... */
  double u_turb = std::pow(0.09, 0.25) * std::sqrt(in.k[n_l]) ;
  double beta = u_turb / dv0 ;
  double epsilon_turb = 0.09 * in.k[n_l]*in.k[n_l] /in.nut[n_l] ;
  double L = std::pow(0.09, 0.75)  * std::pow(in.k[n_l], 1.5) / epsilon_turb ;
  double tau_rel = dv0/(2*norm_g);
  double lambda_turb = std::sqrt(10*in.mu[n_l]/in.rho[n_l]*in.k[n_l]/epsilon_turb);
  double mu_turb = lambda_turb/tau_rel ;
  double diffu_Spelt_Biesheuvel = 0.5*beta*u_turb*L + 3./8.*beta*u_turb*L*(tau_rel/mu_turb)*(tau_rel/mu_turb)*lambda_turb/L;

  for (int d = 0; d < D; d++) vg0(d) +=  - diffu_Spelt_Biesheuvel * in.gradAlpha(n_g, d)/std::max(in.alpha[n_g], 1.e-4) ;

  for (int d = 0; d < D; d++) vg0(d) *=  (1.0 - C0 * in.alpha[n_g]) ;
}
