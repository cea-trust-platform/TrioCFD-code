/****************************************************************************
* Copyright (c) 2023, CEA
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

#include <Portance_interfaciale_Sugrue_Modifiee.h>
#include <Frottement_interfacial_base.h>
#include <Pb_Multiphase.h>
#include <math.h>

Implemente_instanciable(Portance_interfaciale_Sugrue_Modifiee, "Portance_interfaciale_Sugrue_Modifiee", Portance_interfaciale_base);

Sortie& Portance_interfaciale_Sugrue_Modifiee::printOn(Sortie& os) const
{
  return os;
}

Entree& Portance_interfaciale_Sugrue_Modifiee::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("constante_gravitation", &g_);
  param.lire_avec_accolades_depuis(is);

  const Pb_Multiphase *pbm = sub_type(Pb_Multiphase, pb_.valeur()) ? &ref_cast(Pb_Multiphase, pb_.valeur()) : nullptr;

  if (!pbm || pbm->nb_phases() == 1) Process::exit(que_suis_je() + " : not needed for single-phase flow!");
  for (int n = 0; n < pbm->nb_phases(); n++) //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    if (pbm->nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))  n_l = n;

  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");

  return is;
}

void Portance_interfaciale_Sugrue_Modifiee::coefficient(const input_t& in, output_t& out) const
{
  int k, N = out.Cl.dimension(0);
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, pb_.valeur());

  for (k = 0; k < N; k++)
    if (k!=n_l) // k gas phase
      {
        int ind_trav = (k>n_l) ? (n_l*(N-1)-(n_l-1)*(n_l)/2) + (k-n_l-1) : (k*(N-1)-(k-1)*(k)/2) + (n_l-k-1);

        // Newton method to determine dv along gravity
        double dv0 = 0.2, epsilon = 1.e-4; // Initialize dv at random
        int step = 1, iter_max = 20;
        DoubleTab dv(N, N), coeff(N, N, 2), alpha_l(N);
        const Frottement_interfacial_base& correlation_fi = ref_cast(Frottement_interfacial_base, pbm.get_correlation("frottement_interfacial").valeur());
        double sum_alpha = 0;
        for (int n=0; n<N ; n++) alpha_l(n)= std::max(in.alpha(n), alpha_lim_), sum_alpha+=alpha_l(n);
        for (int n=0; n<N ; n++) alpha_l(n)/=sum_alpha;

        do
          {
            dv(k,n_l) = dv0;
            dv(n_l,k) = dv0;
            correlation_fi.coefficient(alpha_l, in.p, in.T, in.rho, in.mu, in.sigma, in.dh, dv, in.d_bulles, coeff);
            dv0 = dv0 - (coeff(n_l, k, 0)*dv0 - g_*alpha_l(k)*(in.rho[n_l]*in.alpha[n_l]+in.rho[k]*in.alpha[k] - in.rho[k])) / (coeff(n_l, k, 1)*dv0 + coeff(n_l, k, 0));
            step = step+1;
            if(step > iter_max) Process::exit(que_suis_je() + " : Newton algorithm not converging to find relative velocity !");
          }
        while(std::abs(coeff(n_l, k, 0)*dv0 - g_*alpha_l(k)*(in.rho[n_l]*in.alpha[n_l]+in.rho[k]*in.alpha[k]- in.rho[k])) > epsilon);

        // Rest of the correlation

        const double Eo = g_ * std::abs(in.rho[n_l]-in.rho[k]) * in.d_bulles[k]*in.d_bulles[k]/in.sigma[ind_trav];

        double Wo = std::min(Eo * in.k_turb[n_l]/std::max(dv0*dv0, 1.e-8) , 6.); // Experimental validation up to Wo=6
        double f_Wo = std::min(0.03, 5.0404 - 5.0781*std::pow(Wo, 0.0108));
        double f_alp= in.alpha[k] < .7 ? 1 : ( in.alpha[k] > .95 ? 0 : (.95-in.alpha[k])/.25 );
        double Cl = f_Wo*f_alp ;

        out.Cl(k, n_l) = Cl * in.rho[n_l] * in.alpha[k] ;
        out.Cl(n_l, k) = out.Cl(k, n_l);
      }
}
