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
// File:        Dispersion_bulles_turbulente_Burns.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Dispersion_bulles_turbulente_Burns.h>
#include <Pb_Multiphase.h>
#include <QDM_Multiphase.h>
#include <TRUSTTrav.h>
#include <Frottement_interfacial_base.h>
#include <math.h>
//#include <Energie_cinetique_turbulente_WIT.h>

Implemente_instanciable(Dispersion_bulles_turbulente_Burns, "Dispersion_bulles_turbulente_Burns", Dispersion_bulles_base);

Sortie& Dispersion_bulles_turbulente_Burns::printOn(Sortie& os) const
{
  return os;
}

Entree& Dispersion_bulles_turbulente_Burns::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("minimum", &minimum_);
  param.ajouter("a_res", &a_res_);
  param.ajouter("g_", &g_);
  param.ajouter("coefBIA_", &coefBIA_);
  param.lire_avec_accolades_depuis(is);

  const Pb_Multiphase *pbm = sub_type(Pb_Multiphase, pb_.valeur()) ? &ref_cast(Pb_Multiphase, pb_.valeur()) : nullptr;

  if (!pbm || pbm->nb_phases() == 1) Process::exit(que_suis_je() + " : not needed for single-phase flow!");
  for (int n = 0; n < pbm->nb_phases(); n++) //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    if (pbm->nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))  n_l = n;
  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");

  if (pbm->has_correlation("frottement_interfacial")) correlation_drag_ = pbm->get_correlation("frottement_interfacial"); //correlation fournie par le bloc correlation
  else Correlation_base::typer_lire_correlation(correlation_drag_, *pbm, "frottement_interfacial", is); //sinon -> on la lit

  return is;
}

void Dispersion_bulles_turbulente_Burns::completer()
{
  if ((a_res_ == -1) && (ref_cast(QDM_Multiphase, pb_->equation(0)).alpha_res()<0.99) ) // Pas en homogene
    a_res_ = std::min(1., std::max(1.e-4, ref_cast(QDM_Multiphase, pb_->equation(0)).alpha_res()*100.));
  else if (a_res_ == -1) a_res_ = 1.e-6;
}

void Dispersion_bulles_turbulente_Burns::coefficient(const input_t& in, output_t& out) const
{
  const Frottement_interfacial_base& corr = ref_cast(Frottement_interfacial_base, correlation_drag_->valeur());
  int N = out.Ctd.dimension(0);

  DoubleTrav coeff_drag(N, N, 2);
  corr.coefficient( in.alpha, in.p, in.T, in.rho, in.mu, in.sigma, in.dh, in.nv, in.d_bulles, coeff_drag);

  out.Ctd = 0;

  for (int k = 0; k < N; k++)
    if (k!=n_l)
      {
        double nuBIA = 0.;

        if (coefBIA_)
          {
            // Calcul de nuBIA = (k_WIT+k_WIF)/omega_WIT
            double u_r = in.nv(k,n_l); // vitesse relative
            double Reb = in.rho[n_l]*in.d_bulles[k]*u_r/in.mu[n_l]; // Reynolds bulle
            int ind_trav = (k>n_l) ? (n_l*(N-1)-(n_l-1)*(n_l)/2) + (k-n_l-1) : (k*(N-1)-(k-1)*(k)/2) + (n_l-k-1);
            double Eo = g_ * std::abs(in.rho[n_l] - in.rho[k]) * in.d_bulles[k] * in.d_bulles[k]/in.sigma[ind_trav]; // Eotvos
            double Cd = (u_r!=0) ? std::max( std::min( 16./Reb*(1.+0.15*std::pow(Reb, 0.687)) , 48./Reb )   , 8.*Eo/(3.*(Eo+4.))) : 0.; // si u_r=0 alors pas de trainée, pas de WIT donc dissipation=0
            double omega_WIT = 2.0 * in.mu[n_l] * Cd * Reb / (C_lambda_*C_lambda_*in.d_bulles[k]*in.d_bulles[k]); // dissipation spécifique de la WIT (definie comme epsilon_WIT/kWIT)
            //double k_WIF = 1./2. * in.alpha[k] * (u_r*u_r * 0.5 + 3./2.*0.25*gamma_*gamma_*gamma_); // energie cinetique turbulente de la composante WIF
            nuBIA = coefBIA_ * (omega_WIT == 0.0 ? 0.0 : in.k_WIT/omega_WIT); // l'équivalent de nu_t pour l'agitation induite par les bulles (BIA)
          }

        // Calcul des coefficients de dispersion turbulente
        out.Ctd(k, n_l) = std::max( minimum_, (in.alpha[k]  >a_res_) ? (nuBIA + in.nut[n_l])/Prt_ * coeff_drag(k, n_l, 0)/in.alpha[k]  : (nuBIA + in.nut[n_l])/Prt_ * coeff_drag(k, n_l, 0)*in.alpha[k]  /(a_res_*a_res_) );
        out.Ctd(n_l, k) = std::max( minimum_, (in.alpha[n_l]>a_res_) ? (nuBIA + in.nut[n_l])/Prt_ * coeff_drag(n_l, k, 0)/in.alpha[n_l]: (nuBIA + in.nut[n_l])/Prt_ * coeff_drag(n_l, k, 0)*in.alpha[n_l]/(a_res_*a_res_) );
      }
}
