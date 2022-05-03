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
// File:        Flux_parietal_Kommajosyula.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Flux_parietal_Kommajosyula.h>
#include <Flux_parietal_adaptatif.h>
#include <Loi_paroi_adaptative.h>
#include <Correlation.h>
#include <Pb_Multiphase.h>
#include <Domaine_dis.h>
#include <Zone_VF.h>
#include <Op_Diff_PolyMAC_base.h>
#include <Milieu_composite.h>
#include <Saturation_base.h>

#include <math.h>

Implemente_instanciable(Flux_parietal_Kommajosyula, "Flux_parietal_Kommajosyula", Flux_parietal_base);

Sortie& Flux_parietal_Kommajosyula::printOn(Sortie& os) const
{
  return os;
}

Entree& Flux_parietal_Kommajosyula::readOn(Entree& is)
{
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, pb_.valeur());
  correlation_monophasique_.typer_lire(pbm, "Flux_parietal", is);
  Cout << que_suis_je() << " : single-phase wall heat flux is " << correlation_monophasique_->que_suis_je() << finl;

  Param param(que_suis_je());
  param.ajouter("contact_angle_deg",&theta_,Param::REQUIRED);
  param.ajouter("molar_mass",&molar_mass_,Param::REQUIRED);
  param.lire_avec_accolades(is);

  if ( !sub_type(Milieu_composite, pb_->milieu())) Process::exit("Flux_parietal_Kommajosyula::readOn : the medium must be composite !");
  if (!pbm.nom_phase(0).debute_par("liquide")) Process::exit("Flux_parietal_Kommajosyula::readOn : the first phase must be liquid !");

  const Milieu_composite& milc = ref_cast(Milieu_composite, pb_->milieu());
  int n_g = -1;
  for (int k = 1 ; k<pbm.nb_phases() ; k++) if (milc.has_saturation(0, k)) n_g += 1;
  if (n_g > 0) Process::exit("Flux_parietal_Kommajosyula::readOn : there can only be one evaporating phase for the carrying liquid for now ! Please feel free to update the code if you need.");

  return is;
}

void Flux_parietal_Kommajosyula::completer()
{
  correlation_monophasique_->completer();
}

void Flux_parietal_Kommajosyula::qp(int N, int f, double D_h, double D_ch,
                                    const double *alpha, const double *T, const double p, const double *v, const double Tp,
                                    const double *lambda, const double *mu, const double *rho, const double *Cp,
                                    DoubleTab *qpk, DoubleTab *da_qpk, DoubleTab *dp_qpk, DoubleTab *dv_qpk, DoubleTab *dTf_qpk, DoubleTab *dTp_qpk,
                                    DoubleTab *qpi, DoubleTab *da_qpi, DoubleTab *dp_qpi, DoubleTab *dv_qpi, DoubleTab *dTf_qpi, DoubleTab *dTp_qpi,
                                    DoubleTab *d_nuc, int& nonlinear) const
{
  // On remplit le monophasique ; pas besoin du flux interfacial normalement
  ref_cast(Flux_parietal_base, correlation_monophasique_.valeur()).qp(N, f, D_h, D_ch, alpha, T, p, v, Tp, lambda, mu, rho, Cp,
                                                                      qpk, da_qpk, dp_qpk, dv_qpk, dTf_qpk, dTp_qpk, NULL, NULL, NULL, NULL, NULL, NULL,
                                                                      d_nuc, nonlinear);

  // Ici la phase liquide est forcement la phase 0 car la correlation monophasique ne remplit que la phase 0
  int n_l = 0 ;
  const Milieu_composite& milc = ref_cast(Milieu_composite, pb_->milieu());

  for (int k = 0 ; k<N ; k++) if (n_l != k) if (milc.has_saturation(n_l, k))
        {
          Saturation_base& sat = milc.get_saturation(n_l, k);

          // Nucleation site density (Hibiki Ishii 2003)
          double N_sites = Hibiki_Ishii_Site_density(rho[k], rho[n_l], T[k], T[n_l], p, Tp, sat.Lvap(p), sat.Tsat(p), sat.sigma(sat.Tsat(p), p), theta_, molar_mass_);

          double Delta_T_sup = Tp - sat.Tsat(p); // Wall superheat
          double Delta_T_sub = sat.Tsat(p) - T[n_l] ; // Subcooling
          double Ja_sup = rho[n_l]*Cp[n_l]*Delta_T_sup/(rho[k] * sat.Lvap(p));// Superheat Jakob number
          double Ja_sub = rho[n_l]*Cp[n_l]*Delta_T_sub/(rho[k] * sat.Lvap(p));// Subcooling Jakob number

          double u_bulk = 0;
          if (sub_type(Flux_parietal_adaptatif, correlation_monophasique_.valeur()))
            {
              const Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, ref_cast(Pb_Multiphase, pb_.valeur()).get_correlation("Loi_paroi").valeur());
              const double u_tau = corr_loi_paroi.get_utau(f);
              u_bulk = 15.*u_tau; // Big approximation...
            }
          else Cerr << "Flux_parietal_Kommajosyula::qp isn't adapted to a single-phase " << correlation_monophasique_->que_suis_je() << " heat flux for now ! " , Process::exit();

          // Single phase heat transfer coefficient
          double h_fc = (*dTp_qpk)(n_l);

          // Departure diameter (page 47 Kommajosyula PhD)
          double D_d = 18.9e-6 * std::pow( (rho[n_l]-rho[k])/rho[k], 0.27 ) * std::pow(Ja_sup, 0.75) * std::pow(1+Ja_sub, -0.3) * std::pow(u_bulk, -0.26);

          // Liftoff diameter (page 47 Kommajosyula PhD)
          double D_lo = 1.2*D_d ;

          // Bubble growth constant (page 42 Kommajosyula PhD) excluding the microlayer contribution ; ! different formulation in the manuscript and algorithmic appendix
          double chi = 1.55 - 0.05*Delta_T_sub/Delta_T_sup;
          double K = chi * 2*std::sqrt(3/M_PI) *Ja_sup*std::sqrt(lambda[n_l]/(rho[n_l]*Cp[n_l]));

          // Bubble growth time (page 44 Kommajosyula PhD)
          double t_g = D_d*D_d/(16*K*K) ;

          // Bubble wait time (page 56 Kommajosyula PhD)
          double t_w = 0.0061*std::pow(Ja_sub, 0.63) /Delta_T_sup;

          // Bubble departure frequency (page 51 Kommajosyula PhD)
          double f_dep = 1/(t_g+t_w) ;

          // Time necessary to reform the termal boundary layer (page 61 Kommajosyula PhD ; error in equation 4.5 of the manuscript)
          double t_star = lambda[n_l] * rho[n_l]*Cp[n_l]/(h_fc*h_fc*M_PI) ;

          // Active nucleation site density (page 66 Kommajosyula PhD)
          double N_0 = f_dep * t_g * M_PI * D_d*D_d/4;
          double N_active = Lambert_W_function(N_0 * N_sites)/N_0;

          // Single bubble sliding area
          double A_sl = 1/std::sqrt(N_active)*(D_d+D_lo)/2;
          // Percentage of the surface occupied by sliding bubbles
          double S_sl = std::min(1., A_sl * N_active * f_dep * t_star) ;

          // We correct the single phase heat flux
          if (qpk)     (*qpk)     *= (1-S_sl);
          if (da_qpk)  (*da_qpk)  *= (1-S_sl);
          if (dp_qpk)  (*dp_qpk)  *= (1-S_sl);
          if (dv_qpk)  (*dv_qpk)  *= (1-S_sl);
          if (dTf_qpk) (*dTf_qpk) *= (1-S_sl);
          if (dTp_qpk) (*dTp_qpk) *= (1-S_sl);

          // Bubble sliding
          if (qpk)     (*qpk)     += S_sl*2*h_fc*(Tp-T[n_l]);
          if (dTf_qpk) (*dTf_qpk) -= S_sl*2*h_fc;
          if (dTp_qpk) (*dTp_qpk) += S_sl*2*h_fc;

          // Evaporation (calculer les derivees/T apres)
          if (qpi)  (*qpi)(n_l, k) = 1./6.*M_PI*D_lo*D_lo*D_lo*rho[k]*sat.Lvap(p)*f_dep*N_active;
          if (qpi)  (*qpi)(k, n_l) = - (*qpi)(n_l, k);

          if (d_nuc) (*d_nuc)(k) = D_lo;
        }
}

double Flux_parietal_Kommajosyula::Lambert_W_function(double x)  const
{
  if (x < 1/M_E) return x;
  else if (x > M_E) return std::log(x) - std::log(log(x));
  else return 0.2689*x+0.2690;
}

// See Hibiki, Ishii 2003
double Flux_parietal_Kommajosyula::Hibiki_Ishii_Site_density(double rho_v, double rho_l, double T_v, double T_l, double p, double Tp, double L_vap, double T_sat, double sigma, double theta, double molar_mass) const
{
  double N_bar = 4.72e5;
  double mu2 = 0.722*180/M_PI;
  mu2 *= mu2;

  double rho_p = std::log( (rho_l-rho_v)/rho_v);
  double f_rho_plus = -0.01064 + 0.48246*rho_p - 0.22712*rho_p*rho_p + 0.05468*rho_p*rho_p*rho_p;

  double lambda_prime = 2.5e-6;
  double R = 8.314462618 / molar_mass ;

  double inv_Rc = p * (-1 + std::exp(L_vap*(T_v-T_sat)/( R * T_v * T_sat ))) / (2 * sigma * (1+rho_v/rho_l));

  return N_bar * (1 - std::exp(-theta*theta/(8*mu2))) * (-1 + std::exp(f_rho_plus * lambda_prime * inv_Rc));
}
