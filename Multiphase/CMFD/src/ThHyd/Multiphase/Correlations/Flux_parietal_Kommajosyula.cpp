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
#include <TRUSTTrav.h>
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
  for (int k = 1 ; k<pbm.nb_phases() ; k++)
    if (milc.has_saturation(0, k)) n_g += 1;
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

  for (int k = 0 ; k<N ; k++)
    if (n_l != k)
      if (milc.has_saturation(n_l, k))
        {
          Saturation_base& sat = milc.get_saturation(n_l, k);

          double Delta_T_sup = Tp - sat.Tsat(p); // Wall superheat

          if (Delta_T_sup > 0) // Else : no wall superheat => no nucleation => single phase heat transfer only
            {
              double Delta_T_sub = std::max(sat.Tsat(p) - T[n_l], 1.e-8) ; // Subcooling ; non negative for numerical reasons
              double dTp_Delta_T_sup = 1.;
              double dTl_Delta_T_sub = -1.;
              double Ja_sup = rho[n_l]*Cp[n_l]*Delta_T_sup/(rho[k] * sat.Lvap(p));// Superheat Jakob number
              double Ja_sub = rho[n_l]*Cp[n_l]*Delta_T_sub/(rho[k] * sat.Lvap(p));// Subcooling Jakob number
              double dTp_Ja_sup = rho[n_l]*Cp[n_l]/(rho[k] * sat.Lvap(p));
              double dTl_Ja_sub =-rho[n_l]*Cp[n_l]/(rho[k] * sat.Lvap(p));

              // Nucleation site density (Hibiki Ishii 2003)
              double N_sites =              Hibiki_Ishii_Site_density(rho[k], rho[n_l], T[k], T[n_l], p, Tp,         sat.Lvap(p), sat.Tsat(p), sat.sigma(sat.Tsat(p), p), theta_, molar_mass_);
              double dTp_N_sites = 1.e6*(   Hibiki_Ishii_Site_density(rho[k], rho[n_l], T[k], T[n_l], p, Tp + .5e-7, sat.Lvap(p), sat.Tsat(p), sat.sigma(sat.Tsat(p), p), theta_, molar_mass_)
                                            - Hibiki_Ishii_Site_density(rho[k], rho[n_l], T[k], T[n_l], p, Tp - .5e-7, sat.Lvap(p), sat.Tsat(p), sat.sigma(sat.Tsat(p), p), theta_, molar_mass_));
              double dTl_N_sites = 1.e6*(   Hibiki_Ishii_Site_density(rho[k], rho[n_l], T[k], T[n_l] + .5e-7, p, Tp, sat.Lvap(p), sat.Tsat(p), sat.sigma(sat.Tsat(p), p), theta_, molar_mass_)
                                            - Hibiki_Ishii_Site_density(rho[k], rho[n_l], T[k], T[n_l] - .5e-7, p, Tp, sat.Lvap(p), sat.Tsat(p), sat.sigma(sat.Tsat(p), p), theta_, molar_mass_));
              double dTk_N_sites = 1.e6*(   Hibiki_Ishii_Site_density(rho[k], rho[n_l], T[k] + .5e-7, T[n_l], p, Tp, sat.Lvap(p), sat.Tsat(p), sat.sigma(sat.Tsat(p), p), theta_, molar_mass_)
                                            - Hibiki_Ishii_Site_density(rho[k], rho[n_l], T[k] - .5e-7, T[n_l], p, Tp, sat.Lvap(p), sat.Tsat(p), sat.sigma(sat.Tsat(p), p), theta_, molar_mass_));

              double u_bulk = 0;
              if (sub_type(Flux_parietal_adaptatif, correlation_monophasique_.valeur()))
                {
                  const Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, ref_cast(Pb_Multiphase, pb_.valeur()).get_correlation("Loi_paroi").valeur());
                  const double u_tau = corr_loi_paroi.get_utau(f);
                  u_bulk = 20.*u_tau; // Big approximation...
                }
              else Cerr << "Flux_parietal_Kommajosyula::qp isn't adapted to a single-phase " << correlation_monophasique_->que_suis_je() << " heat flux for now ! " , Process::exit();

              // Single phase heat transfer coefficient
              double h_fc = (*dTp_qpk)(n_l);

              // Departure diameter (page 47 Kommajosyula PhD)
              double D_d     = 18.9e-6 * std::pow( (rho[n_l]-rho[k])/rho[k], 0.27 ) *                     std::pow(Ja_sup, 0.75) *                     std::pow(1+Ja_sub, -0.3) * std::pow(u_bulk, -0.26);
              double dTp_D_d = 18.9e-6 * std::pow( (rho[n_l]-rho[k])/rho[k], 0.27 ) * 0.75 * dTp_Ja_sup * std::pow(Ja_sup,-0.25) *                     std::pow(1+Ja_sub, -0.3) * std::pow(u_bulk, -0.26);
              double dTl_D_d = 18.9e-6 * std::pow( (rho[n_l]-rho[k])/rho[k], 0.27 ) *                     std::pow(Ja_sup, 0.75) * -0.3 * dTl_Ja_sub * std::pow(1+Ja_sub, -1.3) * std::pow(u_bulk, -0.26);

              // Liftoff diameter (page 47 Kommajosyula PhD)
              double D_lo     = 1.2*D_d ;
              double dTp_D_lo = 1.2*dTp_D_d ;
              double dTl_D_lo = 1.2*dTl_D_d ;

              // Bubble growth constant (page 42 Kommajosyula PhD) excluding the microlayer contribution ; ! different formulation in the manuscript and algorithmic appendix
              double chi = 1.55 - 0.05 * Delta_T_sub     / Delta_T_sup;
              double dTp_chi =  - 0.05 * Delta_T_sub     * -dTp_Delta_T_sup / std::pow(Delta_T_sup, 2);
              double dTl_chi =  - 0.05 * dTl_Delta_T_sub / Delta_T_sup;

              double K     = chi     * 2*std::sqrt(3/M_PI) *Ja_sup    *lambda[n_l]/(rho[n_l]*Cp[n_l]);
              double dTp_K = dTp_chi * 2*std::sqrt(3/M_PI) *Ja_sup    *lambda[n_l]/(rho[n_l]*Cp[n_l])
                             + chi   * 2*std::sqrt(3/M_PI) *dTp_Ja_sup*lambda[n_l]/(rho[n_l]*Cp[n_l]);
              double dTl_K = dTl_chi * 2*std::sqrt(3/M_PI) *Ja_sup    *lambda[n_l]/(rho[n_l]*Cp[n_l]);

              // Bubble growth time (page 44 Kommajosyula PhD)
              double t_g     = D_d*D_d       / (16*K*K) ;
              double dTp_t_g = 2*dTp_D_d*D_d / (16*K*K)
                               + D_d*D_d     * -2*dTp_K / (16*K*K*K) ;
              double dTl_t_g = 2*dTl_D_d*D_d / (16*K*K)
                               + D_d*D_d     * -2*dTl_K / (16*K*K*K) ;

              // Bubble wait time (page 56 Kommajosyula PhD)
              double t_w     = 0.0061 *                std::pow(Ja_sub, 0.63) / Delta_T_sup;
              double dTp_t_w = 0.0061 *                std::pow(Ja_sub, 0.63) * -dTp_Delta_T_sup / std::pow(Delta_T_sup,2);
              double dTl_t_w = 0.0061 * .63*dTl_Ja_sub*std::pow(Ja_sub, -.37) / Delta_T_sup;

              // Bubble departure frequency (page 51 Kommajosyula PhD)
              double f_dep     =  1/(t_g+t_w) ;
              double dTp_f_dep = -(dTp_t_g+dTp_t_w)/std::pow(t_g+t_w, 2) ;
              double dTl_f_dep = -(dTl_t_g+dTl_t_w)/std::pow(t_g+t_w, 2) ;

              // Time necessary to reform the termal boundary layer (page 61 Kommajosyula PhD ; error in equation 4.5 of the manuscript)
              double t_star = lambda[n_l] * rho[n_l]*Cp[n_l]/(h_fc*h_fc*M_PI) ;

              // Active nucleation site density (page 66 Kommajosyula PhD)
              double N_0     = f_dep * t_g * M_PI * D_d*D_d/4;
              double dTp_N_0 = dTp_f_dep * t_g     * M_PI *       D_d*D_d/4
                               + f_dep   * dTp_t_g * M_PI *       D_d*D_d/4
                               + f_dep   * t_g     * M_PI * 2*dTp_D_d*D_d/4;
              double dTl_N_0 = dTl_f_dep * t_g     * M_PI *       D_d*D_d/4
                               + f_dep   * dTl_t_g * M_PI *       D_d*D_d/4
                               + f_dep   * t_g     * M_PI * 2*dTl_D_d*D_d/4;

              double N_active     = Lambert_W_function(N_0 * N_sites)/N_0;
              double dLoc_Lambert = 1e6*(Lambert_W_function(N_0*N_sites + 5e-7)-Lambert_W_function(N_0*N_sites-5e-7)); // Super elegant !
              double dTp_N_active = (dTp_N_0*N_sites + N_0*dTp_N_sites)*dLoc_Lambert / N_0
                                    +                Lambert_W_function(N_0*N_sites) * -dTp_N_0/std::pow(N_0,2);
              double dTl_N_active = (dTl_N_0*N_sites + N_0*dTl_N_sites)*dLoc_Lambert / N_0
                                    +                Lambert_W_function(N_0*N_sites) * -dTl_N_0/std::pow(N_0,2);
              double dTk_N_active =                     N_0*dTk_N_sites*dLoc_Lambert / N_0;

              // Single bubble sliding area
              double A_sl     =                  std::pow(N_active, -.5)*(D_d+D_lo)        /2;
              double dTp_A_sl = -.5*dTp_N_active*std::pow(N_active,-1.5)*(D_d+D_lo)        /2
                                +                  std::pow(N_active, -.5)*(dTp_D_d+dTp_D_lo)/2;
              double dTl_A_sl = -.5*dTl_N_active*std::pow(N_active,-1.5)*(D_d+D_lo)        /2
                                +                  std::pow(N_active, -.5)*(dTl_D_d+dTl_D_lo)/2;
              double dTk_A_sl = -.5*dTk_N_active*std::pow(N_active,-1.5)*(D_d+D_lo)        /2;

              // Percentage of the surface occupied by sliding bubbles
              double S_sl     = std::min(1., A_sl * N_active * f_dep * t_star) ;
              double dTp_S_sl = (1. <  A_sl * N_active * f_dep * t_star) ? 0 :
                                dTp_A_sl    * N_active     * f_dep     * t_star
                                + A_sl      * dTp_N_active * f_dep     * t_star
                                + A_sl      * N_active     * dTp_f_dep * t_star  ;
              double dTl_S_sl = (1. <  A_sl * N_active * f_dep * t_star) ? 0 :
                                dTl_A_sl    * N_active     * f_dep     * t_star
                                + A_sl      * dTl_N_active * f_dep     * t_star
                                + A_sl      * N_active     * dTl_f_dep * t_star  ;
              double dTk_S_sl = (1. <  A_sl * N_active * f_dep * t_star) ? 0 :
                                dTk_A_sl    * N_active     * f_dep     * t_star
                                + A_sl      * dTk_N_active * f_dep     * t_star;

              DoubleTrav qpk_loc, da_qpk_loc, dp_qpk_loc, dv_qpk_loc, dTf_qpk_loc, dTp_qpk_loc;

              if (qpk)     qpk_loc    = *qpk;
              if (da_qpk)  da_qpk_loc = *da_qpk;
              if (dp_qpk)  dp_qpk_loc = *dp_qpk;
              if (dv_qpk)  dv_qpk_loc = *dv_qpk;
              if (dTf_qpk) dTf_qpk_loc= *dTf_qpk;
              if (dTp_qpk) dTp_qpk_loc= *dTp_qpk;

              // We correct the single phase heat flux
              if (qpk)     (*qpk)    *= (1-S_sl);
              if (da_qpk)  (*da_qpk) *= (1-S_sl);
              if (dp_qpk)  (*dp_qpk) *= (1-S_sl);
              if (dv_qpk)  (*dv_qpk) *= (1-S_sl);
              if (dTf_qpk)
                {
                  (*dTf_qpk)         *= (1-S_sl);
                  (*dTf_qpk)(n_l,n_l)+= -dTl_S_sl*qpk_loc(n_l,n_l);
                  (*dTf_qpk)(n_l, k )+= -dTk_S_sl*qpk_loc(n_l, k );
                }
              if (dTp_qpk) (*dTp_qpk) = (1-S_sl); // Non !
              {
                (*dTp_qpk)         *= (1-S_sl);
                (*dTp_qpk)(n_l)    += -dTp_S_sl*qpk_loc(n_l,n_l);
              }

              // Bubble sliding
              if (qpk)     (*qpk)(n_l)          += S_sl*2*h_fc*(Tp-T[n_l]);
              if (dTf_qpk) (*dTf_qpk)(n_l, n_l) +=-S_sl*2*h_fc + dTl_S_sl*2*h_fc*(Tp-T[n_l]) ;
              if (dTf_qpk) (*dTf_qpk)(n_l, k)   +=               dTk_S_sl*2*h_fc*(Tp-T[n_l]) ;
              if (dTp_qpk) (*dTp_qpk)(n_l, k)   += S_sl*2*h_fc + dTp_S_sl*2*h_fc*(Tp-T[n_l]) ;

              // Evaporation (calculer les derivees/T apres)
              if (qpi)         (*qpi)(n_l, k)      = 1./6.*M_PI * rho[k] * sat.Lvap(p) *              std::pow(D_lo,3) *     f_dep *     N_active;
              if (dTp_qpi) (*dTp_qpi)(n_l, k)      = 1./6.*M_PI * rho[k] * sat.Lvap(p) * ( 3*dTp_D_lo*std::pow(D_lo,2) *     f_dep *     N_active
                                                                                             +        std::pow(D_lo,3) * dTp_f_dep *     N_active
                                                                                             +        std::pow(D_lo,3) *     f_dep * dTp_N_active);
              if (dTf_qpi) (*dTf_qpi)(n_l, k, n_l) = 1./6.*M_PI * rho[k] * sat.Lvap(p) * ( 3*dTl_D_lo*std::pow(D_lo,2) *     f_dep *     N_active
                                                                                             +        std::pow(D_lo,3) * dTl_f_dep *     N_active
                                                                                             +        std::pow(D_lo,3) *     f_dep * dTl_N_active);
              if (dTf_qpi) (*dTf_qpi)(n_l, k,   k) = 1./6.*M_PI * rho[k] * sat.Lvap(p) * (            std::pow(D_lo,3) *     f_dep * dTk_N_active);

              if (d_nuc) (*d_nuc)(k) = D_lo;
            }
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
