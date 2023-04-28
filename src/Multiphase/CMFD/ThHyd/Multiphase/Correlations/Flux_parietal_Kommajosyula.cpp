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
#include <Domaine_VF.h>
#include <TRUSTTrav.h>
#include <Milieu_composite.h>
#include <Saturation_base.h>

#include <math.h>

Implemente_instanciable(Flux_parietal_Kommajosyula, "Flux_parietal_Kommajosyula", Flux_parietal_base);

Sortie& Flux_parietal_Kommajosyula::printOn(Sortie& os) const { return Flux_parietal_base::printOn(os); }

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

void Flux_parietal_Kommajosyula::qp(const input_t& in, output_t& out) const
{
  // On met tout a 0 a tout hasard
  if (out.qpk)     (*out.qpk)    = 0.;
  if (out.da_qpk)  (*out.da_qpk) = 0.;
  if (out.dp_qpk)  (*out.dp_qpk) = 0.;
  if (out.dv_qpk)  (*out.dv_qpk) = 0.;
  if (out.dTf_qpk) (*out.dTf_qpk)= 0.;
  if (out.dTp_qpk) (*out.dTp_qpk)= 0.;
  if (out.qpi)     (*out.qpi)    = 0.;
  if (out.da_qpi)  (*out.da_qpi) = 0.;
  if (out.dp_qpi)  (*out.dp_qpi) = 0.;
  if (out.dv_qpi)  (*out.dv_qpi) = 0.;
  if (out.dTf_qpi) (*out.dTf_qpi)= 0.;
  if (out.dTp_qpi) (*out.dTp_qpi)= 0.;
  if (out.nonlinear) (*out.nonlinear) = 1;

  // On remplit le monophasique ; pas besoin du flux interfacial normalement
  ref_cast(Flux_parietal_base, correlation_monophasique_.valeur()).qp(in, out);

  // Ici la phase liquide est forcement la phase 0 car la correlation monophasique ne remplit que la phase 0
  int n_l = 0 ;
  const Milieu_composite& milc = ref_cast(Milieu_composite, pb_->milieu());

  for (int k = 0 ; k < in.N ; k++)
    if (n_l != k)
      if (milc.has_saturation(n_l, k))
        {
          int ind_sat = k<n_l ? ( k *(in.N-1)-( k -1)*( k )/2) + (n_l- k -1) :
                        (n_l*(in.N-1)-(n_l-1)*(n_l)/2) + ( k -n_l-1);

          double Delta_T_sup = in.Tp - in.Tsat[ind_sat]; // Wall superheat

          if (Delta_T_sup > 0) // Else : no wall superheat => no nucleation => single phase heat transfer only
            {

              double Delta_T_sub = std::max(in.Tsat[ind_sat] - in.T[n_l], 1.e-8) ; // Subcooling ; non negative for numerical reasons
              double dTp_Delta_T_sup = 1.;
              double dTl_Delta_T_sub = -1.;
              double Ja_sup = in.rho[n_l]*in.Cp[n_l]*Delta_T_sup/(in.rho[k] * in.Lvap[ind_sat]);// Superheat Jakob number
              double Ja_sub = in.rho[n_l]*in.Cp[n_l]*Delta_T_sub/(in.rho[k] * in.Lvap[ind_sat]);// Subcooling Jakob number
              double dTp_Ja_sup = in.rho[n_l]*in.Cp[n_l]/(in.rho[k] * in.Lvap[ind_sat]);
              double dTl_Ja_sub =-in.rho[n_l]*in.Cp[n_l]/(in.rho[k] * in.Lvap[ind_sat]);

              // Nucleation site density (Hibiki Ishii 2003)
              double N_sites, dTp_N_sites, dTl_N_sites, dTk_N_sites;
              N_sites     =     Hibiki_Ishii_Site_density(in.rho[k], in.rho[n_l], in.T[k], in.T[n_l], in.p, in.Tp, in.Lvap[ind_sat], in.Tsat[ind_sat], in.Sigma[ind_sat], theta_, molar_mass_);
              dTp_N_sites = dTp_Hibiki_Ishii_Site_density(in.rho[k], in.rho[n_l], in.T[k], in.T[n_l], in.p, in.Tp, in.Lvap[ind_sat], in.Tsat[ind_sat], in.Sigma[ind_sat], theta_, molar_mass_);
              dTl_N_sites = dTl_Hibiki_Ishii_Site_density(in.rho[k], in.rho[n_l], in.T[k], in.T[n_l], in.p, in.Tp, in.Lvap[ind_sat], in.Tsat[ind_sat], in.Sigma[ind_sat], theta_, molar_mass_);
              dTk_N_sites = dTk_Hibiki_Ishii_Site_density(in.rho[k], in.rho[n_l], in.T[k], in.T[n_l], in.p, in.Tp, in.Lvap[ind_sat], in.Tsat[ind_sat], in.Sigma[ind_sat], theta_, molar_mass_);

              double u_bulk = 0;
              if (sub_type(Flux_parietal_adaptatif, correlation_monophasique_.valeur()))
                {
                  const Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, ref_cast(Pb_Multiphase, pb_.valeur()).get_correlation("Loi_paroi").valeur());
                  const double u_tau = corr_loi_paroi.get_utau(in.f);
                  u_bulk = 20.*u_tau; // Big approximation...
                }
              else Cerr << "Flux_parietal_Kommajosyula::qp isn't adapted to a single-phase " << correlation_monophasique_->que_suis_je() << " heat flux for now ! " , Process::exit();

              // Single phase heat transfer coefficient
              double h_fc = (*out.dTp_qpk)(n_l);

              // Departure diameter (page 47 Kommajosyula PhD)
              double D_d     = 18.9e-6 * std::pow( (in.rho[n_l]-in.rho[k])/in.rho[k], 0.27 ) *                     std::pow(Ja_sup, 0.75) *                     std::pow(1+Ja_sub, -0.3) * std::pow(u_bulk, -0.26);
              double dTp_D_d = 18.9e-6 * std::pow( (in.rho[n_l]-in.rho[k])/in.rho[k], 0.27 ) * 0.75 * dTp_Ja_sup * std::pow(Ja_sup,-0.25) *                     std::pow(1+Ja_sub, -0.3) * std::pow(u_bulk, -0.26);
              double dTl_D_d = 18.9e-6 * std::pow( (in.rho[n_l]-in.rho[k])/in.rho[k], 0.27 ) *                     std::pow(Ja_sup, 0.75) * -0.3 * dTl_Ja_sub * std::pow(1+Ja_sub, -1.3) * std::pow(u_bulk, -0.26);

              // Liftoff diameter (page 47 Kommajosyula PhD)
              double D_lo     = 1.2*D_d ;
              double dTp_D_lo = 1.2*dTp_D_d ;
              double dTl_D_lo = 1.2*dTl_D_d ;

              // Bubble growth constant (page 42 Kommajosyula PhD) excluding the microlayer contribution ; ! different formulation in the manuscript and algorithmic appendix
              double chi = 1.55 - 0.05 * Delta_T_sub     / Delta_T_sup;
              double dTp_chi =  - 0.05 * Delta_T_sub     * -dTp_Delta_T_sup / std::pow(Delta_T_sup, 2);
              double dTl_chi =  - 0.05 * dTl_Delta_T_sub / Delta_T_sup;

              double K     = chi     * 2.*std::sqrt(3./M_PI) *Ja_sup    *in.lambda[n_l]/(in.rho[n_l]*in.Cp[n_l]);
              double dTp_K = dTp_chi * 2.*std::sqrt(3./M_PI) *Ja_sup    *in.lambda[n_l]/(in.rho[n_l]*in.Cp[n_l])
                             + chi   * 2.*std::sqrt(3./M_PI) *dTp_Ja_sup*in.lambda[n_l]/(in.rho[n_l]*in.Cp[n_l]);
              double dTl_K = dTl_chi * 2.*std::sqrt(3./M_PI) *Ja_sup    *in.lambda[n_l]/(in.rho[n_l]*in.Cp[n_l]);

              // Bubble growth time (page 44 Kommajosyula PhD)
              double t_g     = D_d*D_d       / (16.*K*K) ;
              double dTp_t_g = 2*dTp_D_d*D_d / (16.*K*K)
                               + D_d*D_d     * -2.*dTp_K / (16.*K*K*K) ;
              double dTl_t_g = 2*dTl_D_d*D_d / (16*K*K)
                               + D_d*D_d     * -2.*dTl_K / (16.*K*K*K) ;

              // Bubble wait time (page 56 Kommajosyula PhD)
              double t_w     = 0.0061 *                std::pow(Ja_sub, 0.63) / Delta_T_sup;
              double dTp_t_w = 0.0061 *                std::pow(Ja_sub, 0.63) * -dTp_Delta_T_sup / std::pow(Delta_T_sup,2);
              double dTl_t_w = 0.0061 * .63*dTl_Ja_sub*std::pow(Ja_sub, -.37) / Delta_T_sup;

              // Bubble departure frequency (page 51 Kommajosyula PhD)
              double f_dep     =  1/(t_g+t_w) ;
              double dTp_f_dep = -(dTp_t_g+dTp_t_w)/std::pow(t_g+t_w, 2) ;
              double dTl_f_dep = -(dTl_t_g+dTl_t_w)/std::pow(t_g+t_w, 2) ;

              // Time necessary to reform the termal boundary layer (page 61 Kommajosyula PhD ; error in equation 4.5 of the manuscript)
              double t_star = in.lambda[n_l] * in.rho[n_l]*in.Cp[n_l]/(h_fc*h_fc*M_PI) ;

              // Active nucleation site density (page 66 Kommajosyula PhD)
              double N_0     =   f_dep   * t_g     * M_PI *        D_d*D_d/4.;
              double dTp_N_0 = dTp_f_dep * t_g     * M_PI *        D_d*D_d/4.
                               + f_dep   * dTp_t_g * M_PI *        D_d*D_d/4.
                               + f_dep   * t_g     * M_PI * 2.*dTp_D_d*D_d/4.;
              double dTl_N_0 = dTl_f_dep * t_g     * M_PI *        D_d*D_d/4.
                               + f_dep   * dTl_t_g * M_PI *        D_d*D_d/4.
                               + f_dep   * t_g     * M_PI * 2.*dTl_D_d*D_d/4.;

              double N_active     =    Lambert_W_function(N_0 * N_sites)/N_0;
              double dLoc_Lambert = dx_Lambert_W_function(N_0 * N_sites);
              double dTp_N_active = (dTp_N_0*N_sites + N_0*dTp_N_sites)*dLoc_Lambert / N_0
                                    +                Lambert_W_function(N_0*N_sites) * -dTp_N_0/std::pow(N_0,2);
              double dTl_N_active = (dTl_N_0*N_sites + N_0*dTl_N_sites)*dLoc_Lambert / N_0
                                    +                Lambert_W_function(N_0*N_sites) * -dTl_N_0/std::pow(N_0,2);
              double dTk_N_active =                     N_0*dTk_N_sites*dLoc_Lambert / N_0;

              if (N_active < 0) Cerr << "nsite " << N_sites << "N0 " << N_0 ;

              // Single bubble sliding area
              double A_sl     =                  std::pow(N_active, -.5)*(D_d+D_lo)        /2.;
              double dTp_A_sl = -.5*dTp_N_active*std::pow(N_active,-1.5)*(D_d+D_lo)        /2.
                                +                std::pow(N_active, -.5)*(dTp_D_d+dTp_D_lo)/2.;
              double dTl_A_sl = -.5*dTl_N_active*std::pow(N_active,-1.5)*(D_d+D_lo)        /2.
                                +                std::pow(N_active, -.5)*(dTl_D_d+dTl_D_lo)/2.;
              double dTk_A_sl = -.5*dTk_N_active*std::pow(N_active,-1.5)*(D_d+D_lo)        /2.;

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

              if (out.qpk)     qpk_loc    = *out.qpk;
              if (out.da_qpk)  da_qpk_loc = *out.da_qpk;
              if (out.dp_qpk)  dp_qpk_loc = *out.dp_qpk;
              if (out.dv_qpk)  dv_qpk_loc = *out.dv_qpk;
              if (out.dTf_qpk) dTf_qpk_loc= *out.dTf_qpk;
              if (out.dTp_qpk) dTp_qpk_loc= *out.dTp_qpk;

              // We correct the single phase heat flux
              if (out.qpk)     (*out.qpk)    *= (1-S_sl);
              if (out.da_qpk)  (*out.da_qpk) *= (1-S_sl);
              if (out.dp_qpk)  (*out.dp_qpk) *= (1-S_sl);
              if (out.dv_qpk)  (*out.dv_qpk) *= (1-S_sl);
              if (out.dTf_qpk)
                {
                  (*out.dTf_qpk)         *= (1-S_sl);
                  (*out.dTf_qpk)(n_l,n_l)+= -dTl_S_sl*qpk_loc(n_l);
                  (*out.dTf_qpk)(n_l, k )+= -dTk_S_sl*qpk_loc( k );
                }
              if (out.dTp_qpk)
                {
                  (*out.dTp_qpk)         *= (1-S_sl);
                  (*out.dTp_qpk)(n_l)    += -dTp_S_sl*qpk_loc(n_l);
                }

              // Bubble sliding
              if (out.qpk)     (*out.qpk)(n_l)          += S_sl*2.*h_fc*(in.Tp-in.T[n_l]);
              if (out.dTf_qpk) (*out.dTf_qpk)(n_l, n_l) +=-S_sl*2.*h_fc + dTl_S_sl*2.*h_fc*(in.Tp-in.T[n_l]) ;
              if (out.dTf_qpk) (*out.dTf_qpk)(n_l, k)   +=                dTk_S_sl*2.*h_fc*(in.Tp-in.T[n_l]) ;
              if (out.dTp_qpk) (*out.dTp_qpk)( k )      += S_sl*2.*h_fc + dTp_S_sl*2.*h_fc*(in.Tp-in.T[n_l]) ;

              // Evaporation (calculer les derivees/T apres)
              if (out.qpi)         (*out.qpi)(n_l, k)      = 1./6.*M_PI * in.rho[k] * in.Lvap[ind_sat] *              std::pow(D_lo,3.) *     f_dep *     N_active;
              if (out.dTp_qpi) (*out.dTp_qpi)(n_l, k)      = 1./6.*M_PI * in.rho[k] * in.Lvap[ind_sat] * (3.*dTp_D_lo*std::pow(D_lo,2.) *     f_dep *     N_active
                                                                                                            +        std::pow(D_lo,3.) * dTp_f_dep *     N_active
                                                                                                            +        std::pow(D_lo,3.) *     f_dep * dTp_N_active);
              if (out.dTf_qpi) (*out.dTf_qpi)(n_l, k, n_l) = 1./6.*M_PI * in.rho[k] * in.Lvap[ind_sat] * (3.*dTl_D_lo*std::pow(D_lo,2.) *     f_dep *     N_active
                                                                                                            +        std::pow(D_lo,3.) * dTl_f_dep *     N_active
                                                                                                            +        std::pow(D_lo,3.) *     f_dep * dTl_N_active);
              if (out.dTf_qpi) (*out.dTf_qpi)(n_l, k,   k) = 1./6.*M_PI * in.rho[k] * in.Lvap[ind_sat] * (            std::pow(D_lo,3.) *     f_dep * dTk_N_active);

              if (out.d_nuc) (*out.d_nuc)(k) = D_lo;

            }
        }
}

double Flux_parietal_Kommajosyula::Lambert_W_function(double x)  const
{
  if (x < 1/M_E) return x;
  else if (x > M_E) return std::log(x) - std::log(std::log(x));
  else return 0.2689*x+0.2690;
}

double Flux_parietal_Kommajosyula::dx_Lambert_W_function(double x)  const
{
  if (x < 1/M_E) return 1.;
  else if (x > M_E) return 1./x - 1./(x*std::log(x));
  else return 0.2689;
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

  double inv_Rc = p * (-1 + std::exp(L_vap*std::max(Tp-T_sat, 0.)/( R * Tp * T_sat ))) / (2 * sigma * (1+rho_v/rho_l)); // T_p is put here rather than T_v as we need T_v at the wall, not in first element.

  return N_bar * (1 - std::exp(-theta*theta/(8.*mu2))) * (-1. + std::exp(f_rho_plus * lambda_prime * inv_Rc));
}

double Flux_parietal_Kommajosyula::dTp_Hibiki_Ishii_Site_density(double rho_v, double rho_l, double T_v, double T_l, double p, double Tp, double L_vap, double T_sat, double sigma, double theta, double molar_mass) const
{
  double N_bar = 4.72e5;
  double mu2 = 0.722*180/M_PI;
  mu2 *= mu2;

  double rho_p = std::log( (rho_l-rho_v)/rho_v);
  double f_rho_plus = -0.01064 + 0.48246*rho_p - 0.22712*rho_p*rho_p + 0.05468*rho_p*rho_p*rho_p;

  double lambda_prime = 2.5e-6;
  double R = 8.314462618 / molar_mass ;

  double inv_Rc = p * (-1 + std::exp(L_vap*std::max(Tp-T_sat, 0.)/( R * Tp * T_sat ))) / (2. * sigma * (1+rho_v/rho_l)); // T_p is put here rather than T_v as we need T_v at the wall, not in first element.
  double dTp_inv_Rc = 0.;
  if (Tp-T_sat> 0.)
    dTp_inv_Rc = p * L_vap / ( R * Tp*Tp ) * std::exp(L_vap*(Tp-T_sat)/(R*Tp*T_sat)) / (2. * sigma * (1+rho_v/rho_l)); // T_p is put here rather than T_v as we need T_v at the wall, not in first element.

  return N_bar * (1 - std::exp(-theta*theta/(8.*mu2))) * (0. + f_rho_plus * lambda_prime * dTp_inv_Rc* std::exp(f_rho_plus * lambda_prime * inv_Rc));
}

double Flux_parietal_Kommajosyula::dTl_Hibiki_Ishii_Site_density(double rho_v, double rho_l, double T_v, double T_l, double p, double Tp, double L_vap, double T_sat, double sigma, double theta, double molar_mass) const
{
  return 0.;
}

double Flux_parietal_Kommajosyula::dTk_Hibiki_Ishii_Site_density(double rho_v, double rho_l, double T_v, double T_l, double p, double Tp, double L_vap, double T_sat, double sigma, double theta, double molar_mass) const
{
  return 0.;
}
