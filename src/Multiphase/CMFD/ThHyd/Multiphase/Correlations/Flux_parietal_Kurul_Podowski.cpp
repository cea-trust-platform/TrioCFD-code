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
// File:        Flux_parietal_Kurul_Podowski.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Flux_parietal_Kurul_Podowski.h>
#include <Flux_parietal_adaptatif.h>
#include <Loi_paroi_adaptative.h>
#include <Pb_Multiphase.h>
#include <Domaine_dis.h>
#include <TRUSTTrav.h>
#include <Milieu_composite.h>
#include <Saturation_base.h>
#include <math.h>

Implemente_instanciable(Flux_parietal_Kurul_Podowski, "Flux_parietal_Kurul_Podowski", Flux_parietal_base);

Sortie& Flux_parietal_Kurul_Podowski::printOn(Sortie& os) const { return Flux_parietal_base::printOn(os); }

Entree& Flux_parietal_Kurul_Podowski::readOn(Entree& is)
{
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, pb_.valeur());
  Correlation_base::typer_lire_correlation(correlation_monophasique_, pbm, "Flux_parietal", is);
  Cout << que_suis_je() << " : single-phase wall heat flux is " << correlation_monophasique_->que_suis_je() << finl;

  Param param(que_suis_je());
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

void Flux_parietal_Kurul_Podowski::completer()
{
  correlation_monophasique_->completer();
}

void Flux_parietal_Kurul_Podowski::qp(const input_t& in, output_t& out) const
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
          // on part de ind_sat = (k*(in.N-1)-(k-1)*(k)/2) + (l-k-1) // avec k<l

          int ind_sat = k<n_l ? ( k *(in.N-1)-( k -1)*( k )/2) + (n_l- k -1) :
                        (n_l*(in.N-1)-(n_l-1)*(n_l)/2) + ( k -n_l-1);

          if (in.Tp - in.Tsat[ind_sat] > 0) // Else : no wall superheat => no nucleation => single phase heat transfer only
            {

              double dd = 1.e-4*(in.Tp - in.Tsat[ind_sat])+0.0014;
              double dTp_dd = 1.e-4;
              // NB: dans Neptune, on utilise un correlation plus complexe:
              /*double b = (sat.Tsat(in.p)-in.T[n_l])/(2-2*in.rho[k]/in.rho[n_l])
              double a = (in.Tp - sat.Tsat(in.p))*in.lambda[n_l]/(2*in.rho[k]*)
              double phi = std::max(1., std::pow(v[n_l]/.61, .47))
              double dd = 2.42e-5*std::pow(in.p, 0.709)*a/std::sqrt(b*phi);
              // NB: dans Neptune, le b est un peu different mais il faudrait un newton pour le calculer...
              */

              double     N_sites =          std::pow(210.*(in.Tp - in.Tsat[ind_sat]), 1.8);
              double dTp_N_sites = 210.*1.8*std::pow(210.*(in.Tp - in.Tsat[ind_sat]),  .8);

              double     A_bubbles = std::min(1., 3.1415/4.*N_sites*dd*dd);
              double dTp_A_bubbles = (3.1415/4.*N_sites*dd*dd>1) ? 0. : 3.1415/4.*dTp_N_sites*dd*dd+3.1415/4.*N_sites*2.*dTp_dd*dd ;

              double     f_dep = std::sqrt(4./3*9.81*(in.rho[n_l]-in.rho[k])/(in.rho[n_l]))   *            std::pow(dd, -0.5);
              double dTp_f_dep = std::sqrt(4./3*9.81*(in.rho[n_l]-in.rho[k])/(in.rho[n_l]))   *-0.5*dTp_dd*std::pow(dd, -1.5);

              double   q_evap  =       f_dep * 3.1415/6. * dd*dd*dd       * in.rho[k] * in.Lvap[ind_sat] * N_sites ;
              double dTp_q_evap=   dTp_f_dep * 3.1415/6. * dd*dd*dd       * in.rho[k] * in.Lvap[ind_sat] * N_sites
                                   +f_dep    * 3.1415/6. *3.*dTp_dd*dd*dd * in.rho[k] * in.Lvap[ind_sat] * N_sites
                                   +f_dep    * 3.1415/6. * dd*dd*dd       * in.rho[k] * in.Lvap[ind_sat] * dTp_N_sites ;


              double     q_quench= A_bubbles     *2. * in.lambda[n_l] * (in.Tp - in.T[n_l]) / std::sqrt(3.1415*in.lambda[n_l]/(in.rho[n_l]*in.Cp[n_l])) * std::sqrt(f_dep);
              double dTl_q_quench= A_bubbles     *2. * in.lambda[n_l] * (-1.)               / std::sqrt(3.1415*in.lambda[n_l]/(in.rho[n_l]*in.Cp[n_l])) * std::sqrt(f_dep);
              double dTp_q_quench= dTp_A_bubbles *2. * in.lambda[n_l] * (in.Tp - in.T[n_l]) / std::sqrt(3.1415*in.lambda[n_l]/(in.rho[n_l]*in.Cp[n_l])) * std::sqrt(f_dep)
                                   +A_bubbles    *2. * in.lambda[n_l] * (1.)                / std::sqrt(3.1415*in.lambda[n_l]/(in.rho[n_l]*in.Cp[n_l])) * std::sqrt(f_dep)
                                   +A_bubbles    *2. * in.lambda[n_l] * (in.Tp - in.T[n_l]) / std::sqrt(3.1415*in.lambda[n_l]/(in.rho[n_l]*in.Cp[n_l])) * .5*dTp_f_dep/std::sqrt(f_dep);

              // We correct the single phase heat flux
              double qpk_nl_loc = -1;
              if (out.qpk)        qpk_nl_loc  = (*out.qpk)(n_l);
              if (out.qpk)     (*out.qpk)    *= (1-A_bubbles);
              if (out.da_qpk)  (*out.da_qpk) *= (1-A_bubbles);
              if (out.dp_qpk)  (*out.dp_qpk) *= (1-A_bubbles);
              if (out.dv_qpk)  (*out.dv_qpk) *= (1-A_bubbles);
              if (out.dTf_qpk) (*out.dTf_qpk)*= (1-A_bubbles);
              if (out.dTp_qpk)
                {
                  (*out.dTp_qpk)         *= (1-A_bubbles);
                  (*out.dTp_qpk)(n_l)    += -dTp_A_bubbles*qpk_nl_loc;
                }

              // Quenching
              if (out.qpk)     (*out.qpk)(n_l)        += q_quench;
              if (out.dTp_qpk) (*out.dTp_qpk)(n_l)    += dTp_q_quench;
              if (out.dTf_qpk) (*out.dTf_qpk)(n_l,n_l)+= dTl_q_quench;

              // Evaporation
              if (out.qpi)     (*out.qpi)(n_l, k)     += q_evap;
              if (out.dTp_qpi) (*out.dTp_qpi)(n_l, k) += dTp_q_evap;

            }
        }
}

