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
// File:        Flux_parietal_adaptatif.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Flux_parietal_adaptatif.h>
#include <Loi_paroi_adaptative.h>
#include <Pb_Multiphase.h>
#include <Correlation_base.h>

#include <math.h>

Implemente_instanciable(Flux_parietal_adaptatif, "Flux_parietal_adaptatif", Flux_parietal_base);

Sortie& Flux_parietal_adaptatif::printOn(Sortie& os) const { return Flux_parietal_base::printOn(os); }
Entree& Flux_parietal_adaptatif::readOn(Entree& is) { return Flux_parietal_base::readOn(is); }

void Flux_parietal_adaptatif::completer()
{
  correlation_loi_paroi_ = pb_->get_correlation("Loi_paroi");
}

void Flux_parietal_adaptatif::qp(const input_t& in, output_t& out) const
{
  const Loi_paroi_base& corr_loi_paroi = ref_cast(Loi_paroi_base, correlation_loi_paroi_->valeur());
  const double u_tau = corr_loi_paroi.get_utau(in.f);

  double theta_plus = calc_theta_plus(in.y, u_tau, in.mu[0], in.lambda[0], in.rho[0], in.Cp[0], in.D_h),
         fac = in.rho[0] * in.Cp[0] * u_tau / theta_plus ;
  if (out.qpk) (*out.qpk)(0) = fac * (in.Tp - in.T[0]);
  if (out.dTf_qpk) (*out.dTf_qpk)(0,0) = -fac;
  if (out.dTp_qpk) (*out.dTp_qpk)(0)   = fac;
}

double Flux_parietal_adaptatif::calc_theta_plus(double y, double u_tau, double mu, double lambda, double rho, double Cp,double Diam_hyd_) const
{
  double Prandtl = mu*Cp/lambda;
  double visc = mu/rho;
  double y_plus = y * u_tau/visc;
  double beta = std::pow(3.85*std::pow(Prandtl, 1./3)-1.3, 2) + 2.12*std::log(Prandtl);
  double gamma = 0.01*(Prandtl*y_plus)*(Prandtl*y_plus)*(Prandtl*y_plus)*(Prandtl*y_plus)/(1+5*Prandtl*Prandtl*Prandtl*y_plus);
  double y_on_D_h = 0 ; //(Diam_hyd_>0)? y/Diam_hyd_:0;
  return Prandtl*y_plus*std::exp(-gamma) + (2.12*std::log((1+y_plus)*1.5*(2-y_on_D_h)/(1+2*(1-y_on_D_h)*(1-y_on_D_h)))+beta)*std::exp(-1/gamma);
}

