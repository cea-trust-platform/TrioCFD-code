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
#include <Correlation.h>
#include <Pb_Multiphase.h>
#include <Domaine_dis.h>
#include <Zone_VF.h>
#include <Op_Diff_PolyMAC_base.h>

#include <math.h>

Implemente_instanciable(Flux_parietal_adaptatif, "Flux_parietal_adaptatif", Flux_parietal_base);

Sortie& Flux_parietal_adaptatif::printOn(Sortie& os) const
{
  return os;
}

Entree& Flux_parietal_adaptatif::readOn(Entree& is)
{
  return is;
}

void Flux_parietal_adaptatif::completer()
{
  correlation_loi_paroi_ = ref_cast(Pb_Multiphase, pb_.valeur()).get_correlation("Loi_paroi");
}

void Flux_parietal_adaptatif::qp(int N, int f, double D_h, double D_ch,
                                 const double *alpha, const double *T, const double p, const double *v, const double Tp,
                                 const double *lambda, const double *mu, const double *rho, const double *Cp,
                                 DoubleTab *qpk, DoubleTab *da_qpk, DoubleTab *dp_qpk, DoubleTab *dv_qpk, DoubleTab *dTf_qpk, DoubleTab *dTp_qpk,
                                 DoubleTab *qpi, DoubleTab *da_qpi, DoubleTab *dp_qpi, DoubleTab *dv_qpi, DoubleTab *dTf_qpi, DoubleTab *dTp_qpi,
                                 DoubleTab *d_nuc, int& nonlinear) const
{
  const Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, correlation_loi_paroi_.valeur().valeur());
  const double y = corr_loi_paroi.get_y(f);
  const double u_tau = corr_loi_paroi.get_utau(f);
  int e = ref_cast(Zone_VF, pb_->domaine_dis().zone_dis(0).valeur()).face_voisins()(f,0);
  const DoubleTab& diffu_th = ref_cast(Op_Diff_PolyMAC_base, pb_->equation(2).operateur(0).l_op_base()).nu();

  double theta_plus = calc_theta_plus(y, u_tau, mu[0], lambda[0], rho[0], Cp[0], D_h),
         fac = alpha[0] *rho[0] *Cp[0] * u_tau / theta_plus   *   1./(1-rho[0]*Cp[0]*y/ diffu_th(e, 0) * u_tau/theta_plus);
  if (qpk) (*qpk)(0) = fac * (Tp - T[0]);
  if (dTf_qpk) (*dTf_qpk)(0,0) = -fac;
  if (dTp_qpk) (*dTp_qpk)(0)   = fac;
  return;
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

