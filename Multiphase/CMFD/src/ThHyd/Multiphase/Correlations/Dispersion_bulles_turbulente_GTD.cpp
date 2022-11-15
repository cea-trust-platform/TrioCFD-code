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
// File:        Dispersion_bulles_turbulente_GTD.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Dispersion_bulles_turbulente_GTD.h>
#include <Pb_Multiphase.h>
#include <TRUSTTrav.h>
#include <Frottement_interfacial_base.h>
#include <Masse_ajoutee_base.h>
#include <math.h>

Implemente_instanciable(Dispersion_bulles_turbulente_GTD, "Dispersion_bulles_turbulente_GTD", Dispersion_bulles_base);

Sortie& Dispersion_bulles_turbulente_GTD::printOn(Sortie& os) const
{
  return os;
}

Entree& Dispersion_bulles_turbulente_GTD::readOn(Entree& is)
{
  const Pb_Multiphase *pbm = sub_type(Pb_Multiphase, pb_.valeur()) ? &ref_cast(Pb_Multiphase, pb_.valeur()) : NULL;

  if (!pbm || pbm->nb_phases() == 1) Process::exit(que_suis_je() + " : not needed for single-phase flow!");
  for (int n = 0; n < pbm->nb_phases(); n++) //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    if (pbm->nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))  n_l = n;
  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");

  if (pbm->has_correlation("frottement_interfacial")) correlation_drag_ = pbm->get_correlation("frottement_interfacial"); //correlation fournie par le bloc correlation
  else correlation_drag_.typer_lire(*pbm, "frottement_interfacial", is); //sinon -> on la lit

  if (pbm->has_correlation("masse_ajoutee")) correlation_drag_ = pbm->get_correlation("masse_ajoutee"); //correlation fournie par le bloc correlation
  else correlation_drag_.typer_lire(*pbm, "masse_ajoutee", is); //sinon -> on la lit

  return is;
}


void Dispersion_bulles_turbulente_GTD::coefficient( const DoubleTab& alpha, const DoubleTab& p, const DoubleTab& T,
                                                    const DoubleTab& rho, const DoubleTab& mu, const DoubleTab& sigma,
                                                    const DoubleTab& nut, const DoubleTab& k_turb, const DoubleTab& d_bulles,
                                                    const DoubleTab& ndv, DoubleTab& coeff) const
{
  int N = ndv.dimension(0);

  const Frottement_interfacial_base& corr_drag = ref_cast(Frottement_interfacial_base, correlation_drag_.valeur());
  DoubleTrav coeff_CD(N, N);
  corr_drag.coefficient_CD( alpha, p, T, rho, mu, sigma, 0., ndv, d_bulles, coeff_CD);

  const Masse_ajoutee_base& corr_MA = ref_cast(Masse_ajoutee_base, correlation_MA_.valeur());
  DoubleTrav coeff_MA(N);
  corr_MA.coefficient( &alpha(0), &rho(0), coeff_MA);

  for (int k = 0; k < N; k++)
    if (k!=n_l)
      for (int i = 0 ; i<2 ; i++) // k gas phase
        {
          double f_D = 3./4.*coeff_CD(k, n_l)*ndv(n_l,k)/d_bulles(k);
          double t_12_f = 1/f_D * (rho(k)/rho(n_l) + coeff_MA(k)) ;
          double t_12_t = 3./2. * nut(n_l)/k_turb(n_l) * std::pow(1+2.7 * ndv(n_l,k)*ndv(n_l,k)/k_turb(n_l), .5);
          double b = rho(n_l) * (1+coeff_MA(k)) / (rho(k) + coeff_MA(k)*rho(n_l));
          double eta_loc = t_12_t / t_12_f ;
          double GTD = (f_D * t_12_t - 1)*(b + eta_loc)/(1+eta_loc) + coeff_MA(k)*(b*b + eta_loc)/(1+eta_loc);
          coeff(k, n_l, i) = GTD * rho(n_l)*k_turb(n_l) ;
        }
}