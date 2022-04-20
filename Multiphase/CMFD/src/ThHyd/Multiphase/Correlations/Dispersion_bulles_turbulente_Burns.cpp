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
#include <DoubleTrav.h>
#include <Frottement_interfacial_base.h>
#include <math.h>

Implemente_instanciable(Dispersion_bulles_turbulente_Burns, "Dispersion_turbulente_Burns", Dispersion_bulles_base);

Sortie& Dispersion_bulles_turbulente_Burns::printOn(Sortie& os) const
{
  return os;
}

Entree& Dispersion_bulles_turbulente_Burns::readOn(Entree& is)
{
  const Pb_Multiphase *pbm = sub_type(Pb_Multiphase, pb_.valeur()) ? &ref_cast(Pb_Multiphase, pb_.valeur()) : NULL;

  if (!pbm || pbm->nb_phases() == 1) Process::exit(que_suis_je() + " : not needed for single-phase flow!");
  for (int n = 0; n < pbm->nb_phases(); n++) //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    if (pbm->nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))  n_l = n;
  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");

  if (pbm->has_correlation("frottement_interfacial")) correlation_drag_ = pbm->get_correlation("frottement_interfacial"); //correlation fournie par le bloc correlation
  else correlation_drag_.typer_lire(*pbm, "frottement_interfacial", is); //sinon -> on la lit

  return is;
}


void Dispersion_bulles_turbulente_Burns::coefficient( const DoubleTab& alpha, const DoubleTab& p, const DoubleTab& T,
                                                      const DoubleTab& rho, const DoubleTab& mu, const DoubleTab& sigma,
                                                      const DoubleTab& nut, const DoubleTab& k_turb, const DoubleTab& d_bulles,
                                                      const DoubleTab& ndv, DoubleTab& coeff) const
{
  const Frottement_interfacial_base& corr = ref_cast(Frottement_interfacial_base, correlation_drag_);
  int N = ndv.dimension(0);

  DoubleTrav coeff_drag(N, N, 2);
  corr.coefficient( alpha, p, T, rho, mu, sigma, 0, ndv, d_bulles, coeff_drag);

  coeff = 0;

  for (int k = 0; k < N; k++) if (k!=n_l) for (int i = 0 ; i<2 ; i++) // k gas phase
        {
          coeff(k, n_l, i) = (alpha(k) > 1.e-6) ? nut(n_l)/Prt_ * coeff_drag(k, n_l, i)/alpha(k)  : 0 ;
          coeff(n_l, k, i) = (alpha(n_l)>1.e-6) ? nut(n_l)/Prt_ * coeff_drag(n_l, k, i)/alpha(n_l): 0 ;
        }
}