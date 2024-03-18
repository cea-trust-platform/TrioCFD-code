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

#include <Dispersion_bulles_turbulente_Bertodano.h>
#include <Pb_Multiphase.h>
#include <TRUSTTrav.h>
#include <Frottement_interfacial_base.h>
#include <math.h>

Implemente_instanciable(Dispersion_bulles_turbulente_Bertodano, "Dispersion_bulles_turbulente_Bertodano", Dispersion_bulles_base);

Sortie& Dispersion_bulles_turbulente_Bertodano::printOn(Sortie& os) const
{
  return os;
}

Entree& Dispersion_bulles_turbulente_Bertodano::readOn(Entree& is)
{
  const Pb_Multiphase *pbm = sub_type(Pb_Multiphase, pb_.valeur()) ? &ref_cast(Pb_Multiphase, pb_.valeur()) : nullptr;

  if (!pbm || pbm->nb_phases() == 1) Process::exit(que_suis_je() + " : not needed for single-phase flow!");
  for (int n = 0; n < pbm->nb_phases(); n++) //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    if (pbm->nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))  n_l = n;
  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");

  if (pbm->has_correlation("frottement_interfacial")) correlation_drag_ = pbm->get_correlation("frottement_interfacial"); //correlation fournie par le bloc correlation
  else correlation_drag_.typer_lire(*pbm, "frottement_interfacial", is); //sinon -> on la lit

  return is;
}


void Dispersion_bulles_turbulente_Bertodano::coefficient(const input_t& in, output_t& out) const
{
  const Frottement_interfacial_base& corr = ref_cast(Frottement_interfacial_base, correlation_drag_.valeur());
  int N = out.Ctd.dimension(0);

  DoubleTrav coeff_CD(N, N);
  corr.coefficient_CD( in.alpha, in.p, in.T, in.rho, in.mu, in.sigma, in.dh, in.nv, in.d_bulles, coeff_CD);

  for (int k = 0; k < N; k++)
    if (k!=n_l)
      {
        double t_c = in.nut[n_l] / in.k_turb[n_l] ;
        double t_d = 4./3.*in.rho[k]*in.d_bulles[k]/(coeff_CD(k, n_l)*in.rho[n_l]*in.nv(n_l,k));
        out.Ctd(k, n_l) = 2.*in.rho[n_l]*in.k_turb[n_l] * t_c*t_c/(t_d*(t_c+t_d));
      }
}
