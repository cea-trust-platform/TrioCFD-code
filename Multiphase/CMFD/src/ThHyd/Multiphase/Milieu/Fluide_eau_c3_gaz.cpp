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
// File:        Fluide_eau_c3_gaz.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Milieu
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Fluide_eau_c3_gaz.h>
#include <Lois_eau_c3.h>

Implemente_instanciable(Fluide_eau_c3_gaz, "Fluide_eau_c3_gaz", Fluide_reel_base);

Sortie& Fluide_eau_c3_gaz::printOn(Sortie& os) const
{
  return os;
}

Entree& Fluide_eau_c3_gaz::readOn(Entree& is)
{
  Fluide_reel_base::readOn(is);
  return is;
}

double Fluide_eau_c3_gaz::rho_(const double T, const double P) const
{
  /* calcul a saturation */
  int un = 1, ienc = 0, ier, itest;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  /* calcul a T */
  double vapa, vapb, vapc, vapdb, vapdc, hg, dP_hg, cpg, dP_cpg, dT_cpg, rhog, dP_rhog, dT_rhog, delta_h;
  F77NAME(FTVAP)(&un, &ienc, &ier, &itest, &P, &T, &Ts, &dP_Ts, &hgs, &dP_hgs, &vapa, &vapb, &vapc, &vapdb, &vapdc,
                 &hg, &dP_hg, &cpg, &dP_cpg, &dT_cpg, &rhog, &dP_rhog, &dT_rhog, &delta_h);
  return rhog;
}

double Fluide_eau_c3_gaz::dT_rho_(const double T, const double P) const
{
  /* calcul a saturation */
  int un = 1, ienc = 0, ier, itest;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  /* calcul a T */
  double vapa, vapb, vapc, vapdb, vapdc, hg, dP_hg, cpg, dP_cpg, dT_cpg, rhog, dP_rhog, dT_rhog, delta_h;
  F77NAME(FTVAP)(&un, &ienc, &ier, &itest, &P, &T, &Ts, &dP_Ts, &hgs, &dP_hgs, &vapa, &vapb, &vapc, &vapdb, &vapdc,
                 &hg, &dP_hg, &cpg, &dP_cpg, &dT_cpg, &rhog, &dP_rhog, &dT_rhog, &delta_h);
  return dT_rhog;
}

double Fluide_eau_c3_gaz::dP_rho_(const double T, const double P) const
{
  /* calcul a saturation */
  int un = 1, ienc = 0, ier, itest;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  /* calcul a T */
  double vapa, vapb, vapc, vapdb, vapdc, hg, dP_hg, cpg, dP_cpg, dT_cpg, rhog, dP_rhog, dT_rhog, delta_h;
  F77NAME(FTVAP)(&un, &ienc, &ier, &itest, &P, &T, &Ts, &dP_Ts, &hgs, &dP_hgs, &vapa, &vapb, &vapc, &vapdb, &vapdc,
                 &hg, &dP_hg, &cpg, &dP_cpg, &dT_cpg, &rhog, &dP_rhog, &dT_rhog, &delta_h);
  return dP_rhog;
}

double Fluide_eau_c3_gaz::h_(const double T, const double P) const
{
  /* calcul a saturation */
  int un = 1, ienc = 0, ier, itest;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  /* calcul a T */
  double vapa, vapb, vapc, vapdb, vapdc, hg, dP_hg, cpg, dP_cpg, dT_cpg, rhog, dP_rhog, dT_rhog, delta_h;
  F77NAME(FTVAP)(&un, &ienc, &ier, &itest, &P, &T, &Ts, &dP_Ts, &hgs, &dP_hgs, &vapa, &vapb, &vapc, &vapdb, &vapdc,
                 &hg, &dP_hg, &cpg, &dP_cpg, &dT_cpg, &rhog, &dP_rhog, &dT_rhog, &delta_h);
  return hg;
}

double Fluide_eau_c3_gaz::dT_h_(const double T, const double P) const
{
  /* calcul a saturation */
  int un = 1, ienc = 0, ier, itest;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  /* calcul a T */
  double vapa, vapb, vapc, vapdb, vapdc, hg, dP_hg, cpg, dP_cpg, dT_cpg, rhog, dP_rhog, dT_rhog, delta_h;
  F77NAME(FTVAP)(&un, &ienc, &ier, &itest, &P, &T, &Ts, &dP_Ts, &hgs, &dP_hgs, &vapa, &vapb, &vapc, &vapdb, &vapdc,
                 &hg, &dP_hg, &cpg, &dP_cpg, &dT_cpg, &rhog, &dP_rhog, &dT_rhog, &delta_h);
  return cpg;
}

double Fluide_eau_c3_gaz::dP_h_(const double T, const double P) const
{
  /* calcul a saturation */
  int un = 1, ienc = 0, ier, itest;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  /* calcul a T */
  double vapa, vapb, vapc, vapdb, vapdc, hg, dP_hg, cpg, dP_cpg, dT_cpg, rhog, dP_rhog, dT_rhog, delta_h;
  F77NAME(FTVAP)(&un, &ienc, &ier, &itest, &P, &T, &Ts, &dP_Ts, &hgs, &dP_hgs, &vapa, &vapb, &vapc, &vapdb, &vapdc,
                 &hg, &dP_hg, &cpg, &dP_cpg, &dT_cpg, &rhog, &dP_rhog, &dT_rhog, &delta_h);
  return dP_hg;
}

double Fluide_eau_c3_gaz::cp_(const double T, const double P) const
{
  /* calcul a saturation */
  int un = 1, ienc = 0, ier, itest;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  /* calcul a T */
  double vapa, vapb, vapc, vapdb, vapdc, hg, dP_hg, cpg, dP_cpg, dT_cpg, rhog, dP_rhog, dT_rhog, delta_h;
  F77NAME(FTVAP)(&un, &ienc, &ier, &itest, &P, &T, &Ts, &dP_Ts, &hgs, &dP_hgs, &vapa, &vapb, &vapc, &vapdb, &vapdc,
                 &hg, &dP_hg, &cpg, &dP_cpg, &dT_cpg, &rhog, &dP_rhog, &dT_rhog, &delta_h);
  return cpg;
}

double Fluide_eau_c3_gaz::mu_(const double T, const double P) const
{
  /* calcul a saturation */
  int un = 1, ienc = 0;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  /* calcul a T */
  double cond, dT_cond, dP_cond, visc, dP_visc, dT_visc, sigma, dP_sigma;
  F77NAME(FHVAPA)(&un, &ienc, &P, &T, &Ts, &dP_Ts, &cond, &dP_cond, &dT_cond, &visc, &dP_visc, &dT_visc, &sigma, &dP_sigma);
  return visc;
}

double Fluide_eau_c3_gaz::lambda_(const double T, const double P) const
{
  /* calcul a saturation */
  int un = 1, ienc = 0;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  /* calcul a T */
  double cond, dT_cond, dP_cond, visc, dP_visc, dT_visc, sigma, dP_sigma;
  F77NAME(FHVAPA)(&un, &ienc, &P, &T, &Ts, &dP_Ts, &cond, &dP_cond, &dT_cond, &visc, &dP_visc, &dT_visc, &sigma, &dP_sigma);
  return cond;
}

double Fluide_eau_c3_gaz::beta_(const double T, const double P) const
{
  return dT_rho_(T, P) / rho_(T, P);
}
