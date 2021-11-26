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
// File:        Saturation_eau_c3.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Milieu
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Saturation_eau_c3.h>
#include <Lois_eau_c3.h>

Implemente_instanciable(Saturation_eau_c3, "Saturation_eau_c3", Saturation_base);

Sortie& Saturation_eau_c3::printOn(Sortie& os) const
{
  return os;
}

Entree& Saturation_eau_c3::readOn(Entree& is)
{
  return Saturation_base::readOn(is);
}

double Saturation_eau_c3::Tsat_(const double P) const
{
  int un = 1, ienc = 0;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  return Ts;
}
double Saturation_eau_c3::dP_Tsat_(const double P) const
{
  int un = 1, ienc = 0;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  return dP_Ts;
}
double Saturation_eau_c3::Psat_(const double T) const
{
  int ienc = 0;
  double Ps, dT_Ps, d2T_Ps, hls, dT_hls, hgs, dT_hgs, cpls, dT_cpls, cpgs, dT_cpgs, rhols, dT_rhols, rhogs, dT_rhogs;
  F77NAME(FPSATT)(&ienc, &T, &Ps, &dT_Ps, &d2T_Ps, &hls, &dT_hls, &hgs, &dT_hgs, &cpls, &dT_cpls, &cpgs, &dT_cpgs, &rhols, &dT_rhols, &rhogs, &dT_rhogs);
  return Ps;
}
double Saturation_eau_c3::dT_Psat_(const double T) const
{
  int ienc = 0;
  double Ps, dT_Ps, d2T_Ps, hls, dT_hls, hgs, dT_hgs, cpls, dT_cpls, cpgs, dT_cpgs, rhols, dT_rhols, rhogs, dT_rhogs;
  F77NAME(FPSATT)(&ienc, &T, &Ps, &dT_Ps, &d2T_Ps, &hls, &dT_hls, &hgs, &dT_hgs, &cpls, &dT_cpls, &cpgs, &dT_cpgs, &rhols, &dT_rhols, &rhogs, &dT_rhogs);
  return dT_Ps;
}
double Saturation_eau_c3::Lvap_(const double P) const
{
  int un = 1, ienc = 0;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  return hgs - hls;
}
double Saturation_eau_c3::dP_Lvap_(const double P) const
{
  int un = 1, ienc = 0;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  return dP_hgs - dP_hls;
}
double Saturation_eau_c3::Hls_(const double P) const
{
  int un = 1, ienc = 0;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  return hls;
}
double Saturation_eau_c3::dP_Hls_(const double P) const
{
  int un = 1, ienc = 0;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  return dP_hls;
}
double Saturation_eau_c3::Hvs_(const double P) const
{
  int un = 1, ienc = 0;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  return hgs;
}
double Saturation_eau_c3::dP_Hvs_(const double P) const
{
  int un = 1, ienc = 0;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  return dP_hgs;
}

double Saturation_eau_c3::sigma_(const double T, const double P) const
{
  /* calcul a saturation */
  int un = 1, ienc = 0;
  double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
  F77NAME(FTSATP)(&un, &ienc, &P, &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
  /* calcul a T */
  double cond, dT_cond, dP_cond, visc, dP_visc, dT_visc, sigma, dP_sigma;
  F77NAME(FHVAPA)(&un, &ienc, &P, &T, &Ts, &dP_Ts, &cond, &dP_cond, &dT_cond, &visc, &dP_visc, &dT_visc, &sigma, &dP_sigma);
  return sigma;
}
