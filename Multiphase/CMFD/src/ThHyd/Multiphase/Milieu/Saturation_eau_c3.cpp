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
#if HAVE_LIBC3
  Saturation_base::readOn(is);
#else
  Process::exit(que_suis_je() + " : this binary was not compiled with C3 water laws!");
#endif
  return is;
}

void Saturation_eau_c3::Tsat_(const SpanD P, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      int un = 1, ienc = 0;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[i], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      res[i * ncomp + ind] = Ts;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}

void Saturation_eau_c3::dP_Tsat_(const SpanD P, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      int un = 1, ienc = 0;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[i], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      res[i * ncomp + ind] =  dP_Ts;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}

void Saturation_eau_c3::Psat_(const SpanD T, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)T.size() / ncomp; i++)
    {
      int ienc = 0;
      double Ps, dT_Ps, d2T_Ps, hls, dT_hls, hgs, dT_hgs, cpls, dT_cpls, cpgs, dT_cpgs, rhols, dT_rhols, rhogs, dT_rhogs;
      F77NAME(FPSATT)(&ienc, &T[i * ncomp + ind], &Ps, &dT_Ps, &d2T_Ps, &hls, &dT_hls, &hgs, &dT_hgs, &cpls, &dT_cpls, &cpgs, &dT_cpgs, &rhols, &dT_rhols, &rhogs, &dT_rhogs);
      res[i * ncomp + ind] = Ps;
    }
#else
  for (int i =0; i < (int)T.size() / ncomp; i++) res[i] = 0;
#endif
}

void Saturation_eau_c3::dT_Psat_(const SpanD T, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)T.size() / ncomp; i++)
    {
      int ienc = 0;
      double Ps, dT_Ps, d2T_Ps, hls, dT_hls, hgs, dT_hgs, cpls, dT_cpls, cpgs, dT_cpgs, rhols, dT_rhols, rhogs, dT_rhogs;
      F77NAME(FPSATT)(&ienc, &T[i * ncomp + ind], &Ps, &dT_Ps, &d2T_Ps, &hls, &dT_hls, &hgs, &dT_hgs, &cpls, &dT_cpls, &cpgs, &dT_cpgs, &rhols, &dT_rhols, &rhogs, &dT_rhogs);
      res[i * ncomp + ind] = dT_Ps;
    }
#else
  for (int i =0; i < (int)T.size() / ncomp; i++) res[i] = 0;
#endif
}

void Saturation_eau_c3::Lvap_(const SpanD P, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      int un = 1, ienc = 0;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[i], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      res[i * ncomp + ind] =  hgs - hls;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}

void Saturation_eau_c3::dP_Lvap_(const SpanD P, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      int un = 1, ienc = 0;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[i], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      res[i * ncomp + ind] = dP_hgs - dP_hls;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}

void Saturation_eau_c3::Hls_(const SpanD P, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      int un = 1, ienc = 0;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[i], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      res[i * ncomp + ind] = hls;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}

void Saturation_eau_c3::dP_Hls_(const SpanD P, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      int un = 1, ienc = 0;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[i], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      res[i * ncomp + ind] = dP_hls;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}

void Saturation_eau_c3::Hvs_(const SpanD P, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      int un = 1, ienc = 0;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[i], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      res[i * ncomp + ind] = hgs;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}

void Saturation_eau_c3::dP_Hvs_(const SpanD P, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      int un = 1, ienc = 0;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[i], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      res[i * ncomp + ind] = dP_hgs;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}

void Saturation_eau_c3::sigma_(const SpanD T, const SpanD P, SpanD res, int ncomp, int ind) const
{
  //TODO : FIXME : verifier indices ...
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      /* calcul a saturation */
      int un = 1, ienc = 0;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[i], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      /* calcul a T */
      double cond, dT_cond, dP_cond, visc, dP_visc, dT_visc, surfaceTension, dP_surfaceTension;
      F77NAME(FHVAPA)(&un, &ienc, &P[i], &T[i * ncomp + ind], &Ts, &dP_Ts, &cond, &dP_cond, &dT_cond, &visc, &dP_visc, &dT_visc, &surfaceTension, &dP_surfaceTension);
      res[i] = surfaceTension;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}
