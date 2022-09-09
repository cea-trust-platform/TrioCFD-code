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

#include <Fluide_eau_c3_gaz.h>
#include <Lois_eau_c3.h>

Implemente_instanciable(Fluide_eau_c3_gaz, "Fluide_eau_c3_gaz", Fluide_reel_base);

Sortie& Fluide_eau_c3_gaz::printOn(Sortie& os) const { return os; }

Entree& Fluide_eau_c3_gaz::readOn(Entree& is)
{
#if HAVE_LIBC3
  Fluide_reel_base::readOn(is);
#else
  Process::exit(que_suis_je() + " : this binary was not compiled with C3 water laws!");
#endif
  return is;
}

#define ind std::distance(res.begin(), &val)

void Fluide_eau_c3_gaz::rho_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1, ienc = 0, ier, itest;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[ind], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      /* calcul a T */
      double vapa, vapb, vapc, vapdb, vapdc, hg, dP_hg, cpg, dP_cpg, dT_cpg, rhog, dP_rhog, dT_rhog, delta_h;
      F77NAME(FTVAP)(&un, &ienc, &ier, &itest, &P[ind], &T[ind * ncomp + id], &Ts, &dP_Ts, &hgs, &dP_hgs, &vapa, &vapb, &vapc, &vapdb, &vapdc,
                     &hg, &dP_hg, &cpg, &dP_cpg, &dT_cpg, &rhog, &dP_rhog, &dT_rhog, &delta_h);
      val = rhog;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_eau_c3_gaz::dP_rho_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1, ienc = 0, ier, itest;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[ind], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      /* calcul a T */
      double vapa, vapb, vapc, vapdb, vapdc, hg, dP_hg, cpg, dP_cpg, dT_cpg, rhog, dP_rhog, dT_rhog, delta_h;
      F77NAME(FTVAP)(&un, &ienc, &ier, &itest, &P[ind], &T[ind * ncomp + id], &Ts, &dP_Ts, &hgs, &dP_hgs, &vapa, &vapb, &vapc, &vapdb, &vapdc,
                     &hg, &dP_hg, &cpg, &dP_cpg, &dT_cpg, &rhog, &dP_rhog, &dT_rhog, &delta_h);
      val = dP_rhog;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_eau_c3_gaz::dT_rho_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1, ienc = 0, ier, itest;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[ind], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      /* calcul a T */
      double vapa, vapb, vapc, vapdb, vapdc, hg, dP_hg, cpg, dP_cpg, dT_cpg, rhog, dP_rhog, dT_rhog, delta_h;
      F77NAME(FTVAP)(&un, &ienc, &ier, &itest, &P[ind], &T[ind * ncomp + id], &Ts, &dP_Ts, &hgs, &dP_hgs, &vapa, &vapb, &vapc, &vapdb, &vapdc,
                     &hg, &dP_hg, &cpg, &dP_cpg, &dT_cpg, &rhog, &dP_rhog, &dT_rhog, &delta_h);
      val = dT_rhog;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_eau_c3_gaz::h_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1, ienc = 0, ier, itest;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[ind], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      /* calcul a T */
      double vapa, vapb, vapc, vapdb, vapdc, hg, dP_hg, cpg, dP_cpg, dT_cpg, rhog, dP_rhog, dT_rhog, delta_h;
      F77NAME(FTVAP)(&un, &ienc, &ier, &itest, &P[ind], &T[ind * ncomp + id], &Ts, &dP_Ts, &hgs, &dP_hgs, &vapa, &vapb, &vapc, &vapdb, &vapdc,
                     &hg, &dP_hg, &cpg, &dP_cpg, &dT_cpg, &rhog, &dP_rhog, &dT_rhog, &delta_h);
      val = hg;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_eau_c3_gaz::dP_h_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1, ienc = 0, ier, itest;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[ind], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      /* calcul a T */
      double vapa, vapb, vapc, vapdb, vapdc, hg, dP_hg, cpg, dP_cpg, dT_cpg, rhog, dP_rhog, dT_rhog, delta_h;
      F77NAME(FTVAP)(&un, &ienc, &ier, &itest, &P[ind], &T[ind * ncomp + id], &Ts, &dP_Ts, &hgs, &dP_hgs, &vapa, &vapb, &vapc, &vapdb, &vapdc,
                     &hg, &dP_hg, &cpg, &dP_cpg, &dT_cpg, &rhog, &dP_rhog, &dT_rhog, &delta_h);
      val = dP_hg;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_eau_c3_gaz::dT_h_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1, ienc = 0, ier, itest;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[ind], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      /* calcul a T */
      double vapa, vapb, vapc, vapdb, vapdc, hg, dP_hg, cpg, dP_cpg, dT_cpg, rhog, dP_rhog, dT_rhog, delta_h;
      F77NAME(FTVAP)(&un, &ienc, &ier, &itest, &P[ind], &T[ind * ncomp + id], &Ts, &dP_Ts, &hgs, &dP_hgs, &vapa, &vapb, &vapc, &vapdb, &vapdc,
                     &hg, &dP_hg, &cpg, &dP_cpg, &dT_cpg, &rhog, &dP_rhog, &dT_rhog, &delta_h);
      val = cpg;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_eau_c3_gaz::cp_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1, ienc = 0, ier, itest;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[ind], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      /* calcul a T */
      double vapa, vapb, vapc, vapdb, vapdc, hg, dP_hg, cpg, dP_cpg, dT_cpg, rhog, dP_rhog, dT_rhog, delta_h;
      F77NAME(FTVAP)(&un, &ienc, &ier, &itest, &P[ind], &T[ind * ncomp + id], &Ts, &dP_Ts, &hgs, &dP_hgs, &vapa, &vapb, &vapc, &vapdb, &vapdc,
                     &hg, &dP_hg, &cpg, &dP_cpg, &dT_cpg, &rhog, &dP_rhog, &dT_rhog, &delta_h);
      val = cpg;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_eau_c3_gaz::beta_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
  VectorD dT_rho___((int )res.size()), rho___((int )res.size());
  dT_rho_(T,P,SpanD(dT_rho___),ncomp,id);
  rho_(T,P,SpanD(rho___),ncomp,id);
  for (auto& val : res) val = dT_rho___[ind] / rho___[ind];
}

void Fluide_eau_c3_gaz::mu_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1, ienc = 0;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[ind], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      /* calcul a T */
      double cond, dT_cond, dP_cond, visc, dP_visc, dT_visc, sigma, dP_sigma;
      F77NAME(FHVAPA)(&un, &ienc, &P[ind], &T[ind * ncomp + id], &Ts, &dP_Ts, &cond, &dP_cond, &dT_cond, &visc, &dP_visc, &dT_visc, &sigma, &dP_sigma);
      val = visc;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_eau_c3_gaz::lambda_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1, ienc = 0;
      double Ts, dP_Ts, d2P_Ts, hls, dP_hls, hgs, dP_hgs, cpls, dP_cpls, cpgs, dP_cpgs, rhols, dP_rhols, rhogs, dP_rhogs;
      F77NAME(FTSATP)(&un, &ienc, &P[ind], &Ts, &dP_Ts, &d2P_Ts, &hls, &dP_hls, &hgs, &dP_hgs, &cpls, &dP_cpls, &cpgs, &dP_cpgs, &rhols, &dP_rhols, &rhogs, &dP_rhogs);
      /* calcul a T */
      double cond, dT_cond, dP_cond, visc, dP_visc, dT_visc, sigma, dP_sigma;
      F77NAME(FHVAPA)(&un, &ienc, &P[ind], &T[ind * ncomp + id], &Ts, &dP_Ts, &cond, &dP_cond, &dT_cond, &visc, &dP_visc, &dT_visc, &sigma, &dP_sigma);
      val = cond;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

#undef ind
