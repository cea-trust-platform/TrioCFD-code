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

#include <Fluide_eau_c3_liquide.h>
#include <Lois_eau_c3.h>

Implemente_instanciable(Fluide_eau_c3_liquide, "Fluide_eau_c3_liquide", Fluide_reel_base);

Sortie& Fluide_eau_c3_liquide::printOn(Sortie& os) const { return os; }

Entree& Fluide_eau_c3_liquide::readOn(Entree& is)
{
#if HAVE_LIBC3
  Fluide_reel_base::readOn(is);
#else
  Process::exit(que_suis_je() + " : this binary was not compiled with C3 water laws!");
#endif
  return is;
}

#define ind std::distance(res.begin(), &val)

void Fluide_eau_c3_liquide::rho_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1;
      double hl, dP_hl, dT_hl, cpl, dT_cpl, dP_cpl, rhol, dT_rhol, dP_rhol; //sorties
      F77NAME(FTLIQ)(&un, &P[ind], &T[ind * ncomp + id], &hl, &dP_hl, &dT_hl, &cpl, &dP_cpl, &dT_cpl, &rhol, &dP_rhol, &dT_rhol);
      return rhol;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_eau_c3_liquide::dP_rho_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1;
      double hl, dP_hl, dT_hl, cpl, dT_cpl, dP_cpl, rhol, dT_rhol, dP_rhol; //sorties
      F77NAME(FTLIQ)(&un, &P[ind], &T[ind * ncomp + id], &hl, &dP_hl, &dT_hl, &cpl, &dP_cpl, &dT_cpl, &rhol, &dP_rhol, &dT_rhol);
      return dP_rhol;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_eau_c3_liquide::dT_rho_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1;
      double hl, dP_hl, dT_hl, cpl, dT_cpl, dP_cpl, rhol, dT_rhol, dP_rhol; //sorties
      F77NAME(FTLIQ)(&un, &P[ind], &T[ind * ncomp + id], &hl, &dP_hl, &dT_hl, &cpl, &dP_cpl, &dT_cpl, &rhol, &dP_rhol, &dT_rhol);
      return dT_rhol;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_eau_c3_liquide::h_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1;
      double hl, dP_hl, dT_hl, cpl, dT_cpl, dP_cpl, rhol, dT_rhol, dP_rhol; //sorties
      F77NAME(FTLIQ)(&un, &P[ind], &T[ind * ncomp + id], &hl, &dP_hl, &dT_hl, &cpl, &dP_cpl, &dT_cpl, &rhol, &dP_rhol, &dT_rhol);
      return hl;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_eau_c3_liquide::dP_h_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1;
      double hl, dP_hl, dT_hl, cpl, dT_cpl, dP_cpl, rhol, dT_rhol, dP_rhol; //sorties
      F77NAME(FTLIQ)(&un, &P[ind], &T[ind * ncomp + id], &hl, &dP_hl, &dT_hl, &cpl, &dP_cpl, &dT_cpl, &rhol, &dP_rhol, &dT_rhol);
      return dP_hl;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_eau_c3_liquide::dT_h_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1;
      double hl, dP_hl, dT_hl, cpl, dT_cpl, dP_cpl, rhol, dT_rhol, dP_rhol; //sorties
      F77NAME(FTLIQ)(&un, &P[ind], &T[ind * ncomp + id], &hl, &dP_hl, &dT_hl, &cpl, &dP_cpl, &dT_cpl, &rhol, &dP_rhol, &dT_rhol);
      return dT_hl;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_eau_c3_liquide::cp_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1;
      double hl, dP_hl, dT_hl, cpl, dT_cpl, dP_cpl, rhol, dT_rhol, dP_rhol; //sorties
      F77NAME(FTLIQ)(&un, &P[ind], &T[ind * ncomp + id], &hl, &dP_hl, &dT_hl, &cpl, &dP_cpl, &dT_cpl, &rhol, &dP_rhol, &dT_rhol);
      return cpl;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_eau_c3_liquide::beta_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
  VectorD dT_rho___((int )res.size()), rho___((int )res.size());
  dT_rho_(T,P,SpanD(dT_rho___),ncomp,id);
  rho_(T,P,SpanD(rho___),ncomp,id);
  for (auto& val : res) val = dT_rho___[ind] / rho___[ind];
}

void Fluide_eau_c3_liquide::mu_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1;
      double hl = h_(T, P), zero = 0, cond, dcond1, dcond2, visc, dvisc1, dvisc2;
      F77NAME(FHLIQA)(&un, &P[ind], &hl, &T[ind * ncomp + id], &zero, &zero, &cond, &dcond1, &dcond2, &visc, &dvisc1, &dvisc2);
      return visc;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_eau_c3_liquide::lambda_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1;
      double hl = h_(T, P), zero = 0, cond, dcond1, dcond2, visc, dvisc1, dvisc2;
      F77NAME(FHLIQA)(&un, &P[ind], &hl, &T[ind * ncomp + id], &zero, &zero, &cond, &dcond1, &dcond2, &visc, &dvisc1, &dvisc2);
      return cond;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

#undef ind
