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
// File:        Saturation_R12_c1.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Milieu
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Saturation_R12_c1.h>
#include <Lois_R12_c1.h>

Implemente_instanciable(Saturation_R12_c1, "Saturation_R12_c1", Saturation_base);

Sortie& Saturation_R12_c1::printOn(Sortie& os) const
{
  return os;
}

Entree& Saturation_R12_c1::readOn(Entree& is)
{
#if HAVE_LIBC3
  Saturation_base::readOn(is);
#else
  Process::exit(que_suis_je() + " : this binary was not compiled with C3 water laws!");
#endif
  return is;
}

void Saturation_R12_c1::Tsat_(const SpanD P, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      int un=1, ill, ivstat, ierrth ;
      double tsp, dtsp1 ;
      F77NAME(FPSATR12)(&un, &P[i], &tsp, &dtsp1, &ill, &ivstat, &ierrth);
      res[i * ncomp + ind] = tsp;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}

void Saturation_R12_c1::dP_Tsat_(const SpanD P, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      int un=1, ill, ivstat, ierrth ;
      double tsp, dtsp1 ;
      F77NAME(FPSATR12)(&un, &P[i], &tsp, &dtsp1, &ill, &ivstat, &ierrth);
      res[i * ncomp + ind] = dtsp1;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}

void Saturation_R12_c1::Lvap_(const SpanD P, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      // Calcule Tsat
      int un=1, ill, ivstat, ierrth ;
      double tsp, dtsp1 ;
      F77NAME(FPSATR12)(&un, &P[i], &tsp, &dtsp1, &ill, &ivstat, &ierrth);

      // Calcule masses volumiques (avec tl, tg = tsat)
      double rg, rl, drl1, drg1 ;
      double DP_Tl = 0.; // Renvoye a 0. dans FCPLR12, on gange un appel
      F77NAME(FROVLR12)(&un, &P[i], &tsp, &tsp, &rg, &rl, &DP_Tl, &drl1, &drg1,  &ill, &ivstat, &ierrth);

      // Calcule les enthalpies
      double hvsp, hlsp, dhvsp1, dhlsp1  ;
      F77NAME(FHSATR12)(&un, &tsp, &rg, &rl, &P[i], &hvsp, &hlsp, &dtsp1, &drg1, &drl1, &dhvsp1, &dhlsp1);

      // Calcule (enfin...) Lvap
      res[i * ncomp + ind] =  hvsp - hlsp;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}

void Saturation_R12_c1::dP_Lvap_(const SpanD P, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      // Calcule Tsat
      int un=1, ill, ivstat, ierrth ;
      double tsp, dtsp1 ;
      F77NAME(FPSATR12)(&un, &P[i], &tsp, &dtsp1, &ill, &ivstat, &ierrth);

      // Calcule masses volumiques (avec tl, tg = tsat)
      double rg, rl, drl1, drg1 ;
      double DP_Tl = 0.; // Renvoye a 0. dans FCPLR12, on gange un appel
      F77NAME(FROVLR12)(&un, &P[i], &tsp, &tsp, &rg, &rl, &DP_Tl, &drl1, &drg1,    &ill, &ivstat, &ierrth);

      // Calcule les derivees des enthalpies
      double hvsp, hlsp, dhvsp1, dhlsp1 ;
      F77NAME(FHSATR12)(&un, &tsp, &rg, &rl, &P[i], &hvsp, &hlsp, &dtsp1, &drg1, &drl1, &dhvsp1, &dhlsp1);

      // Calcule (enfin...) dP_Lvap
      res[i * ncomp + ind] = dhvsp1 - dhlsp1;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}

void Saturation_R12_c1::Hls_(const SpanD P, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      // Calcule Tsat
      int un=1, ill, ivstat, ierrth ;
      double tsp, dtsp1 ;
      F77NAME(FPSATR12)(&un, &P[i], &tsp, &dtsp1, &ill, &ivstat, &ierrth);

      // Calcule masses volumiques (avec tl, tg = tsat)
      double rg, rl, drl1, drg1 ;
      double DP_Tl = 0.; // Renvoye a 0. dans FCPLR12, on gange un appel
      F77NAME(FROVLR12)(&un, &P[i], &tsp, &tsp, &rg, &rl, &DP_Tl, &drl1, &drg1,    &ill, &ivstat, &ierrth);

      // Calcule les enthalpies
      double hvsp, hlsp, dhvsp1, dhlsp1 ;
      F77NAME(FHSATR12)(&un, &tsp, &rg, &rl, &P[i], &hvsp, &hlsp, &dtsp1, &drg1, &drl1, &dhvsp1, &dhlsp1);
      res[i * ncomp + ind] = hlsp;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}

void Saturation_R12_c1::dP_Hls_(const SpanD P, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      // Calcule Tsat
      int un=1, ill, ivstat, ierrth ;
      double tsp, dtsp1 ;
      F77NAME(FPSATR12)(&un, &P[i], &tsp, &dtsp1, &ill, &ivstat, &ierrth);

      // Calcule masses volumiques (avec tl, tg = tsat)
      double rg, rl, drl1, drg1 ;
      double DP_Tl = 0.; // Renvoye a 0. dans FCPLR12, on gange un appel
      F77NAME(FROVLR12)(&un, &P[i], &tsp, &tsp, &rg, &rl, &DP_Tl, &drl1, &drg1,    &ill, &ivstat, &ierrth);

      // Calcule les enthalpies
      double hvsp, hlsp, dhvsp1, dhlsp1;
      F77NAME(FHSATR12)(&un, &tsp, &rg, &rl, &P[i], &hvsp, &hlsp, &dtsp1, &drg1, &drl1, &dhvsp1, &dhlsp1);
      res[i * ncomp + ind] = dhlsp1;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}

void Saturation_R12_c1::Hvs_(const SpanD P, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      // Calcule Tsat
      int un=1, ill, ivstat, ierrth ;
      double tsp, dtsp1 ;
      F77NAME(FPSATR12)(&un, &P[i], &tsp, &dtsp1, &ill, &ivstat, &ierrth);

      // Calcule masses volumiques (avec tl, tg = tsat)
      double rg, rl, drl1, drg1 ;
      double DP_Tl = 0.; // Renvoye a 0. dans FCPLR12, on gange un appel
      F77NAME(FROVLR12)(&un, &P[i], &tsp, &tsp, &rg, &rl, &DP_Tl, &drl1, &drg1,    &ill, &ivstat, &ierrth);

      // Calcule les enthalpies
      double hvsp, hlsp, dhvsp1, dhlsp1 ;
      F77NAME(FHSATR12)(&un, &tsp, &rg, &rl, &P[i], &hvsp, &hlsp, &dtsp1, &drg1, &drl1, &dhvsp1, &dhlsp1);
      res[i * ncomp + ind] = hvsp;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}

void Saturation_R12_c1::dP_Hvs_(const SpanD P, SpanD res, int ncomp, int ind) const
{
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      // Calcule Tsat
      int un=1, ill, ivstat, ierrth ;
      double tsp, dtsp1 ;
      F77NAME(FPSATR12)(&un, &P[i], &tsp, &dtsp1, &ill, &ivstat, &ierrth);

      // Calcule masses volumiques (avec tl, tg = tsat)
      double rg, rl, drl1, drg1 ;
      double DP_Tl = 0.; // Renvoye a 0. dans FCPLR12, on gange un appel
      F77NAME(FROVLR12)(&un, &P[i], &tsp, &tsp, &rg, &rl, &DP_Tl, &drl1, &drg1,    &ill, &ivstat, &ierrth);

      // Calcule les enthalpies
      double hvsp, hlsp, dhvsp1, dhlsp1 ;
      F77NAME(FHSATR12)(&un, &tsp, &rg, &rl, &P[i], &hvsp, &hlsp, &dtsp1, &drg1, &drl1, &dhvsp1, &dhlsp1);
      res[i * ncomp + ind] = dhvsp1;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}

void Saturation_R12_c1::sigma_(const SpanD T, const SpanD P, SpanD res, int ncomp, int ind) const
{
  //TODO : FIXME : verifier indices ...
#if HAVE_LIBC3
  for (int i =0; i < (int)P.size(); i++)
    {
      // Calcule Tsat
      int un=1, ill, ivstat, ierrth ;
      double tsp, dtsp1 ;
      F77NAME(FPSATR12)(&un, &P[i], &tsp, &dtsp1, &ill, &ivstat, &ierrth);

      // Calcule sigma
      double sigma, dsig1;
      F77NAME(FSIGMAR12)(&un, &tsp, &dtsp1, &sigma, &dsig1);
      res[i*ncomp+ind] = sigma;
    }
#else
  for (int i =0; i < (int)P.size(); i++) res[i] = 0;
#endif
}
