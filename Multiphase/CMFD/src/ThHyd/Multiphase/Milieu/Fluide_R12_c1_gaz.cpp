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

#include <Fluide_R12_c1_gaz.h>
#include <Lois_R12_c1.h>

Implemente_instanciable(Fluide_R12_c1_gaz, "Fluide_R12_c1_gaz", Fluide_reel_base);

Sortie& Fluide_R12_c1_gaz::printOn(Sortie& os) const { return os; }

Entree& Fluide_R12_c1_gaz::readOn(Entree& is)
{
#if HAVE_LIBC3
  Param param(que_suis_je());
  param.ajouter("like_eos",&like_eos_);
  param.lire_avec_accolades_depuis(is);
#else
  Process::exit(que_suis_je() + " : this binary was not compiled with C3 water laws!");
#endif
  return is;
}

#define ind std::distance(res.begin(), &val)

void Fluide_R12_c1_gaz::rho_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  if (like_eos_)
    for (auto& val : res)
      {
        // On fait comme dans EOS : des trucs chelous
        // 1. Enthalpie gaz
        int un = 1;
        double Hg;
        int ill, ivstat, ierrth;
        F77NAME(FPTHGR12)(&un, &P[ind], &T[ind * ncomp + id], &Hg, &ill, &ivstat, &ierrth);

        // 2. Saturation
        double tsat, dtsp1 ;
        F77NAME(FPSATR12)(&un, &P[ind], &tsat, &dtsp1, &ill, &ivstat, &ierrth);

        // 3. Masse volumique a saturation
        double DP_Tl = 0.; // Renvoye a 0. dans FCPLR12, on gange un appel
        double rhog_sat, rhol_sat, DP_rhoL_sat, DP_rhoG_sat;
        F77NAME(FROVLR12)(&un, &P[ind], &tsat, &tsat , &rhog_sat, &rhol_sat , &DP_Tl, &DP_rhoL_sat, &DP_rhoG_sat, &ill, &ivstat, &ierrth);

        // 4. Enthalpie a saturation
        double hvsp, hlsp, dhvsp1, dhlsp1;
        F77NAME(FHSATR12)(&un, &tsat, &rhog_sat, &rhol_sat, &P[ind], &hvsp, &hlsp, &dtsp1, &DP_rhoG_sat, &DP_rhoL_sat, &dhvsp1, &dhlsp1);

        // 5. Inversion des inconnues, sort la masse volumique...
        double Tg, dP_Tg, dHg_Tg, rhog, dP_rhog, dHg_rhog;
        F77NAME(FTGR12)(&un, &P[ind], &Hg, &tsat, &rhog_sat, &Tg , &dP_Tg, &dHg_Tg, &rhog, &dP_rhog, &dHg_rhog , &ill, &ivstat, &ierrth);

        val = rhog; // NB : only a function of T
      }
  else   for (auto& val : res)
      {
        int un = 1;
        int ill, ivstat, ierrth;

        // 1. Saturation
        double tsat, dtsp1 ;
        F77NAME(FPSATR12)(&un, &P[ind], &tsat, &dtsp1, &ill, &ivstat, &ierrth);

        // 2. finex
        double DP_Tl = 0.; // Renvoye a 0. dans FCPLR12, on gange un appel
        double rhog, rhol, DP_rhoL, DP_rhoG;
        F77NAME(FROVLR12)(&un, &P[ind], &T[ind * ncomp + id], &tsat ,  &rhog , &rhol, &DP_Tl, &DP_rhoL, &DP_rhoG, &ill, &ivstat, &ierrth);

        val = rhog; // NB : only a function of T
      }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_R12_c1_gaz::dP_rho_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      // On fait comme dans EOS : des trucs chelous
      // 1. Enthalpie gaz
      int un = 1;
      double Hg;
      int ill, ivstat, ierrth;
      F77NAME(FPTHGR12)(&un, &P[ind], &T[ind * ncomp + id], &Hg, &ill, &ivstat, &ierrth);

      // 2. Saturation
      double tsat, dtsp1 ;
      F77NAME(FPSATR12)(&un, &P[ind], &tsat, &dtsp1, &ill, &ivstat, &ierrth);

      // 3. Masse volumique a saturation
      double DP_Tl = 0.; // Renvoye a 0. dans FCPLR12, on gange un appel
      double rhog_sat, rhol_sat, DP_rhoL_sat, DP_rhoG_sat;
      F77NAME(FROVLR12)(&un, &P[ind], &tsat, &tsat , &rhog_sat, &rhol_sat , &DP_Tl, &DP_rhoL_sat, &DP_rhoG_sat, &ill, &ivstat, &ierrth);

      // 4. Enthalpie a saturation
      double hvsp, hlsp, dhvsp1, dhlsp1;
      F77NAME(FHSATR12)(&un, &tsat, &rhog_sat, &rhol_sat, &P[ind], &hvsp, &hlsp, &dtsp1, &DP_rhoG_sat, &DP_rhoL_sat, &dhvsp1, &dhlsp1);

      // 5. Inversion des inconnues, sort la masse volumique...
      double Tg, dP_Tg, dHg_Tg, rhog, dP_rhog, dHg_rhog;
      F77NAME(FTGR12)(&un, &P[ind], &Hg, &tsat, &rhog_sat, &Tg , &dP_Tg, &dHg_Tg, &rhog, &dP_rhog, &dHg_rhog ,  &ill,  &ivstat, &ierrth);

      val = dP_rhog - dHg_rhog* dP_Tg/dHg_Tg; // NB : this is derivative with T fixed, not with H fixed
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_R12_c1_gaz::dT_rho_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      // On fait comme dans EOS : des trucs chelous
      // 1. Enthalpie gaz
      int un = 1;
      double Hg;
      int ill, ivstat, ierrth;
      F77NAME(FPTHGR12)(&un, &P[ind], &T[ind * ncomp + id], &Hg, &ill, &ivstat, &ierrth);

      // 2. Saturation
      double tsat, dtsp1 ;
      F77NAME(FPSATR12)(&un, &P[ind], &tsat, &dtsp1, &ill, &ivstat, &ierrth);

      // 3. Masse volumique a saturation
      double DP_Tl = 0.; // Renvoye a 0. dans FCPLR12, on gange un appel
      double rhog_sat, rhol_sat, DP_rhoL_sat, DP_rhoG_sat;
      F77NAME(FROVLR12)(&un, &P[ind], &tsat, &tsat , &rhog_sat, &rhol_sat , &DP_Tl, &DP_rhoL_sat, &DP_rhoG_sat, &ill, &ivstat, &ierrth);

      // 4. Enthalpie a saturation
      double hvsp, hlsp, dhvsp1, dhlsp1;
      F77NAME(FHSATR12)(&un, &tsat, &rhog_sat, &rhol_sat, &P[ind], &hvsp, &hlsp, &dtsp1, &DP_rhoG_sat, &DP_rhoL_sat, &dhvsp1, &dhlsp1);

      // 5. Inversion des inconnues, sort la masse volumique...
      double Tg, dP_Tg, dHg_Tg, rhog, dP_rhog, dHg_rhog;
      F77NAME(FTGR12)(&un, &P[ind], &Hg, &tsat, &rhog_sat, &Tg , &dP_Tg, &dHg_Tg, &rhog, &dP_rhog, &dHg_rhog , &ill,  &ivstat, &ierrth);

      val = dHg_rhog/dHg_Tg;  // Dans la formule de FPTHGR, H ne depend que de rho et de T.
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_R12_c1_gaz::h_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1;
      double Hg;
      int ill, ivstat, ierrth;
      F77NAME(FPTHGR12)(&un, &P[ind], &T[ind * ncomp + id], &Hg, &ill, &ivstat, &ierrth);
      val = Hg;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_R12_c1_gaz::dP_h_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      // On fait comme dans EOS : des trucs chelous
      // 1. Enthalpie gaz
      int un = 1;
      double Hg;
      int ill, ivstat, ierrth;
      F77NAME(FPTHGR12)(&un, &P[ind], &T[ind * ncomp + id], &Hg, &ill, &ivstat, &ierrth);

      // 2. Saturation
      double tsat, dtsp1 ;
      F77NAME(FPSATR12)(&un, &P[ind], &tsat, &dtsp1, &ill, &ivstat, &ierrth);

      // 3. Masse volumique a saturation
      double DP_Tl = 0.; // Renvoye a 0. dans FCPLR12, on gange un appel
      double rhog_sat, rhol_sat, DP_rhoL_sat, DP_rhoG_sat;
      F77NAME(FROVLR12)(&un, &P[ind], &tsat, &tsat , &rhog_sat, &rhol_sat , &DP_Tl, &DP_rhoL_sat, &DP_rhoG_sat, &ill, &ivstat, &ierrth);

      // 4. Enthalpie a saturation
      double hvsp, hlsp, dhvsp1, dhlsp1;
      F77NAME(FHSATR12)(&un, &tsat, &rhog_sat, &rhol_sat, &P[ind], &hvsp, &hlsp, &dtsp1, &DP_rhoG_sat, &DP_rhoL_sat, &dhvsp1, &dhlsp1);

      // 5. Inversion des inconnues, sort la masse volumique...
      double Tg, dP_Tg, dHg_Tg, rhog, dP_rhog, dHg_rhog;
      F77NAME(FTGR12)(&un, &P[ind], &Hg, &tsat, &rhog_sat, &Tg , &dP_Tg, &dHg_Tg, &rhog, &dP_rhog, &dHg_rhog , &ill,  &ivstat, &ierrth);

      val = -dP_Tg/dHg_Tg;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_R12_c1_gaz::dT_h_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      // On fait comme dans EOS : des trucs chelous
      // 1. Enthalpie gaz
      int un = 1;
      double Hg;
      int ill, ivstat, ierrth;
      F77NAME(FPTHGR12)(&un, &P[ind], &T[ind * ncomp + id], &Hg, &ill, &ivstat, &ierrth);

      // 2. Saturation
      double tsat, dtsp1 ;
      F77NAME(FPSATR12)(&un, &P[ind], &tsat, &dtsp1, &ill, &ivstat, &ierrth);

      // 3. Masse volumique a saturation
      double DP_Tl = 0.; // Renvoye a 0. dans FCPLR12, on gange un appel
      double rhog_sat, rhol_sat, DP_rhoL_sat, DP_rhoG_sat;
      F77NAME(FROVLR12)(&un, &P[ind], &tsat, &tsat , &rhog_sat, &rhol_sat , &DP_Tl, &DP_rhoL_sat, &DP_rhoG_sat, &ill, &ivstat, &ierrth);

      // 4. Enthalpie a saturation
      double hvsp, hlsp, dhvsp1, dhlsp1;
      F77NAME(FHSATR12)(&un, &tsat, &rhog_sat, &rhol_sat, &P[ind], &hvsp, &hlsp, &dtsp1, &DP_rhoG_sat, &DP_rhoL_sat, &dhvsp1, &dhlsp1);

      // 5. Inversion des inconnues, sort la masse volumique...
      double Tg, dP_Tg, dHg_Tg, rhog, dP_rhog, dHg_rhog;
      F77NAME(FTGR12)(&un, &P[ind], &Hg, &tsat, &rhog_sat, &Tg , &dP_Tg, &dHg_Tg, &rhog, &dP_rhog, &dHg_rhog , &ill,  &ivstat, &ierrth);

      val = 1./dHg_Tg;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_R12_c1_gaz::cp_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      // On fait comme dans EOS : des trucs chelous
      // 1. Enthalpie gaz
      int un = 1;
      double Hg;
      int ill, ivstat, ierrth;
      F77NAME(FPTHGR12)(&un, &P[ind], &T[ind * ncomp + id], &Hg, &ill, &ivstat, &ierrth);

      // 2. Saturation
      double tsat, dtsp1 ;
      F77NAME(FPSATR12)(&un, &P[ind], &tsat, &dtsp1, &ill, &ivstat, &ierrth);

      // 3. Masse volumique a saturation
      double DP_Tl = 0.; // Renvoye a 0. dans FCPLR12, on gange un appel
      double rhog_sat, rhol_sat, DP_rhoL_sat, DP_rhoG_sat;
      F77NAME(FROVLR12)(&un, &P[ind], &tsat, &tsat , &rhog_sat, &rhol_sat , &DP_Tl, &DP_rhoL_sat, &DP_rhoG_sat, &ill, &ivstat, &ierrth);

      // 4. Enthalpie a saturation
      double hvsp, hlsp, dhvsp1, dhlsp1;
      F77NAME(FHSATR12)(&un, &tsat, &rhog_sat, &rhol_sat, &P[ind], &hvsp, &hlsp, &dtsp1, &DP_rhoG_sat, &DP_rhoL_sat, &dhvsp1, &dhlsp1);

      // 5. Inversion des inconnues, sort la masse volumique...
      double Tg, dP_Tg, dHg_Tg, rhog, dP_rhog, dHg_rhog;
      F77NAME(FTGR12)(&un, &P[ind], &Hg, &tsat, &rhog_sat, &Tg , &dP_Tg, &dHg_Tg, &rhog, &dP_rhog, &dHg_rhog , &ill,  &ivstat, &ierrth);

      // 5. On appelle la fonction qu'on cherche
      double cpg, dP_cpg, dhg_cpg;
      F77NAME(FCPVR12)(&un, &Tg, &dP_Tg, &dHg_Tg, &rhog, &dP_rhog, &dHg_rhog, &cpg, &dP_cpg, &dhg_cpg);

      val = cpg;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_R12_c1_gaz::beta_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
  VectorD dT_rho___((int )res.size()), rho___((int )res.size());
  dT_rho_(T,P,SpanD(dT_rho___),ncomp,id);
  rho_(T,P,SpanD(rho___),ncomp,id);
  for (auto& val : res) val = dT_rho___[ind] / rho___[ind];
}

void Fluide_R12_c1_gaz::mu_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      // 1. Enthalpie gaz
      int un = 1;
      double Hg;
      int ill, ivstat, ierrth;
      F77NAME(FPTHGR12)(&un, &P[ind], &T[ind * ncomp + id], &Hg, &ill, &ivstat, &ierrth);

      // 2. Saturation
      double tsat, dtsp1 ;
      F77NAME(FPSATR12)(&un, &P[ind], &tsat, &dtsp1, &ill, &ivstat, &ierrth);

      // 3. Masse volumique du gaz a saturation
      double DP_Tl = 0.; // Renvoye a 0. dans FCPLR12, on gange un appel
      double rhog_sat, rhol_sat, DP_rhoL_sat, DP_rhoG_sat;
      F77NAME(FROVLR12)(&un, &P[ind], &tsat, &tsat , &rhog_sat, &rhol_sat , &DP_Tl, &DP_rhoL_sat, &DP_rhoG_sat, &ill, &ivstat, &ierrth);

      // 4. On appelle la fonction chelou
      double Tg, dP_Tg, dHg_Tg, rhog, dP_rhog, dHg_rhog;
      F77NAME(FTGR12)(&un, &P[ind], &Hg, &tsat, &rhog_sat, &Tg, &dP_Tg, &dHg_Tg, &rhog, &dP_rhog, &dHg_rhog, &ill, &ivstat, &ierrth);

      // 5. On appelle la fonction qu'on cherche
      double mug, dP_mug, DHg_mug;
      F77NAME(FMUVR12)(&un, &P[ind], &T[ind * ncomp + id], &dP_Tg, &dHg_Tg, &mug, &dP_mug, &DHg_mug);

      val = mug;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_R12_c1_gaz::lambda_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      // 1. Enthalpie gaz
      int un = 1;
      double Hg;
      int ill, ivstat, ierrth;
      F77NAME(FPTHGR12)(&un, &P[ind], &T[ind * ncomp + id], &Hg, &ill, &ivstat, &ierrth);

      // 2. Saturation
      double tsat, dtsp1 ;
      F77NAME(FPSATR12)(&un, &P[ind], &tsat, &dtsp1, &ill, &ivstat, &ierrth);

      // 3. Masse volumique du gaz a saturation
      double DP_Tl = 0.; // Renvoye a 0. dans FCPLR12, on gange un appel
      double rhog_sat, rhol_sat, DP_rhoL_sat, DP_rhoG_sat;
      F77NAME(FROVLR12)(&un, &P[ind], &tsat, &tsat , &rhog_sat, &rhol_sat , &DP_Tl, &DP_rhoL_sat, &DP_rhoG_sat, &ill, &ivstat, &ierrth);

      // 4. On appelle la fonction chelou
      double Tg, dP_Tg, dHg_Tg, rhog, dP_rhog, dHg_rhog;
      F77NAME(FTGR12)(&un, &P[ind], &Hg, &tsat, &rhog_sat, &Tg, &dP_Tg, &dHg_Tg, &rhog, &dP_rhog, &dHg_rhog, &ill, &ivstat, &ierrth);

      // 5. On appelle la fonction qu'on cherche
      double lambdag, dP_lambdag, DHg_lambdag;
      F77NAME(FMUVR12)(&un, &P[ind], &T[ind * ncomp + id], &dP_Tg, &dHg_Tg, &lambdag, &dP_lambdag, &DHg_lambdag);

      val = lambdag;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

#undef ind
