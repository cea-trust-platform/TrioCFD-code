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

#include <Fluide_R12_c1_liquide.h>
#include <Lois_R12_c1.h>

Implemente_instanciable(Fluide_R12_c1_liquide, "Fluide_R12_c1_liquide", Fluide_reel_base);

Sortie& Fluide_R12_c1_liquide::printOn(Sortie& os) const { return os; }

Entree& Fluide_R12_c1_liquide::readOn(Entree& is)
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

void Fluide_R12_c1_liquide::rho_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  if (like_eos_)
    for (auto& val : res)
      {
        // On fait comme dans EOS : des trucs chelous
        // 1. Enthalpie liquide
        int un = 1;
        double Hl;
        int ill, ivstat, ierrth;
        F77NAME(FPTHLR12)(&un, &T[ind * ncomp + id], &Hl, &ill, &ivstat, &ierrth);

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
        double Tl, rhol;
        F77NAME(FTLR12)(&un, &P[ind], &tsat, &rhog_sat, &hvsp , &Hl, &Tl, &rhol, &ill, &ivstat, &ierrth);

        val = rhol; // NB : only a function of T
      }
  else
    for (auto& val : res)
      {
        int un = 1;
        int ill, ivstat, ierrth;

        // 1. Saturation
        double tsat, dtsp1 ;
        F77NAME(FPSATR12)(&un, &P[ind], &tsat, &dtsp1, &ill, &ivstat, &ierrth);

        // 2. Allez hop
        double rhog, rhol, DP_rhoL, DP_rhoG;
        double DP_Tl = 0.; // Renvoye a 0. dans FCPLR12, on gagne un appel
        F77NAME(FROVLR12)(&un, &P[ind], &tsat, &T[ind * ncomp + id] , &rhog, &rhol , &DP_Tl, &DP_rhoL, &DP_rhoG, &ill, &ivstat, &ierrth);

        val = rhol; // NB : only a function of T
      }

#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_R12_c1_liquide::dP_rho_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  if (like_eos_)
    for (auto& val : res)
      {
        // On fait comme dans EOS : des trucs chelous
        // 1. Enthalpie liquide
        int un = 1;
        double Hl;
        int ill, ivstat, ierrth;
        F77NAME(FPTHLR12)(&un, &T[ind * ncomp + id], &Hl, &ill, &ivstat, &ierrth);

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

        // 5. Inversion des inconnues...
        double Tl, rhol;
        F77NAME(FTLR12)(&un, &P[ind], &tsat, &rhog_sat, &hvsp , &Hl, &Tl, &rhol, &ill, &ivstat, &ierrth);

        // 6. Va chercher les derivees
        double cpl, dP_Tl, dHl_Tl, dP_rhol, dHl_rhol, dP_cpl, dHl_cpl;
        F77NAME(FCPLR12)(&un, &Tl, &cpl, &dP_Tl, &dHl_Tl, &dP_rhol, &dHl_rhol, &dP_cpl, &dHl_cpl, &ill, &ivstat, &ierrth);

        val = dP_rhol - dHl_rhol* dP_Tl/dHl_Tl; // NB : this is derivative with T fixed, not with H fixed
      }
  else
    for (auto& val : res)
      {
        int un = 1;
        int ill, ivstat, ierrth;

        // 2. Direct !!
        double cpl, dP_Tl, dHl_Tl, dP_rhol, dHl_rhol, dP_cpl, dHl_cpl;
        F77NAME(FCPLR12)(&un, &T[ind * ncomp + id], &cpl, &dP_Tl, &dHl_Tl, &dP_rhol, &dHl_rhol, &dP_cpl, &dHl_cpl, &ill, &ivstat, &ierrth);

//        val = dP_rhol - dHl_rhol* dP_Tl/dHl_Tl ;
      }

#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_R12_c1_liquide::dT_rho_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  if (like_eos_)
    for (auto& val : res)
      {
        // On fait comme dans EOS : des trucs chelous
        // 1. Enthalpie liquide
        int un = 1;
        double Hl;
        int ill, ivstat, ierrth;
        F77NAME(FPTHLR12)(&un, &T[ind * ncomp + id], &Hl, &ill, &ivstat, &ierrth);

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

        // 5. Inversion des inconnues...
        double Tl, rhol;
        F77NAME(FTLR12)(&un, &P[ind], &tsat, &rhog_sat, &hvsp , &Hl, &Tl, &rhol, &ill, &ivstat, &ierrth);

        // 6. Va chercher les derivees
        double cpl, dP_Tl, dHl_Tl, dP_rhol, dHl_rhol, dP_cpl, dHl_cpl;
        F77NAME(FCPLR12)(&un, &Tl, &cpl, &dP_Tl, &dHl_Tl, &dP_rhol, &dHl_rhol, &dP_cpl, &dHl_cpl, &ill, &ivstat, &ierrth);

        val = dHl_rhol/dHl_Tl;
      }
  else for (auto& val : res)
      {
        // 1. Si facile
        int un = 1;
        int ill, ivstat, ierrth;
        double cpl, dP_Tl, dHl_Tl, dP_rhol, dHl_rhol, dP_cpl, dHl_cpl;
        F77NAME(FCPLR12)(&un, &T[ind * ncomp + id], &cpl, &dP_Tl, &dHl_Tl, &dP_rhol, &dHl_rhol, &dP_cpl, &dHl_cpl, &ill, &ivstat, &ierrth);

        val = dHl_rhol/dHl_Tl;
      }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_R12_c1_liquide::h_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1;
      double hl; //sortie
      int ill, ivstat, ierrth;
      F77NAME(FPTHLR12)(&un, &T[ind * ncomp + id], &hl, &ill, &ivstat, &ierrth);
      val =  hl;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_R12_c1_liquide::dP_h_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  if (like_eos_)
    for (auto& val : res)
      {
        // On fait comme dans EOS : des trucs chelous
        // 1. Enthalpie liquide
        int un = 1;
        double Hl;
        int ill, ivstat, ierrth;
        F77NAME(FPTHLR12)(&un, &T[ind * ncomp + id], &Hl, &ill, &ivstat, &ierrth);

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

        // 5. Inversion des inconnues...
        double Tl, rhol;
        F77NAME(FTLR12)(&un, &P[ind], &tsat, &rhog_sat, &hvsp , &Hl, &Tl, &rhol, &ill, &ivstat, &ierrth);

        // 6. Va chercher les derivees
        double cpl, dP_Tl, dHl_Tl, dP_rhol, dHl_rhol, dP_cpl, dHl_cpl;
        F77NAME(FCPLR12)(&un, &Tl, &cpl, &dP_Tl, &dHl_Tl, &dP_rhol, &dHl_rhol, &dP_cpl, &dHl_cpl, &ill, &ivstat, &ierrth);

        val = -dP_Tl/dHl_Tl;
      }
  else for (auto& val : res)
      {
        int un = 1;
        int ill, ivstat, ierrth;
        // 1. Si beau
        double cpl, dP_Tl, dHl_Tl, dP_rhol, dHl_rhol, dP_cpl, dHl_cpl;
        F77NAME(FCPLR12)(&un, &T[ind * ncomp + id], &cpl, &dP_Tl, &dHl_Tl, &dP_rhol, &dHl_rhol, &dP_cpl, &dHl_cpl, &ill, &ivstat, &ierrth);

        val = -dP_Tl/dHl_Tl;
      }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_R12_c1_liquide::dT_h_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  if (like_eos_)
    for (auto& val : res)
      {
        // On fait comme dans EOS : des trucs chelous
        // 1. Enthalpie liquide
        int un = 1;
        double Hl;
        int ill, ivstat, ierrth;
        F77NAME(FPTHLR12)(&un, &T[ind * ncomp + id], &Hl, &ill, &ivstat, &ierrth);

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

        // 5. Inversion des inconnues...
        double Tl, rhol;
        F77NAME(FTLR12)(&un, &P[ind], &tsat, &rhog_sat, &hvsp , &Hl, &Tl, &rhol, &ill, &ivstat, &ierrth);

        // 6. Va chercher les derivees
        double cpl, dP_Tl, dHl_Tl, dP_rhol, dHl_rhol, dP_cpl, dHl_cpl;
        F77NAME(FCPLR12)(&un, &Tl, &cpl, &dP_Tl, &dHl_Tl, &dP_rhol, &dHl_rhol, &dP_cpl, &dHl_cpl, &ill, &ivstat, &ierrth);

        val = 1./dHl_Tl; // =cpl normalement... mais ici non !


        // 7. On teste la version plus rapide
        double cpl2, dP_Tl2, dHl_Tl2, dP_rhol2, dHl_rhol2, dP_cpl2, dHl_cpl2;
        F77NAME(FCPLR12)(&un, &T[ind * ncomp + id], &cpl2, &dP_Tl2, &dHl_Tl2, &dP_rhol2, &dHl_rhol2, &dP_cpl2, &dHl_cpl2, &ill, &ivstat, &ierrth);
        if  ((abs(1./dHl_Tl)>1.e-10) && (abs((1./dHl_Tl-1./dHl_Tl2)/(1./dHl_Tl)) > 1.e-4 )) Cout << "dT_h2 : " << 1./dHl_Tl2 << " dT_h : " << 1./dHl_Tl << finl ;
      }
  else for (auto& val : res)
      {
        int un = 1;
        int ill, ivstat, ierrth;
        double cpl, dP_Tl, dHl_Tl, dP_rhol, dHl_rhol, dP_cpl, dHl_cpl;
        F77NAME(FCPLR12)(&un, &T[ind * ncomp + id], &cpl, &dP_Tl, &dHl_Tl, &dP_rhol, &dHl_rhol, &dP_cpl, &dHl_cpl, &ill, &ivstat, &ierrth);
        val = 1./dHl_Tl; // =cpl normalement... mais ici non !
      }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_R12_c1_liquide::cp_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1;
      double cpl, dP_Tl, DHl_Tl, DP_Rhol, DHl_Rhol, DP_cpl, DHl_cpl ; //sorties
      int ill, ivstat, ierrth;
      F77NAME(FCPLR12)(&un, &T[ind * ncomp + id], &cpl, &dP_Tl, &DHl_Tl, &DP_Rhol, &DHl_Rhol, &DP_cpl, &DHl_cpl, &ill, &ivstat, &ierrth);
      val =  cpl;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_R12_c1_liquide::beta_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
  VectorD dT_rho___((int )res.size()), rho___((int )res.size());
  dT_rho_(T,P,SpanD(dT_rho___),ncomp,id);
  rho_(T,P,SpanD(rho___),ncomp,id);
  for (auto& val : res) val = dT_rho___[ind] / rho___[ind];
}

void Fluide_R12_c1_liquide::mu_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1;
      double dtl1 = 0., dtl2 = 0.; // We don't need the derivative of lambda
      double mul, dmu1, dmu2;
      F77NAME(FMULR12)(&un, &T[ind * ncomp + id], &dtl1, &dtl2, &mul, &dmu1, &dmu2);
      val =  mul;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

void Fluide_R12_c1_liquide::lambda_(const SpanD T, const SpanD P, SpanD res, int ncomp, int id) const
{
  assert((int )T.size() == ncomp * (int )P.size() && (int )T.size() == ncomp * (int )res.size());
#if HAVE_LIBC3
  /* calcul a saturation */
  for (auto& val : res)
    {
      int un = 1;
      double dtl1 = 0., dtl2 = 0.; // We don't need the derivative of lambda
      double lambdal, dlambl1, dlambl2;
      F77NAME(FCONLR12)(&un, &T[ind * ncomp + id], &dtl1, &dtl2, &lambdal, &dlambl1, &dlambl2);
      val =  lambdal;
    }
#else
  for (auto& val : res) val = 0;
#endif
}

#undef ind
