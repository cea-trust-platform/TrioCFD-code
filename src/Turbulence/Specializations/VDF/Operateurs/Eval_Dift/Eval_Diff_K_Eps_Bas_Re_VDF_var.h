/****************************************************************************
* Copyright (c) 2019, CEA
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
// File:        Eval_Diff_K_Eps_Bas_Re_VDF_var.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Operateurs/Eval_Dift
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Eval_Diff_K_Eps_Bas_Re_VDF_var_included
#define Eval_Diff_K_Eps_Bas_Re_VDF_var_included

#include <Eval_Diff_VDF.h>
#include <Champ_Fonc.h>
#include <TRUST_Ref.h>

class Champ_Uniforme;

class  Eval_Diff_K_Eps_Bas_Re_VDF_var : public Eval_Diff_VDF
{
public:
  Eval_Diff_K_Eps_Bas_Re_VDF_var(double Prandt_K = PRDT_K_DEFAUT , double Prandt_Eps =PRDT_EPS_DEFAUT ) : Prdt_K(Prandt_K) , Prdt_Eps(Prandt_Eps)
  {
    Prdt[0]=Prandt_K;
    Prdt[1]=Prandt_Eps;
  }

  inline void associer_diff_turb(const Champ_Fonc& diffu) { diffusivite_turbulente_ = diffu; }
  inline const Champ_Fonc& diffusivite_turbulente() const { return diffusivite_turbulente_.valeur(); }
  inline void mettre_a_jour() override { dv_diffusivite_turbulente.ref(diffusivite_turbulente_->valeur().valeurs()); }

  // pour CRTP
  inline double nu_1_impl(int i, int compo) const
  {
    return tab_diffusivite_(i) + dv_diffusivite_turbulente(i)/Prdt[compo];
  }

  inline double compute_heq_impl(double d0, int i, double d1, int j, int compo) const
  {
    return 0.5*(nu_1_impl(i,compo) + nu_1_impl(j,compo))/(d0+d1);
  }

protected:
  static constexpr double PRDT_K_DEFAUT = 1.0, PRDT_EPS_DEFAUT = 1.3;
  double Prdt_K, Prdt_Eps, Prdt[2];
  REF(Champ_Fonc) diffusivite_turbulente_;
  DoubleVect dv_diffusivite_turbulente;
};

#endif /* Eval_Diff_K_Eps_Bas_Re_VDF_var_included */
