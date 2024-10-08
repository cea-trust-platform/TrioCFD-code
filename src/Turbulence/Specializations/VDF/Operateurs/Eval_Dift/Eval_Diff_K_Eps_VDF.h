/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
// File:        Eval_Diff_K_Eps_VDF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Operateurs/Eval_Dift
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Eval_Diff_K_Eps_VDF_included
#define Eval_Diff_K_Eps_VDF_included

#include <Champ_Fonc_base.h>
#include <Champ_Uniforme.h>
#include <Champ_base.h>
#include <TRUST_Ref.h>


class Eval_Diff_K_Eps_VDF
{
public:
  virtual ~Eval_Diff_K_Eps_VDF() { }

  Eval_Diff_K_Eps_VDF(double Prandt_K = PRDT_K_DEFAUT, double Prandt_Eps = PRDT_EPS_DEFAUT ) : Prdt_K(Prandt_K) , Prdt_Eps(Prandt_Eps), db_diffusivite(-123.)
  {
    Prdt[0]=Prandt_K;
    Prdt[1]=Prandt_Eps;
  }

  inline void associer_diff_turb(const Champ_Fonc_base& diffu) { diffusivite_turbulente_ = diffu; }

  inline void associer_mvolumique(const Champ_base& mvol)
  {
    masse_volumique_ = mvol;
    dv_mvol.ref(mvol.valeurs());
  }

  inline void associer_Pr_K_Eps(double Pr_K,double Pr_Eps)
  {
    Prdt_K = Pr_K;
    Prdt_Eps = Pr_Eps;
    Prdt[0] = Pr_K;
    Prdt[1] = Pr_Eps;
  }

  inline void associer(const Champ_base& diffu)
  {
    diffusivite_ =  diffu;
    if (sub_type(Champ_Uniforme, diffu))  db_diffusivite = diffu.valeurs()(0,0);
  }

  inline virtual void mettre_a_jour()
  {
    dv_diffusivite_turbulente.ref(diffusivite_turbulente_->valeurs());
    if (sub_type(Champ_Uniforme, diffusivite_.valeur())) db_diffusivite = diffusivite_->valeurs()(0,0);
  }

  inline const Champ_Fonc_base& diffusivite_turbulente() const { return diffusivite_turbulente_.valeur(); }
  inline const Champ_base& diffusivite() const { return diffusivite_.valeur(); }

  // Pour CRTP !
  inline int get_ind_Fluctu_Term() const { throw; }
  inline double get_equivalent_distance(int boundary_index,int local_face) const { throw; }
  inline double get_dv_mvol(const int i) const { return dv_mvol[i]; }

protected:
  static constexpr double PRDT_K_DEFAUT = 1.0, PRDT_EPS_DEFAUT = 1.3;
  double Prdt_K, Prdt_Eps, db_diffusivite, Prdt[2];
  DoubleVect dv_diffusivite_turbulente, dv_mvol;
  OBS_PTR(Champ_Fonc_base) diffusivite_turbulente_;
  OBS_PTR(Champ_base) masse_volumique_, diffusivite_;
};

#endif /* Eval_Diff_K_Eps_VDF_included */
