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
// File:        Op_Diff_K_Eps_Bas_Re_VDF_leaves.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Operateurs/Op_Dift
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Op_Diff_K_Eps_Bas_Re_VDF_leaves_included
#define Op_Diff_K_Eps_Bas_Re_VDF_leaves_included

#include <Op_Diff_K_Eps_Bas_Re_VDF_base.h>
#include <Op_Diff_K_Eps_VDF_base.h>

/// \cond DO_NOT_DOCUMENT
class Op_Diff_K_Eps_Bas_Re_VDF_leaves { };
/// \endcond

class Op_Diff_K_Eps_Bas_Re_VDF_Elem_Axi : public Op_Diff_K_Eps_Bas_Re_VDF_base, public Op_Diff_K_Eps_VDF_Generique<Op_Diff_K_Eps_Bas_Re_VDF_Elem_Axi>
{
  Declare_instanciable_sans_constructeur(Op_Diff_K_Eps_Bas_Re_VDF_Elem_Axi);
public:
  Op_Diff_K_Eps_Bas_Re_VDF_Elem_Axi();
  inline void associer_diffusivite_turbulente() override { associer_diffusivite_turbulente_impl<Eval_Diff_K_Eps_Bas_Re_VDF_const_Elem_Axi>(); }
  inline void associer(const Domaine_dis_base& zd, const Domaine_Cl_dis_base& zcd, const Champ_Inc_base& ch) override { associer_impl<Eval_Diff_K_Eps_Bas_Re_VDF_const_Elem_Axi>(zd,zcd,ch); }
  inline void associer_diffusivite(const Champ_base& ch) override { associer_diffusivite_impl<Eval_Diff_K_Eps_Bas_Re_VDF_const_Elem_Axi>(ch); }
  inline const Champ_Fonc_base& diffusivite_turbulente() const { return diffusivite_turbulente_impl<Eval_Diff_K_Eps_Bas_Re_VDF_const_Elem_Axi>(); }
  inline const Champ_base& diffusivite() const override { return diffusivite_impl<Eval_Diff_K_Eps_Bas_Re_VDF_const_Elem_Axi>(); }
};

class Op_Diff_K_Eps_Bas_Re_VDF_Elem : public Op_Diff_K_Eps_Bas_Re_VDF_base, public Op_Diff_K_Eps_VDF_Generique<Op_Diff_K_Eps_Bas_Re_VDF_Elem>
{
  Declare_instanciable_sans_constructeur(Op_Diff_K_Eps_Bas_Re_VDF_Elem);
public:
  Op_Diff_K_Eps_Bas_Re_VDF_Elem();
  inline void associer_diffusivite_turbulente() override { associer_diffusivite_turbulente_impl<Eval_Diff_K_Eps_Bas_Re_VDF_const_Elem>(); }
  inline void associer(const Domaine_dis_base& zd, const Domaine_Cl_dis_base& zcd, const Champ_Inc_base& ch) override { associer_impl<Eval_Diff_K_Eps_Bas_Re_VDF_const_Elem>(zd,zcd,ch); }
  inline void associer_diffusivite(const Champ_base& ch) override { associer_diffusivite_impl<Eval_Diff_K_Eps_Bas_Re_VDF_const_Elem>(ch); }
  inline const Champ_Fonc_base& diffusivite_turbulente() const { return diffusivite_turbulente_impl<Eval_Diff_K_Eps_Bas_Re_VDF_const_Elem>(); }
  inline const Champ_base& diffusivite() const override { return diffusivite_impl<Eval_Diff_K_Eps_Bas_Re_VDF_const_Elem>(); }
};

///////////////////////////////////
// VAR

class Op_Diff_K_Eps_Bas_Re_VDF_var_Elem : public Op_Diff_K_Eps_Bas_Re_VDF_base, public Op_Diff_K_Eps_VDF_Generique<Op_Diff_K_Eps_Bas_Re_VDF_var_Elem>
{
  Declare_instanciable_sans_constructeur(Op_Diff_K_Eps_Bas_Re_VDF_var_Elem);
public:
  Op_Diff_K_Eps_Bas_Re_VDF_var_Elem();
  double calculer_dt_stab() const override;
  inline void associer_diffusivite_turbulente() override { associer_diffusivite_turbulente_impl<Eval_Diff_K_Eps_Bas_Re_VDF_var_Elem,Eval_Diff_K_Eps_Bas_Re_VDF_const_Elem>(); }
  inline void associer(const Domaine_dis_base& zd, const Domaine_Cl_dis_base& zcd, const Champ_Inc_base& ch) override { associer_impl<Eval_Diff_K_Eps_Bas_Re_VDF_var_Elem>(zd,zcd,ch); }
  inline void associer_diffusivite(const Champ_base& ch) override { associer_diffusivite_impl<Eval_Diff_K_Eps_Bas_Re_VDF_var_Elem>(ch); }
  inline const Champ_Fonc_base& diffusivite_turbulente() const { return diffusivite_turbulente_impl<Eval_Diff_K_Eps_Bas_Re_VDF_var_Elem>(); }
  inline const Champ_base& diffusivite() const override { return diffusivite_impl<Eval_Diff_K_Eps_Bas_Re_VDF_var_Elem>(); }
  inline void completer() override
  {
    Op_Diff_K_Eps_Bas_Re_VDF_base::completer();
    Eval_Diff_K_Eps_Bas_Re_VDF_var_Elem& eval_diff_turb = static_cast<Eval_Diff_K_Eps_Bas_Re_VDF_var_Elem&>(iter->evaluateur());
    eval_diff_turb.associer_pb(equation().probleme());
  }
};

#endif /* Op_Diff_K_Eps_Bas_Re_VDF_leaves_included */
