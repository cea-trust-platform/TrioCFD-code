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
// File:        Op_Diff_K_Eps_V2_VDF_Elem.h
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Op_Diff_K_Eps_V2_VDF_Elem_included
#define Op_Diff_K_Eps_V2_VDF_Elem_included

#include <Eval_Diff_K_Eps_V2_VDF_const_Elem.h>
#include <Op_Diff_K_Eps_Bas_Re_VDF_base.h>

declare_It_VDF_Elem(Eval_Diff_K_Eps_V2_VDF_const_Elem)
class Op_Diff_K_Eps_V2_VDF_Elem : public Op_Diff_K_Eps_Bas_Re_VDF_base, public Op_Diff_K_Eps_Bas_Re_VDF_Generique<Op_Diff_K_Eps_V2_VDF_Elem>
{
  Declare_instanciable_sans_constructeur(Op_Diff_K_Eps_V2_VDF_Elem);
public:
  Op_Diff_K_Eps_V2_VDF_Elem();
  void associer_diffusivite_turbulente();

  inline void associer(const Zone_dis& zd, const Zone_Cl_dis& zcd, const Champ_Inc& ch) { associer_impl<Eval_Diff_K_Eps_V2_VDF_const_Elem>(zd,zcd,ch); }
  inline void associer_diffusivite(const Champ_base& ch) { associer_diffusivite_impl<Eval_Diff_K_Eps_V2_VDF_const_Elem>(ch); }
  inline const Champ_Fonc& diffusivite_turbulente() const { return diffusivite_turbulente_impl<Eval_Diff_K_Eps_V2_VDF_const_Elem>(); }
  inline const Champ_base& diffusivite() const { return diffusivite_impl<Eval_Diff_K_Eps_V2_VDF_const_Elem>(); }
};

#endif /* Op_Diff_K_Eps_V2_VDF_Elem_included */
