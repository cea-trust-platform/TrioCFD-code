/****************************************************************************
* Copyright (c) 2023, CEA
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
// File:        Eval_Diff_K_Omega_VDF_leaves.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Operateurs/Eval_Dift
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Eval_Diff_K_Omega_VDF_leaves_included
#define Eval_Diff_K_Omega_VDF_leaves_included

#include <Eval_Diff_K_Omega_VDF_const.h>
#include <Eval_Diff_K_Omega_VDF_var.h>
#include <Eval_Diff_VDF_Elem_Gen.h>

/// \cond DO_NOT_DOCUMENT
class Eval_Diff_K_Omega_VDF_leaves { };
/// \endcond

// const
class Eval_Diff_K_Omega_VDF_Elem : public Eval_Diff_VDF_Elem_Gen<Eval_Diff_K_Omega_VDF_Elem>, public Eval_Diff_K_Omega_VDF_const
{
public:
  static constexpr bool IS_MULTD = false, IS_RANS = true;
};

// class Eval_Diff_K_Omega_VDF_Elem_Axi : public Eval_Diff_VDF_Elem_Gen<Eval_Diff_K_Omega_VDF_Elem_Axi>, public Eval_Diff_K_ou_Omega_VDF_const
// {
// public:
//   static constexpr bool IS_MULTD = false, IS_AXI = true, IS_RANS = true;
// };

// class Eval_Diff_K_Omega_QC_VDF_Elem : public Eval_Diff_VDF_Elem_Gen<Eval_Diff_K_Omega_QC_VDF_Elem>, public Eval_Diff_K_Omega_VDF_const
// {
// public:
//   static constexpr bool IS_MULTD = false, IS_QUASI = true, IS_RANS = true;
// };

// var
class Eval_Diff_K_Omega_var_VDF_Elem : public Eval_Diff_VDF_Elem_Gen<Eval_Diff_K_Omega_var_VDF_Elem>, public Eval_Diff_K_Omega_VDF_var
{
public:
  static constexpr bool IS_RANS = true;
};

// class Eval_Diff_K_Omega_QC_var_VDF_Elem : public Eval_Diff_VDF_Elem_Gen<Eval_Diff_K_Omega_QC_var_VDF_Elem>, public Eval_Diff_K_Omega_VDF_var
// {
// public:
//   static constexpr bool IS_MULTD = false, IS_QUASI = true, IS_RANS = true;
// };

// Bicephale
// class Eval_Diff_K_VDF_Elem : public Eval_Diff_VDF_Elem_Gen<Eval_Diff_K_VDF_Elem>, public Eval_Diff_K_ou_Omega_VDF_const
// {
// public:
//   static constexpr bool IS_MULTD = false, IS_RANS = true;
// };

// class Eval_Diff_Omega_VDF_Elem : public Eval_Diff_VDF_Elem_Gen<Eval_Diff_Omega_VDF_Elem>, public Eval_Diff_K_ou_Omega_VDF_const
// {
// public:
//   static constexpr bool IS_MULTD = false, IS_RANS = true;
// };

#endif /* Eval_Diff_K_Omega_VDF_leaves_included */
