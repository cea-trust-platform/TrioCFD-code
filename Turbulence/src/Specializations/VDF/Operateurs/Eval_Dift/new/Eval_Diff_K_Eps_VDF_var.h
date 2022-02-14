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
// File:        Eval_Diff_K_Eps_VDF_var.h
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Eval_Diff_K_Eps_VDF_var_included
#define Eval_Diff_K_Eps_VDF_var_included

#include <Eval_Diff_K_Eps_VDF.h>
#include <Champ_Don.h>

class Eval_Diff_K_Eps_VDF_var : public Eval_Diff_K_Eps_VDF
{
public:

  inline void associer(const Champ_Don& diffu)
  {
    diffusivite_ = diffu.valeur();
    dv_diffusivite.ref(diffu.valeurs());
  }

  inline virtual void mettre_a_jour()
  {
    (diffusivite_->valeurs().echange_espace_virtuel());
    dv_diffusivite.ref(diffusivite_->valeurs());
    dv_diffusivite_turbulente.ref(diffusivite_turbulente_->valeurs());
  }

  // Methods used by the flux computation in template class:
  inline double nu_1_impl(int i, int compo) const
  {
    return dv_diffusivite(i) + dv_diffusivite_turbulente(i)/Prdt[compo];
  }

  inline double nu_2_impl(int i, int compo) const { return nu_1_impl(i,compo); }

  inline double compute_heq_impl(double d0, int i, double d1, int j, int compo) const
  {
    return 0.5*(nu_1_impl(i,compo)+nu_1_impl(j,compo))/(d0+d1);
  }

protected:
  DoubleVect dv_diffusivite;
};

#endif /* Eval_Diff_K_Eps_VDF_var_included */
