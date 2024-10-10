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
// File:        Op_Diff_K_Omega_VDF_base.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Operateurs/Op_Dift
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Op_Diff_K_Omega_VDF_base_included
#define Op_Diff_K_Omega_VDF_base_included

#include <Eval_Diff_K_Omega_VDF_leaves.h>
#include <Op_Diff_K_Omega_base.h>
#include <Iterateur_VDF_Elem.h>
#include <Op_VDF_Elem.h>




class Op_Diff_K_Omega_VDF_base : public Op_Diff_K_Omega_base, public Op_VDF_Elem
{
  Declare_base(Op_Diff_K_Omega_VDF_base);
public:
  Op_Diff_K_Omega_VDF_base(const Iterateur_VDF_base& iter_base) { iter = iter_base; }

  void completer() override;
  void associer_diffusivite(const Champ_base& ch_diff) override;
  void associer_diffusivite_turbulente() override;
  void modifier_pour_Cl(Matrice_Morse&, DoubleTab&) const override;
  const Champ_Fonc_base& diffusivite_turbulente() const;
  const Champ_base& diffusivite() const override;

  virtual inline void mettre_a_jour_diffusivite() const { assert(mon_equation.non_nul()); }

  inline DoubleTab& calculer(const DoubleTab& inco, DoubleTab& resu) const override
  {
    mettre_a_jour_diffusivite();
    return iter->calculer(inco, resu);
  }

  inline OWN_PTR(Iterateur_VDF_base)& get_iter() { return iter; }

  inline int has_interface_blocs() const override { return 1; }

  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const override;
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const override;

protected:
  OWN_PTR(Iterateur_VDF_base) iter;
};

// one for all !
template<typename OP_TYPE>
class Op_Diff_K_Omega_VDF_Generique
{
protected:

  template <typename EVAL_TYPE>
  inline void associer_impl(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_cl_dis, const Champ_Inc_base& ch_diffuse)
  {
    const Champ_P0_VDF& inco = ref_cast(Champ_P0_VDF,ch_diffuse);
    const Domaine_VDF& zvdf = ref_cast(Domaine_VDF,domaine_dis);
    const Domaine_Cl_VDF& zclvdf = ref_cast(Domaine_Cl_VDF,domaine_cl_dis);
    iter_()->associer(zvdf, zclvdf,static_cast<OP_TYPE&>(*this));
    EVAL_TYPE& eval_diff = static_cast<EVAL_TYPE&> (iter_()->evaluateur());
    eval_diff.associer_domaines(zvdf, zclvdf );
    eval_diff.associer_inconnue(inco );
  }

private:
  inline OWN_PTR(Iterateur_VDF_base)& iter_() { return static_cast<OP_TYPE *>(this)->get_iter(); } // CRTP pour recuperer l'iter
};

#endif /* Op_Diff_K_Omega_VDF_base_included */
