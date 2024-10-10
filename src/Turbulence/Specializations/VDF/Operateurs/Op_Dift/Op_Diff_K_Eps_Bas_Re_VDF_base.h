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
// File:        Op_Diff_K_Eps_Bas_Re_VDF_base.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Operateurs/Op_Dift
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Op_Diff_K_Eps_Bas_Re_VDF_base_included
#define Op_Diff_K_Eps_Bas_Re_VDF_base_included

#include <Modele_turbulence_hyd_K_Eps_Bas_Reynolds.h>
#include <Modele_turbulence_hyd_K_Eps_2_Couches.h>
#include <Eval_Diff_K_Eps_Bas_Re_VDF_leaves.h>
#include <Op_Diff_K_Eps_Bas_Re_base.h>
#include <Iterateur_VDF_Elem.h>
#include <Statistiques.h>
#include <Op_VDF_Elem.h>




extern Stat_Counter_Id diffusion_counter_;

class Op_Diff_K_Eps_Bas_Re_VDF_base : public Op_Diff_K_Eps_Bas_Re_base, public Op_VDF_Elem
{
  Declare_base(Op_Diff_K_Eps_Bas_Re_VDF_base);
public:
  Op_Diff_K_Eps_Bas_Re_VDF_base(const Iterateur_VDF_base& iter_base) { iter = iter_base; }

  void completer() override;

  inline void mettre_a_jour_diffusivite() const  { /* do nothing */ }
  inline void modifier_pour_Cl(Matrice_Morse& matrice, DoubleTab& secmem) const  override { Op_VDF_Elem::modifier_pour_Cl(iter->domaine(), iter->domaine_Cl(), matrice, secmem); }
  inline DoubleTab& calculer(const DoubleTab& inco, DoubleTab& resu) const override { return iter->calculer(inco, resu); }

  inline void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const override
  {
    const std::string& nom_inco = equation().inconnue().le_nom().getString();
    Matrice_Morse *mat = matrices.count(nom_inco) ? matrices.at(nom_inco) : nullptr, mat2;
    if (!mat)
      return;
    Op_VDF_Elem::dimensionner(iter->domaine(), iter->domaine_Cl(), mat2);
    mat->nb_colonnes() ? *mat += mat2 : *mat = mat2;
  }

  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const override
  {
    statistiques().begin_count(diffusion_counter_);
    iter->ajouter_blocs(matrices,secmem,semi_impl);
    statistiques().end_count(diffusion_counter_);
  }
  inline int has_interface_blocs() const override { return 1; }

  inline OWN_PTR(Iterateur_VDF_base)& get_iter() { return iter; }

  template <typename EVAL_TYPE>
  void associer_diffusivite_impl(const Champ_base& ch_diff)
  {
    EVAL_TYPE& eval_diff_turb = static_cast<EVAL_TYPE&> (iter->evaluateur());
    eval_diff_turb.associer(ch_diff);
  }

  template <typename EVAL_TYPE>
  const Champ_Fonc_base& diffusivite_turbulente_impl() const
  {
    const EVAL_TYPE& eval_diff = static_cast<const EVAL_TYPE&> (iter->evaluateur());
    return eval_diff.diffusivite_turbulente();
  }

  template <typename EVAL_TYPE>
  const Champ_base& diffusivite_impl() const
  {
    const EVAL_TYPE& eval_diff_turb = static_cast<const EVAL_TYPE&> (iter->evaluateur());
    return eval_diff_turb.diffusivite();
  }

  template <typename EVAL_TYPE, typename EVAL_TYPE2 = EVAL_TYPE /* for var only */>
  void associer_diffusivite_turbulente_impl();

protected:
  OWN_PTR(Iterateur_VDF_base) iter;
};

template <typename EVAL_TYPE, typename EVAL_TYPE2>
void Op_Diff_K_Eps_Bas_Re_VDF_base::associer_diffusivite_turbulente_impl()
{
  assert(mon_equation.non_nul());
  if(sub_type(Transport_K_KEps,mon_equation.valeur()))
    {
      const Transport_K_KEps& eqn_transport = ref_cast(Transport_K_KEps,mon_equation.valeur());
      const Modele_turbulence_hyd_K_Eps_2_Couches& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_2_Couches,eqn_transport.modele_turbulence());
      const Champ_Fonc_base& diff_turb = mod_turb.viscosite_turbulente();
      EVAL_TYPE& eval_diff = static_cast<EVAL_TYPE&> (iter->evaluateur());
      eval_diff.associer_diff_turb(diff_turb);
    }
  else if(sub_type(Transport_K_Eps_Bas_Reynolds,mon_equation.valeur()))
    {
      const Transport_K_Eps_Bas_Reynolds& eqn_transport = ref_cast(Transport_K_Eps_Bas_Reynolds,mon_equation.valeur());
      const Modele_turbulence_hyd_K_Eps_Bas_Reynolds& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Bas_Reynolds,eqn_transport.modele_turbulence());
      const Champ_Fonc_base& diff_turb = mod_turb.viscosite_turbulente();
      EVAL_TYPE& eval_diff = static_cast<EVAL_TYPE&> (iter->evaluateur());
      eval_diff.associer_diff_turb(diff_turb);
    }
  else if(sub_type(Transport_K_Eps,mon_equation.valeur()))
    {
      const Transport_K_Eps& eqn_transport = ref_cast(Transport_K_Eps,mon_equation.valeur());
      const Modele_turbulence_hyd_K_Eps& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps,eqn_transport.modele_turbulence());
      const Champ_Fonc_base& diff_turb = mod_turb.viscosite_turbulente();
      EVAL_TYPE2& eval_diff = static_cast<EVAL_TYPE2&> (iter->evaluateur());
      eval_diff.associer_diff_turb(diff_turb);
    }
}

#endif /* Op_Diff_K_Eps_Bas_Re_VDF_base_included */
