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
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Op_Diff_K_Eps_Bas_Re_VDF_base_included
#define Op_Diff_K_Eps_Bas_Re_VDF_base_included

#include <Eval_Diff_K_Eps_Bas_Re_VDF_leaves.h>
#include <Op_Diff_K_Eps_Bas_Re_base.h>
#include <Op_VDF_Elem.h>
#include <ItVDFEl.h>

class Zone_Cl_dis;
class Zone_dis;
class Champ_Inc;

class Op_Diff_K_Eps_Bas_Re_VDF_base : public Op_Diff_K_Eps_Bas_Re_base, public Op_VDF_Elem
{
  Declare_base(Op_Diff_K_Eps_Bas_Re_VDF_base);
public:
  Op_Diff_K_Eps_Bas_Re_VDF_base(const Iterateur_VDF_base& iter_base) : iter(iter_base) { }

  void completer();

  inline void mettre_a_jour_diffusivite() const { /* do nothing */ }
  inline void contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& matrice) const { iter->ajouter_contribution(inco, matrice); }
  inline void contribuer_au_second_membre(DoubleTab& resu) const { iter->contribuer_au_second_membre(resu); }
  inline void dimensionner(Matrice_Morse& matrice) const { Op_VDF_Elem::dimensionner(iter->zone(), iter->zone_Cl(), matrice); }
  inline void modifier_pour_Cl(Matrice_Morse& matrice, DoubleTab& secmem) const { Op_VDF_Elem::modifier_pour_Cl(iter->zone(), iter->zone_Cl(), matrice, secmem); }
  inline DoubleTab& ajouter(const DoubleTab& inco, DoubleTab& resu) const { return iter->ajouter(inco, resu); }
  inline DoubleTab& calculer(const DoubleTab& inco, DoubleTab& resu) const { return iter->calculer(inco, resu); }

  inline Iterateur_VDF& get_iter() { return iter; }
  inline const Iterateur_VDF& get_iter() const { return iter; }

protected:
  Iterateur_VDF iter;
};

// one for all !
template<typename OP_TYPE>
class Op_Diff_K_Eps_Bas_Re_VDF_Generique
{
protected:

  template <typename EVAL_TYPE>
  inline void associer_impl(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_cl_dis, const Champ_Inc& ch_diffuse)
  {
    const Champ_P0_VDF& inco = ref_cast(Champ_P0_VDF,ch_diffuse.valeur());
    const Zone_VDF& zvdf = ref_cast(Zone_VDF,zone_dis.valeur());
    const Zone_Cl_VDF& zclvdf = ref_cast(Zone_Cl_VDF,zone_cl_dis.valeur());
    iter_()->associer(zvdf, zclvdf,static_cast<OP_TYPE&>(*this));
    EVAL_TYPE& eval_diff = static_cast<EVAL_TYPE&> (iter_()->evaluateur());
    eval_diff.associer_zones(zvdf, zclvdf );
    eval_diff.associer_inconnue(inco );
  }

  template <typename EVAL_TYPE>
  void associer_diffusivite_impl(const Champ_base& ch_diff)
  {
    EVAL_TYPE& eval_diff_turb = static_cast<EVAL_TYPE&> (iter_()->evaluateur());
    eval_diff_turb.associer(ch_diff);
  }

  template <typename EVAL_TYPE>
  const Champ_Fonc& diffusivite_turbulente_impl() const
  {
    const EVAL_TYPE& eval_diff = static_cast<const EVAL_TYPE&> (iter_()->evaluateur());
    return eval_diff.diffusivite_turbulente();
  }

  template <typename EVAL_TYPE>
  const Champ_base& diffusivite_impl() const
  {
    const EVAL_TYPE& eval_diff_turb = static_cast<const EVAL_TYPE&> (iter_()->evaluateur());
    return eval_diff_turb.diffusivite();
  }

private:
  inline Iterateur_VDF& iter_() { return static_cast<OP_TYPE *>(this)->get_iter(); } // CRTP pour recuperer l'iter
  inline const Iterateur_VDF& iter_() const { return static_cast<const OP_TYPE *>(this)->get_iter(); } // CRTP pour recuperer l'iter
};


#endif /* Op_Diff_K_Eps_Bas_Re_VDF_base_included */
