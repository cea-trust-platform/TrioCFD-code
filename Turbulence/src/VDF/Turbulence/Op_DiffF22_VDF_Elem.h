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
// File:        Op_DiffF22_VDF_Elem.h
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Op_DiffF22_VDF_Elem_included
#define Op_DiffF22_VDF_Elem_included

#include <Eval_Diff_VDF_Elem_leaves.h>
#include <Op_Diff_VDF_base.h>
#include <Op_VDF_Elem.h>
#include <ItVDFEl.h>

//
// .DESCRIPTION class Op_DiffF22_VDF_Elem
//
//  Cette classe represente l'operateur de diffusion associe a l'equation de la quantite F22 dans le modele V2F
//  La discretisation est VDF
//  Le champ diffuse est scalaire
//  Le champ de diffusivite est uniforme
//  L'iterateur associe est de type Iterateur_VDF_Elem
//  L'evaluateur associe est de type Eval_DiffF22_VDF_const_Elem

//
// .SECTION voir aussi
//
//

declare_It_VDF_Elem(Eval_DiffF22_VDF_const_Elem)

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Op_DiffF22_VDF_Elem
//
//////////////////////////////////////////////////////////////////////////////

class Op_DiffF22_VDF_Elem : public Op_Diff_VDF_base, public Op_VDF_Elem
{

  Declare_instanciable_sans_constructeur(Op_DiffF22_VDF_Elem);

public:

  Op_DiffF22_VDF_Elem();
  void associer(const Zone_dis& , const Zone_Cl_dis& ,const Champ_Inc& ) override;
  void associer_diffusivite(const Champ_base& ) override;
  const Champ_base& diffusivite() const override;
  inline  void dimensionner(Matrice_Morse& ) const override;
  inline void modifier_pour_Cl(Matrice_Morse&, DoubleTab&) const override;
  double calculer_dt_stab() const override;

};

// Description:
// on dimensionne notre matrice.
inline  void Op_DiffF22_VDF_Elem::dimensionner(Matrice_Morse& matrice) const
{
  Op_VDF_Elem::dimensionner(iter.zone(), iter.zone_Cl(), matrice);
}

inline void Op_DiffF22_VDF_Elem::modifier_pour_Cl(Matrice_Morse& matrice, DoubleTab& secmem) const
{
  Op_VDF_Elem::modifier_pour_Cl(iter.zone(), iter.zone_Cl(), matrice, secmem);
}

#endif
