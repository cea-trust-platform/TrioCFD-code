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
/////////////////////////////////////////////////////////////////////////////
//
// File      : OpVDFElQuickb.h
// Directory : $IJK_ROOT/src/IJK/OpVDF
//
/////////////////////////////////////////////////////////////////////////////


#ifndef Op_VDF_Elem_Quickb_H
#define Op_VDF_Elem_Quickb_H

#include <ItVDFEl.h>
#include <EvQuickbVDFEl.h>
#include <Op_VDF_Elem.h>
/*! @brief class Op_Conv_Quickb_VDF_Elem
 *
 *   Cette classe represente l'operateur de convection associe a une equation de
 *   transport d'un scalaire.
 *   La discretisation est VDF
 *   Le champ convecte est scalaire
 *   Le schema de convection est du type Quickb
 *   L'iterateur associe est de type Iterateur_VDF_Elem
 *   L'evaluateur associe est de type Eval_Quickb_VDF_Elem
 *
 *
 */
declare_It_VDF_Elem(Eval_Quickb_VDF_Elem)

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Op_Conv_Quickb_VDF_Elem
//
//////////////////////////////////////////////////////////////////////////////

class Op_Conv_Quickb_VDF_Elem : public Op_Conv_VDF_base, Op_VDF_Elem
{
  Declare_instanciable_sans_constructeur(Op_Conv_Quickb_VDF_Elem);

public:

  Op_Conv_Quickb_VDF_Elem();
  void associer(const Domaine_dis& , const Domaine_Cl_dis& , const Champ_Inc&  );
  inline void associer_vitesse(const Champ_base& );
  const Champ_Inc_base& vitesse() const;
  Champ_Inc_base& vitesse();
  inline void dimensionner(Matrice_Morse& ) const;
  inline void modifier_pour_Cl(Matrice_Morse&, DoubleTab&) const;

protected:
  inline Op_Conv_Quickb_VDF_Elem(const Iterateur_VDF_base&);
};

//
// Fonctions inline de la classe Op_Conv_Quickb_VDF_Elem
//


// Ce constructeur permet de creer des classes filles (exemple : front_tracking)
inline Op_Conv_Quickb_VDF_Elem::Op_Conv_Quickb_VDF_Elem(const Iterateur_VDF_base& it)
  : Op_Conv_VDF_base(it)
{
}

/*! @brief associe le champ de vitesse a l'evaluateur
 *
 */
inline void Op_Conv_Quickb_VDF_Elem::associer_vitesse(const Champ_base& ch_vit)
{
  const Champ_Face_VDF& vit = (Champ_Face_VDF&) ch_vit;

  Eval_Quickb_VDF_Elem& eval_conv = (Eval_Quickb_VDF_Elem&) iter.evaluateur();
  eval_conv.associer(vit );                // Eval_Conv_VDF::associer
}

inline void Op_Conv_Quickb_VDF_Elem::dimensionner(Matrice_Morse& matrice) const
{
  Op_VDF_Elem::dimensionner(iter.domaine(), iter.domaine_Cl(), matrice);
}

inline void Op_Conv_Quickb_VDF_Elem::modifier_pour_Cl(Matrice_Morse& matrice, DoubleTab& secmem) const
{
  Op_VDF_Elem::modifier_pour_Cl(iter.domaine(), iter.domaine_Cl(), matrice, secmem);
}

#endif

