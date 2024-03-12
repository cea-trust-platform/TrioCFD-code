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
// File      : OpVDFFacentre4b.h
// Directory : $IJK_ROOT/src/IJK/OpVDF
//
/////////////////////////////////////////////////////////////////////////////


#ifndef Op_VDF_Face_centre4b_H
#define Op_VDF_Face_centre4b_H

#include <Op_Conv_VDF_base.h>
#include <ItVDFFa.h>
#include <Evcentre4bVDFFa.h>

/*! @brief class Op_Conv_centre4b_VDF_Face
 *
 *   Cette classe represente l'operateur de convection associe a une equation de
 *   la quantite de mouvement.
 *   La discretisation est VDF
 *   Le champ convecte est de type Champ_Face
 *   Le schema de convection est du type centre4b (centre sur 4 points)
 *   L'iterateur associe est de type Iterateur_VDF_Face
 *   L'evaluateur associe est de type Eval_centre4b_VDF_Face
 *
 *
 */
declare_It_VDF_Face(Eval_centre4b_VDF_Face)

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Op_Conv_centre4b_VDF_Face
//
//////////////////////////////////////////////////////////////////////////////


class Op_Conv_centre4b_VDF_Face : public Op_Conv_VDF_base
{

  Declare_instanciable_sans_constructeur(Op_Conv_centre4b_VDF_Face);

public:

  Op_Conv_centre4b_VDF_Face();
  void associer(const Domaine_dis& , const Domaine_Cl_dis& , const Champ_Inc& );
  inline void associer_vitesse(const Champ_base&) ;
  const Champ_Inc_base& vitesse() const;
  Champ_Inc_base& vitesse();
};

/*! @brief associe le champ de vitesse a l'evaluateur
 *
 */
inline void Op_Conv_centre4b_VDF_Face::associer_vitesse(const Champ_base& ch_vit)
{
  const Champ_Face_VDF& vit = (Champ_Face_VDF&) ch_vit;

  Eval_centre4b_VDF_Face& eval_conv = (Eval_centre4b_VDF_Face&) iter.evaluateur();
  eval_conv.associer(vit );                // Eval_Conv_VDF::associer
}

#endif
