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
// File      : Op_Conv_QuickIJK_VDF_Face.h
// Directory : $IJK_ROOT/src/IJK/OpVDF
//
/////////////////////////////////////////////////////////////////////////////
#ifndef Op_Conv_QuickIJK_VDF_Face_H
#define Op_Conv_QuickIJK_VDF_Face_H

#include <Op_Conv_Quick_VDF_Face.h>
#include <OpQuickIJK.h>
#include <VDF_to_IJK.h>

class Op_Conv_QuickIJK_VDF_Face : public Op_Conv_Quick_VDF_Face
{
  Declare_instanciable(Op_Conv_QuickIJK_VDF_Face);
public:
  //void associer(const Domaine_dis& , const Domaine_Cl_dis& , const Champ_Inc& );
  void completer();
  void associer_vitesse(const Champ_base&) ;
  //const Champ_Inc_base& vitesse() const;
  //Champ_Inc_base& vitesse();

  DoubleTab& ajouter(const DoubleTab& inco, DoubleTab& resu) const;
  DoubleTab& calculer(const DoubleTab& inco, DoubleTab& resu ) const;
  //double calculer_dt_stab() const;

protected:
  mutable OpConvQuickIJK_double op_ijk_;

  REF(Champ_base) vconv_;
};
#endif
