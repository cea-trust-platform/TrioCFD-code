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
// File      : Op_Dift_IJK_VDF_Face.h
// Directory : $IJK_ROOT/src/IJK/OpVDF
//
/////////////////////////////////////////////////////////////////////////////
#ifndef Op_Dift_IJK_VDF_Face_H
#define Op_Dift_IJK_VDF_Face_H

#include <Op_Dift_VDF_var_Face.h>
#include <OpDiffTurbIJK.h>
#include <VDF_to_IJK.h>
#include <Boundary_Conditions.h>


class Op_Dift_IJK_VDF_Face : public Op_Dift_VDF_var_Face
{
  Declare_instanciable(Op_Dift_IJK_VDF_Face);
public:
  void associer_diffusivite(const Champ_base& nu);
  void completer();
  DoubleTab& ajouter(const DoubleTab& inco, DoubleTab& resu) const;

protected:
  mutable OpDiffTurbIJK_double op_ijk_;

  // Reference to diffusivity and k_energy in vdf storage format
  REF(Champ_base) vdf_nu_;
  REF(Champ_base) vdf_nut_;
  REF(Champ_base) vdf_kenergy_;

  Boundary_Conditions bc_;

  mutable IJK_Field_double vx_, vy_, vz_, nu_, nut_, kenergy_, dvx_, dvy_, dvz_;
};
#endif
