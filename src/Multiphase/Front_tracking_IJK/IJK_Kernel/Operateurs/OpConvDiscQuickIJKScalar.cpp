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

#include <OpConvDiscQuickIJKScalar.h>

Implemente_instanciable_sans_constructeur(OpConvDiscQuickIJKScalar_double, "OpConvDiscQuickIJKScalar_double", Operateur_IJK_elem_conv_base_double);

Sortie& OpConvDiscQuickIJKScalar_double::printOn(Sortie& os) const
{
  return os;
}

Entree& OpConvDiscQuickIJKScalar_double::readOn(Entree& is)
{
  return is;
}

void OpConvDiscQuickIJKScalar_double::calculer(const IJK_Field_double& field,
                                               const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                               IJK_Field_double& result)
{
  Operateur_IJK_elem_conv_base_double::calculer(field, vx, vy, vz, result);
  input_indicatrice_ = 0;

}

void OpConvDiscQuickIJKScalar_double::ajouter(const IJK_Field_double& field,
                                              const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                              IJK_Field_double& result)
{
  Operateur_IJK_elem_conv_base_double::ajouter(field, vx, vy, vz, result);
  input_indicatrice_ = 0;

}

//void OpConvDiscQuickIJKScalar_double::initialize(const IJK_Splitting& splitting)
//{
//  Operateur_IJK_elem_conv_base_double::initialize(splitting, indicatrice);
//  input_indicatrice_ = 0;
//}
