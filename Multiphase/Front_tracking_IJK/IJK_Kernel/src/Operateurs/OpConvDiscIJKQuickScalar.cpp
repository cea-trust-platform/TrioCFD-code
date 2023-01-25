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

#include <OpConvDiscIJKQuickScalar.h>


void OpConvDiscIJKQuickScalar_double::initialize(const IJK_Splitting& splitting)
{
  OpConvIJKElemCommon_double::initialize(splitting);
  input_indicatrice_ = 0;
}

void OpConvDiscIJKQuickScalar_double::calculer(const IJK_Field_double& field,
                                               const IJK_Field_double& ref_ijk,
                                               const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                               IJK_Field_double& result)
{
  // Si ce test plante, c'est qu'on a oublie d'appeler la methode initialize() !!!
  assert(channel_data_.get_delta_z().size() == field.nk());

  input_velocity_x_ = &vx;
  input_velocity_y_ = &vy;
  input_velocity_z_ = &vz;
  input_indicatrice_ = &ref_ijk;
  input_field_ = &field;
  stored_curv_fram_layer_z_ = -1000; // put a non-existant layer index: curv_fram will be computed at first call
  // Storage for curvature and fram limiter. We need 1 ghost layer:
  //  flux at left of the leftmost field data requires curv and fram on the element at left
  //  flux at the right of the rightmost field data requires on the element at right
  // We need 4 layers of temporary storage: curv and fram, and, for each, two consecutive layers in z
  const int ni = field.ni();
  const int nj = field.nj();
  tmp_curv_fram_.allocate(ni, nj, 4, 1);
  compute_set(result);
  input_field_ = 0;
  input_velocity_x_ = 0;
  input_velocity_y_ = 0;
  input_velocity_z_ = 0;
  input_indicatrice_ = 0;

}

void OpConvDiscIJKQuickScalar_double::ajouter(const IJK_Field_double& field,
                                              const IJK_Field_double& ref_ijk,
                                              const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                              IJK_Field_double& result)
{
  // Si ce test plante, c'est qu'on a oublie d'appeler la methode initialize() !!!
  assert(channel_data_.get_delta_z().size() == field.nk());

  input_velocity_x_ = &vx;
  input_velocity_y_ = &vy;
  input_velocity_z_ = &vz;
  input_field_ = &field;
  input_indicatrice_ = &ref_ijk;
  stored_curv_fram_layer_z_ = -1000; // put a non-existant layer index: curv_fram will be computed at first call
  // Storage for curvature and fram limiter. We need 1 ghost layer:
  //  flux at left of the leftmost field data requires curv and fram on the element at left
  //  flux at the right of the rightmost field data requires on the element at right
  // We need 4 layers of temporary storage: curv and fram, and, for each, two consecutive layers in z
  const int ni = field.ni();
  const int nj = field.nj();
  tmp_curv_fram_.allocate(ni, nj, 4, 1);
  compute_add(result);
  input_field_ = 0;
  input_velocity_x_ = 0;
  input_velocity_y_ = 0;
  input_velocity_z_ = 0;
  input_indicatrice_ = 0;

}
