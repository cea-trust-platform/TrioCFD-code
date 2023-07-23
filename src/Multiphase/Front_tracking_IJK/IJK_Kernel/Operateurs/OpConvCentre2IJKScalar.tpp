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

#ifndef OpConvCentre2IJKScalar_TPP_included
#define OpConvCentre2IJKScalar_TPP_included
#include <IJK_Field.h>

template <DIRECTION _DIR_>
void OpConvCentre2IJKScalar_double::compute_flux_(IJK_Field_local_double& resu, const int k_layer)
{
  ConstIJK_double_ptr velocity_dir(get_input_velocity(_DIR_), 0, 0, k_layer);
  ConstIJK_double_ptr input_field(*input_field_, 0, 0, k_layer);
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);
  const int nx = _DIR_==DIRECTION::X ? input_field_->ni() + 1 : input_field_->ni();
  const int ny = _DIR_==DIRECTION::Y ? input_field_->nj() + 1 : input_field_->nj();
  const double surface = channel_data_.get_delta_y() * channel_data_.get_delta_z()[k_layer];

  if(_DIR_ == DIRECTION::Z)
    {
      // Are we on the wall ?
      const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
      // global index of the layer of flux of the wall
      //  (assume one walls at zmin and zmax)
      const int first_global_k_layer = 0;
      const int last_global_k_layer = channel_data_.nb_elem_k_tot();

      // GB 21/12/2020 : Similarly to velocity, we make the same adjustments.
      // if (global_k_layer == first_global_k_layer || global_k_layer == last_global_k_layer) {
      // ie (i) replace the former "==" by "<=" and ">=" to be identic to the condition in OpCentre4IJK
      // and (ii) add the condition on perio_k
      if (!perio_k_ && (global_k_layer <= first_global_k_layer || global_k_layer >= last_global_k_layer))
        {
          // We are on (or worse inside) the wall, zero flux
          putzero(resu);
          return;
        }
    }

  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double velocity;
          velocity_dir.get_center(i, velocity);
          Simd_double T0, T1; // scalar value at left and at right of the computed flux
          input_field.get_left_center(_DIR_,i, T0, T1);
          Simd_double flux = (T0 + T1) * 0.5;
          flux = flux * velocity * surface;
          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      input_field.next_j();
      velocity_dir.next_j();
      resu_ptr.next_j();
    }

  if(_DIR_ == DIRECTION::Z)
    {
      // store curv and fram for next layer of fluxes in z direction
      shift_curv_fram(tmp_curv_fram_);
      stored_curv_fram_layer_z_ = k_layer;
    }
}

#endif


