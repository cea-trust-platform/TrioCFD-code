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

#ifndef OpConvIJKQuickScalar_H_TPP
#define OpConvIJKQuickScalar_H_TPP

#include <IJK_Field.h>

template <DIRECTION _DIR_>

void OpConvQuickIJKScalar_double::compute_flux_(IJK_Field_local_double& resu, const int k_layer)
{
  if(_DIR_==DIRECTION::Z)
    {
      // The previous layer of curv and fram values might have been computed already
      if (stored_curv_fram_layer_z_ != k_layer-1)
        {
          compute_curv_fram(_DIR_, k_layer-1);
          shift_curv_fram(tmp_curv_fram_);
        }
    }
  compute_curv_fram(_DIR_, k_layer);

  ConstIJK_double_ptr velocity_dir(get_input_velocity(_DIR_), 0, 0, k_layer);
  ConstIJK_double_ptr input_field(*input_field_, 0, 0, k_layer);
  ConstIJK_double_ptr curv_values(tmp_curv_fram_, 0, 0, 1); /* if z direction, "left" will be in layer 0 */
  ConstIJK_double_ptr fram_values(tmp_curv_fram_, 0, 0, 3); /* if z direction, "left" is in layer 2 */
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);
  const int nx = _DIR_==DIRECTION::X ? input_field_->ni() + 1 : input_field_->ni();
  const int ny = _DIR_==DIRECTION::Y ? input_field_->nj() + 1 : input_field_->nj();
  const double dx = channel_data_.get_delta_x();
  const double surface = channel_data_.get_delta_y() * channel_data_.get_delta_z()[k_layer];
  if(_DIR_==DIRECTION::Z)
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

  const double dx_squared_over_8 = dx * dx * 0.125;
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
          Simd_double fram0, fram1;
          Simd_double curv0, curv1;
          input_field.get_left_center(_DIR_, i, T0, T1);
          fram_values.get_left_center(_DIR_, i, fram0, fram1);
          curv_values.get_left_center(_DIR_, i, curv0, curv1);
          Simd_double fram = max(fram0, fram1);
          Simd_double curv	   = select_double(velocity, 0., curv1, curv0);
          Simd_double T_amont = select_double(velocity, 0., T1 /* if velocity < 0 */, T0 /* if velocity > 0 */);
          Simd_double flux	   = (T0 + T1) * 0.5 - dx_squared_over_8 * curv;
          flux		   = ((1. - fram) * flux + fram * T_amont) * velocity * surface;

          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      input_field.next_j();
      velocity_dir.next_j();
      fram_values.next_j();
      curv_values.next_j();
      resu_ptr.next_j();
    }
}
#endif
