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

#ifndef OpConvIJKAmont_H_TPP
#define OpConvIJKAmont_H_TPP

#include <IJK_Field_simd_tools.h>

/*! @brief compute fluxes in direction _DIR_ for velocity component _VCOMPO_ for the layer of fluxes k_layer
 *
 *  4-th order centered convection scheme
 *
 */
template <DIRECTION _DIR_, DIRECTION _VCOMPO_>
void OpConvAmontIJK_double::compute_flux_(IJK_Field_local_double& resu, const int k_layer)
{
  // convected field
  const IJK_Field_local_double& src = get_input(_VCOMPO_);
  // Convected vector field:
  ConstIJK_double_ptr src_ptr(src, 0, 0, k_layer);
  // Velocity in x direction (convecting velocity)
  ConstIJK_double_ptr vconv_ptr(get_v(_DIR_), 0, 0, k_layer);

  const int idir = (int)_DIR_;
  const int nx = _DIR_ == DIRECTION::X ? src.ni() + 1 : src.ni();
  const int ny = _DIR_ == DIRECTION::Y ? src.nj() + 1 : src.nj();
  const int icompo = (int)_VCOMPO_;

  // Result (fluxes in direction x for component x of the convected field)
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);

  const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
  // global index of the layer of flux of the wall
  //  (assume one walls at zmin and zmax)
  const int first_global_k_layer = channel_data_.first_global_k_layer_flux(icompo, idir);
  const int last_global_k_layer = channel_data_.last_global_k_layer_flux(icompo, idir);

  // GB 21/12/2020 : Similarly to velocity,
  // if (global_k_layer < first_global_k_layer || global_k_layer > last_global_k_layer) {
  if (!perio_k_ && (global_k_layer <= first_global_k_layer || global_k_layer >= last_global_k_layer))
    {
      // We are in the wall
      putzero(resu);
      return;
    }

  double half_surface = channel_data_.get_surface(k_layer, icompo, idir) * 0.5;
  {
    const int imax = nx;
    const int jmax = ny;
    const int vsize = Simd_double::size();
    for (int j = 0; ; j++)
      {
        for (int i = 0; i < imax; i += vsize)     // specific coding for uniform mesh in x and y: surfaces are constant on an xy plane
          {
            Simd_double vit_0,vit_1; // 2 adjacent velocity values
            src_ptr.get_left_center(_DIR_,i,vit_0,vit_1);
            // get convecting velocity
            Simd_double vconv0, vconv1;
            vconv_ptr.get_left_center(_VCOMPO_, i, vconv0, vconv1);

            // Average of the convecting velocity (copied from Eval_Amont_VDF_Face : not weighted)
            Simd_double psc = (vconv0 + vconv1) * half_surface;
            Simd_double upwind_velocity = select_double(psc, 0., vit_1 /* if psc < 0 */, vit_0 /* if psc > 0 */);
            Simd_double flux = psc * upwind_velocity;
            resu_ptr.put_val(i, flux);
          }
        // do not execute end_iloop at last iteration (because of assert on valid j+1)
        if (j+1==jmax)
          break;
        // instructions to perform to jump to next row
        src_ptr.next_j();
        resu_ptr.next_j();
        vconv_ptr.next_j();
      }
  }
}

#endif
