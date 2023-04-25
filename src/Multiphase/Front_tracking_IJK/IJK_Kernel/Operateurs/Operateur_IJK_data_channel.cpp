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

#include <Operateur_IJK_data_channel.h>

void Operateur_IJK_data_channel::initialize(const IJK_Splitting& splitting)
{
  const IJK_Grid_Geometry& grid_geom = splitting.get_grid_geometry();

  nb_elem_k_tot_ = grid_geom.get_nb_elem_tot(2);
  offset_to_global_k_layer_ = splitting.get_offset_local(2);
  delta_x_ = grid_geom.get_constant_delta(0);
  delta_y_ = grid_geom.get_constant_delta(1);

  splitting.get_local_mesh_delta(2 /*direction*/, 2 /* number of ghost cells to get */, delta_z_);

  inv_delta_x_ = 1. / delta_x_;
  inv_delta_y_ = 1. / delta_y_;
  {
    const int n = delta_z_.size();
    inv_elem_size_z_.resize(n, 1);
    for (int i = -1; i < n+1; i++)
      inv_elem_size_z_[i] = 1. / delta_z_[i];
    inv_dist_z_elemcenter_.resize(n+1,0); // neither
    for (int i = 0; i < n+1; i++)
      {
        const int global_k = i + offset_to_global_k_layer_;
        if (global_k == first_global_k_layer_flux(0, 2) - 1)
          {
            // bottom wall
            inv_dist_z_elemcenter_[i] = 2. / delta_z_[i];
          }
        else if (global_k == last_global_k_layer_flux(0, 2) + 1)
          {
            // top wall
            inv_dist_z_elemcenter_[i] = 2. / delta_z_[i-1];
          }
        else
          {
            inv_dist_z_elemcenter_[i] = 2. / (delta_z_[i-1] + delta_z_[i]);
          }
      }
  }
}
