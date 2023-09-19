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

#include <Operateur_IJK_elem_conv_base.h>

Implemente_base_sans_constructeur(Operateur_IJK_elem_conv_base_double,"Operateur_IJK_elem_conv_base_double",Operateur_IJK_elem_base_double);

Operateur_IJK_elem_conv_base_double::Operateur_IJK_elem_conv_base_double()
{

  stored_curv_fram_layer_z_ = -1000;
  input_field_ = nullptr;
  input_velocity_x_ = nullptr;
  input_velocity_y_ = nullptr;
  input_velocity_z_ = nullptr;
  perio_k_ = false;
  stored_curv_fram_layer_z_ = 0; // which (local) layer is currently stored in layer 0 of the tmp array ?

  corrige_flux_ = nullptr;
  indicatrice_ = nullptr;

  is_corrected_ = false;
  is_grad_ = false;
}

Sortie& Operateur_IJK_elem_conv_base_double::printOn(Sortie& os) const
{
  return os;
}

Entree& Operateur_IJK_elem_conv_base_double::readOn(Entree& is)
{
  return is;
}

void Operateur_IJK_elem_conv_base_double::initialize(const IJK_Splitting& splitting)
{
  perio_k_= splitting.get_grid_geometry().get_periodic_flag(DIRECTION_K);
  channel_data_.initialize(splitting);
  input_field_ = nullptr;
  input_velocity_x_ = nullptr;
  input_velocity_y_ = nullptr;
  input_velocity_z_ = nullptr;
}

void Operateur_IJK_elem_conv_base_double::calculer(const IJK_Field_double& field,
                                                   const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                   IJK_Field_double& result)
{
  // Si ce test plante, c'est qu'on a oublie d'appeler la methode initialize() !!!
  assert(channel_data_.get_delta_z().size() == field.nk());

  input_velocity_x_ = &vx;
  input_velocity_y_ = &vy;
  input_velocity_z_ = &vz;
  input_field_ = &field;
  stored_curv_fram_layer_z_ = -1000; // put a non-existant layer index: curv_fram will be computed at first call
  // Storage for curvature and fram limiter. We need 1 ghost layer:
  // flux at left of the leftmost field data requires curv and fram on the element at left
  // flux at the right of the rightmost field data requires on the element at right
  // We need 4 layers of temporary storage: curv and fram, and, for each, two consecutive layers in z
  const int ni = field.ni();
  const int nj = field.nj();
  tmp_curv_fram_.allocate(ni, nj, 4, 1);
  compute_set(result);
  input_field_ = nullptr;
  input_velocity_x_ = nullptr;
  input_velocity_y_ = nullptr;
  input_velocity_z_ = nullptr;

}

void Operateur_IJK_elem_conv_base_double::ajouter(const IJK_Field_double& field,
                                                  const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                  IJK_Field_double& result)
{
  // Si ce test plante, c'est qu'on a oublie d'appeler la methode initialize() !!!
  assert(channel_data_.get_delta_z().size() == field.nk());

  input_velocity_x_ = &vx;
  input_velocity_y_ = &vy;
  input_velocity_z_ = &vz;
  input_field_ = &field;
  stored_curv_fram_layer_z_ = -1000; // put a non-existant layer index: curv_fram will be computed at first call
  // Storage for curvature and fram limiter. We need 1 ghost layer:
  //  flux at left of the leftmost field data requires curv and fram on the element at left
  //  flux at the right of the rightmost field data requires on the element at right
  // We need 4 layers of temporary storage: curv and fram, and, for each, two consecutive layers in z
  const int ni = field.ni();
  const int nj = field.nj();
  tmp_curv_fram_.allocate(ni, nj, 4, 1);
  compute_add(result);
  input_field_ = nullptr;
  input_velocity_x_ = nullptr;
  input_velocity_y_ = nullptr;
  input_velocity_z_ = nullptr;

}


// Copy curv_fram values from layers 1 and 3 to layers 0 and 2
// (called at end of the computation of fluxes in z direction to keep the values
//  for the next layer).
void Operateur_IJK_elem_conv_base_double::shift_curv_fram(IJK_Field_local_double& tmp_curv_fram)
{
  const int ni = tmp_curv_fram.ni();
  const int nj = tmp_curv_fram.nj();
  int i, j;

  for (j = 0; j < nj; j++)
    for (i = 0; i < ni; i++)
      tmp_curv_fram(i,j,0) = tmp_curv_fram(i,j,1);

  for (j = 0; j < nj; j++)
    for (i = 0; i < ni; i++)
      tmp_curv_fram(i,j,2) = tmp_curv_fram(i,j,3);

}


// compute the "curv" value and part of the fram limiter (in the Fram routine from Eval_Quick_VDF_Elem.h) for 3 adjacent T values
// store result in tmp_curv_fram_, z layer=1 (curv) and z layer=3 (fram)
void Operateur_IJK_elem_conv_base_double::compute_curv_fram(DIRECTION _DIR_, int k_layer)
{
  int index = _DIR_ == DIRECTION::Y ? -1 : 0;
  ConstIJK_double_ptr input_field(*input_field_, 0, index, k_layer);
  // Where to store result:
  IJK_double_ptr curv_values(tmp_curv_fram_, 0, index, 1);
  IJK_double_ptr fram_values(tmp_curv_fram_, 0, index, 3);
  const int ni = tmp_curv_fram_.ni();
  const int nj = tmp_curv_fram_.nj();

  if(_DIR_ == DIRECTION::Z)
    {
      // Are we on the wall ?
      const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
      // global index of the layer of flux of the wall
      //  (assume one walls at zmin and zmax)
      const int first_global_k_layer = 0;
      const int last_global_k_layer = channel_data_.nb_elem_k_tot();

      if (!perio_k_ && (global_k_layer <= first_global_k_layer || global_k_layer >= last_global_k_layer))
        {
          for (int j = 0; j < nj; j++)
            for (int i = 0; i < ni; i++)
              tmp_curv_fram_(i, j, 1) = 0.;
          for (int j = 0; j < nj; j++)
            for (int i = 0; i < ni; i++)
              tmp_curv_fram_(i, j, 3) = 1.; // fram = 1 => upwind scheme on first layer of cells in z direction.

          return;
        }
    }

  double factor01 = -1.0;
  double factor12 = -1.0;
  if(_DIR_==DIRECTION::X)
    {
      // Compute invd1 and invd2 factors:
      const double inv_dx = 1./channel_data_.get_delta_x();
      factor01 = inv_dx * inv_dx;
      factor12 = factor01;
    }
  if(_DIR_==DIRECTION::Y)
    {
      // Compute invd1 and invd2 factors:
      const double inv_dy = 1. / channel_data_.get_delta_y();
      factor01 = inv_dy * inv_dy;
      factor12 = factor01;
    }
  if(_DIR_==DIRECTION::Z)
    {
      // Compute invd1 and invd2 factors:
      const double dz0 = channel_data_.get_delta_z()[k_layer-1];
      const double dz1 = channel_data_.get_delta_z()[k_layer];
      const double dz2 = channel_data_.get_delta_z()[k_layer+1];
      factor01 = 1. / (dz1 * (dz0 + dz1) * 0.5);
      factor12 = 1. / (dz1 * (dz1 + dz2) * 0.5);
    }

  const int vsize = Simd_double::size();
  const int imax = _DIR_==DIRECTION::X ? ni + 1 : ni;
  const int imax1 = imax - (vsize-1); // test to check for end of vectorizable part
  const int jmax =  _DIR_==DIRECTION::Y ? nj + 1 : nj;
  for (int j = 0; ; j++)
    {
      int i;
      for (i = -1; i < 0; i++)
        compute_curv_fram_loop_(_DIR_, i, factor12, factor01, input_field, curv_values, fram_values );

      for (i = 0; i < imax1; i += vsize)
        {
          Simd_double t0 = 0., t1 = 0., t2 = 0.;
          input_field.get_left_center_right(_DIR_, i, t0, t1, t2);
          Simd_double curv = (t2 - t1) * factor12 - (t1 - t0) * factor01;
          curv_values.put_val(i, curv);
          Simd_double smin = min(t0, t2);
          Simd_double smax = max(t0, t2);
          // Compared to original code (Eval_Quick_VDF_Elem.h), we first compute the 4th power,
          // then take the max (dabs is then useless)
          Simd_double zeroVec = 0.;
          Simd_double oneVec = 1.;
          Simd_double minVec = DMINFLOAT;
          Simd_double dsabs = SimdSelect(zeroVec, smax - smin, smax - smin, smin - smax);
          Simd_double ds = SimdSelect(dsabs, minVec, oneVec, smax - smin);
          Simd_double sr = SimdSelect(dsabs, minVec, zeroVec, ((t1 - smin) / ds - 0.5) * 2.);
          sr *= sr;
          sr *= sr;
          sr = min(sr, oneVec);

          fram_values.put_val(i, sr);
        }
      for (; i < imax; i++)
        compute_curv_fram_loop_(_DIR_, i, factor12, factor01, input_field, curv_values, fram_values );

      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      input_field.next_j();
      curv_values.next_j();
      fram_values.next_j();

    }
}

void Operateur_IJK_elem_conv_base_double::compute_curv_fram_loop_(DIRECTION _DIR_, int iter, double factor12, double factor01, const ConstIJK_double_ptr& input_field, IJK_double_ptr& curv_values, IJK_double_ptr& fram_values )
{
  double t0 = 0., t1 = 0., t2 = 0.;
  input_field.get_left_center_right(_DIR_, iter, t0, t1, t2);
  double curv = (t2 - t1) * factor12 - (t1 - t0) * factor01;
  curv_values.put_val(iter, curv);
  double smin = std::min(t0, t2);
  double smax = std::max(t0, t2);
  // Compared to original code (Eval_Quick_VDF_Elem.h), we first compute the 4th power,
  // then take the max (dabs is then useless)
  double sr = std::fabs(smax - smin)<DMINFLOAT ? 0. : ((t1 - smin) / (smax - smin) - 0.5) * 2.;
  sr *= sr;
  sr *= sr;
  sr = std::min(sr, 1.);

  fram_values.put_val(iter, sr);
}

