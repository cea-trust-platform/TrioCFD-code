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

#ifndef OpDiffTurbIJKScalar_H_TPP
#define OpDiffTurbIJKScalar_H_TPP

static void copy_boundary_condition(const IJK_Field_local_double& boundary_flux, IJK_Field_local_double& resu)
{
  assert(resu.ni() >= boundary_flux.ni());
  assert(resu.nj() >= boundary_flux.nj());
  // resu is the temporary array where all fluxes are stored before computing divergence,
  // they might have more place than ni and nj because there 1 more flux value dans velocity values
  // to compute divergence
  const int ni = boundary_flux.ni();
  const int nj = boundary_flux.nj();
  for (int j = 0; j < nj; j++)
    for (int i = 0; i < ni; i++)
      resu(i,j,0) = boundary_flux(i,j,0);
}

template <DIRECTION _DIR_>
void OpDiffIJKScalarGeneric_double::compute_flux_(IJK_Field_local_double& resu, const int k_layer)
{
  const int nx = _DIR_ == DIRECTION::X ? input_field_->ni() + 1 : input_field_->ni() ;
  const int ny = _DIR_ == DIRECTION::Y ? input_field_->nj() + 1 : input_field_->nj();

  ConstIJK_double_ptr input_field(*input_field_, 0, 0, k_layer);
//  const IJK_Field_local_double dummy_field = *input_field_;
  Simd_double uniform_lambda(1.);
  Simd_double avg_lambda;
  if (is_uniform_ and uniform_lambda_!=0)
    {
      Simd_double uniform_lambda_tmp(*uniform_lambda_);
      uniform_lambda *= uniform_lambda_tmp;
    }
  if (is_uniform_)
    {
      lambda_=input_field_;
    }
  ConstIJK_double_ptr lambda(is_vectorial_? get_model(_DIR_) : *lambda_, 0, 0, k_layer);

  ConstIJK_double_ptr structural_model(is_structural_ ? get_model(_DIR_) : *lambda_, 0, 0, k_layer);

  IJK_double_ptr resu_ptr(resu, 0, 0, 0);

  if(_DIR_ == DIRECTION::Z)
    {
      // Are we on the wall ?
      const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
      // global index of the layer of flux of the wall
      //  (assume one walls at zmin and zmax)
      const int first_global_k_layer = 0; // index of k_layer when we are on the wall
      // Fluxes in direction k are on the faces, hence, index of last flux is equal to number of elements
      const int last_global_k_layer =  channel_data_.nb_elem_k_tot();

      if (!perio_k_)
        {
          if (global_k_layer == first_global_k_layer)
            {
              // We are on wall at kmin, copy boundary condition fluxes to "resu"
              if (boundary_flux_kmin_) // boundary condition is not zero flux
                copy_boundary_condition(*boundary_flux_kmin_, resu);
              else
                putzero(resu);
              return;
            }
          else if (global_k_layer == last_global_k_layer)
            {
              if (boundary_flux_kmax_) // boundary condition is not zero flux
                copy_boundary_condition(*boundary_flux_kmax_, resu);
              else
                putzero(resu);
              return;
            }
        }
    }

  double d0 = 0., d1 = 0.;
  double surface = 0.;
  if(_DIR_ == DIRECTION::X)
    {
      d0 = channel_data_.get_delta_x() * 0.5;
      d1 = d0;
      surface = channel_data_.get_delta_y() * channel_data_.get_delta_z()[k_layer];
    }
  if(_DIR_ == DIRECTION::Y)
    {
      d0 = channel_data_.get_delta_y() * 0.5;
      d1 = d0;
      surface = channel_data_.get_delta_x() * channel_data_.get_delta_z()[k_layer];
    }
  if(_DIR_ == DIRECTION::Z)
    {
      d0 = channel_data_.get_delta_z()[k_layer-1] * 0.5;
      d1 = channel_data_.get_delta_z()[k_layer] * 0.5;
      surface = channel_data_.get_delta_x() * channel_data_.get_delta_y();
    }

  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double flux = 0.;
          if (is_structural_)
            {
              Simd_double s_mo, s_mo_dummy;
              structural_model.get_left_center(_DIR_, i, s_mo_dummy, s_mo);
              flux = (-1.) * s_mo * surface;
            }
          else
            {
              Simd_double left_val, right_val;
              input_field.get_left_center(_DIR_, i, left_val, right_val);
              Simd_double zeroVec = 0.;
              Simd_double oneVec = 1.;
              Simd_double minVec = DMINFLOAT;
              Simd_double d = 1.;
              Simd_double lambda_m1(uniform_lambda);
              Simd_double lambda_m2(uniform_lambda);
              // Fetch conductivity on neighbour cells:
              if (!is_uniform_)
                {
                  lambda.get_left_center(_DIR_, i, lambda_m1, lambda_m2);
                }
              // geometric avg: (d0+d1) / ( d0 / lambda_m1 + d1 / lambda_m2 ), optimized with only 1 division:
              Simd_double dsabs = SimdSelect(zeroVec, d0 * lambda_m2 + d1 * lambda_m1, d0 * lambda_m2 + d1 * lambda_m1, (-1) * (d0 * lambda_m2 + d1 * lambda_m1));
              Simd_double ds = SimdSelect(dsabs, minVec, oneVec, d0 * lambda_m2 + d1 * lambda_m1);
              if(is_anisotropic_)
                d = d0 + d1;
              avg_lambda = SimdSelect(dsabs, minVec, zeroVec, SimdDivideMed(d * lambda_m1 * lambda_m2, ds));
//                  // thermal flux is positive if going from left to right => -grad(T)
              flux = (left_val - right_val) * avg_lambda * surface;
//              flux = (left_val - right_val) * 0.1 * surface;
//              if (flux.data_ != 0.)
//                Cerr << "Diffusive flux : " << i << j << k_layer << "=" << flux.data_ << finl;
            }
          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      input_field.next_j();
      if (!is_uniform_)
        { lambda.next_j(); }
      if(is_structural_)
        structural_model.next_j();
      resu_ptr.next_j();
    }

}

#endif

