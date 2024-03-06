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

#ifndef OpConvQuickSharpIJK_TPP_included
#define OpConvQuickSharpIJK_TPP_included

inline double sharp(const double utc)
{
  double cf;
  if (utc > 0)
    {
      if (utc > 0.25)
        {
          if (utc <= 1.)
            {
              cf = 0.25 * (1. - utc);
            }
          else if (utc < 1.5)
            {
              cf = 0.25 * (utc - 1.);
            }
          else
            {
              cf = 0.125;
            }
        }
      else
        {
          cf = 0.5 - 0.625*sqrt(utc);
        }
    }
  else
    {
      if (utc > -1)
        {
          cf = 0.5 + 0.375*utc;
        }
      else
        {
          cf = 0.125;
        }
    }
  return cf;
}

inline double conv_quick_sharp_plus(const double psc,const double vit_0, const double vit_1,
                                    const double vit_0_0, const double dx,
                                    const double dm, const double dxam)
{
  double cf;
  double curv;
  double delta_0 = vit_0 - vit_0_0;
  double delta = vit_1 - vit_0;
  double dd1,utc;
  double delta_delta;

  curv = (delta/dx - delta_0/dxam)/dm ;

  // Calcul de cf:

  delta_delta = delta_0+delta;
  dd1 = std::fabs(delta_delta);
  if (dd1 < 1.e-5)
    cf = 0.125;
  else
    {
      utc = delta_0/delta_delta;
      cf = sharp(utc);
    }

  return (0.5*(vit_0 + vit_1) - cf*(dx*dx)*curv)*psc;

}

inline double conv_quick_sharp_moins(const double psc,const double vit_0,const double vit_1,
                                     const double vit_1_1,const double dx,
                                     const double dm, const double dxam)
{
  double cf;
  double curv;
  double delta_1 = vit_1_1 - vit_1;
  double delta = vit_1 - vit_0;
  double dd1,utc;
  double delta_delta;

  curv = ( delta_1/dxam - delta/dx )/dm ;

  // Calcul de cf:

  delta_delta = delta_1+delta;
  dd1 = std::fabs(delta_delta);
  if (dd1 < 1.e-5)
    cf = 0.125;
  else
    {
      utc = delta_1/delta_delta;
      cf = sharp(utc);
    }

  return (0.5*(vit_0 + vit_1) - cf*(dx*dx)*curv)*psc;

}

/*! @brief compute fluxes in direction x for velocity component x for the layer of fluxes k_layer
 *
 *  4-th order centered convection scheme
 *
 */
template <DIRECTION _DIR_, DIRECTION _VCOMPO_>
void OpConvQuickSharpIJK_double::compute_flux_(IJK_Field_local_double& resu, const int k_layer)
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

  // GB 23/04/2020. Before, it was not working for tri-periodic cases :
  // if (global_k_layer < first_global_k_layer || global_k_layer > last_global_k_layer) {
  // BUT ALSO, it was strict inf or sup to.
  // It was replaced by "<=" and ">=" to be identic to the condition in OpCentre4IJK
  // I'm not sure which one is correct, or if it makes a difference on a wall
  // (maybe the convection is zero anyway at the wall, and besides, the contribution
  //  of convection at the wal is set to zero anyway)?
  if (!perio_k_ && (global_k_layer <= first_global_k_layer || global_k_layer >= last_global_k_layer))
    {
      // We are in the wall
      putzero(resu);
      return;
    }

  const double half_surface = channel_data_.get_surface(k_layer, icompo, idir) * 0.5;
  const double delta = get_delta(_DIR_);
  // Non vectorized coding
  for (int j = 0; j < ny; j++)
    {
      for (int i = 0; i < nx; i++)
        {
          // get convecting velocity
          double vconv0, vconv1;
          vconv_ptr.get_left_center(_VCOMPO_,i, vconv0, vconv1);
          double psc = (vconv0 + vconv1) * half_surface;
          double vit_0_0,vit_0,vit_1,vit_1_1; // 4 adjacent velocity values
          src_ptr.get_leftleft_left_center_right(_DIR_,i,vit_0_0,vit_0,vit_1,vit_1_1);
          double flux;
          if (psc > 0)
            {
              flux = conv_quick_sharp_plus(psc, vit_0, vit_1, vit_0_0, delta, delta, delta);
            }
          else
            {
              flux = conv_quick_sharp_moins(psc, vit_0, vit_1, vit_1_1, delta, delta, delta);
            }
          resu_ptr.put_val(i, flux);
        }
      if (j < ny-1)
        {
          src_ptr.next_j();
          resu_ptr.next_j();
          vconv_ptr.next_j();
        }
    }
}

#endif

