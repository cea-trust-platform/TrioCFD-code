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

#include <Operateur_IJK_base.h>
// Definition of the fluxes
//
//  control volumes for velocity in x direction:
//  this face_i(i,j,k) is located at left of element(i,j,k)
//
//     -------------------------
//     .           |           .
//     .           |           .
//     .           |           .
//     .           |           .
//     .           |----->     .
//     .           |           .
//     .           |           .
//     .           |           .
//     -------------------------
//  face(i,j,k) receives fluxes from all neighbour faces:
//   flux_i(i,j,k) at left
//   flux_j(i,j,k)
//   flux_k(i,j,k)
//   flux_i(i+1,j,k) at right
//   ...

static inline void swap_data(IJK_Field_local_double*& a, IJK_Field_local_double*& b)
{
  IJK_Field_local_double *tmp=a;
  a=b;
  b=tmp;
}

// Compute the "- divergence" of the given flux, integrated over control volume
// We assume than flux contain the integral of the flux on the faces of the control volumes
// For "div_set" method:
//  resu(i,j,k) = flux_x(i,j) - flux_x(i+1,j) + flux_y(i,j) - flux_y(i,j+1) + flux_zmin(i,j) - flux_zmax(i,j)
// For "div_add" method: same but
//  resu(i,j,k) += ...
static void Operator_IJK_div(const IJK_Field_local_double& flux_x, const IJK_Field_local_double& flux_y,
                             const IJK_Field_local_double& flux_zmin, const IJK_Field_local_double& flux_zmax,
                             IJK_Field_local_double& resu, int k_layer, bool add)
{
  ConstIJK_double_ptr fx(flux_x, 0, 0, 0);
  ConstIJK_double_ptr fy(flux_y, 0, 0, 0);
  ConstIJK_double_ptr fzmin(flux_zmin, 0, 0, 0);
  ConstIJK_double_ptr fzmax(flux_zmax, 0, 0, 0);
  IJK_double_ptr resu_ptr(resu, 0, 0, k_layer);


  const int imax = resu.ni();
  const int jmax = resu.nj();
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double r, a, b;
          fx.get_center_right(DIRECTION::X, i, a, b);
          r = a-b;

          fy.get_center_right(DIRECTION::Y,i, a, b);
          r += a-b;

          fzmin.get_center(i, a);
          fzmax.get_center(i, b);
          r += a-b;

          if(add)
            {
              resu_ptr.get_center(i, a);
              r += a;
            }
          resu_ptr.put_val(i, r);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      fx.next_j();
      fy.next_j();
      fzmin.next_j();
      fzmax.next_j();
      resu_ptr.next_j();
    }

}


/*
 * Operateur_IJK_faces_base
 */

Implemente_base(Operateur_IJK_faces_base_double,"Operateur_IJK_faces_base_double", Objet_U);

Sortie& Operateur_IJK_faces_base_double::printOn(Sortie& os) const
{
  //  Objet_U::printOn(os);
  return os;
}

Entree& Operateur_IJK_faces_base_double::readOn(Entree& is)
{
  //  Objet_U::readOn(is);
  return is;
}

// Compute the maximum stable timestep for a convection scheme on the local processor.
// The global stability timestep is the minimum of all values on all processors.
// The ghost velocities must be uptodate on the right each subdomain.
double Operateur_IJK_faces_base_double::compute_dtstab_convection_local(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz)
{
  double inverse_dt_conv = 0.;
  const IJK_Splitting& split = vx.get_splitting();
  const int ni = split.get_nb_elem_local(DIRECTION_I);
  const int nj = split.get_nb_elem_local(DIRECTION_J);
  const int nk = split.get_nb_elem_local(DIRECTION_K);
  const IJK_Grid_Geometry& geom = split.get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const ArrOfDouble& delta_z = geom.get_delta(DIRECTION_K);
  const int k_offset = split.get_offset_local(DIRECTION_K);

  for (int k = 0 ; k < nk ; k++)
    {
      // pointeur pour optimisation du calcul des voisins ( marche pour les derivees)
      ConstIJK_double_ptr ptr_vx(vx, 0, 0, k);
      ConstIJK_double_ptr ptr_vy(vy, 0, 0, k);
      ConstIJK_double_ptr ptr_vz(vz, 0, 0, k);
      const double dz = delta_z[k + k_offset];
      const double inverse_dx = 1. / dx;
      const double inverse_dy = 1. / dy;
      const double inverse_dz = 1. / dz;

      for (int j = 0 ; j < nj ; j++)
        {
          for (int i = 0 ; i < ni ; i++)
            {
              // pas de temps de convection
              double vx0 , vx1 , vy0 , vy1 , vz0 , vz1 ; // vitesses face gauche en 0 et face droite en 1
              ptr_vx.get_center_right(DIRECTION::X, i, vx0, vx1);
              ptr_vy.get_center_right(DIRECTION::Y, i, vy0, vy1);
              ptr_vz.get_center_right(DIRECTION::Z, i, vz0, vz1);
              double inverse_dt_conv_loc =
                (std::max(vx0, 0.) - std::min(vx1, 0.)) * inverse_dx
                + (std::max(vy0, 0.) - std::min(vy1, 0.)) * inverse_dy
                + (std::max(vz0, 0.) - std::min(vz1, 0.)) * inverse_dz;

              inverse_dt_conv = std::max(inverse_dt_conv, inverse_dt_conv_loc);
            }
          if ( j < nj-1 )
            {
              ptr_vx.next_j();
              ptr_vy.next_j();
              ptr_vz.next_j();
            }
        }
    }
  inverse_dt_conv = std::max(1e-20, inverse_dt_conv);
  return 1. / inverse_dt_conv;
}

void Operateur_IJK_faces_base_double::compute_set(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  compute_(dvx, dvy, dvz, 0 /* setting */);
}

void Operateur_IJK_faces_base_double::compute_add(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  compute_(dvx, dvy, dvz, 1 /* adding */);
}

void Operateur_IJK_faces_base_double::compute_(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz, bool add)
{
  IJK_Field_local_double storage;
  storage.allocate(dvx.ni()+1, dvy.nj()+1, 6, 0);
  IJK_Field_local_double tmp[6];
  for (int i = 0; i < 6; i++)
    tmp[i].ref_ij(storage, i);
  IJK_Field_local_double *flux_vx_zmin = &tmp[0]; // flux pour la composante x a travers la face zmin des volumes de controles
  IJK_Field_local_double *flux_vy_zmin = &tmp[1]; // flux pour la composante y
  IJK_Field_local_double *flux_vz_zmin = &tmp[2]; // flux pour la composante z
  IJK_Field_local_double *flux_zmax = &tmp[3];    // flux a travers la face zmax des volumes de controles
  IJK_Field_local_double * const flux_x = &tmp[4]; // pointer will not change, content will change
  IJK_Field_local_double * const flux_y = &tmp[5];

  compute_flux_z_vx(*flux_vx_zmin, 0);
  compute_flux_z_vy(*flux_vy_zmin, 0);
  compute_flux_z_vz(*flux_vz_zmin, 0);

//  correct_flux(flux_zmin, -1, 0);
//  correct_flux(flux_zmin, -1, 1);
//  correct_flux(flux_zmin, -1, 2);

  // At the bottom end of the domain, there is more layer of unknowns for vz:
  const int kmax = std::max(std::max(dvx.nk(), dvy.nk()), dvz.nk());
  for (int k = 0; k < kmax ; k++)
    {
      if (k < dvz.nk())
        {

          compute_flux_x_vz(*flux_x, k);
          compute_flux_y_vz(*flux_y, k);
          compute_flux_z_vz(*flux_zmax, k+1);

//          correct_flux(flux_x, k, 2);
//          correct_flux(flux_y, k, 2);
//          correct_flux(flux_zmax, k, 2);

          Operator_IJK_div(*flux_x, *flux_y, *flux_vz_zmin, *flux_zmax, dvz, k, add);

          exec_after_divergence_flux_z(dvz, k);
          swap_data(flux_vz_zmin, flux_zmax); // conserve le flux en z pour la couche z suivante
        }
      if (k < dvx.nk())
        {
          compute_flux_x_vx(*flux_x, k);
          compute_flux_y_vx(*flux_y, k);
          compute_flux_z_vx(*flux_zmax, k+1);

//          correct_flux(flux_x, k, 0);
//          correct_flux(flux_y, k, 0);
//          correct_flux(flux_zmax, k, 0);

          Operator_IJK_div(*flux_x, *flux_y, *flux_vx_zmin, *flux_zmax, dvx, k, add);

          exec_after_divergence_flux_x(dvx, k);
          swap_data(flux_vx_zmin, flux_zmax); // conserve le flux en z pour la couche z suivante
        }
      if (k < dvy.nk())
        {
          compute_flux_x_vy(*flux_x, k);
          compute_flux_y_vy(*flux_y, k);
          compute_flux_z_vy(*flux_zmax, k+1);

//          correct_flux(flux_x, k, 1);
//          correct_flux(flux_y, k, 1);
//          correct_flux(flux_zmax, k, 1);

          Operator_IJK_div(*flux_x, *flux_y, *flux_vy_zmin, *flux_zmax, dvy, k, add);

          exec_after_divergence_flux_y(dvy, k);
          swap_data(flux_vy_zmin, flux_zmax); // conserve le flux en z pour la couche z suivante
        }
    }
}


/*
 * Operateur_IJK_elem_base
 */

Implemente_base(Operateur_IJK_elem_base_double,"Operateur_IJK_elem_base_double",Objet_U);

Sortie& Operateur_IJK_elem_base_double::printOn(Sortie& os) const
{
  //  Objet_U::printOn(os);
  return os;
}

Entree& Operateur_IJK_elem_base_double::readOn(Entree& is)
{
  //  Objet_U::readOn(is);
  return is;
}

void Operateur_IJK_elem_base_double::compute_set(IJK_Field_double& dx)
{
  compute_(dx, 0);
}

void Operateur_IJK_elem_base_double::compute_add(IJK_Field_double& dx)
{
  compute_(dx, 1);
}


void Operateur_IJK_elem_base_double::compute_(IJK_Field_double& dx, bool add)
{
  IJK_Field_local_double storage;
  storage.allocate(dx.ni()+1, dx.nj()+1, 4, 0);
  IJK_Field_local_double tmp[4];
  for (int i = 0; i < 4; i++)
    tmp[i].ref_ij(storage, i);
  IJK_Field_local_double *flux_zmin = &tmp[0];
  IJK_Field_local_double *flux_zmax = &tmp[1];    // flux a travers la face zmax des volumes de controles
  IJK_Field_local_double *const flux_x = &tmp[2]; // pointer will not change, content will change
  IJK_Field_local_double *const flux_y = &tmp[3];

  compute_flux_z(*flux_zmin, 0);
  correct_flux(flux_zmin, -1, 2);
  /*
   * TODO: correct_flux_z_min();
   * fluxes_to_correct ?
   * Critere fabs(fluxes_to_correct)<DMAX_FLOAT ?
   * Zero valeur admissible
   */
  const int kmax = dx.nk();
  for (int k = 0; k < kmax; k++)
    {
      compute_flux_x(*flux_x, k);
      compute_flux_y(*flux_y, k);
      compute_flux_z(*flux_zmax, k+1);
      /*
       * Correct the fluxes before it is summed
       * with Operator_IJK_div
       */
      correct_flux(flux_x, k, 0);
      correct_flux(flux_y, k, 1);
      correct_flux(flux_zmin, k, 2);
      correct_flux(flux_zmax, k, 2);
      Operator_IJK_div(*flux_x, *flux_y, *flux_zmin, *flux_zmax, dx, k, add);
      swap_data(flux_zmin, flux_zmax); // conserve le flux en z pour la couche z suivante
    }
}

void Operateur_IJK_elem_base_double::compute_grad(FixedVector<IJK_Field_double, 3>& dx)
{
  IJK_Field_local_double storage;
  storage.allocate(dx[0].ni()+1, dx[0].nj()+1, 4, 0);
  IJK_Field_local_double tmp[4];
  for (int i = 0; i < 4; i++)
    tmp[i].ref_ij(storage, i);
  IJK_Field_local_double *flux_zmin = &tmp[0];
  IJK_Field_local_double *flux_zmax = &tmp[1];    // flux a travers la face zmax des volumes de controles
  IJK_Field_local_double *const flux_x = &tmp[2]; // pointer will not change, content will change
  IJK_Field_local_double *const flux_y = &tmp[3];

  compute_flux_z(*flux_zmin, 0);

  const int kmax = dx[0].nk();
  for (int k = 0; k < kmax; k++)
    {
      compute_flux_x(*flux_x, k);
      compute_flux_y(*flux_y, k);
      compute_flux_z(*flux_zmax, k+1);
      fill_grad_field_x_(*flux_x, dx, k);
      fill_grad_field_y_(*flux_y, dx, k);
      fill_grad_field_z_(*flux_zmin, *flux_zmax, dx, k);
      swap_data(flux_zmin, flux_zmax); // conserve le flux en z pour la couche z suivante
    }
}

void Operateur_IJK_elem_base_double::compute_grad_x(IJK_Field_double& dx)
{
  IJK_Field_local_double storage;
  storage.allocate(dx.ni()+1, dx.nj()+1, 1, 0);
  IJK_Field_local_double tmp;
  tmp.ref_ij(storage, 0);
  IJK_Field_local_double *const flux_x = &tmp;
  const int kmax = dx.nk();
  for (int k = 0; k < kmax; k++)
    {
      compute_flux_x(*flux_x, k);
      fill_grad_field_x_y_(*flux_x, dx, k, 0);
    }
}

void Operateur_IJK_elem_base_double::compute_grad_y(IJK_Field_double& dx)
{
  IJK_Field_local_double storage;
  storage.allocate(dx.ni()+1, dx.nj()+1, 1, 0);
  IJK_Field_local_double tmp;
  tmp.ref_ij(storage, 0);
  IJK_Field_local_double *const flux_y = &tmp;
  const int kmax = dx.nk();
  for (int k = 0; k < kmax; k++)
    {
      compute_flux_y(*flux_y, k);
      fill_grad_field_x_y_(*flux_y, dx, k, 1);
    }
}

void Operateur_IJK_elem_base_double::compute_grad_z(IJK_Field_double& dx)
{
  IJK_Field_local_double storage;
  storage.allocate(dx.ni()+1, dx.nj()+1, 2, 0);
  IJK_Field_local_double tmp[2];
  for (int i = 0; i < 2; i++)
    tmp[i].ref_ij(storage, i);
  IJK_Field_local_double *flux_zmin = &tmp[0];
  IJK_Field_local_double *flux_zmax = &tmp[1];
  compute_flux_z(*flux_zmin, 0);
  const int kmax = dx.nk();
  for (int k = 0; k < kmax; k++)
    {
      compute_flux_z(*flux_zmax, k+1);
      fill_grad_field_z_(*flux_zmin, *flux_zmax, dx, k);
      swap_data(flux_zmin, flux_zmax); // conserve le flux en z pour la couche z suivante
    }
}

void Operateur_IJK_elem_base_double::fill_grad_field_x_(IJK_Field_local_double& flux, FixedVector<IJK_Field_double, 3>& resu, int k)
{
  fill_grad_field_x_y_(flux, resu[0], k, 0);
}

void Operateur_IJK_elem_base_double::fill_grad_field_y_(IJK_Field_local_double& flux, FixedVector<IJK_Field_double, 3>& resu, int k)
{
  fill_grad_field_x_y_(flux, resu[1], k, 1);
}

void Operateur_IJK_elem_base_double::fill_grad_field_z_(IJK_Field_local_double& flux_min, IJK_Field_local_double& flux_max, FixedVector<IJK_Field_double, 3>& resu, int k)
{
  fill_grad_field_z_(flux_min, flux_max, resu[2], k);
}
