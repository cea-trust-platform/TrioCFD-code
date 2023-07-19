/****************************************************************************
 * Copyright (c) 2015 - 2016, CEA
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 *modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *this list of conditions and the following disclaimer in the documentation
 *and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *may be used to endorse or promote products derived from this software without
 *specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *****************************************************************************/
/////////////////////////////////////////////////////////////////////////////
//
// File      : OpConvIJKQuickScalar_interface.cpp
// Directory : $IJK_ROOT/src/IJK/OpVDF
//
/////////////////////////////////////////////////////////////////////////////
#include <Corrige_flux_FT.h>
#include <IJK_Field_simd_tools.h>
#include <OpConvIJKQuickScalar_interface.h>
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

Implemente_instanciable_sans_constructeur(OpConvQuickInterfaceIJKScalar_double, "OpConvQuickInterfaceIJKScalar_double", OpConvQuickIJKScalar_double);

Sortie& OpConvQuickInterfaceIJKScalar_double::printOn(Sortie& os) const
{
  //  OpConvIJKQuickScalar_double::printOn(os);
  return os;
}

Entree& OpConvQuickInterfaceIJKScalar_double::readOn(Entree& is)
{
  //  OpConvIJKQuickScalar_double::readOn(is);
  return is;
}

static inline void swap_data(IJK_Field_local_double *& a,
                             IJK_Field_local_double *& b)
{
  IJK_Field_local_double *tmp = a;
  a = b;
  b = tmp;
}

// Compute the "- divergence" of the given flux, integrated over control volume
// We assume than flux contain the integral of the flux on the faces of the
// control volumes For "div_set" method:
//  resu(i,j,k) = flux_x(i,j) - flux_x(i+1,j) + flux_y(i,j) - flux_y(i,j+1) +
//  flux_zmin(i,j) - flux_zmax(i,j)
// For "div_add" method: same but
//  resu(i,j,k) += ...
static void Operator_IJK_div_set(const IJK_Field_local_double& flux_x,
                                 const IJK_Field_local_double& flux_y,
                                 const IJK_Field_local_double& flux_zmin,
                                 const IJK_Field_local_double& flux_zmax,
                                 IJK_Field_local_double& resu, int k_layer)
{
  ConstIJK_double_ptr fx(flux_x, 0, 0, 0);
  ConstIJK_double_ptr fy(flux_y, 0, 0, 0);
  ConstIJK_double_ptr fzmin(flux_zmin, 0, 0, 0);
  ConstIJK_double_ptr fzmax(flux_zmax, 0, 0, 0);
  IJK_double_ptr resu_ptr(resu, 0, 0, k_layer);

  {
    const int imax = resu.ni();
    const int jmax = resu.nj();
    const int vsize = Simd_double::size();
    for (int j = 0;; j++)
      {
        for (int i = 0; i < imax; i += vsize)
          {
            Simd_double r, a, b;
            fx.get_center_right(DIRECTION::X, i, a, b);
            r = a - b;
            fy.get_center_right(DIRECTION::Y, i, a, b);
            r += a - b;
            fzmin.get_center(i, a);
            fzmax.get_center(i, b);
            r += a - b;
            resu_ptr.put_val(i, r);
          }
        // do not execute end_iloop at last iteration (because of assert on valid
        // j+1)
        if (j + 1 == jmax)
          break;
        fx.next_j();
        fy.next_j();
        fzmin.next_j();
        fzmax.next_j();
        resu_ptr.next_j();
      }
  }
}

// Finalement on va faire ces manipulations dans la correction de flux, comme ça
// il n'y a pas d'erreur possible. void
// OpConvQuickInterfaceIJKScalar_double::calculer(const IJK_Field_double&
// convect, const IJK_Field_double& disc_prop, const IJK_Field_double& vx, const
// IJK_Field_double& vy, const IJK_Field_double& vz, IJK_Field_double& res)
// {
//   // TODO: on calcule l'interpolation de convect aux faces
//   // on fait le produit avec la vitesse
//   // on fait le produit avec les propriétés de la manière suivante :
//   //   - si un des 2 voisin est pur, on utilise cette propriété
//   //   - sinon on fait le produit avec la moyenne arith avec l'indicatrice
//   moyenne (arith aussi)
//   // ensuite on corrige le flux avec la methode correction flux
// }

void OpConvQuickInterfaceIJKScalar_double::compute_set(IJK_Field_double& dx)
{
  IJK_Field_local_double storage;
  storage.allocate(dx.ni() + 1, dx.nj() + 1, 4, 0);
  IJK_Field_local_double tmp[4];
  for (int i = 0; i < 4; i++)
    tmp[i].ref_ij(storage, i);
  IJK_Field_local_double *flux_zmin = &tmp[0];
  IJK_Field_local_double *flux_zmax =
    &tmp[1]; // flux a travers la face zmax des volumes de controles
  IJK_Field_local_double *const flux_x =
    &tmp[2]; // pointer will not change, content will change
  IJK_Field_local_double *const flux_y = &tmp[3];

  // Une fois pour tout le domaine
  corrige_flux_->update();

  compute_flux_z(*flux_zmin, 0);
  corrige_flux_->corrige_flux_faceIJ(flux_zmin, -1, 2);

  const int kmax = dx.nk();
  for (int k = 0; k < kmax; k++)
    {
      compute_flux_x(*flux_x, k);
      compute_flux_y(*flux_y, k);
      compute_flux_z(*flux_zmax, k + 1);
      corrige_flux_->corrige_flux_faceIJ(flux_x, k, 0);
      corrige_flux_->corrige_flux_faceIJ(flux_y, k, 1);
      // correction_->corrige_flux_faceIJ(flux_zmin, k, 2);
      corrige_flux_->corrige_flux_faceIJ(flux_zmax, k, 2);
      Operator_IJK_div_set(*flux_x, *flux_y, *flux_zmin, *flux_zmax, dx, k);
      swap_data(flux_zmin,
                flux_zmax); // conserve le flux en z pour la couche z suivante
    }
}

// void OpConvQuickInterfaceIJKScalar_double::compute_add(IJK_Field_double& dx)
// {
// }

// OpConvQuickInterfaceIJKScalar_double::~OpConvQuickInterfaceIJKScalar_double()
// {
//   delete correction_;
// }
