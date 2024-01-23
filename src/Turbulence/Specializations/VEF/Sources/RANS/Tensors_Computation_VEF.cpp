/****************************************************************************
 * Copyright (c) 2023, CEA
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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Tensors_Computation_VEF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#include <Tensors_Computation_VEF.h>
#include <TRUSTTab.h>
#include <Domaine_VEF.h>
#include <Domaine_Cl_VEF.h>
#include <Champ_P1NC.h>

using TCV = Tensors_Computation_VEF;

// Compute the antisymetric tensor
void TCV::compute_enstrophy(const Domaine_VEF& dom_VEF,
                            const Domaine_Cl_VEF& dom_BC_VEF,
                            const DoubleTab& velocity,
                            DoubleTab& enstrophy) const
{
  // Initialisation
  int nb_elem_tot = dom_VEF.nb_elem_tot();
  int dimension = Objet_U::dimension;
  DoubleTab gradient_elem(nb_elem_tot, dimension, dimension);
  gradient_elem = 0;

  // Computation of gradient
  Champ_P1NC::calcul_gradient(velocity, gradient_elem, dom_BC_VEF);

  // Loop on edge faces
  antisym_loop_edge_faces(dom_VEF, dom_BC_VEF, gradient_elem, enstrophy);

  // Loop on internal faces
  antisym_loop_internal_faces(dom_VEF, gradient_elem, enstrophy);

}

void TCV::antisym_loop_edge_faces(const Domaine_VEF& dom_VEF,
                                  const Domaine_Cl_VEF& dom_BC_VEF,
                                  const DoubleTab& gradient_elem,
                                  DoubleTab& enstrophy) const
{
  for (int n_edge = 0; n_edge < dom_VEF.nb_front_Cl(); ++n_edge)
    {
      const Cond_lim& current_BC = dom_BC_VEF.les_conditions_limites(n_edge);
      const Front_VF& the_edge = ref_cast(Front_VF, current_BC.frontiere_dis());

      if (sub_type(Periodique, current_BC.valeur()))
        antisym_loop_edges_periodiqueBC(dom_VEF, the_edge, gradient_elem, enstrophy);
      else
        antisym_loop_edges_general(dom_VEF, the_edge, gradient_elem, enstrophy);
    }
}

void TCV::antisym_loop_edges_periodiqueBC(const Domaine_VEF& dom_VEF,
                                          const Front_VF& the_edge,
                                          const DoubleTab& gradient_elem,
                                          DoubleTab& enstrophy) const
{
  const IntTab& neighbour_face = dom_VEF.face_voisins();
  const DoubleVect& volumes = dom_VEF.volumes();

  const int nbegin = the_edge.num_premiere_face();
  const int nend = nbegin + the_edge.nb_faces();

  for (int nface = nbegin; nface < nend; ++nface)
    {
      const int f0 = neighbour_face(nface, 0);
      const int f1 = neighbour_face(nface, 1);

      const double inv_vol = 1./(volumes(f0) + volumes(f1));
      const double vol0 = volumes(f0)*inv_vol;
      const double vol1 = volumes(f1)*inv_vol;

      // gradient_elem(elem_num, dimension, dimension)
      // const double du_dx = vol0*gradient_elem(f0, 0, 0) + vol1*gradient_elem(f1, 0, 0);
      const double du_dy = vol0*gradient_elem(f0, 0, 1) + vol1*gradient_elem(f1, 0, 1);
      const double dv_dx = vol0*gradient_elem(f0, 1, 0) + vol1*gradient_elem(f1, 1, 0);
      // const double dv_dy = vol0*gradient_elem(f0, 1, 1) + vol1*gradient_elem(f1, 1, 1);

      double tmp = 0.5*((du_dy - dv_dx)*(du_dy - dv_dx) +
                        (dv_dx - du_dy)*(dv_dx - du_dy));

      if (Objet_U::dimension == 3)
        {
          const double du_dz = vol0*gradient_elem(f0, 0, 2) + vol1*gradient_elem(f1, 0, 2);
          const double dv_dz = vol0*gradient_elem(f0, 1, 2) + vol1*gradient_elem(f1, 1, 2);
          const double dw_dx = vol0*gradient_elem(f0, 2, 0) + vol1*gradient_elem(f1, 2, 0);
          const double dw_dy = vol0*gradient_elem(f0, 2, 1) + vol1*gradient_elem(f1, 2, 1);
          // const double dw_dz = vol0*gradient_elem(f0, 2, 2) + vol1*gradient_elem(f1, 2, 2);

          tmp += 0.5*((du_dz - dw_dx)*(du_dz - dw_dx) +
                      (dv_dz - dw_dy)*(dv_dz - dw_dy) +
                      (dw_dx - du_dz)*(dw_dx - du_dz) +
                      (dw_dy - dv_dz)*(dw_dy - dv_dz));
        }

      enstrophy(nface) = sqrt(tmp);
    }
}


void TCV::antisym_loop_edges_general(const Domaine_VEF& dom_VEF,
                                     const Front_VF& the_edge,
                                     const DoubleTab& gradient_elem,
                                     DoubleTab& enstrophy) const
{
  const IntTab& neighbour_face = dom_VEF.face_voisins();

  const int nbegin = the_edge.num_premiere_face();
  const int nend = nbegin + the_edge.nb_faces();

  for (int nface = nbegin; nface < nend; ++nface)
    {
      const int f0 = neighbour_face(nface, 0);

      // gradient_elem(elem_num, dimension, dimension)
      // const double du_dx = gradient_elem(f0, 0, 0);
      const double du_dy = gradient_elem(f0, 0, 1);
      const double dv_dx = gradient_elem(f0, 1, 0);
      // const double dv_dy = gradient_elem(f0, 1, 1);

      double tmp = 0.5*((du_dy - dv_dx)*(du_dy - dv_dx) +
                        (dv_dx - du_dy)*(dv_dx - du_dy));

      if (Objet_U::dimension == 3)
        {
          const double du_dz = gradient_elem(f0, 0, 2);
          const double dv_dz = gradient_elem(f0, 1, 2);
          const double dw_dx = gradient_elem(f0, 2, 0);
          const double dw_dy = gradient_elem(f0, 2, 1);
          // const double dw_dz = gradient_elem(f0, 2, 2);

          tmp += 0.5*((du_dz - dw_dx)*(du_dz - dw_dx) +
                      (dv_dz - dw_dy)*(dv_dz - dw_dy) +
                      (dw_dx - du_dz)*(dw_dx - du_dz) +
                      (dw_dy - dv_dz)*(dw_dy - dv_dz));
        }

      enstrophy(nface) = sqrt(tmp);
    }
}

// Même fonction que périodique !
void TCV::antisym_loop_internal_faces(const Domaine_VEF& dom_VEF,
                                      const DoubleTab& grad_elem,
                                      DoubleTab& enstrophy) const
{
  const IntTab& neighbour_face = dom_VEF.face_voisins();
  const DoubleVect& volumes = dom_VEF.volumes();
  int first_internal_face = dom_VEF.premiere_face_int();
  int nb_faces = dom_VEF.nb_faces();

  for (int nface = first_internal_face; nface < nb_faces; ++nface)
    {
      const int f0 = neighbour_face(nface, 0);
      const int f1 = neighbour_face(nface, 1);

      const double inv_vol = 1./(volumes(f0) + volumes(f1));
      const double vol0 = volumes(f0)*inv_vol;
      const double vol1 = volumes(f1)*inv_vol;

      // format: grad_elem(elem_num, dimension, dimension)
      // const double du_dx = vol0*grad_elem(f0, 0, 0) + vol1*grad_elem(f1, 0, 0);
      const double du_dy = vol0*grad_elem(f0, 0, 1) + vol1*grad_elem(f1, 0, 1);
      const double dv_dx = vol0*grad_elem(f0, 1, 0) + vol1*grad_elem(f1, 1, 0);
      // const double dv_dy = vol0*grad_elem(f0, 1, 1) + vol1*grad_elem(f1, 1, 1);

      double tmp = 0.5*((du_dy - dv_dx)*(du_dy - dv_dx) +
                        (dv_dx - du_dy)*(dv_dx - du_dy));

      if (Objet_U::dimension == 3)
        {
          const double du_dz = vol0*grad_elem(f0, 0, 2) + vol1*grad_elem(f1, 0, 2);
          const double dv_dz = vol0*grad_elem(f0, 1, 2) + vol1*grad_elem(f1, 1, 2);
          const double dw_dx = vol0*grad_elem(f0, 2, 0) + vol1*grad_elem(f1, 2, 0);
          const double dw_dy = vol0*grad_elem(f0, 2, 1) + vol1*grad_elem(f1, 2, 1);
          // const double dw_dz = vol0*grad_elem(f0, 2, 2) + vol1*grad_elem(f1, 2, 2);

          tmp += 0.5*((du_dz - dw_dx)*(du_dz - dw_dx) +
                      (dv_dz - dw_dy)*(dv_dz - dw_dy) +
                      (dw_dx - du_dz)*(dw_dx - du_dz) +
                      (dw_dy - dv_dz)*(dw_dy - dv_dz));
        }

      enstrophy(nface) = sqrt(tmp);
    }
}
