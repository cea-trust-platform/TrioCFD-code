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
/////////////////////////////////////////////////////////////////////////////
//
// File      : IJK_One_Dimensional_Subproblem.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_One_Dimensional_Subproblem.h>

Implemente_instanciable_sans_constructeur( IJK_One_Dimensional_Subproblem, "IJK_One_Dimensional_Subproblem", Objet_U ) ;

// IJK_One_Dimensional_Subproblem::IJK_One_Dimensional_Subproblem() {}

IJK_One_Dimensional_Subproblem::IJK_One_Dimensional_Subproblem()
{
  points_per_thermal_subproblem_ = nullptr;
  points_per_thermal_subproblem_base_ = nullptr;
  alpha_ = nullptr;
  coeff_distance_diagonal_ = nullptr;
  cell_diagonal_ = nullptr;

  interfaces_ = nullptr;
  eulerian_distance_ = nullptr;
  eulerian_curvature_ = nullptr;
  eulerian_interfacial_area_ = nullptr;
  eulerian_normal_vect_ = nullptr;
  eulerian_facets_barycentre_ = nullptr;

  temperature_ = nullptr;
  temperature_ft_ = nullptr;
  velocity_ = nullptr;
  velocity_ft_ = nullptr;

  grad_T_elem_ = nullptr;
  hess_diag_T_elem_ = nullptr;
  hess_cross_T_elem_ = nullptr;

  radial_coordinates_ = nullptr;
  radial_coordinates_base_ = nullptr;

  radial_first_order_operator_raw_base_ = nullptr;
  radial_second_order_operator_raw_base_ = nullptr;
  radial_first_order_operator_base_ = nullptr;
  radial_second_order_operator_base_ = nullptr;
  radial_diffusion_matrix_base_ = nullptr;
  radial_convection_matrix_base_ = nullptr;
  radial_velocity_convection_matrix_base_=nullptr;
}

Sortie& IJK_One_Dimensional_Subproblem::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& IJK_One_Dimensional_Subproblem::readOn( Entree& is )
{
  Objet_U::readOn( is );
  return is;
}

void IJK_One_Dimensional_Subproblem::associate_eulerian_fields_references(const IJK_Interfaces& interfaces,
                                                                          const IJK_Field_local_double& eulerian_distance,
                                                                          const IJK_Field_local_double& eulerian_curvature,
                                                                          const IJK_Field_local_double& eulerian_interfacial_area,
                                                                          const FixedVector<IJK_Field_double, 3>& eulerian_normal_vect,
                                                                          const FixedVector<IJK_Field_double, 3>& eulerian_facets_barycentre,
                                                                          const IJK_Field_local_double& temperature,
                                                                          const IJK_Field_local_double& temperature_ft,
                                                                          const FixedVector<IJK_Field_double, 3>& velocity,
                                                                          const FixedVector<IJK_Field_double, 3>& velocity_ft,
                                                                          const FixedVector<IJK_Field_double, 3>& grad_T_elem,
                                                                          const FixedVector<IJK_Field_double, 3>& hess_diag_T_elem,
                                                                          const FixedVector<IJK_Field_double, 3>& hess_cross_T_elem)
{
  interfaces_ = &interfaces;
  eulerian_distance_ = &eulerian_distance;
  eulerian_curvature_ = &eulerian_curvature;
  eulerian_interfacial_area_ = &eulerian_interfacial_area;
  eulerian_normal_vect_ = &eulerian_normal_vect;
  eulerian_facets_barycentre_ = &eulerian_facets_barycentre;
  temperature_ = &temperature ;
  temperature_ft_ = &temperature_ft ;
  velocity_ = &velocity ;
  velocity_ft_ = &velocity_ft;
  grad_T_elem_ = &grad_T_elem ;
  hess_diag_T_elem_ = &hess_diag_T_elem ;
  hess_cross_T_elem_ = &hess_cross_T_elem ;
}

void IJK_One_Dimensional_Subproblem::associate_sub_problem_to_inputs(int sub_problem_index,
                                                                     int i, int j, int k, int compo_connex,
                                                                     double distance,
                                                                     double curvature,
                                                                     double interfacial_area,
                                                                     ArrOfDouble facet_barycentre,
                                                                     ArrOfDouble normal_vector,
                                                                     double bubble_rising_velocity,
                                                                     ArrOfDouble bubble_rising_vector,
                                                                     const int& points_per_thermal_subproblem,
                                                                     const double& alpha,
                                                                     const double& coeff_distance_diagonal,
                                                                     const double& cell_diagonal,
                                                                     const double& dr_base,
                                                                     const DoubleVect& radial_coordinates,
                                                                     const Matrice& radial_first_order_operator_raw,
                                                                     const Matrice& radial_second_order_operator_raw,
                                                                     const Matrice& radial_first_order_operator,
                                                                     const Matrice& radial_second_order_operator,
                                                                     Matrice& radial_diffusion_matrix,
                                                                     Matrice& radial_convection_matrix,
                                                                     const IJK_Interfaces& interfaces,
                                                                     const IJK_Field_local_double& eulerian_distance,
                                                                     const IJK_Field_local_double& eulerian_curvature,
                                                                     const IJK_Field_local_double& eulerian_interfacial_area,
                                                                     const FixedVector<IJK_Field_double, 3>& eulerian_normal_vect,
                                                                     const FixedVector<IJK_Field_double, 3>& eulerian_facets_barycentre,
                                                                     const IJK_Field_local_double& temperature,
                                                                     const IJK_Field_local_double& temperature_ft,
                                                                     const FixedVector<IJK_Field_double, 3>& velocity,
                                                                     const FixedVector<IJK_Field_double, 3>& velocity_ft,
                                                                     const FixedVector<IJK_Field_double, 3>& grad_T_elem,
                                                                     const FixedVector<IJK_Field_double, 3>& hess_diag_T_elem,
                                                                     const FixedVector<IJK_Field_double, 3>& hess_cross_T_elem)
{
  sub_problem_index_ = sub_problem_index;
  associate_cell_ijk(i, j, k);
  associate_compos(compo_connex);
  associate_interface_related_parameters(distance, curvature, interfacial_area, facet_barycentre, normal_vector);
  associate_rising_velocity(bubble_rising_velocity, bubble_rising_vector);
  associate_eulerian_fields_references(interfaces, eulerian_distance, eulerian_curvature, eulerian_interfacial_area, eulerian_normal_vect,
                                       eulerian_facets_barycentre, temperature, temperature_ft, velocity, velocity_ft, grad_T_elem, hess_diag_T_elem, hess_cross_T_elem);
  associate_probe_parameters(points_per_thermal_subproblem, alpha, coeff_distance_diagonal, cell_diagonal, dr_base, radial_coordinates);
  associate_finite_difference_operators(radial_first_order_operator_raw, radial_second_order_operator_raw,
                                        radial_first_order_operator, radial_second_order_operator,
                                        radial_diffusion_matrix, radial_convection_matrix);
  initialise_thermal_probe();
  // FIXME: What happen with highly deformable bubbles (concave interface portions) ?
}

void IJK_One_Dimensional_Subproblem::associate_probe_parameters(const int& points_per_thermal_subproblem,
                                                                const double& alpha,
                                                                const double& coeff_distance_diagonal,
                                                                const double& cell_diagonal,
                                                                const double& dr_base,
                                                                const DoubleVect& radial_coordinates)
{
  /*
   * If the probe characteristics are the same (don't copy attributes)
   */
  points_per_thermal_subproblem_base_ = &points_per_thermal_subproblem;
  coeff_distance_diagonal_ = &coeff_distance_diagonal;
  alpha_ = &alpha;
  cell_diagonal_ = &cell_diagonal;
  dr_base_ = &dr_base;
  radial_coordinates_base_ = &radial_coordinates;

  /*
   * If each probe differs (create attributes !)
   */
  if (global_probes_characteristics)
    points_per_thermal_subproblem_ = points_per_thermal_subproblem_base_;
  else
    points_per_thermal_subproblem_ = increase_number_of_points(); //copy if modified later
}

void IJK_One_Dimensional_Subproblem::associate_finite_difference_operators(const Matrice& radial_first_order_operator_raw,
                                                                           const Matrice& radial_second_order_operator_raw,
                                                                           const Matrice& radial_first_order_operator,
                                                                           const Matrice& radial_second_order_operator,
                                                                           Matrice& radial_diffusion_matrix,
                                                                           Matrice& radial_convection_matrix)
{
  radial_first_order_operator_raw_base_ = &radial_first_order_operator_raw;
  radial_second_order_operator_raw_base_ = &radial_second_order_operator_raw;
  radial_first_order_operator_base_ = &radial_first_order_operator;
  radial_second_order_operator_base_ = &radial_second_order_operator;
  radial_diffusion_matrix_base_ = &radial_diffusion_matrix;
  radial_convection_matrix_base_ = &radial_convection_matrix;
}

void IJK_One_Dimensional_Subproblem::initialise_thermal_probe()
{
  /*
   *  Curvature is negative for a convex bubble
   *  but R should be positive in that case
   */
  int i;
  if (fabs(curvature_) > DMINFLOAT)
    osculating_radius_ = fabs(2 / curvature_);
  if (global_probes_characteristics)
    {
      radial_coordinates_ = radial_coordinates_base_;
    }
  else
    {
      /*
       * coeff_distance_diagonal_ as well as
       * points_per_thermal_subproblem_ could be adapted
       */
      radial_coordinates_modified_ = DoubleVect(*points_per_thermal_subproblem_);
      radial_coordinates_ = &radial_coordinates_modified_;
      probe_length_ = (*coeff_distance_diagonal_) * (*cell_diagonal_);
      dr_ = probe_length_ / (*points_per_thermal_subproblem_ - 1);
      for (i=0; i < *points_per_thermal_subproblem_; i++)
        radial_coordinates_modified_(i) = i * dr_;
      // (*radial_coordinates_)(i) = i * dr_;
    }
  /*
   * Following attributes differ anyway !
   */
  osculating_radial_coordinates_ = DoubleVect(*points_per_thermal_subproblem_);
  radial_coordinates_cartesian_compo_ = DoubleTab(*points_per_thermal_subproblem_, 3);
  osculating_radial_coordinates_cartesian_compo_ = DoubleTab(*points_per_thermal_subproblem_, 3);
  osculating_radial_coordinates_ += osculating_radius_;
  osculating_radial_coordinates_inv_ = DoubleVect(*points_per_thermal_subproblem_);
  for (i=0; i < *points_per_thermal_subproblem_; i++)
    {
      osculating_radial_coordinates_inv_[i] = 1 / osculating_radial_coordinates_[i];
      for (int dir=0; dir<3; dir++)
        {
          radial_coordinates_cartesian_compo_(i, dir) = (*radial_coordinates_)(i) * normal_vector_compo_[dir];
          osculating_radial_coordinates_cartesian_compo_(i, dir) = osculating_radial_coordinates_(i) * normal_vector_compo_[dir];
        }
    }
}

const int * IJK_One_Dimensional_Subproblem::increase_number_of_points()
{
  increased_point_numbers_ = *points_per_thermal_subproblem_base_;
  return &increased_point_numbers_;
}

/*
 * TODO: Associate a basis to each subproblem
 * Use Rodrigues' rotation formula to determine ephi
 * Needs an axis of (rotation gravity_dir x relative_vectors)
 * and an angle (gravity_dir dot relative_vectors) / (norm(gravity_dir)*norm(relative_vectors))
 * ephi is determined in the gravity_align rising direction
 * 		 | gravity_dir
 * 		 |
 *   *****
 * ***   ***
 * **     **
 * ***   ***
 *   *****
 *     |
 *     |
 */

void IJK_One_Dimensional_Subproblem::interpolate_velocity_on_probes()
{

}

void IJK_One_Dimensional_Subproblem::reinitialise_local_finite_difference_operators()
{

}

void IJK_One_Dimensional_Subproblem::compute_radial_convection_diffusion_operators(IJK_Finite_Difference_One_Dimensional_Matrix_Assembler& finite_difference_assembler,
                                                                                   DoubleVect& thermal_subproblems_rhs_assembly,
                                                                                   const int& boundary_condition_interface,
                                                                                   const double& interfacial_boundary_condition_value,
                                                                                   const int& boundary_condition_end,
                                                                                   const double& impose_user_boundary_condition_end_value)
{
  interpolate_velocity_on_probes();
  const double alpha_inv = 1 / *alpha_;
  DoubleVect osculating_radial_coefficient = osculating_radial_coordinates_inv_;
  osculating_radial_coefficient *= 2;
  radial_convection_prefactor_ = radial_velocity_;
  radial_convection_prefactor_ *= alpha_inv;
  radial_convection_prefactor_ +=	osculating_radial_coefficient;
  const int boundary_conditions = 1;
  Matrice radial_convection_matrix_base_test;
  finite_difference_assembler.scale_matrix_subproblem_by_vector(radial_convection_matrix_base_,
                                                                radial_convection_prefactor_,
                                                                sub_problem_index_,
                                                                boundary_conditions);
}

