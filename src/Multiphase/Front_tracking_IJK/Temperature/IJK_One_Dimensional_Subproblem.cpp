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
#include <IJK_Navier_Stokes_tools.h>
#include <Ouvrir_fichier.h>
#include <IJK_FT.h>

Implemente_instanciable_sans_constructeur( IJK_One_Dimensional_Subproblem, "IJK_One_Dimensional_Subproblem", Objet_U ) ;

IJK_One_Dimensional_Subproblem::IJK_One_Dimensional_Subproblem()
{
  points_per_thermal_subproblem_ = nullptr;
  points_per_thermal_subproblem_base_ = nullptr;
  alpha_ = nullptr;
  lambda_ = nullptr;
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

  radial_first_order_operator_=nullptr;
  radial_second_order_operator_=nullptr;
  finite_difference_assembler_=nullptr;
  thermal_subproblems_matrix_assembly_=nullptr;
  thermal_subproblems_rhs_assembly_=nullptr;
  thermal_subproblems_temperature_solution_=nullptr;

  interfacial_boundary_condition_value_ = 0.;
  end_boundary_condition_value_ = 0.;
  start_index_ = 0;
  end_index_ = 0;

  correct_radial_velocity_ = 1;
  correct_tangential_temperature_gradient_ = 0;
  correct_tangential_temperature_hessian_ = 0;

  first_tangential_vector_compo_solver_=nullptr;
  second_tangential_vector_compo_solver_=nullptr;
}

IJK_One_Dimensional_Subproblem::IJK_One_Dimensional_Subproblem(const IJK_FT_double& ijk_ft) : IJK_One_Dimensional_Subproblem()
{
  associer(ijk_ft);
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

void IJK_One_Dimensional_Subproblem::get_ijk_indices(int& i, int& j, int& k) const
{
  i = (int) index_i_;
  j = (int) index_j_;
  k = (int) index_k_;
}

void IJK_One_Dimensional_Subproblem::associate_eulerian_fields_references(const IJK_Interfaces& interfaces,
                                                                          const IJK_Field_double& eulerian_distance,
                                                                          const IJK_Field_double& eulerian_curvature,
                                                                          const IJK_Field_double& eulerian_interfacial_area,
                                                                          const FixedVector<IJK_Field_double, 3>& eulerian_normal_vect,
                                                                          const FixedVector<IJK_Field_double, 3>& eulerian_facets_barycentre,
                                                                          const IJK_Field_double& temperature,
                                                                          const IJK_Field_double& temperature_ft,
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
                                                                     ArrOfDouble bubble_barycentre,
                                                                     const int& points_per_thermal_subproblem,
                                                                     const double& alpha,
                                                                     const double& lambda,
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
                                                                     const IJK_Field_double& eulerian_distance,
                                                                     const IJK_Field_double& eulerian_curvature,
                                                                     const IJK_Field_double& eulerian_interfacial_area,
                                                                     const FixedVector<IJK_Field_double, 3>& eulerian_normal_vect,
                                                                     const FixedVector<IJK_Field_double, 3>& eulerian_facets_barycentre,
                                                                     const IJK_Field_double& temperature,
                                                                     const IJK_Field_double& temperature_ft,
                                                                     const FixedVector<IJK_Field_double, 3>& velocity,
                                                                     const FixedVector<IJK_Field_double, 3>& velocity_ft,
                                                                     const FixedVector<IJK_Field_double, 3>& grad_T_elem,
                                                                     const FixedVector<IJK_Field_double, 3>& hess_diag_T_elem,
                                                                     const FixedVector<IJK_Field_double, 3>& hess_cross_T_elem,
                                                                     IJK_Finite_Difference_One_Dimensional_Matrix_Assembler& finite_difference_assembler,
                                                                     Matrice& thermal_subproblems_matrix_assembly,
                                                                     DoubleVect& thermal_subproblems_rhs_assembly,
                                                                     DoubleVect& thermal_subproblems_temperature_solution,
                                                                     const int& source_terms_type,
                                                                     const int& source_terms_correction)
{
  sub_problem_index_ = sub_problem_index;
  associate_cell_ijk(i, j, k);
  associate_compos(compo_connex);
  associate_interface_related_parameters(distance, curvature, interfacial_area, facet_barycentre, normal_vector);
  associate_rising_velocity(bubble_rising_velocity, bubble_rising_vector, bubble_barycentre);
  associate_eulerian_fields_references(interfaces, eulerian_distance, eulerian_curvature, eulerian_interfacial_area, eulerian_normal_vect,
                                       eulerian_facets_barycentre, temperature, temperature_ft, velocity, velocity_ft, grad_T_elem, hess_diag_T_elem, hess_cross_T_elem);
  associate_probe_parameters(points_per_thermal_subproblem, alpha, lambda, coeff_distance_diagonal, cell_diagonal, dr_base, radial_coordinates);
  associate_finite_difference_operators(radial_first_order_operator_raw, radial_second_order_operator_raw,
                                        radial_first_order_operator, radial_second_order_operator,
                                        radial_diffusion_matrix, radial_convection_matrix);
  finite_difference_assembler_ = &finite_difference_assembler;
  thermal_subproblems_matrix_assembly_ = &thermal_subproblems_matrix_assembly;
  thermal_subproblems_rhs_assembly_ = &thermal_subproblems_rhs_assembly;
  thermal_subproblems_temperature_solution_ = &thermal_subproblems_temperature_solution;
  source_terms_type_ = source_terms_type;
  correct_tangential_temperature_gradient_ = source_terms_correction;
  correct_tangential_temperature_hessian_ = source_terms_correction;
  initialise_thermal_probe();
  recompute_finite_difference_matrices();
  // FIXME: What happen with highly deformable bubbles (concave interface portions) ?
}

void IJK_One_Dimensional_Subproblem::associate_probe_parameters(const int& points_per_thermal_subproblem,
                                                                const double& alpha,
                                                                const double& lambda,
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
  lambda_ = &lambda;
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
  compute_interface_basis_vectors();
  compute_pure_spherical_basis_vectors();
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
      dr_ = *dr_base_;
    }
  else
    {
      /*
       * coeff_distance_diagonal_ as well as
       * points_per_thermal_subproblem_ could be adapted
       */
      radial_coordinates_modified_.resize(*points_per_thermal_subproblem_);
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
  dr_inv_ = 1 / dr_;
  osculating_radial_coordinates_ = (*radial_coordinates_); // DoubleVect(*points_per_thermal_subproblem_);
  radial_coordinates_cartesian_compo_.resize(*points_per_thermal_subproblem_, 3);
  coordinates_cartesian_compo_.resize(*points_per_thermal_subproblem_, 3);
  osculating_radial_coordinates_cartesian_compo_.resize(*points_per_thermal_subproblem_, 3);
  osculating_radial_coordinates_inv_.resize(*points_per_thermal_subproblem_);
  osculating_radial_coordinates_ += osculating_radius_;
  for (i=0; i < *points_per_thermal_subproblem_; i++)
    {
      osculating_radial_coordinates_inv_[i] = 1 / osculating_radial_coordinates_[i];
      for (int dir=0; dir<3; dir++)
        {
          radial_coordinates_cartesian_compo_(i, dir) = (*radial_coordinates_)(i) * normal_vector_compo_[dir];
          osculating_radial_coordinates_cartesian_compo_(i, dir) = osculating_radial_coordinates_(i) * normal_vector_compo_[dir];
          coordinates_cartesian_compo_(i, dir) = radial_coordinates_cartesian_compo_(i, dir) + facet_barycentre_[dir];
        }
    }
  surface_ = (*eulerian_interfacial_area_)(index_i_, index_j_, index_k_);
  rhs_assembly_.resize(*points_per_thermal_subproblem_);
}

void IJK_One_Dimensional_Subproblem::compute_interface_basis_vectors()
{
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
//	bubble_rising_velocity_ = bubble_rising_velocity;
//	bubble_rising_vector_ = Vecteur3(bubble_rising_vector);
//	facet_barycentre_;
//	normal_vector_compo_;

  facet_barycentre_relative_ = facet_barycentre_ - bubble_barycentre_;
  Vecteur3 facet_barycentre_relative_normed = facet_barycentre_relative_;
  const double facet_barycentre_relative_norm = facet_barycentre_relative_normed.length();
  facet_barycentre_relative_normed *= (1 / facet_barycentre_relative_norm);
  Vecteur3 normal_contrib;
  const double normal_vector_compo_norm = normal_vector_compo_.length();
  normal_vector_compo_ *= (1 / normal_vector_compo_norm);
  /*
   * First method with tangential direction of maximum velocity variations
   */

  DoubleTab facet_barycentre(1, 3);
  interfacial_velocity_compo_ = 0.;
  for (int dir=0; dir<3; dir++)
    facet_barycentre(0, dir) = facet_barycentre_[dir];
  for (int dir=0; dir<3; dir++)
    {
      DoubleVect interfacial_velocity_component(1);
      ijk_interpolate_skip_unknown_points((*velocity_)[dir], facet_barycentre, interfacial_velocity_component, INVALID_INTERP);
      interfacial_velocity_compo_[0] = interfacial_velocity_component[0];
    }
  if (interfacial_velocity_compo_.length() < INVALID_VELOCITY)
    {
      normal_contrib = normal_vector_compo_;
      normal_contrib *= Vecteur3::produit_scalaire(facet_barycentre_relative_normed, normal_vector_compo_);
      first_tangential_vector_compo_ = facet_barycentre_relative_normed - normal_contrib;
    }
  else
    {
      normal_contrib = normal_vector_compo_;
      normal_contrib *= Vecteur3::produit_scalaire(interfacial_velocity_compo_, normal_vector_compo_);
      interfacial_tangential_velocity_compo_ = interfacial_velocity_compo_ - normal_contrib;
      first_tangential_vector_compo_ = interfacial_tangential_velocity_compo_;
    }
  const double norm_first_tangential_vector = first_tangential_vector_compo_.length();
  first_tangential_vector_compo_ *= (1 / norm_first_tangential_vector);
  Vecteur3::produit_vectoriel(normal_vector_compo_, first_tangential_vector_compo_, second_tangential_vector_compo_);
  const double norm_second_tangential_vector = second_tangential_vector_compo_.length();
  second_tangential_vector_compo_ *= (1 / norm_second_tangential_vector);
  /*
   * Second method with rising velocity
   */
  Vecteur3::produit_vectoriel(bubble_rising_vector_, facet_barycentre_relative_, azymuthal_vector_compo_raw_);

  const int sign_vector = signbit(Vecteur3::produit_scalaire(bubble_rising_vector_, normal_vector_compo_));
  azymuthal_vector_compo_ = azymuthal_vector_compo_raw_;
  const double norm_azymuthal_vector_compo_raw_ = sqrt(Vecteur3::produit_scalaire(azymuthal_vector_compo_raw_, azymuthal_vector_compo_raw_));
  if (sign_vector)
    azymuthal_vector_compo_ *= -1;
  azymuthal_vector_compo_ *= norm_azymuthal_vector_compo_raw_;

  normal_contrib = normal_vector_compo_;
  normal_contrib *=	Vecteur3::produit_scalaire(azymuthal_vector_compo_, normal_vector_compo_);
  azymuthal_vector_compo_ = azymuthal_vector_compo_ - normal_contrib;
  Vecteur3::produit_vectoriel(azymuthal_vector_compo_, normal_vector_compo_, first_tangential_vector_compo_from_rising_dir_);
  const double norm_second_tangential_vector_from_rising_dir = first_tangential_vector_compo_from_rising_dir_.length();
  first_tangential_vector_compo_from_rising_dir_ *= (1 / norm_second_tangential_vector_from_rising_dir);


  if (tangential_from_rising_vel)
    {
      first_tangential_vector_compo_solver_ = &first_tangential_vector_compo_from_rising_dir_;
      second_tangential_vector_compo_solver_ = &azymuthal_vector_compo_;
    }
  else
    {
      // By default
      first_tangential_vector_compo_solver_= &first_tangential_vector_compo_;
      second_tangential_vector_compo_solver_ = &second_tangential_vector_compo_;
    }
}

void IJK_One_Dimensional_Subproblem::compute_pure_spherical_basis_vectors()
{
  /*
   * FIXME: It is align with gravity z but it should be modified to be align with the gravity dir
   */
  r_sph_ = sqrt(facet_barycentre_relative_[0] * facet_barycentre_relative_[0]
                + facet_barycentre_relative_[1] * facet_barycentre_relative_[1]
                + facet_barycentre_relative_[2] * facet_barycentre_relative_[2]);
  theta_sph_ = atan(sqrt(facet_barycentre_relative_[0] * facet_barycentre_relative_[0]
                         + facet_barycentre_relative_[1] * facet_barycentre_relative_[1])/ facet_barycentre_relative_[2]);
  phi_sph_ = atan(facet_barycentre_relative_[1] / facet_barycentre_relative_[0]);

  for (int dir=0; dir<3; dir++)
    er_sph_[dir] = facet_barycentre_relative_[dir] / r_sph_;

  const double length = sqrt(facet_barycentre_relative_[0] * facet_barycentre_relative_[0]
                             + facet_barycentre_relative_[1] * facet_barycentre_relative_[1]);
  for (int dir=0; dir<2; dir++)
    etheta_sph_[dir] = facet_barycentre_relative_[dir] * facet_barycentre_relative_[2] / (r_sph_ * length);
  etheta_sph_[2] = - facet_barycentre_relative_[2] * length / r_sph_;

  ephi_sph_ = {0., 0., 0.};
  ephi_sph_[0] = - facet_barycentre_relative_[1];
  ephi_sph_[1] = facet_barycentre_relative_[0];
}

const int * IJK_One_Dimensional_Subproblem::increase_number_of_points()
{
  increased_point_numbers_ = *points_per_thermal_subproblem_base_;
  return &increased_point_numbers_;
}

void IJK_One_Dimensional_Subproblem::compute_first_order_operator_local(Matrice& radial_first_order_operator)
{
  int check_nb_elem;
  check_nb_elem = (*finite_difference_assembler_).build(radial_first_order_operator, *points_per_thermal_subproblem_, 0);
  Cerr << "Check_nb_elem" << check_nb_elem << finl;
  radial_first_order_operator_local_ *= dr_inv_;
}

void IJK_One_Dimensional_Subproblem::compute_second_order_operator_local(Matrice& radial_second_order_operator)
{
  int check_nb_elem;
  check_nb_elem = (*finite_difference_assembler_).build(radial_second_order_operator, *points_per_thermal_subproblem_, 0);
  Cerr << "Check_nb_elem" << check_nb_elem << finl;
  const double dr_squared_inv = 1 / pow(dr_, 2);
  radial_second_order_operator_local_ *= dr_squared_inv;
}

void IJK_One_Dimensional_Subproblem::recompute_finite_difference_matrices()
{
  if (!global_probes_characteristics)
    {
      compute_first_order_operator_local(radial_first_order_operator_local_);
      compute_second_order_operator_local(radial_second_order_operator_local_);
    }
}

void IJK_One_Dimensional_Subproblem::interpolate_project_velocities_on_probes()
{
  radial_velocity_.resize(*points_per_thermal_subproblem_);
  first_tangential_velocity_.resize(*points_per_thermal_subproblem_);
  second_tangential_velocity_.resize(*points_per_thermal_subproblem_);
  azymuthal_velocity_.resize(*points_per_thermal_subproblem_);

  x_velocity_.resize(*points_per_thermal_subproblem_);
  y_velocity_.resize(*points_per_thermal_subproblem_);
  z_velocity_.resize(*points_per_thermal_subproblem_);

  interpolate_cartesian_velocities_on_probes();
  project_velocities_on_probes();
}

void IJK_One_Dimensional_Subproblem::interpolate_cartesian_velocities_on_probes()
{
  DoubleVect * vel_compo = nullptr;
  for (int dir = 0; dir < 3; dir++)
    {
      switch (dir)
        {
        case 0:
          vel_compo = &x_velocity_;
          break;
        case 1:
          vel_compo = &y_velocity_;
          break;
        case 2:
          vel_compo = &z_velocity_;
          break;
        default:
          vel_compo = &x_velocity_;
          break;
        }
      ijk_interpolate_skip_unknown_points((*velocity_)[dir], coordinates_cartesian_compo_, *vel_compo, INVALID_INTERP);
    }
}

void IJK_One_Dimensional_Subproblem::project_velocities_on_probes()
{
  project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, normal_vector_compo_, radial_velocity_);
  project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, *first_tangential_vector_compo_solver_, first_tangential_velocity_);
  project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, *second_tangential_vector_compo_solver_, second_tangential_velocity_);
  project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, azymuthal_vector_compo_, azymuthal_velocity_);
  correct_radial_velocity();
}

void IJK_One_Dimensional_Subproblem::correct_radial_velocity()
{
  radial_velocity_corrected_ = radial_velocity_;
  for (int i=0; i<radial_velocity_corrected_.size(); i++)
    radial_velocity_corrected_[i] -= radial_velocity_[0];
}

void IJK_One_Dimensional_Subproblem::project_cartesian_onto_basis_vector(const DoubleVect& compo_x, const DoubleVect& compo_y, const DoubleVect& compo_z, const Vecteur3& basis, DoubleVect& projection)
{
  const int size_x = compo_x.size();
//  const int size_y = compo_y.size();
//  const int size_z = compo_z.size();
//  const int size_projection = projection.size();
//  assert((size_x == size_y) && (size_x == size_z) && (size_x == size_projection));
//  Cerr << "size_x << size_y << size_z << size_projection" <<  finl;
//  Cerr << size_x << size_y << size_z << size_projection <<  finl;
  for (int i=0; i<size_x; i++)
    projection[i] = compo_x[i] * basis[0] + compo_y[i] * basis[1] + compo_z[i] * basis[2];
}

void IJK_One_Dimensional_Subproblem::project_basis_vector_onto_cartesian_dir(const int& dir, const DoubleVect& compo_u, const DoubleVect& compo_v, const DoubleVect& compo_w,
                                                                             const Vecteur3& basis_u, const Vecteur3& basis_v, const Vecteur3& basis_w,
                                                                             DoubleVect& projection)
{
  const int size_u = compo_u.size();
  const int size_v = compo_v.size();
  const int size_w = compo_w.size();
  for (int i=0; i<projection.size(); i++)
    {
      if (i< size_u)
        projection[i] += compo_u[i] * basis_u[dir];
      if (i< size_v)
        projection[i] += compo_v[i] * basis_v[dir];
      if (i< size_w)
        projection[i] += compo_w[i] * basis_w[dir];
    }
}

void IJK_One_Dimensional_Subproblem::interpolate_temperature_on_probe()
{
  temperature_interp_.resize(*points_per_thermal_subproblem_);
  ijk_interpolate_skip_unknown_points(*temperature_, coordinates_cartesian_compo_, temperature_interp_, INVALID_INTERP);
}

void IJK_One_Dimensional_Subproblem::interpolate_temperature_gradient_on_probe()
{
  for (int dir = 0; dir < 3; dir++)
    {
      grad_T_elem_interp_[dir].resize(*points_per_thermal_subproblem_);
      ijk_interpolate_skip_unknown_points((*grad_T_elem_)[dir], coordinates_cartesian_compo_, grad_T_elem_interp_[dir], INVALID_INTERP);
    }
}

void IJK_One_Dimensional_Subproblem::project_temperature_gradient_on_probes()
{
  normal_temperature_gradient_.resize(*points_per_thermal_subproblem_);
  project_cartesian_onto_basis_vector(grad_T_elem_interp_[0], grad_T_elem_interp_[1], grad_T_elem_interp_[2], normal_vector_compo_, normal_temperature_gradient_);

  tangential_temperature_gradient_first_.resize(*points_per_thermal_subproblem_);
  tangential_temperature_gradient_second_.resize(*points_per_thermal_subproblem_);
  project_cartesian_onto_basis_vector(grad_T_elem_interp_[0], grad_T_elem_interp_[1], grad_T_elem_interp_[2], *first_tangential_vector_compo_solver_, tangential_temperature_gradient_first_);
  project_cartesian_onto_basis_vector(grad_T_elem_interp_[0], grad_T_elem_interp_[1], grad_T_elem_interp_[2], *second_tangential_vector_compo_solver_, tangential_temperature_gradient_second_);

  azymuthal_temperature_gradient_.resize(*points_per_thermal_subproblem_);
  project_cartesian_onto_basis_vector(grad_T_elem_interp_[0], grad_T_elem_interp_[1], grad_T_elem_interp_[2], azymuthal_vector_compo_, azymuthal_temperature_gradient_);
}

void IJK_One_Dimensional_Subproblem::interpolate_temperature_hessian_on_probe()
{
  for (int dir = 0; dir < 3; dir++)
    {
      hess_diag_T_elem_interp_[dir].resize(*points_per_thermal_subproblem_);
      hess_cross_T_elem_interp_[dir].resize(*points_per_thermal_subproblem_);
      ijk_interpolate_skip_unknown_points((*hess_diag_T_elem_)[dir], coordinates_cartesian_compo_, hess_diag_T_elem_interp_[dir], INVALID_INTERP);
      ijk_interpolate_skip_unknown_points((*hess_cross_T_elem_)[dir], coordinates_cartesian_compo_, hess_cross_T_elem_interp_[dir], INVALID_INTERP);
    }
}

void IJK_One_Dimensional_Subproblem::project_temperature_hessian_on_probes()
{

}

void IJK_One_Dimensional_Subproblem::initialise_radial_convection_operator_local()
{
  if (!global_probes_characteristics)
    (*finite_difference_assembler_).reinitialise_matrix_subproblem(*radial_convection_matrix_base_, radial_first_order_operator_local_, sub_problem_index_);
}

void IJK_One_Dimensional_Subproblem::initialise_radial_diffusion_operator_local()
{
  if (!global_probes_characteristics)
    (*finite_difference_assembler_).reinitialise_matrix_subproblem(*radial_diffusion_matrix_base_, radial_second_order_operator_local_, sub_problem_index_);
}

void IJK_One_Dimensional_Subproblem::compute_radial_convection_diffusion_operators()
{
  const double alpha_inv = - 1 / *alpha_;
  DoubleVect osculating_radial_coefficient = osculating_radial_coordinates_inv_;
  osculating_radial_coefficient *= 2;
  radial_convection_prefactor_.resize(*points_per_thermal_subproblem_);
  if (source_terms_type_ != linear_diffusion && source_terms_type_ != spherical_diffusion)
    {
      interpolate_project_velocities_on_probes();
      if (correct_radial_velocity_)
        radial_convection_prefactor_ = radial_velocity_corrected_;
      else
        radial_convection_prefactor_ = radial_velocity_;
      radial_convection_prefactor_ *= alpha_inv;
    }
  else
    {
      radial_velocity_.resize(*points_per_thermal_subproblem_);
      radial_velocity_corrected_.resize(*points_per_thermal_subproblem_);
      first_tangential_velocity_.resize(*points_per_thermal_subproblem_);
      second_tangential_velocity_.resize(*points_per_thermal_subproblem_);
      azymuthal_velocity_.resize(*points_per_thermal_subproblem_);
    }
  if (source_terms_type_ != linear_diffusion)
    radial_convection_prefactor_ +=	osculating_radial_coefficient;
  const int boundary_conditions = 0;
  initialise_radial_convection_operator_local();
  initialise_radial_diffusion_operator_local();
  (*finite_difference_assembler_).scale_matrix_subproblem_by_vector(radial_convection_matrix_base_,
                                                                    radial_convection_prefactor_,
                                                                    sub_problem_index_,
                                                                    boundary_conditions);
}

void IJK_One_Dimensional_Subproblem::impose_boundary_conditions(DoubleVect& thermal_subproblems_rhs_assembly,
                                                                const int& boundary_condition_interface,
                                                                const double& interfacial_boundary_condition_value,
                                                                const int& impose_boundary_condition_interface_from_simulation,
                                                                const int& boundary_condition_end,
                                                                const double& end_boundary_condition_value,
                                                                const int& impose_user_boundary_condition_end_value)
{
  interpolate_temperature_on_probe();
  for (int i=0; i<temperature_interp_.size(); i++)
    if (fabs(temperature_interp_[i]) > INVALID_INTERP_TEST)
      Cerr << "Error in the temperature_interpolation" << temperature_interp_[i] << finl;
  if (!impose_boundary_condition_interface_from_simulation)
    interfacial_boundary_condition_value_ = interfacial_boundary_condition_value;
  else
    interfacial_boundary_condition_value_ = temperature_interp_[0];
  if (impose_user_boundary_condition_end_value)
    end_boundary_condition_value_ = end_boundary_condition_value;
  else
    end_boundary_condition_value_ = temperature_interp_[temperature_interp_.size() -1];

  const int thermal_subproblems_rhs_size = (*thermal_subproblems_rhs_assembly_).size();
  start_index_ = thermal_subproblems_rhs_size;
  end_index_ = start_index_ + (*points_per_thermal_subproblem_);

  thermal_subproblems_rhs_assembly.resize(thermal_subproblems_rhs_size + *points_per_thermal_subproblem_);

  (*finite_difference_assembler_).impose_boundary_conditions_subproblem(thermal_subproblems_matrix_assembly_,
                                                                        thermal_subproblems_rhs_assembly_,
                                                                        rhs_assembly_,
                                                                        boundary_condition_interface,
                                                                        interfacial_boundary_condition_value_,
                                                                        boundary_condition_end,
                                                                        end_boundary_condition_value_,
                                                                        sub_problem_index_,
                                                                        dr_inv_);
}

void IJK_One_Dimensional_Subproblem::compute_add_source_terms()
{
  if (source_terms_type_ != linear_diffusion && source_terms_type_ != spherical_diffusion)
    {
      interpolate_temperature_gradient_on_probe();
      project_temperature_gradient_on_probes();

      if (source_terms_type_ == tangential_conv_2D_tangential_diffusion_2D ||
          source_terms_type_ == tangential_conv_3D_tangentual_diffusion_3D)
        {
          interpolate_temperature_hessian_on_probe();
          project_temperature_hessian_on_probes();
        }
      /*
       * Compute at least the tangential convection
       */
      tangential_convection_source_terms_first_ = tangential_temperature_gradient_first_;
      if (correct_tangential_temperature_gradient_)
        correct_tangential_temperature_gradient(tangential_convection_source_terms_first_);
      tangential_convection_source_terms_first_ *= first_tangential_velocity_;
      const double alpha_inv = 1 / *alpha_;
      tangential_convection_source_terms_first_ *= alpha_inv;
      switch (source_terms_type_)
        {
        case 0:
          tangential_convection_source_terms_ = tangential_convection_source_terms_first_;
          source_terms_ = tangential_convection_source_terms_;
          break;
        case 1:
          tangential_convection_source_terms_second_ = tangential_temperature_gradient_second_;
          if (correct_tangential_temperature_gradient_)
            correct_tangential_temperature_gradient(tangential_convection_source_terms_second_);
          tangential_convection_source_terms_second_ *= second_tangential_velocity_;
          tangential_convection_source_terms_second_ *= alpha_inv;
          tangential_convection_source_terms_ = tangential_convection_source_terms_first_;
          tangential_convection_source_terms_ += tangential_convection_source_terms_second_;
          source_terms_ = tangential_convection_source_terms_;
          break;
        case 2:
          tangential_convection_source_terms_ = tangential_convection_source_terms_first_;
          source_terms_ = tangential_convection_source_terms_;
          // tangential_diffusion_source_terms_;
          // if (correct_tangential_temperature_hessian_)
          // 	correct_tangential_temperature_hessian(tangential_diffusion_source_terms_)
          source_terms_ += tangential_diffusion_source_terms_;
          break;
        case 3:
          tangential_convection_source_terms_second_ = tangential_temperature_gradient_second_;
          if (correct_tangential_temperature_gradient_)
            correct_tangential_temperature_gradient(tangential_convection_source_terms_second_);
          tangential_convection_source_terms_second_ *= second_tangential_velocity_;
          tangential_convection_source_terms_second_ *= alpha_inv;
          tangential_convection_source_terms_ = tangential_convection_source_terms_first_;
          tangential_convection_source_terms_ += tangential_convection_source_terms_second_;
          source_terms_ = tangential_convection_source_terms_;
          // tangential_diffusion_source_terms_;
          // if (correct_tangential_temperature_hessian_)
          // 	correct_tangential_temperature_hessian(tangential_diffusion_source_terms_)
          source_terms_ += tangential_diffusion_source_terms_;
          break;
        default:
          tangential_convection_source_terms_ = tangential_convection_source_terms_first_;
          source_terms_ = tangential_convection_source_terms_;
          break;
        }

      (*finite_difference_assembler_).add_source_terms(thermal_subproblems_rhs_assembly_, rhs_assembly_);
    }
}

void IJK_One_Dimensional_Subproblem::correct_tangential_temperature_gradient(DoubleVect& tangential_convection_source_terms)
{
  DoubleVect tangential_convection_source_terms_tmp = tangential_convection_source_terms;
  for (int i=0; i<tangential_convection_source_terms.size(); i++)
    tangential_convection_source_terms[i] -= tangential_convection_source_terms_tmp[0];
}

void IJK_One_Dimensional_Subproblem::correct_tangential_temperature_hessian(DoubleVect& tangential_diffusion_source_terms)
{
  DoubleVect tangential_diffusion_source_terms_tmp = tangential_diffusion_source_terms_tmp;
  for (int i=0; i<tangential_diffusion_source_terms_tmp.size(); i++)
    tangential_diffusion_source_terms_tmp[i] -= tangential_diffusion_source_terms_tmp[0];
}

void IJK_One_Dimensional_Subproblem::retrieve_temperature_solution()
{
  temperature_solution_.resize(rhs_assembly_.size());
  for (int i=start_index_; i<end_index_; i++)
    temperature_solution_[i - start_index_] = (*thermal_subproblems_temperature_solution_)[i];
}

void IJK_One_Dimensional_Subproblem::compute_local_temperature_gradient_solution()
{
  //	for (int dir=0; dir<3; dir++)
  //		temperature_gradient_solution_[dir].resize(temperature_solution_.size());
  normal_temperature_gradient_solution_.resize(temperature_solution_.size());
  (*finite_difference_assembler_).compute_operator(radial_first_order_operator_base_, temperature_solution_, normal_temperature_gradient_solution_);

  normal_temperature_double_derivative_solution_.resize(temperature_solution_.size());
  (*finite_difference_assembler_).compute_operator(radial_second_order_operator_base_, temperature_solution_, normal_temperature_double_derivative_solution_);

  temperature_x_gradient_solution_.resize(normal_temperature_gradient_solution_.size());
  project_basis_vector_onto_cartesian_dir(0, tangential_temperature_gradient_first_, tangential_temperature_gradient_second_, normal_temperature_gradient_solution_,
                                          *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                          temperature_x_gradient_solution_);
  temperature_y_gradient_solution_.resize(normal_temperature_gradient_solution_.size());
  project_basis_vector_onto_cartesian_dir(1, tangential_temperature_gradient_first_, tangential_temperature_gradient_second_, normal_temperature_gradient_solution_,
                                          *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                          temperature_y_gradient_solution_);
  temperature_z_gradient_solution_.resize(normal_temperature_gradient_solution_.size());
  project_basis_vector_onto_cartesian_dir(2, tangential_temperature_gradient_first_, tangential_temperature_gradient_second_, normal_temperature_gradient_solution_,
                                          *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                          temperature_z_gradient_solution_);

  if (source_terms_type_ == linear_diffusion && source_terms_type_ == spherical_diffusion)
    {
      tangential_temperature_gradient_first_.resize(*points_per_thermal_subproblem_);
      tangential_temperature_gradient_second_.resize(*points_per_thermal_subproblem_);
    }

  thermal_flux_ = normal_temperature_gradient_solution_;
  thermal_flux_*= ((*lambda_) * surface_);
}

double IJK_One_Dimensional_Subproblem::get_interfacial_gradient_corrected() const
{
  return normal_temperature_gradient_solution_[0];
}

void IJK_One_Dimensional_Subproblem::compute_local_velocity_gradient()
{
  normal_velocity_normal_gradient_.resize(*points_per_thermal_subproblem_);
  tangential_velocity_normal_gradient_.resize(*points_per_thermal_subproblem_);
  second_tangential_velocity_normal_gradient_.resize(*points_per_thermal_subproblem_);
  azymuthal_velocity_normal_gradient_.resize(*points_per_thermal_subproblem_);

  if (source_terms_type_ != linear_diffusion && source_terms_type_ != spherical_diffusion)
    {
      (*finite_difference_assembler_).compute_operator(radial_first_order_operator_base_,
                                                       radial_velocity_,
                                                       normal_velocity_normal_gradient_);
      (*finite_difference_assembler_).compute_operator(radial_first_order_operator_base_,
                                                       first_tangential_velocity_,
                                                       tangential_velocity_normal_gradient_);
      (*finite_difference_assembler_).compute_operator(radial_first_order_operator_base_,
                                                       second_tangential_velocity_,
                                                       second_tangential_velocity_normal_gradient_);
      (*finite_difference_assembler_).compute_operator(radial_first_order_operator_base_,
                                                       azymuthal_velocity_,
                                                       azymuthal_velocity_normal_gradient_);
    }
}

double IJK_One_Dimensional_Subproblem::get_normal_velocity_normal_gradient() const
{
  return normal_velocity_normal_gradient_[0];
}

double IJK_One_Dimensional_Subproblem::get_tangential_velocity_normal_gradient() const
{
  return tangential_velocity_normal_gradient_[0];
}

double IJK_One_Dimensional_Subproblem::get_second_tangential_velocity_normal_gradient() const
{
  return second_tangential_velocity_normal_gradient_[0];
}

double IJK_One_Dimensional_Subproblem::get_azymuthal_velocity_normal_gradient() const
{
  return azymuthal_velocity_normal_gradient_[0];
}

double IJK_One_Dimensional_Subproblem::get_field_profile_at_point(const double& dist, const DoubleVect& field) const
{
  double field_value = INVALID_TEMPERATURE;// temperature_solution_;
  if (dist >= (*radial_coordinates_)[0] && dist <= (*radial_coordinates_)[*points_per_thermal_subproblem_-1])
    {
      /*
       * Dummy dichotomy and linear interpolation along the probe
       */
      int left_interval = 0;
      int right_interval = *points_per_thermal_subproblem_-1;
      find_interval(dist, left_interval, right_interval);
      const double field_interp = (field[right_interval] - field[left_interval])
                                  / ((*radial_coordinates_)[right_interval] - (*radial_coordinates_)[left_interval]) *
                                  (dist-(*radial_coordinates_)[left_interval]) + field[left_interval];
      field_value = field_interp;
    }
  if (dist < (*radial_coordinates_)[0])
    {
      const double interfacial_temperature_gradient_solution = get_interfacial_gradient_corrected();
      field_value = interfacial_temperature_gradient_solution * dist;
    }
  return field_value;
}

double IJK_One_Dimensional_Subproblem::get_temperature_profile_at_point(const double& dist) const
{
  return get_field_profile_at_point(dist, temperature_solution_);
}

double IJK_One_Dimensional_Subproblem::get_temperature_gradient_profile_at_point(const double& dist, const int& dir) const
{
  double temperature_gradient = 0;
  switch(dir)
    {
    case 0:
      temperature_gradient = get_field_profile_at_point(dist, temperature_x_gradient_solution_);
      break;
    case 1:
      temperature_gradient = get_field_profile_at_point(dist, temperature_y_gradient_solution_);
      break;
    case 2:
      temperature_gradient = get_field_profile_at_point(dist, temperature_z_gradient_solution_);
      break;
    default:
      temperature_gradient = get_field_profile_at_point(dist, temperature_x_gradient_solution_);
      break;
    }
  return temperature_gradient;
}

void IJK_One_Dimensional_Subproblem::find_interval(const double& dist, int& left_interval, int& right_interval) const
{
  int mid_interval = left_interval + (right_interval - left_interval) / 2;
  while ((right_interval - left_interval) != 1)
    {
      if (dist > (*radial_coordinates_)[mid_interval])
        left_interval = mid_interval;
      else
        right_interval = mid_interval;
      mid_interval = left_interval + (right_interval - left_interval) / 2;
    }
}

void IJK_One_Dimensional_Subproblem::thermal_subresolution_outputs()
{
  post_process_interfacial_quantities();
  post_process_radial_quantities();
}

void IJK_One_Dimensional_Subproblem::post_process_interfacial_quantities()
{
  if (Process::je_suis_maitre())
    {
      const int reset = (!ref_ijk_ft_->get_reprise()) && (ref_ijk_ft_->get_tstep()==0);
      Nom probe_name = Nom("thermal_subproblems_interfacial_quantities") + Nom(".out");
      Nom probe_header = Nom("tstep\ttime\tsubproblem\ttemperature_gradient\ttemperature_double_deriv"
                             "ttemperature_gradient_tangential\ttemperature_gradient_tangential2\ttemperature_gradient_azymuthal"
                             "\tsurface\tthermal_flux\tlambda\talpha\tprandtl_liq"
                             "\tu_x\tu_y\tu_z\tu_r\tu_r_corr\tu_theta\tu_theta2\tu_phi\tdu_r_dr\tdu_theta_dr\tdu_theta2_dr\tdu_phi_dr");
      SFichier fic = Ouvrir_fichier(probe_name, probe_header, reset);
      fic << ref_ijk_ft_->get_tstep() << " " << ref_ijk_ft_->get_current_time() << " ";
      fic << sub_problem_index_;
      fic << normal_temperature_gradient_solution_[0] << normal_temperature_double_derivative_solution_[0];
      fic << tangential_temperature_gradient_first_[0];
      fic << tangential_temperature_gradient_second_[0] << azymuthal_temperature_gradient_[0];
      fic << surface_ << thermal_flux_[0] << *lambda_ << *alpha_ << Pr_l_;
      fic << x_velocity_[0] << y_velocity_[0] << z_velocity_[0];
      fic << radial_velocity_[0] << radial_velocity_corrected_[0];
      fic << first_tangential_velocity_[0];
      fic << second_tangential_velocity_[0] << azymuthal_velocity_[0];
      fic << normal_velocity_normal_gradient_[0] << tangential_velocity_normal_gradient_[0];
      fic << second_tangential_velocity_normal_gradient_[0] << azymuthal_velocity_normal_gradient_[0];
      fic << finl;
      fic.close();
    }
}

void IJK_One_Dimensional_Subproblem::post_process_radial_quantities()
{
  if (Process::je_suis_maitre())
    {
      const int reset = (!ref_ijk_ft_->get_reprise()) && (ref_ijk_ft_->get_tstep()==0);
      Nom probe_name = Nom("thermal_subproblem_") + Nom(sub_problem_index_) + ("_radial_quantities_time_index_")
                       + Nom(ref_ijk_ft_->get_tstep()) + Nom(".out");
      Nom probe_header = Nom("tstep\ttime\tsubproblem\tradial_coord\ttemperature_interp\ttemperature_solution"
                             "\ttemperature_gradient\ttemperature_double_deriv"
                             "\tsurface\tthermal_flux\tlambda\talpha\tprandtl_liq"
                             "\tu_x\tu_y\tu_z\tu_r\tu_r_corr\tu_theta\tu_theta2\tu_phi\tdu_r_dr\tdu_theta_dr\tdu_theta2_dr\tdu_phi_dr");
      SFichier fic = Ouvrir_fichier(probe_name, probe_header, reset);
      for (int i=0; i<(*points_per_thermal_subproblem_); i++)
        {
          fic << ref_ijk_ft_->get_tstep() << " " << ref_ijk_ft_->get_current_time() << " ";
          fic << sub_problem_index_;
          fic << (*radial_coordinates_)[i] << temperature_interp_[i] << temperature_solution_[i];
          fic << normal_temperature_gradient_[i] << normal_temperature_double_derivative_solution_[i];
          fic << surface_ << thermal_flux_[i] << *lambda_ << *alpha_ << Pr_l_;
          fic << x_velocity_[i] << y_velocity_[i] << z_velocity_[i];
          fic << radial_velocity_[i] << radial_velocity_corrected_[i];
          fic << first_tangential_velocity_[i];
          fic << second_tangential_velocity_[i] << azymuthal_velocity_[i];
          fic << normal_velocity_normal_gradient_[i] << tangential_velocity_normal_gradient_[i];
          fic << second_tangential_velocity_normal_gradient_[i] << azymuthal_velocity_normal_gradient_[i];
          fic << finl;
        }
      fic.close();
    }
}

//  double temperature_cell_centre = INVALID_TEMPERATURE;// temperature_solution_;
//  if (dist >= (*radial_coordinates_)[0] && dist <= (*radial_coordinates_)[*points_per_thermal_subproblem_])
//    {
//      /*
//       * Dummy dichotomy and linear interpolation along the probe
//       */
//      int left_interval = 0;
//      int right_interval = *points_per_thermal_subproblem_-1;
//      find_interval(dist, left_interval, right_interval);
//      const double temperature_interp = (temperature_solution_[right_interval] - temperature_solution_[left_interval])
//                                        / ((*radial_coordinates_)[right_interval] - (*radial_coordinates_)[left_interval]) *
//                                        (dist-(*radial_coordinates_)[left_interval]) + temperature_solution_[left_interval];
//      temperature_cell_centre = temperature_interp;
//    }
//  if (dist < (*radial_coordinates_)[0])
//    {
//      const double interfacial_temperature_gradient_solution = get_interfacial_gradient_corrected();
//      temperature_cell_centre = interfacial_temperature_gradient_solution * dist;
//    }
//  return temperature_cell_centre;

