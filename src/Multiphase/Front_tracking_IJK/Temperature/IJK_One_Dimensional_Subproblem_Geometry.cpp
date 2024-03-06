/****************************************************************************
* Copyright (c) 2024, CEA
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
// File      : IJK_One_Dimensional_Subproblem_Geometry.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_One_Dimensional_Subproblem_Geometry.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_FT.h>

Implemente_instanciable_sans_constructeur( IJK_One_Dimensional_Subproblem_Geometry, "IJK_One_Dimensional_Subproblem_Geometry", Objet_U ) ;

IJK_One_Dimensional_Subproblem_Geometry::IJK_One_Dimensional_Subproblem_Geometry()
{
  interfaces_ = nullptr;

  points_per_thermal_subproblem_ = nullptr;

  disable_probe_collision_ = 0;
  disable_find_cell_centre_probe_tip_ = 0;
  enable_resize_probe_collision_ = 0;
  resize_probe_collision_index_ = 0;
}

Sortie& IJK_One_Dimensional_Subproblem_Geometry::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& IJK_One_Dimensional_Subproblem_Geometry::readOn( Entree& is )
{
  Objet_U::readOn( is );
  return is;
}

//void IJK_One_Dimensional_Subproblem_Geometry::compute_interface_basis_vectors()
//{
//  /*
//   * TODO: Associate a basis to each subproblem
//   * Use Rodrigues' rotation formula to determine ephi ?
//   * Needs an axis of (rotation gravity_dir x relative_vectors)
//   * and an angle (gravity_dir dot relative_vectors) / (norm(gravity_dir)*norm(relative_vectors))
//   * ephi is determined in the gravity_align rising direction
//   * 		 | gravity_dir
//   * 		 |
//   *   *****
//   * ***   ***
//   * **     **
//   * ***   ***
//   *   *****
//   *     |
//   *     |
//   */
//
//  facet_barycentre_relative_ = facet_barycentre_ - bubble_barycentre_;
//  if (debug_)
//    {
//      Cerr << "bubble_barycentre_"<< bubble_barycentre_[0] << " ; " << bubble_barycentre_[1] << " ; " << bubble_barycentre_[2] << finl;
//      Cerr << "facet_barycentre_"<< facet_barycentre_[0] << " ; " << facet_barycentre_[1] << " ; " << facet_barycentre_[2] << finl;
//      Cerr << "facet_barycentre_relative_"<< facet_barycentre_relative_[0] << " ; " << facet_barycentre_relative_[1] << " ; " << facet_barycentre_relative_[2] << finl;
//    }
//  Vecteur3 facet_barycentre_relative_normed = facet_barycentre_relative_;
//  const double facet_barycentre_relative_norm = facet_barycentre_relative_normed.length();
//  facet_barycentre_relative_normed *= (1 / facet_barycentre_relative_norm);
//  Vecteur3 normal_contrib;
//  const double normal_vector_compo_norm = normal_vector_compo_.length();
//  normal_vector_compo_ *= (1 / normal_vector_compo_norm);
//
//  if (debug_)
//    Cerr << "Normal vector norm:" << normal_vector_compo_norm << finl;
//  /*
//   * First method with tangential direction of maximum velocity variations
//   */
//  DoubleTab facet_barycentre(1, 3);
//  interfacial_velocity_compo_ = 0.;
//  for (int dir=0; dir<3; dir++)
//    facet_barycentre(0, dir) = facet_barycentre_[dir];
//  for (int dir=0; dir<3; dir++)
//    {
//      DoubleVect interfacial_velocity_component(1);
//      ijk_interpolate_skip_unknown_points((*velocity_)[dir], facet_barycentre, interfacial_velocity_component, INVALID_INTERP);
//      interfacial_velocity_compo_[dir] = interfacial_velocity_component[0];
//    }
//  if (interfacial_velocity_compo_.length() < INVALID_VELOCITY)
//    {
//      normal_contrib = normal_vector_compo_;
//      normal_contrib *= Vecteur3::produit_scalaire(facet_barycentre_relative_normed, normal_vector_compo_);
//      first_tangential_vector_compo_ = facet_barycentre_relative_normed - normal_contrib;
//    }
//  else
//    {
//      // Should I remove the rising velocity ?
//      interfacial_velocity_compo_ = interfacial_velocity_compo_ - bubble_rising_velocity_compo_;
//      normal_contrib = normal_vector_compo_;
//      normal_contrib *= Vecteur3::produit_scalaire(interfacial_velocity_compo_, normal_vector_compo_);
//      interfacial_tangential_velocity_compo_ = interfacial_velocity_compo_ - normal_contrib;
//      first_tangential_vector_compo_ = interfacial_tangential_velocity_compo_;
//    }
//  const double norm_first_tangential_vector = first_tangential_vector_compo_.length();
//  first_tangential_vector_compo_ *= (1 / norm_first_tangential_vector);
//  Vecteur3::produit_vectoriel(normal_vector_compo_, first_tangential_vector_compo_, second_tangential_vector_compo_);
//  const double norm_second_tangential_vector = second_tangential_vector_compo_.length();
//  second_tangential_vector_compo_ *= (1 / norm_second_tangential_vector);
//
//  /*
//   * Second method with rising velocity
//   */
//  Vecteur3::produit_vectoriel(bubble_rising_vector_, facet_barycentre_relative_, azymuthal_vector_compo_raw_);
//
//  azymuthal_vector_compo_ = azymuthal_vector_compo_raw_;
//  const double norm_azymuthal_vector_compo_raw_ = azymuthal_vector_compo_raw_.length();
////  const int sign_vector = signbit(Vecteur3::produit_scalaire(bubble_rising_vector_, normal_vector_compo_));
////  if (sign_vector)
////    azymuthal_vector_compo_ *= -1;
//  azymuthal_vector_compo_ *= (1 / norm_azymuthal_vector_compo_raw_);
//
//  normal_contrib = normal_vector_compo_;
//  normal_contrib *=	Vecteur3::produit_scalaire(azymuthal_vector_compo_, normal_vector_compo_);
//  azymuthal_vector_compo_ = azymuthal_vector_compo_ - normal_contrib;
//  Vecteur3::produit_vectoriel(azymuthal_vector_compo_, normal_vector_compo_, first_tangential_vector_compo_from_rising_dir_);
//  const double norm_first_tangential_vector_from_rising_dir = first_tangential_vector_compo_from_rising_dir_.length();
//  first_tangential_vector_compo_from_rising_dir_ *= (1 / norm_first_tangential_vector_from_rising_dir);
//
//
//  if (tangential_from_rising_vel_)
//    {
//      first_tangential_vector_compo_solver_ = &first_tangential_vector_compo_from_rising_dir_;
//      second_tangential_vector_compo_solver_ = &azymuthal_vector_compo_;
//      first_tangential_velocity_solver_ = &first_tangential_velocity_from_rising_dir_corrected_;
//      second_tangential_velocity_solver_ = &azymuthal_velocity_corrected_;
//      tangential_temperature_gradient_first_solver_ = &tangential_temperature_gradient_first_from_rising_dir_;
//      tangential_temperature_gradient_second_solver_ = &azymuthal_temperature_gradient_;
//    }
//  else
//    {
//      // By default
//      first_tangential_vector_compo_solver_= &first_tangential_vector_compo_;
//      second_tangential_vector_compo_solver_ = &second_tangential_vector_compo_;
//      first_tangential_velocity_solver_ = &first_tangential_velocity_corrected_;
//      second_tangential_velocity_solver_ = &second_tangential_velocity_corrected_;
//      tangential_temperature_gradient_first_solver_ = &tangential_temperature_gradient_first_;
//      tangential_temperature_gradient_second_solver_ = &tangential_temperature_gradient_second_;
//    }
//}
//
//void IJK_One_Dimensional_Subproblem_Geometry::compute_pure_spherical_basis_vectors()
//{
//  /*
//   * FIXME: It is align with gravity z but it should be modified to be align with the gravity dir ?
//   */
//  if (debug_)
//    Cerr << "r_sph_ calculation"  << finl;
//  r_sph_ = sqrt(facet_barycentre_relative_[0] * facet_barycentre_relative_[0]
//                + facet_barycentre_relative_[1] * facet_barycentre_relative_[1]
//                + facet_barycentre_relative_[2] * facet_barycentre_relative_[2]);
//  if (debug_)
//    {
//      Cerr << "r_sph_ = " << r_sph_ << finl;
//      Cerr << "theta_sph_ calculation"  << finl;
//    }
//
////  theta_sph_ = atan(sqrt(facet_barycentre_relative_[0] * facet_barycentre_relative_[0]
////                         + facet_barycentre_relative_[1] * facet_barycentre_relative_[1])/ facet_barycentre_relative_[2]);
//  theta_sph_ = atan2(sqrt(facet_barycentre_relative_[0] * facet_barycentre_relative_[0]
//                          + facet_barycentre_relative_[1] * facet_barycentre_relative_[1]), facet_barycentre_relative_[2]);
//  const double atan_theta_incr_ini = M_PI / 2;
//  const double atan_incr_factor = -1;
//  theta_sph_ = (theta_sph_ - atan_theta_incr_ini) * atan_incr_factor;
//
//  if (debug_)
//    {
//      Cerr << "theta_sph_ = " << theta_sph_ << finl;
//      Cerr << "phi_sph_ calculation"  << finl;
//    }
//  phi_sph_ = atan2(facet_barycentre_relative_[1], facet_barycentre_relative_[0]);
//
//  if (debug_)
//    {
//      Cerr << "phi_sph_ = " << phi_sph_ << finl;
//      Cerr << "er_sph_ calculation"  << finl;
//    }
//  for (int dir=0; dir<3; dir++)
//    er_sph_[dir] = facet_barycentre_relative_[dir] / r_sph_;
//
//  const double length = sqrt(facet_barycentre_relative_[0] * facet_barycentre_relative_[0]
//                             + facet_barycentre_relative_[1] * facet_barycentre_relative_[1]);
//
//  if (debug_)
//    {
//      Cerr << "er_sph_ = " << er_sph_[0] << finl;
//      Cerr << "etheta_sph_ calculation"  << finl;
//    }
//  for (int dir=0; dir<2; dir++)
//    etheta_sph_[dir] = facet_barycentre_relative_[dir] * facet_barycentre_relative_[2] / (r_sph_ * length);
//  etheta_sph_[2] = - facet_barycentre_relative_[2] * length / r_sph_;
//
//  ephi_sph_ = {0., 0., 0.};
//  ephi_sph_[0] = - facet_barycentre_relative_[1];
//  ephi_sph_[1] = facet_barycentre_relative_[0];
//}
//
//void IJK_One_Dimensional_Subproblem::compute_local_discretisation()
//{
//  int i;
//  if (global_probes_characteristics_)
//    {
//      if (!probe_variations_enabled_)
//        {
//          if (!velocities_calculation_counter_)
//            {
//              radial_coordinates_ = radial_coordinates_base_;
//              dr_ = *dr_base_;
//            }
//        }
//      else
//        {
//          dr_ = probe_length_ / (*points_per_thermal_subproblem_ - 1);
//          radial_coordinates_modified_.resize(*points_per_thermal_subproblem_);
//          for (i=0; i < *points_per_thermal_subproblem_; i++)
//            radial_coordinates_modified_(i) = i * dr_;
//          radial_coordinates_ = &radial_coordinates_modified_;
//        }
//    }
//  else
//    {
//      /*
//       * coeff_distance_diagonal_ as well as
//       * points_per_thermal_subproblem_ could be adapted
//       */
//      if (!probe_variations_enabled_)
//        radial_coordinates_modified_.resize(*points_per_thermal_subproblem_);
//      dr_ = probe_length_ / (*points_per_thermal_subproblem_ - 1);
//      for (i=0; i < *points_per_thermal_subproblem_; i++)
//        radial_coordinates_modified_(i) = i * dr_;
//      radial_coordinates_ = &radial_coordinates_modified_;
//    }
//  /*
//   * Following attributes differ anyway !
//   */
//  if (!velocities_calculation_counter_ || probe_variations_enabled_)
//    {
//      dr_inv_ = 1 / dr_;
//      osculating_radial_coordinates_ = (*radial_coordinates_);
//      osculating_radial_coordinates_ += osculating_radius_;
//      if (!probe_variations_enabled_)
//        {
//          radial_coordinates_cartesian_compo_.resize(*points_per_thermal_subproblem_, 3);
//          coordinates_cartesian_compo_.resize(*points_per_thermal_subproblem_, 3);
//          osculating_radial_coordinates_cartesian_compo_.resize(*points_per_thermal_subproblem_, 3);
//          osculating_radial_coordinates_inv_.resize(*points_per_thermal_subproblem_);
//        }
//      for (i=0; i < *points_per_thermal_subproblem_; i++)
//        {
//          osculating_radial_coordinates_inv_[i] = 1 / osculating_radial_coordinates_[i];
//          for (int dir=0; dir<3; dir++)
//            {
//              radial_coordinates_cartesian_compo_(i, dir) = (*radial_coordinates_)(i) * normal_vector_compo_[dir];
//              osculating_radial_coordinates_cartesian_compo_(i, dir) = osculating_radial_coordinates_(i) * normal_vector_compo_[dir];
//              coordinates_cartesian_compo_(i, dir) = radial_coordinates_cartesian_compo_(i, dir) + facet_barycentre_[dir];
//            }
//        }
//    }
//}

void IJK_One_Dimensional_Subproblem_Geometry::interpolate_indicator_on_probes()
{
  indicator_interp_.resize(*points_per_thermal_subproblem_);
  const IJK_Field_double& indicator = interfaces_->I();
  ijk_interpolate_skip_unknown_points(indicator, coordinates_cartesian_compo_, indicator_interp_, INVALID_INTERP);
  if (enable_resize_probe_collision_)
    {
      for (int i=(*points_per_thermal_subproblem_)-1; i>=0; i--)
        {
          const double indic_last = find_cell_related_indicator_on_probes(i);
          if (indic_last > LIQUID_INDICATOR_TEST)
            {
              resize_probe_collision_index_ = i;
              return;
            }
          disable_probe_collision_ = (i == 0);
        }
    }
  else
    {
      if (!disable_find_cell_centre_probe_tip_)
        {

          const int last_index = (*points_per_thermal_subproblem_) - 1;
          const double indic_last = find_cell_related_indicator_on_probes(last_index);
          if (indic_last < LIQUID_INDICATOR_TEST)
            {
              disable_probe_collision_ = 1;
              return;
            }
        }
      else
        {
          for (int i=(*points_per_thermal_subproblem_)-1; i>=0; i--)
            {
              const double indicator_val = indicator_interp_(i);
              if (indicator_val < LIQUID_INDICATOR_TEST)
                {
                  disable_probe_collision_ = 1;
                  return;
                }
            }
        }
    }
}

double IJK_One_Dimensional_Subproblem_Geometry::find_cell_related_indicator_on_probes(const int& last_index)
{
  const IJK_Field_double& indicator = interfaces_->I();
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();

  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  Vecteur3 xyz_cart_end = {coordinates_cartesian_compo_(last_index, 0),
                           coordinates_cartesian_compo_(last_index, 1),
                           coordinates_cartesian_compo_(last_index, 2)
                          };
  Vecteur3 centre_elem = ref_ijk_ft_->get_splitting_ns().get_coords_of_dof(index_i_, index_j_, index_k_, IJK_Splitting::ELEM);
  Vecteur3 displacement_centre_probe = centre_elem;
  displacement_centre_probe *= (-1);
  displacement_centre_probe += xyz_cart_end;
  Vecteur3 displacement_factor = {displacement_centre_probe[0] / dx,
                                  displacement_centre_probe[1] / dy,
                                  displacement_centre_probe[2] / dz
                                 };
  const int offset_x = (int) displacement_factor[0];
  const int offset_y = (int) displacement_factor[1];
  const int offset_z = (int) displacement_factor[2];
  const int offset_x_elem = ((abs(displacement_factor[0] - offset_x) >= 0.5) ? 1 : 0);
  const int offset_y_elem = ((abs(displacement_factor[1] - offset_y) >= 0.5) ? 1 : 0);
  const int offset_z_elem = ((abs(displacement_factor[2] - offset_z) >= 0.5) ? 1 : 0);
  const int real_offset_x = offset_x + (signbit(offset_x) ? - offset_x_elem : offset_x_elem);
  const int real_offset_y = offset_y + (signbit(offset_y) ? - offset_y_elem : offset_y_elem);
  const int real_offset_z = offset_z + (signbit(offset_z) ? - offset_z_elem : offset_z_elem);
  const double indic_last = indicator(index_i_ + real_offset_x,
                                      index_j_ + real_offset_y,
                                      index_k_ + real_offset_z);
  return indic_last;
}

void IJK_One_Dimensional_Subproblem_Geometry::compute_distance_cell_centre()
{
  if (!has_computed_cell_centre_distance_)
    {
      Vecteur3 centre = ref_ijk_ft_->get_splitting_ns().get_coords_of_dof(index_i_, index_j_, index_k_, IJK_Splitting::ELEM);

      Vecteur3 facet_to_cell_centre = facet_barycentre_;
      facet_to_cell_centre *= -1;
      facet_to_cell_centre += centre;
      cell_centre_distance_ = Vecteur3::produit_scalaire(facet_to_cell_centre, normal_vector_compo_);

      Vecteur3 normal_contrib = normal_vector_compo_;
      normal_contrib *= cell_centre_distance_;
      Vecteur3 tangential_displacement = normal_contrib;
      tangential_displacement *= (-1);
      tangential_displacement += facet_to_cell_centre;
      cell_centre_tangential_distance_ = tangential_displacement.length();
      tangential_distance_vector_ = tangential_displacement;
      if (cell_centre_tangential_distance_ > 1e-16)
        tangential_distance_vector_ *= (1 / cell_centre_tangential_distance_);
      has_computed_cell_centre_distance_ = true;
    }
  else if (debug_)
    Cerr << "Cell centre distances have already been computed" << finl;
}

void IJK_One_Dimensional_Subproblem_Geometry::compute_distance_faces_centres()
{
  if (!has_computed_cell_faces_distance_)
    {
      Vecteur3 bary_face {0., 0., .0};
      Vecteur3 vector_relative {0., 0., 0.};
      Vecteur3 bary_vertex {0., 0., 0.};
      int neighbours_i[6] = NEIGHBOURS_I;
      int neighbours_j[6] = NEIGHBOURS_J;
      int neighbours_k[6] = NEIGHBOURS_K;
      int neighbours_faces_i[6] = NEIGHBOURS_FACES_I;
      int neighbours_faces_j[6] = NEIGHBOURS_FACES_J;
      int neighbours_faces_k[6] = NEIGHBOURS_FACES_K;
      int face_dir[6] = FACES_DIR;
      int m;
      for (int l=0; l<6; l++)
        {
          const int ii = neighbours_i[l];
          const int jj = neighbours_j[l];
          const int kk = neighbours_k[l];
          const double indic_neighbour = ref_ijk_ft_->itfce().I()(index_i_+ii, index_j_+jj, index_k_+kk);
          if (fabs(indic_neighbour) > LIQUID_INDICATOR_TEST)
            {
              const int ii_f = neighbours_faces_i[l];
              const int jj_f = neighbours_faces_j[l];
              const int kk_f = neighbours_faces_k[l];
              pure_liquid_neighbours_[l] = 1;
              if (ii)
                bary_face = ref_ijk_ft_->get_splitting_ns().get_coords_of_dof(index_i_+ii_f, index_j_+jj_f, index_k_+kk_f, IJK_Splitting::FACES_I);
              if (jj)
                bary_face = ref_ijk_ft_->get_splitting_ns().get_coords_of_dof(index_i_+ii_f, index_j_+jj_f, index_k_+kk_f, IJK_Splitting::FACES_J);
              if (kk)
                bary_face = ref_ijk_ft_->get_splitting_ns().get_coords_of_dof(index_i_+ii_f, index_j_+jj_f, index_k_+kk_f, IJK_Splitting::FACES_K);

              // Normal distance
              vector_relative = facet_barycentre_;
              vector_relative *= (-1);
              vector_relative += bary_face;
              {
                const double distance_face_centre = Vecteur3::produit_scalaire(vector_relative, normal_vector_compo_);
                face_centres_distance_[l] = distance_face_centre;
                // Tangential distance
                Vecteur3 normal_contrib = normal_vector_compo_;
                normal_contrib *= distance_face_centre;
                Vecteur3 tangential_displacement = normal_contrib;
                tangential_displacement *= (-1);
                tangential_displacement += vector_relative;
                face_centres_tangential_distance_[l] = tangential_displacement.length();
                if (face_centres_tangential_distance_[l] > 1e-16)
                  tangential_displacement *= (1 / face_centres_tangential_distance_[l]);
                face_tangential_distance_vector_[l] = tangential_displacement;
              }
              // Distance to vertex
              for (m=0; m<4; m++)
                {
                  double distance_vertex_centre = 0.;
                  double tangential_distance_vertex_centre = 0.;
                  Vecteur3 tangential_distance_vector_vertex_centre = {0., 0., 0.};
                  bary_vertex = vector_relative;
                  compute_vertex_position(m,
                                          face_dir[l],
                                          bary_vertex,
                                          distance_vertex_centre,
                                          tangential_distance_vertex_centre,
                                          tangential_distance_vector_vertex_centre);
                  vertices_centres_distance_[l][m] = distance_vertex_centre;
                  vertices_centres_tangential_distance_[l][m] = tangential_distance_vertex_centre;
                  vertices_tangential_distance_vector_[l][m] = tangential_distance_vector_vertex_centre;
                }
            }
          else
            {
              pure_liquid_neighbours_[l] = 0;
              face_centres_distance_[l] = 0.;
              face_centres_tangential_distance_[l] = 0.;
              face_tangential_distance_vector_[l] = {0., 0., 0.};
              for (m=0; m<4; m++)
                {
                  vertices_centres_distance_[l][m] = 0.;
                  vertices_centres_tangential_distance_[l][m] = 0.;
                  vertices_tangential_distance_vector_[l][m] = {0., 0., 0.};
                }
            }
        }
      has_computed_cell_faces_distance_ = true;
    }
  else if (debug_)
    Cerr << "Cell face distances have already been computed" << finl;
}

void IJK_One_Dimensional_Subproblem_Geometry::compute_vertex_position(const int& vertex_number,
                                                                      const int& face_dir,
                                                                      Vecteur3& bary_vertex,
                                                                      double& distance_vertex_centre,
                                                                      double& tangential_distance_vertex_centre,
                                                                      Vecteur3& tangential_distance_vector_vertex_centre)
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  const double neighbours_first_dir[4] = NEIGHBOURS_FIRST_DIR;
  const double neighbours_second_dir[4] = NEIGHBOURS_SECOND_DIR;
  Vecteur3 point_coords {0., 0., 0.};
  double dl1;
  double dl2;
  switch(face_dir)
    {
    case 0:
      dl1 = dy / 2.;
      dl2 = dz / 2.;
      point_coords[1] = dl1 * neighbours_first_dir[vertex_number];
      point_coords[2] = dl2 * neighbours_second_dir[vertex_number];
      break;
    case 1:
      dl1 = dx / 2.;
      dl2 = dz / 2.;
      point_coords[0] = dl1 * neighbours_first_dir[vertex_number];
      point_coords[2] = dl2 * neighbours_second_dir[vertex_number];
      break;
    case 2:
      dl1 = dx / 2.;
      dl2 = dy / 2.;
      point_coords[0] = dl1 * neighbours_first_dir[vertex_number];
      point_coords[1] = dl2 * neighbours_second_dir[vertex_number];
      break;
    default:
      dl1 = dx / 2.;
      dl2 = dy / 2.;
      point_coords[0] = dl1 * neighbours_first_dir[vertex_number];
      point_coords[1] = dl2 * neighbours_second_dir[vertex_number];
      break;
    }
  bary_vertex += point_coords;
  distance_vertex_centre = Vecteur3::produit_scalaire(bary_vertex, normal_vector_compo_);
  Vecteur3 tangential_distance_vector = normal_vector_compo_;
  tangential_distance_vector *= distance_vertex_centre;
  tangential_distance_vector *= (-1);
  tangential_distance_vector += bary_vertex;
  tangential_distance_vertex_centre = tangential_distance_vector.length();
  if (tangential_distance_vertex_centre > 1e-16)
    tangential_distance_vector *= (1 / tangential_distance_vertex_centre);
  tangential_distance_vector_vertex_centre = tangential_distance_vector;
}

void IJK_One_Dimensional_Subproblem_Geometry::compute_distance_cell_centres_neighbours()
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  int l, m, n;

  int dxyz_increment_max = get_dxyz_increment_max();
  /*
   * 8-1 values for one neighbour in each dir... (Too much, enhance later)
   * Positive OR Negative dir depending on the normal vector
   */
  pure_neighbours_to_correct_.resize(dxyz_increment_max + 1);
  pure_neighbours_corrected_distance_.resize(dxyz_increment_max + 1);
  if (neighbours_weighting_)
    pure_neighbours_corrected_colinearity_.resize(dxyz_increment_max + 1);
  for (l=dxyz_increment_max; l>=0; l--)
    {
      pure_neighbours_to_correct_[l].resize(dxyz_increment_max + 1);
      pure_neighbours_corrected_distance_[l].resize(dxyz_increment_max + 1);
      if (neighbours_weighting_)
        pure_neighbours_corrected_colinearity_[l].resize(dxyz_increment_max + 1);
      for (m=dxyz_increment_max; m>=0; m--)
        {
          pure_neighbours_to_correct_[l][m].resize(dxyz_increment_max + 1);
          pure_neighbours_corrected_distance_[l][m].resize(dxyz_increment_max + 1);
          if (neighbours_weighting_)
            pure_neighbours_corrected_colinearity_[l][m].resize(dxyz_increment_max + 1);
          for (n=dxyz_increment_max; n>=0; n--)
            {
              pure_neighbours_to_correct_[l][m][n] = false;
              pure_neighbours_corrected_distance_[l][m][n] = 0.;
              if (neighbours_weighting_)
                pure_neighbours_corrected_colinearity_[l][m][n] = 0.;
            }
        }
    }

  // get_maximum_remaining_distance()
  double remaining_distance_diag = probe_length_ - cell_centre_distance_;
  Vecteur3 remaining_distance_diag_projected = normal_vector_compo_;
  remaining_distance_diag_projected *= remaining_distance_diag;
  for (int i=0; i<3; i++)
    pure_neighbours_corrected_sign_[i] = signbit(normal_vector_compo_[i]);
  int dx_increment = (int) abs(remaining_distance_diag_projected[0] / dx);
  int dy_increment = (int) abs(remaining_distance_diag_projected[1] / dy);
  int dz_increment = (int) abs(remaining_distance_diag_projected[2] / dz);
  if (correct_neighbours_rank_)
    {
      dx_increment = std::min(dx_increment, dxyz_increment_max);
      dy_increment = std::min(dy_increment, dxyz_increment_max);
      dz_increment = std::min(dz_increment, dxyz_increment_max);
    }
  dxyz_increment_bool_ = (dx_increment || dy_increment || dz_increment);
  for (l=dx_increment; l>=0; l--)
    for (m=dy_increment; m>=0; m--)
      for (n=dz_increment; n>=0; n--)
        if (l!=0 || m!=0 || n!=0)
          {
            const int l_dir = (pure_neighbours_corrected_sign_[0]) ? l * (-1) : l;
            const int m_dir = (pure_neighbours_corrected_sign_[1]) ? m * (-1) : m;
            const int n_dir = (pure_neighbours_corrected_sign_[2]) ? n * (-1) : n;
            const double indic_neighbour = ref_ijk_ft_->itfce().I()(index_i_ + l_dir, index_j_ + m_dir, index_k_ + n_dir);
            if (indic_neighbour > LIQUID_INDICATOR_TEST)
              {
                pure_neighbours_to_correct_[l][m][n] = true;
                const double dx_contrib = l_dir * dx;
                const double dy_contrib = m_dir * dy;
                const double dz_contrib = n_dir * dz;
                Vecteur3 distance_contrib = {dx_contrib, dy_contrib, dz_contrib};
                pure_neighbours_corrected_distance_[l][m][n] = cell_centre_distance_
                                                               + Vecteur3::produit_scalaire(normal_vector_compo_, distance_contrib);
                if (neighbours_weighting_)
                  {
                    const double colinearity = compute_cell_weighting(dx_contrib, dy_contrib, dz_contrib);
                    pure_neighbours_corrected_colinearity_[l][m][n] = colinearity;
                  }
              }
          }
}

double IJK_One_Dimensional_Subproblem_Geometry::compute_cell_weighting(const double& dx_contrib,
                                                                       const double& dy_contrib,
                                                                       const double& dz_contrib)
{
  if (neighbours_colinearity_weighting_)
    return compute_colinearity(dx_contrib, dy_contrib, dz_contrib);
  if (neighbours_distance_weighting_)
    return compute_distance_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  if (neighbours_colinearity_distance_weighting_)
    return compute_colinearity(dx_contrib, dy_contrib, dz_contrib) * compute_distance_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  return 1;
}

void IJK_One_Dimensional_Subproblem_Geometry::compute_distance_last_cell_faces_neighbours()
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  const double dx_over_two = dx / 2.;
  const double dy_over_two = dy / 2.;
  const double dz_over_two = dz / 2.;
  int l, m, n;
  int l_cell, m_cell, n_cell;

  int dxyz_increment_max = get_dxyz_increment_max();
  int dxyz_over_two_increment_max = get_dxyz_over_two_increment_max();
  const int first_increment[3] = {dxyz_over_two_increment_max + 1, dxyz_increment_max, dxyz_increment_max};
  const int second_increment[3] = {dxyz_increment_max, dxyz_over_two_increment_max + 1, dxyz_increment_max};
  const int third_increment[3] = {dxyz_increment_max, dxyz_increment_max, dxyz_over_two_increment_max + 1};

  /*
   * 8-1 values for one neighbour in each dir... (Too much, enhance later)
   * Positive OR Negative dir depending on the normal vector
   */
  pure_neighbours_last_faces_to_correct_.resize(3);
  pure_neighbours_last_faces_corrected_distance_.resize(3);
  if (neighbours_last_faces_weighting_)
    pure_neighbours_last_faces_corrected_colinearity_.resize(3);
  for (int c=0; c<3; c++)
    {
      const int first_incr = first_increment[c];
      const int second_incr = second_increment[c];
      const int third_incr = third_increment[c];
      pure_neighbours_last_faces_to_correct_[c].resize(first_incr + 1);
      pure_neighbours_last_faces_corrected_distance_[c].resize(first_incr + 1);
      if (neighbours_last_faces_weighting_)
        pure_neighbours_last_faces_corrected_colinearity_[c].resize(first_incr + 1);
      for (l=first_incr; l>=0; l--)
        {
          pure_neighbours_last_faces_to_correct_[c][l].resize(second_incr + 1);
          pure_neighbours_last_faces_corrected_distance_[c][l].resize(second_incr + 1);
          if (neighbours_last_faces_weighting_)
            pure_neighbours_last_faces_corrected_colinearity_[c][l].resize(second_incr + 1);
          for (m=second_incr; m>=0; m--)
            {
              pure_neighbours_last_faces_to_correct_[c][l][m].resize(third_incr + 1);
              pure_neighbours_last_faces_corrected_distance_[c][l][m].resize(third_incr + 1);
              if (neighbours_last_faces_weighting_)
                pure_neighbours_last_faces_corrected_colinearity_[c][l][m].resize(third_incr + 1);
              for (n=third_incr; n>=0; n--)
                {
                  pure_neighbours_last_faces_to_correct_[c][l][m][n] = false;
                  pure_neighbours_last_faces_corrected_distance_[c][l][m][n] = 0.;
                  if (neighbours_last_faces_weighting_)
                    pure_neighbours_last_faces_corrected_colinearity_[c][l][m][n] = 0.;
                }
            }
        }
    }

  double remaining_distance_diag = probe_length_ - cell_centre_distance_;
  Vecteur3 remaining_distance_diag_projected = normal_vector_compo_;
  remaining_distance_diag_projected *= remaining_distance_diag;
  for (int i=0; i<3; i++)
    pure_neighbours_corrected_sign_[i] = signbit(normal_vector_compo_[i]);
  int dx_over_two_increment = (int) abs(remaining_distance_diag_projected[0] / dx_over_two);
  int dy_over_two_increment = (int) abs(remaining_distance_diag_projected[1] / dy_over_two);
  int dz_over_two_increment = (int) abs(remaining_distance_diag_projected[2] / dz_over_two);
  int dx_increment = (int) (dx_over_two_increment / 2);
  int dy_increment = (int) (dy_over_two_increment / 2);
  int dz_increment = (int) (dz_over_two_increment / 2);
  dx_over_two_increment -= dx_increment;
  dy_over_two_increment -= dy_increment;
  dz_over_two_increment -= dz_increment;
  if (correct_neighbours_rank_)
    {
      dx_over_two_increment = std::min(dx_over_two_increment, dxyz_over_two_increment_max);
      dy_over_two_increment = std::min(dy_over_two_increment, dxyz_over_two_increment_max);
      dz_over_two_increment = std::min(dz_over_two_increment, dxyz_over_two_increment_max);
      dx_increment = std::min(dx_increment, dxyz_increment_max);
      dy_increment = std::min(dy_increment, dxyz_increment_max);
      dz_increment = std::min(dz_increment, dxyz_increment_max);
    }
  dxyz_over_two_increment_bool_ = (dx_over_two_increment >0 || dy_over_two_increment>0 || dz_over_two_increment >0 );
  if (dx_over_two_increment>0)
    dx_over_two_increment--;
  if (dy_over_two_increment>0)
    dy_over_two_increment--;
  if (dz_over_two_increment>0)
    dz_over_two_increment--;

  /*
   * TODO: Should we look to cells faces that are not in the normal vector direction ?
   */
  for (l=dx_over_two_increment + 1; l>=0; l--)
    for (m_cell=dy_increment; m_cell>=0; m_cell--)
      for (n_cell=dz_increment; n_cell>=0; n_cell--)
        {
          const int l_dir = (pure_neighbours_corrected_sign_[0]) ? l * (-1) + 1 : l;
          const int m_dir = (pure_neighbours_corrected_sign_[1]) ? m_cell * (-1) : m_cell;
          const int n_dir = (pure_neighbours_corrected_sign_[2]) ? n_cell * (-1) : n_cell;
          const int l_dir_elem = (pure_neighbours_corrected_sign_[0]) ? (l + 1) * (-1) + 1 : l;
          const double indic_neighbour = ref_ijk_ft_->itfce().I()(index_i_ + l_dir_elem, index_j_ + m_dir, index_k_ + n_dir);
          if (indic_neighbour > LIQUID_INDICATOR_TEST)
            {
              pure_neighbours_last_faces_to_correct_[0][l][m_cell][n_cell] = true;
              const double lmn_zero = (l > 0) ? 1. : 0.;
              const double contrib_factor = (pure_neighbours_corrected_sign_[0]) ? (lmn_zero * (2 * abs(l_dir) + 1) - (1. - lmn_zero)) * (-1):
                                            lmn_zero * (2 * (l_dir - 1) + 1) - (1. - lmn_zero);
              const double dx_contrib = contrib_factor * dx_over_two;
              const double dy_contrib = m_dir * dy;
              const double dz_contrib = n_dir * dz;
              Vecteur3 distance_contrib = {dx_contrib, dy_contrib, dz_contrib};
              pure_neighbours_last_faces_corrected_distance_[0][l][m_cell][n_cell] = cell_centre_distance_
                                                                                     + Vecteur3::produit_scalaire(normal_vector_compo_, distance_contrib);
              if (pure_neighbours_last_faces_corrected_distance_[0][l][m_cell][n_cell] < 0)
                pure_neighbours_last_faces_to_correct_[0][l][m_cell][n_cell] = false;
              if (neighbours_last_faces_weighting_)
                {
                  const double colinearity = compute_cell_faces_weighting(dx_contrib, dy_contrib, dz_contrib, 0);
                  pure_neighbours_last_faces_corrected_colinearity_[0][l][m_cell][n_cell] = colinearity;
                }
            }
        }

  for (l_cell=dx_increment; l_cell>=0; l_cell--)
    for (m=dy_over_two_increment + 1; m>=0; m--)
      for (n_cell=dz_increment; n_cell>=0; n_cell--)
        {
          const int l_dir = (pure_neighbours_corrected_sign_[0]) ? l_cell * (-1) : l_cell;
          const int m_dir = (pure_neighbours_corrected_sign_[1]) ? m * (-1) + 1: m;
          const int n_dir = (pure_neighbours_corrected_sign_[2]) ? n_cell * (-1) : n_cell;
          const int m_dir_elem = (pure_neighbours_corrected_sign_[1]) ? (m + 1) * (-1) + 1 : m;
          const double indic_neighbour = ref_ijk_ft_->itfce().I()(index_i_ + l_dir, index_j_ + m_dir_elem, index_k_ + n_dir);
          if (indic_neighbour > LIQUID_INDICATOR_TEST)
            {
              pure_neighbours_last_faces_to_correct_[1][l_cell][m][n_cell] = true;
              const double lmn_zero = (m > 0) ? 1. : 0.;
              const double contrib_factor = (pure_neighbours_corrected_sign_[1]) ? (lmn_zero * (2 * abs(m_dir) + 1) - (1. - lmn_zero)) * (-1):
                                            lmn_zero * (2 * (m_dir - 1) + 1) - (1. - lmn_zero);
              const double dx_contrib = l_dir * dx;
              const double dy_contrib = contrib_factor * dy_over_two;
              const double dz_contrib = n_dir * dz;
              Vecteur3 distance_contrib = {dx_contrib, dy_contrib, dz_contrib};
              pure_neighbours_last_faces_corrected_distance_[1][l_cell][m][n_cell] = cell_centre_distance_
                                                                                     + Vecteur3::produit_scalaire(normal_vector_compo_, distance_contrib);
              if (pure_neighbours_last_faces_corrected_distance_[1][l_cell][m][n_cell] < 0)
                pure_neighbours_last_faces_to_correct_[1][l_cell][m][n_cell] = false;
              if (neighbours_last_faces_weighting_)
                {
                  const double colinearity = compute_cell_faces_weighting(dx_contrib, dy_contrib, dz_contrib, 1);
                  pure_neighbours_last_faces_corrected_colinearity_[1][l_cell][m][n_cell] = colinearity;
                }
            }
        }

  for (l_cell=dx_increment; l_cell>=0; l_cell--)
    for (m_cell=dy_increment; m_cell>=0; m_cell--)
      for (n=dz_over_two_increment + 1; n>=0; n--)
        {
          const int l_dir = (pure_neighbours_corrected_sign_[0]) ? l_cell * (-1) : l_cell;
          const int m_dir = (pure_neighbours_corrected_sign_[1]) ? m_cell * (-1) : m_cell;
          const int n_dir = (pure_neighbours_corrected_sign_[2]) ? n * (-1) + 1: n;
          const int n_dir_elem = (pure_neighbours_corrected_sign_[2]) ? (n + 1) * (-1) + 1 : n;
          const double indic_neighbour = ref_ijk_ft_->itfce().I()(index_i_ + l_dir, index_j_ + m_dir, index_k_ + n_dir_elem);
          if (indic_neighbour > LIQUID_INDICATOR_TEST)
            {
              pure_neighbours_last_faces_to_correct_[2][l_cell][m_cell][n] = true;
              const double lmn_zero = (n > 0) ? 1. : 0.;
              const double contrib_factor = (pure_neighbours_corrected_sign_[2]) ? (lmn_zero * (2 * abs(n_dir) + 1) - (1. - lmn_zero)) * (-1):
                                            lmn_zero * (2 * (n_dir - 1) + 1) - (1. - lmn_zero);
              const double dx_contrib = l_dir * dx;
              const double dy_contrib = m_dir * dy;
              const double dz_contrib = contrib_factor * dz_over_two;
              Vecteur3 distance_contrib = {dx_contrib, dy_contrib, dz_contrib};
              pure_neighbours_last_faces_corrected_distance_[2][l_cell][m_cell][n] = cell_centre_distance_ +
                                                                                     Vecteur3::produit_scalaire(normal_vector_compo_, distance_contrib);
              if (pure_neighbours_last_faces_corrected_distance_[2][l_cell][m_cell][n] < 0)
                pure_neighbours_last_faces_to_correct_[2][l_cell][m_cell][n] = false;
              if (neighbours_last_faces_weighting_)
                {
                  const double colinearity = compute_cell_faces_weighting(dx_contrib, dy_contrib, dz_contrib, 2);
                  pure_neighbours_last_faces_corrected_colinearity_[2][l_cell][m_cell][n] = colinearity;
                }
            }
        }
}

double IJK_One_Dimensional_Subproblem_Geometry::compute_cell_faces_weighting(const double& dx_contrib,
                                                                             const double& dy_contrib,
                                                                             const double& dz_contrib,
                                                                             const int& dir)
{
  if (neighbours_last_faces_colinearity_weighting_)
    return compute_colinearity(dx_contrib, dy_contrib, dz_contrib);
  if (neighbours_last_faces_colinearity_face_weighting_)
    return compute_colinearity_cell_faces(dx_contrib, dy_contrib, dz_contrib, dir);
  if (neighbours_last_faces_distance_weighting_)
    return compute_distance_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  if (neighbours_last_faces_distance_colinearity_weighting_)
    return compute_colinearity(dx_contrib, dy_contrib, dz_contrib) * compute_distance_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  if (neighbours_last_faces_distance_colinearity_face_weighting_)
    return compute_colinearity_cell_faces(dx_contrib, dy_contrib, dz_contrib, dir) * compute_distance_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  return 1;
}

Vecteur3 IJK_One_Dimensional_Subproblem_Geometry::compute_relative_vector_cell_faces(const double& dx_contrib,
                                                                                     const double& dy_contrib,
                                                                                     const double& dz_contrib)
{
  Vecteur3 relative_vector = normal_vector_compo_;
  relative_vector *= cell_centre_distance_;
  Vecteur3 tangential_relative_vector = tangential_distance_vector_;
  tangential_relative_vector *= cell_centre_tangential_distance_;
  relative_vector += tangential_relative_vector;
  Vecteur3 dxyz_contrib = {dx_contrib, dy_contrib, dz_contrib};
  relative_vector += dxyz_contrib;
  return relative_vector;
}

double IJK_One_Dimensional_Subproblem_Geometry::compute_colinearity(const double& dx_contrib,
                                                                    const double& dy_contrib,
                                                                    const double& dz_contrib)
{
  Vecteur3 relative_vector = compute_relative_vector_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  const double relative_vector_norm = relative_vector.length();
  relative_vector *= (1 / relative_vector_norm);
  const double colinearity = Vecteur3::produit_scalaire(normal_vector_compo_, relative_vector);
  return abs(colinearity);
}

double IJK_One_Dimensional_Subproblem_Geometry::compute_colinearity_cell_faces(const double& dx_contrib,
                                                                               const double& dy_contrib,
                                                                               const double& dz_contrib,
                                                                               const int& dir)
{
  Vecteur3 relative_vector = compute_relative_vector_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  const double relative_vector_norm = relative_vector.length();
  relative_vector *= (1 / relative_vector_norm);
  return abs(relative_vector[dir]);
}

double IJK_One_Dimensional_Subproblem_Geometry::compute_distance_cell_faces(const double& dx_contrib,
                                                                            const double& dy_contrib,
                                                                            const double& dz_contrib)
{
  Vecteur3 relative_vector = compute_relative_vector_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  const double distance = Vecteur3::produit_scalaire(tangential_distance_vector_, relative_vector);
  return abs(1 / (distance + 1e-16));
}

int IJK_One_Dimensional_Subproblem_Geometry::get_dxyz_increment_max()
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);

  int dxyz_increment_max;
  if (!correct_neighbours_rank_)
    {
      const int dx_increment_max = (int) ((dx + probe_length_) / dx);
      const int dy_increment_max = (int) ((dy + probe_length_) / dy);
      const int dz_increment_max = (int) ((dz + probe_length_) / dz);
      dxyz_increment_max = std::max(std::max(dx_increment_max, dy_increment_max), dz_increment_max);
    }
  else
    {
      dxyz_increment_max = neighbours_corrected_rank_;
    }
  return dxyz_increment_max;
}

int IJK_One_Dimensional_Subproblem_Geometry::get_dxyz_over_two_increment_max()
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);

  int dxyz_over_two_increment_max;
  if (!correct_neighbours_rank_)
    {
      const int dx_increment_max = (int) ((probe_length_ + dx) / (dx / 2.));
      const int dy_increment_max = (int) ((probe_length_ + dy) / (dy / 2.));
      const int dz_increment_max = (int) ((probe_length_ + dz) / (dz / 2.));
      dxyz_over_two_increment_max = std::max(std::max(dx_increment_max, dy_increment_max), dz_increment_max);
    }
  else
    {
      dxyz_over_two_increment_max = neighbours_face_corrected_rank_;
    }
  return dxyz_over_two_increment_max;
}

void IJK_One_Dimensional_Subproblem_Geometry::get_maximum_remaining_distance(int& dx_remaining,
                                                                             int& dy_remaining,
                                                                             int& dz_remaining)
{
  double remaining_distance_diag = probe_length_ - cell_centre_distance_;
  Vecteur3 remaining_distance_diag_projected = normal_vector_compo_;
  dx_remaining = (int) (remaining_distance_diag / remaining_distance_diag_projected[0]);
  dy_remaining = (int) (remaining_distance_diag / remaining_distance_diag_projected[1]);
  dz_remaining = (int) (remaining_distance_diag / remaining_distance_diag_projected[2]);
}

