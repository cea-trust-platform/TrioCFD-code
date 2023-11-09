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
// File      : IJK_One_Dimensional_Subproblems.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_One_Dimensional_Subproblems.h>
#include <IJK_FT.h>
#include <IJK_switch_FT.h>

Implemente_instanciable_sans_constructeur(IJK_One_Dimensional_Subproblems, "IJK_One_Dimensional_Subproblems", LIST(IJK_One_Dimensional_Subproblem));

IJK_One_Dimensional_Subproblems::IJK_One_Dimensional_Subproblems()
{
  interfacial_thermal_flux_per_bubble_.set_smart_resize(1);
  total_surface_per_bubble_.set_smart_resize(1);
  overall_nusselt_number_per_bubble_.set_smart_resize(1);
  overall_shear_stress_per_bubble_.set_smart_resize(1);
  overall_shear_force_per_bubble_.set_smart_resize(1);
}

IJK_One_Dimensional_Subproblems::IJK_One_Dimensional_Subproblems(const IJK_FT_double& ijk_ft) : IJK_One_Dimensional_Subproblems()
{
  ref_ijk_ft_ = ijk_ft;
}

Sortie& IJK_One_Dimensional_Subproblems::printOn( Sortie& os ) const
{
  return os;
}

Entree& IJK_One_Dimensional_Subproblems::readOn( Entree& is )
{
  LIST(IJK_One_Dimensional_Subproblem)::readOn( is );
  return is;
}

void IJK_One_Dimensional_Subproblems::initialise_thermal_subproblems_list_params(const int& pre_initialise_thermal_subproblems_list,
                                                                                 const double& pre_factor_subproblems_number)
{
  pre_initialise_thermal_subproblems_list_ = pre_initialise_thermal_subproblems_list;
  pre_factor_subproblems_number_ = pre_factor_subproblems_number;
}

void IJK_One_Dimensional_Subproblems::clean()
{
  clean(pre_initialise_thermal_subproblems_list_);
}

void IJK_One_Dimensional_Subproblems::clean(int add)
{
  if (add)
    clean_add();
  else
    clean_remove();
  subproblems_counter_ = 0;
}

void IJK_One_Dimensional_Subproblems::clean_remove()
{
  vide();
}

void IJK_One_Dimensional_Subproblems::clean_add()
{
  complete_subproblems();
}

void IJK_One_Dimensional_Subproblems::complete_subproblems()
{
  if (init_)
    {
      max_subproblems_ = (int) (pre_factor_subproblems_number_ * subproblems_counter_);
      if (subproblems_counter_ < max_subproblems_)//
        {
          const int sub_problem_end_index = subproblems_counter_ - 1;
          IJK_One_Dimensional_Subproblem subproblem = (*this)[sub_problem_end_index];
          while(subproblems_counter_ < max_subproblems_)
            {
              (*this).add(subproblem);
              subproblems_counter_ ++;
            }
          init_ = 0;
        }
    }
}

void IJK_One_Dimensional_Subproblems::add_subproblems(int n)
{
  for (int i=0; i<n; i++)
    {
      IJK_One_Dimensional_Subproblem subproblem;
      (*this).add(subproblem);
    }
  max_subproblems_ = n;
}

void IJK_One_Dimensional_Subproblems::associate_sub_problem_to_inputs(int debug,
                                                                      int i, int j, int k,
                                                                      double global_time_step,
                                                                      double current_time,
                                                                      const IJK_Field_double& eulerian_compo_connex,
                                                                      const IJK_Field_double& eulerian_distance,
                                                                      const IJK_Field_double& eulerian_curvature,
                                                                      const IJK_Field_double& eulerian_interfacial_area,
                                                                      const FixedVector<IJK_Field_double, 3>& eulerian_facets_barycentre,
                                                                      const FixedVector<IJK_Field_double, 3>& eulerian_normal_vectors,
                                                                      const ArrOfDouble& rising_velocities,
                                                                      const DoubleTab& rising_vectors,
                                                                      const DoubleTab& bubbles_barycentre,
                                                                      int advected_frame_of_reference,
                                                                      int neglect_frame_of_reference_radial_advection,
                                                                      const int& points_per_thermal_subproblem,
                                                                      const double& alpha,
                                                                      const double& lambda,
                                                                      const double& coeff_distance_diagonal,
                                                                      const double& cell_diagonal,
                                                                      const double& dr_base,
                                                                      const DoubleVect& radial_coordinates,
                                                                      const Matrice& identity_matrix_explicit_implicit,
                                                                      const Matrice& radial_first_order_operator_raw,
                                                                      const Matrice& radial_second_order_operator_raw,
                                                                      const Matrice& radial_first_order_operator,
                                                                      const Matrice& radial_second_order_operator,
                                                                      Matrice& identity_matrix_subproblems,
                                                                      Matrice& radial_diffusion_matrix,
                                                                      Matrice& radial_convection_matrix,
                                                                      const IJK_Interfaces& interfaces,
                                                                      const double& indicator,
                                                                      const IJK_Field_double& temperature,
                                                                      const IJK_Field_double& temperature_ft,
                                                                      const IJK_Field_double& temperature_before_extrapolation,
                                                                      const FixedVector<IJK_Field_double, 3>& velocity,
                                                                      const FixedVector<IJK_Field_double, 3>& velocity_ft,
                                                                      const IJK_Field_double& pressure,
                                                                      const FixedVector<IJK_Field_double, 3>& grad_T_elem,
                                                                      const FixedVector<IJK_Field_double, 3>& hess_diag_T_elem,
                                                                      const FixedVector<IJK_Field_double, 3>& hess_cross_T_elem,
                                                                      IJK_Finite_Difference_One_Dimensional_Matrix_Assembler& finite_difference_assembler,
                                                                      Matrice& thermal_subproblems_matrix_assembly,
                                                                      DoubleVect& thermal_subproblems_rhs_assembly,
                                                                      DoubleVect& thermal_subproblems_temperature_solution_ini,
                                                                      DoubleVect& thermal_subproblems_temperature_solution,
                                                                      const int& source_terms_type,
                                                                      const int& source_terms_correction,
                                                                      const bool& is_first_time_step,
                                                                      int& first_time_step_temporal,
                                                                      const int& first_time_step_explicit,
                                                                      const double& local_fourier,
                                                                      const double& local_cfl,
                                                                      const double& min_delta_xyz,
                                                                      const double& delta_T_subcooled_overheated,
                                                                      const int& first_time_step_varying_probes,
                                                                      const int& probe_variations_priority,
                                                                      const int& disable_interpolation_in_mixed_cells,
                                                                      const int& max_u_radial,
                                                                      const int& correct_fluxes,
                                                                      const int& distance_cell_faces_from_lrs,
                                                                      const int& pre_initialise_thermal_subproblems_list,
                                                                      const int& correct_temperature_cell_neighbours,
                                                                      const int& correct_neighbours_rank,
                                                                      const int& neighbours_corrected_rank,
                                                                      const int& neighbours_colinearity_weighting,
                                                                      const int& compute_reachable_fluxes,
                                                                      const int& find_cell_neighbours_for_fluxes_spherical_correction,
                                                                      const int& n_iter_distance)
{
  if (!init_ && subproblems_counter_ > max_subproblems_)
    {
      Cerr << max_subproblems_ << "subproblems were expected but" << subproblems_counter_ << "subproblem try to be associated with the list of subproblems";
      Process::exit();
    }

  debug_ = debug;
  ArrOfDouble bubble_rising_vector(3);
  ArrOfDouble normal_vector(3);
  ArrOfDouble facet_barycentre(3);
  ArrOfDouble bubble_barycentre(3);

  if (debug_)
    Cerr << "Mixed cell indices (i,j,k) : (" << i << ";" << j << ";" << k << ")" << finl;
  const int compo_connex = int(eulerian_compo_connex(i, j, k));
  if (debug_)
    {
      Cerr << "compo_connex : " << compo_connex << finl;
      Cerr << "bubbles_barycentre : " << bubbles_barycentre << finl;
    }

  // Need for a Navier-Stokes field (NOT FT)
  const double distance = eulerian_distance(i, j ,k);
  const double curvature = eulerian_curvature(i, j ,k);
  const double interfacial_area = eulerian_interfacial_area(i, j ,k);

  IJK_Splitting splitting = eulerian_compo_connex.get_splitting();
  const double bubble_rising_velocity = rising_velocities(compo_connex);
  for (int dir=0; dir < 3; dir++)
    {
      facet_barycentre(dir) = eulerian_facets_barycentre[dir](i, j, k);
      normal_vector(dir) = eulerian_normal_vectors[dir](i, j, k);
      bubble_rising_vector(dir) = rising_vectors(compo_connex, dir);
      bubble_barycentre(dir) = bubbles_barycentre(compo_connex, dir);
    }

  if (init_)
    {
      IJK_One_Dimensional_Subproblem subproblem(ref_ijk_ft_);
      (*this).add(subproblem);
    }
  (*this)[subproblems_counter_].associate_sub_problem_to_inputs(init_,
                                                                debug,
                                                                subproblems_counter_,
                                                                i, j, k,
                                                                global_time_step,
                                                                current_time,
                                                                compo_connex,
                                                                distance,
                                                                curvature,
                                                                interfacial_area,
                                                                facet_barycentre,
                                                                normal_vector,
                                                                bubble_rising_velocity,
                                                                bubble_rising_vector,
                                                                bubble_barycentre,
                                                                advected_frame_of_reference,
                                                                neglect_frame_of_reference_radial_advection,
                                                                points_per_thermal_subproblem,
                                                                alpha,
                                                                lambda,
                                                                coeff_distance_diagonal,
                                                                cell_diagonal,
                                                                dr_base,
                                                                radial_coordinates,
                                                                identity_matrix_explicit_implicit,
                                                                radial_first_order_operator_raw,
                                                                radial_second_order_operator_raw,
                                                                radial_first_order_operator,
                                                                radial_second_order_operator,
                                                                identity_matrix_subproblems,
                                                                radial_diffusion_matrix,
                                                                radial_convection_matrix,
                                                                interfaces,
                                                                indicator,
                                                                eulerian_distance,
                                                                eulerian_curvature,
                                                                eulerian_interfacial_area,
                                                                eulerian_normal_vectors,
                                                                eulerian_facets_barycentre,
                                                                temperature,
                                                                temperature_ft,
                                                                temperature_before_extrapolation,
                                                                velocity,
                                                                velocity_ft,
                                                                pressure,
                                                                grad_T_elem,
                                                                hess_diag_T_elem,
                                                                hess_cross_T_elem,
                                                                finite_difference_assembler,
                                                                thermal_subproblems_matrix_assembly,
                                                                thermal_subproblems_rhs_assembly,
                                                                thermal_subproblems_temperature_solution_ini,
                                                                thermal_subproblems_temperature_solution,
                                                                source_terms_type,
                                                                source_terms_correction,
                                                                is_first_time_step,
                                                                first_time_step_temporal,
                                                                first_time_step_explicit,
                                                                local_fourier,
                                                                local_cfl,
                                                                min_delta_xyz,
                                                                delta_T_subcooled_overheated,
                                                                first_time_step_varying_probes,
                                                                probe_variations_priority,
                                                                disable_interpolation_in_mixed_cells,
                                                                max_u_radial,
                                                                correct_fluxes,
                                                                distance_cell_faces_from_lrs,
                                                                pre_initialise_thermal_subproblems_list,
                                                                correct_temperature_cell_neighbours,
                                                                correct_neighbours_rank,
                                                                neighbours_corrected_rank,
                                                                neighbours_colinearity_weighting,
                                                                compute_reachable_fluxes,
                                                                find_cell_neighbours_for_fluxes_spherical_correction,
                                                                n_iter_distance);

  subproblems_counter_++;
}

void IJK_One_Dimensional_Subproblems::interpolate_project_velocities_on_probes()
{
  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].interpolate_project_velocities_on_probes();
}

void IJK_One_Dimensional_Subproblems::reajust_probes_length()
{
  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].reajust_probe_length();
}

int IJK_One_Dimensional_Subproblems::get_probe_variations_enabled(const int& probe_variations_priority)
{
  if (probe_variations_priority)
    return get_probe_variations_enabled_priority();
  else
    return get_probe_variations_enabled_non_priority();
}

int IJK_One_Dimensional_Subproblems::get_probe_variations_enabled_priority()
{
  int probe_variations_enabled = 0;
  for (int itr=0; itr < subproblems_counter_; itr++)
    probe_variations_enabled = (probe_variations_enabled || (*this)[itr].get_probe_variations_enabled());
  return probe_variations_enabled;
}

int IJK_One_Dimensional_Subproblems::get_probe_variations_enabled_non_priority()
{
  int probe_variations_enabled = 1;
  for (int itr=0; itr < subproblems_counter_; itr++)
    probe_variations_enabled = (probe_variations_enabled && (*this)[itr].get_probe_variations_enabled());
  return probe_variations_enabled;
}

void IJK_One_Dimensional_Subproblems::compute_modified_probe_length(const int& probe_variations_enabled)
{
  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].compute_modified_probe_length(probe_variations_enabled);
}

void IJK_One_Dimensional_Subproblems::compute_radial_convection_diffusion_operators()
{
  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].compute_radial_convection_diffusion_operators();
}

void IJK_One_Dimensional_Subproblems::compute_source_terms_impose_boundary_conditions(DoubleVect& thermal_subproblems_rhs_assembly,
                                                                                      DoubleVect& thermal_subproblems_temperature_solution_ini,
                                                                                      const int& boundary_condition_interface,
                                                                                      const double& interfacial_boundary_condition_value,
                                                                                      const int& impose_boundary_condition_interface_from_simulation,
                                                                                      const int& boundary_condition_end,
                                                                                      const double& end_boundary_condition_value,
                                                                                      const int& impose_user_boundary_condition_end_value)
{
  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].compute_source_terms_impose_boundary_conditions(thermal_subproblems_rhs_assembly,
                                                                 thermal_subproblems_temperature_solution_ini,
                                                                 boundary_condition_interface,
                                                                 interfacial_boundary_condition_value,
                                                                 impose_boundary_condition_interface_from_simulation,
                                                                 boundary_condition_end,
                                                                 end_boundary_condition_value,
                                                                 impose_user_boundary_condition_end_value);
}

void IJK_One_Dimensional_Subproblems::approximate_temperature_increment_material_derivative()
{
  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].approximate_temperature_increment_material_derivative();
}

void IJK_One_Dimensional_Subproblems::retrieve_radial_quantities()
{
  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].retrieve_radial_quantities();
}

void IJK_One_Dimensional_Subproblems::retrieve_temperature_solutions()
{
  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].retrieve_temperature_solution();
}

void IJK_One_Dimensional_Subproblems::compute_local_temperature_gradient_solutions()
{
  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].compute_local_temperature_gradient_solution();
}

void IJK_One_Dimensional_Subproblems::compute_local_velocity_gradient()
{
  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].compute_local_velocity_gradient();
}

void IJK_One_Dimensional_Subproblems::get_subproblem_ijk_indices(int& i, int& j, int& k, int& subproblem_index) const
{
  (*this)[subproblem_index].get_ijk_indices(i,j,k);
}

const int& IJK_One_Dimensional_Subproblems::get_dxyz_increment_bool(const int& subproblem_index) const
{
  return (*this)[subproblem_index].get_dxyz_increment_bool();
}

const int& IJK_One_Dimensional_Subproblems::get_dxyz_over_two_increment_bool(const int& subproblem_index) const
{
  return (*this)[subproblem_index].get_dxyz_over_two_increment_bool();
}

const FixedVector<int,3>& IJK_One_Dimensional_Subproblems::get_pure_neighbours_corrected_sign(const int& subproblem_index) const
{
  return (*this)[subproblem_index].get_pure_neighbours_corrected_sign();
}

const std::vector<std::vector<std::vector<bool>>>& IJK_One_Dimensional_Subproblems::get_pure_neighbours_to_correct(const int& subproblem_index) const
{
  return (*this)[subproblem_index].get_pure_neighbours_to_correct();
}

const std::vector<std::vector<std::vector<double>>>& IJK_One_Dimensional_Subproblems::get_pure_neighbours_corrected_distance(const int& subproblem_index) const
{
  return (*this)[subproblem_index].get_pure_neighbours_corrected_distance();
}

const std::vector<std::vector<std::vector<double>>>& IJK_One_Dimensional_Subproblems::get_pure_neighbours_corrected_colinearity(const int& subproblem_index) const
{
  return (*this)[subproblem_index].get_pure_neighbours_corrected_colinearity();
}

const std::vector<std::vector<std::vector<std::vector<bool>>>> IJK_One_Dimensional_Subproblems::get_pure_neighbours_last_faces_to_correct(const int& subproblem_index) const
{
  return (*this)[subproblem_index].get_pure_neighbours_last_faces_to_correct();
}
const std::vector<std::vector<std::vector<std::vector<double>>>> IJK_One_Dimensional_Subproblems::get_pure_neighbours_last_faces_corrected_distance(const int& subproblem_index) const
{
  return (*this)[subproblem_index].get_pure_neighbours_last_faces_corrected_distance();
}
const std::vector<std::vector<std::vector<std::vector<double>>>> IJK_One_Dimensional_Subproblems::get_pure_neighbours_last_faces_corrected_colinearity(const int& subproblem_index) const
{
  return (*this)[subproblem_index].get_pure_neighbours_last_faces_corrected_colinearity();
}

double IJK_One_Dimensional_Subproblems::get_interfacial_gradient_corrected(int i)
{
  return (*this)[i].get_interfacial_gradient_corrected();
}

double IJK_One_Dimensional_Subproblems::get_temperature_profile_at_point(const int& i, const double& dist) const
{
  return (*this)[i].get_temperature_profile_at_point(dist);
}

const Vecteur3& IJK_One_Dimensional_Subproblems::get_bary_facet(const int& i) const
{
  return (*this)[i].get_bary_facet();
}

const double& IJK_One_Dimensional_Subproblems::get_dist_cell_interface(const int& i) const
{
  return (*this)[i].get_dist_cell();
}

const FixedVector<double,6>& IJK_One_Dimensional_Subproblems::get_dist_faces_interface(const int& i) const
{
  return (*this)[i].get_dist_faces();
}

double IJK_One_Dimensional_Subproblems::get_temperature_times_velocity_profile_at_point(const int& i, const double& dist, const int& dir) const
{
  return (*this)[i].get_temperature_times_velocity_profile_at_point(dist, dir);
}

DoubleVect IJK_One_Dimensional_Subproblems::get_temperature_profile_discrete_integral_at_point(const int& i, const double& dist, const int& level, const int& dir) const
{
  return (*this)[i].get_temperature_profile_discrete_integral_at_point(dist, level, dir);
}

DoubleVect IJK_One_Dimensional_Subproblems::get_temperature_times_velocity_profile_discrete_integral_at_point(const int& i, const double& dist, const int& level, const int& dir) const
{
  return (*this)[i].get_temperature_times_velocity_profile_discrete_integral_at_point(dist, level, dir);
}

double IJK_One_Dimensional_Subproblems::get_temperature_gradient_profile_at_point(const int& i, const double& dist, const int& dir) const
{
  return (*this)[i].get_temperature_gradient_profile_at_point(dist, dir);
}

double IJK_One_Dimensional_Subproblems::get_temperature_gradient_times_conductivity_profile_at_point(const int& i, const double& dist, const int& dir) const
{
  return (*this)[i].get_temperature_gradient_times_conductivity_profile_at_point(dist, dir);
}

DoubleVect IJK_One_Dimensional_Subproblems::get_temperature_gradient_profile_discrete_integral_at_point(const int& i, const double& dist, const int& level, const int& dir) const
{
  return (*this)[i].get_temperature_gradient_profile_discrete_integral_at_point(dist, level, dir);
}

DoubleVect IJK_One_Dimensional_Subproblems::get_temperature_gradient_times_conductivity_profile_discrete_integral_at_point(const int& i, const double& dist, const int& level, const int& dir) const
{
  return (*this)[i].get_temperature_gradient_times_conductivity_profile_discrete_integral_at_point(dist, level, dir);
}

void IJK_One_Dimensional_Subproblems::thermal_subresolution_outputs(SFichier& fic, const int rank)
{
  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].thermal_subresolution_outputs(fic, rank);
  post_process_overall_bubbles_quantities(rank);
}

double IJK_One_Dimensional_Subproblems::get_min_temperature() const
{
  return (*this)[0].get_min_temperature();
}

double IJK_One_Dimensional_Subproblems::get_max_temperature() const
{
  return (*this)[0].get_max_temperature();
}

double IJK_One_Dimensional_Subproblems::get_min_temperature_domain_ends() const
{
  return (*this)[0].get_min_temperature_domain_ends();
}

double IJK_One_Dimensional_Subproblems::get_max_temperature_domain_ends() const
{
  return (*this)[0].get_max_temperature_domain_ends();
}

double IJK_One_Dimensional_Subproblems::get_min_euler_time_step(int& nb_iter_explicit)
{
  double min_euler_time_step = 1e20;
  nb_iter_explicit = 1;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      min_euler_time_step = std::min(min_euler_time_step, (*this)[itr].get_local_time_step_round());
      nb_iter_explicit = std::max(nb_iter_explicit, (*this)[itr].get_nb_iter_explicit());
    }
  return min_euler_time_step;
}

double IJK_One_Dimensional_Subproblems::get_local_max_fourier_time_step_probe_length()
{
  double max_local_fourier_time_step_probe_length = 0.;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      max_local_fourier_time_step_probe_length = std::max(max_local_fourier_time_step_probe_length,
                                                          (*this)[itr].get_local_fourier_time_step_probe_length());
    }
  return max_local_fourier_time_step_probe_length;
}


double IJK_One_Dimensional_Subproblems::get_local_max_cfl_time_step_probe_length()
{
  double max_local_cfl_time_step_probe_length = 0.;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      max_local_cfl_time_step_probe_length = std::max(max_local_cfl_time_step_probe_length,
                                                      (*this)[itr].get_local_cfl_time_step_probe_length());
    }
  return max_local_cfl_time_step_probe_length;
}

double IJK_One_Dimensional_Subproblems::get_local_min_fourier_time_step_probe_length()
{
  double min_local_fourier_time_step_probe_length = 1.e20;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      min_local_fourier_time_step_probe_length = std::min(min_local_fourier_time_step_probe_length,
                                                          (*this)[itr].get_local_fourier_time_step_probe_length());
    }
  return min_local_fourier_time_step_probe_length;
}


double IJK_One_Dimensional_Subproblems::get_local_min_cfl_time_step_probe_length()
{
  double min_local_cfl_time_step_probe_length = 1.e20;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      min_local_cfl_time_step_probe_length = std::min(min_local_cfl_time_step_probe_length,
                                                      (*this)[itr].get_local_cfl_time_step_probe_length());
    }
  return min_local_cfl_time_step_probe_length;
}

double IJK_One_Dimensional_Subproblems::get_local_dt_cfl()
{
  double min_local_dt_cfl = 1.e20;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      min_local_dt_cfl = std::min(min_local_dt_cfl, (*this)[itr].get_local_dt_cfl());
    }
  return min_local_dt_cfl;
}

double IJK_One_Dimensional_Subproblems::get_local_dt_cfl_min_delta_xyz()
{
  double min_local_dt_cfl_min_delta_xyz = 1.e20;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      min_local_dt_cfl_min_delta_xyz = std::min(min_local_dt_cfl_min_delta_xyz, (*this)[itr].get_local_dt_cfl_min_delta_xyz());
    }
  return min_local_dt_cfl_min_delta_xyz;
}

void IJK_One_Dimensional_Subproblems::set_local_time_step(const double& local_time_step)
{
  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].set_local_time_step(local_time_step);
}

void IJK_One_Dimensional_Subproblems::prepare_temporal_schemes()
{
  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].prepare_temporal_schemes();
}

const int& IJK_One_Dimensional_Subproblems::get_end_index_subproblem(const int index) const
{
  return (*this)[index].get_end_index_subproblem();
}

//static std::vector<int> arg_min(ArrOfDouble theta_phi_scope)
//{
//  const int n = theta_phi_scope.size();
//  // IntVect indices(n);
//  std::vector<int> indices(n);
//  for (int j=0; j<n; j++)
//    indices[j]=j;
//  std::sort(indices.begin(), indices.end(), [&theta_phi_scope](int i, int j) {return theta_phi_scope[i] < theta_phi_scope[j];});
//  return indices;
//}

void IJK_One_Dimensional_Subproblems::post_processed_all_probes()
{
  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].set_post_processing_theta_phi_scope();
}

void IJK_One_Dimensional_Subproblems::sort_limited_probes_spherical_coords_post_processing(const int nb_theta, const int nb_phi,
                                                                                           const int theta_diag_val, const int phi_diag_val)
{
  ArrOfDouble r_sph(subproblems_counter_);
  ArrOfDouble theta_sph(subproblems_counter_);
  ArrOfDouble phi_sph(subproblems_counter_);
  ArrOfDouble theta_scope;
  ArrOfDouble phi_scope;
  int nb_theta_even = nb_theta;
  int nb_phi_even = nb_phi;
  if (nb_theta_even % 2)
    nb_theta_even += 1;
  if (nb_phi_even % 2)
    nb_phi_even += 1;
  theta_scope.set_smart_resize(1);
  phi_scope.set_smart_resize(1);
  int i, j, k;
  for (i=0; i<subproblems_counter_; i++)
    {
      r_sph(i) = (*this)[i].get_radius_spherical_coords();
      theta_sph(i) = (*this)[i].get_theta_spherical_coords();
      phi_sph(i) = (*this)[i].get_phi_spherical_coords();
    }
  const int nb_outputs = nb_phi_even * nb_theta_even;
  int theta_incr;
  int phi_incr;
  theta_incr = (int) (M_PI / (double) nb_theta_even);
  phi_incr = (int) ((2 * M_PI) / (double) nb_phi_even);
  if (theta_diag_val)
    for (i=0; i<nb_theta_even; i++)
      theta_scope.append_array(theta_incr * (i + 0.5));
  else
    for (i=0; i<nb_theta_even; i++)
      theta_scope.append_array(theta_incr * i);
  if (phi_diag_val)
    for (i=0; i<nb_phi_even; i++)
      phi_scope.append_array(phi_incr * (i + 0.5));
  else
    for (i=0; i<nb_phi_even; i++)
      phi_scope.append_array(phi_incr * i);
  /*
   * Sort by phi and theta simultaneously
   */
  ArrOfDouble radius_outputs(nb_outputs);
  ArrOfDouble theta_outputs(nb_outputs);
  ArrOfDouble phi_outputs(nb_outputs);
  int phi_theta_counter = 0;
  for (j=0; j<nb_phi_even; j++)
    for (i=0; i<nb_theta_even; i++)
      {
        ArrOfDouble theta_diff = theta_sph;
        ArrOfDouble phi_diff = phi_sph;
        std::vector<double> sum_errors_theta_phi_scope;
        theta_diff -= theta_scope(i);
        phi_diff -= phi_scope(j);
        for (k=0; k<subproblems_counter_; k++)
          {
            theta_diff(k) = abs(theta_diff(k));
            phi_diff(k) = abs(phi_diff(k));
            sum_errors_theta_phi_scope.push_back(theta_diff(k));
            sum_errors_theta_phi_scope[k] += phi_diff(k);
          }
        const int theta_phi_scope_index = (int) std::distance(sum_errors_theta_phi_scope.begin(), std::min_element(sum_errors_theta_phi_scope.begin(), sum_errors_theta_phi_scope.end()));
        (*this)[theta_phi_scope_index].set_post_processing_theta_phi_scope();
        radius_outputs(phi_theta_counter) = r_sph(theta_phi_scope_index);
        theta_outputs(phi_theta_counter) = theta_sph(theta_phi_scope_index);
        phi_outputs(phi_theta_counter) = phi_sph(theta_phi_scope_index);
        phi_theta_counter ++;
      }
}

void IJK_One_Dimensional_Subproblems::compute_overall_quantities_per_bubbles(const double& delta_temperature,
                                                                             const double& lambda)
{
  std::vector<int> compo_found;
  int local_compo;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      local_compo = (*this)[itr].get_compo();
      if (std::find(compo_found.begin(), compo_found.end(), local_compo) == compo_found.end())
        compo_found.push_back(local_compo);
    }
  nb_bubbles_ = (int) compo_found.size();

  interfacial_thermal_flux_per_bubble_.resize(nb_bubbles_);
  total_surface_per_bubble_.resize(nb_bubbles_);
  overall_nusselt_number_per_bubble_.resize(nb_bubbles_);
  overall_shear_stress_per_bubble_.resize(nb_bubbles_);
  overall_shear_force_per_bubble_.resize(nb_bubbles_);
  compute_nusselt_numbers_per_bubbles(delta_temperature, lambda);
  compute_shear_per_bubbles();
  compo_found.clear();
}

void IJK_One_Dimensional_Subproblems::compute_overall_bubbles_quantities(const double& delta_temperature,
                                                                         const double& lambda)
{
  compute_overall_quantities_per_bubbles(delta_temperature, lambda);
  compute_overall_nusselt_numbers_bubbles();
  is_updated_ = true;
}

void IJK_One_Dimensional_Subproblems::compute_nusselt_numbers_per_bubbles(const double& delta_temperature,
                                                                          const double& lambda)
{
  overall_nusselt_number_ = 0.;
  interfacial_thermal_flux_ = 0.;
  total_surface_ = 0.;
  lambda_ = lambda;
  delta_temperature_ = delta_temperature;
  int local_compo;
  interfacial_thermal_flux_per_bubble_ *= 0.;
  total_surface_per_bubble_ *= 0.;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      local_compo = (*this)[itr].get_compo();
      interfacial_thermal_flux_per_bubble_(local_compo) += (*this)[itr].get_interfacial_thermal_flux();
      total_surface_per_bubble_(local_compo) += (*this)[itr].get_local_surface_area();
    }
  for (int i=0; i < nb_bubbles_; i++)
    overall_nusselt_number_per_bubble_(i) = (interfacial_thermal_flux_per_bubble_(i) * caracteristic_length_) / (total_surface_per_bubble_(i) * delta_temperature_ * lambda_);
}

void IJK_One_Dimensional_Subproblems::compute_shear_per_bubbles()
{
  int local_compo;
  overall_shear_force_per_bubble_ *= 0.;
  overall_shear_force_ = 0.;
  overall_shear_stress_ = 0.;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      local_compo = (*this)[itr].get_compo();
      overall_shear_force_per_bubble_(local_compo) += (*this)[itr].get_shear_force();
    }
  overall_shear_stress_per_bubble_ = overall_shear_force_per_bubble_;
  for (int i=0; i < nb_bubbles_; i++)
    overall_shear_stress_per_bubble_(i) *= (1 / total_surface_per_bubble_(i));
  for (int i=0; i < nb_bubbles_; i++)
    overall_shear_force_ += overall_shear_force_per_bubble_(i);
  overall_shear_stress_ = overall_shear_force_ / total_surface_;

}

void IJK_One_Dimensional_Subproblems::compute_overall_nusselt_numbers_bubbles()
{
  overall_nusselt_number_ = 0.;
  interfacial_thermal_flux_ = 0.;
  total_surface_ = 0.;
  for (int i=0; i < nb_bubbles_; i++)
    {
      interfacial_thermal_flux_ += interfacial_thermal_flux_per_bubble_(i);
      total_surface_ += total_surface_per_bubble_(i);
      overall_nusselt_number_ += overall_nusselt_number_per_bubble_(i);
    }
}

void IJK_One_Dimensional_Subproblems::post_process_overall_bubbles_quantities(const int rank)
{
  if (Process::je_suis_maitre())
    {
      const int reset = 1;
      Nom probe_name = Nom("_thermal_rank_") + Nom(rank) + Nom("_thermal_subproblems") + ("_overall_bubbles_quantities_")
                       + Nom(ref_ijk_ft_->get_tstep()) + Nom(".out");
      Nom probe_header = Nom("tstep\ttime\tthermalrank"
                             "\tnusseltoverall");
      SFichier fic = Ouvrir_fichier(probe_name, probe_header, reset);
      if (is_updated_)
        {
          for (int i=0; i<subproblems_counter_; i++)
            {
              fic << ref_ijk_ft_->get_tstep() << " " << ref_ijk_ft_->get_current_time() << " ";
              fic << rank << " ";
              fic << overall_nusselt_number_ << " ";
              fic << finl;
            }
        }
      fic.close();
    }
}


