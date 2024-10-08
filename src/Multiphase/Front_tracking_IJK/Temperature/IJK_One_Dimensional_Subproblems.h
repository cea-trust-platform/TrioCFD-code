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
// File      : IJK_One_Dimensional_Subproblems.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_One_Dimensional_Subproblems_included
#define IJK_One_Dimensional_Subproblems_included

#include <IJK_One_Dimensional_Subproblem.h>
#include <IJK_Finite_Difference_One_Dimensional_Matrix_Assembler.h>
#include <TRUSTList.h>
#include <TRUST_List.h>


/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_One_Dimensional_Subproblems
//
// <Description of class IJK_One_Dimensional_Subproblems>
//
/////////////////////////////////////////////////////////////////////////////

class IJK_FT_double;
class Switch_FT_double;
class IJK_Thermal_base;
class IJK_Thermal_Subresolution;

class IJK_One_Dimensional_Subproblems : public LIST(IJK_One_Dimensional_Subproblem)
{

  Declare_instanciable(IJK_One_Dimensional_Subproblems);
  friend class IJK_One_Dimensional_Subproblem;

public :
  IJK_One_Dimensional_Subproblems(const IJK_FT_double& ijk_ft);
  void associer(const IJK_FT_double& ijk_ft) { ref_ijk_ft_ = ijk_ft; };
  void initialise_thermal_subproblems_list_params(const int& pre_initialise_thermal_subproblems_list,
                                                  const double& pre_factor_subproblems_number,
                                                  const int& remove_append_subproblems,
                                                  const int& use_sparse_matrix);
  void associate_variables_for_post_processing(IJK_Thermal_Subresolution& ref_thermal_subresolution);
  void set_max_subproblems(const int max_subproblems) { max_subproblems_ = max_subproblems; };
  void clean();
  void clean(int add, int append);
  void clean_add();
  void clean_remove();
  void clean_append();
  void complete_subproblems();
  void shorten_subproblems();
  void add_subproblems(int n);
  void compute_global_indices();
  void set_global_index();
  void associate_subproblem_to_ijk_indices(const int& i,
                                           const int& j,
                                           const int& k);
  void get_ijk_indices_from_subproblems(const int& rank,
                                        int& i,
                                        int& j,
                                        int& k);
  void associate_ijk_indices_to_subproblem(const int& rank,
                                           const int& i,
                                           const int& j,
                                           const int& k);
  int get_subproblem_index_from_ijk_indices(const int& i,
                                            const int& j,
                                            const int& k) const;
  void associate_sub_problem_to_inputs(IJK_Thermal_Subresolution& ref_thermal_subresolution,
                                       int i, int j, int k,
                                       const double& indicator,
                                       double global_time_step,
                                       double current_time,
                                       const IJK_Interfaces& interfaces,
                                       const FixedVector<IJK_Field_double, 3>& velocity,
                                       const FixedVector<IJK_Field_double, 3>& velocity_ft,
                                       const IJK_Field_double& pressure);
  void store_previous_temperature_indicator_velocities();
  void share_previous_temperature_indicator_velocities();
  void retrieve_boundary_previous_values();
  void share_boundary_previous_values();
  void complete_boundary_previous_values();
  int is_in_map_index_ijk(const std::map<int, std::map<int, std::map<int, int>>>& subproblem_to_ijk_indices,
                          const int& index_i,
                          const int& index_j,
                          const int& index_k);
  void set_effective_subproblems(const int& enable_probe_collision_detection);
  void interpolate_indicator_on_probes();
  void clear_sort_problems_colliding_bubbles();
  void interpolate_project_velocities_on_probes();
  void reajust_probes_length();
  void reajust_probes_length(const int probe_length_condition);
  void compute_modified_probe_length(const int& probe_variations_enabled);
  void compute_radial_convection_diffusion_operators();
  void compute_source_terms_impose_boundary_conditions(const int& boundary_condition_interface,
                                                       const double& interfacial_boundary_condition_value,
                                                       const int& impose_boundary_condition_interface_from_simulation,
                                                       const int& boundary_condition_end,
                                                       const double& end_boundary_condition_value,
                                                       const int& impose_user_boundary_condition_end_value);
  void approximate_temperature_increment_material_derivative();
  void retrieve_radial_quantities();
  void retrieve_temperature_solutions();
  void compute_local_temperature_gradient_solutions();
  void compute_local_velocity_gradient();
  void get_subproblem_ijk_indices(int& i, int& j, int& k, int& subproblem_index) const;
  const int& get_dxyz_increment_bool(const int& subproblem_index) const;
  const int& get_dxyz_over_two_increment_bool(const int& subproblem_index) const;
  const FixedVector<int,3>& get_pure_neighbours_corrected_sign(const int& subproblem_index) const;
  const std::vector<std::vector<std::vector<bool>>>& get_pure_neighbours_to_correct(const int& subproblem_index) const;
  const std::vector<std::vector<std::vector<double>>>& get_pure_neighbours_corrected_distance(const int& subproblem_index) const;
  const std::vector<std::vector<std::vector<double>>>& get_pure_neighbours_corrected_colinearity(const int& subproblem_index) const;
  const std::vector<std::vector<std::vector<std::vector<bool>>>> get_pure_neighbours_last_faces_to_correct(const int& subproblem_index) const;
  const std::vector<std::vector<std::vector<std::vector<double>>>> get_pure_neighbours_last_faces_corrected_distance(const int& subproblem_index) const;
  const std::vector<std::vector<std::vector<std::vector<double>>>> get_pure_neighbours_last_faces_corrected_colinearity(const int& subproblem_index) const;
  double get_interfacial_gradient_corrected(int i);
  double get_temperature_profile_at_point(const int& i, const double& dist) const;
  const double& get_dist_cell_interface(const int& i) const;
  const FixedVector<double,6>& get_dist_faces_interface(const int& i) const;
  const Vecteur3& get_bary_facet(const int& i) const;
  double get_temperature_times_velocity_profile_at_point(const int& i, const double& dist, const int& dir) const;
  DoubleVect get_temperature_profile_discrete_integral_at_point(const int& i, const double& dist, const int& level, const int& dir) const;
  DoubleVect get_temperature_times_velocity_profile_discrete_integral_at_point(const int& i, const double& dist, const int& level, const int& dir) const;
  DoubleVect get_temperature_gradient_profile_discrete_integral_at_point(const int& i, const double& dist, const int& level, const int& dir) const;
  DoubleVect get_temperature_gradient_times_conductivity_profile_discrete_integral_at_point(const int& i, const double& dist, const int& level, const int& dir) const;
  double get_temperature_gradient_profile_at_point(const int& i, const double& dist, const int& dir) const;
  double get_temperature_gradient_times_conductivity_profile_at_point(const int& i, const double& dist, const int& dir) const;

  Nom get_header_from_string_lists(const std::vector<std::string>& key_results_int,
                                   const std::vector<std::string>& key_results_double);
  void set_results_probes_size(const std::vector<std::string>& key_results_int,
                               const std::vector<std::string>& key_results_double,
                               std::map<std::string, ArrOfInt>& results_probes_int,
                               std::map<std::string, ArrOfDouble>& results_probes_double);
  void set_results_probes_fic(SFichier& fic,
                              const std::vector<std::string>& key_results_int,
                              const std::vector<std::string>& key_results_double,
                              std::map<std::string, ArrOfInt>& results_probes_int,
                              std::map<std::string, ArrOfDouble>& results_probes_double);
  void thermal_subresolution_outputs_parallel(const int& rank,
                                              const Nom& interfacial_quantities_thermal_probes,
                                              const Nom& overall_bubbles_quantities,
                                              const Nom& local_quantities_thermal_probes_time_index_folder);
  void thermal_subresolution_outputs(const int& rank,
                                     const Nom& interfacial_quantities_thermal_probes,
                                     const Nom& overall_bubbles_quantities,
                                     const Nom& local_quantities_thermal_probes_time_index_folder);

  const int& get_subproblems_counter() const
  {
    return subproblems_counter_;
  }
  const int& get_effective_subproblems_counter() const
  {
    return effective_subproblems_counter_;
  }

  double get_min_temperature() const;
  double get_max_temperature() const;
  double get_min_temperature_domain_ends() const;
  double get_max_temperature_domain_ends() const;
  double get_min_euler_time_step(int& nb_iter_explicit);
  double get_local_max_fourier_time_step_probe_length();
  double get_local_min_fourier_time_step_probe_length();
  double get_local_max_cfl_time_step_probe_length();
  double get_local_min_cfl_time_step_probe_length();
  double get_local_dt_cfl();
  double get_local_dt_cfl_min_delta_xyz();
  int get_probe_variations_enabled(const int& probe_variations_priority);
  int get_probe_variations_enabled_priority();
  int get_probe_variations_enabled_non_priority();
  void set_local_time_step(const double& local_time_step);
  void prepare_temporal_schemes();
  const int& get_end_index_subproblem(const int index) const;
  void post_processed_all_probes();
  void sort_limited_probes_spherical_coords_post_processing(const int& post_process_all_probes,
                                                            const int& nb_theta, const int& nb_phi,
                                                            const int theta_diag_val, const int phi_diag_val);
  void compute_overall_quantities_per_bubbles(const IJK_Field_double& temperature_ghost,
                                              const double& delta_temperature,
                                              const double& lambda);
  void compute_nusselt_numbers_per_bubbles(const IJK_Field_double& temperature_gradient_ghost,
                                           const double& delta_temperature,
                                           const double& lambda);
  void compute_shear_per_bubbles();
  void compute_overall_bubbles_quantities(IJK_Thermal_Subresolution& ref_thermal_subresolution);
  void compute_dynamics_per_bubbles();
  void compute_overall_quantities();
  void post_process_overall_bubbles_quantities(const int rank, const Nom& overall_bubbles_quantities);


protected :
  std::vector<IJK_One_Dimensional_Subproblem*> one_dimensional_effective_subproblems_;
  std::vector<IJK_One_Dimensional_Subproblem*> one_dimensional_disabled_subproblems_;
  int init_ = 1;
  int debug_ = 0;
  int max_subproblems_ = 0;
  int subproblems_counter_ = 0;
  int effective_subproblems_counter_ = 0;
  int disabled_subproblems_counter_ = 0;
  int effective_and_disabled_subproblems_counter_ = 0;
  int global_subproblems_counter_ = 0;
  int index_ini_=0;
  int index_end_=0;
  int reallocate_subproblems_ = 1;
  bool is_updated_ = 0;
  OBS_PTR(IJK_FT_double) ref_ijk_ft_;

  FixedVector<ArrOfInt, 3> ijk_indices_to_subproblem_;
  std::map<int, std::map<int, std::map<int, int>>> subproblem_to_ijk_indices_previous_;
  std::map<int, std::map<int, std::map<int, int>>> subproblem_to_ijk_indices_;
  std::vector<DoubleVect> temperature_probes_previous_;
  std::vector<double> indicator_probes_previous_;
  std::vector<Vecteur3> velocities_probes_previous_;
  std::vector<Vecteur3> normal_vector_compo_probes_previous_;

  int * points_per_thermal_subproblem_;
  std::map<int, std::map<int, std::map<int, int>>> subproblem_to_ijk_indices_previous_global_;
  std::map<int, std::map<int, std::map<int, int>>> subproblem_to_ijk_indices_global_;

  FixedVector<std::vector<ArrOfInt>,2> index_ij_subproblems_local_perio_;
  std::vector<std::vector<ArrOfDouble>> temperature_probes_previous_local_perio_;
  std::vector<ArrOfDouble> indicator_probes_previous_local_perio_;
  std::vector<std::vector<ArrOfDouble>> velocities_probes_previous_local_perio_;
  std::vector<std::vector<ArrOfDouble>> normal_vector_compo_probes_previous_local_perio_;

  FixedVector<std::vector<ArrOfInt>,2> index_ij_subproblems_global_;
  std::vector<std::vector<ArrOfDouble>> temperature_probes_previous_global_;
  std::vector<ArrOfDouble> indicator_probes_previous_global_;
  std::vector<std::vector<ArrOfDouble>> velocities_probes_previous_global_;
  std::vector<std::vector<ArrOfDouble>> normal_vector_compo_probes_previous_global_;

  int reference_gfm_on_probes_ = 0;
  int pre_initialise_thermal_subproblems_list_ = 0;
  double pre_factor_subproblems_number_ = 1.;
  int remove_append_subproblems_ = 0;
  int use_sparse_matrix_ = 0;

  double spherical_nusselt_ = 2.;
  double spherical_nusselt_liquid_ = 2.;
  double heat_flux_spherical_ = 0.;
  double overall_shear_stress_ = 0.;
  double overall_shear_force_ = 0.;

  double overall_nusselt_number_ = 0.;
  double overall_nusselt_number_gfm_ = 0.;
  double overall_nusselt_number_spherical_ = 0.;
  double overall_nusselt_number_liquid_ = 0.;
  double overall_nusselt_number_gfm_liquid_ = 0.;
  double overall_nusselt_number_spherical_liquid_ = 0.;

  double overall_nusselt_number_error_ = 0.;
  double overall_nusselt_number_gfm_error_ = 0.;
  double overall_nusselt_number_liquid_error_ = 0.;
  double overall_nusselt_number_gfm_liquid_error_ = 0.;

  double overall_nusselt_number_error_rel_ = 0.;
  double overall_nusselt_number_gfm_error_rel_ = 0.;
  double overall_nusselt_number_liquid_error_rel_ = 0.;
  double overall_nusselt_number_gfm_liquid_error_rel_ = 0.;

  double caracteristic_length_ = 0.;
  double caracteristic_length_from_surfaces_ = 0.;
  double caracteristic_length_from_volumes_ = 0.;

  double delta_temperature_ = -1;
  double mean_liquid_temperature_= -1.;

  double interfacial_thermal_flux_ = 0.;
  double interfacial_thermal_flux_gfm_ = 0.;

  double total_surface_ = 0.;
  double total_volume_ = 0.;
  double lambda_=0.;
  double * prandtl_number_;

  double radius_from_surfaces_ = 0.;
  double radius_from_volumes_ = 0.;

  double uniform_alpha_ = 0.;
  double error_temperature_ana_total_ = 0.;
  double error_temperature_ana_squared_total_ = 0.;
  double error_temperature_ana_rel_total_ = 0.;

  const ArrOfDouble * bubbles_volume_;
  const ArrOfDouble * bubbles_rising_velocities_;
  const DoubleTab * bubbles_rising_vectors_per_bubble_;
  ArrOfDouble bubbles_rising_relative_velocities_;
  ArrOfDouble bubbles_rising_relative_velocities_upstream_;
  ArrOfDouble bubbles_reynolds_numbers_per_bubble_;
  ArrOfDouble bubbles_reynolds_numbers_from_surface_per_bubble_;
  ArrOfDouble bubbles_reynolds_numbers_from_volume_per_bubble_;
  ArrOfDouble bubbles_reynolds_numbers_per_bubble_upstream_;
  ArrOfDouble bubbles_reynolds_numbers_from_surface_per_bubble_upstream_;
  ArrOfDouble bubbles_reynolds_numbers_from_volume_per_bubble_upstream_;
  FixedVector<ArrOfDouble, 3> relative_rising_dir_compo_;
  FixedVector<ArrOfDouble, 3> relative_rising_dir_compo_upstream_;

  double velocity_upstream_ = 0.;
  int gravity_dir_ = 0;
  Vecteur3 total_relative_rising_dir_compo_;
  Vecteur3 total_relative_rising_dir_compo_upstream_;
  Vecteur3 total_rising_dir_compo_;
  double bubbles_total_rising_velocities_ = 0.;
  double bubbles_total_rising_relative_velocities_ = 0.;
  double bubbles_total_rising_relative_velocities_upstream_ = 0.;
  double bubbles_total_reynolds_numbers_ = 0.;
  double bubbles_total_reynolds_numbers_from_surface_per_bubble_ = 0.;
  double bubbles_total_reynolds_numbers_from_volume_per_bubble_ = 0.;
  double bubbles_total_reynolds_numbers_upstream_ = 0.;
  double bubbles_total_reynolds_numbers_from_surface_per_bubble_upstream_ = 0.;
  double bubbles_total_reynolds_numbers_from_volume_per_bubble_upstream_ = 0.;

  ArrOfDouble radius_outputs_;
  ArrOfDouble theta_outputs_;
  ArrOfDouble phi_outputs_;
  ArrOfInt global_indices_post_processed_;

  ArrOfDouble interfacial_thermal_flux_per_bubble_ ;
  ArrOfDouble interfacial_thermal_flux_per_bubble_gfm_;
  ArrOfDouble interfacial_thermal_flux_per_bubble_spherical_;

  ArrOfDouble total_surface_per_bubble_;

  ArrOfDouble overall_nusselt_number_per_bubble_;
  ArrOfDouble overall_nusselt_number_per_bubble_gfm_;
  ArrOfDouble overall_nusselt_number_per_bubble_spherical_;
  ArrOfDouble overall_nusselt_number_per_bubble_liquid_;
  ArrOfDouble overall_nusselt_number_per_bubble_gfm_liquid_;
  ArrOfDouble overall_nusselt_number_per_bubble_spherical_liquid_;

  ArrOfDouble overall_shear_stress_per_bubble_;
  ArrOfDouble overall_shear_force_per_bubble_;

  ArrOfDouble radius_from_surfaces_per_bubble_;
  ArrOfDouble radius_from_volumes_per_bubble_;

  ArrOfDouble caracteristic_length_from_surfaces_per_bubble_;
  ArrOfDouble caracteristic_length_from_volumes_per_bubble_;

  ArrOfDouble overall_nusselt_number_per_bubble_error_;
  ArrOfDouble overall_nusselt_number_per_bubble_gfm_error_;
  ArrOfDouble overall_nusselt_number_per_bubble_liquid_error_;
  ArrOfDouble overall_nusselt_number_per_bubble_gfm_liquid_error_;

  ArrOfDouble overall_nusselt_number_per_bubble_error_rel_;
  ArrOfDouble overall_nusselt_number_per_bubble_gfm_error_rel_;
  ArrOfDouble overall_nusselt_number_per_bubble_liquid_error_rel_;
  ArrOfDouble overall_nusselt_number_per_bubble_gfm_liquid_error_rel_;

  ArrOfDouble bubbles_peclet_numbers_per_bubble_;
  ArrOfDouble bubbles_peclet_numbers_from_surface_per_bubble_;
  ArrOfDouble bubbles_peclet_numbers_from_volume_per_bubble_;
  ArrOfDouble bubbles_peclet_numbers_per_bubble_upstream_;
  ArrOfDouble bubbles_peclet_numbers_from_surface_per_bubble_upstream_;
  ArrOfDouble bubbles_peclet_numbers_from_volume_per_bubble_upstream_;
  double bubbles_total_peclet_numbers_ = 0.;
  double bubbles_total_peclet_numbers_from_surface_per_bubble_ = 0.;
  double bubbles_total_peclet_numbers_from_volume_per_bubble_ = 0.;
  double bubbles_total_peclet_numbers_upstream_  = 0.;
  double bubbles_total_peclet_numbers_from_surface_per_bubble_upstream_ = 0.;
  double bubbles_total_peclet_numbers_from_volume_per_bubble_upstream_  = 0.;

  const Vecteur3 * liquid_velocity_;

  int nb_bubbles_ = 0;
  int * latastep_reprise_;

  ArrOfInt collision_indices_;

};

#endif /* IJK_One_Dimensional_Subproblems_included */
