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
#include <IJK_Bubble_tools.h>
#include <IJK_Thermal_base.h>
#include <IJK_Thermal_Subresolution.h>

Implemente_instanciable_sans_constructeur(IJK_One_Dimensional_Subproblems, "IJK_One_Dimensional_Subproblems", LIST(IJK_One_Dimensional_Subproblem));

IJK_One_Dimensional_Subproblems::IJK_One_Dimensional_Subproblems()
{
  interfacial_thermal_flux_per_bubble_.set_smart_resize(1);
  interfacial_thermal_flux_per_bubble_gfm_.set_smart_resize(1);
  interfacial_thermal_flux_per_bubble_spherical_.set_smart_resize(1);
  total_surface_per_bubble_.set_smart_resize(1);
  overall_nusselt_number_per_bubble_.set_smart_resize(1);
  overall_nusselt_number_per_bubble_gfm_.set_smart_resize(1);
  overall_nusselt_number_per_bubble_spherical_.set_smart_resize(1);
  overall_nusselt_number_per_bubble_liquid_.set_smart_resize(1);
  overall_nusselt_number_per_bubble_gfm_liquid_.set_smart_resize(1);
  overall_nusselt_number_per_bubble_spherical_liquid_.set_smart_resize(1);
  overall_shear_stress_per_bubble_.set_smart_resize(1);
  overall_shear_force_per_bubble_.set_smart_resize(1);
  radius_outputs_.set_smart_resize(1);
  theta_outputs_.set_smart_resize(1);
  phi_outputs_.set_smart_resize(1);
  global_indices_post_processed_.set_smart_resize(1);
  radius_from_surfaces_per_bubble_.set_smart_resize(1);
  radius_from_volumes_per_bubble_.set_smart_resize(1);
  caracteristic_length_from_surfaces_per_bubble_.set_smart_resize(1);
  caracteristic_length_from_volumes_per_bubble_.set_smart_resize(1);
  bubbles_volume_ = nullptr;
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
                                                                                 const double& pre_factor_subproblems_number,
                                                                                 const int& remove_append_subproblems)
{
  pre_initialise_thermal_subproblems_list_ = pre_initialise_thermal_subproblems_list;
  pre_factor_subproblems_number_ = pre_factor_subproblems_number;
  remove_append_subproblems_ = remove_append_subproblems;
}

void IJK_One_Dimensional_Subproblems::clean()
{
  clean(pre_initialise_thermal_subproblems_list_, remove_append_subproblems_);
}

void IJK_One_Dimensional_Subproblems::clean(int add, int append)
{

  if (append)
    clean_append();
  else if (add)
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

void IJK_One_Dimensional_Subproblems::clean_append()
{
  shorten_subproblems();
}

void IJK_One_Dimensional_Subproblems::complete_subproblems()
{
  if (init_)
    {
      int total_subproblems = subproblems_counter_;
      total_subproblems = Process::mp_sum(total_subproblems);
      max_subproblems_ = (int) (pre_factor_subproblems_number_ * total_subproblems);
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

void IJK_One_Dimensional_Subproblems::shorten_subproblems()
{
  const int active_subproblems = subproblems_counter_;
  const int total_subproblems = (*this).size();
  for (int i=0; i<(total_subproblems - active_subproblems); i++)
    (*this).suppr((*this).dernier());
  assert((*this).size() == subproblems_counter_);
  init_ = 0;
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

void IJK_One_Dimensional_Subproblems::compute_global_indices()
{
  global_subproblems_counter_ = subproblems_counter_;
  global_subproblems_counter_ = Process::mp_sum(global_subproblems_counter_);
  const int proc_number = Process::nproc();
  ArrOfInt indices(proc_number);
  const int my_process_number = Process::me();
  indices(my_process_number) = subproblems_counter_;
  mp_sum_for_each_item(indices);
  ArrOfInt indices_ini(proc_number);
  ArrOfInt indices_max(proc_number);
  indices_ini(0) = 0;
  indices_max(0) = indices(0);
  for (int i=1; i<proc_number; i++)
    {
      indices_ini(i) = indices_ini(i-1) + indices(i-1);
      indices_max(i) = indices_max(i-1) + indices(i);
    }
  assert(indices_max(proc_number - 1) == global_subproblems_counter_);
  index_ini_ = indices_ini(my_process_number);
  index_end_ = indices_max(my_process_number);
  set_global_index();
}

void IJK_One_Dimensional_Subproblems::set_global_index()
{
  int index_counter = index_ini_;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      (*this)[itr].set_global_index(index_counter);
      index_counter++;
    }
}

void IJK_One_Dimensional_Subproblems::associate_variables_for_post_processing(IJK_Thermal_Subresolution& ref_thermal_subresolution)
{
  bubbles_volume_ = &ref_thermal_subresolution.bubbles_volume_;
}

void IJK_One_Dimensional_Subproblems::associate_sub_problem_to_inputs(IJK_Thermal_Subresolution& ref_thermal_subresolution,
                                                                      int i, int j, int k,
                                                                      const double& indicator,
                                                                      double global_time_step,
                                                                      double current_time,
                                                                      const IJK_Interfaces& interfaces,
                                                                      const FixedVector<IJK_Field_double, 3>& velocity,
                                                                      const FixedVector<IJK_Field_double, 3>& velocity_ft,
                                                                      const IJK_Field_double& pressure)
{
  if (!init_ && subproblems_counter_ > max_subproblems_ && (pre_initialise_thermal_subproblems_list_ && !remove_append_subproblems_))
    {
      Cerr << max_subproblems_ << "subproblems were expected but" << subproblems_counter_ << "subproblem try to be associated with the list of subproblems";
      Process::exit();
    }

  debug_ = ref_thermal_subresolution.debug_;

  bubbles_volume_ = &ref_thermal_subresolution.bubbles_volume_;

  ArrOfDouble bubble_rising_vector(3);
  ArrOfDouble normal_vector(3);
  ArrOfDouble facet_barycentre(3);
  ArrOfDouble bubble_barycentre(3);

  if (debug_)
    Cerr << "Mixed cell indices (i,j,k) : (" << i << ";" << j << ";" << k << ")" << finl;
  const int compo_connex = (*(ref_thermal_subresolution.eulerian_compo_connex_from_interface_int_ns_))(i, j, k);
  const int compo_connex_ghost = (*(ref_thermal_subresolution.eulerian_compo_connex_from_interface_int_ns_))(i,j,k);
  if (debug_)
    {
      Cerr << "compo_connex : " << compo_connex << finl;
      Cerr << "compo_connex_ghost : " << compo_connex_ghost << finl;
      Cerr << "bubbles_barycentre : " << ref_thermal_subresolution.bubbles_barycentre_ << finl;
    }

  // Need for a Navier-Stokes field (NOT FT)
  const double distance = ref_thermal_subresolution.eulerian_distance_ns_(i, j ,k);
  const double curvature = ref_thermal_subresolution.eulerian_curvature_ns_(i, j ,k);
  const double interfacial_area = ref_thermal_subresolution.eulerian_interfacial_area_ns_(i, j ,k);

  IJK_Splitting splitting = (*(ref_thermal_subresolution.eulerian_compo_connex_from_interface_int_ns_)).get_splitting();
  const double bubble_rising_velocity = ref_thermal_subresolution.rising_velocities_(compo_connex);
  //  const double bubble_rising_velocity = rising_velocities(compo_connex);
  for (int dir=0; dir < 3; dir++)
    {
      facet_barycentre(dir) = ref_thermal_subresolution.eulerian_facets_barycentre_ns_[dir](i, j, k);
      normal_vector(dir) = ref_thermal_subresolution.eulerian_normal_vectors_ns_[dir](i, j, k);
      bubble_rising_vector(dir) = ref_thermal_subresolution.rising_vectors_(compo_connex, dir);
      bubble_barycentre(dir) = ref_thermal_subresolution.bubbles_barycentre_(compo_connex_ghost, dir);
    }

  const int total_subproblems = (*this).size();
  if (subproblems_counter_ >= total_subproblems)
    {
      IJK_One_Dimensional_Subproblem subproblem(ref_ijk_ft_);
      (*this).add(subproblem);
      init_ = 1;
    }
  (*this)[subproblems_counter_].associate_sub_problem_to_inputs(ref_thermal_subresolution,
                                                                i, j, k,
                                                                init_,
                                                                subproblems_counter_,
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
                                                                indicator,
                                                                interfaces,
                                                                velocity,
                                                                velocity_ft,
                                                                pressure);

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

Nom IJK_One_Dimensional_Subproblems::get_header_from_string_lists(const std::vector<std::string>& key_results_int,
                                                                  const std::vector<std::string>& key_results_double)
{
  Nom probe_header;
  int i;
  probe_header = key_results_int[0];
  const int size_int = (int) key_results_int.size();
  const int size_double = (int) key_results_double.size();
  for (i=1; i<size_int; i++)
    probe_header += ("\t" + key_results_int[i]);
  for (i=0; i<size_double; i++)
    probe_header += ("\t" + key_results_double[i]);
  return probe_header;
}

void IJK_One_Dimensional_Subproblems::set_results_probes_size(const std::vector<std::string>& key_results_int,
                                                              const std::vector<std::string>& key_results_double,
                                                              std::map<std::string, ArrOfInt>& results_probes_int,
                                                              std::map<std::string, ArrOfDouble>& results_probes_double)
{
  // Below 1000 per simu ? * nb_procs ? Is it a problem ?
  const int size_outputs = global_subproblems_counter_;
  const int size_int = (int) key_results_int.size();
  const int size_double = (int) key_results_double.size();
  int i;
  for (i=0; i<size_int; i++)
    results_probes_int[key_results_int[i]] = ArrOfInt(size_outputs);
  for (i=0; i<size_double; i++)
    results_probes_double[key_results_double[i]] = ArrOfDouble(size_outputs);
}

void IJK_One_Dimensional_Subproblems::set_results_probes_fic(SFichier& fic,
                                                             const std::vector<std::string>& key_results_int,
                                                             const std::vector<std::string>& key_results_double,
                                                             std::map<std::string, ArrOfInt>& results_probes_int,
                                                             std::map<std::string, ArrOfDouble>& results_probes_double)
{
  const int size_outputs = global_subproblems_counter_;
  const int size_int = (int) results_probes_int.size();
  const int size_double = (int) results_probes_double.size();
  int i,j;
  for (j=0; j<size_outputs; j++)
    {
      for (i=0; i<size_int; i++)
        fic << results_probes_int[key_results_int[i]](j) << " ";
      for (i=0; i<size_double-1; i++)
        fic << results_probes_double[key_results_double[i]](j) << " ";
      fic << results_probes_double[key_results_double[size_double-1]](j) << finl;
    }
}

void IJK_One_Dimensional_Subproblems::thermal_subresolution_outputs_parallel(const int& rank,
                                                                             const Nom& interfacial_quantities_thermal_probes,
                                                                             const Nom& overall_bubbles_quantities,
                                                                             const Nom& local_quantities_thermal_probes_time_index_folder)
{

  std::map<std::string, ArrOfDouble> results_probes_double;
  std::map<std::string, ArrOfInt> results_probes_int;

  std::vector<std::string> key_results_int = {"tstep", "thermalrank", "postproindex", "globalsubproblem", "localsubproblem"};
  std::vector<std::string> key_results_double = {"time",
                                                 "nx", "ny", "nz",
                                                 "t1x", "t1y", "t1z", "t2x", "t2y", "t2z",
                                                 "s1x", "s1y", "s1z", "s2x", "s2y", "s2z",
                                                 "rsph", "thetasph", "phisph",
                                                 "temperature_interp","temperature_solution","temperature_gradient", "temperature_gradient_sol",
                                                 "temperature_double_deriv_sol",
                                                 "temperature_gradient_tangential","temperature_gradient_tangential2",
                                                 "temperature_gradient_tangential_rise","temperature_gradient_azymuthal",
                                                 "temperature_diffusion_hessian_cartesian_trace",
                                                 "temperature_diffusion_hessian_trace",
                                                 "radial_temperature_diffusion",
                                                 "tangential_temperature_diffusion",
                                                 "surface","thermal_flux","lambda","alpha","prandtl_liq",
                                                 "shear","force",
                                                 "pressure",
                                                 "u_x","u_y","u_z",
                                                 "u_r","u_r_corr","u_r_static","u_r_advected",
                                                 "u_theta","u_theta_corr","u_theta_static","u_theta_advected",
                                                 "u_theta2","u_theta2_corr","u_theta2_static","u_theta2_advected",
                                                 "u_theta_rise","u_theta_rise_corr","u_theta_rise_static","u_theta_rise_advected",
                                                 "u_phi","u_phi_corr","u_phi_static","u_phi_advected",
                                                 "du_r_dr","du_theta_dr","du_theta2_dr","du_theta_rise_dr","du_phi_dr"
                                                };

  Nom probe_header = get_header_from_string_lists(key_results_int, key_results_double);
  set_results_probes_size(key_results_int,
                          key_results_double,
                          results_probes_int,
                          results_probes_double);
  // Common to all procs
  int i;
  const int size_outputs = global_subproblems_counter_;
  const int size_int = (int) key_results_int.size();
  const int size_double = (int) key_results_double.size();


  for (int itr=0; itr < size_outputs; itr++)
    {
      // const int global_index_post_processed = global_indices_post_processed_(itr);
      if (itr >= index_ini_ && itr < index_end_)
        {

          Cerr << "Post-process this probe on proc:" << Process::me() << finl;
          (*this)[itr - index_ini_].retrieve_interfacial_quantities(rank,
                                                                    itr,
                                                                    key_results_int,
                                                                    key_results_double,
                                                                    results_probes_int,
                                                                    results_probes_double);
        }
    }
  for (i=0; i<size_int; i++)
    {
      ArrOfInt& array_int_tmp = results_probes_int[key_results_int[i]];
      mp_sum_for_each_item(array_int_tmp);
    }
  for (i=0; i<size_double; i++)
    {
      ArrOfDouble& array_double_tmp = results_probes_double[key_results_double[i]];
      mp_sum_for_each_item(array_double_tmp);
    }
  /*
   * Post-process all probes for interfacial quantities
   */
  const int reset = 1;
  const int last_time_index = ref_ijk_ft_->get_tstep();
  const int max_digit = 3;
  const int max_digit_time = 8;
  const int max_rank_digit = rank < 1 ? 1 : (int) (log10(rank) + 1);
  const int nb_digit_tstep = last_time_index < 1 ? 1 : (int) (log10(last_time_index) + 1);

  Nom probe_name = Nom("_thermal_rank_") +  Nom(std::string(max_digit - max_rank_digit, '0'))  + Nom(rank)
                   + Nom("_thermal_subproblems_interfacial_quantities_time_index_")
                   + Nom(std::string(max_digit_time - nb_digit_tstep, '0')) + Nom(last_time_index) + Nom(".out");

  if (Process::je_suis_maitre())
    {
      SFichier fic = Open_file_folder(interfacial_quantities_thermal_probes, probe_name, probe_header, reset);
      set_results_probes_fic(fic,
                             key_results_int,
                             key_results_double,
                             results_probes_int,
                             results_probes_double);
      fic.close();
    }

  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].thermal_subresolution_outputs_parallel(rank, local_quantities_thermal_probes_time_index_folder);

  post_process_overall_bubbles_quantities(rank, overall_bubbles_quantities);
}

void IJK_One_Dimensional_Subproblems::thermal_subresolution_outputs(const int& rank,
                                                                    const Nom& interfacial_quantities_thermal_probes,
                                                                    const Nom& overall_bubbles_quantities,
                                                                    const Nom& local_quantities_thermal_probes_time_index_folder)
{
  /*
   * Replace routines for parallel calculation
   */
  Cerr << "Post-processing on the probes" << finl;
  //if (Process::je_suis_maitre())
  // {
  const int reset = 1;
  const int last_time_index = ref_ijk_ft_->get_tstep();
  Nom probe_header = Nom("tstep\tthermalrank\tpostproindex\tglobalsubproblem\tlocalsubproblem\ttime"
                         "\tnx\tny\tnz\tt1x\tt1y\tt2z\tt2x\tt2y\tt2z\ts1x\ts1y\ts1z\ts2x\ts2y\ts2z"
                         "\trsph\tthetasph\tphisph"
                         "\ttemperature_interp\ttemperature_solution"
                         "\ttemperature_gradient\ttemperature_gradient_sol"
                         "\ttemperature_double_deriv_sol"
                         "\ttemperature_gradient_tangential\ttemperature_gradient_tangential2"
                         "\ttemperature_gradient_tangential_rise\ttemperature_gradient_azymuthal"
                         "\ttemperature_diffusion_hessian_cartesian_trace"
                         "\ttemperature_diffusion_hessian_trace"
                         "\tradial_temperature_diffusion"
                         "\ttangential_temperature_diffusion"
                         "\tsurface\tthermal_flux\tlambda\talpha\tprandtl_liq"
                         "\tshear\tforce"
                         "\tpressure"
                         "\tu_x\tu_y\tu_z"
                         "\tu_r\tu_r_corr\tu_r_static\tu_r_advected"
                         "\tu_theta\tu_theta_corr\tu_theta_static\tu_theta_advected"
                         "\tu_theta2\tu_theta2_corr\tu_theta2_static\tu_theta2_advected"
                         "\tu_theta_rise\tu_theta_rise_corr\tu_theta_rise_static\tu_theta_rise_advected"
                         "\tu_phi\tu_phi_corr\tu_phi_static\tu_phi_advected"
                         "\tdu_r_dr\tdu_theta_dr\tdu_theta2_dr\tdu_theta_rise_dr\tdu_phi_dr");

  const int max_digit = 3;
  const int max_digit_time = 8;
  const int max_rank_digit = rank < 1 ? 1 : (int) (log10(rank) + 1);
  const int nb_digit_tstep = last_time_index < 1 ? 1 : (int) (log10(last_time_index) + 1);

  Nom probe_name = Nom("_thermal_rank_") +  Nom(std::string(max_digit - max_rank_digit, '0'))  + Nom(rank)
                   + Nom("_thermal_subproblems_interfacial_quantities_time_index_")
                   + Nom(std::string(max_digit_time - nb_digit_tstep, '0')) + Nom(last_time_index) + Nom(".out");

  const int proc_number = Process::nproc();
  if (proc_number != 1)
    {
      const int my_process_number = Process::me();
      Nom my_process_string = Nom(".processor_") + Nom(my_process_number);
      probe_name += my_process_string;
    }

  /*
   * Post-process all probes for interfacial quantities
   */
  SFichier fic = Open_file_folder(interfacial_quantities_thermal_probes, probe_name, probe_header, reset);

  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].thermal_subresolution_outputs(fic, rank, local_quantities_thermal_probes_time_index_folder);
  fic.close();

  post_process_overall_bubbles_quantities(rank, overall_bubbles_quantities);
}

double IJK_One_Dimensional_Subproblems::get_min_temperature() const
{
  double min_temperature = 1e20;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      min_temperature = std::min(min_temperature, (*this)[itr].get_min_temperature());
    }
  min_temperature = Process::mp_min(min_temperature);
  return min_temperature;
}

double IJK_One_Dimensional_Subproblems::get_max_temperature() const
{
  double max_temperature = -1e20;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      max_temperature = std::max(max_temperature, (*this)[itr].get_max_temperature());
    }
  max_temperature = Process::mp_max(max_temperature);
  return max_temperature;
}

double IJK_One_Dimensional_Subproblems::get_min_temperature_domain_ends() const
{
  double min_temperature = -1e20;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      min_temperature = std::min(min_temperature, (*this)[itr].get_min_temperature_domain_ends());
    }
  min_temperature = Process::mp_min(min_temperature);
  return min_temperature;
}

double IJK_One_Dimensional_Subproblems::get_max_temperature_domain_ends() const
{
  double max_temperature = -1e20;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      max_temperature = std::max(max_temperature, (*this)[itr].get_max_temperature_domain_ends());
    }
  max_temperature = Process::mp_max(max_temperature);
  return max_temperature;
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
  min_euler_time_step = Process::mp_min(min_euler_time_step);
  nb_iter_explicit = Process::mp_max(nb_iter_explicit);
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
  max_local_fourier_time_step_probe_length = Process::mp_max(max_local_fourier_time_step_probe_length);
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
  max_local_cfl_time_step_probe_length = Process::mp_max(max_local_cfl_time_step_probe_length);
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
  min_local_fourier_time_step_probe_length = Process::mp_min(min_local_fourier_time_step_probe_length);
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
  min_local_cfl_time_step_probe_length = Process::mp_min(min_local_cfl_time_step_probe_length);
  return min_local_cfl_time_step_probe_length;
}

double IJK_One_Dimensional_Subproblems::get_local_dt_cfl()
{
  double min_local_dt_cfl = 1.e20;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      min_local_dt_cfl = std::min(min_local_dt_cfl, (*this)[itr].get_local_dt_cfl());
    }
  min_local_dt_cfl = Process::mp_min(min_local_dt_cfl);
  return min_local_dt_cfl;
}

double IJK_One_Dimensional_Subproblems::get_local_dt_cfl_min_delta_xyz()
{
  double min_local_dt_cfl_min_delta_xyz = 1.e20;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      min_local_dt_cfl_min_delta_xyz = std::min(min_local_dt_cfl_min_delta_xyz, (*this)[itr].get_local_dt_cfl_min_delta_xyz());
    }
  min_local_dt_cfl_min_delta_xyz = Process::mp_min(min_local_dt_cfl_min_delta_xyz);
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
    (*this)[itr].set_post_processing_theta_phi_scope(0);
}

void IJK_One_Dimensional_Subproblems::sort_limited_probes_spherical_coords_post_processing(const int& post_process_all_probes,
                                                                                           const int& nb_theta, const int& nb_phi,
                                                                                           const int theta_diag_val, const int phi_diag_val)
{
  if (post_process_all_probes)
    post_processed_all_probes();
  else
    {
      const int nb_subproblems_total = global_subproblems_counter_;
      Cerr << global_subproblems_counter_ << finl;

      ArrOfDouble r_sph(nb_subproblems_total);
      ArrOfDouble theta_sph(nb_subproblems_total);
      ArrOfDouble phi_sph(nb_subproblems_total);
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
          r_sph(i + index_ini_) = (*this)[i].get_radius_spherical_coords();
          theta_sph(i + index_ini_) = (*this)[i].get_theta_spherical_coords();
          phi_sph(i + index_ini_) = (*this)[i].get_phi_spherical_coords();
        }

      mp_sum_for_each_item(r_sph);
      mp_sum_for_each_item(theta_sph);
      mp_sum_for_each_item(phi_sph);

      int nb_outputs = ((nb_phi_even * nb_theta_even) > nb_subproblems_total) ? nb_subproblems_total : (nb_phi_even * nb_theta_even);
      nb_theta_even = (nb_outputs == nb_subproblems_total) ? (int) sqrt(2 * nb_subproblems_total) / 2 : nb_theta_even;
      nb_phi_even = (nb_outputs == nb_subproblems_total) ? (int) sqrt(2 * nb_subproblems_total) : nb_phi_even;
      nb_outputs = nb_phi_even * nb_theta_even;
      double theta_incr, phi_incr;
      theta_incr = (double) (M_PI / (double) nb_theta_even);
      phi_incr = (double) ((2 * M_PI) / (double) nb_phi_even);
      const double atan_theta_incr_ini = M_PI / 2;
      const double atan_incr_factor = -1;
      const double atan_phi_incr_ini = M_PI;
      // PI/2 -> -PI/2
      if (theta_diag_val)
        for (i=0; i<nb_theta_even; i++)
          theta_scope.append_array((theta_incr * (i + 0.5) - atan_theta_incr_ini) * atan_incr_factor);
      else
        for (i=0; i<nb_theta_even; i++)
          theta_scope.append_array((theta_incr * i - atan_theta_incr_ini) * atan_incr_factor);
      if (phi_diag_val)
        for (i=0; i<nb_phi_even; i++)
          phi_scope.append_array(phi_incr * (i + 0.5) - atan_phi_incr_ini);
      else
        for (i=0; i<nb_phi_even; i++)
          phi_scope.append_array(phi_incr * i - atan_phi_incr_ini);
      /*
       * Sort by phi and theta simultaneously
       */
      radius_outputs_.resize(nb_outputs);
      theta_outputs_.resize(nb_outputs);
      phi_outputs_.resize(nb_outputs);
      global_indices_post_processed_.resize(nb_outputs);
      int phi_theta_counter = 0;
      for (j=0; j<nb_phi_even; j++)
        for (i=0; i<nb_theta_even; i++)
          {
            ArrOfDouble theta_diff = theta_sph;
            ArrOfDouble phi_diff = phi_sph;
            std::vector<double> sum_errors_theta_phi_scope;
            theta_diff -= theta_scope(i);
            phi_diff -= phi_scope(j);
            for (k=0; k<nb_subproblems_total; k++)
              {
                theta_diff(k) = abs(theta_diff(k));
                phi_diff(k) = abs(phi_diff(k));
                sum_errors_theta_phi_scope.push_back(theta_diff(k));
                sum_errors_theta_phi_scope[k] += phi_diff(k);
              }


            const int theta_phi_scope_index = (int) std::distance(sum_errors_theta_phi_scope.begin(),
                                                                  std::min_element(sum_errors_theta_phi_scope.begin(),
                                                                                   sum_errors_theta_phi_scope.end()));

//            if (theta_phi_scope_index >= index_ini_ && theta_phi_scope_index < index_end_)
//              (*this)[theta_phi_scope_index - index_ini_].set_post_processing_theta_phi_scope(phi_theta_counter);
            radius_outputs_(phi_theta_counter) = r_sph(theta_phi_scope_index);
            theta_outputs_(phi_theta_counter) = theta_sph(theta_phi_scope_index);
            phi_outputs_(phi_theta_counter) = phi_sph(theta_phi_scope_index);
            global_indices_post_processed_(phi_theta_counter) = theta_phi_scope_index;
            phi_theta_counter ++;
          }

      // Cerr << "test" << finl;
      //      mp_sum_for_each_item(radius_outputs_);
      //      mp_sum_for_each_item(theta_outputs_);
      //      mp_sum_for_each_item(phi_outputs_);
      //      mp_sum_for_each_item(global_indices_post_processed_);

      Cerr << radius_outputs_ << finl;
      Cerr << theta_outputs_ << finl;
      Cerr << phi_outputs_ << finl;
      Cerr << global_indices_post_processed_ << finl;

      ArrOfDouble radius_outputs_tmp = radius_outputs_;
      ArrOfDouble theta_outputs_tmp = theta_outputs_;
      ArrOfDouble phi_outputs_tmp = phi_outputs_;
      ArrOfInt global_indices_post_processed_tmp = global_indices_post_processed_;

      int ii;
      std::vector<int> indices_theta_sorted = arg_sort_array(theta_outputs_);
      for (ii=0; ii<nb_outputs; ii++)
        {
          radius_outputs_(ii) = radius_outputs_tmp(indices_theta_sorted[ii]);
          theta_outputs_(ii) = theta_outputs_tmp(indices_theta_sorted[ii]);
          phi_outputs_(ii) = phi_outputs_tmp(indices_theta_sorted[ii]);
          global_indices_post_processed_(ii) = global_indices_post_processed_tmp(indices_theta_sorted[ii]);
        }

      std::vector<int> indices_phi_sorted = arg_sort_array_phi(theta_scope, theta_outputs_, phi_outputs_);
      radius_outputs_tmp = radius_outputs_;
      theta_outputs_tmp = theta_outputs_;
      phi_outputs_tmp = phi_outputs_;
      global_indices_post_processed_tmp = global_indices_post_processed_;
      for (ii=0; ii<nb_outputs; ii++)
        {
          radius_outputs_(ii) = radius_outputs_tmp(indices_phi_sorted[ii]);
          theta_outputs_(ii) = theta_outputs_tmp(indices_phi_sorted[ii]);
          phi_outputs_(ii) = phi_outputs_tmp(indices_phi_sorted[ii]);
          global_indices_post_processed_(ii) = global_indices_post_processed_tmp(indices_phi_sorted[ii]);
        }

      phi_theta_counter = 0;
      for (ii=0; ii<nb_outputs; ii++)
        {
          const int global_index = global_indices_post_processed_(ii);
          if (global_index >= index_ini_ && global_index < index_end_)
            (*this)[global_index - index_ini_].set_post_processing_theta_phi_scope(phi_theta_counter);
          phi_theta_counter++;
        }
    }
}

void IJK_One_Dimensional_Subproblems::compute_overall_quantities_per_bubbles(const IJK_Field_double& temperature_gradient_ghost,
                                                                             const double& delta_temperature,
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
  // First bubble at index zero
  nb_bubbles_ = (int) std::distance(compo_found.begin(), std::max_element(compo_found.begin(), compo_found.end()));
  nb_bubbles_ += 1;
  nb_bubbles_ = Process::mp_max(nb_bubbles_);
  // nb_bubbles_ = (int) compo_found.size();

  /*
   * Should be the same size on each processor
   */
  interfacial_thermal_flux_per_bubble_.resize(nb_bubbles_);
  interfacial_thermal_flux_per_bubble_gfm_.resize(nb_bubbles_);
  interfacial_thermal_flux_per_bubble_spherical_.resize(nb_bubbles_);
  total_surface_per_bubble_.resize(nb_bubbles_);
  overall_nusselt_number_per_bubble_.resize(nb_bubbles_);
  overall_nusselt_number_per_bubble_gfm_.resize(nb_bubbles_);
  overall_nusselt_number_per_bubble_spherical_.resize(nb_bubbles_);
  overall_nusselt_number_per_bubble_liquid_.resize(nb_bubbles_);
  overall_nusselt_number_per_bubble_gfm_liquid_.resize(nb_bubbles_);
  overall_nusselt_number_per_bubble_spherical_liquid_.resize(nb_bubbles_);
  overall_shear_stress_per_bubble_.resize(nb_bubbles_);
  overall_shear_force_per_bubble_.resize(nb_bubbles_);
  radius_from_surfaces_per_bubble_.resize(nb_bubbles_);
  radius_from_volumes_per_bubble_.resize(nb_bubbles_);
  caracteristic_length_from_surfaces_per_bubble_.resize(nb_bubbles_);
  caracteristic_length_from_volumes_per_bubble_.resize(nb_bubbles_);

  compute_nusselt_numbers_per_bubbles(temperature_gradient_ghost, delta_temperature, lambda);
  compute_shear_per_bubbles();
  compo_found.clear();
}

void IJK_One_Dimensional_Subproblems::compute_overall_bubbles_quantities(IJK_Thermal_Subresolution& ref_thermal_subresolution)
{
  caracteristic_length_ = (ref_thermal_subresolution.single_centred_bubble_radius_ini_) * 2;
  spherical_nusselt_ = ref_thermal_subresolution.nusselt_spherical_diffusion_;
  spherical_nusselt_liquid_ = ref_thermal_subresolution.nusselt_spherical_diffusion_liquid_;
  heat_flux_spherical_ = ref_thermal_subresolution.heat_flux_spherical_;
  mean_liquid_temperature_ = ref_thermal_subresolution.mean_liquid_temperature_;
  compute_overall_quantities_per_bubbles(ref_thermal_subresolution.eulerian_grad_T_interface_ns_,
                                         ref_thermal_subresolution.delta_T_subcooled_overheated_,
                                         ref_thermal_subresolution.uniform_lambda_);
  compute_overall_quantities();
  is_updated_ = true;
}

void IJK_One_Dimensional_Subproblems::compute_nusselt_numbers_per_bubbles(const IJK_Field_double& temperature_gradient_ghost,
                                                                          const double& delta_temperature,
                                                                          const double& lambda)
{
  lambda_ = lambda;
  delta_temperature_ = delta_temperature;
  int local_compo;
  interfacial_thermal_flux_per_bubble_ *= 0.;
  interfacial_thermal_flux_per_bubble_gfm_ *= 0.;
  overall_nusselt_number_per_bubble_ *= 0.;
  overall_nusselt_number_per_bubble_gfm_ *= 0.;
  total_surface_per_bubble_ *= 0.;
  int index_i, index_j, index_k;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      local_compo = (*this)[itr].get_compo();
      (*this)[itr].get_ijk_indices(index_i, index_j, index_k);
      const double local_temperature_gradient_gfm = temperature_gradient_ghost(index_i, index_j, index_k);
      const double ai = (*this)[itr].get_local_surface_area();
      interfacial_thermal_flux_per_bubble_gfm_(local_compo) += (local_temperature_gradient_gfm * ai * lambda_);
      interfacial_thermal_flux_per_bubble_(local_compo) += (*this)[itr].get_interfacial_thermal_flux();
      total_surface_per_bubble_(local_compo) += (*this)[itr].get_local_surface_area();
    }

  mp_sum_for_each_item(interfacial_thermal_flux_per_bubble_gfm_);
  mp_sum_for_each_item(interfacial_thermal_flux_per_bubble_);
  mp_sum_for_each_item(total_surface_per_bubble_);

  // Same on each proc
  for (int i=0; i < nb_bubbles_; i++)
    {
      overall_nusselt_number_per_bubble_(i) = abs((interfacial_thermal_flux_per_bubble_(i) * caracteristic_length_)
                                                  / (total_surface_per_bubble_(i) * delta_temperature_ * lambda_));
      overall_nusselt_number_per_bubble_gfm_(i) = abs((interfacial_thermal_flux_per_bubble_gfm_(i) * caracteristic_length_)
                                                      / (total_surface_per_bubble_(i) * delta_temperature_ * lambda_));
      overall_nusselt_number_per_bubble_liquid_(i) = abs((interfacial_thermal_flux_per_bubble_(i) * caracteristic_length_)
                                                         / (total_surface_per_bubble_(i) * mean_liquid_temperature_ * lambda_));
      overall_nusselt_number_per_bubble_gfm_liquid_(i) = abs((interfacial_thermal_flux_per_bubble_gfm_(i) * caracteristic_length_)
                                                             / (total_surface_per_bubble_(i) * mean_liquid_temperature_ * lambda_));
    }
}

void IJK_One_Dimensional_Subproblems::compute_shear_per_bubbles()
{
  int local_compo;
  overall_shear_force_per_bubble_ *= 0.;
  overall_shear_stress_per_bubble_ *= 0.;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      local_compo = (*this)[itr].get_compo();
      overall_shear_force_per_bubble_(local_compo) += (*this)[itr].get_shear_force();
    }
  mp_sum_for_each_item(overall_shear_force_per_bubble_);

  // Same on each proc
  overall_shear_stress_per_bubble_ = overall_shear_force_per_bubble_;
  for (int i=0; i < nb_bubbles_; i++)
    overall_shear_stress_per_bubble_(i) *= (1 / total_surface_per_bubble_(i));
}

void IJK_One_Dimensional_Subproblems::compute_overall_quantities()
{
  overall_nusselt_number_ = 0.;
  overall_nusselt_number_gfm_ = 0.;
  overall_nusselt_number_spherical_ = 0.;
  overall_nusselt_number_liquid_ = 0.;
  overall_nusselt_number_gfm_liquid_ = 0.;
  overall_nusselt_number_spherical_liquid_ = 0.;
  interfacial_thermal_flux_ = 0.;
  interfacial_thermal_flux_gfm_ = 0.;
  total_surface_ = 0.;
  total_volume_ = 0.;
  overall_shear_force_ = 0.;
  overall_shear_stress_ = 0.;

  for (int i=0; i < nb_bubbles_; i++)
    {
      interfacial_thermal_flux_ += interfacial_thermal_flux_per_bubble_(i);
      interfacial_thermal_flux_gfm_ += interfacial_thermal_flux_per_bubble_gfm_(i);
      total_surface_ += total_surface_per_bubble_(i);
      overall_nusselt_number_ += overall_nusselt_number_per_bubble_(i);
      overall_nusselt_number_gfm_ += overall_nusselt_number_per_bubble_gfm_(i);
      overall_nusselt_number_liquid_ += overall_nusselt_number_per_bubble_liquid_(i);
      overall_nusselt_number_gfm_liquid_ += overall_nusselt_number_per_bubble_gfm_liquid_(i);
      interfacial_thermal_flux_per_bubble_spherical_(i) = heat_flux_spherical_ * total_surface_per_bubble_(i);
      overall_shear_force_ += overall_shear_force_per_bubble_(i);
      total_volume_ += (*bubbles_volume_)(i);
      radius_from_surfaces_per_bubble_(i) = sqrt(total_surface_per_bubble_(i) / (4. * M_PI));
      radius_from_volumes_per_bubble_(i) = pow((*bubbles_volume_)(i) * 3. / (4. * M_PI), (1. / 3.));
      caracteristic_length_from_surfaces_per_bubble_(i) = radius_from_surfaces_per_bubble_(i) * 2;
      caracteristic_length_from_volumes_per_bubble_(i) = radius_from_volumes_per_bubble_(i) * 2;
      overall_nusselt_number_per_bubble_spherical_(i) = (heat_flux_spherical_ * caracteristic_length_from_surfaces_per_bubble_(i))
                                                        / (delta_temperature_ * lambda_);
      overall_nusselt_number_per_bubble_spherical_liquid_(i) = (heat_flux_spherical_ * caracteristic_length_from_surfaces_per_bubble_(i))
                                                               / (mean_liquid_temperature_ * lambda_);
    }
  overall_shear_stress_ = overall_shear_force_ / total_surface_;
  radius_from_surfaces_ = sqrt((total_surface_ / nb_bubbles_) / (4. * M_PI));
  radius_from_volumes_ = pow((total_volume_ / nb_bubbles_) * 3. / (4. * M_PI), (1. / 3.));
  caracteristic_length_from_surfaces_ = radius_from_surfaces_ * 2;
  caracteristic_length_from_volumes_ = radius_from_volumes_ * 2;
  heat_flux_spherical_ *= total_surface_;
  overall_nusselt_number_spherical_ = (heat_flux_spherical_ * caracteristic_length_from_surfaces_)
                                      / (delta_temperature_ * lambda_);
  overall_nusselt_number_spherical_liquid_ = (heat_flux_spherical_ * caracteristic_length_from_surfaces_)
                                             / (mean_liquid_temperature_ * lambda_);
//  interfacial_thermal_flux_ = Process::mp_sum(interfacial_thermal_flux_);
//  interfacial_thermal_flux_gfm_ = Process::mp_sum(interfacial_thermal_flux_gfm_);
//  total_surface_ = Process::mp_sum(total_surface_);
//  overall_nusselt_number_ = Process::mp_sum(overall_nusselt_number_);
//  overall_nusselt_number_gfm_ = Process::mp_sum(overall_nusselt_number_gfm_);
}

void IJK_One_Dimensional_Subproblems::post_process_overall_bubbles_quantities(const int rank, const Nom& overall_bubbles_quantities)
{
  if (Process::je_suis_maitre())
    {
      const int reset = 1;
      const int last_time_index = ref_ijk_ft_->get_tstep();
      const int max_digit = 3;
      const int max_digit_time = 8;
      const int max_rank_digit = rank < 1 ? 1 : (int) (log10(rank) + 1);
      const int nb_digit_tstep = last_time_index < 1 ? 1 : (int) (log10(last_time_index) + 1);

      Nom probe_name = Nom("_thermal_rank_") + Nom(std::string(max_digit - max_rank_digit,'0')) + Nom(rank)
                       + Nom("_thermal_subproblems") + ("_overall_bubbles_quantities_")
                       + Nom(std::string(max_digit_time - nb_digit_tstep, '0')) + Nom(last_time_index) + Nom(".out");
      Nom probe_header = Nom("tstep\ttime\tthermalrank\tbubbleindex"
                             "\tnusseltoverall\tnusseltoverallgfm\tnusseltspherical\tnusseltsphericalth"
                             "\tnusseltoverallliq\tnusseltoverallgfmliq\tnusseltsphericalliq\tnusseltsphericalthliq"
                             "\theatflux\theatfluxgfm\theatfluxspherical"
                             "\ttotalsurface\ttotalvolume"
                             "\tradiussurface\tradiusvolume");
      SFichier fic = Open_file_folder(overall_bubbles_quantities, probe_name, probe_header, reset);
      int max_counter = nb_bubbles_;
      const double last_time = ref_ijk_ft_->get_current_time() - ref_ijk_ft_->get_timestep();
      /*
       * TODO: fill the Array in parallel
       */
      for (int i=0; i<max_counter; i++)
        {
          fic << ref_ijk_ft_->get_tstep() << " " << last_time << " ";
          fic << rank << " ";
          fic << i << " ";
          fic << overall_nusselt_number_per_bubble_(i) << " ";
          fic << overall_nusselt_number_per_bubble_gfm_(i) << " ";
          fic << overall_nusselt_number_per_bubble_spherical_(i) << " ";
          fic << spherical_nusselt_ << " ";
          fic << overall_nusselt_number_per_bubble_liquid_(i) << " ";
          fic << overall_nusselt_number_per_bubble_gfm_liquid_(i) << " ";
          fic << overall_nusselt_number_per_bubble_spherical_liquid_(i) << " ";
          fic << spherical_nusselt_liquid_ << " ";
          fic << interfacial_thermal_flux_per_bubble_(i) << " ";
          fic << interfacial_thermal_flux_per_bubble_gfm_(i) << " ";
          fic << interfacial_thermal_flux_per_bubble_spherical_(i) << " ";
          fic << total_surface_per_bubble_(i) << " ";
          fic << (*bubbles_volume_)(i) << " ";
          fic << radius_from_surfaces_per_bubble_(i) << " ";
          fic << radius_from_volumes_per_bubble_(i) << " ";
          fic << finl;
        }
      /*
       * Should be good for parallel
       */
      if(max_counter > 1)
        {
          fic << ref_ijk_ft_->get_tstep() << " " << last_time << " ";
          fic << rank << " ";
          fic << nb_bubbles_ << " ";
          fic << overall_nusselt_number_ << " ";
          fic << overall_nusselt_number_gfm_ << " ";
          fic << overall_nusselt_number_spherical_ << " ";
          fic << spherical_nusselt_ << " ";
          fic << overall_nusselt_number_liquid_ << " ";
          fic << overall_nusselt_number_gfm_liquid_ << " ";
          fic << overall_nusselt_number_spherical_liquid_ << " ";
          fic << spherical_nusselt_liquid_ << " ";
          fic << interfacial_thermal_flux_ << " ";
          fic << interfacial_thermal_flux_gfm_ << " ";
          fic << heat_flux_spherical_ << " ";
          fic << total_surface_ << " ";
          fic << total_volume_ << " ";
          fic << radius_from_surfaces_ << " ";
          fic << radius_from_volumes_ << " ";
          fic << finl;
        }
      fic.close();
    }
}


