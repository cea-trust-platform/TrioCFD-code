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
  for (int dir=0; dir<3; dir++)
    ijk_indices_to_subproblem_[dir].set_smart_resize(1);
  subproblem_to_ijk_indices_.clear();

  interfacial_thermal_flux_per_bubble_.set_smart_resize(1);
  interfacial_thermal_flux_per_bubble_gfm_.set_smart_resize(1);
  interfacial_thermal_flux_per_bubble_spherical_.set_smart_resize(1);

  total_surface_per_bubble_.set_smart_resize(1);
  bubbles_rising_velocities_ = nullptr;
  bubbles_rising_vectors_per_bubble_ = nullptr;
  bubbles_rising_relative_velocities_.set_smart_resize(1);
  bubbles_rising_relative_velocities_upstream_.set_smart_resize(1);
  bubbles_reynolds_numbers_per_bubble_.set_smart_resize(1);
  bubbles_reynolds_numbers_from_surface_per_bubble_.set_smart_resize(1);
  bubbles_reynolds_numbers_from_surface_per_bubble_.set_smart_resize(1);
  bubbles_reynolds_numbers_per_bubble_upstream_.set_smart_resize(1);
  bubbles_reynolds_numbers_from_surface_per_bubble_upstream_.set_smart_resize(1);
  bubbles_reynolds_numbers_from_surface_per_bubble_upstream_.set_smart_resize(1);
  total_rising_dir_compo_ = 0.;
  total_relative_rising_dir_compo_ = 0.;
  total_relative_rising_dir_compo_upstream_ =0.;
  for (int dir=0; dir<3; dir++)
    {
      relative_rising_dir_compo_[dir].set_smart_resize(1);
      relative_rising_dir_compo_upstream_[dir].set_smart_resize(1);
    }

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

  overall_nusselt_number_per_bubble_error_.set_smart_resize(1);
  overall_nusselt_number_per_bubble_gfm_error_.set_smart_resize(1);
  overall_nusselt_number_per_bubble_liquid_error_.set_smart_resize(1);
  overall_nusselt_number_per_bubble_gfm_liquid_error_.set_smart_resize(1);

  overall_nusselt_number_per_bubble_error_rel_.set_smart_resize(1);
  overall_nusselt_number_per_bubble_gfm_error_rel_.set_smart_resize(1);
  overall_nusselt_number_per_bubble_liquid_error_rel_.set_smart_resize(1);
  overall_nusselt_number_per_bubble_gfm_liquid_error_rel_.set_smart_resize(1);

  bubbles_peclet_numbers_per_bubble_.set_smart_resize(1);
  bubbles_peclet_numbers_from_surface_per_bubble_.set_smart_resize(1);
  bubbles_peclet_numbers_from_volume_per_bubble_.set_smart_resize(1);

  bubbles_peclet_numbers_per_bubble_upstream_.set_smart_resize(1);
  bubbles_peclet_numbers_from_surface_per_bubble_upstream_.set_smart_resize(1);
  bubbles_peclet_numbers_from_volume_per_bubble_upstream_.set_smart_resize(1);

  liquid_velocity_ = nullptr;
  bubbles_volume_ = nullptr;
  prandtl_number_ = nullptr;
  latastep_reprise_ = nullptr;

  points_per_thermal_subproblem_ = nullptr;

  collision_indices_.set_smart_resize(1);
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
                                                                                 const int& remove_append_subproblems,
                                                                                 const int& use_sparse_matrix)
{
  pre_initialise_thermal_subproblems_list_ = pre_initialise_thermal_subproblems_list;
  pre_factor_subproblems_number_ = pre_factor_subproblems_number;
  remove_append_subproblems_ = remove_append_subproblems;
  use_sparse_matrix_ = use_sparse_matrix;
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
  effective_subproblems_counter_ = 0;
  effective_and_disabled_subproblems_counter_ = 0;
  for (int dir=0; dir<3; dir++)
    ijk_indices_to_subproblem_[dir].reset();
  subproblem_to_ijk_indices_previous_ = subproblem_to_ijk_indices_;
  subproblem_to_ijk_indices_.clear();
  one_dimensional_effective_subproblems_.clear();
  one_dimensional_disabled_subproblems_.clear();
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
      // int total_subproblems = subproblems_counter_;
      int total_subproblems = subproblems_counter_;
      if (!(ref_ijk_ft_->get_disable_convection_qdm() && ref_ijk_ft_->get_disable_diffusion_qdm()))
        total_subproblems = Process::mp_sum(total_subproblems);
      max_subproblems_ = (int) (pre_factor_subproblems_number_ * total_subproblems);
      if (effective_subproblems_counter_ < max_subproblems_)//
        {
          const int sub_problem_end_index = effective_subproblems_counter_ - 1;
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

void IJK_One_Dimensional_Subproblems::associate_subproblem_to_ijk_indices(const int& i,
                                                                          const int& j,
                                                                          const int& k)
{
  ijk_indices_to_subproblem_[0].append_array(i);
  ijk_indices_to_subproblem_[1].append_array(j);
  ijk_indices_to_subproblem_[2].append_array(k);
}

void IJK_One_Dimensional_Subproblems::get_ijk_indices_from_subproblems(const int& rank,
                                                                       int& i,
                                                                       int& j,
                                                                       int& k)
{
  i = ijk_indices_to_subproblem_[0][rank];
  j = ijk_indices_to_subproblem_[1][rank];
  k = ijk_indices_to_subproblem_[2][rank];
}

void IJK_One_Dimensional_Subproblems::associate_ijk_indices_to_subproblem(const int& rank,
                                                                          const int& i,
                                                                          const int& j,
                                                                          const int& k)
{
  subproblem_to_ijk_indices_[i][j][k] = rank;
  //	void print_map(std::string_view comment, const std::map<std::string, std::string>& m)
  //	{
  //	    std::cout << comment;
  //	    // Iterate using C++17 facilities
  //	    for (const auto& [key, value] : m)
  //	        std::cout << '[' << key << "] = " << value << "; ";
  //
  //	// C++11 alternative:
  //	//  for (const auto& n : m)
  //	//      std::cout << n.first << " = " << n.second << "; ";
  //	//
  //	// C++98 alternative:
  //	//  for (std::map<std::string, int>::const_iterator it = m.begin(); it != m.end(); ++it)
  //	//      std::cout << it->first << " = " << it->second << "; ";
  //
  //	    std::cout << '\n';
  //	}
  //	std::map<std::string, std::string> inner;
  //	inner.insert(std::make_pair("key2", "value2"));
  //	someStorage.insert(std::make_pair("key", inner));
}

int IJK_One_Dimensional_Subproblems::get_subproblem_index_from_ijk_indices(const int& i,
                                                                           const int& j,
                                                                           const int& k) const
{
  return ((subproblem_to_ijk_indices_.at(i)).at(j)).at(k);
}

void IJK_One_Dimensional_Subproblems::associate_variables_for_post_processing(IJK_Thermal_Subresolution& ref_thermal_subresolution)
{
  if (init_)
    {
      debug_ = ref_thermal_subresolution.debug_;
      reference_gfm_on_probes_ = ref_thermal_subresolution.reference_gfm_on_probes_;
      bubbles_volume_ = ref_thermal_subresolution.bubbles_volume_;
      bubbles_rising_velocities_ = ref_thermal_subresolution.rising_velocities_;
      bubbles_rising_vectors_per_bubble_ = ref_thermal_subresolution.rising_vectors_;
      liquid_velocity_ = ref_thermal_subresolution.liquid_velocity_;
      prandtl_number_ = &ref_thermal_subresolution.prandtl_number_;
      latastep_reprise_= &ref_thermal_subresolution.latastep_reprise_ini_;
      points_per_thermal_subproblem_ = &ref_thermal_subresolution.points_per_thermal_subproblem_;
    }
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

  associate_variables_for_post_processing(ref_thermal_subresolution);
  associate_ijk_indices_to_subproblem(subproblems_counter_, i,j,k);

  ArrOfDouble bubble_rising_vector(3);
  ArrOfDouble normal_vector(3);
  ArrOfDouble facet_barycentre(3);
  ArrOfDouble bubble_barycentre(3);


  if (debug_)
    Cerr << "Mixed cell indices (i,j,k) : (" << i << ";" << j << ";" << k << ")" << finl;
  const int compo_connex = (*ref_thermal_subresolution.eulerian_compo_connex_from_interface_int_ns_)(i, j, k);
  const int compo_connex_ghost = (*ref_thermal_subresolution.eulerian_compo_connex_from_interface_int_ns_)(i,j,k);
  if (debug_)
    {
      Cerr << "compo_connex : " << compo_connex << finl;
      Cerr << "compo_connex_ghost : " << compo_connex_ghost << finl;
      Cerr << "bubbles_barycentre : " << *ref_thermal_subresolution.bubbles_barycentre_ << finl;
    }

  // Need for a Navier-Stokes field (NOT FT)
  const double distance = (*ref_thermal_subresolution.eulerian_distance_ns_)(i, j ,k);
  const double curvature = (*ref_thermal_subresolution.eulerian_curvature_ns_)(i, j ,k);
  const double interfacial_area = (*ref_thermal_subresolution.eulerian_interfacial_area_ns_)(i, j ,k);

  //IJK_Splitting splitting = (*ref_thermal_subresolution.eulerian_compo_connex_from_interface_int_ns_).get_splitting();
  const double bubble_rising_velocity = (*ref_thermal_subresolution.rising_velocities_)(compo_connex);

  for (int dir=0; dir < 3; dir++)
    {
      facet_barycentre(dir) = (*ref_thermal_subresolution.eulerian_facets_barycentre_ns_)[dir](i, j, k);
      normal_vector(dir) = (*ref_thermal_subresolution.eulerian_normal_vectors_ns_)[dir](i, j, k);
      bubble_rising_vector(dir) = (*ref_thermal_subresolution.rising_vectors_)(compo_connex, dir);
      bubble_barycentre(dir) = (*ref_thermal_subresolution.bubbles_barycentre_)(compo_connex_ghost, dir);
    }

  const int total_subproblems = (*this).size();
  if (subproblems_counter_ >= total_subproblems)
    {
      IJK_One_Dimensional_Subproblem subproblem(ref_ijk_ft_);
      (*this).add(subproblem);
      init_ = 1;
    }

  (*this)[subproblems_counter_].associate_sub_problem_to_inputs(ref_thermal_subresolution,
                                                                (*this),
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
  effective_subproblems_counter_++;
}

void IJK_One_Dimensional_Subproblems::store_previous_temperature_indicator_velocities()
{
  temperature_probes_previous_.clear();
  indicator_probes_previous_.clear();
  velocities_probes_previous_.clear();
  normal_vector_compo_probes_previous_.clear();
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      temperature_probes_previous_.push_back((*this)[itr].get_current_temperature_solution());
      indicator_probes_previous_.push_back((*this)[itr].get_current_indicator());
      velocities_probes_previous_.push_back((*this)[itr].get_current_cell_xyz_velocities());
      normal_vector_compo_probes_previous_.push_back((*this)[itr].get_normal_vector_compo());
    }
  share_previous_temperature_indicator_velocities();
}

void IJK_One_Dimensional_Subproblems::share_previous_temperature_indicator_velocities()
{
  const int nk_tot = ref_ijk_ft_->get_splitting_ns().get_grid_geometry().get_nb_elem_tot(2);

  if (ref_ijk_ft_->get_tstep() == 0)
    {
      for (int ij = 0; ij<2; ij++)
        index_ij_subproblems_local_perio_[ij].resize(nk_tot);
      temperature_probes_previous_local_perio_.resize(nk_tot);
      indicator_probes_previous_local_perio_.resize(nk_tot);
      velocities_probes_previous_local_perio_.resize(nk_tot);
      normal_vector_compo_probes_previous_local_perio_.resize(nk_tot);

      for (int ij = 0; ij<2; ij++)
        index_ij_subproblems_global_[ij].resize(nk_tot);
      temperature_probes_previous_global_.resize(nk_tot);
      indicator_probes_previous_global_.resize(nk_tot);
      velocities_probes_previous_global_.resize(nk_tot);
      normal_vector_compo_probes_previous_global_.resize(nk_tot);

      for (int k=0; k<nk_tot; k++)
        {
          for (int ij = 0; ij<2; ij++)
            index_ij_subproblems_local_perio_[ij][k].set_smart_resize(1);
          temperature_probes_previous_local_perio_[k].resize(*points_per_thermal_subproblem_);
          for (int point=0; point<(*points_per_thermal_subproblem_); point++)
            temperature_probes_previous_local_perio_[k][point].set_smart_resize(1);
          indicator_probes_previous_local_perio_[k].set_smart_resize(1);
          velocities_probes_previous_local_perio_[k].resize(3);
          normal_vector_compo_probes_previous_local_perio_[k].resize(3);
          for (int compo=0; compo<3; compo++)
            {
              velocities_probes_previous_local_perio_[k][compo].set_smart_resize(1);
              normal_vector_compo_probes_previous_local_perio_[k][compo].set_smart_resize(1);
            }
          for (int ij = 0; ij<2; ij++)
            index_ij_subproblems_global_[ij][k].set_smart_resize(1);
          temperature_probes_previous_global_[k].resize(*points_per_thermal_subproblem_);
          for (int point=0; point<(*points_per_thermal_subproblem_); point++)
            temperature_probes_previous_global_[k][point].set_smart_resize(1);
          indicator_probes_previous_global_[k].set_smart_resize(1);
          velocities_probes_previous_global_[k].resize(3);
          normal_vector_compo_probes_previous_global_[k].resize(3);
          for (int compo=0; compo<3; compo++)
            {
              velocities_probes_previous_global_[k][compo].set_smart_resize(1);
              normal_vector_compo_probes_previous_global_[k][compo].set_smart_resize(1);
            }
        }
    }
  else
    {
      for (int k=0; k<nk_tot; k++)
        {
          for (int ij = 0; ij<2; ij++)
            index_ij_subproblems_local_perio_[ij][k].reset();
          for (int point=0; point<(*points_per_thermal_subproblem_); point++)
            temperature_probes_previous_local_perio_[k][point].reset();
          indicator_probes_previous_local_perio_[k].reset();
          for (int compo=0; compo<3; compo++)
            {
              velocities_probes_previous_local_perio_[k][compo].reset();
              normal_vector_compo_probes_previous_local_perio_[k][compo].reset();
            }

          for (int ij = 0; ij<2; ij++)
            index_ij_subproblems_global_[ij][k].reset();
          for (int point=0; point<(*points_per_thermal_subproblem_); point++)
            temperature_probes_previous_global_[k][point].reset();
          indicator_probes_previous_global_[k].reset();
          for (int compo=0; compo<3; compo++)
            {
              velocities_probes_previous_global_[k][compo].reset();
              normal_vector_compo_probes_previous_global_[k][compo].reset();
            }
        }
    }
  retrieve_boundary_previous_values();
  share_boundary_previous_values();
  complete_boundary_previous_values();
}

void IJK_One_Dimensional_Subproblems::retrieve_boundary_previous_values()
{
  const int ni = ref_ijk_ft_->itfce().I().ni();
  const int nj = ref_ijk_ft_->itfce().I().nj();
  const int nk = ref_ijk_ft_->itfce().I().nk();

  const int offset_i = ref_ijk_ft_->get_splitting_ns().get_offset_local(0);
  const int offset_j = ref_ijk_ft_->get_splitting_ns().get_offset_local(1);
  const int offset_k = ref_ijk_ft_->get_splitting_ns().get_offset_local(2);

  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      int i, j, k;
      (*this)[itr].get_ijk_indices(i, j, k);
      if ((i==ni-1 || i==0) || (j==nj-1 || j==0) || (k==nk-1 || k==0))
        {
          const int i_global = i + offset_i;
          const int j_global = j + offset_j;
          const int k_global = k + offset_k;
          index_ij_subproblems_local_perio_[0][k_global].append_array(i_global);
          index_ij_subproblems_local_perio_[1][k_global].append_array(j_global);
          DoubleVect& temperature_vect = temperature_probes_previous_[itr];
          Vecteur3& velocities_xyz = velocities_probes_previous_[itr];
          Vecteur3& normal_vector_compo = normal_vector_compo_probes_previous_[itr];
          for (int point=0; point<(*points_per_thermal_subproblem_); point++)
            temperature_probes_previous_local_perio_[k_global][point].append_array(temperature_vect[point]);
          indicator_probes_previous_local_perio_[k_global].append_array(indicator_probes_previous_[itr]);
          for (int compo=0; compo<3; compo++)
            {
              velocities_probes_previous_local_perio_[k_global][compo].append_array(velocities_xyz[compo]);
              normal_vector_compo_probes_previous_local_perio_[k_global][compo].append_array(normal_vector_compo[compo]);
            }
        }
    }
}

void IJK_One_Dimensional_Subproblems::share_boundary_previous_values()
{
  const int nb_procs = Process::nproc();
  const int proc_num = Process::me();
  if (nb_procs)
    {
      const int size_k_layers = (int) index_ij_subproblems_local_perio_[0].size();
      for (int k=0; k<size_k_layers; k++)
        {
          const int size_array = index_ij_subproblems_local_perio_[0][k].size_array();
          int size_array_global = size_array;
          size_array_global = mp_sum(size_array_global);
          ArrOfInt overall_numerotation(nb_procs);
          ArrOfInt start_indices(nb_procs);
          overall_numerotation(proc_num) = size_array;
          mp_sum_for_each_item(overall_numerotation);
          int l;
          for (l=1; l<overall_numerotation.size_array(); l++)
            start_indices(l) = start_indices(l-1) + overall_numerotation(l-1);

          const ArrOfInt& local_indices_i_tmp = index_ij_subproblems_local_perio_[0][k];
          const ArrOfInt& local_indices_j_tmp = index_ij_subproblems_local_perio_[1][k];
          const ArrOfDouble& local_indicator = indicator_probes_previous_local_perio_[k];


          ArrOfInt& global_indices_i_tmp = index_ij_subproblems_global_[0][k];
          ArrOfInt& global_indices_j_tmp = index_ij_subproblems_global_[1][k];
          ArrOfDouble& global_indicator = indicator_probes_previous_global_[k];

          global_indices_i_tmp.resize(size_array_global);
          global_indices_j_tmp.resize(size_array_global);
          global_indicator.resize(size_array_global);

          global_indices_i_tmp *= 0;
          global_indices_j_tmp *= 0;
          global_indicator *= 0.;

          for (l=0; l<local_indices_i_tmp.size_array(); l++)
            {
              global_indices_i_tmp(start_indices(proc_num) + l) = local_indices_i_tmp(l);
              global_indices_j_tmp(start_indices(proc_num) + l) = local_indices_j_tmp(l);
              global_indicator(start_indices(proc_num) + l) = local_indicator(l);
            }
          mp_sum_for_each_item(global_indices_i_tmp);
          mp_sum_for_each_item(global_indices_j_tmp);
          mp_sum_for_each_item(global_indicator);

          for (int point=0; point<(*points_per_thermal_subproblem_); point++)
            {
              const ArrOfDouble& local_temperature = temperature_probes_previous_local_perio_[k][point];
              ArrOfDouble& global_temperature = temperature_probes_previous_global_[k][point];
              global_temperature.resize(size_array_global);
              global_temperature *= 0.;
              for (l=0; l<local_temperature.size_array(); l++)
                global_temperature(start_indices(proc_num) + l) = local_temperature(l);
              mp_sum_for_each_item(global_temperature);
            }

          for (int compo=0; compo<3; compo++)
            {
              const ArrOfDouble& local_velocity = velocities_probes_previous_local_perio_[k][compo];
              const ArrOfDouble& local_normal_compo = normal_vector_compo_probes_previous_local_perio_[k][compo];
              ArrOfDouble& global_velocity = velocities_probes_previous_global_[k][compo];
              ArrOfDouble& global_normal_compo = normal_vector_compo_probes_previous_global_[k][compo];
              global_velocity.resize(size_array_global);
              global_normal_compo.resize(size_array_global);
              global_velocity *= 0.;
              global_normal_compo *= 0.;
              for (l=0; l<local_velocity.size_array(); l++)
                {
                  global_velocity(start_indices(proc_num) + l) = local_velocity(l);
                  global_normal_compo(start_indices(proc_num) + l) = local_normal_compo(l);
                }
              mp_sum_for_each_item(global_velocity);
              mp_sum_for_each_item(global_normal_compo);
            }
        }
    }
}

void IJK_One_Dimensional_Subproblems::complete_boundary_previous_values()
{
  const int ni = ref_ijk_ft_->itfce().I().ni();
  const int nj = ref_ijk_ft_->itfce().I().nj();
  const int nk = ref_ijk_ft_->itfce().I().nk();

  const int offset_i = ref_ijk_ft_->get_splitting_ns().get_offset_local(0);
  const int offset_j = ref_ijk_ft_->get_splitting_ns().get_offset_local(1);
  const int offset_k = ref_ijk_ft_->get_splitting_ns().get_offset_local(2);

  const int ni_tot = ref_ijk_ft_->get_splitting_ns().get_grid_geometry().get_nb_elem_tot(0);
  const int nj_tot = ref_ijk_ft_->get_splitting_ns().get_grid_geometry().get_nb_elem_tot(1);
  const int nk_tot = ref_ijk_ft_->get_splitting_ns().get_grid_geometry().get_nb_elem_tot(2);

  const int nb_procs = Process::nproc();
  int counter_prev = subproblems_counter_;
  if (nb_procs)
    {
      const int size_k_layers = (int) index_ij_subproblems_global_[0].size();
      for (int k_global=0; k_global<size_k_layers; k_global++)
        {
          const int size_val = (int) index_ij_subproblems_global_[0][k_global].size_array();
          for (int ival=0; ival<size_val; ival++)
            {
              const int i_global = index_ij_subproblems_global_[0][k_global][ival];
              const int j_global = index_ij_subproblems_global_[1][k_global][ival];
              /*
               * Global perio
               */
              int i = i_global - offset_i;
              int j = j_global - offset_j;
              int k = k_global - offset_k;
              if (i_global == ni_tot-1 && offset_i==0)
                i = -1;
              if (i_global == 0 && (ni + offset_i==ni_tot))
                i = ni;
              if (j_global == nj_tot-1 && offset_j==0)
                j = -1;
              if (j_global == 0 && (nj + offset_j==nj_tot))
                j = nj;
              if (k_global == nk_tot-1 && offset_k==0)
                k = -1;
              if (k_global == 0 && (nk + offset_k==nk_tot))
                k = nk;
              /*
               * Neighbours
               */
              const int i_dir = (i==-1 || i==ni) && (0 <= j && j < nj) && (0 <= k && k < nk);
              const int j_dir = (0 <= i && i < ni) && (j==-1 || j==nj) && (0 <= k && k < nk);
              const int k_dir = (0 <= i && i < ni) && (0 <= j && j < nj) && (k==-1 || k==nk);
              if (i_dir || j_dir || k_dir)
                if (!is_in_map_index_ijk(subproblem_to_ijk_indices_, i, j, k))
                  {
                    subproblem_to_ijk_indices_[i][j][k] = counter_prev;

                    DoubleVect temperature_prev;
                    temperature_prev.resize((*points_per_thermal_subproblem_));
                    for (int point=0; point<(*points_per_thermal_subproblem_); point++)
                      temperature_prev[point]	= temperature_probes_previous_global_[k_global][point][ival];
                    temperature_probes_previous_.push_back(temperature_prev);

                    indicator_probes_previous_.push_back(indicator_probes_previous_global_[k_global][ival]);

                    Vecteur3 velocities_compo;
                    Vecteur3 normal_vector_compo;
                    for (int compo=0; compo<3; compo++)
                      {
                        velocities_compo[compo] = velocities_probes_previous_global_[k_global][compo][ival];
                        normal_vector_compo[compo] = normal_vector_compo_probes_previous_global_[k_global][compo][ival];
                      }
                    velocities_probes_previous_.push_back(velocities_compo);
                    normal_vector_compo_probes_previous_.push_back(normal_vector_compo);

                    counter_prev++;
                  }
            }
        }
    }
  else
    {
      const int size_k_layers = (int) index_ij_subproblems_local_perio_[0].size();
      for (int k_global=0; k_global<size_k_layers; k_global++)
        {
          const int size_val = (int) index_ij_subproblems_global_[0][k_global].size_array();
          for (int ival=0; ival<size_val; ival++)
            {
              /*
               * Global perio
               */
              int i = index_ij_subproblems_local_perio_[0][k_global][ival];
              int j = index_ij_subproblems_local_perio_[1][k_global][ival];
              int k = k_global - offset_k;
              if (i == ni-1)
                i = 0;
              if (i == 0)
                i = ni-1;
              if (j == nj-1)
                j = 0;
              if (j == 0)
                j = nj-1;
              if (k == nk-1)
                k = 0;
              if (k == 0)
                k = nk-1;
              const int i_dir = (i == ni-1 || i == 0);
              const int j_dir = (j == nj-1 || j == 0);
              const int k_dir = (k == nk-1 || k == 0);
              if (i_dir || j_dir || k_dir)
                if (!is_in_map_index_ijk(subproblem_to_ijk_indices_, i, j, k))
                  {
                    subproblem_to_ijk_indices_[i][j][k] = counter_prev;

                    DoubleVect temperature_prev;
                    temperature_prev.resize((*points_per_thermal_subproblem_));
                    for (int point=0; point<(*points_per_thermal_subproblem_); point++)
                      temperature_prev[point]	= temperature_probes_previous_local_perio_[k_global][point][ival];
                    temperature_probes_previous_.push_back(temperature_prev);

                    indicator_probes_previous_.push_back(indicator_probes_previous_local_perio_[k_global][ival]);

                    Vecteur3 velocities_compo;
                    Vecteur3 normal_vector_compo;
                    for (int compo=0; compo<3; compo++)
                      {
                        velocities_compo[compo] = velocities_probes_previous_local_perio_[k_global][compo][ival];
                        normal_vector_compo[compo] = normal_vector_compo_probes_previous_local_perio_[k_global][compo][ival];
                      }
                    velocities_probes_previous_.push_back(velocities_compo);
                    normal_vector_compo_probes_previous_.push_back(normal_vector_compo);

                    counter_prev++;
                  }
            }
        }
    }
}

int IJK_One_Dimensional_Subproblems::is_in_map_index_ijk(const std::map<int, std::map<int, std::map<int, int>>>& subproblem_to_ijk_indices,
                                                         const int& index_i,
                                                         const int& index_j,
                                                         const int& index_k)
{
  const int count_index_i = (int) subproblem_to_ijk_indices.count(index_i);
  int count_index_j = 0;
  int count_index_k = 0;
  if (count_index_i)
    {
      count_index_j = (int) subproblem_to_ijk_indices.at(index_i).count(index_j);
      if (count_index_j)
        count_index_k = (int) subproblem_to_ijk_indices.at(index_i).at(index_j).count(index_k);
    }
  return (count_index_i && count_index_j && count_index_k);
}

void IJK_One_Dimensional_Subproblems::set_effective_subproblems(const int& enable_probe_collision_detection)
{
  if (!enable_probe_collision_detection)
    {
      effective_subproblems_counter_ = subproblems_counter_;
      effective_and_disabled_subproblems_counter_ = subproblems_counter_;
      for (int itr=0; itr < effective_subproblems_counter_; itr++)
        one_dimensional_effective_subproblems_.push_back(&((*this)[itr]));
    }
}

void IJK_One_Dimensional_Subproblems::interpolate_indicator_on_probes()
{
  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].interpolate_indicator_on_probes();
}

void IJK_One_Dimensional_Subproblems::clear_sort_problems_colliding_bubbles()
{
  effective_subproblems_counter_ = subproblems_counter_;
  collision_indices_.reset();
  //  if ((*bubbles_rising_velocities_).size_array() == 1)
  //    return;
  int disable_probe_collision = 0;
  int counter_enabled = 0;
  for (int itr=0; itr < subproblems_counter_; itr++)
    {
      disable_probe_collision = (*this)[itr].get_disable_probe_collision();
      if (disable_probe_collision)
        {
          collision_indices_.append_array(itr);
          (*this)[itr].set_reference_gfm_on_probes(1);
        }
      else
        {
          (*this)[itr].set_subproblem_index(counter_enabled);
          one_dimensional_effective_subproblems_.push_back(&((*this)[itr]));
          counter_enabled++;
        }
    }
  effective_subproblems_counter_ -= (collision_indices_.size_array());
  assert(effective_subproblems_counter_==counter_enabled);
  int counter_disabled = effective_subproblems_counter_;
  for (int itr=0; itr < collision_indices_.size_array(); itr++)
    {
      (*this)[collision_indices_[itr]].set_subproblem_index(counter_disabled);
      one_dimensional_disabled_subproblems_.push_back(&((*this)[collision_indices_[itr]]));
      counter_disabled++;
    }
  effective_and_disabled_subproblems_counter_ = counter_disabled;
  assert(effective_and_disabled_subproblems_counter_ == subproblems_counter_);
  disabled_subproblems_counter_ = effective_and_disabled_subproblems_counter_ - effective_subproblems_counter_;
}

void IJK_One_Dimensional_Subproblems::interpolate_project_velocities_on_probes()
{
  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].interpolate_project_velocities_on_probes();
}

void IJK_One_Dimensional_Subproblems::reajust_probes_length()
{
  for (int itr=0; itr < effective_subproblems_counter_; itr++)
    (*this)[itr].reajust_probe_length();
}

void IJK_One_Dimensional_Subproblems::reajust_probes_length(const int probe_length_condition)
{
  for (int itr=0; itr < effective_subproblems_counter_; itr++)
    (*this)[itr].compute_modified_probe_length_condition(probe_length_condition);
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
  for (int itr=0; itr < effective_subproblems_counter_; itr++)
    probe_variations_enabled = (probe_variations_enabled || one_dimensional_effective_subproblems_[itr]->get_probe_variations_enabled());
  return probe_variations_enabled;
}

int IJK_One_Dimensional_Subproblems::get_probe_variations_enabled_non_priority()
{
  int probe_variations_enabled = 1;
  for (int itr=0; itr < effective_subproblems_counter_; itr++)
    probe_variations_enabled = (probe_variations_enabled && one_dimensional_effective_subproblems_[itr]->get_probe_variations_enabled());
  return probe_variations_enabled;
}

void IJK_One_Dimensional_Subproblems::compute_modified_probe_length(const int& probe_variations_enabled)
{
  for (int itr=0; itr < effective_subproblems_counter_; itr++)
    (*this)[itr].compute_modified_probe_length(probe_variations_enabled);
}

void IJK_One_Dimensional_Subproblems::compute_radial_convection_diffusion_operators()
{
  for (int itr=0; itr < effective_subproblems_counter_; itr++)
    one_dimensional_effective_subproblems_[itr]->compute_radial_convection_diffusion_operators();
}

void IJK_One_Dimensional_Subproblems::compute_source_terms_impose_boundary_conditions(const int& boundary_condition_interface,
                                                                                      const double& interfacial_boundary_condition_value,
                                                                                      const int& impose_boundary_condition_interface_from_simulation,
                                                                                      const int& boundary_condition_end,
                                                                                      const double& end_boundary_condition_value,
                                                                                      const int& impose_user_boundary_condition_end_value)
{
  for (int itr=0; itr < effective_subproblems_counter_; itr++)
    one_dimensional_effective_subproblems_[itr]->compute_source_terms_impose_boundary_conditions(boundary_condition_interface,
                                                                                                 interfacial_boundary_condition_value,
                                                                                                 impose_boundary_condition_interface_from_simulation,
                                                                                                 boundary_condition_end,
                                                                                                 end_boundary_condition_value,
                                                                                                 impose_user_boundary_condition_end_value);
}

void IJK_One_Dimensional_Subproblems::approximate_temperature_increment_material_derivative()
{
  for (int itr=0; itr < effective_subproblems_counter_; itr++)
    one_dimensional_effective_subproblems_[itr]->approximate_temperature_increment_material_derivative();
  for (int itr=0; itr < disabled_subproblems_counter_; itr++)
    one_dimensional_disabled_subproblems_[itr]->approximate_temperature_increment_material_derivative();
}

void IJK_One_Dimensional_Subproblems::retrieve_radial_quantities()
{
  for (int itr=0; itr < effective_subproblems_counter_; itr++)
    one_dimensional_effective_subproblems_[itr]->retrieve_radial_quantities();
  for (int itr=0; itr < disabled_subproblems_counter_; itr++)
    one_dimensional_disabled_subproblems_[itr]->retrieve_radial_quantities();
}

void IJK_One_Dimensional_Subproblems::retrieve_temperature_solutions()
{
  for (int itr=0; itr < effective_subproblems_counter_; itr++)
    one_dimensional_effective_subproblems_[itr]->retrieve_temperature_solution();
}

void IJK_One_Dimensional_Subproblems::compute_local_temperature_gradient_solutions()
{
  for (int itr=0; itr < effective_subproblems_counter_; itr++)
    one_dimensional_effective_subproblems_[itr]->compute_local_temperature_gradient_solution();
}

void IJK_One_Dimensional_Subproblems::compute_local_velocity_gradient()
{
  for (int itr=0; itr < effective_subproblems_counter_; itr++)
    one_dimensional_effective_subproblems_[itr]->compute_local_velocity_gradient();
}

void IJK_One_Dimensional_Subproblems::get_subproblem_ijk_indices(int& i, int& j, int& k, int& subproblem_index) const
{
  one_dimensional_effective_subproblems_[subproblem_index]->get_ijk_indices(i,j,k);
}

const int& IJK_One_Dimensional_Subproblems::get_dxyz_increment_bool(const int& subproblem_index) const
{
  return one_dimensional_effective_subproblems_[subproblem_index]->get_dxyz_increment_bool();
}

const int& IJK_One_Dimensional_Subproblems::get_dxyz_over_two_increment_bool(const int& subproblem_index) const
{
  return one_dimensional_effective_subproblems_[subproblem_index]->get_dxyz_over_two_increment_bool();
}

const FixedVector<int,3>& IJK_One_Dimensional_Subproblems::get_pure_neighbours_corrected_sign(const int& subproblem_index) const
{
  return one_dimensional_effective_subproblems_[subproblem_index]->get_pure_neighbours_corrected_sign();
}

const std::vector<std::vector<std::vector<bool>>>& IJK_One_Dimensional_Subproblems::get_pure_neighbours_to_correct(const int& subproblem_index) const
{
  return one_dimensional_effective_subproblems_[subproblem_index]->get_pure_neighbours_to_correct();
}

const std::vector<std::vector<std::vector<double>>>& IJK_One_Dimensional_Subproblems::get_pure_neighbours_corrected_distance(const int& subproblem_index) const
{
  return one_dimensional_effective_subproblems_[subproblem_index]->get_pure_neighbours_corrected_distance();
}

const std::vector<std::vector<std::vector<double>>>& IJK_One_Dimensional_Subproblems::get_pure_neighbours_corrected_colinearity(const int& subproblem_index) const
{
  return one_dimensional_effective_subproblems_[subproblem_index]->get_pure_neighbours_corrected_colinearity();
}

const std::vector<std::vector<std::vector<std::vector<bool>>>> IJK_One_Dimensional_Subproblems::get_pure_neighbours_last_faces_to_correct(const int& subproblem_index) const
{
  return one_dimensional_effective_subproblems_[subproblem_index]->get_pure_neighbours_last_faces_to_correct();
}
const std::vector<std::vector<std::vector<std::vector<double>>>> IJK_One_Dimensional_Subproblems::get_pure_neighbours_last_faces_corrected_distance(const int& subproblem_index) const
{
  return one_dimensional_effective_subproblems_[subproblem_index]->get_pure_neighbours_last_faces_corrected_distance();
}
const std::vector<std::vector<std::vector<std::vector<double>>>> IJK_One_Dimensional_Subproblems::get_pure_neighbours_last_faces_corrected_colinearity(const int& subproblem_index) const
{
  return one_dimensional_effective_subproblems_[subproblem_index]->get_pure_neighbours_last_faces_corrected_colinearity();
}

double IJK_One_Dimensional_Subproblems::get_interfacial_gradient_corrected(int i)
{
  return one_dimensional_effective_subproblems_[i]->get_interfacial_gradient_corrected();
}

double IJK_One_Dimensional_Subproblems::get_temperature_profile_at_point(const int& i, const double& dist) const
{
  return one_dimensional_effective_subproblems_[i]->get_temperature_profile_at_point(dist);
}

const Vecteur3& IJK_One_Dimensional_Subproblems::get_bary_facet(const int& i) const
{
  return one_dimensional_effective_subproblems_[i]->get_bary_facet();
}

const double& IJK_One_Dimensional_Subproblems::get_dist_cell_interface(const int& i) const
{
  return one_dimensional_effective_subproblems_[i]->get_dist_cell();
}

const FixedVector<double,6>& IJK_One_Dimensional_Subproblems::get_dist_faces_interface(const int& i) const
{
  return one_dimensional_effective_subproblems_[i]->get_dist_faces();
}

double IJK_One_Dimensional_Subproblems::get_temperature_times_velocity_profile_at_point(const int& i, const double& dist, const int& dir) const
{
  return one_dimensional_effective_subproblems_[i]->get_temperature_times_velocity_profile_at_point(dist, dir);
}

DoubleVect IJK_One_Dimensional_Subproblems::get_temperature_profile_discrete_integral_at_point(const int& i, const double& dist, const int& level, const int& dir) const
{
  return one_dimensional_effective_subproblems_[i]->get_temperature_profile_discrete_integral_at_point(dist, level, dir);
}

DoubleVect IJK_One_Dimensional_Subproblems::get_temperature_times_velocity_profile_discrete_integral_at_point(const int& i, const double& dist, const int& level, const int& dir) const
{
  return one_dimensional_effective_subproblems_[i]->get_temperature_times_velocity_profile_discrete_integral_at_point(dist, level, dir);
}

double IJK_One_Dimensional_Subproblems::get_temperature_gradient_profile_at_point(const int& i, const double& dist, const int& dir) const
{
  return one_dimensional_effective_subproblems_[i]->get_temperature_gradient_profile_at_point(dist, dir);
}

double IJK_One_Dimensional_Subproblems::get_temperature_gradient_times_conductivity_profile_at_point(const int& i, const double& dist, const int& dir) const
{
  return one_dimensional_effective_subproblems_[i]->get_temperature_gradient_times_conductivity_profile_at_point(dist, dir);
}

DoubleVect IJK_One_Dimensional_Subproblems::get_temperature_gradient_profile_discrete_integral_at_point(const int& i, const double& dist, const int& level, const int& dir) const
{
  return one_dimensional_effective_subproblems_[i]->get_temperature_gradient_profile_discrete_integral_at_point(dist, level, dir);
}

DoubleVect IJK_One_Dimensional_Subproblems::get_temperature_gradient_times_conductivity_profile_discrete_integral_at_point(const int& i, const double& dist, const int& level, const int& dir) const
{
  return one_dimensional_effective_subproblems_[i]->get_temperature_gradient_times_conductivity_profile_discrete_integral_at_point(dist, level, dir);
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

  std::vector<std::string> key_results_int = {"tstep", "thermal_rank", "post_pro_index", "global_subproblem", "local_subproblem"};
  std::vector<std::string> key_results_double = {"time",
                                                 "nx", "ny", "nz",
                                                 "t1x", "t1y", "t1z", "t2x", "t2y", "t2z",
                                                 "s1x", "s1y", "s1z", "s2x", "s2y", "s2z",
                                                 "r_sph", "theta_sph", "phi_sph",
                                                 "temperature_interp","temperature_sol", "temperature_prev",
                                                 "temperature_gradient", "temperature_gradient_sol",
                                                 "temperature_double_deriv", "temperature_double_deriv_sol",
                                                 "temperature_gradient_tangential","temperature_gradient_tangential2",
                                                 "temperature_gradient_tangential_rise","temperature_gradient_azymuthal",
                                                 "temperature_diffusion_hessian_cartesian_trace",
                                                 "temperature_diffusion_hessian_trace",
                                                 "radial_temperature_diffusion",
                                                 "radial_temperature_diffusion_sol",
                                                 "tangential_temperature_diffusion",
                                                 "radial_scale_factor_inerp", "radial_scale_factor_sol",
                                                 "radial_convection_interp", "radial_convection_sol",
                                                 "tangential_convection_first", "tangential_convection_second",
                                                 "surface","thermal_flux","lambda_liq","alpha_liq","prandtl_liq",
                                                 "nusselt_number","nusselt_number_liquid_temperature",
                                                 "nusselt_number_integrand","nusselt_number_liquid_temperature_integrand",
                                                 "shear","force",
                                                 "pressure",
                                                 "u_x","u_y","u_z",
                                                 "u_r","u_r_corr","u_r_static","u_r_advected",
                                                 "u_theta","u_theta_corr","u_theta_static","u_theta_advected",
                                                 "u_theta2","u_theta2_corr","u_theta2_static","u_theta2_advected",
                                                 "u_theta_rise","u_theta_rise_corr","u_theta_rise_static","u_theta_rise_advected",
                                                 "u_phi","u_phi_corr","u_phi_static","u_phi_advected",
                                                 "du_r_dr","du_theta_dr","du_theta2_dr","du_theta_rise_dr","du_phi_dr",
                                                 "total_surface", "total_volume", "radius_from_surface", "radius_from_volume",
                                                 "delta_temperature", "mean_liquid_temperature",
                                                 "rising_dir_x", "rising_dir_y", "rising_dir_z"
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
  const int last_time_index = ref_ijk_ft_->get_tstep() + (*latastep_reprise_);
  const int max_digit = 3;
  const int max_digit_time = 8;
  const int max_rank_digit = rank < 1 ? 1 : (int) (log10(rank) + 1);
  const int nb_digit_tstep = last_time_index < 1 ? 1 : (int) (log10(last_time_index) + 1);

  if (debug_)
    Cerr << "Post-process interfacial quantities" << finl;
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

  if (debug_)
    Cerr << "Post-process local radial quantities" << finl;
  for (int itr=0; itr < subproblems_counter_; itr++)
    (*this)[itr].thermal_subresolution_outputs_parallel(rank, local_quantities_thermal_probes_time_index_folder);

  if (debug_)
    Cerr << "Post-process overall bubble quantities" << finl;
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
  const int reset = 1;
  const int last_time_index = ref_ijk_ft_->get_tstep() + (*latastep_reprise_);
  Nom probe_header = Nom("tstep\tthermal_rank\tpost_pro_index\tglobal_subproblem\tlocal_subproblem\ttime"
                         "\tnx\tny\tnz\tt1x\tt1y\tt2z\tt2x\tt2y\tt2z\ts1x\ts1y\ts1z\ts2x\ts2y\ts2z"
                         "\tr_sph\ttheta_sph\tphi_sph"
                         "\ttemperature_interp\ttemperature_sol\ttemperature_prev"
                         "\ttemperature_gradient\ttemperature_gradient_sol"
                         "\ttemperature_double_deriv\ttemperature_double_deriv_sol"
                         "\ttemperature_gradient_tangential\ttemperature_gradient_tangential2"
                         "\ttemperature_gradient_tangential_rise\ttemperature_gradient_azymuthal"
                         "\ttemperature_diffusion_hessian_cartesian_trace"
                         "\ttemperature_diffusion_hessian_trace"
                         "\tradial_temperature_diffusion"
                         "\tradial_temperature_diffusion_sol"
                         "\ttangential_temperature_diffusion"
                         "\tradial_scale_factor_interp\tradial_scale_factor_sol"
                         "\tradial_convection_interp\tradial_convection_sol"
                         "\ttangential_convection_first\ttangential_convection_second"
                         "\tsurface\tthermal_flux\tlambda_liq\talpha_liq\tprandtl_liq"
                         "\tnusselt_number\tnusselt_number_liquid_temperature"
                         "\tnusselt_number_integrand\tnusselt_number_liquid_temperature_integrand"
                         "\tshear\tforce"
                         "\tpressure"
                         "\tu_x\tu_y\tu_z"
                         "\tu_r\tu_r_corr\tu_r_static\tu_r_advected"
                         "\tu_theta\tu_theta_corr\tu_theta_static\tu_theta_advected"
                         "\tu_theta2\tu_theta2_corr\tu_theta2_static\tu_theta2_advected"
                         "\tu_theta_rise\tu_theta_rise_corr\tu_theta_rise_static\tu_theta_rise_advected"
                         "\tu_phi\tu_phi_corr\tu_phi_static\tu_phi_advected"
                         "\tdu_r_dr\tdu_theta_dr\tdu_theta2_dr\tdu_theta_rise_dr\tdu_phi_dr"
                         "\ttotal_surface\ttotal_volume\tradius_from_surface\tradius_from_volume"
                         "\tdelta_temperature\tmean_liquid_temperature"
                         "\trising_dir_x\trising_dir_y\trising_dir_z");

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

            radius_outputs_(phi_theta_counter) = r_sph(theta_phi_scope_index);
            theta_outputs_(phi_theta_counter) = theta_sph(theta_phi_scope_index);
            phi_outputs_(phi_theta_counter) = phi_sph(theta_phi_scope_index);
            global_indices_post_processed_(phi_theta_counter) = theta_phi_scope_index;
            phi_theta_counter ++;
          }

      if (debug_)
        {
          Cerr << radius_outputs_ << finl;
          Cerr << theta_outputs_ << finl;
          Cerr << phi_outputs_ << finl;
          Cerr << global_indices_post_processed_ << finl;
        }

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

  overall_nusselt_number_per_bubble_error_.resize(nb_bubbles_);
  overall_nusselt_number_per_bubble_gfm_error_.resize(nb_bubbles_);
  overall_nusselt_number_per_bubble_liquid_error_.resize(nb_bubbles_);
  overall_nusselt_number_per_bubble_gfm_liquid_error_.resize(nb_bubbles_);

  overall_nusselt_number_per_bubble_error_rel_.resize(nb_bubbles_);
  overall_nusselt_number_per_bubble_gfm_error_rel_.resize(nb_bubbles_);
  overall_nusselt_number_per_bubble_liquid_error_rel_.resize(nb_bubbles_);
  overall_nusselt_number_per_bubble_gfm_liquid_error_rel_.resize(nb_bubbles_);

  bubbles_rising_relative_velocities_.resize(nb_bubbles_);
  bubbles_rising_relative_velocities_upstream_.resize(nb_bubbles_);

  bubbles_reynolds_numbers_per_bubble_.resize(nb_bubbles_);
  bubbles_reynolds_numbers_from_surface_per_bubble_.resize(nb_bubbles_);
  bubbles_reynolds_numbers_from_volume_per_bubble_.resize(nb_bubbles_);

  bubbles_reynolds_numbers_per_bubble_upstream_.resize(nb_bubbles_);
  bubbles_reynolds_numbers_from_surface_per_bubble_upstream_.resize(nb_bubbles_);
  bubbles_reynolds_numbers_from_volume_per_bubble_upstream_.resize(nb_bubbles_);

  bubbles_peclet_numbers_per_bubble_.resize(nb_bubbles_);
  bubbles_peclet_numbers_from_surface_per_bubble_.resize(nb_bubbles_);
  bubbles_peclet_numbers_from_volume_per_bubble_.resize(nb_bubbles_);

  bubbles_peclet_numbers_per_bubble_upstream_.resize(nb_bubbles_);
  bubbles_peclet_numbers_from_surface_per_bubble_upstream_.resize(nb_bubbles_);
  bubbles_peclet_numbers_from_volume_per_bubble_upstream_.resize(nb_bubbles_);

  for (int dir=0; dir<3; dir++)
    {
      relative_rising_dir_compo_[dir].resize(nb_bubbles_);
      relative_rising_dir_compo_upstream_[dir].resize(nb_bubbles_);
    }

  assert((*bubbles_rising_velocities_).size_array() == nb_bubbles_);
  compute_dynamics_per_bubbles();

  compute_nusselt_numbers_per_bubbles(temperature_gradient_ghost, delta_temperature, lambda);
  compute_shear_per_bubbles();
  compo_found.clear();
}

void IJK_One_Dimensional_Subproblems::compute_overall_bubbles_quantities(IJK_Thermal_Subresolution& ref_thermal_subresolution)
{
  uniform_alpha_ = ref_thermal_subresolution.uniform_alpha_;
  error_temperature_ana_total_ = ref_thermal_subresolution.error_temperature_ana_total_;
  error_temperature_ana_squared_total_ = ref_thermal_subresolution.error_temperature_ana_squared_total_;
  error_temperature_ana_rel_total_ = ref_thermal_subresolution.error_temperature_ana_rel_total_;
  caracteristic_length_ = (ref_thermal_subresolution.single_centred_bubble_radius_ini_) * 2;

  spherical_nusselt_ = ref_thermal_subresolution.nusselt_spherical_diffusion_;
  spherical_nusselt_liquid_ = ref_thermal_subresolution.nusselt_spherical_diffusion_liquid_;

  heat_flux_spherical_ = ref_thermal_subresolution.heat_flux_spherical_;
  mean_liquid_temperature_ = ref_thermal_subresolution.mean_liquid_temperature_;

  velocity_upstream_ = ref_ijk_ft_->get_vitesse_upstream();
  if (velocity_upstream_ < -1.e20)
    velocity_upstream_ = 0.;
  gravity_dir_ = ref_ijk_ft_->get_direction_gravite();
  compute_overall_quantities_per_bubbles(ref_thermal_subresolution.eulerian_grad_T_interface_ns_,
                                         ref_thermal_subresolution.delta_T_subcooled_overheated_,
                                         ref_thermal_subresolution.uniform_lambda_);
  compute_overall_quantities();
  for (auto& itr : *this)
    itr.compute_bubble_related_quantities();
  is_updated_ = true;
}

void IJK_One_Dimensional_Subproblems::compute_dynamics_per_bubbles()
{
  bubbles_rising_relative_velocities_ = (*bubbles_rising_velocities_);
  bubbles_rising_relative_velocities_upstream_ = (*bubbles_rising_velocities_);

  Vecteur3 upstream_velocity_vector;
  upstream_velocity_vector = 0.;
  upstream_velocity_vector[gravity_dir_] = 1.;
  upstream_velocity_vector *= velocity_upstream_;

  for (int i=0; i<(*bubbles_rising_vectors_per_bubble_).dimension(0); i++)
    {
      Vecteur3 rising_velocity_vector;
      rising_velocity_vector = 1.;
      for (int dir=0; dir<3; dir++)
        rising_velocity_vector[dir] = (*bubbles_rising_vectors_per_bubble_)(i, dir);
      rising_velocity_vector *= (*bubbles_rising_velocities_)(i);

      Vecteur3 rising_relative_velocity_vector = (*liquid_velocity_);
      rising_relative_velocity_vector *= (-1);
      rising_relative_velocity_vector += rising_velocity_vector;
      double rising_relative_velocity = rising_relative_velocity_vector.length();
      rising_relative_velocity_vector *= (1 / (rising_relative_velocity + 1e-30));
      for (int dir=0; dir<3; dir++)
        relative_rising_dir_compo_[dir](i) = rising_relative_velocity_vector[dir];

      Vecteur3 rising_relative_velocity_vector_upstream = upstream_velocity_vector;
      rising_relative_velocity_vector_upstream *= (-1);
      rising_relative_velocity_vector_upstream += rising_velocity_vector;
      double rising_relative_velocity_upstream = rising_relative_velocity_vector_upstream.length();
      rising_relative_velocity_vector_upstream *= (1 / (rising_relative_velocity_upstream + 1e-30));
      for (int dir=0; dir<3; dir++)
        relative_rising_dir_compo_upstream_[dir](i) = rising_relative_velocity_vector_upstream[dir];

      bubbles_rising_relative_velocities_(i) = rising_relative_velocity;
      bubbles_reynolds_numbers_per_bubble_(i) = (rising_relative_velocity * caracteristic_length_ * ref_ijk_ft_->get_rho_l()) / ref_ijk_ft_->get_mu_liquid();

      bubbles_rising_relative_velocities_upstream_(i) = rising_relative_velocity_upstream;
      bubbles_reynolds_numbers_per_bubble_upstream_(i) = (rising_relative_velocity_upstream * caracteristic_length_ * ref_ijk_ft_->get_rho_l()) / ref_ijk_ft_->get_mu_liquid();

    }
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

  overall_nusselt_number_error_ = 0.;
  overall_nusselt_number_gfm_error_ = 0.;
  overall_nusselt_number_liquid_error_ = 0.;
  overall_nusselt_number_gfm_liquid_error_ = 0.;

  overall_nusselt_number_error_rel_ = 0.;
  overall_nusselt_number_gfm_error_rel_ = 0.;
  overall_nusselt_number_liquid_error_rel_ = 0.;
  overall_nusselt_number_gfm_liquid_error_rel_ = 0.;

  bubbles_total_rising_velocities_ = 0.;
  bubbles_total_rising_relative_velocities_ = 0.;
  bubbles_total_rising_relative_velocities_upstream_ = 0.;
  bubbles_total_reynolds_numbers_ = 0.;
  bubbles_total_reynolds_numbers_from_surface_per_bubble_ = 0.;
  bubbles_total_reynolds_numbers_from_volume_per_bubble_ = 0.;
  bubbles_total_reynolds_numbers_upstream_ = 0.;
  bubbles_total_reynolds_numbers_from_surface_per_bubble_upstream_ = 0.;
  bubbles_total_reynolds_numbers_from_volume_per_bubble_upstream_ = 0.;

  bubbles_total_peclet_numbers_ = 0.;
  bubbles_total_peclet_numbers_from_surface_per_bubble_ = 0.;
  bubbles_total_peclet_numbers_from_volume_per_bubble_ = 0.;
  bubbles_total_peclet_numbers_upstream_  = 0.;
  bubbles_total_peclet_numbers_from_surface_per_bubble_upstream_ = 0.;
  bubbles_total_peclet_numbers_from_volume_per_bubble_upstream_  = 0.;

  total_relative_rising_dir_compo_ = 0.;
  total_relative_rising_dir_compo_upstream_ = 0.;
  for (int i=0; i < nb_bubbles_; i++)
    {
      interfacial_thermal_flux_ += interfacial_thermal_flux_per_bubble_(i);
      interfacial_thermal_flux_gfm_ += interfacial_thermal_flux_per_bubble_gfm_(i);
      interfacial_thermal_flux_per_bubble_spherical_(i) = heat_flux_spherical_ * total_surface_per_bubble_(i);


      overall_nusselt_number_ += overall_nusselt_number_per_bubble_(i);
      overall_nusselt_number_gfm_ += overall_nusselt_number_per_bubble_gfm_(i);
      overall_nusselt_number_liquid_ += overall_nusselt_number_per_bubble_liquid_(i);
      overall_nusselt_number_gfm_liquid_ += overall_nusselt_number_per_bubble_gfm_liquid_(i);

      overall_shear_force_ += overall_shear_force_per_bubble_(i);

      total_surface_ += total_surface_per_bubble_(i);
      total_volume_ += (*bubbles_volume_)(i);

      radius_from_surfaces_per_bubble_(i) = sqrt(total_surface_per_bubble_(i) / (4. * M_PI));
      radius_from_volumes_per_bubble_(i) = pow((*bubbles_volume_)(i) * 3. / (4. * M_PI), (1. / 3.));

      caracteristic_length_from_surfaces_per_bubble_(i) = radius_from_surfaces_per_bubble_(i) * 2;
      caracteristic_length_from_volumes_per_bubble_(i) = radius_from_volumes_per_bubble_(i) * 2;

      overall_nusselt_number_per_bubble_spherical_(i) = (heat_flux_spherical_ * caracteristic_length_from_surfaces_per_bubble_(i))
                                                        / (delta_temperature_ * lambda_);
      overall_nusselt_number_per_bubble_spherical_liquid_(i) = (heat_flux_spherical_ * caracteristic_length_from_surfaces_per_bubble_(i))
                                                               / (mean_liquid_temperature_ * lambda_);

      overall_nusselt_number_per_bubble_error_(i) = overall_nusselt_number_per_bubble_(i) - spherical_nusselt_;
      overall_nusselt_number_per_bubble_gfm_error_(i) = overall_nusselt_number_per_bubble_gfm_(i) - spherical_nusselt_;
      overall_nusselt_number_per_bubble_liquid_error_(i) = overall_nusselt_number_per_bubble_liquid_(i) - spherical_nusselt_;
      overall_nusselt_number_per_bubble_gfm_liquid_error_(i) = overall_nusselt_number_per_bubble_gfm_liquid_(i) - spherical_nusselt_;

      overall_nusselt_number_per_bubble_error_rel_(i) = overall_nusselt_number_per_bubble_error_(i) / (1e-16 + spherical_nusselt_);
      overall_nusselt_number_per_bubble_gfm_error_rel_(i) = overall_nusselt_number_per_bubble_gfm_error_(i) / (1e-16 + spherical_nusselt_);
      overall_nusselt_number_per_bubble_liquid_error_rel_(i) = overall_nusselt_number_per_bubble_liquid_error_(i) / (1e-16 + spherical_nusselt_);
      overall_nusselt_number_per_bubble_gfm_liquid_error_rel_(i) = overall_nusselt_number_per_bubble_gfm_liquid_error_(i) / (1e-16 + spherical_nusselt_);

      bubbles_total_reynolds_numbers_ += bubbles_reynolds_numbers_per_bubble_(i) * (*bubbles_volume_)(i);
      bubbles_total_reynolds_numbers_upstream_ += bubbles_reynolds_numbers_per_bubble_upstream_(i) * (*bubbles_volume_)(i);

      bubbles_reynolds_numbers_from_surface_per_bubble_(i) = (bubbles_rising_relative_velocities_(i)
                                                              * caracteristic_length_from_surfaces_per_bubble_(i) * ref_ijk_ft_->get_rho_l()) / ref_ijk_ft_->get_mu_liquid();
      bubbles_reynolds_numbers_from_volume_per_bubble_(i) = (bubbles_rising_relative_velocities_(i)
                                                             * caracteristic_length_from_volumes_per_bubble_(i) * ref_ijk_ft_->get_rho_l()) / ref_ijk_ft_->get_mu_liquid();
      bubbles_reynolds_numbers_from_surface_per_bubble_ += bubbles_reynolds_numbers_from_surface_per_bubble_(i) * (*bubbles_volume_)(i);
      bubbles_reynolds_numbers_from_volume_per_bubble_ += bubbles_reynolds_numbers_from_volume_per_bubble_(i) * (*bubbles_volume_)(i);

      bubbles_reynolds_numbers_from_surface_per_bubble_upstream_(i) = (bubbles_rising_relative_velocities_upstream_(i)
                                                                       * caracteristic_length_from_surfaces_per_bubble_(i) * ref_ijk_ft_->get_rho_l()) / ref_ijk_ft_->get_mu_liquid();
      bubbles_reynolds_numbers_from_volume_per_bubble_upstream_(i) = (bubbles_rising_relative_velocities_upstream_(i)
                                                                      * caracteristic_length_from_volumes_per_bubble_(i) * ref_ijk_ft_->get_rho_l()) / ref_ijk_ft_->get_mu_liquid();
      bubbles_reynolds_numbers_from_surface_per_bubble_upstream_ += bubbles_reynolds_numbers_from_surface_per_bubble_upstream_(i) * (*bubbles_volume_)(i);
      bubbles_reynolds_numbers_from_volume_per_bubble_upstream_ += bubbles_reynolds_numbers_from_volume_per_bubble_upstream_(i) * (*bubbles_volume_)(i);

      bubbles_total_rising_velocities_ += (*bubbles_rising_velocities_)(i) * (*bubbles_volume_)(i);
      bubbles_total_rising_relative_velocities_ += bubbles_rising_relative_velocities_(i) * (*bubbles_volume_)(i);
      bubbles_total_rising_relative_velocities_upstream_ += bubbles_rising_relative_velocities_upstream_(i) * (*bubbles_volume_)(i);
      for (int dir=0; dir<3; dir++)
        {
          total_rising_dir_compo_[dir] = (*bubbles_rising_vectors_per_bubble_)(i, dir) * (*bubbles_volume_)(i);
          total_relative_rising_dir_compo_[dir] = relative_rising_dir_compo_[dir](i) * (*bubbles_volume_)(i);
          total_relative_rising_dir_compo_upstream_[dir] = relative_rising_dir_compo_upstream_[dir](i) * (*bubbles_volume_)(i);
        }
      bubbles_peclet_numbers_per_bubble_(i) = bubbles_reynolds_numbers_per_bubble_(i) * (*prandtl_number_);
      bubbles_peclet_numbers_from_surface_per_bubble_(i) = bubbles_reynolds_numbers_from_surface_per_bubble_(i) * (*prandtl_number_);
      bubbles_peclet_numbers_from_volume_per_bubble_(i) = bubbles_reynolds_numbers_from_volume_per_bubble_(i) * (*prandtl_number_);

      bubbles_total_peclet_numbers_ += bubbles_peclet_numbers_per_bubble_(i) * total_surface_per_bubble_(i);
      bubbles_total_peclet_numbers_from_surface_per_bubble_ += bubbles_peclet_numbers_from_surface_per_bubble_(i) * total_surface_per_bubble_(i);
      bubbles_total_peclet_numbers_from_volume_per_bubble_ += bubbles_peclet_numbers_from_volume_per_bubble_(i) * total_surface_per_bubble_(i);

      bubbles_peclet_numbers_per_bubble_upstream_(i) = bubbles_reynolds_numbers_per_bubble_upstream_(i) * (*prandtl_number_);
      bubbles_peclet_numbers_from_surface_per_bubble_upstream_(i) = bubbles_reynolds_numbers_from_surface_per_bubble_upstream_(i) * (*prandtl_number_);
      bubbles_peclet_numbers_from_volume_per_bubble_upstream_(i) = bubbles_reynolds_numbers_from_volume_per_bubble_upstream_(i) * (*prandtl_number_);

      bubbles_total_peclet_numbers_upstream_ += bubbles_peclet_numbers_per_bubble_upstream_(i) * total_surface_per_bubble_(i);
      bubbles_total_peclet_numbers_from_surface_per_bubble_upstream_ += bubbles_peclet_numbers_from_surface_per_bubble_upstream_(i) * total_surface_per_bubble_(i);
      bubbles_total_peclet_numbers_from_volume_per_bubble_upstream_ += bubbles_peclet_numbers_from_volume_per_bubble_upstream_(i) * total_surface_per_bubble_(i);

    }
  bubbles_total_rising_velocities_ /= total_volume_;
  bubbles_total_rising_relative_velocities_ /= total_volume_;
  bubbles_total_rising_relative_velocities_upstream_ /= total_volume_;

  bubbles_total_reynolds_numbers_ /= total_volume_;
  bubbles_total_reynolds_numbers_from_surface_per_bubble_ /= total_volume_;
  bubbles_total_reynolds_numbers_from_volume_per_bubble_ /= total_volume_;

  bubbles_total_reynolds_numbers_upstream_ /= total_volume_;
  bubbles_total_reynolds_numbers_from_surface_per_bubble_upstream_ /= total_volume_;
  bubbles_total_reynolds_numbers_from_volume_per_bubble_upstream_ /= total_volume_;

  const double total_rising_dir_norm = total_rising_dir_compo_.length();
  total_rising_dir_compo_ *= (1. / (total_rising_dir_norm + 1e-30));
  const double total_relative_rising_dir_norm = total_relative_rising_dir_compo_.length();
  total_relative_rising_dir_compo_ *= (1. / (total_relative_rising_dir_norm + 1e-30));
  const double total_relative_rising_dir_norm_upstream = total_relative_rising_dir_compo_upstream_.length();
  total_relative_rising_dir_compo_upstream_ *= (1. / (total_relative_rising_dir_norm_upstream + 1e-30));

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

  overall_nusselt_number_error_ = (overall_nusselt_number_ - spherical_nusselt_);
  overall_nusselt_number_gfm_error_ = (overall_nusselt_number_gfm_ - spherical_nusselt_);
  overall_nusselt_number_liquid_error_ = (overall_nusselt_number_liquid_ - spherical_nusselt_);
  overall_nusselt_number_gfm_liquid_error_ = (overall_nusselt_number_gfm_liquid_ - spherical_nusselt_);

  overall_nusselt_number_error_rel_ = overall_nusselt_number_error_ / (1e-16 + spherical_nusselt_);
  overall_nusselt_number_gfm_error_rel_ = overall_nusselt_number_gfm_error_ / (1e-16 + spherical_nusselt_);
  overall_nusselt_number_liquid_error_rel_ = overall_nusselt_number_liquid_error_ / (1e-16 + spherical_nusselt_);
  overall_nusselt_number_gfm_liquid_error_rel_ = overall_nusselt_number_gfm_liquid_error_ / (1e-16 + spherical_nusselt_);

  bubbles_total_peclet_numbers_ /= total_surface_;
  bubbles_total_peclet_numbers_from_surface_per_bubble_ /= total_surface_;
  bubbles_total_peclet_numbers_from_volume_per_bubble_ /= total_surface_;

  bubbles_total_peclet_numbers_upstream_ /= total_surface_;
  bubbles_total_peclet_numbers_from_surface_per_bubble_upstream_ /= total_surface_;
  bubbles_total_peclet_numbers_from_volume_per_bubble_upstream_ /= total_surface_;
}

void IJK_One_Dimensional_Subproblems::post_process_overall_bubbles_quantities(const int rank, const Nom& overall_bubbles_quantities)
{
  if (Process::je_suis_maitre())
    {
      const int reset = 1;
      const int last_time_index = (*latastep_reprise_) + ref_ijk_ft_->get_tstep();
      const int max_digit = 3;
      const int max_digit_time = 8;
      const int max_rank_digit = rank < 1 ? 1 : (int) (log10(rank) + 1);
      const int nb_digit_tstep = last_time_index < 1 ? 1 : (int) (log10(last_time_index) + 1);

      Nom probe_name = Nom("_thermal_rank_") + Nom(std::string(max_digit - max_rank_digit,'0')) + Nom(rank)
                       + Nom("_thermal_subproblems") + ("_overall_bubbles_quantities_")
                       + Nom(std::string(max_digit_time - nb_digit_tstep, '0')) + Nom(last_time_index) + Nom(".out");
      Nom probe_header = Nom("tstep\ttime\tthermal_rank\tbubble_index"
                             "\ttime_dimensionless"
                             "\tnusselt_overall\tnusselt_overall_gfm\tnusselt_spherical\tnusselt_spherical_th"
                             "\tnusselt_overall_liq\tnusselt_overall_gfm_liq\tnusselt_spherical_liq\tnusselt_spherical_th_liq"
                             "\tnusselt_overall_error\tnusselt_overall_gfm_error\tnusselt_overall_liq_error\tnusselt_overall_gfm_liq_error"
                             "\tnusselt_overall_error_rel\tnusselt_overall_gfm_error_rel\tnusselt_overall_liq_error_rel\tnusselt_overall_gfm_liq_error_rel"
                             "\theat_flux\theat_flux_gfm\theat_flux_spherical"
                             "\ttotal_surface\ttotal_volume"
                             "\tradius_surface\tradius_volume"
                             "\terror_temperature_ana\terror_temperature_ana_norm\terror_temperature_ana_rel"
                             "\tliquid_velocity"
                             "\tliquid_velocity_x\tliquid_velocity_y\tliquid_velocity_z"
                             "\tgravity_dir\tupstream_velocity"
                             "\trising_dir_x\trising_dir_y\trising_dir_z"
                             "\trelative_rising_dir_x\trelative_rising_dir_y\trelative_rising_dir_z"
                             "\trelative_rising_dir_upstream_x\trelative_rising_dir_upstream_y\trelative_rising_dir_upstream_z"
                             "\trising_velocity\trising_relative_velocity\trising_relative_velocity_upstream"
                             "\treynolds_number\treynolds_number_surface\treynolds_number_volume"
                             "\treynolds_number_upstream\treynolds_number_surface_upstream\treynolds_number_volume_upstream"
                             "\tprandtl_liq"
                             "\tpeclet_number\tpeclet_number_surface\tpeclet_number_volume"
                             "\tpeclet_number_upstream\tpeclet_number_surface_upstream\tpeclet_number_volume_upstream");
      SFichier fic = Open_file_folder(overall_bubbles_quantities, probe_name, probe_header, reset);
      int max_counter = nb_bubbles_;
      const double last_time = ref_ijk_ft_->get_current_time() - ref_ijk_ft_->get_timestep();
      const double dimensionless_time = caracteristic_length_ / (1e-16 + sqrt(M_PI * last_time * uniform_alpha_));
      /*
       * TODO: fill the Array in parallel
       */
      for (int i=0; i<max_counter; i++)
        {
          fic << last_time_index << " " << last_time << " ";
          fic << rank << " ";
          fic << i << " ";
          fic << dimensionless_time << " ";
          fic << overall_nusselt_number_per_bubble_(i) << " ";
          fic << overall_nusselt_number_per_bubble_gfm_(i) << " ";
          fic << overall_nusselt_number_per_bubble_spherical_(i) << " ";
          fic << spherical_nusselt_ << " ";
          fic << overall_nusselt_number_per_bubble_liquid_(i) << " ";
          fic << overall_nusselt_number_per_bubble_gfm_liquid_(i) << " ";
          fic << overall_nusselt_number_per_bubble_spherical_liquid_(i) << " ";
          fic << spherical_nusselt_liquid_ << " ";
          fic << overall_nusselt_number_per_bubble_error_(i) << " ";
          fic << overall_nusselt_number_per_bubble_gfm_error_(i) << " ";
          fic << overall_nusselt_number_per_bubble_liquid_error_(i)  << " ";
          fic << overall_nusselt_number_per_bubble_gfm_liquid_error_(i) << " ";
          fic << overall_nusselt_number_per_bubble_error_rel_(i)  << " ";
          fic << overall_nusselt_number_per_bubble_gfm_error_rel_(i) << " ";
          fic << overall_nusselt_number_per_bubble_liquid_error_rel_(i)  << " ";
          fic << overall_nusselt_number_per_bubble_gfm_liquid_error_rel_(i)  << " ";
          fic << interfacial_thermal_flux_per_bubble_(i) << " ";
          fic << interfacial_thermal_flux_per_bubble_gfm_(i) << " ";
          fic << interfacial_thermal_flux_per_bubble_spherical_(i) << " ";
          fic << total_surface_per_bubble_(i) << " ";
          fic << (*bubbles_volume_)(i) << " ";
          fic << radius_from_surfaces_per_bubble_(i) << " ";
          fic << radius_from_volumes_per_bubble_(i) << " ";
          fic << error_temperature_ana_total_ << " ";
          fic << error_temperature_ana_squared_total_ << " ";
          fic << error_temperature_ana_rel_total_ << " ";
          fic << (*liquid_velocity_).length() << " ";
          fic << (*liquid_velocity_)[0] << " ";
          fic << (*liquid_velocity_)[1] << " ";
          fic << (*liquid_velocity_)[2] << " ";
          fic << gravity_dir_ << " ";
          fic << velocity_upstream_ << " ";
          fic << (*bubbles_rising_vectors_per_bubble_)(i, 0) << " ";
          fic << (*bubbles_rising_vectors_per_bubble_)(i, 1) << " ";
          fic << (*bubbles_rising_vectors_per_bubble_)(i, 2) << " ";
          fic << relative_rising_dir_compo_[0](i) << " ";
          fic << relative_rising_dir_compo_[1](i) << " ";
          fic << relative_rising_dir_compo_[2](i) << " ";
          fic << relative_rising_dir_compo_upstream_[0](i) << " ";
          fic << relative_rising_dir_compo_upstream_[1](i) << " ";
          fic << relative_rising_dir_compo_upstream_[2](i) << " ";
          fic << (*bubbles_rising_velocities_)(i) << " ";
          fic << bubbles_rising_relative_velocities_(i) << " ";
          fic << bubbles_rising_relative_velocities_upstream_(i) << " ";
          fic << bubbles_reynolds_numbers_per_bubble_(i) << " ";
          fic << bubbles_reynolds_numbers_from_surface_per_bubble_(i) << " ";
          fic << bubbles_reynolds_numbers_from_volume_per_bubble_(i) << " ";
          fic << bubbles_reynolds_numbers_per_bubble_upstream_(i) << " ";
          fic << bubbles_reynolds_numbers_from_surface_per_bubble_upstream_(i) << " ";
          fic << bubbles_reynolds_numbers_from_volume_per_bubble_upstream_(i) << " ";
          fic << (*prandtl_number_) << " ";
          fic << bubbles_peclet_numbers_per_bubble_(i) << " ";
          fic << bubbles_peclet_numbers_from_surface_per_bubble_(i) << " ";
          fic << bubbles_peclet_numbers_from_volume_per_bubble_(i) << " ";
          fic << bubbles_peclet_numbers_per_bubble_upstream_(i) << " ";
          fic << bubbles_peclet_numbers_from_surface_per_bubble_upstream_(i) << " ";
          fic << bubbles_peclet_numbers_from_volume_per_bubble_upstream_(i) << " ";
          fic << finl;
        }
      /*
       * Should be good for parallel
       */
      if(max_counter > 1)
        {
          fic << last_time_index << " " << last_time << " ";
          fic << rank << " ";
          fic << nb_bubbles_ << " ";
          fic << dimensionless_time << " ";
          fic << overall_nusselt_number_ << " ";
          fic << overall_nusselt_number_gfm_ << " ";
          fic << overall_nusselt_number_spherical_ << " ";
          fic << spherical_nusselt_ << " ";
          fic << overall_nusselt_number_liquid_ << " ";
          fic << overall_nusselt_number_gfm_liquid_ << " ";
          fic << overall_nusselt_number_spherical_liquid_ << " ";
          fic << spherical_nusselt_liquid_ << " ";
          fic << overall_nusselt_number_error_ << " ";
          fic << overall_nusselt_number_gfm_error_ << " ";
          fic << overall_nusselt_number_liquid_error_  << " ";
          fic << overall_nusselt_number_gfm_liquid_error_  << " ";
          fic << overall_nusselt_number_error_rel_  << " ";
          fic << overall_nusselt_number_gfm_error_rel_  << " ";
          fic << overall_nusselt_number_liquid_error_rel_  << " ";
          fic << overall_nusselt_number_gfm_liquid_error_rel_  << " ";
          fic << interfacial_thermal_flux_ << " ";
          fic << interfacial_thermal_flux_gfm_ << " ";
          fic << heat_flux_spherical_ << " ";
          fic << total_surface_ << " ";
          fic << total_volume_ << " ";
          fic << radius_from_surfaces_ << " ";
          fic << radius_from_volumes_ << " ";
          fic << error_temperature_ana_total_ << " ";
          fic << error_temperature_ana_squared_total_ << " ";
          fic << error_temperature_ana_rel_total_ << " ";
          fic << (*liquid_velocity_).length() << " ";
          fic << (*liquid_velocity_)[0] << " ";
          fic << (*liquid_velocity_)[1] << " ";
          fic << (*liquid_velocity_)[2] << " ";
          fic << gravity_dir_ << " ";
          fic << velocity_upstream_ << " ";
          fic << total_rising_dir_compo_[0] << " ";
          fic << total_rising_dir_compo_[1] << " ";
          fic << total_rising_dir_compo_[2] << " ";
          fic << total_relative_rising_dir_compo_[0] << " ";
          fic << total_relative_rising_dir_compo_[1] << " ";
          fic << total_relative_rising_dir_compo_[2] << " ";
          fic << total_relative_rising_dir_compo_upstream_[0] << " ";
          fic << total_relative_rising_dir_compo_upstream_[1] << " ";
          fic << total_relative_rising_dir_compo_upstream_[2] << " ";
          fic << bubbles_total_rising_velocities_ << " ";
          fic << bubbles_total_rising_relative_velocities_ << " ";
          fic << bubbles_total_rising_relative_velocities_upstream_ << " ";
          fic << bubbles_total_reynolds_numbers_ << " ";
          fic << bubbles_reynolds_numbers_from_surface_per_bubble_ << " ";
          fic << bubbles_reynolds_numbers_from_volume_per_bubble_ << " ";
          fic << bubbles_total_reynolds_numbers_upstream_ << " ";
          fic << bubbles_reynolds_numbers_from_surface_per_bubble_upstream_ << " ";
          fic << bubbles_reynolds_numbers_from_volume_per_bubble_upstream_ << " ";
          fic << (*prandtl_number_) << finl;
          fic << bubbles_total_peclet_numbers_ << " ";
          fic << bubbles_peclet_numbers_from_surface_per_bubble_ << " ";
          fic << bubbles_peclet_numbers_from_volume_per_bubble_ << " ";
          fic << bubbles_total_peclet_numbers_upstream_ << " ";
          fic << bubbles_peclet_numbers_from_surface_per_bubble_upstream_ << " ";
          fic << bubbles_peclet_numbers_from_volume_per_bubble_upstream_ << " ";
          fic << finl;
        }
      fic.close();
    }
}


