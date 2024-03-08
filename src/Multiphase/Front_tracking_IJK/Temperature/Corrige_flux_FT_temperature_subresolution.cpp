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
// File      : Corrige_flux_FT_temperature_subresolution.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <Corrige_flux_FT_temperature_subresolution.h>
#include <IJK_FT.h>

Implemente_instanciable_sans_constructeur( Corrige_flux_FT_temperature_subresolution, "Corrige_flux_FT_temperature_subresolution", Corrige_flux_FT_base ) ;

Corrige_flux_FT_temperature_subresolution::Corrige_flux_FT_temperature_subresolution()
{
  thermal_subproblems_ = nullptr;
  has_checked_consistency_ = false;
  distance_cell_faces_from_lrs_= 0;
  find_temperature_cell_neighbours_ = 0;
  neighbours_colinearity_weighting_ = 0;
  find_cell_neighbours_for_fluxes_spherical_correction_ = 0;
  use_cell_neighbours_for_fluxes_spherical_correction_ = 0;
  find_reachable_fluxes_ = 0;
  use_reachable_fluxes_ = 0;
  keep_first_reachable_fluxes_ = 0.;
  cell_faces_neighbours_corrected_bool_ = nullptr;
  eulerian_normal_vectors_ns_normed_ = nullptr;

  convection_negligible_ = 0;
  diffusion_negligible_ = 0;
  debug_=0;
  levels_=0;
  discrete_integral_=0;
  flux_init_ = 0;
  convective_flux_correction_ = 0;
  diffusive_flux_correction_ = 0;
  smooth_temperature_field_=0;

  copy_fluxes_on_every_procs_ = 1;
  copy_temperature_on_every_procs_ = 1;

  for (int c=0; c<3; c++)
    indices_temperature_neighbours_on_procs_[c].set_smart_resize(1);
  temperature_neighbours_on_procs_.set_smart_resize(1);
  neighbours_weighting_colinearity_on_procs_.set_smart_resize(1);
}

void Corrige_flux_FT_temperature_subresolution::associate_thermal_problems(const IJK_One_Dimensional_Subproblems& thermal_subproblems)
{
  thermal_subproblems_ = &thermal_subproblems;
}

Sortie& Corrige_flux_FT_temperature_subresolution::printOn( Sortie& os ) const
{
  Corrige_flux_FT_base::printOn( os );
  return os;
}

Entree& Corrige_flux_FT_temperature_subresolution::readOn( Entree& is )
{
  Corrige_flux_FT_base::readOn( is );
  return is;
}

void Corrige_flux_FT_temperature_subresolution::clear_vectors()
{
  int c, l;

  for (c=0; c<3; c++)
    indices_temperature_neighbours_on_procs_[c].reset();
  temperature_neighbours_on_procs_.reset();
  if (neighbours_colinearity_weighting_)
    neighbours_weighting_colinearity_on_procs_.reset();

  for (l=0; l<2; l++)
    for (c=0; c<3; c++)
      clear_std_vectors_array_of_int(index_face_ij_flux_xyz_sorted_[l][c]);

  for (c=0; c<3; c++)
    {
      if (convective_flux_correction_)
        clear_std_vectors_array_of_double(convective_diffusive_flux_xyz_sorted_[0][c]);
      if (diffusive_flux_correction_)
        clear_std_vectors_array_of_double(convective_diffusive_flux_xyz_sorted_[1][c]);
    }

  for (l=0; l<2; l++)
    for (c=0; c<3; c++)
      clear_std_vectors_array_of_int(index_face_ij_flux_xyz_remaining_global_sorted_[l][c]);

  for (c=0; c<3; c++)
    {
      if (convective_flux_correction_)
        clear_std_vectors_array_of_double(convective_diffusive_flux_xyz_remaining_global_sorted_[0][c]);
      if (diffusive_flux_correction_)
        clear_std_vectors_array_of_double(convective_diffusive_flux_xyz_remaining_global_sorted_[1][c]);
    }

  for (c=0; c<3; c++)
    flux_frontier_map_[c].clear();

  /*
   * M.G (07/12/23) All fluxes may be useless; Diagonal fluxes are not relevant for the moment
   */

  for (l=0; l<2; l++)
    for (c=0; c<3; c++)
      clear_std_vectors_array_of_int(index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[l][c]);

  for (l=0; l<2; l++)
    for (c=0; c<3; c++)
      clear_std_vectors_array_of_int(index_face_ij_flux_xyz_neighbours_all_faces_sorted_[l][c]);

  /*
   * Face fluxes on a reconstructed convex shell around the bubble !
   */

  for (l=0; l<2; l++)
    for (c=0; c<3; c++)
      clear_std_vectors_array_of_int(index_face_ij_flux_xyz_neighbours_min_max_faces_sorted_[l][c]);

  for (c=0; c<3; c++)
    {
      if (convective_flux_correction_)
        clear_std_vectors_array_of_double(convective_diffusive_flux_xyz_min_max_faces_sorted_[0][c]);
      if (diffusive_flux_correction_)
        clear_std_vectors_array_of_double(convective_diffusive_flux_xyz_min_max_faces_sorted_[1][c]);
    }

  /*
   * Face fluxes on a reconstructed convex shell around the bubble (Parallel) !
   */

  for (l=0; l<2; l++)
    for (c=0; c<3; c++)
      clear_std_vectors_array_of_int(index_face_ij_flux_xyz_neighbours_all_faces_remaining_global_sorted_[l][c]);

  for (c=0; c<3; c++)
    clear_std_vectors_array_of_int(weighting_flux_xyz_neighbours_all_faces_remaining_global_sorted_[c]);


  for (c=0; c<3; c++)
    clear_std_vectors_array_of_double(colinearity_flux_xyz_neighbours_all_faces_remaining_global_sorted_[c]);

  for (c=0; c<3; c++)
    {
      if (convective_flux_correction_)
        clear_std_vectors_array_of_double(convective_diffusive_flux_all_faces_remaining_global_sorted_[0][c]);
      if (diffusive_flux_correction_)
        clear_std_vectors_array_of_double(convective_diffusive_flux_all_faces_remaining_global_sorted_[1][c]);
    }

  for (c=0; c<3; c++)
    {
      flux_outside_frontier_all_map_[c].clear();
      flux_frontier_all_map_[c].clear();
    }
}

void Corrige_flux_FT_temperature_subresolution::clear_std_vectors_array_of_int(std::vector<ArrOfInt>& indices_to_clear)
{
  for (int k=0; k<(int) indices_to_clear.size(); k++)
    indices_to_clear[k].reset();
}
void Corrige_flux_FT_temperature_subresolution::clear_std_vectors_array_of_double(std::vector<ArrOfDouble>& values_to_clear)
{
  for (int k=0; k<(int) values_to_clear.size(); k++)
    values_to_clear[k].reset();
}

int compute_periodic_index(const int index, const int n)
{
  return (n + index % n) % n;
}

void Corrige_flux_FT_temperature_subresolution::initialize_with_subproblems(const IJK_Splitting& splitting,
                                                                            const IJK_Field_double& field,
                                                                            const IJK_Interfaces& interfaces,
                                                                            const IJK_FT_double& ijk_ft,
                                                                            Intersection_Interface_ijk_face& intersection_ijk_face,
                                                                            Intersection_Interface_ijk_cell& intersection_ijk_cell,
                                                                            const IJK_One_Dimensional_Subproblems& thermal_subproblems)
{
  Corrige_flux_FT_base::initialize_with_subproblems(splitting, field, interfaces, ijk_ft, intersection_ijk_face, intersection_ijk_cell, thermal_subproblems);
  associate_thermal_problems(thermal_subproblems);
}

void Corrige_flux_FT_temperature_subresolution::update_intersections()
{

  // if (!distance_cell_faces_from_lrs_)
  {
    // On commence par calculer les temperatures aux faces mouillées
    intersection_ijk_cell_->update_interpolations_cell_centres_on_interface();
    /*
     * TODO update with face cell centres positions
     */
    if (!convection_negligible_ || !diffusion_negligible_)
      intersection_ijk_cell_->update_interpolations_cell_faces_on_interface();

    Cerr << "The intersections have been updated" << finl;
  }
}

void Corrige_flux_FT_temperature_subresolution::update()
{
  associate_indices_and_check_subproblems_consistency();
}

void Corrige_flux_FT_temperature_subresolution::associate_indices_and_check_subproblems_consistency()
{
  if (!has_checked_consistency_)
    {
      const int nb_diph = intersection_ijk_cell_->get_nb_diph();
      const int nb_subproblems = thermal_subproblems_->get_subproblems_counter();
      const int nb_effective_subproblems = thermal_subproblems_->get_effective_subproblems_counter();
      has_checked_consistency_ = (nb_diph==nb_subproblems);
      assert(has_checked_consistency_);
      ijk_intersections_subproblems_indices_.reset();
      ijk_intersections_subproblems_indices_.resize(nb_effective_subproblems);
      int index_i_problem = 0;
      int index_j_problem = 0;
      int index_k_problem = 0;
      int problem_index;
      for (problem_index=0; problem_index<nb_effective_subproblems; problem_index++)
        {
          thermal_subproblems_->get_subproblem_ijk_indices(index_i_problem, index_j_problem, index_k_problem, problem_index);
          const int ijk_intersections_index = (*intersection_ijk_cell_)(index_i_problem, index_j_problem, index_k_problem);
          if (ijk_intersections_index > -1)
            {
              ijk_intersections_subproblems_indices_[problem_index] = ijk_intersections_index;
              has_checked_consistency_ = has_checked_consistency_ && true;
            }
          else
            {
              Cerr << "Inconsistency between intersection_ijk_cell and the thermal_subproblem (index " << problem_index << ")" << finl;
              has_checked_consistency_ = has_checked_consistency_ && false;
            }
        }
      assert(has_checked_consistency_);
      if (!has_checked_consistency_)
        {
          Cerr << "There is an inconsistency between the LRS sub-problem and the intersection_ijk_elem" << finl;
          Process::exit();
        }
    }
  else
    Cerr << "Inconsistency has already be checked" << finl;
}

void Corrige_flux_FT_temperature_subresolution::clean()
{
  // if (!distance_cell_faces_from_lrs_)
  {
    has_checked_consistency_=false;
    intersection_ijk_cell_->set_pas_a_jour();
  }
}


void Corrige_flux_FT_temperature_subresolution::compute_temperature_cell_centre(IJK_Field_double& temperature) const
{
  /*
   * For each subproblem fill the right interfacial_cell
   */
  DoubleTab dist_interf;
  if (!distance_cell_faces_from_lrs_)
    dist_interf = intersection_ijk_cell_->dist_interf();

  const double min_temperature = thermal_subproblems_->get_min_temperature_domain_ends();
  const double max_temperature = thermal_subproblems_->get_max_temperature_domain_ends();

  for (int i=0; i<thermal_subproblems_->get_effective_subproblems_counter(); i++)
    {
      double dist = 0;
      double dist_sub_res = 0;
      double temperature_ghost = 0.;
      int intersection_ijk_cell_index = 0;
      int ijk_indices_i = 0;
      int ijk_indices_j = 0;
      int ijk_indices_k = 0;
      dist_sub_res = thermal_subproblems_->get_dist_cell_interface(i);
      if (!distance_cell_faces_from_lrs_)
        {
          intersection_ijk_cell_index = ijk_intersections_subproblems_indices_[i];
          dist = dist_interf(intersection_ijk_cell_index, 0);
          temperature_ghost = thermal_subproblems_->get_temperature_profile_at_point(i, dist);
          ijk_indices_i = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 0);
          ijk_indices_j = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 1);
          ijk_indices_k = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 2);
        }
      else
        {
          temperature_ghost = thermal_subproblems_->get_temperature_profile_at_point(i, dist_sub_res);
          thermal_subproblems_->get_subproblem_ijk_indices(ijk_indices_i, ijk_indices_j, ijk_indices_k, i);
        }

      temperature(ijk_indices_i, ijk_indices_j, ijk_indices_k) = temperature_ghost;

      if (debug_)
        {
          Cerr << "Time-step: " << ref_ijk_ft_->get_tstep() << "--" << ref_ijk_ft_->get_timestep() << " s" << finl;
          if (!distance_cell_faces_from_lrs_)
            {
              Cerr << "Distance at cell : " << intersection_ijk_cell_index <<
                   ", subproblem " << i << "."<<
                   " -- (" << ijk_indices_i << ", " << ijk_indices_j << ", " << ijk_indices_k << ")" << finl;
              Cerr << "Distance from intersection_ijk_cell: " << dist << finl;
            }
          Cerr << "Distance from sub-resolution: " << dist_sub_res << finl;
          Vecteur3 bary_facet_debug = thermal_subproblems_->get_bary_facet(i);
          Cerr << "Facet barycentre: " << bary_facet_debug[0] << ";"
               << bary_facet_debug[1] << ";"
               << bary_facet_debug[2] << finl;

          const IJK_Field_double& indicator = ref_ijk_ft_->itfce().I();
          const double indic = indicator(ijk_indices_i, ijk_indices_j, ijk_indices_k);
          if (temperature_ghost < min_temperature && indic > 0.5)
            Cerr << "Ghost temperature: " << temperature_ghost << " is lower than the minimum temperature:" << min_temperature << finl;
          if (temperature_ghost > max_temperature && indic > 0.5)
            Cerr << "Ghost temperature: " << temperature_ghost << " is higher than the maximum temperature:" << max_temperature << finl;
        }
    }
}

void Corrige_flux_FT_temperature_subresolution::compute_temperature_cell_centre_neighbours(IJK_Field_double& temperature_neighbours,
                                                                                           IJK_Field_int& neighbours_weighting,
                                                                                           IJK_Field_double& neighbours_weighting_colinearity)
{
  if (distance_cell_faces_from_lrs_ && find_temperature_cell_neighbours_)
    {

      const int ni = neighbours_weighting.ni();
      const int nj = neighbours_weighting.nj();
      const int nk = neighbours_weighting.nk();

      const int offset_i = neighbours_weighting.get_splitting().get_offset_local(0);
      const int offset_j = neighbours_weighting.get_splitting().get_offset_local(1);
      const int offset_k = neighbours_weighting.get_splitting().get_offset_local(2);

      const int ni_tot = neighbours_weighting.get_splitting().get_grid_geometry().get_nb_elem_tot(0);
      const int nj_tot = neighbours_weighting.get_splitting().get_grid_geometry().get_nb_elem_tot(1);
      const int nk_tot = neighbours_weighting.get_splitting().get_grid_geometry().get_nb_elem_tot(2);

      neighbours_weighting.data() = 0;
      neighbours_weighting.echange_espace_virtuel(neighbours_weighting.ghost());
      temperature_neighbours.data() = 0.;
      temperature_neighbours.echange_espace_virtuel(temperature_neighbours.ghost());
      if (neighbours_colinearity_weighting_)
        {
          neighbours_weighting_colinearity.data() = 0.;
          neighbours_weighting_colinearity.echange_espace_virtuel(neighbours_weighting_colinearity.ghost());
        }
      std::vector<std::vector<std::vector<double>>> pure_neighbours_corrected_colinearity;
      int index_i_problem, index_j_problem, index_k_problem;
      int index_i_neighbour, index_j_neighbour, index_k_neighbour;
      int index_i_procs, index_j_procs, index_k_procs;
      int index_i_neighbour_global, index_j_neighbour_global, index_k_neighbour_global;
      int m,l,n;
      double neighbours_colinearity = 0.;
      for (int i=0; i<ijk_intersections_subproblems_indices_.size_array(); i++)
        {
          if (thermal_subproblems_->get_dxyz_increment_bool(i))
            {
              thermal_subproblems_->get_subproblem_ijk_indices(index_i_problem, index_j_problem, index_k_problem, i);
              const FixedVector<int,3>& pure_neighbours_corrected_sign = thermal_subproblems_->get_pure_neighbours_corrected_sign(i);
              const std::vector<std::vector<std::vector<bool>>>& pure_neighbours_to_correct = thermal_subproblems_->get_pure_neighbours_to_correct(i);
              const std::vector<std::vector<std::vector<double>>>& pure_neighbours_corrected_distance = thermal_subproblems_->get_pure_neighbours_corrected_distance(i);
              if (neighbours_colinearity_weighting_)
                pure_neighbours_corrected_colinearity = thermal_subproblems_->get_pure_neighbours_corrected_colinearity(i);
              const int l_dir_size = (int) pure_neighbours_to_correct.size() -1;
              const int m_dir_size = (int) pure_neighbours_to_correct[0].size() -1;
              const int n_dir_size = (int) pure_neighbours_to_correct[0][0].size() -1;
              for (l=l_dir_size; l>=0; l--)
                for (m=m_dir_size; m>=0; m--)
                  for (n=n_dir_size; n>=0; n--)
                    if (pure_neighbours_to_correct[l][m][n])
                      {
                        const double dist_sub_res = pure_neighbours_corrected_distance[l][m][n];
                        const double temperature_ghost = thermal_subproblems_->get_temperature_profile_at_point(i, dist_sub_res);
                        index_i_neighbour = index_i_problem + ((pure_neighbours_corrected_sign[0]) ?  l * (-1) : l);
                        index_j_neighbour = index_j_problem + ((pure_neighbours_corrected_sign[1]) ?  m * (-1) : m);
                        index_k_neighbour = index_k_problem + ((pure_neighbours_corrected_sign[2]) ?  n * (-1) : n);
                        /*
                         * Handle both positive and negative periodicity !
                         */
                        index_i_neighbour_global = compute_periodic_index((index_i_neighbour + offset_i), ni_tot);
                        index_j_neighbour_global = compute_periodic_index((index_j_neighbour + offset_j), nj_tot);
                        index_k_neighbour_global = compute_periodic_index((index_k_neighbour + offset_k), nk_tot);
                        index_i_procs = compute_periodic_index(index_i_neighbour, ni);
                        index_j_procs = compute_periodic_index(index_j_neighbour, nj);
                        index_k_procs = compute_periodic_index(index_k_neighbour, nk);
                        if (neighbours_colinearity_weighting_)
                          neighbours_colinearity = pure_neighbours_corrected_colinearity[l][m][n];
                        if (index_i_procs == index_i_neighbour
                            && index_j_procs == index_j_neighbour
                            && index_k_procs == index_k_neighbour)
                          {
                            if (neighbours_colinearity_weighting_)
                              {
                                neighbours_weighting_colinearity(index_i_neighbour, index_j_neighbour, index_k_neighbour)
                                += neighbours_colinearity;
                                temperature_neighbours(index_i_neighbour, index_j_neighbour, index_k_neighbour)
                                += temperature_ghost * neighbours_colinearity;
                              }
                            else
                              temperature_neighbours(index_i_neighbour, index_j_neighbour, index_k_neighbour) += temperature_ghost;
                            neighbours_weighting(index_i_neighbour, index_j_neighbour, index_k_neighbour) += 1;
                          }
                        else
                          {
                            compute_temperature_cell_centre_neighbours_on_procs(temperature_ghost,
                                                                                neighbours_colinearity,
                                                                                index_i_neighbour_global,
                                                                                index_j_neighbour_global,
                                                                                index_k_neighbour_global);
                          }
                      }
            }
        }
      receive_temperature_cell_centre_neighbours_from_procs();
      combine_temperature_cell_centre_neighbours_from_procs(temperature_neighbours,
                                                            neighbours_weighting,
                                                            neighbours_weighting_colinearity,
                                                            ni,
                                                            nj,
                                                            nk,
                                                            offset_i,
                                                            offset_j,
                                                            offset_k);
      if (neighbours_colinearity_weighting_)
        neighbours_weighting_colinearity.echange_espace_virtuel(neighbours_weighting_colinearity.ghost());
      neighbours_weighting.echange_espace_virtuel(neighbours_weighting.ghost());
      temperature_neighbours.echange_espace_virtuel(temperature_neighbours.ghost());
    }
}


void Corrige_flux_FT_temperature_subresolution::compute_temperature_cell_centre_neighbours_on_procs(const double& temperature_neighbours,
                                                                                                    const double& neighbours_weighting_colinearity,
                                                                                                    const int& index_i_neighbour_global,
                                                                                                    const int& index_j_neighbour_global,
                                                                                                    const int& index_k_neighbour_global)
{
  indices_temperature_neighbours_on_procs_[0].append_array(index_i_neighbour_global);
  indices_temperature_neighbours_on_procs_[1].append_array(index_j_neighbour_global);
  indices_temperature_neighbours_on_procs_[2].append_array(index_k_neighbour_global);
  temperature_neighbours_on_procs_.append_array(temperature_neighbours);
  if (neighbours_colinearity_weighting_)
    neighbours_weighting_colinearity_on_procs_.append_array(neighbours_weighting_colinearity);
}

void Corrige_flux_FT_temperature_subresolution::receive_temperature_cell_centre_neighbours_from_procs()
{
  Cerr << "Copy temperature on every processors" << finl;
  if (copy_temperature_on_every_procs_)
    {
      const int nb_procs = Process::nproc();
      const int proc_num = Process::me();
      if (nb_procs > 1)
        {
          ArrOfInt local_indices_tmp;
          ArrOfDouble local_values_tmp;

          const int size_vector = indices_temperature_neighbours_on_procs_[0].size_array();
          int size_vector_total = size_vector;
          size_vector_total = Process::mp_sum(size_vector_total);

          ArrOfInt overall_numerotation(nb_procs);
          ArrOfInt start_indices(nb_procs);
          overall_numerotation(proc_num) = size_vector;
          mp_sum_for_each_item(overall_numerotation);
          int l;
          for (l=1; l<overall_numerotation.size_array(); l++)
            start_indices(l) = start_indices(l-1) + overall_numerotation(l-1);

          for (int c=0; c<3; c++)
            {
              local_indices_tmp = indices_temperature_neighbours_on_procs_[c];
              indices_temperature_neighbours_on_procs_[c].resize(size_vector_total);
              indices_temperature_neighbours_on_procs_[c] *= 0;
              for (l=0; l<local_indices_tmp.size_array(); l++)
                indices_temperature_neighbours_on_procs_[c](start_indices(proc_num) + l) = local_indices_tmp(l);
            }
          local_indices_tmp.reset();

          local_values_tmp = temperature_neighbours_on_procs_;
          temperature_neighbours_on_procs_.resize(size_vector_total);
          temperature_neighbours_on_procs_*= 0.;
          for (l=0; l<local_values_tmp.size_array(); l++)
            temperature_neighbours_on_procs_(start_indices(proc_num) + l) = local_values_tmp(l);
          if (neighbours_colinearity_weighting_)
            {
              local_values_tmp = neighbours_weighting_colinearity_on_procs_;
              neighbours_weighting_colinearity_on_procs_.resize(size_vector_total);
              neighbours_weighting_colinearity_on_procs_*= 0.;
              for (l=0; l<local_values_tmp.size_array(); l++)
                neighbours_weighting_colinearity_on_procs_(start_indices(proc_num) + l) = local_values_tmp(l);
            }

          for (int c=0; c<3; c++)
            {
              ArrOfInt& indices_dir = indices_temperature_neighbours_on_procs_[c];
              mp_sum_for_each_item(indices_dir);
            }
          mp_sum_for_each_item(temperature_neighbours_on_procs_);
          if (neighbours_colinearity_weighting_)
            mp_sum_for_each_item(neighbours_weighting_colinearity_on_procs_);
          assert(indices_temperature_neighbours_on_procs_[0].size_array() == temperature_neighbours_on_procs_.size_array());
        }
    }
}

void Corrige_flux_FT_temperature_subresolution::combine_temperature_cell_centre_neighbours_from_procs(IJK_Field_double& temperature_neighbours,
                                                                                                      IJK_Field_int& neighbours_weighting,
                                                                                                      IJK_Field_double& neighbours_weighting_colinearity,
                                                                                                      const int& ni,
                                                                                                      const int& nj,
                                                                                                      const int& nk,
                                                                                                      const int& offset_i,
                                                                                                      const int& offset_j,
                                                                                                      const int& offset_k)
{
  Cerr << "Combine temperature on every processors" << finl;
  if (copy_temperature_on_every_procs_)
    {
      const int size_array = indices_temperature_neighbours_on_procs_[0].size_array();
      int index_i_local, index_j_local, index_k_local;
      for (int l=0; l<size_array; l++)
        {
          index_i_local = indices_temperature_neighbours_on_procs_[0](l) - offset_i;
          index_j_local = indices_temperature_neighbours_on_procs_[1](l) - offset_j;
          index_k_local = indices_temperature_neighbours_on_procs_[2](l) - offset_k;
          if ((0 <= index_i_local && index_i_local < ni) &&
              (0 <= index_j_local && index_j_local < nj) &&
              (0 <= index_k_local && index_k_local < nk))
            {
              const double temperature_ghost = temperature_neighbours_on_procs_(l);
              neighbours_weighting(index_i_local, index_j_local, index_k_local) += 1;
              if (neighbours_colinearity_weighting_)
                {
                  const double colinearity = neighbours_weighting_colinearity_on_procs_(l);
                  neighbours_weighting_colinearity(index_i_local, index_j_local, index_k_local) += colinearity;
                  temperature_neighbours(index_i_local, index_j_local, index_k_local) += colinearity * temperature_ghost;
                }
              else
                temperature_neighbours(index_i_local, index_j_local, index_k_local) += temperature_ghost;
            }
        }
    }
}

void Corrige_flux_FT_temperature_subresolution::initialise_fixed_vectors(FixedVector<FixedVector<std::vector<ArrOfInt>,3>,2>& fixed_vectors,
                                                                         const int nb_k_layer)
{
  for (int l=0; l<2; l++)
    {
      for (int dir=0; dir<2; dir++)
        fixed_vectors[l][dir].resize(nb_k_layer);
      fixed_vectors[l][2].resize(nb_k_layer + 1);
    }
}

void Corrige_flux_FT_temperature_subresolution::initialise_fixed_vector(FixedVector<std::vector<ArrOfInt>,3>& fixed_vector,
                                                                        const int nb_k_layer)
{
  for (int dir=0; dir<2; dir++)
    fixed_vector[dir].resize(nb_k_layer);
  fixed_vector[2].resize(nb_k_layer + 1);
}

void Corrige_flux_FT_temperature_subresolution::initialise_fixed_vector_values(FixedVector<std::vector<ArrOfDouble>,3>& fixed_vector_values,
                                                                               const int nb_k_layer)
{
  for (int dir=0; dir<2; dir++)
    fixed_vector_values[dir].resize(nb_k_layer);
  fixed_vector_values[2].resize(nb_k_layer + 1);
}

void Corrige_flux_FT_temperature_subresolution::initialise_any_cell_neighbours_indices_to_correct(FixedVector<FixedVector<std::vector<ArrOfInt>,3>,2>& index_face_ij_flux_xyz_faces_sorted,
                                                                                                  const int global_indices)
{
  int nb_k_layer;
  if (global_indices)
    nb_k_layer = ref_ijk_ft_->itfce().I().get_splitting().get_grid_geometry().get_nb_elem_tot(2);
  else
    nb_k_layer = ref_ijk_ft_->itfce().I().nk();

  const int first_iter = !(ref_ijk_ft_->get_tstep());

  int dir, l;
  if (first_iter)
    initialise_fixed_vectors(index_face_ij_flux_xyz_faces_sorted, nb_k_layer);

  for (l=0; l<2; l++)
    for (dir=0; dir<3; dir++)
      for (int k_layer=0; k_layer<nb_k_layer+1; k_layer++)
        {
          if ((dir==DIRECTION_I || dir==DIRECTION_J) && k_layer==nb_k_layer)
            break;
          index_face_ij_flux_xyz_faces_sorted[l][dir][k_layer].reset();
          if (first_iter)
            index_face_ij_flux_xyz_faces_sorted[l][dir][k_layer].set_smart_resize(1);
        }
}

void Corrige_flux_FT_temperature_subresolution::initialise_any_cell_neighbours_indices_to_correct_with_flux(FixedVector<FixedVector<std::vector<ArrOfInt>,3>,2>& index_face_ij_flux_xyz_faces_sorted,
                                                                                                            FixedVector<std::vector<ArrOfDouble>,3>& fluxes,
                                                                                                            FixedVector<std::vector<ArrOfInt>,3>& weighting_flux_xyz_faces_sorted,
                                                                                                            FixedVector<std::vector<ArrOfDouble>,3>& colinearity_flux_xyz_faces_sorted,
                                                                                                            const bool& ini_index,
                                                                                                            const int global_indices,
                                                                                                            const int weighting_colinearity)
{
  if (!ini_index)
    initialise_any_cell_neighbours_indices_to_correct(index_face_ij_flux_xyz_faces_sorted,
                                                      global_indices);

  int nb_k_layer;
  if (global_indices)
    nb_k_layer = ref_ijk_ft_->itfce().I().get_splitting().get_grid_geometry().get_nb_elem_tot(2);
  else
    nb_k_layer = ref_ijk_ft_->itfce().I().nk();

  const int first_iter = !(ref_ijk_ft_->get_tstep());

  if (first_iter)
    {
      initialise_fixed_vector_values(fluxes, nb_k_layer);
      if (weighting_colinearity && !ini_index)
        {
          initialise_fixed_vector(weighting_flux_xyz_faces_sorted, nb_k_layer);
          initialise_fixed_vector_values(colinearity_flux_xyz_faces_sorted, nb_k_layer);
        }
    }

  int dir;
  for (dir=0; dir<3; dir++)
    for (int k_layer=0; k_layer<nb_k_layer+1; k_layer++)
      {
        if ((dir==DIRECTION_I || dir==DIRECTION_J) && k_layer==nb_k_layer)
          break;
        fluxes[dir][k_layer].reset();
        if (first_iter)
          fluxes[dir][k_layer].set_smart_resize(1);
        if (weighting_colinearity && !ini_index)
          {
            weighting_flux_xyz_faces_sorted[dir][k_layer].reset();
            colinearity_flux_xyz_faces_sorted[dir][k_layer].reset();
            if (first_iter)
              {
                weighting_flux_xyz_faces_sorted[dir][k_layer].set_smart_resize(1);
                colinearity_flux_xyz_faces_sorted[dir][k_layer].set_smart_resize(1);
              }
          }
      }
}

void Corrige_flux_FT_temperature_subresolution::initialise_cell_neighbours_indices_to_correct()
{
  if (distance_cell_faces_from_lrs_ && find_temperature_cell_neighbours_ && find_cell_neighbours_for_fluxes_spherical_correction_)
    {
      /*
       * Test to correct flux in the diagonal (strongly not aligned with the cartesian grid !)
       */
      initialise_any_cell_neighbours_indices_to_correct(index_face_ij_flux_xyz_neighbours_diag_faces_sorted_);
    }

  if (distance_cell_faces_from_lrs_ && find_reachable_fluxes_)
    {
      const int seq = (Process::nproc() == 1);
      if (!seq)
        {
          if (!convection_negligible_)
            {
              initialise_any_cell_neighbours_indices_to_correct_with_flux(index_face_ij_flux_xyz_neighbours_all_faces_remaining_global_sorted_,
                                                                          convective_diffusive_flux_all_faces_remaining_global_sorted_[0],
                                                                          weighting_flux_xyz_neighbours_all_faces_remaining_global_sorted_,
                                                                          colinearity_flux_xyz_neighbours_all_faces_remaining_global_sorted_,
                                                                          flux_init_,
                                                                          1,
                                                                          1);
              Cerr << "Thermal Sub-resolutions all reachable convective fluxes variables are now initialised" << finl;
              flux_init_ = 1;
            }
          if (!diffusion_negligible_)
            {
              initialise_any_cell_neighbours_indices_to_correct_with_flux(index_face_ij_flux_xyz_neighbours_all_faces_remaining_global_sorted_,
                                                                          convective_diffusive_flux_all_faces_remaining_global_sorted_[1],
                                                                          weighting_flux_xyz_neighbours_all_faces_remaining_global_sorted_,
                                                                          colinearity_flux_xyz_neighbours_all_faces_remaining_global_sorted_,
                                                                          flux_init_,
                                                                          1,
                                                                          1);
              Cerr << "Thermal Sub-resolutions all reachable diffusive fluxes variables are now initialised" << finl;
              flux_init_ = 1;
            }
          flux_init_ = 0;
        }
    }

  if (distance_cell_faces_from_lrs_ && use_reachable_fluxes_)
    {
      if (!convection_negligible_)
        {
          Cerr << "Sort the thermal min-max reachable convective fluxes" << finl;
          FixedVector<std::vector<ArrOfInt>,3>& dummy_int_array = index_face_ij_flux_xyz_neighbours_min_max_faces_sorted_[0];
          FixedVector<std::vector<ArrOfDouble>,3>& dummy_double_array = convective_diffusive_flux_xyz_min_max_faces_sorted_[0];
          initialise_any_cell_neighbours_indices_to_correct_with_flux(index_face_ij_flux_xyz_neighbours_min_max_faces_sorted_,
                                                                      convective_diffusive_flux_xyz_min_max_faces_sorted_[0],
                                                                      dummy_int_array,
                                                                      dummy_double_array,
                                                                      flux_init_);
          Cerr << "Thermal Sub-resolutions min-max reachable convective fluxes variables are now initialised" << finl;
          flux_init_ = 1;
        }
      if (!diffusion_negligible_)
        {
          Cerr << "Sort the thermal diffusive min-max reachable fluxes" << finl;
          FixedVector<std::vector<ArrOfInt>,3>& dummy_int_array = index_face_ij_flux_xyz_neighbours_min_max_faces_sorted_[0];
          FixedVector<std::vector<ArrOfDouble>,3>& dummy_double_array = convective_diffusive_flux_xyz_min_max_faces_sorted_[0];
          initialise_any_cell_neighbours_indices_to_correct_with_flux(index_face_ij_flux_xyz_neighbours_min_max_faces_sorted_,
                                                                      convective_diffusive_flux_xyz_min_max_faces_sorted_[1],
                                                                      dummy_int_array,
                                                                      dummy_double_array,
                                                                      flux_init_);
          Cerr << "Thermal Sub-resolutions min-max reachable diffusive fluxes variables are now initialised" << finl;
          flux_init_ = 1;
        }
      flux_init_ = 0;
    }
}


void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_faces_indices_for_spherical_correction(const int& n_iter_distance)
{
  /*
   * TODO: Works only in sequential
   */
  if (distance_cell_faces_from_lrs_ && find_cell_neighbours_for_fluxes_spherical_correction_)
    {
      for (int c=0; c<3; c++)
        (*cell_faces_neighbours_corrected_bool_)[c].data() = 0;
      (*cell_faces_neighbours_corrected_bool_).echange_espace_virtuel();
      const int nb_i_layer = ref_ijk_ft_->itfce().I().ni();
      const int nb_j_layer = ref_ijk_ft_->itfce().I().nj();
      const int nb_k_layer = ref_ijk_ft_->itfce().I().nk();

      // index_face_i_sorted[0] = &index_face_i_flux_x_neighbours_diag_faces_sorted_;

      int m,l,n;
      int index_i_problem, index_j_problem, index_k_problem;
      for (int i=0; i<ijk_intersections_subproblems_indices_.size_array(); i++)
        {
          if (thermal_subproblems_->get_dxyz_increment_bool(i))
            {
              thermal_subproblems_->get_subproblem_ijk_indices(index_i_problem, index_j_problem, index_k_problem, i);
              const FixedVector<int,3>& pure_neighbours_corrected_sign = thermal_subproblems_->get_pure_neighbours_corrected_sign(i);
              const std::vector<std::vector<std::vector<bool>>>& pure_neighbours_to_correct = thermal_subproblems_->get_pure_neighbours_to_correct(i);
              const int l_dir_size = (int) pure_neighbours_to_correct.size() -1;
              const int m_dir_size = (int) pure_neighbours_to_correct[0].size() -1;
              const int n_dir_size = (int) pure_neighbours_to_correct[0][0].size() -1;
              for (l=l_dir_size; l>=0; l--)
                for (m=m_dir_size; m>=0; m--)
                  for (n=n_dir_size; n>=0; n--)
                    if (l <= n_iter_distance && m <= n_iter_distance && n <= n_iter_distance) // Need an underlying normal vector for correction
                      if (pure_neighbours_to_correct[l][m][n] && ((l!=0 && m!=0) || (l!=0 && n!=0) || (m!=0 && n!=0))) // (l!=0 || m!=0 || n!=0) &&
                        {
                          const int i_offset = ((pure_neighbours_corrected_sign[0]) ?  l * (-1) : l);
                          const int j_offset = ((pure_neighbours_corrected_sign[1]) ?  m * (-1) : m);
                          const int k_offset = ((pure_neighbours_corrected_sign[2]) ?  n * (-1) : n);
                          const int index_i_neighbour = index_i_problem + i_offset;
                          const int index_j_neighbour = index_j_problem + j_offset;
                          const int index_k_neighbour = index_k_problem + k_offset;
                          const int i_offset_face = ((signbit(i_offset)) ?  (index_i_neighbour + 1) : index_i_neighbour);
                          const int j_offset_face = ((signbit(j_offset)) ?  (index_j_neighbour + 1) : index_j_neighbour);
                          const int k_offset_face = ((signbit(k_offset)) ?  (index_k_neighbour + 1) : index_k_neighbour);
                          if (l !=0 )
                            {
                              if (!(*cell_faces_neighbours_corrected_bool_)[0](i_offset_face, index_j_neighbour, index_k_neighbour))
                                {
                                  index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[0][0][index_k_neighbour].append_array(i_offset_face);
                                  index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[1][0][index_k_neighbour].append_array(index_j_neighbour);
                                }
                              (*cell_faces_neighbours_corrected_bool_)[0](i_offset_face, index_j_neighbour, index_k_neighbour) += 1;
                              if (i_offset_face == 0)
                                {
                                  if(!(*cell_faces_neighbours_corrected_bool_)[0](nb_i_layer, index_j_neighbour, index_k_neighbour))
                                    {
                                      index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[0][0][index_k_neighbour].append_array(nb_i_layer);
                                      index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[1][0][index_k_neighbour].append_array(index_j_neighbour);
                                    }
                                  (*cell_faces_neighbours_corrected_bool_)[0](nb_i_layer, index_j_neighbour, index_k_neighbour) += 1;
                                }
                              if (i_offset_face == nb_i_layer)
                                {
                                  if (!(*cell_faces_neighbours_corrected_bool_)[0](0, index_j_neighbour, index_k_neighbour))
                                    {
                                      index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[0][0][index_k_neighbour].append_array(0);
                                      index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[1][0][index_k_neighbour].append_array(index_j_neighbour);
                                    }
                                  (*cell_faces_neighbours_corrected_bool_)[0](0, index_j_neighbour, index_k_neighbour) += 1;
                                }
                            }
                          if (m != 0)
                            {
                              if (!(*cell_faces_neighbours_corrected_bool_)[1](index_i_neighbour, j_offset_face, index_k_neighbour))
                                {
                                  index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[0][1][index_k_neighbour].append_array(index_i_neighbour);
                                  index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[1][1][index_k_neighbour].append_array(j_offset_face);
                                }
                              (*cell_faces_neighbours_corrected_bool_)[1](index_i_neighbour, j_offset_face, index_k_neighbour) += 1;
                              if (j_offset_face == 0)
                                {
                                  if (!(*cell_faces_neighbours_corrected_bool_)[1](index_i_neighbour, nb_j_layer, index_k_neighbour))
                                    {
                                      index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[0][1][index_k_neighbour].append_array(index_i_neighbour);
                                      index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[1][1][index_k_neighbour].append_array(nb_j_layer);
                                    }
                                  (*cell_faces_neighbours_corrected_bool_)[1](index_i_neighbour, nb_j_layer, index_k_neighbour) += 1;
                                }
                              if (i_offset_face == nb_j_layer)
                                {
                                  if (!(*cell_faces_neighbours_corrected_bool_)[1](index_i_neighbour, 0, index_k_neighbour))
                                    {
                                      index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[0][1][index_k_neighbour].append_array(index_i_neighbour);
                                      index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[1][1][index_k_neighbour].append_array(0);
                                    }
                                  (*cell_faces_neighbours_corrected_bool_)[1](index_i_neighbour, 0, index_k_neighbour) += 1;
                                }
                            }
                          if (n != 0)
                            {
                              if (!(*cell_faces_neighbours_corrected_bool_)[2](index_i_neighbour, index_j_neighbour, k_offset_face))
                                {
                                  index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[0][2][k_offset_face].append_array(index_i_neighbour);
                                  index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[1][2][k_offset_face].append_array(index_j_neighbour);
                                }
                              (*cell_faces_neighbours_corrected_bool_)[2](index_i_neighbour, index_j_neighbour, k_offset_face) += 1;
                              if (j_offset_face == 0)
                                {
                                  if (!(*cell_faces_neighbours_corrected_bool_)[2](index_i_neighbour, index_j_neighbour, nb_k_layer))
                                    {
                                      index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[0][2][nb_k_layer].append_array(index_i_neighbour);
                                      index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[1][2][nb_k_layer].append_array(index_j_neighbour);
                                    }
                                  (*cell_faces_neighbours_corrected_bool_)[2](index_i_neighbour, index_j_neighbour, nb_k_layer) += 1;
                                }
                              if (i_offset_face == nb_k_layer)
                                {
                                  if (!(*cell_faces_neighbours_corrected_bool_)[2](index_i_neighbour, index_j_neighbour, 0))
                                    {
                                      index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[0][2][0].append_array(index_i_neighbour);
                                      index_face_ij_flux_xyz_neighbours_diag_faces_sorted_[1][2][0].append_array(index_j_neighbour);
                                    }
                                  (*cell_faces_neighbours_corrected_bool_)[2](index_i_neighbour, index_j_neighbour, 0) += 1;
                                }
                            }
                        }
            }
        }
      (*cell_faces_neighbours_corrected_bool_).echange_espace_virtuel();
    }
}

void Corrige_flux_FT_temperature_subresolution::correct_flux_spherical(Simd_double& a,
                                                                       Simd_double& b,
                                                                       const int& i,
                                                                       const int& j,
                                                                       const int& k_layer,
                                                                       const int dir)
{
  /*
   * I wanted to use index_face_i_sorted, index_face_j_sorted to avoid reading cell_faces_neighbours_corrected_bool again
   */
  if (use_cell_neighbours_for_fluxes_spherical_correction_)
    {
      int bool_a, bool_b;
      bool_a = (*cell_faces_neighbours_corrected_bool_)[dir](i, j, k_layer);
      bool_a = bool_a || (*cell_faces_neighbours_corrected_bool_)[dir](i, j, k_layer);
      bool_a = bool_a || (*cell_faces_neighbours_corrected_bool_)[dir](i, j, k_layer);
      bool_b = (*cell_faces_neighbours_corrected_bool_)[dir](i+1, j, k_layer);
      bool_b = bool_b || (*cell_faces_neighbours_corrected_bool_)[dir](i, j+1, k_layer);
      bool_b = bool_b || (*cell_faces_neighbours_corrected_bool_)[dir](i, j, k_layer + 1);
      if (bool_a || bool_b)
        {
          double n;
          n = abs((*eulerian_normal_vectors_ns_normed_)[dir](i, j, k_layer));
          a *= n;
          b *= n;
        }
    }
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_mixed_cell_faces_indices_to_correct(FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_bool_mixed_cell,
                                                                                                            FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_convective_mixed_cell,
                                                                                                            FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_diffusive_mixed_cell,
                                                                                                            FixedVector<IJK_Field_double, 3>& neighbours_weighting_colinearity_mixed_cell)
{

  if (distance_cell_faces_from_lrs_ && find_reachable_fluxes_ && keep_first_reachable_fluxes_)
    {

      const int ni = ref_ijk_ft_->itfce().I().ni();
      const int nj = ref_ijk_ft_->itfce().I().nj();
      const int nk = ref_ijk_ft_->itfce().I().nk();

      IJK_Field_local_int cell_faces_neighbours_corrected_bool_tmp;
      cell_faces_neighbours_corrected_bool_tmp.allocate(ni, nj, nk, 1);

      const int neighbours_i[6] = NEIGHBOURS_I;
      const int neighbours_j[6] = NEIGHBOURS_J;
      const int neighbours_k[6] = NEIGHBOURS_K;
      const int neighbours_faces_i[6] = NEIGHBOURS_FACES_I;
      const int neighbours_faces_j[6] = NEIGHBOURS_FACES_J;
      const int neighbours_faces_k[6] = NEIGHBOURS_FACES_K;

      for (int c=0; c<3; c++)
        {
          const int index_ini = 2*c;
          cell_faces_neighbours_corrected_bool_tmp.data() = 0;
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  const double indic = ref_ijk_ft_->itfce().I()(i,j,k);
                  if (indic < LIQUID_INDICATOR_TEST && indic > VAPOUR_INDICATOR_TEST)
                    {
                      for (int l=index_ini; l<index_ini + 2; l++)
                        {
                          const int i_neighbour = neighbours_i[l];
                          const int j_neighbour = neighbours_j[l];
                          const int k_neighbour = neighbours_k[l];
                          const int ii = neighbours_faces_i[l];
                          const int jj = neighbours_faces_j[l];
                          const int kk = neighbours_faces_k[l];
                          const int cell_faces_neighbours_ijk = cell_faces_neighbours_corrected_bool_mixed_cell[c](i + ii,j + jj, k + kk);
                          const double indic_neighbour = ref_ijk_ft_->itfce().I()(i+i_neighbour,j+j_neighbour,k+k_neighbour);
                          if (cell_faces_neighbours_ijk && indic_neighbour > LIQUID_INDICATOR_TEST)
                            cell_faces_neighbours_corrected_bool_tmp(i+ii,j+jj,k+kk) = cell_faces_neighbours_ijk;
                        }
                    }
                  if (indic > LIQUID_INDICATOR_TEST)
                    {
                      for (int l=index_ini; l<index_ini + 2; l++)
                        {
                          const int i_neighbour = neighbours_i[l];
                          const int j_neighbour = neighbours_j[l];
                          const int k_neighbour = neighbours_k[l];
                          const int ii = neighbours_faces_i[l];
                          const int jj = neighbours_faces_j[l];
                          const int kk = neighbours_faces_k[l];
                          const int cell_faces_neighbours_ijk = cell_faces_neighbours_corrected_bool_mixed_cell[c](i + ii,j + jj, k + kk);
                          const double indic_neighbour = ref_ijk_ft_->itfce().I()(i+i_neighbour,j+j_neighbour,k+k_neighbour);
                          if (cell_faces_neighbours_ijk && (indic_neighbour < LIQUID_INDICATOR_TEST && indic_neighbour > VAPOUR_INDICATOR_TEST))
                            cell_faces_neighbours_corrected_bool_tmp(i+ii,j+jj,k+kk) = cell_faces_neighbours_ijk;
                        }
                    }
                }
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                cell_faces_neighbours_corrected_bool_mixed_cell[c](i,j,k) = cell_faces_neighbours_corrected_bool_tmp(i,j,k);
        }
      // cell_faces_neighbours_corrected_bool_tmp.~IJK_Field_local_int();
      if (use_reachable_fluxes_)
        {
          IJK_Field_local_double cell_faces_neighbours_corrected_field_tmp;
          if (convective_flux_correction_ || diffusive_flux_correction_ || neighbours_colinearity_weighting_)
            cell_faces_neighbours_corrected_field_tmp.allocate(ni, nj, nk, 0);

          if (convective_flux_correction_)
            compute_cell_neighbours_mixed_cell_faces_any_field(cell_faces_neighbours_corrected_bool_mixed_cell,
                                                               cell_faces_neighbours_corrected_field_tmp,
                                                               cell_faces_neighbours_corrected_convective_mixed_cell);


          if (diffusive_flux_correction_)
            compute_cell_neighbours_mixed_cell_faces_any_field(cell_faces_neighbours_corrected_bool_mixed_cell,
                                                               cell_faces_neighbours_corrected_field_tmp,
                                                               cell_faces_neighbours_corrected_diffusive_mixed_cell);

          if (neighbours_colinearity_weighting_)
            compute_cell_neighbours_mixed_cell_faces_any_field(cell_faces_neighbours_corrected_bool_mixed_cell,
                                                               cell_faces_neighbours_corrected_field_tmp,
                                                               neighbours_weighting_colinearity_mixed_cell);
        }
    }
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_mixed_cell_faces_any_field(FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_bool,
                                                                                                   IJK_Field_local_double& cell_faces_neighbours_corrected_field,
                                                                                                   FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_field_mixed_cell)
{
  const int ni = cell_faces_neighbours_corrected_field_mixed_cell[0].ni();
  const int nj = cell_faces_neighbours_corrected_field_mixed_cell[1].nj();
  const int nk = cell_faces_neighbours_corrected_field_mixed_cell[2].nk();

  for (int c=0; c<3; c++)
    {
      cell_faces_neighbours_corrected_field.data() = 0;
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            if (cell_faces_neighbours_corrected_bool[c](i,j,k))
              cell_faces_neighbours_corrected_field(i,j,k) = cell_faces_neighbours_corrected_field_mixed_cell[c](i,j,k);
      cell_faces_neighbours_corrected_field_mixed_cell[c].data() = 0.;
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            cell_faces_neighbours_corrected_field_mixed_cell[c](i,j,k) = cell_faces_neighbours_corrected_field(i,j,k);

    }
  cell_faces_neighbours_corrected_field_mixed_cell.echange_espace_virtuel();
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_faces_indices_to_correct(FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_bool,
                                                                                                 FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_convective,
                                                                                                 FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_diffusive,
                                                                                                 FixedVector<IJK_Field_double, 3>& neighbours_weighting_colinearity)
{
  if (distance_cell_faces_from_lrs_ && find_reachable_fluxes_)
    {
      for (int c=0; c<3; c++)
        cell_faces_neighbours_corrected_bool[c].data() = 0;
      cell_faces_neighbours_corrected_bool.echange_espace_virtuel();

      if (neighbours_colinearity_weighting_ && use_reachable_fluxes_)
        {
          for (int c=0; c<3; c++)
            neighbours_weighting_colinearity[c].data() = 0;
          neighbours_weighting_colinearity.echange_espace_virtuel();
        }

      if (convective_flux_correction_ && use_reachable_fluxes_)
        {
          for (int c=0; c<3; c++)
            cell_faces_neighbours_corrected_convective[c].data() = 0;
          cell_faces_neighbours_corrected_convective.echange_espace_virtuel();
        }

      if (diffusive_flux_correction_ && use_reachable_fluxes_)
        {
          for (int c=0; c<3; c++)
            cell_faces_neighbours_corrected_diffusive[c].data() = 0;
          cell_faces_neighbours_corrected_diffusive.echange_espace_virtuel();
        }

      const int seq = (Process::nproc() == 1);

      const int nb_i_layer = cell_faces_neighbours_corrected_bool[0].ni();
      const int nb_j_layer = cell_faces_neighbours_corrected_bool[0].nj();
      const int nb_k_layer = cell_faces_neighbours_corrected_bool[0].nk();

      const IJK_Splitting& splitting_ns = ref_ijk_ft_->itfce().I().get_splitting();
      const int offset_i = splitting_ns.get_offset_local(0);
      const int offset_j = splitting_ns.get_offset_local(1);
      const int offset_k = splitting_ns.get_offset_local(2);

      const IJK_Grid_Geometry& geometry = splitting_ns.get_grid_geometry();
      const int nb_i_layer_tot = geometry.get_nb_elem_tot(0);
      const int nb_j_layer_tot = geometry.get_nb_elem_tot(1);
      const int nb_k_layer_tot = geometry.get_nb_elem_tot(2);

      int m,l,n;
      int index_i_problem, index_j_problem, index_k_problem;
      int index_i_neighbour, index_j_neighbour, index_k_neighbour;
      int index_i_procs, index_j_procs, index_k_procs;
      int index_i_neighbour_global, index_j_neighbour_global, index_k_neighbour_global;
      int l_dir, m_dir, n_dir;
      double distance, colinearity, convective_flux, diffusive_flux;
      for (int i=0; i<ijk_intersections_subproblems_indices_.size_array(); i++)
        {
          if (thermal_subproblems_->get_dxyz_over_two_increment_bool(i))
            {
              thermal_subproblems_->get_subproblem_ijk_indices(index_i_problem, index_j_problem, index_k_problem, i);
              const FixedVector<int,3>& pure_neighbours_corrected_sign = thermal_subproblems_->get_pure_neighbours_corrected_sign(i);
              const std::vector<std::vector<std::vector<std::vector<bool>>>>& pure_neighbours_to_correct = thermal_subproblems_->get_pure_neighbours_last_faces_to_correct(i);
              const std::vector<std::vector<std::vector<std::vector<double>>>>& pure_neighbours_distance_to_correct = thermal_subproblems_->get_pure_neighbours_last_faces_corrected_distance(i);
              const std::vector<std::vector<std::vector<std::vector<double>>>>& pure_neighbours_colinearity_to_correct = thermal_subproblems_->get_pure_neighbours_last_faces_corrected_colinearity(i);
              for (int c=0; c<3; c++)
                {
                  const int l_dir_size = (int) pure_neighbours_to_correct[c].size() -1;
                  const int m_dir_size = (int) pure_neighbours_to_correct[c][0].size() -1;
                  const int n_dir_size = (int) pure_neighbours_to_correct[c][0][0].size() -1;
                  for (l=l_dir_size; l>=0; l--)
                    for (m=m_dir_size; m>=0; m--)
                      for (n=n_dir_size; n>=0; n--)
                        if (pure_neighbours_to_correct[c][l][m][n])
                          {
                            switch(c)
                              {
                              case 0:
                                // l_dir = (pure_neighbours_corrected_sign[0]) ? l * (-1) : l + 1;
                                l_dir = (pure_neighbours_corrected_sign[0]) ? l * (-1) + 1 : l;
                                m_dir = (pure_neighbours_corrected_sign[1]) ? m * (-1) : m;
                                n_dir = (pure_neighbours_corrected_sign[2]) ? n * (-1) : n;
                                break;
                              case 1:
                                l_dir = (pure_neighbours_corrected_sign[0]) ? l * (-1) : l;
                                m_dir = (pure_neighbours_corrected_sign[1]) ? m * (-1) + 1 : m;
                                n_dir = (pure_neighbours_corrected_sign[2]) ? n * (-1) : n;
                                break;
                              case 2:
                                l_dir = (pure_neighbours_corrected_sign[0]) ? l * (-1) : l;
                                m_dir = (pure_neighbours_corrected_sign[1]) ? m * (-1) : m;
                                n_dir = (pure_neighbours_corrected_sign[2]) ? n * (-1) + 1 : n;
                                break;
                              default:
                                l_dir = (pure_neighbours_corrected_sign[0]) ? l * (-1) + 1 : l;
                                m_dir = (pure_neighbours_corrected_sign[1]) ? m * (-1) : m;
                                n_dir = (pure_neighbours_corrected_sign[2]) ? n * (-1) : n;
                                break;
                              }
                            /*
                             * TODO: Handle the periodicity and check if it works
                             */
                            index_i_neighbour = (index_i_problem + l_dir);
                            index_j_neighbour = (index_j_problem + m_dir);
                            index_k_neighbour = (index_k_problem + n_dir);
                            index_i_neighbour_global = compute_periodic_index((index_i_neighbour + offset_i), nb_i_layer_tot);
                            index_j_neighbour_global = compute_periodic_index((index_j_neighbour + offset_j), nb_j_layer_tot);
                            index_k_neighbour_global = compute_periodic_index((index_k_neighbour + offset_k), nb_k_layer_tot);
                            index_i_procs = compute_periodic_index(index_i_neighbour, nb_i_layer);
                            index_j_procs = compute_periodic_index(index_j_neighbour, nb_j_layer);
                            index_k_procs = compute_periodic_index(index_k_neighbour, nb_k_layer);
                            colinearity = 1.;
                            distance = pure_neighbours_distance_to_correct[c][l][m][n];
                            if (neighbours_colinearity_weighting_ && use_reachable_fluxes_)
                              colinearity = pure_neighbours_colinearity_to_correct[c][l][m][n];
                            if ((index_i_procs == index_i_neighbour
                                 && index_j_procs == index_j_neighbour
                                 && index_k_procs == index_k_neighbour)
                                || seq)
                              {
                                cell_faces_neighbours_corrected_bool[c](index_i_procs, index_j_neighbour, index_k_procs) += 1;
                                if (neighbours_colinearity_weighting_ && use_reachable_fluxes_)
                                  neighbours_weighting_colinearity[c](index_i_procs, index_j_procs, index_k_procs) += colinearity;

                                compute_cell_neighbours_fluxes_to_correct(cell_faces_neighbours_corrected_convective,
                                                                          cell_faces_neighbours_corrected_diffusive,
                                                                          i,
                                                                          index_i_procs, index_j_procs, index_k_procs,
                                                                          distance,
                                                                          c,
                                                                          colinearity,
                                                                          use_reachable_fluxes_,
                                                                          convective_flux,
                                                                          diffusive_flux);
                              }
                            else
                              {
                                compute_flux_neighbours_on_procs(index_i_neighbour_global,
                                                                 index_j_neighbour_global,
                                                                 index_k_neighbour_global,
                                                                 i,
                                                                 distance,
                                                                 c,
                                                                 colinearity,
                                                                 use_reachable_fluxes_);
                              }
                          }
                }
            }
        }
      receive_all_fluxes_from_outisde_frontier_on_procs();
      if (debug_)
        Cerr << "Fluxes have been received from procs" << finl;

      combine_all_fluxes_from_outisde_frontier_on_procs(cell_faces_neighbours_corrected_bool,
                                                        cell_faces_neighbours_corrected_convective,
                                                        cell_faces_neighbours_corrected_diffusive,
                                                        neighbours_weighting_colinearity);
      if (debug_)
        Cerr << "Fluxes have been combined on procs" << finl;


      complete_neighbours_and_weighting_colinearity(cell_faces_neighbours_corrected_bool,
                                                    cell_faces_neighbours_corrected_convective,
                                                    cell_faces_neighbours_corrected_diffusive,
                                                    neighbours_weighting_colinearity,
                                                    use_reachable_fluxes_);
      if (debug_)
        Cerr << "Weighted calculation of all fluxes has been performed" << finl;


      compute_cell_neighbours_mixed_cell_faces_indices_to_correct(cell_faces_neighbours_corrected_bool,
                                                                  cell_faces_neighbours_corrected_convective,
                                                                  cell_faces_neighbours_corrected_diffusive,
                                                                  neighbours_weighting_colinearity);
      if (debug_)
        Cerr << "Only the fluxes in the immediate interface vicinity have been eventually kept" << finl;
    }
}


void Corrige_flux_FT_temperature_subresolution::compute_flux_neighbours_on_procs(const int& index_i_neighbour_global,
                                                                                 const int& index_j_neighbour_global,
                                                                                 const int& index_k_neighbour_global,
                                                                                 const int& subproblem_index,
                                                                                 const double& dist,
                                                                                 const int& dir,
                                                                                 const double& colinearity,
                                                                                 const double& convective_flux_computed,
                                                                                 const double& diffusive_flux_computed)
{
  double convective_flux = convective_flux_computed;
  double diffusive_flux = diffusive_flux_computed;
  if (use_reachable_fluxes_)
    {
      if (convective_flux_correction_)
        compute_cell_neighbours_convective_fluxes_to_correct(convective_flux, subproblem_index, dist, dir, colinearity);

      if (diffusive_flux_correction_)
        compute_cell_neighbours_diffusive_fluxes_to_correct(diffusive_flux, subproblem_index, dist, dir, colinearity);
    }

  const int linear_index_global = get_linear_index_global(index_i_neighbour_global, index_j_neighbour_global, index_k_neighbour_global, dir);
  const int count_val = (int) flux_outside_frontier_all_map_[dir].count(linear_index_global);

  if (!count_val)
    {
      const int size_array = (int) index_face_ij_flux_xyz_neighbours_all_faces_remaining_global_sorted_[0][dir][index_k_neighbour_global].size_array();
      flux_outside_frontier_all_map_[dir][linear_index_global] = size_array;
      index_face_ij_flux_xyz_neighbours_all_faces_remaining_global_sorted_[0][dir][index_k_neighbour_global].append_array(index_i_neighbour_global);
      index_face_ij_flux_xyz_neighbours_all_faces_remaining_global_sorted_[1][dir][index_k_neighbour_global].append_array(index_j_neighbour_global);
      weighting_flux_xyz_neighbours_all_faces_remaining_global_sorted_[dir][index_k_neighbour_global].append_array(1);
      if (use_reachable_fluxes_)
        {
          if(neighbours_colinearity_weighting_)
            colinearity_flux_xyz_neighbours_all_faces_remaining_global_sorted_[dir][index_k_neighbour_global].append_array(colinearity);
          if (convective_flux_correction_)
            convective_diffusive_flux_all_faces_remaining_global_sorted_[0][dir][index_k_neighbour_global].append_array(convective_flux);
          if (diffusive_flux_correction_)
            convective_diffusive_flux_all_faces_remaining_global_sorted_[1][dir][index_k_neighbour_global].append_array(diffusive_flux);
        }
    }
  else
    {
      const int index_array = flux_outside_frontier_all_map_[dir][linear_index_global];
      weighting_flux_xyz_neighbours_all_faces_remaining_global_sorted_[dir][index_k_neighbour_global](index_array) += 1;
      if (use_reachable_fluxes_)
        {
          if(neighbours_colinearity_weighting_)
            colinearity_flux_xyz_neighbours_all_faces_remaining_global_sorted_[dir][index_k_neighbour_global](index_array) += colinearity;
          if (convective_flux_correction_)
            convective_diffusive_flux_all_faces_remaining_global_sorted_[0][dir][index_k_neighbour_global](index_array) += convective_flux;
          if (diffusive_flux_correction_)
            convective_diffusive_flux_all_faces_remaining_global_sorted_[1][dir][index_k_neighbour_global](index_array) += diffusive_flux;
        }
    }
}


void Corrige_flux_FT_temperature_subresolution::receive_all_fluxes_from_outisde_frontier_on_procs()
{
  const int nb_procs = Process::nproc();
  const int proc_num = Process::me();
  if (copy_fluxes_on_every_procs_)
    {
      if (nb_procs > 1)
        {
          for (int c=0; c<3; c++)
            {
              const int size_k_layers = (int) index_face_ij_flux_xyz_neighbours_all_faces_remaining_global_sorted_[0][c].size();
              for (int k=0; k<size_k_layers; k++)
                {
                  const int size_array = index_face_ij_flux_xyz_neighbours_all_faces_remaining_global_sorted_[0][c][k].size_array();
                  int size_array_global = size_array;
                  size_array_global = mp_sum(size_array_global);
                  ArrOfInt overall_numerotation(nb_procs);
                  ArrOfInt start_indices(nb_procs);
                  overall_numerotation(proc_num) = size_array;
                  mp_sum_for_each_item(overall_numerotation);
                  int l;
                  for (l=1; l<overall_numerotation.size_array(); l++)
                    start_indices(l) = start_indices(l-1) + overall_numerotation(l-1);

                  Cerr << "Size array" << size_array << finl;
                  Cerr << "Size array global" << size_array_global << finl;
                  Cerr << "Overall_numerotation" << overall_numerotation(0) << "-" << overall_numerotation(1) << finl;

                  ArrOfInt local_indices_i_tmp;
                  ArrOfInt local_indices_j_tmp;
                  ArrOfInt local_weighting_tmp;

                  ArrOfInt& global_indices_i_tmp = index_face_ij_flux_xyz_neighbours_all_faces_remaining_global_sorted_[0][c][k];
                  ArrOfInt& global_indices_j_tmp = index_face_ij_flux_xyz_neighbours_all_faces_remaining_global_sorted_[1][c][k];
                  ArrOfInt& global_weighting_tmp = weighting_flux_xyz_neighbours_all_faces_remaining_global_sorted_[c][k];

                  local_indices_i_tmp = global_indices_i_tmp;
                  local_indices_j_tmp = global_indices_j_tmp;
                  local_weighting_tmp = global_weighting_tmp;

                  global_indices_i_tmp.resize(size_array_global);
                  global_indices_j_tmp.resize(size_array_global);
                  global_weighting_tmp.resize(size_array_global);

                  global_indices_i_tmp *= 0;
                  global_indices_j_tmp *= 0;
                  global_weighting_tmp *= 0;

                  for (l=0; l<local_indices_i_tmp.size_array(); l++)
                    {
                      global_indices_i_tmp(start_indices(proc_num) + l) = local_indices_i_tmp(l);
                      global_indices_j_tmp(start_indices(proc_num) + l) = local_indices_j_tmp(l);
                      global_weighting_tmp(start_indices(proc_num) + l) = local_weighting_tmp(l);
                    }

                  mp_sum_for_each_item(global_indices_i_tmp);
                  mp_sum_for_each_item(global_indices_j_tmp);
                  mp_sum_for_each_item(global_weighting_tmp);

                  if (use_reachable_fluxes_)
                    {
                      if(neighbours_colinearity_weighting_)
                        {
                          ArrOfDouble local_colinearity_tmp;
                          ArrOfDouble& global_colinearity_tmp = colinearity_flux_xyz_neighbours_all_faces_remaining_global_sorted_[c][k];
                          local_colinearity_tmp = global_colinearity_tmp;
                          global_colinearity_tmp.resize(size_array_global);
                          global_colinearity_tmp *= 0.;
                          for (l=0; l<local_indices_i_tmp.size_array(); l++)
                            global_colinearity_tmp(start_indices(proc_num) + l) = local_colinearity_tmp(l);
                          mp_sum_for_each_item(global_colinearity_tmp);
                        }


                      if (convective_flux_correction_)
                        {
                          ArrOfDouble local_convective_flux_values_tmp;
                          ArrOfDouble& global_convective_flux_values_tmp = convective_diffusive_flux_all_faces_remaining_global_sorted_[0][c][k];
                          local_convective_flux_values_tmp = global_convective_flux_values_tmp;
                          global_convective_flux_values_tmp.resize(size_array_global);
                          global_convective_flux_values_tmp *= 0.;
                          for (l=0; l<local_indices_i_tmp.size_array(); l++)
                            global_convective_flux_values_tmp(start_indices(proc_num) + l) = local_convective_flux_values_tmp(l);
                          mp_sum_for_each_item(global_convective_flux_values_tmp);
                        }

                      if (diffusive_flux_correction_)
                        {
                          ArrOfDouble local_diffusive_flux_values_tmp;
                          ArrOfDouble& global_diffusive_flux_values_tmp = convective_diffusive_flux_all_faces_remaining_global_sorted_[1][c][k];
                          local_diffusive_flux_values_tmp = global_diffusive_flux_values_tmp;
                          global_diffusive_flux_values_tmp.resize(size_array_global);
                          global_diffusive_flux_values_tmp *= 0.;
                          for (l=0; l<local_indices_i_tmp.size_array(); l++)
                            global_diffusive_flux_values_tmp(start_indices(proc_num) + l) = local_diffusive_flux_values_tmp(l);
                          mp_sum_for_each_item(global_diffusive_flux_values_tmp);
                        }
                    }
                }
            }
        }
    }
}

void Corrige_flux_FT_temperature_subresolution::combine_all_fluxes_from_outisde_frontier_on_procs(FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_bool,
                                                                                                  FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_convective,
                                                                                                  FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_diffusive,
                                                                                                  FixedVector<IJK_Field_double, 3>& neighbours_weighting_colinearity)
{
  const IJK_Field_double& indicator = ref_ijk_ft_->itfce().I();
  const int ni = indicator.ni();
  const int nj = indicator.nj();
  const int nk = indicator.nk();

  const IJK_Splitting& splitting_ns = ref_ijk_ft_->itfce().I().get_splitting();
  const int offset_i = splitting_ns.get_offset_local(0);
  const int offset_j = splitting_ns.get_offset_local(1);
  const int offset_k = splitting_ns.get_offset_local(2);

  double convective_flux = 0.;
  double diffusive_flux = 0.;
  double colinearity = 0.;
  for (int dir=0; dir<3; dir++)
    {
      const int size_k_layers = (int) index_face_ij_flux_xyz_neighbours_all_faces_remaining_global_sorted_[0][dir].size();
      for (int k_global=0; k_global<size_k_layers; k_global++)
        {
          const int global_size_array = index_face_ij_flux_xyz_neighbours_all_faces_remaining_global_sorted_[0][dir][k_global].size_array();
          for (int l=0; l<global_size_array; l++)
            {
              const int i_global = index_face_ij_flux_xyz_neighbours_all_faces_remaining_global_sorted_[0][dir][k_global](l);
              const int j_global = index_face_ij_flux_xyz_neighbours_all_faces_remaining_global_sorted_[1][dir][k_global](l);
              if (use_reachable_fluxes_)
                {
                  if (convective_flux_correction_)
                    convective_flux = convective_diffusive_flux_all_faces_remaining_global_sorted_[0][dir][k_global](l);
                  if (diffusive_flux_correction_)
                    diffusive_flux = convective_diffusive_flux_all_faces_remaining_global_sorted_[1][dir][k_global](l);
                  if (neighbours_colinearity_weighting_)
                    colinearity = colinearity_flux_xyz_neighbours_all_faces_remaining_global_sorted_[dir][k_global](l);
                }
              const int local_weighting = weighting_flux_xyz_neighbours_all_faces_remaining_global_sorted_[dir][k_global](l);
              const int i = i_global - offset_i;
              const int j = j_global - offset_j;
              const int k = k_global - offset_k;
              if ((0 <= i && i < ni) && (0 <= j && j < nj) && (0 <= k && k < nk))
                {
                  if (use_reachable_fluxes_)
                    {
                      if (convective_flux_correction_)
                        cell_faces_neighbours_corrected_convective[dir](i,j,k) += convective_flux;
                      if (diffusive_flux_correction_)
                        cell_faces_neighbours_corrected_diffusive[dir](i,j,k) += diffusive_flux;
                      if (neighbours_colinearity_weighting_)
                        neighbours_weighting_colinearity[dir](i,j,k) += colinearity;
                    }
                  cell_faces_neighbours_corrected_bool[dir](i,j,k) += local_weighting;
                }
            }
        }
    }
}

void Corrige_flux_FT_temperature_subresolution::complete_neighbours_and_weighting_colinearity(FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_bool,
                                                                                              FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_convective,
                                                                                              FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_diffusive,
                                                                                              FixedVector<IJK_Field_double, 3>& neighbours_weighting_colinearity,
                                                                                              const int& compute_fluxes_values)
{
  cell_faces_neighbours_corrected_bool.echange_espace_virtuel();
  if (compute_fluxes_values)
    {
      //      if (convective_flux_correction_)
      //        cell_faces_neighbours_corrected_convective.echange_espace_virtuel();
      //      if (diffusive_flux_correction_)
      //        cell_faces_neighbours_corrected_diffusive.echange_espace_virtuel();
      //			if (neighbours_colinearity_weighting_)
      //				neighbours_weighting_colinearity.echange_espace_virtuel();
      const int ni = cell_faces_neighbours_corrected_bool[0].ni();
      const int nj = cell_faces_neighbours_corrected_bool[0].nj();
      const int nk = cell_faces_neighbours_corrected_bool[0].nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            for (int c=0; c<3; c++)
              if(cell_faces_neighbours_corrected_bool[c](i,j,k))
                {
                  if (neighbours_colinearity_weighting_)
                    {
                      const double colinearity = neighbours_weighting_colinearity[c](i,j,k);
                      if (colinearity > 1e-12)
                        {
                          if (convective_flux_correction_)
                            cell_faces_neighbours_corrected_convective[c](i,j,k) /= colinearity;
                          if (diffusive_flux_correction_)
                            cell_faces_neighbours_corrected_diffusive[c](i,j,k) /= colinearity;
                        }
                    }
                  else
                    {
                      const double weighting = (double) (cell_faces_neighbours_corrected_bool[c](i,j,k));
                      if (convective_flux_correction_)
                        cell_faces_neighbours_corrected_convective[c](i,j,k) /= weighting;
                      if (diffusive_flux_correction_)
                        cell_faces_neighbours_corrected_diffusive[c](i,j,k) /= weighting;
                    }
                }
      if (convective_flux_correction_)
        cell_faces_neighbours_corrected_convective.echange_espace_virtuel();
      if (diffusive_flux_correction_)
        cell_faces_neighbours_corrected_diffusive.echange_espace_virtuel();
      if (neighbours_colinearity_weighting_)
        neighbours_weighting_colinearity.echange_espace_virtuel();
    }
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_fluxes_to_correct(FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_convective,
                                                                                          FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_diffusive,
                                                                                          const int& subproblem_index,
                                                                                          const int& index_i, const int& index_j, const int& index_k,
                                                                                          const double& dist,
                                                                                          const int& dir,
                                                                                          const double& colinearity,
                                                                                          const int& compute_fluxes_values,
                                                                                          double& convective_flux,
                                                                                          double& diffusive_flux)
{
  if (compute_fluxes_values)
    {
      if (convective_flux_correction_)
        {
//        compute_cell_neighbours_convective_fluxes_to_correct(cell_faces_neighbours_corrected_convective);
          compute_cell_neighbours_convective_fluxes_to_correct(convective_flux, subproblem_index, dist, dir, colinearity);
          cell_faces_neighbours_corrected_convective[dir](index_i, index_j, index_k) += convective_flux;
        }
      if (diffusive_flux_correction_)
        {
          compute_cell_neighbours_diffusive_fluxes_to_correct(diffusive_flux, subproblem_index, dist, dir, colinearity);
          cell_faces_neighbours_corrected_diffusive[dir](index_i, index_j, index_k) += diffusive_flux;
        }
    }
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_convective_fluxes_to_correct(double& convective_flux,
                                                                                                     const int& subproblem_index,
                                                                                                     const double& dist,
                                                                                                     const int& dir,
                                                                                                     const double& colinearity)
{
  if (!discrete_integral_)
    compute_cell_neighbours_thermal_convective_fluxes_face_centre(convective_flux,
                                                                  subproblem_index,
                                                                  dist,
                                                                  dir,
                                                                  colinearity);
  else
    compute_cell_neighbours_thermal_convective_fluxes_face_centre_discrete_integral(convective_flux,
                                                                                    subproblem_index,
                                                                                    dist,
                                                                                    dir,
                                                                                    colinearity);
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_thermal_convective_fluxes_face_centre(double& convective_flux,
                                                                                                              const int& subproblem_index,
                                                                                                              const double& dist,
                                                                                                              const int& dir,
                                                                                                              const double& colinearity)
{
  compute_cell_neighbours_thermal_fluxes_face_centre(convective_flux, convection,
                                                     subproblem_index,
                                                     dist,
                                                     dir,
                                                     colinearity);
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_thermal_convective_fluxes_face_centre_discrete_integral(double& convective_flux,
    const int& subproblem_index,
    const double& dist,
    const int& dir,
    const double& colinearity)
{
  compute_cell_neighbours_thermal_fluxes_face_centre_discrete_integral(convective_flux, convection,
                                                                       subproblem_index,
                                                                       dist,
                                                                       dir,
                                                                       colinearity);
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_diffusive_fluxes_to_correct(double& diffusive_flux,
                                                                                                    const int& subproblem_index,
                                                                                                    const double& dist,
                                                                                                    const int& dir,
                                                                                                    const double& colinearity)
{
  if (!discrete_integral_)
    compute_cell_neighbours_thermal_diffusive_fluxes_face_centre(diffusive_flux,
                                                                 subproblem_index,
                                                                 dist,
                                                                 dir,
                                                                 colinearity);
  else
    compute_cell_neighbours_thermal_diffusive_fluxes_face_centre_discrete_integral(diffusive_flux,
                                                                                   subproblem_index,
                                                                                   dist,
                                                                                   dir,
                                                                                   colinearity);
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_thermal_diffusive_fluxes_face_centre(double& diffusive_flux,
                                                                                                             const int& subproblem_index,
                                                                                                             const double& dist,
                                                                                                             const int& dir,
                                                                                                             const double& colinearity)
{
  compute_cell_neighbours_thermal_fluxes_face_centre(diffusive_flux, diffusion,
                                                     subproblem_index,
                                                     dist,
                                                     dir,
                                                     colinearity);
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_thermal_diffusive_fluxes_face_centre_discrete_integral(double& diffusive_flux,
    const int& subproblem_index,
    const double& dist,
    const int& dir,
    const double& colinearity)
{
  compute_cell_neighbours_thermal_fluxes_face_centre_discrete_integral(diffusive_flux, diffusion,
                                                                       subproblem_index,
                                                                       dist,
                                                                       dir,
                                                                       colinearity);
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_thermal_fluxes_face_centre(double& flux,
                                                                                                   const int fluxes_type,
                                                                                                   const int& subproblem_index,
                                                                                                   const double& dist,
                                                                                                   const int& dir,
                                                                                                   const double& colinearity)
{
  int flux_sign[2][6] = FLUX_SIGN;
  double surf_face = 1.;
  for (int c = 0; c < 3; c++)
    if (c!=dir)
      surf_face *= splitting_->get_grid_geometry().get_constant_delta(c);
  surf_face *= flux_sign[fluxes_type][(int) (2*dir)];
  flux = compute_thermal_flux_face_centre(fluxes_type, subproblem_index, dist, dir);
  flux *= surf_face;
  if (neighbours_colinearity_weighting_)
    flux *= colinearity;
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_thermal_fluxes_face_centre_discrete_integral(double& flux,
                                                                                                                     const int fluxes_type,
                                                                                                                     const int& subproblem_index,
                                                                                                                     const double& dist,
                                                                                                                     const int& dir,
                                                                                                                     const double& colinearity)
{
  int flux_sign[2][6] = FLUX_SIGN;
  double surf_face = 1.;
  surf_face *= flux_sign[fluxes_type][(int) (2*dir)];
  DoubleVect discrete_flux_integral;
  discrete_flux_integral = compute_thermal_flux_face_centre_discrete_integral(fluxes_type, subproblem_index, dist, dir);
  for (int val=0; val < discrete_flux_integral.size(); val++)
    flux += discrete_flux_integral[val];
  flux *= surf_face;
  if (neighbours_colinearity_weighting_)
    flux *= colinearity;
}

void Corrige_flux_FT_temperature_subresolution::replace_cell_neighbours_thermal_convective_diffusive_fluxes_faces(const FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_min_max_bool,
                                                                                                                  const FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_fluxes_corrected,
                                                                                                                  const int& fluxes_type)
{
  if (!convection_negligible_ && use_reachable_fluxes_ && fluxes_type == convection)
    {
      replace_cell_neighbours_thermal_fluxes_faces(cell_faces_neighbours_corrected_min_max_bool,
                                                   cell_faces_neighbours_fluxes_corrected,
                                                   convective_diffusive_flux_xyz_min_max_faces_sorted_[0]);
    }
  if (!diffusion_negligible_ && use_reachable_fluxes_ && fluxes_type == diffusion)
    {
      replace_cell_neighbours_thermal_fluxes_faces(cell_faces_neighbours_corrected_min_max_bool,
                                                   cell_faces_neighbours_fluxes_corrected,
                                                   convective_diffusive_flux_xyz_min_max_faces_sorted_[1]);
    }
}

void Corrige_flux_FT_temperature_subresolution::replace_cell_neighbours_thermal_fluxes_faces(const FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_min_max_bool,
                                                                                             const FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_fluxes_corrected,
                                                                                             FixedVector<std::vector<ArrOfDouble>,3>& flux_xyz)
{
  const int ni = cell_faces_neighbours_corrected_min_max_bool[0].ni();
  const int nj = cell_faces_neighbours_corrected_min_max_bool[0].nj();
  const int nk = cell_faces_neighbours_corrected_min_max_bool[0].nk();

  for (int c=0; c<3; c++)
    {
      const int ni_max = (c == 0 ? ni + 1 : ni);
      const int nj_max = (c == 1 ? nj + 1 : nj);
      const int nk_max = (c == 2 ? nk + 1 : nk);
      for (int k = 0; k < nk_max; k++)
        for (int j = 0; j < nj_max; j++)
          for (int i = 0; i < ni_max; i++)
            if (cell_faces_neighbours_corrected_min_max_bool[c](i,j,k))
              {
                flux_xyz[c][k].append_array(cell_faces_neighbours_fluxes_corrected[c](i,j,k));
                index_face_ij_flux_xyz_neighbours_min_max_faces_sorted_[0][c][k].append_array(i);
                index_face_ij_flux_xyz_neighbours_min_max_faces_sorted_[1][c][k].append_array(j);
              }
    }
}

void Corrige_flux_FT_temperature_subresolution::replace_temperature_cell_centre_neighbours(IJK_Field_double& temperature,
                                                                                           IJK_Field_double& temperature_neighbours,
                                                                                           IJK_Field_int& neighbours_weighting,
                                                                                           IJK_Field_double& neighbours_weighting_colinearity) const
{
  if (debug_)
    Cerr << "Corrige flux - Replace temperature cell neighbours - INI" << finl;
  if (distance_cell_faces_from_lrs_ && find_temperature_cell_neighbours_)
    {
      temperature_neighbours.echange_espace_virtuel(temperature_neighbours.ghost());
      neighbours_weighting.echange_espace_virtuel(neighbours_weighting.ghost());
      if (neighbours_colinearity_weighting_)
        neighbours_weighting_colinearity.echange_espace_virtuel(neighbours_weighting_colinearity.ghost());
      const int ni = temperature.ni();
      const int nj = temperature.nj();
      const int nk = temperature.nk();

      const IJK_Splitting& splitting = temperature.get_splitting();
      ArrOfInt corrected_values;
      ArrOfInt out_of_bounds_corrected_values;
      ArrOfDouble out_of_bounds_values;
      corrected_values.set_smart_resize(1);
      out_of_bounds_values.set_smart_resize(1);
      out_of_bounds_corrected_values.set_smart_resize(1);

      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              double neighbours_weighting_ijk;
              const double temperature_old = temperature(i,j,k);
              if (neighbours_colinearity_weighting_)
                neighbours_weighting_ijk = neighbours_weighting_colinearity(i,j,k);
              else
                neighbours_weighting_ijk = (double) neighbours_weighting(i,j,k);
              const int elem = splitting.convert_ijk_cell_to_packed(i,j,k);
              if (neighbours_weighting(i,j,k))
                {
                  temperature(i,j,k) = temperature_neighbours(i,j,k) / neighbours_weighting_ijk;
                  if (smooth_temperature_field_)
                    corrected_values.append_array(elem);

                }
              if (smooth_temperature_field_)
                {
                  const double local_temperature = temperature(i,j,k);
                  const double indic = ref_ijk_ft_->itfce().I()(i,j,k);
                  if (local_temperature > 0)
                    {
                      if (indic < LIQUID_INDICATOR_TEST && indic > VAPOUR_INDICATOR_TEST)
                        if (indic > 0.5)
                          {
                            const int rank = thermal_subproblems_->get_subproblem_index_from_ijk_indices(i,j,k);
                            const double dist_sub_res = thermal_subproblems_->get_dist_cell_interface(rank);
                            if (dist_sub_res > 0)
                              {
                                if (debug_)
                                  {
                                    Cerr << "Thermal subproblem is: " << rank << finl;
                                    Cerr << "Temperature value has the wrong sign at i: " << i << ", j: " << j << ", k: " << k << finl;
                                    Cerr << "Temperature value is: " << local_temperature << finl;
                                    Cerr << "Distance to cell centre is: " << dist_sub_res << finl;
                                    Cerr << "Enforce zero value" << finl;
                                  }
                                // temperature(i,j,k) = 0.;
                                out_of_bounds_corrected_values.append_array(elem);
                                out_of_bounds_values.append_array(local_temperature);
                                temperature(i,j,k) = temperature_old;
                              }
                          }
                      if (indic > LIQUID_INDICATOR_TEST)
                        {
                          if (debug_)
                            {
                              Cerr << "Temperature value has the wrong sign at i: " << i << ", j: " << j << ", k: " << k << finl;
                              Cerr << "Temperature value is: " << local_temperature << finl;
                              Cerr << "Enforce zero value" << finl;
                            }
                          // temperature(i,j,k) = 0.;
                          out_of_bounds_corrected_values.append_array(elem);
                          out_of_bounds_values.append_array(local_temperature);
                          temperature(i,j,k) = temperature_old;
                        }
                    }
                }
            }
      temperature.echange_espace_virtuel(temperature.ghost());
      smooth_temperature_cell_centre_neighbours(temperature,
                                                corrected_values,
                                                out_of_bounds_corrected_values,
                                                out_of_bounds_values,
                                                temperature);
    }

  if (debug_)
    Cerr << "Corrige flux - Replace temperature cell neighbours - END" << finl;
}

void Corrige_flux_FT_temperature_subresolution::smooth_temperature_cell_centre_neighbours(IJK_Field_double& temperature,
                                                                                          ArrOfInt& corrected_values,
                                                                                          ArrOfInt& out_of_bounds_corrected_values,
                                                                                          ArrOfDouble& out_of_bounds_values,
                                                                                          IJK_Field_double& distance) const
{
  // const IJK_Splitting& splitting = temperature.get_splitting();
  //	const Int3 num_elem_ijk = splitting.convert_packed_to_ijk_cell(elem);
  //	(num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
  if (smooth_temperature_field_)
    {
      const IJK_Splitting& splitting = temperature.get_splitting();
      ArrOfDouble temperature_smoothed;
      temperature_smoothed.set_smart_resize(1);
      int counter;
      double temperature_neighbour;
      for (int ielem=0; ielem<out_of_bounds_corrected_values.size_array(); ielem++)
        {
          const int elem = out_of_bounds_corrected_values[ielem];
          const Int3 num_elem_ijk = splitting.convert_packed_to_ijk_cell(elem);
          const int i = num_elem_ijk[DIRECTION_I];
          const int j = num_elem_ijk[DIRECTION_J];
          const int k = num_elem_ijk[DIRECTION_K];
          const int neighbours_i[6] = NEIGHBOURS_I;
          const int neighbours_j[6] = NEIGHBOURS_J;
          const int neighbours_k[6] = NEIGHBOURS_K;
          // const double temperature_old = temperature(i,j,k);
          const double temperature_old_subres = out_of_bounds_values[ielem];
          const double temperature_old = temperature(i,j,k);
          counter=0;
          temperature_neighbour=0;
          for (int l=0; l<6; l++)
            {
              const int ii = neighbours_i[l];
              const int jj = neighbours_j[l];
              const int kk = neighbours_k[l];
              const double indic = ref_ijk_ft_->itfce().I()(i+ii, j+jj, k+kk);
              // if (indic > VAPOUR_INDICATOR_TEST)
              if (indic > LIQUID_INDICATOR_TEST)
                {
                  temperature_neighbour += temperature(i+ii, j+jj, k+kk);
                  counter++;
                }
            }
          double mean_temperature = 1e20;
          if (counter != 0)
            mean_temperature = (temperature_neighbour / counter + temperature_old_subres) * 0.5;
          if (counter != 0)
            {
              if (mean_temperature < 0.)
                temperature_smoothed.append_array(mean_temperature);
              else
                temperature_smoothed.append_array(temperature_neighbour/counter);
            }
          else
            temperature_smoothed.append_array(temperature_old);
        }

      for (int ielem=0; ielem<out_of_bounds_corrected_values.size_array(); ielem++)
        {
          const int elem = out_of_bounds_corrected_values[ielem];
          const Int3 num_elem_ijk = splitting.convert_packed_to_ijk_cell(elem);
          const int i = num_elem_ijk[DIRECTION_I];
          const int j = num_elem_ijk[DIRECTION_J];
          const int k = num_elem_ijk[DIRECTION_K];
          temperature(i,j,k) = temperature_smoothed(ielem);
          if (debug_)
            Cerr << "Smoothed temperature value is:" << temperature_smoothed(ielem) << finl;
        }
    }
}

void Corrige_flux_FT_temperature_subresolution::compute_thermal_convective_fluxes()
{
  if (!discrete_integral_)
    compute_thermal_convective_fluxes_face_centre();
  else
    compute_thermal_convective_fluxes_face_centre_discrete_integral();
}

void Corrige_flux_FT_temperature_subresolution::compute_thermal_diffusive_fluxes()
{
  if (!discrete_integral_)
    compute_thermal_diffusive_fluxes_face_centre();
  else
    compute_thermal_diffusive_fluxes_face_centre_discrete_integral();
}

void Corrige_flux_FT_temperature_subresolution::set_zero_temperature_increment(IJK_Field_double& d_temperature) const
{
  for (int i=0; i<thermal_subproblems_->get_effective_subproblems_counter(); i++)
    {
      int ijk_indices_i, ijk_indices_j, ijk_indices_k;
      thermal_subproblems_->get_subproblem_ijk_indices(ijk_indices_i, ijk_indices_j, ijk_indices_k, i);
      d_temperature(ijk_indices_i, ijk_indices_j, ijk_indices_k) = 0.;
    }
  //  for (int i=0; i<ijk_intersections_subproblems_indices_.size_array(); i++)
  //    {
  //      const int intersection_ijk_cell_index = ijk_intersections_subproblems_indices_[i];
  //      const int ijk_indices_i = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 0);
  //      const int ijk_indices_j = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 1);
  //      const int ijk_indices_k = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 2);
  //      d_temperature(ijk_indices_i, ijk_indices_j, ijk_indices_k) = 0.;
  //    }
}

void Corrige_flux_FT_temperature_subresolution::compute_thermal_convective_fluxes_face_centre()
{
  compute_thermal_fluxes_face_centre(convective_fluxes_, convection);
}

void Corrige_flux_FT_temperature_subresolution::compute_thermal_diffusive_fluxes_face_centre()
{
  compute_thermal_fluxes_face_centre(diffusive_fluxes_, diffusion);
}

void Corrige_flux_FT_temperature_subresolution::compute_thermal_fluxes_face_centre(DoubleVect& fluxes, const int fluxes_type)
{
  const int faces_dir[6] = FACES_DIR;
  int flux_sign[2][6] = FLUX_SIGN;
  const int nb_faces_to_correct = intersection_ijk_cell_->get_nb_faces_to_correct();
  fluxes.reset();
  fluxes.resize(nb_faces_to_correct);
  dist_.reset();
  dist_.resize(nb_faces_to_correct);
  int counter_faces = 0;
  const DoubleTab& dist_interf = intersection_ijk_cell_->dist_pure_faces_interf();
  for (int i=0; i<ijk_intersections_subproblems_indices_.size_array(); i++)
    {
      const int intersection_ijk_cell_index = ijk_intersections_subproblems_indices_[i];
      for (int l=0; l<6; l++)
        {
          const int is_neighbour_pure_liquid = intersection_ijk_cell_->get_ijk_pure_face_neighbours(intersection_ijk_cell_index, l);
          double surf_face = 1.;
          const int dir = faces_dir[l];
          if (is_neighbour_pure_liquid)
            {
              for (int c = 0; c < 3; c++)
                if (c!= faces_dir[l])
                  surf_face *= splitting_->get_grid_geometry().get_constant_delta(c);
              surf_face *= flux_sign[fluxes_type][l];
              const double dist = dist_interf(intersection_ijk_cell_index, l);
              const double dist_sub_res = thermal_subproblems_->get_dist_faces_interface(i)[l];
              double local_flux_face = 0.;
              if (distance_cell_faces_from_lrs_)
                local_flux_face = compute_thermal_flux_face_centre(fluxes_type, i, dist_sub_res, dir);
              else
                local_flux_face = compute_thermal_flux_face_centre(fluxes_type, i, dist, dir);
              const double flux_face = local_flux_face * surf_face;
              fluxes[counter_faces] = flux_face;
              dist_[counter_faces] = dist;
              counter_faces++;
            }
        }
    }
  /*
   * Useless if a treat a sub-problem per mixed cells
   * May be useful if a treat one subproblem per interface portion
   */
  // check_pure_fluxes_duplicates(convective_fluxes_, convective_fluxes_unique_, pure_face_unique_, 1);
}

double Corrige_flux_FT_temperature_subresolution::compute_thermal_flux_face_centre(const int fluxes_type, const int& index_subproblem, const double& dist, const int& dir)
{
  switch(fluxes_type)
    {
    case convection:
      return thermal_subproblems_->get_temperature_times_velocity_profile_at_point(index_subproblem, dist, dir);
    case diffusion:
      return thermal_subproblems_->get_temperature_gradient_times_conductivity_profile_at_point(index_subproblem, dist, dir);
    default:
      return thermal_subproblems_->get_temperature_times_velocity_profile_at_point(index_subproblem, dist, dir);
    }
}

void Corrige_flux_FT_temperature_subresolution::compute_thermal_convective_fluxes_face_centre_discrete_integral()
{
  compute_thermal_fluxes_face_centre_discrete_integral(convective_fluxes_, convection);
}

void Corrige_flux_FT_temperature_subresolution::compute_thermal_diffusive_fluxes_face_centre_discrete_integral()
{
  compute_thermal_fluxes_face_centre_discrete_integral(diffusive_fluxes_, diffusion);
}

void Corrige_flux_FT_temperature_subresolution::compute_thermal_fluxes_face_centre_discrete_integral(DoubleVect& fluxes, const int fluxes_type)
{
  int faces_dir[6] = FACES_DIR;
  int flux_sign[2][6] = FLUX_SIGN;
  const int nb_faces_to_correct = intersection_ijk_cell_->get_nb_faces_to_correct();
  fluxes.reset();
  fluxes.resize(nb_faces_to_correct);
  dist_.reset();
  dist_.resize(nb_faces_to_correct);
  const DoubleTab& dist_interf = intersection_ijk_cell_->dist_pure_faces_interf();
  int counter_faces = 0;
  // const DoubleTab& pos_pure_faces_interf = intersection_ijk_cell_->pos_pure_faces_interf();
  // positions.resize(nb_diph, 3, 6);
  for (int i=0; i<ijk_intersections_subproblems_indices_.size_array(); i++)
    {
      // double surface = 0.;
      const int intersection_ijk_cell_index = ijk_intersections_subproblems_indices_[i];
      for (int l=0; l<6; l++)
        {
          const int dir = faces_dir[l];
          const int is_neighbour_pure_liquid = intersection_ijk_cell_->get_ijk_pure_face_neighbours(intersection_ijk_cell_index, l);
          if (is_neighbour_pure_liquid)
            {
              const double dist = dist_interf(intersection_ijk_cell_index, l);
              const double dist_sub_res = thermal_subproblems_->get_dist_faces_interface(i)[l];
              DoubleVect discrete_flux_integral;
              if (distance_cell_faces_from_lrs_)
                discrete_flux_integral = compute_thermal_flux_face_centre_discrete_integral(fluxes_type, i, dist_sub_res, dir);
              else
                discrete_flux_integral = compute_thermal_flux_face_centre_discrete_integral(fluxes_type, i, dist, dir);
              double flux_face = 0;
              for (int val=0; val < discrete_flux_integral.size(); val++)
                flux_face += discrete_flux_integral[val];
              flux_face *= flux_sign[fluxes_type][l];
              fluxes[counter_faces] = flux_face;
              dist_[counter_faces] = dist;
              counter_faces++;
            }
        }
    }
  /*
   * Useless if a treat a sub-problem per mixed cells
   * May be useful if a treat one subproblem per interface portion
   */
  // check_pure_fluxes_duplicates(convective_fluxes_, convective_fluxes_unique_, pure_face_unique_, 0);
}

DoubleVect Corrige_flux_FT_temperature_subresolution::compute_thermal_flux_face_centre_discrete_integral(const int fluxes_type, const int& index_subproblem, const double& dist, const int& dir)
{
  switch(fluxes_type)
    {
    case convection:
      return thermal_subproblems_->get_temperature_times_velocity_profile_discrete_integral_at_point(index_subproblem, dist, levels_, dir);
    case diffusion:
      return thermal_subproblems_->get_temperature_gradient_times_conductivity_profile_discrete_integral_at_point(index_subproblem, dist, levels_, dir);
    default:
      return thermal_subproblems_->get_temperature_times_velocity_profile_discrete_integral_at_point(index_subproblem, dist, levels_, dir);
    }
}

/* void Corrige_flux_FT_temperature_subresolution::get_discrete_surface_at_level(const int& dir, const int& level)
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  double surface = 0.;
  switch(dir)
    {
    case 0:
      surface = (dy*dz);
      break;
    case 1:
      surface = (dx*dz);
      break;
    case 2:
      surface = (dx*dy);
      break;
    default:
      surface = (dx*dy);
      break;
    }
  surface /= pow(pow(2, level), 2);
} */

void Corrige_flux_FT_temperature_subresolution::compute_min_max_ijk_reachable_fluxes(const FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_all_bool,
                                                                                     const IJK_Field_int& neighbours_temperature_to_correct,
                                                                                     FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_min_max_bool,
                                                                                     const int& max_flux_per_dir,
                                                                                     const int& check_cell_center_neighbour,
                                                                                     const int& remove_external_neighbour_values,
                                                                                     IJK_Field_int& neighbours_temperature_to_correct_trimmed)
{
  /*
   * TODO: Work only for static bubbles
   */
  if (distance_cell_faces_from_lrs_ && find_reachable_fluxes_)
    {

      if (!keep_first_reachable_fluxes_)
        {
          FixedVector<IJK_Field_local_int, 3> cell_faces_neighbours_corrected_min_max_bool_local;

          remove_min_max_ijk_reachable_fluxes_discontinuous(cell_faces_neighbours_corrected_all_bool,
                                                            cell_faces_neighbours_corrected_min_max_bool_local);

          const int ni = cell_faces_neighbours_corrected_all_bool[0].ni();
          const int nj = cell_faces_neighbours_corrected_all_bool[0].nj();
          const int nk = cell_faces_neighbours_corrected_all_bool[0].nk();
          int i,j,k,c;
          for (c = 0; c < 3; c++)
            cell_faces_neighbours_corrected_min_max_bool[c].data() = 0.;
          int c_ini = 0;
          int c_end = 3;

          if (remove_external_neighbour_values)
            {
              for (k = 0; k < nk; k++)
                for (j = 0; j < nj; j++)
                  for (i = 0; i < ni; i++)
                    neighbours_temperature_to_correct_trimmed(i,j,k) = neighbours_temperature_to_correct(i,j,k);
              neighbours_temperature_to_correct_trimmed.echange_espace_virtuel(neighbours_temperature_to_correct_trimmed.ghost());
            }
          /*
           * i_min, i_max
           */
          if (max_flux_per_dir)
            {
              c_ini = 0;
              c_end = 1;
            }
          int index_neighbour_cell = 0;
          bool index_found = false;
          bool centre_neighbour_bool = false;
          for (c = c_ini; c < c_end; c++)
            for (k = 0; k < nk; k++)
              for (j = 0; j < nj; j++)
                {
                  for (i = 0; i < ni; i++)
                    {
                      if (neighbours_temperature_to_correct(i,j,k) && !index_found)
                        {
                          index_neighbour_cell = i;
                          index_found = true;
                        }
                      if(cell_faces_neighbours_corrected_min_max_bool_local[c](i,j,k))
                        {
                          if (neighbours_temperature_to_correct(i - 1,j,k) && check_cell_center_neighbour && !remove_external_neighbour_values)
                            centre_neighbour_bool = true;
                          if (!centre_neighbour_bool)
                            cell_faces_neighbours_corrected_min_max_bool[c](i,j,k) = cell_faces_neighbours_corrected_min_max_bool_local[c](i,j,k);
                          centre_neighbour_bool = false;
                          index_found = false;
                          for(int ii = index_neighbour_cell; ii < i; ii++)
                            neighbours_temperature_to_correct_trimmed(ii,j,k) = 0;
                          break;
                        }
                    }
                  for (i = ni - 1; i >=0; i--)
                    {
                      if (neighbours_temperature_to_correct(i,j,k) && !index_found)
                        {
                          index_neighbour_cell = i;
                          index_found = true;
                        }
                      if(cell_faces_neighbours_corrected_min_max_bool_local[c](i,j,k))
                        {
                          if (neighbours_temperature_to_correct(i,j,k) && check_cell_center_neighbour && !remove_external_neighbour_values)
                            centre_neighbour_bool = true;
                          if (!centre_neighbour_bool)
                            cell_faces_neighbours_corrected_min_max_bool[c](i,j,k) = cell_faces_neighbours_corrected_min_max_bool_local[c](i,j,k);
                          centre_neighbour_bool = false;
                          index_found = false;
                          for(int ii = index_neighbour_cell; ii >= i; ii--)
                            neighbours_temperature_to_correct_trimmed(ii,j,k) = 0;
                          break;
                        }
                    }
                }
          /*
           * j_min, j_max
           */
          if (max_flux_per_dir)
            {
              c_ini = 1;
              c_end = 2;
            }
          for (c = c_ini; c < c_end; c++)
            for (k = 0; k < nk; k++)
              for (i = 0; i < ni; i++)
                {
                  for (j = 0; j < nj; j++)
                    {
                      if (neighbours_temperature_to_correct(i,j,k) && !index_found)
                        {
                          index_neighbour_cell = j;
                          index_found = true;
                        }
                      if(cell_faces_neighbours_corrected_min_max_bool_local[c](i,j,k))
                        {
                          if (neighbours_temperature_to_correct(i,j - 1,k) && check_cell_center_neighbour && !remove_external_neighbour_values)
                            centre_neighbour_bool = true;
                          if (!centre_neighbour_bool)
                            cell_faces_neighbours_corrected_min_max_bool[c](i,j,k) = cell_faces_neighbours_corrected_min_max_bool_local[c](i,j,k);
                          centre_neighbour_bool = false;
                          index_found = false;
                          for(int jj = index_neighbour_cell; jj < j; jj++)
                            neighbours_temperature_to_correct_trimmed(i,jj,k) = 0;
                          break;
                        }
                    }
                  for (j = nj - 1; j >=0; j--)
                    {
                      if (neighbours_temperature_to_correct(i,j,k) && !index_found)
                        {
                          index_neighbour_cell = j;
                          index_found = true;
                        }
                      if(cell_faces_neighbours_corrected_min_max_bool_local[c](i,j,k))
                        {
                          if (neighbours_temperature_to_correct(i,j,k) && check_cell_center_neighbour && !remove_external_neighbour_values)
                            centre_neighbour_bool = true;
                          if (!centre_neighbour_bool)
                            cell_faces_neighbours_corrected_min_max_bool[c](i,j,k) = cell_faces_neighbours_corrected_min_max_bool_local[c](i,j,k);
                          centre_neighbour_bool = false;
                          index_found = false;
                          for(int jj = index_neighbour_cell; jj >= j; jj--)
                            neighbours_temperature_to_correct_trimmed(i,jj,k) = 0;
                          break;
                        }
                    }
                }
          /*
           * k_min, k_max
           */
          if (max_flux_per_dir)
            {
              c_ini = 2;
              c_end = 3;
            }
          for (c = c_ini; c < c_end; c++)
            for (i = 0; i < ni; i++)
              for (j = 0; j < nj; j++)
                {
                  for (k = 0; k < nk; k++)
                    {
                      if (neighbours_temperature_to_correct(i,j,k) && !index_found)
                        {
                          index_neighbour_cell = k;
                          index_found = true;
                        }
                      if(cell_faces_neighbours_corrected_min_max_bool_local[c](i,j,k))
                        {
                          if (neighbours_temperature_to_correct(i,j,k - 1) && check_cell_center_neighbour && !remove_external_neighbour_values)
                            centre_neighbour_bool = true;
                          if (!centre_neighbour_bool)
                            cell_faces_neighbours_corrected_min_max_bool[c](i,j,k) = cell_faces_neighbours_corrected_min_max_bool_local[c](i,j,k);
                          centre_neighbour_bool = false;
                          index_found = false;
                          for(int kk = index_neighbour_cell; kk < k; kk++)
                            neighbours_temperature_to_correct_trimmed(i,j,kk) = 0;
                          break;
                        }
                    }
                  for (k = nk - 1; k >=0; k--)
                    {
                      if (neighbours_temperature_to_correct(i,j,k) && !index_found)
                        {
                          index_neighbour_cell = k;
                          index_found = true;
                        }
                      if(cell_faces_neighbours_corrected_min_max_bool_local[c](i,j,k))
                        {
                          if (neighbours_temperature_to_correct(i,j,k) && check_cell_center_neighbour && !remove_external_neighbour_values)
                            centre_neighbour_bool = true;
                          if (!centre_neighbour_bool)
                            cell_faces_neighbours_corrected_min_max_bool[c](i,j,k) = cell_faces_neighbours_corrected_min_max_bool_local[c](i,j,k);
                          centre_neighbour_bool = false;
                          index_found = false;
                          for(int kk = index_neighbour_cell; kk >= k; kk--)
                            neighbours_temperature_to_correct_trimmed(i,j,kk) = 0;
                          break;
                        }
                    }
                }
          cell_faces_neighbours_corrected_min_max_bool.echange_espace_virtuel();
        }
    }

}

void Corrige_flux_FT_temperature_subresolution::compute_min_max_ijk_any_reachable_fluxes(const FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_all_bool,
                                                                                         const IJK_Field_int& neighbours_temperature_to_correct,
                                                                                         FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_min_max_bool,
                                                                                         const int& max_flux_per_dir,
                                                                                         const int& check_cell_center_neighbour,
                                                                                         const int& remove_external_neighbour_values,
                                                                                         IJK_Field_int& neighbours_temperature_to_correct_trimmed)
{
  if (distance_cell_faces_from_lrs_ && find_reachable_fluxes_)
    {

      if (!keep_first_reachable_fluxes_)
        {
          FixedVector<IJK_Field_local_int, 3> cell_faces_neighbours_corrected_min_max_bool_local;

          remove_min_max_ijk_reachable_fluxes_discontinuous(cell_faces_neighbours_corrected_all_bool,
                                                            cell_faces_neighbours_corrected_min_max_bool_local);

          int c;
          for (c = 0; c < 3; c++)
            cell_faces_neighbours_corrected_min_max_bool[c].data() = 0.;

          const int ni = cell_faces_neighbours_corrected_all_bool[0].ni();
          const int nj = cell_faces_neighbours_corrected_all_bool[0].nj();
          const int nk = cell_faces_neighbours_corrected_all_bool[0].nk();
          int i,j,k;
          int c_ini = 0;
          int c_end = 3;
          if (remove_external_neighbour_values)
            {
              for (k = 0; k < nk; k++)
                for (j = 0; j < nj; j++)
                  for (i = 0; i < ni; i++)
                    neighbours_temperature_to_correct_trimmed(i,j,k) = neighbours_temperature_to_correct(i,j,k);
              neighbours_temperature_to_correct_trimmed.echange_espace_virtuel(neighbours_temperature_to_correct_trimmed.ghost());
              /*
              * i_min, i_max, j_min, j_max, k_min, k_max
              */
              ArrOfInt indices_found_ini, indices_found_end;
              ArrOfInt indices_flux_found_ini, indices_flux_found_end;
              ArrOfInt indices_to_remove;
              ArrOfInt indices_fluxes_to_remove;
              FixedVector<ArrOfInt,2> indices_sorted;
              FixedVector<ArrOfInt,2> indices_fluxes_sorted;
              indices_found_ini.set_smart_resize(1);
              indices_found_end.set_smart_resize(1);
              indices_flux_found_ini.set_smart_resize(1);
              indices_flux_found_end.set_smart_resize(1);

//              if (max_flux_per_dir)
//                {
//                  c_ini = 0;
//                  c_end = 1;
//                }
              const int * n_first = nullptr;
              const int * n_bis = nullptr;
              const int * n_ter = nullptr;
              int * i_ref = nullptr;
              int * j_ref = nullptr;
              int * k_ref = nullptr;
              int * i_ref_remove = nullptr;
              int * j_ref_remove = nullptr;
              int * k_ref_remove = nullptr;
              int index_first = 0, index_bis = 0, index_ter = 0;
              int val = 0;
              const int factor_pos[3][3] = {{1, 0, 0},{0, 1, 0},{0, 0, 1}};
              for (c = c_ini; c < c_end; c++)
                {
                  switch(c)
                    {
                    case 0:
                      n_first = &nk;
                      n_bis = &nj;
                      n_ter = &ni;
                      k_ref = &index_first;
                      j_ref = &index_bis;
                      i_ref = &index_ter;
                      k_ref_remove = &index_first;
                      j_ref_remove = &index_bis;
                      i_ref_remove = &val;
                      break;
                    case 1:
                      n_first = &nk;
                      n_bis = &ni;
                      n_ter = &nj;
                      k_ref = &index_first;
                      i_ref = &index_bis;
                      j_ref = &index_ter;
                      k_ref_remove = &index_first;
                      i_ref_remove = &index_bis;
                      j_ref_remove = &val;
                      break;
                    case 2:
                      n_first = &ni;
                      n_bis = &nj;
                      n_ter = &nk;
                      i_ref = &index_first;
                      j_ref = &index_bis;
                      k_ref = &index_ter;
                      i_ref_remove = &index_first;
                      j_ref_remove = &index_bis;
                      k_ref_remove = &val;
                      break;
                    default:
                      break;
                    }
                  for (index_first = 0; index_first < *n_first; index_first++)
                    for (index_bis = 0; index_bis < *n_bis; index_bis++)
                      {
                        indices_found_ini.reset();
                        indices_found_end.reset();
                        indices_flux_found_ini.reset();
                        indices_flux_found_end.reset();
                        indices_to_remove.reset();
                        indices_fluxes_to_remove.reset();
                        for (int l=0; l<2; l++)
                          {
                            indices_sorted[l].reset();
                            indices_fluxes_sorted[l].reset();
                          }
                        for (index_ter = 0; index_ter < *n_ter; index_ter++)
                          {
                            if (neighbours_temperature_to_correct(*i_ref,*j_ref,*k_ref))
                              {
                                if (!neighbours_temperature_to_correct(*i_ref + factor_pos[c][0],*j_ref + factor_pos[c][1],*k_ref + factor_pos[c][2]))
                                  indices_found_end.append_array(index_ter);
                                if (!neighbours_temperature_to_correct(*i_ref - factor_pos[c][0],*j_ref - factor_pos[c][1],*k_ref - factor_pos[c][2]))
                                  indices_found_ini.append_array(index_ter);
                              }
                            if(cell_faces_neighbours_corrected_min_max_bool_local[c](*i_ref,*j_ref,*k_ref))
                              {
                                if (!cell_faces_neighbours_corrected_min_max_bool_local[c](*i_ref + factor_pos[c][0],*j_ref + factor_pos[c][1],*k_ref + factor_pos[c][2]))
                                  indices_flux_found_end.append_array(index_ter);
                                if (!cell_faces_neighbours_corrected_min_max_bool_local[c](*i_ref - factor_pos[c][0],*j_ref - factor_pos[c][1],*k_ref - factor_pos[c][2]))
                                  indices_flux_found_ini.append_array(index_ter);
                              }
                          }

                        const int indices_test = indices_found_ini.size_array() || indices_found_end.size_array();
                        if (indices_test)
                          {
                            sort_ini_end_arrays(indices_found_ini,
                                                indices_found_end,
                                                indices_sorted,
                                                *n_ter);
                          }

                        const int indices_flux_test = indices_flux_found_ini.size_array() || indices_flux_found_end.size_array();
                        if (indices_flux_test)
                          {
//                            sort_ini_end_arrays(indices_flux_found_ini,
//                                                indices_flux_found_end,
//                                                indices_fluxes_sorted,
//                                                *n_ter);
                            sort_ini_end_arrays(indices_flux_found_ini,
                                                indices_flux_found_end,
                                                indices_fluxes_sorted,
                                                (*n_ter) + 1);
                          }
                        if (indices_test || indices_flux_test)
                          remove_non_overlapping_fluxes_values(indices_sorted,
                                                               indices_fluxes_sorted,
                                                               indices_to_remove,
                                                               indices_fluxes_to_remove,
                                                               index_first, index_bis, c);

//                        for (index_ter = 0; index_ter < indices_fluxes_sorted[0].size_array(); index_ter++)
//                          {
//                            val = indices_fluxes_sorted[0](index_ter);
//                            cell_faces_neighbours_corrected_min_max_bool[c](*i_ref_remove, *j_ref_remove, *k_ref_remove) =
//                              cell_faces_neighbours_corrected_min_max_bool_local[c](*i_ref_remove, *j_ref_remove, *k_ref_remove);
//                            val = indices_fluxes_sorted[1](index_ter);
//                            cell_faces_neighbours_corrected_min_max_bool[c](*i_ref_remove, *j_ref_remove, *k_ref_remove) =
//                              cell_faces_neighbours_corrected_min_max_bool_local[c](*i_ref_remove, *j_ref_remove, *k_ref_remove);
//                          }
//
//                        /*
//                         * Keep max min flux values to build a convex shape
//                         */
                        for (index_ter = 0; index_ter < indices_flux_found_ini.size_array(); index_ter++)
                          {
                            val = indices_flux_found_ini(index_ter);
                            const double indic_ini = ref_ijk_ft_->itfce().I()(*i_ref_remove,
                                                                              *j_ref_remove,
                                                                              *k_ref_remove);
                            const double indic_ini_prev = ref_ijk_ft_->itfce().I()(*i_ref_remove - factor_pos[c][0],
                                                                                   *j_ref_remove - factor_pos[c][1],
                                                                                   *k_ref_remove - factor_pos[c][2]);
                            if(indic_ini > LIQUID_INDICATOR_TEST && indic_ini_prev > LIQUID_INDICATOR_TEST)
                              {
                                cell_faces_neighbours_corrected_min_max_bool[c](*i_ref_remove, *j_ref_remove, *k_ref_remove) =
                                  cell_faces_neighbours_corrected_min_max_bool_local[c](*i_ref_remove, *j_ref_remove, *k_ref_remove);
                                // break;
                              }
                          }
                        for (index_ter = 0; index_ter < indices_flux_found_end.size_array(); index_ter++)
                          {
                            val = indices_flux_found_end(index_ter);
                            const double indic_end = ref_ijk_ft_->itfce().I()(*i_ref_remove,
                                                                              *j_ref_remove,
                                                                              *k_ref_remove);
                            const double indic_end_prev = ref_ijk_ft_->itfce().I()(*i_ref_remove - factor_pos[c][0],
                                                                                   *j_ref_remove - factor_pos[c][1],
                                                                                   *k_ref_remove - factor_pos[c][2]);
                            if(indic_end > LIQUID_INDICATOR_TEST && indic_end_prev > LIQUID_INDICATOR_TEST)
                              {
                                cell_faces_neighbours_corrected_min_max_bool[c](*i_ref_remove, *j_ref_remove, *k_ref_remove) =
                                  cell_faces_neighbours_corrected_min_max_bool_local[c](*i_ref_remove, *j_ref_remove, *k_ref_remove);
                                //break;
                              }
                          }

                        /*
                         * Remove temperature and non-consistent fluxes values
                         */
                        for (index_ter = 0; index_ter < indices_to_remove.size_array(); index_ter++)
                          {
                            val = indices_to_remove(index_ter);
                            neighbours_temperature_to_correct_trimmed(*i_ref_remove, *j_ref_remove, *k_ref_remove) = 0;
                          }

                        for (index_ter = 0; index_ter < indices_fluxes_to_remove.size_array(); index_ter++)
                          {
                            val = indices_fluxes_to_remove(index_ter);
                            cell_faces_neighbours_corrected_min_max_bool[c](*i_ref_remove, *j_ref_remove, *k_ref_remove) = 0;
                          }
                      }
                }
            }
        }
      Cerr << "Echange espace virtuel" << finl;
      cell_faces_neighbours_corrected_min_max_bool.echange_espace_virtuel();
    }
}

void Corrige_flux_FT_temperature_subresolution::sort_ini_end_arrays(ArrOfInt& indices_found_ini,
                                                                    ArrOfInt& indices_found_end,
                                                                    FixedVector<ArrOfInt,2>& indices_sorted,
                                                                    const int& max_n_layer)
{
  int counter_ini = 0;
  int counter_end = 0;
  FixedVector<int*, 2> counters;
  counters[0] = &counter_ini;
  counters[1] = &counter_end;
  FixedVector<ArrOfInt*, 2> indices_references;
  indices_references[0] = &indices_found_ini;
  indices_references[1] = &indices_found_end;
  int size_ini = indices_found_ini.size_array();
  int size_end = indices_found_end.size_array();
  FixedVector<int*, 2> size_references;
  size_references[0] = &size_ini;
  size_references[1] = &size_end;
  int max_counter = 0;
  const int max_index = size_ini + size_end;
  for (int l=0; l<2; l++)
    indices_sorted[l].set_smart_resize(1);
  std::vector<int> index_tmp;
  while(counter_ini < size_ini || counter_end < size_end)
    {
      for (int l=0; l<2; l++)
        {
          const int index_int_tmp = (*counters[l] < *size_references[l]) ? (*indices_references[l])(*counters[l]) : (int) 1e20;
          index_tmp.push_back(index_int_tmp);
        }
      const int index_min = (int) std::distance(index_tmp.begin(), std::min_element(index_tmp.begin(), index_tmp.end()));
      if (!index_min)
        {
          indices_sorted[0].append_array((*indices_references[index_min])(*counters[index_min]));
          if (max_counter == max_index - 1)
            indices_sorted[1].append_array(max_n_layer - 1);
        }
      else
        {
          indices_sorted[1].append_array((*indices_references[index_min])(*counters[index_min]));
          if (!indices_sorted[0].size_array())
            indices_sorted[0].append_array(0);
        }
      (*(counters[index_min]))++;
      max_counter++;
      index_tmp.clear();
    }
  assert(indices_sorted[0].size_array() == indices_sorted[1].size_array());
}

void Corrige_flux_FT_temperature_subresolution::remove_non_overlapping_fluxes_values(const FixedVector<ArrOfInt,2>& indices_sorted,
                                                                                     const FixedVector<ArrOfInt,2>& indices_fluxes_sorted,
                                                                                     ArrOfInt& indices_to_remove,
                                                                                     ArrOfInt& indices_fluxes_to_remove,
                                                                                     int& index_bis,
                                                                                     int& index_ter,
                                                                                     const int& dir)
{
  int * i = nullptr;
  int * j = nullptr;
  int * k = nullptr;
  int val = 0;
  switch(dir)
    {
    case 0:
      j = &index_ter;
      k = &index_bis;
      i = &val;
      break;
    case 1:
      i = &index_ter;
      k = &index_bis;
      j = &val;
      break;
    case 2:
      i = &index_bis;
      j = &index_ter;
      k = &val;
      break;
    default:
      break;
    }
  const int factor_pos[3][3] = {{1, 0, 0},{0, 1, 0},{0, 0, 1}};
  std::vector<int> is_in_interval = {0, 0};
  indices_to_remove.set_smart_resize(1);
  indices_fluxes_to_remove.set_smart_resize(1);
  int ll;
  int m = 0;
  for(int l=0; l<indices_sorted[0].size_array(); l++)
    if (m < indices_fluxes_sorted[0].size_array())
      // for (int m=0; m<indices_fluxes_sorted[0].size_array(); m++)
      {
        for (int n=0; n<2; n++)
          {
            // use indicator
            val = indices_fluxes_sorted[n][m];
            const double indic = ref_ijk_ft_->itfce().I()(*i, *j, *k);
            const double indic_prev = ref_ijk_ft_->itfce().I()(*i - factor_pos[dir][0], *j - factor_pos[dir][1], *k - factor_pos[dir][2]);
            const int indic_test = (indic > LIQUID_INDICATOR_TEST && indic_prev > LIQUID_INDICATOR_TEST);
            if (indic_test)
              switch(n)
                {
                case 0:
                  is_in_interval[n] = (indices_fluxes_sorted[n][m] > indices_sorted[n][l]) ? 1 : 0;
                  break;
                case 1:
                  is_in_interval[n] = (indices_fluxes_sorted[n][m] <= indices_sorted[n][l]) ? 1 : 0;
                  break;
                default:
                  break;
                }
          }
        if (is_in_interval[0] || is_in_interval[1])
          {
            if (is_in_interval[0])
              for (ll=indices_sorted[0][l]; ll<indices_fluxes_sorted[0][m]; ll++)
                indices_to_remove.append_array(ll);
            if (is_in_interval[1])
              for (ll=indices_fluxes_sorted[1][m]; ll<=indices_sorted[1][l]; ll++)
                indices_to_remove.append_array(ll);
            m++;
          }
        is_in_interval = {0, 0};
      }
  for (m=0; m<indices_fluxes_sorted[0].size_array(); m++)
    if (indices_fluxes_sorted[0][m] == indices_fluxes_sorted[1][m])
      indices_fluxes_to_remove.append_array(indices_fluxes_sorted[0][m]);
}

void Corrige_flux_FT_temperature_subresolution::remove_min_max_ijk_reachable_fluxes_discontinuous(const FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_all_bool,
                                                                                                  FixedVector<IJK_Field_local_int, 3>& cell_faces_neighbours_corrected_min_max_bool)
{
  const int ni = cell_faces_neighbours_corrected_all_bool[0].ni();
  const int nj = cell_faces_neighbours_corrected_all_bool[0].nj();
  const int nk = cell_faces_neighbours_corrected_all_bool[0].nk();

  int i,j,k,c;
  ArrOfInt neighbour_left_right(3);
  const int nb_ghost = 1;

  for (c = 0; c < 3; c++)
    cell_faces_neighbours_corrected_min_max_bool[c].allocate(ni,nj,nk,nb_ghost);

  for (c = 0; c < 3; c++)
    for (k = -nb_ghost; k < nk + nb_ghost; k++)
      for (j = -nb_ghost; j < nj + nb_ghost; j++)
        for (i = -nb_ghost; i < ni + nb_ghost; i++)
          cell_faces_neighbours_corrected_min_max_bool[c](i,j,k) = cell_faces_neighbours_corrected_all_bool[c](i,j,k);

  const int neighbouring_face_index_from_sign[2] = {1, -1};
  const int neighbouring_cell_index_from_sign[2] = {0, -1};
  const int neighbouring_face_first_second_index_from_sign[2] = {1, 0};
  int neighbouring_face_index;
  int n_sign_mean, n_sign_first_dir, n_sign_second_dir;

  for (c = 0; c < 3; c++)
    for (k = 0; k < nk; k++)
      for (j = 0; j < nj; j++)
        for (i = 0; i < ni; i++)
          {
            if (cell_faces_neighbours_corrected_all_bool[c](i,j,k))
              {
                switch(c)
                  {
                  case 0:
                    n_sign_mean = signbit((*eulerian_normal_vectors_ns_normed_)[0](i-1,j,k) + (*eulerian_normal_vectors_ns_normed_)[0](i,j,k));
                    neighbouring_face_index = neighbouring_face_index_from_sign[n_sign_mean];
                    if (!cell_faces_neighbours_corrected_min_max_bool[0](i + neighbouring_face_index,j,k))
                      {
                        const int neighbouring_cell_index = neighbouring_cell_index_from_sign[n_sign_mean];
                        n_sign_first_dir = signbit((*eulerian_normal_vectors_ns_normed_)[1](i+neighbouring_cell_index,j,k));
                        n_sign_second_dir = signbit((*eulerian_normal_vectors_ns_normed_)[2](i+neighbouring_cell_index,j,k));

                        const int j_face = neighbouring_face_first_second_index_from_sign[n_sign_first_dir];
                        cell_faces_neighbours_corrected_min_max_bool[1](i+neighbouring_cell_index, j+j_face, k) = 0;
                        const int k_face = neighbouring_face_first_second_index_from_sign[n_sign_second_dir];
                        cell_faces_neighbours_corrected_min_max_bool[2](i+neighbouring_cell_index, j, k+k_face) = 0;
                      }
                    break;
                  case 1:
                    n_sign_mean = signbit((*eulerian_normal_vectors_ns_normed_)[1](i,j-1,k) + (*eulerian_normal_vectors_ns_normed_)[1](i,j,k));
                    neighbouring_face_index = neighbouring_face_index_from_sign[n_sign_mean];
                    if (!cell_faces_neighbours_corrected_min_max_bool[1](i,j + neighbouring_face_index,k))
                      {
                        const int neighbouring_cell_index = neighbouring_cell_index_from_sign[n_sign_mean];
                        n_sign_first_dir = signbit((*eulerian_normal_vectors_ns_normed_)[0](i, j+neighbouring_cell_index, k));
                        n_sign_second_dir = signbit((*eulerian_normal_vectors_ns_normed_)[2](i, j+neighbouring_cell_index, k));

                        const int i_face = neighbouring_face_first_second_index_from_sign[n_sign_first_dir];
                        cell_faces_neighbours_corrected_min_max_bool[0](i+i_face, j+neighbouring_cell_index, k) = 0;
                        const int k_face = neighbouring_face_first_second_index_from_sign[n_sign_second_dir];
                        cell_faces_neighbours_corrected_min_max_bool[2](i, j+neighbouring_cell_index, k+k_face) = 0;

                      }
                    break;
                  case 2:
                    n_sign_mean = signbit((*eulerian_normal_vectors_ns_normed_)[2](i,j,k-1) + (*eulerian_normal_vectors_ns_normed_)[2](i,j,k));
                    neighbouring_face_index = neighbouring_face_index_from_sign[n_sign_mean];
                    if (!cell_faces_neighbours_corrected_min_max_bool[2](i,j,k + neighbouring_face_index))
                      {
                        const int neighbouring_cell_index = neighbouring_cell_index_from_sign[n_sign_mean];
                        n_sign_first_dir = signbit((*eulerian_normal_vectors_ns_normed_)[0](i, j, k+neighbouring_cell_index));
                        n_sign_second_dir = signbit((*eulerian_normal_vectors_ns_normed_)[1](i, j, k+neighbouring_cell_index));

                        const int i_face = neighbouring_face_first_second_index_from_sign[n_sign_first_dir];
                        cell_faces_neighbours_corrected_min_max_bool[0](i+i_face, j, k+neighbouring_cell_index) = 0;
                        const int j_face = neighbouring_face_first_second_index_from_sign[n_sign_second_dir];
                        cell_faces_neighbours_corrected_min_max_bool[1](i, j+j_face, k+neighbouring_cell_index) = 0;
                      }
                    break;
                  default:
                    break;
                  }
              }
          }
}


void Corrige_flux_FT_temperature_subresolution::compute_ijk_pure_faces_indices()
{
  /*
   * Be careful, the ijk_intersection class is not sorting the faces the same way
   */
  const int faces_dir[6] = FACES_DIR;
  const int neighbours_faces_i[6] = NEIGHBOURS_FACES_I;
  const int neighbours_faces_j[6] = NEIGHBOURS_FACES_J;
  const int neighbours_faces_k[6] = NEIGHBOURS_FACES_K;

  int nb_faces_to_correct = 0;
  const int nb_faces_to_correct_from_ijk = intersection_ijk_cell_->get_nb_faces_to_correct();

  IntVect& i_pure_face_to_correct = ijk_faces_to_correct_[0];
  IntVect& j_pure_face_to_correct = ijk_faces_to_correct_[1];
  IntVect& k_pure_face_to_correct = ijk_faces_to_correct_[2];
  IntVect& dir_pure_face_to_correct = ijk_faces_to_correct_[3];

  i_pure_face_to_correct.resize(nb_faces_to_correct_from_ijk);
  j_pure_face_to_correct.resize(nb_faces_to_correct_from_ijk);
  k_pure_face_to_correct.resize(nb_faces_to_correct_from_ijk);
  dir_pure_face_to_correct.resize(nb_faces_to_correct_from_ijk);
  /*
   * TODO: Make a clever loop
   */
  for (int i=0; i<ijk_intersections_subproblems_indices_.size_array(); i++)
    {
      const int intersection_ijk_cell_index = ijk_intersections_subproblems_indices_[i];
      const int ijk_indices_i = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 0);
      const int ijk_indices_j = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 1);
      const int ijk_indices_k = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 2);
      for (int l=0; l<6; l++)
        {
          const int is_neighbour_pure_liquid = intersection_ijk_cell_->get_ijk_pure_face_neighbours(intersection_ijk_cell_index, l);
          if (is_neighbour_pure_liquid)
            {
              const int ii_f = neighbours_faces_i[l];
              const int jj_f = neighbours_faces_j[l];
              const int kk_f = neighbours_faces_k[l];
              i_pure_face_to_correct[nb_faces_to_correct] = (ijk_indices_i + ii_f);
              j_pure_face_to_correct[nb_faces_to_correct] = (ijk_indices_j + jj_f);
              k_pure_face_to_correct[nb_faces_to_correct] = (ijk_indices_k + kk_f);
              dir_pure_face_to_correct[nb_faces_to_correct] = faces_dir[l];
              nb_faces_to_correct++;
            }
        }
    }
  assert(nb_faces_to_correct== nb_faces_to_correct_from_ijk);
  const int convective_fluxes_size = convective_fluxes_.size();
  const int diffusive_fluxes_size = diffusive_fluxes_.size();
  if (convective_fluxes_size > 0)
    assert(convective_fluxes_size ==  nb_faces_to_correct);
  if (diffusive_fluxes_size > 0)
    assert(diffusive_fluxes_size == nb_faces_to_correct);
}

void Corrige_flux_FT_temperature_subresolution::sort_ijk_intersections_subproblems_indices_fluxes_by_k_layers(FixedVector<FixedVector<std::vector<ArrOfInt>,3>,2>& index_face_ij_flux_xyz,
                                                                                                              FixedVector<FixedVector<std::vector<ArrOfInt>,3>,2>& index_face_ij_flux_xyz_remaining_global,
                                                                                                              FixedVector<std::vector<ArrOfDouble>, 3>& flux_xyz,
                                                                                                              FixedVector<std::vector<ArrOfDouble>, 3>& flux_xyz_remaining_global,
                                                                                                              FixedVector<std::map<int, int>, 3>& flux_frontier_map,
                                                                                                              const DoubleVect& fluxes_subgrid,
                                                                                                              const int ini_index)

{

  const IJK_Field_double& indicator = ref_ijk_ft_->itfce().I();
  const int nb_i_layer = indicator.ni();
  const int nb_j_layer = indicator.nj();
  const int nb_k_layer = indicator.nk();

  // Make it parallel
  const IJK_Splitting& splitting_ns = ref_ijk_ft_->itfce().I().get_splitting();
  const int offset_i = splitting_ns.get_offset_local(0);
  const int offset_j = splitting_ns.get_offset_local(1);
  const int offset_k = splitting_ns.get_offset_local(2);

  const IJK_Grid_Geometry& geometry = splitting_ns.get_grid_geometry();
  const int nb_i_layer_tot = geometry.get_nb_elem_tot(0);
  const int nb_j_layer_tot = geometry.get_nb_elem_tot(1);
  const int nb_k_layer_tot = geometry.get_nb_elem_tot(2);

  const int seq = (Process::nproc() == 1);

  FixedVector<std::vector<ArrOfInt>,3>& dummy_int_array = index_face_ij_flux_xyz[0];
  FixedVector<std::vector<ArrOfDouble>,3>& dummy_double_array = flux_xyz;

  initialise_any_cell_neighbours_indices_to_correct_with_flux(index_face_ij_flux_xyz,
                                                              flux_xyz,
                                                              dummy_int_array,
                                                              dummy_double_array,
                                                              ini_index);
  if (!seq)
    {
      initialise_any_cell_neighbours_indices_to_correct_with_flux(index_face_ij_flux_xyz_remaining_global,
                                                                  flux_xyz_remaining_global,
                                                                  dummy_int_array,
                                                                  dummy_double_array,
                                                                  ini_index,
                                                                  1);
    }

  IntVect& i_pure_face_to_correct = ijk_faces_to_correct_[0];
  IntVect& j_pure_face_to_correct = ijk_faces_to_correct_[1];
  IntVect& k_pure_face_to_correct = ijk_faces_to_correct_[2];
  IntVect& dir_pure_face_to_correct = ijk_faces_to_correct_[3];
  const int nb_fluxes = ijk_faces_to_correct_[0].size();

  for (int i_flux=0; i_flux < nb_fluxes; i_flux++)
    {
      const int k = k_pure_face_to_correct[i_flux];
      const int dir = dir_pure_face_to_correct[i_flux];
      const double flux = fluxes_subgrid[i_flux];
      flux_xyz[dir][k].append_array(flux);
      const int i = i_pure_face_to_correct[i_flux];
      const int j = j_pure_face_to_correct[i_flux];
      if (!ini_index)
        {
          index_face_ij_flux_xyz[0][dir][k].append_array(i);
          index_face_ij_flux_xyz[1][dir][k].append_array(j);
        }
      switch(dir)
        {
        case DIRECTION_I:
          if (seq)
            {
              if (i == 0)
                {
                  if (!ini_index)
                    {
                      index_face_ij_flux_xyz[0][dir][k].append_array(nb_i_layer);
                      index_face_ij_flux_xyz[1][dir][k].append_array(j);
                    }
                  flux_xyz[dir][k].append_array(fluxes_subgrid[i_flux]);

                }
              if (i == nb_i_layer)
                {
                  if (!ini_index)
                    {
                      index_face_ij_flux_xyz[0][dir][k].append_array(0);
                      index_face_ij_flux_xyz[1][dir][k].append_array(j);
                    }
                  flux_xyz[dir][k].append_array(fluxes_subgrid[i_flux]);
                }
            }
          else
            {
              if (i == 0 || i == nb_i_layer)
                {
                  const int i_global = (i + offset_i); // % (ni_tot + 1);
                  const int j_global = (j + offset_j); // % nj_tot;
                  const int k_global = (k + offset_k); // % nj_tot;
                  const int index_flux = flux_xyz[dir][k].size_array() - 1;

                  if (!ini_index)
                    {
                      const int linear_index = get_linear_index_local(i, j, k, 0);
                      flux_frontier_map[dir][linear_index] = index_flux;
                      index_face_ij_flux_xyz_remaining_global[0][dir][k_global].append_array(i_global);
                      index_face_ij_flux_xyz_remaining_global[1][dir][k_global].append_array(j_global);
                    }
                  /*
                   * Overall periodicity
                   */
                  if (i_global == 0)
                    {
                      index_face_ij_flux_xyz_remaining_global[0][dir][k_global].append_array(nb_i_layer_tot);
                      index_face_ij_flux_xyz_remaining_global[1][dir][k_global].append_array(j_global);
                    }
                  if (i_global == nb_i_layer_tot)
                    {
                      index_face_ij_flux_xyz_remaining_global[0][dir][k_global].append_array(i_global);
                      index_face_ij_flux_xyz_remaining_global[1][dir][k_global].append_array(j_global);
                    }
                  flux_xyz_remaining_global[dir][k_global].append_array(flux);
                  if (i_global == 0 || i_global == nb_i_layer_tot)
                    flux_xyz_remaining_global[dir][k_global].append_array(flux);
                }
            }
          break;
        case DIRECTION_J:
          if (seq)
            {
              if (j == 0)
                {
                  if (!ini_index)
                    {
                      index_face_ij_flux_xyz[0][dir][k].append_array(i);
                      index_face_ij_flux_xyz[1][dir][k].append_array(nb_j_layer);
                    }
                  flux_xyz[dir][k].append_array(fluxes_subgrid[i_flux]);

                }
              if (j  == nb_j_layer)
                {
                  if (!ini_index)
                    {
                      index_face_ij_flux_xyz[0][dir][k].append_array(i);
                      index_face_ij_flux_xyz[1][dir][k].append_array(0);
                    }
                  flux_xyz[dir][k].append_array(fluxes_subgrid[i_flux]);
                }
            }
          else
            {
              if (j == 0 || j == nb_j_layer)
                {
                  const int i_global = (i + offset_i); // % (ni_tot + 1);
                  const int j_global = (j + offset_j); // % nj_tot;
                  const int k_global = (k + offset_k); // % nj_tot;
                  const int index_flux = flux_xyz[dir][k].size_array() - 1;

                  if (!ini_index)
                    {
                      const int linear_index = get_linear_index_local(i, j, k, 1);
                      flux_frontier_map[dir][linear_index] = index_flux;
                      index_face_ij_flux_xyz_remaining_global[0][dir][k_global].append_array(i_global);
                      index_face_ij_flux_xyz_remaining_global[1][dir][k_global].append_array(j_global);
                    }
                  /*
                   * Overall periodicity
                   */
                  if (j_global == 0)
                    {
                      index_face_ij_flux_xyz_remaining_global[0][dir][k_global].append_array(i_global);
                      index_face_ij_flux_xyz_remaining_global[1][dir][k_global].append_array(nb_j_layer_tot);
                    }
                  if (j_global == nb_j_layer_tot)
                    {
                      index_face_ij_flux_xyz_remaining_global[0][dir][k_global].append_array(i_global);
                      index_face_ij_flux_xyz_remaining_global[1][dir][k_global].append_array(0);
                    }

                  flux_xyz_remaining_global[dir][k_global].append_array(flux);
                  if (j_global == 0 || j_global == nb_j_layer_tot)
                    flux_xyz_remaining_global[dir][k_global].append_array(flux);
                }
            }
          break;
        case DIRECTION_K:
          if (seq)
            {
              if (k == 0)
                {
                  if (!ini_index)
                    {
                      index_face_ij_flux_xyz[0][dir][nb_k_layer_tot].append_array(i);
                      index_face_ij_flux_xyz[1][dir][nb_k_layer_tot].append_array(j);
                    }
                  flux_xyz[dir][nb_k_layer].append_array(fluxes_subgrid[i_flux]);
                }
              if (k == nb_k_layer)
                {
                  if (!ini_index)
                    {
                      index_face_ij_flux_xyz[0][dir][0].append_array(i);
                      index_face_ij_flux_xyz[1][dir][0].append_array(j);
                    }
                  flux_xyz[dir][0].append_array(fluxes_subgrid[i_flux]);
                }
            }
          else
            {
              if (k == 0 || k == nb_k_layer)
                {
                  const int i_global = (i + offset_i); // % (ni_tot + 1);
                  const int j_global = (j + offset_j); // % nj_tot;
                  const int k_global = (k + offset_k); // % nj_tot;
                  const int index_flux = flux_xyz[dir][k].size_array() - 1;

                  if (!ini_index)
                    {
                      const int linear_index = get_linear_index_local(i, j, k, 2);
                      flux_frontier_map[dir][linear_index] = index_flux;
                      index_face_ij_flux_xyz_remaining_global[0][dir][k_global].append_array(i_global);
                      index_face_ij_flux_xyz_remaining_global[1][dir][k_global].append_array(j_global);
                    }
                  /*
                   * Overall periodicity
                   */
                  if (k_global == 0)
                    {
                      index_face_ij_flux_xyz_remaining_global[0][dir][nb_k_layer_tot].append_array(i_global);
                      index_face_ij_flux_xyz_remaining_global[1][dir][nb_k_layer_tot].append_array(j_global);
                    }
                  if (k_global == nb_k_layer_tot)
                    {
                      index_face_ij_flux_xyz_remaining_global[0][dir][0].append_array(i_global);
                      index_face_ij_flux_xyz_remaining_global[1][dir][0].append_array(0);
                    }
                  flux_xyz_remaining_global[dir][k_global].append_array(flux);
                  if (k_global == 0)
                    flux_xyz_remaining_global[dir][nb_k_layer_tot].append_array(flux);
                  if (k_global == nb_k_layer_tot)
                    flux_xyz_remaining_global[dir][0].append_array(flux);
                }
            }
          break;
        }
    }

  receive_fluxes_from_frontier_on_procs(index_face_ij_flux_xyz_remaining_global,
                                        flux_xyz_remaining_global,
                                        ini_index);
  if (debug_)
    Cerr << "Fluxes have been received from procs" << finl;
  combine_fluxes_from_frontier_on_procs(index_face_ij_flux_xyz,
                                        index_face_ij_flux_xyz_remaining_global,
                                        flux_xyz,
                                        flux_xyz_remaining_global,
                                        flux_frontier_map,
                                        ini_index);
  if (debug_)
    Cerr << "Fluxes have been combined on procs" << finl;
}

int Corrige_flux_FT_temperature_subresolution::get_linear_index_local(const int& i, const int& j, const int& k, const int& dir)
{
  int nb_i_layer_tot = ref_ijk_ft_->itfce().I().ni();
  int nb_j_layer_tot = ref_ijk_ft_->itfce().I().nj();
  if (dir == 0)
    nb_i_layer_tot += 1;
  if (dir == 1)
    nb_j_layer_tot += 1;
  return (i + j * (nb_i_layer_tot) + k * (nb_i_layer_tot * nb_j_layer_tot));
}

int Corrige_flux_FT_temperature_subresolution::get_linear_index_global(const int& i, const int& j, const int& k, const int& dir)
{
  const IJK_Grid_Geometry& geometry = ref_ijk_ft_->itfce().I().get_splitting().get_grid_geometry();
  int nb_i_layer_tot = geometry.get_nb_elem_tot(0);
  int nb_j_layer_tot = geometry.get_nb_elem_tot(1);
  if (dir == 0)
    nb_i_layer_tot += 1;
  if (dir == 1)
    nb_j_layer_tot += 1;
  return (i + j * (nb_i_layer_tot) + k * (nb_i_layer_tot * nb_j_layer_tot));
}


void Corrige_flux_FT_temperature_subresolution::receive_fluxes_from_frontier_on_procs(FixedVector<FixedVector<std::vector<ArrOfInt>,3>,2>& index_face_ij_flux_xyz_remaining_global,
                                                                                      FixedVector<std::vector<ArrOfDouble>,3>& flux_xyz_remaining_global,
                                                                                      const int ini_index)
{
  const int nb_procs = Process::nproc();
  const int proc_num = Process::me();
  if (copy_fluxes_on_every_procs_)
    {
      if (nb_procs > 1)
        {
          FixedVector<std::vector<ArrOfInt>, 3> overall_numerotation_array;
          FixedVector<std::vector<ArrOfInt>, 3> start_indices_array;
          FixedVector<ArrOfInt, 3> size_array_global_array;

          int c, k, l;
          for (c=0; c<3; c++)
            {
              const int size_k_layers = (int) flux_xyz_remaining_global[c].size();
              overall_numerotation_array[c].resize(size_k_layers);
              start_indices_array[c].resize(size_k_layers);
              size_array_global_array[c] = ArrOfInt(size_k_layers);
              for (k=0; k<size_k_layers; k++)
                {
                  const int size_array = flux_xyz_remaining_global[c][k].size_array();
                  int size_array_global = size_array;
                  size_array_global = mp_sum(size_array_global);
                  size_array_global_array[c](k) = size_array_global;
                  overall_numerotation_array[c][k] = ArrOfInt(nb_procs);
                  start_indices_array[c][k] = ArrOfInt(nb_procs);
                  overall_numerotation_array[c][k](proc_num) = size_array;
                  mp_sum_for_each_item(overall_numerotation_array[c][k]);
                  for (l=1; l<overall_numerotation_array[c][k].size_array(); l++)
                    start_indices_array[c][k](l) = start_indices_array[c][k](l-1) + overall_numerotation_array[c][k](l-1);
                  // start_indices_array[c][k](l) = start_indices_array[c][k](l-1) + start_indices_array[c][k](l-1);

                  if (debug_)
                    {
                      Cerr << "Size array" << size_array << finl;
                      Cerr << "Size array global" << size_array_global << finl;
                      Cerr << "Overall_numerotation" << overall_numerotation_array[c][k](0) << "-" << overall_numerotation_array[c][k](1) << finl;
                    }

                  ArrOfDouble local_flux_values_tmp;
                  ArrOfDouble& global_flux_values_tmp = flux_xyz_remaining_global[c][k];

                  local_flux_values_tmp = global_flux_values_tmp;
                  global_flux_values_tmp.resize(size_array_global);
                  global_flux_values_tmp *= 0.;

                  for (l=0; l<local_flux_values_tmp.size_array(); l++)
                    global_flux_values_tmp(start_indices_array[c][k](proc_num) + l) = local_flux_values_tmp(l);

                  mp_sum_for_each_item(global_flux_values_tmp);

                }
            }

          if (!ini_index)
            {
              for (c=0; c<3; c++)
                {
                  const int size_k_layers = (int) index_face_ij_flux_xyz_remaining_global[0][c].size();
                  for (k=0; k<size_k_layers; k++)
                    {
                      const int size_array_global = size_array_global_array[c](k);
                      ArrOfInt& start_indices = start_indices_array[c][k];
                      ArrOfInt local_indices_i_tmp;
                      ArrOfInt local_indices_j_tmp;

                      ArrOfInt& global_indices_i_tmp = index_face_ij_flux_xyz_remaining_global[0][c][k];
                      ArrOfInt& global_indices_j_tmp = index_face_ij_flux_xyz_remaining_global[1][c][k];

                      local_indices_i_tmp = global_indices_i_tmp;
                      local_indices_j_tmp = global_indices_j_tmp;

                      global_indices_i_tmp.resize(size_array_global);
                      global_indices_j_tmp.resize(size_array_global);

                      global_indices_i_tmp *= 0;
                      global_indices_j_tmp *= 0;

                      for (l=0; l<local_indices_i_tmp.size_array(); l++)
                        {
                          global_indices_i_tmp(start_indices(proc_num) + l) = local_indices_i_tmp(l);
                          global_indices_j_tmp(start_indices(proc_num) + l) = local_indices_j_tmp(l);
                        }
                      mp_sum_for_each_item(global_indices_i_tmp);
                      mp_sum_for_each_item(global_indices_j_tmp);
                    }
                }
            }
        }
    }
}

void Corrige_flux_FT_temperature_subresolution::combine_fluxes_from_frontier_on_procs(FixedVector<FixedVector<std::vector<ArrOfInt>,3>,2>& index_face_ij_flux_xyz,
                                                                                      FixedVector<FixedVector<std::vector<ArrOfInt>,3>,2>& index_face_ij_flux_xyz_remaining_global,
                                                                                      FixedVector<std::vector<ArrOfDouble>,3>& flux_xyz,
                                                                                      FixedVector<std::vector<ArrOfDouble>,3>& flux_xyz_remaining_global,
                                                                                      FixedVector<std::map<int, int>, 3>& flux_frontier_map,
                                                                                      const int ini_index)
{

  const IJK_Field_double& indicator = ref_ijk_ft_->itfce().I();
  const int ni = indicator.ni();
  const int nj = indicator.nj();
  const int nk = indicator.nk();

  const IJK_Splitting& splitting_ns = ref_ijk_ft_->itfce().I().get_splitting();
  const int offset_i = splitting_ns.get_offset_local(0);
  const int offset_j = splitting_ns.get_offset_local(1);
  const int offset_k = splitting_ns.get_offset_local(2);

  FixedVector<std::map<int, int>, 3> multiple_flux_values_k;
  FixedVector<std::map<int, int>, 3> multiple_flux_values_count;
  FixedVector<std::map<int, double>, 3> multiple_flux_values_sum;

  FixedVector<std::map<int, int>, 3> flux_frontier_map_tmp = flux_frontier_map;

  for (int dir=0; dir<3; dir++)
    {
      const int size_k_layers = (int) index_face_ij_flux_xyz_remaining_global[0][dir].size();
      for (int k_global=0; k_global<size_k_layers; k_global++)
        {
          const int global_size_array = index_face_ij_flux_xyz_remaining_global[0][dir][k_global].size_array();
          for (int l=0; l<global_size_array; l++)
            {
              const int i_global = index_face_ij_flux_xyz_remaining_global[0][dir][k_global](l);
              const int j_global = index_face_ij_flux_xyz_remaining_global[1][dir][k_global](l);
              const double flux = flux_xyz_remaining_global[dir][k_global](l);
              const int i = i_global - offset_i;
              const int j = j_global - offset_j;
              const int k = k_global - offset_k;
              const int ni_max = (dir == 0 ? ni + 1 : ni);
              const int nj_max = (dir == 1 ? nj + 1 : nj);
              const int nk_max = (dir == 2 ? nk + 1 : nk);
              if ((0 <= i && i < ni_max) && (0 <= j && j < nj_max) && (0 <= k && k < nk_max))
                {
                  const int linear_local_index = get_linear_index_local(i, j, k, dir);
                  const int non_zero_value_local = (int) flux_frontier_map_tmp[dir].count(linear_local_index);
                  if (!non_zero_value_local)
                    {
                      // Add flux at the end of the list if new
                      flux_xyz[dir][k].append_array(flux);
                      if (!ini_index)
                        {
                          index_face_ij_flux_xyz[0][dir][k].append_array(i);
                          index_face_ij_flux_xyz[1][dir][k].append_array(j);
                        }
                      const int local_size_array = flux_xyz[dir][k].size_array() - 1;
                      flux_frontier_map_tmp[dir][linear_local_index] = local_size_array;
                      multiple_flux_values_count[dir][linear_local_index] = 1;
                      multiple_flux_values_sum[dir][linear_local_index] = flux;
                      multiple_flux_values_k[dir][linear_local_index] = k;
                    }
                  else
                    {
                      const int non_zero_local_count = (int) multiple_flux_values_count[dir].count(linear_local_index);
                      // const int array_index = (int) flux_frontier_map_[dir][linear_local_index];
                      if (!non_zero_local_count)
                        {
                          /*
                           * The processor will see its value twice
                           */
                          //                      multiple_flux_values_count[dir][linear_local_index] = 1;
                          //                      multiple_flux_values_sum[dir][linear_local_index] = (*(fluxes[dir]))[k](array_index);
                          multiple_flux_values_count[dir][linear_local_index] = 0;
                          multiple_flux_values_sum[dir][linear_local_index] = 0.;
                          multiple_flux_values_k[dir][linear_local_index] = k;
                        }
                      multiple_flux_values_count[dir][linear_local_index] += 1;
                      multiple_flux_values_sum[dir][linear_local_index] += flux;
                    }
                }
            }
        }
      for(std::map<int,double>::iterator it=multiple_flux_values_sum[dir].begin(); it!=multiple_flux_values_sum[dir].end(); ++it)
        {
          const int key = it->first;
          const double val_flux_sum = it->second;
          const double count_val = (double) multiple_flux_values_count[dir][key];
          const int array_index = (int) flux_frontier_map_tmp[dir][key];
          const int k_local = (int) multiple_flux_values_k[dir][key];
          flux_xyz[dir][k_local](array_index) =  val_flux_sum / count_val;
        }
    }
  if (ini_index)
    flux_frontier_map = flux_frontier_map_tmp;
}

void Corrige_flux_FT_temperature_subresolution::initialise_any_cell_neighbours_indices_to_correct_on_processors(FixedVector<FixedVector<std::vector<std::vector<ArrOfInt>>,3>,2>& index_face_ij_flux_xyz,
                                                                                                                FixedVector<std::vector<std::vector<ArrOfDouble>>,3>& flux_xyz,
                                                                                                                const int ini_index)
{
  const int nb_proc = Process::nproc();
  int c;
  if(!ini_index)
    {
      for(int l=0; l<2; l++)
        for (c=0; c<3; c++)
          index_face_ij_flux_xyz[l][c].resize(nb_proc);
    }
  for (c=0; c<3; c++)
    flux_xyz[c].resize(nb_proc);
}

void Corrige_flux_FT_temperature_subresolution::redistribute_indices_fluxes_by_k_layers(FixedVector<FixedVector<std::vector<ArrOfInt>,3>,2>& index_face_i_flux_x,
                                                                                        FixedVector<FixedVector<std::vector<ArrOfInt>,3>,2>& index_face_i_flux_x_remaining_global,
                                                                                        FixedVector<std::vector<ArrOfDouble>,3>& flux_xyz,
                                                                                        FixedVector<std::vector<ArrOfDouble>,3>& flux_xyz_remaining_global,
                                                                                        const int ini_index)
{

}

void Corrige_flux_FT_temperature_subresolution::store_cell_faces_corrected(FixedVector<IJK_Field_int,3>& cell_faces_corrected_bool,
                                                                           FixedVector<IJK_Field_double,3>& cell_faces_corrected_convective,
                                                                           FixedVector<IJK_Field_double,3>& cell_faces_corrected_diffusive)
{
  int counter = 0;
  if (!convection_negligible_)
    {
      store_any_cell_faces_corrected(cell_faces_corrected_bool,
                                     cell_faces_corrected_convective,
                                     convective_fluxes_,
                                     convective_diffusive_flux_xyz_sorted_[0],
                                     counter);
      counter ++;
    }
  if (!diffusion_negligible_)
    store_any_cell_faces_corrected(cell_faces_corrected_bool,
                                   cell_faces_corrected_diffusive,
                                   diffusive_fluxes_,
                                   convective_diffusive_flux_xyz_sorted_[1],
                                   counter);
}

void Corrige_flux_FT_temperature_subresolution::store_any_cell_faces_corrected(FixedVector<IJK_Field_int,3>& cell_faces_corrected_bool,
                                                                               FixedVector<IJK_Field_double,3>& cell_faces_corrected,
                                                                               const DoubleVect& fluxes,
                                                                               FixedVector<std::vector<ArrOfDouble>,3>& flux_xyz,
                                                                               const int counter)
{
  for (int c=0; c<3; c++)
    {
      cell_faces_corrected_bool[c].data() = 0.;
      cell_faces_corrected[c].data() = 0.;
    }
  if (Process::nproc() == 1)
    {
      IntVect& i_pure_face_to_correct = ijk_faces_to_correct_[0];
      IntVect& j_pure_face_to_correct = ijk_faces_to_correct_[1];
      IntVect& k_pure_face_to_correct = ijk_faces_to_correct_[2];
      IntVect& dir_pure_face_to_correct = ijk_faces_to_correct_[3];
      const int nb_fluxes = ijk_faces_to_correct_[0].size();
      int i,j,k;
      for (int i_flux=0; i_flux < nb_fluxes; i_flux++)
        {
          i = i_pure_face_to_correct[i_flux];
          j = j_pure_face_to_correct[i_flux];
          k = k_pure_face_to_correct[i_flux];
          const int dir = dir_pure_face_to_correct[i_flux];
          if (!counter)
            cell_faces_corrected_bool[dir](i,j,k) += 1;
          cell_faces_corrected[dir](i,j,k) += fluxes[i_flux];
        }
    }
  else
    {
      const IJK_Field_double& indicator = ref_ijk_ft_->itfce().I();
      const int ni = indicator.ni();
      const int nj = indicator.nj();
      const int nk = indicator.nk();

      const int size_k = (int) index_face_ij_flux_xyz_sorted_[0][0].size();
      for (int dir=0; dir<3; dir++)
        for (int k=0; k<size_k; k++)
          {
            const int size_array = index_face_ij_flux_xyz_sorted_[0][dir][k].size_array();
            for (int l=0; l<size_array; l++)
              {
                const int i = index_face_ij_flux_xyz_sorted_[0][dir][k](l);
                const int j = index_face_ij_flux_xyz_sorted_[1][dir][k](l);
                const double flux = flux_xyz[dir][k](l);
                /*
                 * Offset for verification (not at the exact places of the fluxes)!
                 */
                const int i_display = (i >= ni ? i-1: i);
                const int j_display = (j >= nj ? j-1: j);
                const int k_display = (k >= nk ? k-1: k);
                if (!counter)
                  cell_faces_corrected_bool[dir](i_display, j_display, k_display) += 1;
                cell_faces_corrected[dir](i_display,j_display,k_display) += flux;
              }
          }
    }
}

void Corrige_flux_FT_temperature_subresolution::sort_ijk_intersections_subproblems_indices_by_k_layers()
{
  if (!convection_negligible_)
    {
      Cerr << "Sort the thermal convective fluxes" << finl;
      sort_ijk_intersections_subproblems_indices_fluxes_by_k_layers(index_face_ij_flux_xyz_sorted_,
                                                                    index_face_ij_flux_xyz_remaining_global_sorted_,
                                                                    convective_diffusive_flux_xyz_sorted_[0],
                                                                    convective_diffusive_flux_xyz_remaining_global_sorted_[0],
                                                                    flux_frontier_map_,
                                                                    convective_fluxes_,
                                                                    flux_init_);
      Cerr << "Thermal Sub-resolutions convective fluxes are now sorted" << finl;
      flux_init_ = 1;
    }
  if (!diffusion_negligible_)
    {
      Cerr << "Sort the thermal diffusive fluxes" << finl;
      sort_ijk_intersections_subproblems_indices_fluxes_by_k_layers(index_face_ij_flux_xyz_sorted_,
                                                                    index_face_ij_flux_xyz_remaining_global_sorted_,
                                                                    convective_diffusive_flux_xyz_sorted_[1],
                                                                    convective_diffusive_flux_xyz_remaining_global_sorted_[1],
                                                                    flux_frontier_map_,
                                                                    diffusive_fluxes_,
                                                                    flux_init_);
      Cerr << "Thermal Sub-resolutions diffusive fluxes are now sorted" << finl;
      flux_init_ = 1;
    }
  flux_init_ = 0;
}

void Corrige_flux_FT_temperature_subresolution::corrige_flux_faceIJ(IJK_Field_local_double *const flux,
                                                                    const int k_layer,
                                                                    const int dir)
{
  if (convective_flux_correction_)
    {
      if(!use_reachable_fluxes_)
        corrige_flux_faceIJ_any_flux(flux,
                                     index_face_ij_flux_xyz_sorted_,
                                     convective_diffusive_flux_xyz_sorted_[0],
                                     k_layer,
                                     dir);
      else
        corrige_flux_faceIJ_any_flux(flux,
                                     index_face_ij_flux_xyz_neighbours_min_max_faces_sorted_,
                                     convective_diffusive_flux_xyz_min_max_faces_sorted_[0],
                                     k_layer,
                                     dir);
    }
}

void Corrige_flux_FT_temperature_subresolution::corrige_flux_diff_faceIJ(IJK_Field_local_double *const flux,
                                                                         const int k_layer,
                                                                         const int dir)
{
  if (diffusive_flux_correction_)
    {
      if(!use_reachable_fluxes_)
        corrige_flux_faceIJ_any_flux(flux,
                                     index_face_ij_flux_xyz_sorted_,
                                     convective_diffusive_flux_xyz_sorted_[1],
                                     k_layer,
                                     dir);
      else
        corrige_flux_faceIJ_any_flux(flux,
                                     index_face_ij_flux_xyz_neighbours_min_max_faces_sorted_,
                                     convective_diffusive_flux_xyz_min_max_faces_sorted_[1],
                                     k_layer,
                                     dir);
    }
}

void Corrige_flux_FT_temperature_subresolution::corrige_flux_faceIJ_any_flux(IJK_Field_local_double *const flux,
                                                                             FixedVector<FixedVector<std::vector<ArrOfInt>,3>,2>& index_face_ij_flux_xyz_sorted,
                                                                             FixedVector<std::vector<ArrOfDouble>,3>& subgrid_fluxes_xyz,
                                                                             const int k_layer,
                                                                             const int dir)
{
  const int nb_fluxes = subgrid_fluxes_xyz[dir][k_layer].size_array();
  for (int i_flux=0; i_flux<nb_fluxes; i_flux++)
    {
      const int i=index_face_ij_flux_xyz_sorted[0][dir][k_layer][i_flux];
      const int j=index_face_ij_flux_xyz_sorted[1][dir][k_layer][i_flux];
      const double flux_ij = subgrid_fluxes_xyz[dir][k_layer][i_flux];
      (*flux)(i,j,0) = flux_ij;
    }
}

void Corrige_flux_FT_temperature_subresolution::check_pure_fluxes_duplicates(const DoubleVect& fluxes,
                                                                             DoubleVect& fluxes_unique,
                                                                             IntVect& pure_face_unique,
                                                                             const int known_unique)
{
  // FixedVector<DoubleVect, 3>& pure_face_to_correct = intersection_ijk_cell_->get_set_ijk_pure_face_to_correct();
  FixedVector<IntVect, 4>& pure_face_to_correct = ijk_faces_to_correct_;
  const int nb_fluxes = fluxes.size();
  fluxes_unique.set_smart_resize(1);
  fluxes_unique.reset();
  int i;
  if (known_unique)
    {
      const int size_face_unique = pure_face_unique.size();
      for (i=0; i<size_face_unique; i++)
        {
          fluxes_unique.append_array(fluxes(pure_face_unique(i)));
        }
    }
  else
    {
      DoubleVect shared_face;
      shared_face.set_smart_resize(1);
      shared_face.reset();
      for (i=0; i<nb_fluxes; i++)
        {
          const int i_f = pure_face_to_correct[0](i);
          const int j_f = pure_face_to_correct[1](i);
          const int k_f = pure_face_to_correct[2](i);
          for (int j=i; j<nb_fluxes; j++)
            if (i != j)
              {
                const int i_ff = pure_face_to_correct[0](j);
                const int j_ff = pure_face_to_correct[1](j);
                const int k_ff = pure_face_to_correct[2](j);
                if ((i_f==i_ff) && (j_f==j_ff) && (k_f==k_ff))
                  {
                    if (shared_face.size() == 0)
                      shared_face.append_array(i);
                    shared_face.append_array(j);
                  }
              }
          /*
           * Take the closest portion
           */
          const int nb_duplicates = shared_face.size();
          int min_dist_index = 0;
          double min_dist = dist_[i];
          for (int j=1; j<nb_duplicates; j++)
            if (min_dist > dist_[j])
              {
                min_dist = dist_[j];
                min_dist_index = j;
              }
          pure_face_unique.append_array(min_dist_index);
          fluxes_unique.append_array(fluxes(min_dist_index));
        }
    }
}
