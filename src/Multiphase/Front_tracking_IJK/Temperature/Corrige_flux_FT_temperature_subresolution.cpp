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
  index_face_i_flux_x_sorted_.clear();
  index_face_j_flux_x_sorted_.clear();

  index_face_i_flux_y_sorted_.clear();
  index_face_j_flux_y_sorted_.clear();

  index_face_i_flux_z_sorted_.clear();
  index_face_j_flux_z_sorted_.clear();

  convective_flux_x_sorted_.clear();
  convective_flux_y_sorted_.clear();
  convective_flux_z_sorted_.clear();

  diffusive_flux_x_sorted_.clear();
  diffusive_flux_y_sorted_.clear();
  diffusive_flux_z_sorted_.clear();

  index_face_i_flux_x_neighbours_diag_faces_sorted_.clear();
  index_face_j_flux_x_neighbours_diag_faces_sorted_.clear();
  index_face_i_flux_y_neighbours_diag_faces_sorted_.clear();
  index_face_j_flux_y_neighbours_diag_faces_sorted_.clear();
  index_face_i_flux_z_neighbours_diag_faces_sorted_.clear();
  index_face_j_flux_z_neighbours_diag_faces_sorted_.clear();

  index_face_i_flux_x_neighbours_all_faces_sorted_.clear();
  index_face_j_flux_x_neighbours_all_faces_sorted_.clear();
  index_face_i_flux_y_neighbours_all_faces_sorted_.clear();
  index_face_j_flux_y_neighbours_all_faces_sorted_.clear();
  index_face_i_flux_z_neighbours_all_faces_sorted_.clear();
  index_face_j_flux_z_neighbours_all_faces_sorted_.clear();

  index_face_i_flux_x_neighbours_min_max_faces_sorted_.clear();
  index_face_j_flux_x_neighbours_min_max_faces_sorted_.clear();
  index_face_i_flux_y_neighbours_min_max_faces_sorted_.clear();
  index_face_j_flux_y_neighbours_min_max_faces_sorted_.clear();
  index_face_i_flux_z_neighbours_min_max_faces_sorted_.clear();
  index_face_j_flux_z_neighbours_min_max_faces_sorted_.clear();

  convective_flux_x_min_max_faces_sorted_.clear();
  convective_flux_y_min_max_faces_sorted_.clear();
  convective_flux_z_min_max_faces_sorted_.clear();

  diffusive_flux_x_min_max_faces_sorted_.clear();
  diffusive_flux_y_min_max_faces_sorted_.clear();
  diffusive_flux_z_min_max_faces_sorted_.clear();
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

  // On commence par calculer les temperatures aux faces mouillées
  intersection_ijk_cell_->update_interpolations_cell_centres_on_interface();
  /*
   * TODO update with face cell centres positions
   */
  if (!convection_negligible_ || !diffusion_negligible_)
    intersection_ijk_cell_->update_interpolations_cell_faces_on_interface();

  Cerr << "The intersections have been updated" << finl;
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
      has_checked_consistency_ = (nb_diph==nb_subproblems);
      assert(has_checked_consistency_);
      ijk_intersections_subproblems_indices_.reset();
      ijk_intersections_subproblems_indices_.resize(nb_subproblems);
      int index_i_problem = 0;
      int index_j_problem = 0;
      int index_k_problem = 0;
      int problem_index;
      for (problem_index=0; problem_index<nb_subproblems; problem_index++)
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
  has_checked_consistency_=false;
  intersection_ijk_cell_->set_pas_a_jour();
}


void Corrige_flux_FT_temperature_subresolution::compute_temperature_cell_centre(IJK_Field_double& temperature) const
{
  /*
   * For each subproblem fill the right interfacial_cell
   */
  const DoubleTab dist_interf = intersection_ijk_cell_->dist_interf();
  const double min_temperature = thermal_subproblems_->get_min_temperature_domain_ends();
  const double max_temperature = thermal_subproblems_->get_max_temperature_domain_ends();
  for (int i=0; i<ijk_intersections_subproblems_indices_.size_array(); i++)
    {
      const int intersection_ijk_cell_index = ijk_intersections_subproblems_indices_[i];
      const double dist = dist_interf(intersection_ijk_cell_index, 0);
      const double dist_sub_res = thermal_subproblems_->get_dist_cell_interface(i);

      double temperature_ghost = 0.;
      if (distance_cell_faces_from_lrs_)
        temperature_ghost = thermal_subproblems_->get_temperature_profile_at_point(i, dist_sub_res);
      else
        temperature_ghost = thermal_subproblems_->get_temperature_profile_at_point(i, dist);

      const int ijk_indices_i = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 0);
      const int ijk_indices_j = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 1);
      const int ijk_indices_k = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 2);

      if (debug_)
        {
          Cerr << "Time-step: " << ref_ijk_ft_->get_tstep() << "--" << ref_ijk_ft_->get_timestep() << " s" << finl;
          Cerr << "Distance at cell : " << intersection_ijk_cell_index <<
               ", subproblem " << i << "."<<
               " -- (" << ijk_indices_i << ", " << ijk_indices_j << ", " << ijk_indices_k << ")" << finl;
          Cerr << "Distance from intersection_ijk_cell: " << dist << finl;
          Cerr << "Distance from sub-resolution: " << dist_sub_res << finl;
          Vecteur3 bary_facet_debug = thermal_subproblems_->get_bary_facet(i);
          Cerr << "Facet barycentre: " << bary_facet_debug[0] << ";"
               << bary_facet_debug[1] << ";"
               << bary_facet_debug[2] << finl;
        }

      const IJK_Field_double& indicator = ref_ijk_ft_->itfce().I();
      const double indic = indicator(ijk_indices_i, ijk_indices_j, ijk_indices_k);
      if (temperature_ghost < min_temperature && indic > 0.5)
        Cerr << "Ghost temperature: " << temperature_ghost << " is lower than the minimum temperature:" << min_temperature << finl;
      if (temperature_ghost > max_temperature && indic > 0.5)
        Cerr << "Ghost temperature: " << temperature_ghost << " is higher than the maximum temperature:" << max_temperature << finl;

      temperature(ijk_indices_i, ijk_indices_j, ijk_indices_k) = temperature_ghost;
    }
}

void Corrige_flux_FT_temperature_subresolution::compute_temperature_cell_centre_neighbours(IJK_Field_double& temperature_neighbours,
                                                                                           IJK_Field_int& neighbours_weighting,
                                                                                           IJK_Field_double& neighbours_weighting_colinearity) const
{
  if (distance_cell_faces_from_lrs_ && find_temperature_cell_neighbours_)
    {

      const int ni = neighbours_weighting.ni();
      const int nj = neighbours_weighting.nj();
      const int nk = neighbours_weighting.nk();

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
      int index_i_perio, index_j_perio, index_k_perio;
      int m,l,n;
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
                    if (pure_neighbours_to_correct[l][m][n]) // (l!=0 || m!=0 || n!=0) &&
                      {
                        const double dist_sub_res = pure_neighbours_corrected_distance[l][m][n];
                        const double temperature_ghost = thermal_subproblems_->get_temperature_profile_at_point(i, dist_sub_res);
                        index_i_neighbour = index_i_problem + ((pure_neighbours_corrected_sign[0]) ?  l * (-1) : l);
                        index_j_neighbour = index_j_problem + ((pure_neighbours_corrected_sign[1]) ?  m * (-1) : m);
                        index_k_neighbour = index_k_problem + ((pure_neighbours_corrected_sign[2]) ?  n * (-1) : n);
                        index_i_perio = index_i_neighbour % ni;
                        index_j_perio = index_j_neighbour % nj;
                        index_k_perio = index_k_neighbour % nk;
                        if (neighbours_colinearity_weighting_)
                          {
                            // neighbours_weighting_colinearity(index_i_neighbour, index_j_neighbour, index_k_neighbour) += pure_neighbours_corrected_colinearity[l][m][n];
                            // temperature_neighbours(index_i_neighbour, index_j_neighbour, index_k_neighbour) += temperature_ghost * pure_neighbours_corrected_colinearity[l][m][n];
                            neighbours_weighting_colinearity(index_i_perio, index_j_perio, index_k_perio) += pure_neighbours_corrected_colinearity[l][m][n];
                            temperature_neighbours(index_i_perio, index_j_perio, index_k_perio) += temperature_ghost * pure_neighbours_corrected_colinearity[l][m][n];
                          }
                        else
                          temperature_neighbours(index_i_perio, index_j_perio, index_k_perio) += temperature_ghost;
                        // temperature_neighbours(index_i_neighbour, index_j_neighbour, index_k_neighbour) += temperature_ghost;
                        // neighbours_weighting(index_i_neighbour, index_j_neighbour, index_k_neighbour) += 1;
                        neighbours_weighting(index_i_perio, index_j_perio, index_k_perio) += 1;
                      }
            }
        }
      if (neighbours_colinearity_weighting_)
        neighbours_weighting_colinearity.echange_espace_virtuel(neighbours_weighting_colinearity.ghost());
      neighbours_weighting.echange_espace_virtuel(neighbours_weighting.ghost());
      temperature_neighbours.echange_espace_virtuel(temperature_neighbours.ghost());
    }
}

void Corrige_flux_FT_temperature_subresolution::initialise_any_cell_neighbours_indices_to_correct(std::vector<ArrOfInt>& index_face_i_flux_x_faces_sorted,
                                                                                                  std::vector<ArrOfInt>& index_face_j_flux_x_faces_sorted,
                                                                                                  std::vector<ArrOfInt>& index_face_i_flux_y_faces_sorted,
                                                                                                  std::vector<ArrOfInt>& index_face_j_flux_y_faces_sorted,
                                                                                                  std::vector<ArrOfInt>& index_face_i_flux_z_faces_sorted,
                                                                                                  std::vector<ArrOfInt>& index_face_j_flux_z_faces_sorted)
{
  const int nb_k_layer = ref_ijk_ft_->itfce().I().nk();

  index_face_i_flux_x_faces_sorted.resize(nb_k_layer);
  index_face_j_flux_x_faces_sorted.resize(nb_k_layer);
  index_face_i_flux_y_faces_sorted.resize(nb_k_layer);
  index_face_j_flux_y_faces_sorted.resize(nb_k_layer);
  index_face_i_flux_z_faces_sorted.resize(nb_k_layer + 1);
  index_face_j_flux_z_faces_sorted.resize(nb_k_layer + 1);

  FixedVector<std::vector<ArrOfInt>*,3> index_face_i_sorted;
  index_face_i_sorted[0] = &index_face_i_flux_x_faces_sorted;
  index_face_i_sorted[1] = &index_face_i_flux_y_faces_sorted;
  index_face_i_sorted[2] = &index_face_i_flux_z_faces_sorted;

  FixedVector<std::vector<ArrOfInt>*,3> index_face_j_sorted;
  index_face_j_sorted[0] = &index_face_j_flux_x_faces_sorted;
  index_face_j_sorted[1] = &index_face_j_flux_y_faces_sorted;
  index_face_j_sorted[2] = &index_face_j_flux_z_faces_sorted;

  for (int dir=0; dir<3; dir++)
    for (int k_layer=0; k_layer<nb_k_layer+1; k_layer++)
      {
        if ((dir==DIRECTION_I || dir==DIRECTION_J) && k_layer==nb_k_layer)
          break;
        (*(index_face_i_sorted[dir]))[k_layer].reset();
        (*(index_face_j_sorted[dir]))[k_layer].reset();
        (*(index_face_i_sorted[dir]))[k_layer].set_smart_resize(1);
        (*(index_face_j_sorted[dir]))[k_layer].set_smart_resize(1);
      }

}

void Corrige_flux_FT_temperature_subresolution::initialise_any_cell_neighbours_indices_to_correct_with_flux(std::vector<ArrOfInt>& index_face_i_flux_x_faces_sorted,
                                                                                                            std::vector<ArrOfInt>& index_face_j_flux_x_faces_sorted,
                                                                                                            std::vector<ArrOfInt>& index_face_i_flux_y_faces_sorted,
                                                                                                            std::vector<ArrOfInt>& index_face_j_flux_y_faces_sorted,
                                                                                                            std::vector<ArrOfInt>& index_face_i_flux_z_faces_sorted,
                                                                                                            std::vector<ArrOfInt>& index_face_j_flux_z_faces_sorted,
                                                                                                            std::vector<ArrOfDouble>& flux_x,
                                                                                                            std::vector<ArrOfDouble>& flux_y,
                                                                                                            std::vector<ArrOfDouble>& flux_z,
                                                                                                            const bool& ini_index)
{
  if (!ini_index)
    initialise_any_cell_neighbours_indices_to_correct(index_face_i_flux_x_neighbours_min_max_faces_sorted_,
                                                      index_face_j_flux_x_neighbours_min_max_faces_sorted_,
                                                      index_face_i_flux_y_neighbours_min_max_faces_sorted_,
                                                      index_face_j_flux_y_neighbours_min_max_faces_sorted_,
                                                      index_face_i_flux_z_neighbours_min_max_faces_sorted_,
                                                      index_face_j_flux_z_neighbours_min_max_faces_sorted_);

  const int nb_k_layer = ref_ijk_ft_->itfce().I().nk();

  flux_x.resize(nb_k_layer);
  flux_y.resize(nb_k_layer);
  flux_z.resize(nb_k_layer + 1);

  FixedVector<std::vector<ArrOfDouble>*,3> fluxes;
  fluxes[0] = &flux_x;
  fluxes[1] = &flux_y;
  fluxes[2] = &flux_z;

  for (int dir=0; dir<3; dir++)
    for (int k_layer=0; k_layer<nb_k_layer+1; k_layer++)
      {
        if ((dir==DIRECTION_I || dir==DIRECTION_J) && k_layer==nb_k_layer)
          break;
        (*(fluxes[dir]))[k_layer].reset();
        (*(fluxes[dir]))[k_layer].set_smart_resize(1);
      }
}

void Corrige_flux_FT_temperature_subresolution::initialise_cell_neighbours_indices_to_correct()
{
  if (distance_cell_faces_from_lrs_ && find_temperature_cell_neighbours_ && find_cell_neighbours_for_fluxes_spherical_correction_)
    {
      initialise_any_cell_neighbours_indices_to_correct(index_face_i_flux_x_neighbours_diag_faces_sorted_,
                                                        index_face_j_flux_x_neighbours_diag_faces_sorted_,
                                                        index_face_i_flux_y_neighbours_diag_faces_sorted_,
                                                        index_face_j_flux_y_neighbours_diag_faces_sorted_,
                                                        index_face_i_flux_z_neighbours_diag_faces_sorted_,
                                                        index_face_j_flux_z_neighbours_diag_faces_sorted_);
    }
  if (distance_cell_faces_from_lrs_ && use_reachable_fluxes_)
    {
      if (!convection_negligible_)
        {
          Cerr << "Sort the thermal min-max reachable convective fluxes" << finl;
          initialise_any_cell_neighbours_indices_to_correct_with_flux(index_face_i_flux_x_neighbours_min_max_faces_sorted_,
                                                                      index_face_j_flux_x_neighbours_min_max_faces_sorted_,
                                                                      index_face_i_flux_y_neighbours_min_max_faces_sorted_,
                                                                      index_face_j_flux_y_neighbours_min_max_faces_sorted_,
                                                                      index_face_i_flux_z_neighbours_min_max_faces_sorted_,
                                                                      index_face_j_flux_z_neighbours_min_max_faces_sorted_,
                                                                      convective_flux_x_min_max_faces_sorted_,
                                                                      convective_flux_y_min_max_faces_sorted_,
                                                                      convective_flux_z_min_max_faces_sorted_,
                                                                      flux_init_);
          Cerr << "Thermal Sub-resolutions min-max reachable convective fluxes are now sorted" << finl;
          flux_init_ = 1;
        }
      if (!diffusion_negligible_)
        {
          Cerr << "Sort the thermal diffusive min-max reachable fluxes" << finl;
          initialise_any_cell_neighbours_indices_to_correct_with_flux(index_face_i_flux_x_neighbours_min_max_faces_sorted_,
                                                                      index_face_j_flux_x_neighbours_min_max_faces_sorted_,
                                                                      index_face_i_flux_y_neighbours_min_max_faces_sorted_,
                                                                      index_face_j_flux_y_neighbours_min_max_faces_sorted_,
                                                                      index_face_i_flux_z_neighbours_min_max_faces_sorted_,
                                                                      index_face_j_flux_z_neighbours_min_max_faces_sorted_,
                                                                      diffusive_flux_x_min_max_faces_sorted_,
                                                                      diffusive_flux_y_min_max_faces_sorted_,
                                                                      diffusive_flux_z_min_max_faces_sorted_,
                                                                      flux_init_);
          Cerr << "Thermal Sub-resolutions min-max reachable diffusive fluxes are now sorted" << finl;
          flux_init_ = 1;
        }
      flux_init_ = 0;

    }
}


void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_faces_indices_for_spherical_correction(const int& n_iter_distance)
{
  /*
   * TODO: What to do with duplicates ?
   */
  if (distance_cell_faces_from_lrs_ && find_cell_neighbours_for_fluxes_spherical_correction_)
    {
      for (int c=0; c<3; c++)
        (*cell_faces_neighbours_corrected_bool_)[c].data() = 0;
      (*cell_faces_neighbours_corrected_bool_).echange_espace_virtuel();
      const int nb_i_layer = ref_ijk_ft_->itfce().I().ni();
      const int nb_j_layer = ref_ijk_ft_->itfce().I().nj();
      const int nb_k_layer = ref_ijk_ft_->itfce().I().nk();

      FixedVector<std::vector<ArrOfInt>*,3> index_face_i_sorted;
      index_face_i_sorted[0] = &index_face_i_flux_x_neighbours_diag_faces_sorted_;
      index_face_i_sorted[1] = &index_face_i_flux_y_neighbours_diag_faces_sorted_;
      index_face_i_sorted[2] = &index_face_i_flux_z_neighbours_diag_faces_sorted_;

      FixedVector<std::vector<ArrOfInt>*,3> index_face_j_sorted;
      index_face_j_sorted[0] = &index_face_j_flux_x_neighbours_diag_faces_sorted_;
      index_face_j_sorted[1] = &index_face_j_flux_y_neighbours_diag_faces_sorted_;
      index_face_j_sorted[2] = &index_face_j_flux_z_neighbours_diag_faces_sorted_;

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
              //if (pure_neighbours_to_correct[l][m][n] && (l!=0 || m!=0 || n!=0))
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
                                  (*(index_face_i_sorted[0]))[index_k_neighbour].append_array(i_offset_face);
                                  (*(index_face_j_sorted[0]))[index_k_neighbour].append_array(index_j_neighbour);
                                }
                              (*cell_faces_neighbours_corrected_bool_)[0](i_offset_face, index_j_neighbour, index_k_neighbour) += 1;
                              if (i_offset_face == 0)
                                {
                                  if(!(*cell_faces_neighbours_corrected_bool_)[0](nb_i_layer, index_j_neighbour, index_k_neighbour))
                                    {
                                      (*(index_face_i_sorted[0]))[index_k_neighbour].append_array(nb_i_layer);
                                      (*(index_face_j_sorted[0]))[index_k_neighbour].append_array(index_j_neighbour);
                                    }
                                  (*cell_faces_neighbours_corrected_bool_)[0](nb_i_layer, index_j_neighbour, index_k_neighbour) += 1;
                                }
                              if (i_offset_face == nb_i_layer)
                                {
                                  if (!(*cell_faces_neighbours_corrected_bool_)[0](0, index_j_neighbour, index_k_neighbour))
                                    {
                                      (*(index_face_i_sorted[0]))[index_k_neighbour].append_array(0);
                                      (*(index_face_j_sorted[0]))[index_k_neighbour].append_array(index_j_neighbour);
                                    }
                                  (*cell_faces_neighbours_corrected_bool_)[0](0, index_j_neighbour, index_k_neighbour) += 1;
                                }
                            }
                          if (m != 0)
                            {
                              if (!(*cell_faces_neighbours_corrected_bool_)[1](index_i_neighbour, j_offset_face, index_k_neighbour))
                                {
                                  (*(index_face_i_sorted[1]))[index_k_neighbour].append_array(index_i_neighbour);
                                  (*(index_face_j_sorted[1]))[index_k_neighbour].append_array(j_offset_face);
                                }
                              (*cell_faces_neighbours_corrected_bool_)[1](index_i_neighbour, j_offset_face, index_k_neighbour) += 1;
                              if (j_offset_face == 0)
                                {
                                  if (!(*cell_faces_neighbours_corrected_bool_)[1](index_i_neighbour, nb_j_layer, index_k_neighbour))
                                    {
                                      (*(index_face_i_sorted[1]))[index_k_neighbour].append_array(index_i_neighbour);
                                      (*(index_face_j_sorted[1]))[index_k_neighbour].append_array(nb_j_layer);
                                    }
                                  (*cell_faces_neighbours_corrected_bool_)[1](index_i_neighbour, nb_j_layer, index_k_neighbour) += 1;
                                }
                              if (i_offset_face == nb_j_layer)
                                {
                                  if (!(*cell_faces_neighbours_corrected_bool_)[1](index_i_neighbour, 0, index_k_neighbour))
                                    {
                                      (*(index_face_i_sorted[1]))[index_k_neighbour].append_array(index_i_neighbour);
                                      (*(index_face_j_sorted[1]))[index_k_neighbour].append_array(0);
                                    }
                                  (*cell_faces_neighbours_corrected_bool_)[1](index_i_neighbour, 0, index_k_neighbour) += 1;
                                }
                            }
                          if (n != 0)
                            {
                              if (!(*cell_faces_neighbours_corrected_bool_)[2](index_i_neighbour, index_j_neighbour, k_offset_face))
                                {
                                  (*(index_face_i_sorted[2]))[k_offset_face].append_array(index_i_neighbour);
                                  (*(index_face_j_sorted[2]))[k_offset_face].append_array(index_j_neighbour);
                                }
                              (*cell_faces_neighbours_corrected_bool_)[2](index_i_neighbour, index_j_neighbour, k_offset_face) += 1;
                              if (j_offset_face == 0)
                                {
                                  if (!(*cell_faces_neighbours_corrected_bool_)[2](index_i_neighbour, index_j_neighbour, nb_k_layer))
                                    {
                                      (*(index_face_i_sorted[2]))[nb_k_layer].append_array(index_i_neighbour);
                                      (*(index_face_j_sorted[2]))[nb_k_layer].append_array(index_j_neighbour);
                                    }
                                  (*cell_faces_neighbours_corrected_bool_)[2](index_i_neighbour, index_j_neighbour, nb_k_layer) += 1;
                                }
                              if (i_offset_face == nb_k_layer)
                                {
                                  if (!(*cell_faces_neighbours_corrected_bool_)[2](index_i_neighbour, index_j_neighbour, 0))
                                    {
                                      (*(index_face_i_sorted[2]))[0].append_array(index_i_neighbour);
                                      (*(index_face_j_sorted[2]))[0].append_array(index_j_neighbour);
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
      cell_faces_neighbours_corrected_bool_tmp.allocate(ni, nj, nk, 0);

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
                          const int cell_faces_neighbours_ijk = cell_faces_neighbours_corrected_bool_mixed_cell[c](i+ii,j+jj,k+kk);
                          const double indic_neighbour = ref_ijk_ft_->itfce().I()(i+i_neighbour,j+j_neighbour,k+k_neighbour);
                          if (cell_faces_neighbours_ijk && indic_neighbour > LIQUID_INDICATOR_TEST)
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

      const int nb_i_layer = cell_faces_neighbours_corrected_bool[0].ni();
      const int nb_j_layer = cell_faces_neighbours_corrected_bool[0].nj();
      const int nb_k_layer = cell_faces_neighbours_corrected_bool[0].nk();

      int m,l,n;
      int index_i_problem, index_j_problem, index_k_problem;
      int l_dir, m_dir, n_dir;
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
                            const int index_i_perio = (index_i_problem + l_dir) % (nb_i_layer); // + 1);
                            const int index_j_perio = (index_j_problem + m_dir) % (nb_j_layer); // + 1);
                            const int index_k_perio = (index_k_problem + n_dir) % (nb_k_layer); // + 1);
                            const double distance = pure_neighbours_distance_to_correct[c][l][m][n];
                            cell_faces_neighbours_corrected_bool[c](index_i_perio, index_j_perio, index_k_perio) += 1;
                            double colinearity = 1.;
                            if (neighbours_colinearity_weighting_ && use_reachable_fluxes_)
                              {
                                colinearity = pure_neighbours_colinearity_to_correct[c][l][m][n];
                                neighbours_weighting_colinearity[c](index_i_perio, index_j_perio, index_k_perio) += colinearity;
                              }
                            // cell_faces_neighbours_corrected_bool[c](index_i_problem + l_dir, index_j_problem + m_dir, index_k_problem + n_dir) += 1;
                            compute_cell_neighbours_fluxes_to_correct(cell_faces_neighbours_corrected_convective,
                                                                      cell_faces_neighbours_corrected_diffusive,
                                                                      i,
                                                                      index_i_problem, index_j_problem, index_k_problem,
                                                                      index_i_perio, index_j_perio, index_k_perio,
                                                                      distance,
                                                                      c,
                                                                      colinearity,
                                                                      use_reachable_fluxes_);
                          }
                }
            }
        }
      complete_neighbours_and_weighting_colinearity(cell_faces_neighbours_corrected_bool,
                                                    cell_faces_neighbours_corrected_convective,
                                                    cell_faces_neighbours_corrected_diffusive,
                                                    neighbours_weighting_colinearity,
                                                    use_reachable_fluxes_);

      compute_cell_neighbours_mixed_cell_faces_indices_to_correct(cell_faces_neighbours_corrected_bool,
                                                                  cell_faces_neighbours_corrected_convective,
                                                                  cell_faces_neighbours_corrected_diffusive,
                                                                  neighbours_weighting_colinearity);
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
      if (convective_flux_correction_)
        cell_faces_neighbours_corrected_convective.echange_espace_virtuel();
      if (diffusive_flux_correction_)
        cell_faces_neighbours_corrected_diffusive.echange_espace_virtuel();
      // if (neighbours_colinearity_weighting_)
      {
        if (neighbours_colinearity_weighting_)
          neighbours_weighting_colinearity.echange_espace_virtuel();
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
      }
    }
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_fluxes_to_correct(FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_convective,
                                                                                          FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_diffusive,
                                                                                          const int& subproblem_index,
                                                                                          const int& index_i, const int& index_j, const int& index_k,
                                                                                          const int& index_i_perio, const int& index_j_perio, const int& index_k_perio,
                                                                                          const double& dist,
                                                                                          const int& dir,
                                                                                          const double& colinearity,
                                                                                          const int& compute_fluxes_values)
{
  if (compute_fluxes_values)
    {
      if (convective_flux_correction_)
        {
          double convective_flux = 0.;
//        compute_cell_neighbours_convective_fluxes_to_correct(cell_faces_neighbours_corrected_convective);
          compute_cell_neighbours_convective_fluxes_to_correct(convective_flux, subproblem_index, index_i, index_j, index_k, dist, dir, colinearity);
          cell_faces_neighbours_corrected_convective[dir](index_i_perio, index_j_perio, index_k_perio) += convective_flux;
        }
      if (diffusive_flux_correction_)
        {
          double diffusive_flux = 0.;
          compute_cell_neighbours_diffusive_fluxes_to_correct(diffusive_flux, subproblem_index, index_i, index_j, index_k, dist, dir, colinearity);
          cell_faces_neighbours_corrected_diffusive[dir](index_i_perio, index_j_perio, index_k_perio) += diffusive_flux;
        }
    }
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_convective_fluxes_to_correct(double& convective_flux,
                                                                                                     const int& subproblem_index,
                                                                                                     const int& index_i, const int& index_j, const int& index_k,
                                                                                                     const double& dist,
                                                                                                     const int& dir,
                                                                                                     const double& colinearity)
{
  if (!discrete_integral_)
    compute_cell_neighbours_thermal_convective_fluxes_face_centre(convective_flux,
                                                                  subproblem_index,
                                                                  index_i, index_j, index_k,
                                                                  dist,
                                                                  dir,
                                                                  colinearity);
  else
    compute_cell_neighbours_thermal_convective_fluxes_face_centre_discrete_integral(convective_flux,
                                                                                    subproblem_index,
                                                                                    index_i, index_j, index_k,
                                                                                    dist,
                                                                                    dir,
                                                                                    colinearity);
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_thermal_convective_fluxes_face_centre(double& convective_flux,
                                                                                                              const int& subproblem_index,
                                                                                                              const int& index_i, const int& index_j, const int& index_k,
                                                                                                              const double& dist,
                                                                                                              const int& dir,
                                                                                                              const double& colinearity)
{
  compute_cell_neighbours_thermal_fluxes_face_centre(convective_flux, convection,
                                                     subproblem_index,
                                                     index_i, index_j, index_k,
                                                     dist,
                                                     dir,
                                                     colinearity);
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_thermal_convective_fluxes_face_centre_discrete_integral(double& convective_flux,
    const int& subproblem_index,
    const int& index_i,
    const int& index_j,
    const int& index_k,
    const double& dist,
    const int& dir,
    const double& colinearity)
{
  compute_cell_neighbours_thermal_fluxes_face_centre_discrete_integral(convective_flux, convection,
                                                                       subproblem_index,
                                                                       index_i, index_j, index_k,
                                                                       dist,
                                                                       dir,
                                                                       colinearity);
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_diffusive_fluxes_to_correct(double& diffusive_flux,
                                                                                                    const int& subproblem_index,
                                                                                                    const int& index_i,
                                                                                                    const int& index_j,
                                                                                                    const int& index_k,
                                                                                                    const double& dist,
                                                                                                    const int& dir,
                                                                                                    const double& colinearity)
{
  if (!discrete_integral_)
    compute_cell_neighbours_thermal_diffusive_fluxes_face_centre(diffusive_flux,
                                                                 subproblem_index,
                                                                 index_i, index_j, index_k,
                                                                 dist,
                                                                 dir,
                                                                 colinearity);
  else
    compute_cell_neighbours_thermal_diffusive_fluxes_face_centre_discrete_integral(diffusive_flux,
                                                                                   subproblem_index,
                                                                                   index_i, index_j, index_k,
                                                                                   dist,
                                                                                   dir,
                                                                                   colinearity);
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_thermal_diffusive_fluxes_face_centre(double& diffusive_flux,
                                                                                                             const int& subproblem_index,
                                                                                                             const int& index_i,
                                                                                                             const int& index_j,
                                                                                                             const int& index_k,
                                                                                                             const double& dist,
                                                                                                             const int& dir,
                                                                                                             const double& colinearity)
{
  compute_cell_neighbours_thermal_fluxes_face_centre(diffusive_flux, diffusion,
                                                     subproblem_index,
                                                     index_i, index_j, index_k,
                                                     dist,
                                                     dir,
                                                     colinearity);
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_thermal_diffusive_fluxes_face_centre_discrete_integral(double& diffusive_flux,
    const int& subproblem_index,
    const int& index_i,
    const int& index_j,
    const int& index_k,
    const double& dist,
    const int& dir,
    const double& colinearity)
{
  compute_cell_neighbours_thermal_fluxes_face_centre_discrete_integral(diffusive_flux, diffusion,
                                                                       subproblem_index,
                                                                       index_i, index_j, index_k,
                                                                       dist,
                                                                       dir,
                                                                       colinearity);
}

void Corrige_flux_FT_temperature_subresolution::compute_cell_neighbours_thermal_fluxes_face_centre(double& flux,
                                                                                                   const int fluxes_type,
                                                                                                   const int& subproblem_index,
                                                                                                   const int& index_i,
                                                                                                   const int& index_j,
                                                                                                   const int& index_k,
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
                                                                                                                     const int& index_i,
                                                                                                                     const int& index_j,
                                                                                                                     const int& index_k,
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
                                                   convective_flux_x_min_max_faces_sorted_,
                                                   convective_flux_y_min_max_faces_sorted_,
                                                   convective_flux_z_min_max_faces_sorted_);
    }
  if (!diffusion_negligible_ && use_reachable_fluxes_ && fluxes_type == diffusion)
    {
      replace_cell_neighbours_thermal_fluxes_faces(cell_faces_neighbours_corrected_min_max_bool,
                                                   cell_faces_neighbours_fluxes_corrected,
                                                   diffusive_flux_x_min_max_faces_sorted_,
                                                   diffusive_flux_y_min_max_faces_sorted_,
                                                   diffusive_flux_z_min_max_faces_sorted_);
    }
}

void Corrige_flux_FT_temperature_subresolution::replace_cell_neighbours_thermal_fluxes_faces(const FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_min_max_bool,
                                                                                             const FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_fluxes_corrected,
                                                                                             std::vector<ArrOfDouble>& flux_x,
                                                                                             std::vector<ArrOfDouble>& flux_y,
                                                                                             std::vector<ArrOfDouble>& flux_z)
{
  FixedVector<std::vector<ArrOfInt>*,3> index_face_i_sorted;
  index_face_i_sorted[0] = &index_face_i_flux_x_neighbours_min_max_faces_sorted_;
  index_face_i_sorted[1] = &index_face_i_flux_y_neighbours_min_max_faces_sorted_;
  index_face_i_sorted[2] = &index_face_i_flux_z_neighbours_min_max_faces_sorted_;

  FixedVector<std::vector<ArrOfInt>*,3> index_face_j_sorted;
  index_face_j_sorted[0] = &index_face_j_flux_x_neighbours_min_max_faces_sorted_;
  index_face_j_sorted[1] = &index_face_j_flux_y_neighbours_min_max_faces_sorted_;
  index_face_j_sorted[2] = &index_face_j_flux_z_neighbours_min_max_faces_sorted_;

  FixedVector<std::vector<ArrOfDouble>*,3> fluxes;
  fluxes[0] = &flux_x;
  fluxes[1] = &flux_y;
  fluxes[2] = &flux_z;

  const int ni = cell_faces_neighbours_corrected_min_max_bool[0].ni();
  const int nj = cell_faces_neighbours_corrected_min_max_bool[0].nj();
  const int nk = cell_faces_neighbours_corrected_min_max_bool[0].nk();
  for (int c=0; c<3; c++)
    for (int k = 0; k < nk; k++)
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          if (cell_faces_neighbours_corrected_min_max_bool[c](i,j,k))
            {
              (*(fluxes[c]))[k].append_array(cell_faces_neighbours_fluxes_corrected[c](i,j,k));
              (*(index_face_i_sorted[c]))[k].append_array(i);
              (*(index_face_j_sorted[c]))[k].append_array(j);
              /*
               * Handle the periodicity
               */
              switch(c)
                {
                case 0:
                  if (i == 0)
                    {
                      (*(fluxes[c]))[k].append_array(cell_faces_neighbours_fluxes_corrected[c](i,j,k));
                      (*(index_face_i_sorted[c]))[k].append_array(ni);
                      (*(index_face_j_sorted[c]))[k].append_array(j);
                    }
                  break;
                case 1:
                  if (j == 0)
                    {
                      (*(fluxes[c]))[k].append_array(cell_faces_neighbours_fluxes_corrected[c](i,j,k));
                      (*(index_face_i_sorted[c]))[k].append_array(i);
                      (*(index_face_j_sorted[c]))[k].append_array(nj);
                    }
                  break;
                case 2:
                  if (k == 0)
                    {
                      (*(fluxes[c]))[nk].append_array(cell_faces_neighbours_fluxes_corrected[c](i,j,k));
                      (*(index_face_i_sorted[c]))[k].append_array(i);
                      (*(index_face_j_sorted[c]))[k].append_array(j);
                    }
                  break;
                default:
                  break;
                }
            }
}

void Corrige_flux_FT_temperature_subresolution::replace_temperature_cell_centre_neighbours(IJK_Field_double& temperature,
                                                                                           IJK_Field_double& temperature_neighbours,
                                                                                           IJK_Field_int& neighbours_weighting,
                                                                                           IJK_Field_double& neighbours_weighting_colinearity) const
{
  if (distance_cell_faces_from_lrs_ && find_temperature_cell_neighbours_)
    {
      temperature_neighbours.echange_espace_virtuel(temperature_neighbours.ghost());
      neighbours_weighting.echange_espace_virtuel(neighbours_weighting.ghost());
      if (neighbours_colinearity_weighting_)
        neighbours_weighting_colinearity.echange_espace_virtuel(neighbours_weighting_colinearity.ghost());
      const int ni = temperature.ni();
      const int nj = temperature.nj();
      const int nk = temperature.nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              double neighbours_weighting_ijk;
              if (neighbours_colinearity_weighting_)
                neighbours_weighting_ijk = neighbours_weighting_colinearity(i,j,k);
              else
                neighbours_weighting_ijk = (double) neighbours_weighting(i,j,k);
              if (neighbours_weighting(i,j,k))
                temperature(i,j,k) = temperature_neighbours(i,j,k) / neighbours_weighting_ijk;
            }
      temperature.echange_espace_virtuel(temperature.ghost());
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
  for (int i=0; i<ijk_intersections_subproblems_indices_.size_array(); i++)
    {
      const int intersection_ijk_cell_index = ijk_intersections_subproblems_indices_[i];
      const int ijk_indices_i = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 0);
      const int ijk_indices_j = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 1);
      const int ijk_indices_k = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 2);
      d_temperature(ijk_indices_i, ijk_indices_j, ijk_indices_k) = 0.;
    }
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

void Corrige_flux_FT_temperature_subresolution::get_discrete_surface_at_level(const int& dir, const int& level)
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
}

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

//void Corrige_flux_FT_temperature_subresolution::compute_min_max_ijk_any_reachable_fluxes(const FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_all_bool,
//                                                                                         const IJK_Field_int& neighbours_temperature_to_correct,
//                                                                                         FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_min_max_bool,
//                                                                                         const int& max_flux_per_dir,
//                                                                                         const int& check_cell_center_neighbour,
//                                                                                         const int& remove_external_neighbour_values,
//                                                                                         IJK_Field_int& neighbours_temperature_to_correct_trimmed)
//{
//  if (distance_cell_faces_from_lrs_ && find_reachable_fluxes_)
//    {
//
//      if (!keep_first_reachable_fluxes_)
//        {
//          FixedVector<IJK_Field_local_int, 3> cell_faces_neighbours_corrected_min_max_bool_local;
//
//          remove_min_max_ijk_reachable_fluxes_discontinuous(cell_faces_neighbours_corrected_all_bool,
//                                                            cell_faces_neighbours_corrected_min_max_bool_local);
//
//          const int ni = cell_faces_neighbours_corrected_all_bool[0].ni();
//          const int nj = cell_faces_neighbours_corrected_all_bool[0].nj();
//          const int nk = cell_faces_neighbours_corrected_all_bool[0].nk();
//          int i,j,k,c;
//          for (c = 0; c < 3; c++)
//            cell_faces_neighbours_corrected_min_max_bool[c].data() = 0.;
//          int c_ini = 0;
//          int c_end = 3;
//
//          if (remove_external_neighbour_values)
//            {
//              for (k = 0; k < nk; k++)
//                for (j = 0; j < nj; j++)
//                  for (i = 0; i < ni; i++)
//                    neighbours_temperature_to_correct_trimmed(i,j,k) = neighbours_temperature_to_correct(i,j,k);
//              neighbours_temperature_to_correct_trimmed.echange_espace_virtuel(neighbours_temperature_to_correct_trimmed.ghost());
//              /*
//              * i_min, i_max
//              */
//              ArrOfInt indices_found_transition_ini, indices_found_transition_end;
//              indices_found_transition_ini.set_smart_resize(1);
//              indices_found_transition_end.set_smart_resize(1);
//              ArrOfInt indices_found_ini, indices_found_end;
//              indices_found_ini.set_smart_resize(1);
//              indices_found_end.set_smart_resize(1);
//              FixedVector<ArrOfInt,2> indices_sorted;
//              if (max_flux_per_dir)
//                {
//                  c_ini = 0;
//                  c_end = 1;
//                }
//              int was_in_liquid = 0;
//              int is_in_liquid = 0;
////              int index_neighbour_cell = 0;
////						bool indices_found_bool[2] = {false, false};
////						int indices_found[2] = {0, 0};
////              bool centre_neighbour_bool = false;
//              for (c = c_ini; c < c_end; c++)
//                for (k = 0; k < nk; k++)
//                  for (j = 0; j < nj; j++)
//                    {
//                      for (i = 0; i < ni; i++)
//                        {
//                          const double indic = ref_ijk_ft_->itfce().I()(i,j,k);
//                          if (indic > VAPOUR_INDICATOR_TEST)
//                            {
//                              is_in_liquid = 1;
//                              if (i==0)
//                                {
//                                  if(indic > LIQUID_INDICATOR_TEST)
//                                    was_in_liquid = 1;
//                                  else
//                                    was_in_liquid = 0;
//                                }
//                            }
//                          else
//                            {
//                              is_in_liquid = 0;
//                              if (i==0)
//                                {
//                                  was_in_liquid = 0;
////														const double indic_last = ref_ijk_ft_->itfce().I()(i-1,j,k);
////														if(indic_last > VAPOUR_INDICATOR_TEST)
////																was_in_liquid = 1;
////														else
////																was_in_liquid = 0;
//                                }
//                            }
//
//                          if (neighbours_temperature_to_correct(i,j,k))
//                            {
//                              if (is_in_liquid != was_in_liquid)
//                                {
//                                  if (is_in_liquid)
//                                    indices_found_transition_ini.append_array(i);
//                                  else
//                                    indices_found_transition_end.append_array(i-1);
//                                  was_in_liquid = is_in_liquid;
//                                }
//                              else
//                                {
//                                  if (!neighbours_temperature_to_correct(i+1,j,k) && is_in_liquid)
//                                    indices_found_end.append_array(i);
//                                  if (!neighbours_temperature_to_correct(i-1,j,k) && is_in_liquid)
//                                    indices_found_ini.append_array(i);
//                                }
//                            }
//
////												if (indices_found_bool[0] && indices_found_bool[1])
////													for(int ii = index_neighbour_cell; ii < i; ii++)
////														neighbours_temperature_to_correct_trimmed(ii,j,k) = 0;
//                        }
//                      sort_ini_end_arrays(indices_found_transition_ini,
//                                          indices_found_transition_end,
//                                          indices_found_ini,
//                                          indices_found_end,
//                                          indices_sorted,
//                                          ni);
//
//                    }
//            }
//        }
//      cell_faces_neighbours_corrected_min_max_bool.echange_espace_virtuel();
//    }
//}
//
//void Corrige_flux_FT_temperature_subresolution::sort_ini_end_arrays(ArrOfInt& indices_found_transition_ini,
//                                                                    ArrOfInt& indices_found_transition_end,
//                                                                    ArrOfInt& indices_found_ini,
//                                                                    ArrOfInt& indices_found_end,
//                                                                    FixedVector<ArrOfInt,2>& indices_sorted,
//                                                                    const int& max_n_layer)
//{
//  int counter_transition_ini, counter_transition_end, counter_ini, counter_end = 0;
//  FixedVector<int*, 4> counters;
//  counters[0] = &counter_transition_ini;
//  counters[1] = &counter_ini;
//  counters[2] = &counter_transition_end;
//  counters[3] = &counter_end;
//  FixedVector<ArrOfInt*, 4> indices_references;
//  indices_references[0] = &indices_found_transition_ini;
//  indices_references[1] = &indices_found_ini;
//  indices_references[2] = &indices_found_transition_end;
//  indices_references[3] = &indices_found_end;
//  int max_counter = 0;
//  const int size_transition_ini = indices_found_transition_ini.size_array();
//  const int size_transition_end = indices_found_transition_end.size_array();
//  const int size_ini = indices_found_ini.size_array();
//  const int size_end = indices_found_end.size_array();
//  const int max_index = size_transition_ini + size_transition_end + size_ini + size_end - 1;
//  for (int l=0; l<2; l++)
//    indices_sorted[l].set_smart_resize(1);
//  while(counter_transition_ini<size_transition_ini
//        || counter_transition_end<size_transition_end
//        || counter_ini<size_ini
//        || counter_end<size_end)
//    {
//      std::vector<int> index_tmp = {counter_transition_ini, counter_ini, counter_transition_end, counter_end};
//      const int index_min = (int) std::distance(index_tmp.begin(), std::min_element(index_tmp.begin(), index_tmp.end()));
//      if (index_min < 2)
//        {
//          indices_sorted[0].append_array((*indices_references[index_min])(*counters[index_min]));
//          if (max_counter == max_index - 1)
//            indices_sorted[1].append_array(max_n_layer - 1);
//        }
//      else
//        {
//          indices_sorted[1].append_array((*indices_references[index_min])(*counters[index_min]));
//          if (!indices_sorted[0].size_array())
//            indices_sorted[0].append_array(0);
//        }
//      (*(counters[index_min]))++;
//      max_counter++;
//    }
//  assert(indices_sorted[0].size_array() == indices_sorted[1].size_array());
//}

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
              int * i_ref;
              int * j_ref;
              int * k_ref;
              int * i_ref_remove;
              int * j_ref_remove;
              int * k_ref_remove;
              const int * n_first;
              const int * n_bis;
              const int * n_ter;
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
                      j_ref = &index_bis;
                      k_ref = &index_first;
                      i_ref = &index_ter;
                      i_ref_remove = &val;
                      j_ref_remove = &index_bis;
                      k_ref_remove = &index_first;
                      break;
                    case 1:
                      n_first = &nk;
                      n_bis = &ni;
                      n_ter = &nj;
                      k_ref = &index_first;
                      i_ref = &index_bis;
                      j_ref = &index_ter;
                      i_ref_remove = &index_bis;
                      j_ref_remove = &val;
                      k_ref_remove = &index_first;
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
                      n_first = &nk;
                      n_bis = &nj;
                      n_ter = &ni;
                      j_ref = &index_bis;
                      k_ref = &index_first;
                      i_ref = &index_ter;
                      i_ref_remove = &val;
                      j_ref_remove = &index_bis;
                      k_ref_remove = &index_first;
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
                            sort_ini_end_arrays(indices_flux_found_ini,
                                                indices_flux_found_end,
                                                indices_fluxes_sorted,
                                                *n_ter);
                          }
                        if (indices_test || indices_flux_test)
                          remove_non_overlapping_fluxes_values(indices_sorted,
                                                               indices_fluxes_sorted,
                                                               indices_to_remove,
                                                               indices_fluxes_to_remove,
                                                               index_first, index_bis, c);

                        /*
                         * Keep max min flux values to build a convex shape
                         */

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
  while(counter_ini < size_ini || counter_end < size_end)
    {
      std::vector<int> index_tmp; // = {indices_found_ini[counter_ini],
      //	indices_found_end[counter_end]};
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
  int * i;
  int * j;
  int * k;
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
      j = &index_ter;
      k = &index_bis;
      i = &val;
      break;
    }
  const int factor_pos[3][3] = {{1, 0, 0},{0, 1, 0},{0, 0, 1}};
  std::vector<int> is_in_interval = {0, 0};
  indices_to_remove.set_smart_resize(1);
  indices_fluxes_to_remove.set_smart_resize(1);
  int ll;
  for(int l=0; l<indices_sorted[0].size_array(); l++)
    if (indices_fluxes_sorted[0].size_array())
      for (int m=0; m<indices_fluxes_sorted[0].size_array(); m++)
        {
          for (int n=0; n<2; n++)
            {
              // use indicator
              val = indices_fluxes_sorted[n][m];
              const double indic = ref_ijk_ft_->itfce().I()(*i,*j,*k);
              const double indic_prev = ref_ijk_ft_->itfce().I()(*i-factor_pos[dir][0],*j-factor_pos[dir][1],*k-factor_pos[dir][2]);
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
            }
          is_in_interval = {0, 0};
        }
  for (int m=0; m<indices_fluxes_sorted[0].size_array(); m++)
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

void Corrige_flux_FT_temperature_subresolution::sort_ijk_intersections_subproblems_indices_fluxes_by_k_layers(std::vector<ArrOfInt>& index_face_i_flux_x,
                                                                                                              std::vector<ArrOfInt>& index_face_j_flux_x,
                                                                                                              std::vector<ArrOfInt>& index_face_i_flux_y,
                                                                                                              std::vector<ArrOfInt>& index_face_j_flux_y,
                                                                                                              std::vector<ArrOfInt>& index_face_i_flux_z,
                                                                                                              std::vector<ArrOfInt>& index_face_j_flux_z,
                                                                                                              std::vector<ArrOfDouble>& flux_x,
                                                                                                              std::vector<ArrOfDouble>& flux_y,
                                                                                                              std::vector<ArrOfDouble>& flux_z,
                                                                                                              const DoubleVect& fluxes_subgrid,
                                                                                                              const int ini_index)

{
  const int nb_i_layer = ref_ijk_ft_->itfce().I().ni();
  const int nb_j_layer = ref_ijk_ft_->itfce().I().nj();
  const int nb_k_layer = ref_ijk_ft_->itfce().I().nk();
  if (!ini_index)
    {
      index_face_i_flux_x.resize(nb_k_layer);
      index_face_j_flux_x.resize(nb_k_layer);
      index_face_i_flux_y.resize(nb_k_layer);
      index_face_j_flux_y.resize(nb_k_layer);
      index_face_i_flux_z.resize(nb_k_layer + 1);
      index_face_j_flux_z.resize(nb_k_layer + 1);
    }

  flux_x.resize(nb_k_layer);
  flux_y.resize(nb_k_layer);
  flux_z.resize(nb_k_layer + 1);

  FixedVector<std::vector<ArrOfInt>*,3> index_face_i_sorted;
  index_face_i_sorted[0] = &index_face_i_flux_x;
  index_face_i_sorted[1] = &index_face_i_flux_y;
  index_face_i_sorted[2] = &index_face_i_flux_z;

  FixedVector<std::vector<ArrOfInt>*,3> index_face_j_sorted;
  index_face_j_sorted[0] = &index_face_j_flux_x;
  index_face_j_sorted[1] = &index_face_j_flux_y;
  index_face_j_sorted[2] = &index_face_j_flux_z;

  FixedVector<std::vector<ArrOfDouble>*,3> fluxes;
  fluxes[0] = &flux_x;
  fluxes[1] = &flux_y;
  fluxes[2] = &flux_z;

  for (int dir=0; dir<3; dir++)
    for (int k_layer=0; k_layer<nb_k_layer+1; k_layer++)
      {
        if ((dir==DIRECTION_I || dir==DIRECTION_J) && k_layer==nb_k_layer)
          break;
        if (!ini_index)
          {
            (*(index_face_i_sorted[dir]))[k_layer].reset();
            (*(index_face_j_sorted[dir]))[k_layer].reset();
            (*(index_face_i_sorted[dir]))[k_layer].set_smart_resize(1);
            (*(index_face_j_sorted[dir]))[k_layer].set_smart_resize(1);
          }

        (*(fluxes[dir]))[k_layer].reset();
        (*(fluxes[dir]))[k_layer].set_smart_resize(1);
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
      (*(fluxes[dir]))[k].append_array(fluxes_subgrid[i_flux]);
      (*(index_face_i_sorted[dir]))[k].append_array(i_pure_face_to_correct[i_flux]);
      (*(index_face_j_sorted[dir]))[k].append_array(j_pure_face_to_correct[i_flux]);

      switch(dir)
        {
        case DIRECTION_I:
          if (i_pure_face_to_correct[i_flux] == 0)
            {
              (*(fluxes[dir]))[k].append_array(fluxes_subgrid[i_flux]);
              (*(index_face_i_sorted[dir]))[k].append_array(nb_i_layer);
              (*(index_face_j_sorted[dir]))[k].append_array(j_pure_face_to_correct[i_flux]);
            }
          if (i_pure_face_to_correct[i_flux] == nb_i_layer)
            {
              (*(fluxes[dir]))[k].append_array(fluxes_subgrid[i_flux]);
              (*(index_face_i_sorted[dir]))[k].append_array(0);
              (*(index_face_j_sorted[dir]))[k].append_array(j_pure_face_to_correct[i_flux]);
            }
          break;
        case DIRECTION_J:
          if (j_pure_face_to_correct[i_flux] == 0)
            {
              (*(fluxes[dir]))[k].append_array(fluxes_subgrid[i_flux]);
              (*(index_face_i_sorted[dir]))[k].append_array(i_pure_face_to_correct[i_flux]);
              (*(index_face_j_sorted[dir]))[k].append_array(nb_j_layer);
            }
          if (j_pure_face_to_correct[i_flux] == nb_j_layer)
            {
              (*(fluxes[dir]))[k].append_array(fluxes_subgrid[i_flux]);
              (*(index_face_i_sorted[dir]))[k].append_array(i_pure_face_to_correct[i_flux]);
              (*(index_face_j_sorted[dir]))[k].append_array(0);
            }
          break;
        case DIRECTION_K:
          if (k == nb_k_layer)
            {
              (*(fluxes[dir]))[0].append_array(fluxes_subgrid[i_flux]);
              (*(index_face_i_sorted[dir]))[0].append_array(i_pure_face_to_correct[i_flux]);
              (*(index_face_j_sorted[dir]))[0].append_array(j_pure_face_to_correct[i_flux]);
            }
          if (k == 0)
            {
              (*(fluxes[dir]))[nb_k_layer].append_array(fluxes_subgrid[i_flux]);
              (*(index_face_i_sorted[dir]))[nb_k_layer].append_array(i_pure_face_to_correct[i_flux]);
              (*(index_face_j_sorted[dir]))[nb_k_layer].append_array(j_pure_face_to_correct[i_flux]);
            }
          break;
        }
    }
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
                                     convective_fluxes_, counter);
      counter ++;
    }
  if (!diffusion_negligible_)
    store_any_cell_faces_corrected(cell_faces_corrected_bool,
                                   cell_faces_corrected_diffusive,
                                   diffusive_fluxes_, counter);
}

void Corrige_flux_FT_temperature_subresolution::store_any_cell_faces_corrected(FixedVector<IJK_Field_int,3>& cell_faces_corrected_bool,
                                                                               FixedVector<IJK_Field_double,3>& cell_faces_corrected,
                                                                               const DoubleVect& fluxes,
                                                                               const int counter)
{
  for (int c=0; c<3; c++)
    {
      cell_faces_corrected_bool[c].data() = 0.;
      cell_faces_corrected[c].data() = 0.;
    }
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

void Corrige_flux_FT_temperature_subresolution::sort_ijk_intersections_subproblems_indices_by_k_layers()
{
  if (!convection_negligible_)
    {
      Cerr << "Sort the thermal convective fluxes" << finl;
      sort_ijk_intersections_subproblems_indices_fluxes_by_k_layers(index_face_i_flux_x_sorted_,
                                                                    index_face_j_flux_x_sorted_,
                                                                    index_face_i_flux_y_sorted_,
                                                                    index_face_j_flux_y_sorted_,
                                                                    index_face_i_flux_z_sorted_,
                                                                    index_face_j_flux_z_sorted_,
                                                                    convective_flux_x_sorted_,
                                                                    convective_flux_y_sorted_,
                                                                    convective_flux_z_sorted_,
                                                                    convective_fluxes_,
                                                                    flux_init_);
      Cerr << "Thermal Sub-resolutions convective fluxes are now sorted" << finl;
      flux_init_ = 1;
    }
  if (!diffusion_negligible_)
    {
      Cerr << "Sort the thermal diffusive fluxes" << finl;
      sort_ijk_intersections_subproblems_indices_fluxes_by_k_layers(index_face_i_flux_x_sorted_,
                                                                    index_face_j_flux_x_sorted_,
                                                                    index_face_i_flux_y_sorted_,
                                                                    index_face_j_flux_y_sorted_,
                                                                    index_face_i_flux_z_sorted_,
                                                                    index_face_j_flux_z_sorted_,
                                                                    diffusive_flux_x_sorted_,
                                                                    diffusive_flux_y_sorted_,
                                                                    diffusive_flux_z_sorted_,
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
                                     index_face_i_flux_x_sorted_,
                                     index_face_j_flux_x_sorted_,
                                     index_face_i_flux_y_sorted_,
                                     index_face_j_flux_y_sorted_,
                                     index_face_i_flux_z_sorted_,
                                     index_face_j_flux_z_sorted_,
                                     convective_flux_x_sorted_,
                                     convective_flux_y_sorted_,
                                     convective_flux_z_sorted_,
                                     k_layer,
                                     dir);
      else
        corrige_flux_faceIJ_any_flux(flux,
                                     index_face_i_flux_x_neighbours_min_max_faces_sorted_,
                                     index_face_j_flux_x_neighbours_min_max_faces_sorted_,
                                     index_face_i_flux_y_neighbours_min_max_faces_sorted_,
                                     index_face_j_flux_y_neighbours_min_max_faces_sorted_,
                                     index_face_i_flux_z_neighbours_min_max_faces_sorted_,
                                     index_face_j_flux_z_neighbours_min_max_faces_sorted_,
                                     convective_flux_x_min_max_faces_sorted_,
                                     convective_flux_y_min_max_faces_sorted_,
                                     convective_flux_z_min_max_faces_sorted_,
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
                                     index_face_i_flux_x_sorted_,
                                     index_face_j_flux_x_sorted_,
                                     index_face_i_flux_y_sorted_,
                                     index_face_j_flux_y_sorted_,
                                     index_face_i_flux_z_sorted_,
                                     index_face_j_flux_z_sorted_,
                                     diffusive_flux_x_sorted_,
                                     diffusive_flux_y_sorted_,
                                     diffusive_flux_z_sorted_,
                                     k_layer,
                                     dir);
      else
        corrige_flux_faceIJ_any_flux(flux,
                                     index_face_i_flux_x_neighbours_min_max_faces_sorted_,
                                     index_face_j_flux_x_neighbours_min_max_faces_sorted_,
                                     index_face_i_flux_y_neighbours_min_max_faces_sorted_,
                                     index_face_j_flux_y_neighbours_min_max_faces_sorted_,
                                     index_face_i_flux_z_neighbours_min_max_faces_sorted_,
                                     index_face_j_flux_z_neighbours_min_max_faces_sorted_,
                                     diffusive_flux_x_min_max_faces_sorted_,
                                     diffusive_flux_y_min_max_faces_sorted_,
                                     diffusive_flux_z_min_max_faces_sorted_,
                                     k_layer,
                                     dir);
    }
}

void Corrige_flux_FT_temperature_subresolution::corrige_flux_faceIJ_any_flux(IJK_Field_local_double *const flux,
                                                                             std::vector<ArrOfInt>& index_face_i_flux_x_sorted,
                                                                             std::vector<ArrOfInt>& index_face_j_flux_x_sorted,
                                                                             std::vector<ArrOfInt>& index_face_i_flux_y_sorted,
                                                                             std::vector<ArrOfInt>& index_face_j_flux_y_sorted,
                                                                             std::vector<ArrOfInt>& index_face_i_flux_z_sorted,
                                                                             std::vector<ArrOfInt>& index_face_j_flux_z_sorted,
                                                                             std::vector<ArrOfDouble>& subgrid_fluxes_x,
                                                                             std::vector<ArrOfDouble>& subgrid_fluxes_y,
                                                                             std::vector<ArrOfDouble>& subgrid_fluxes_z,
                                                                             const int k_layer,
                                                                             const int dir)
{
  FixedVector<std::vector<ArrOfDouble>*,3> subgrid_fluxes;
  subgrid_fluxes[0] = &subgrid_fluxes_x;
  subgrid_fluxes[1] = &subgrid_fluxes_y;
  subgrid_fluxes[2] = &subgrid_fluxes_z;

  FixedVector<std::vector<ArrOfInt>*,3> index_face_i_sorted;
  index_face_i_sorted[0] = &index_face_i_flux_x_sorted;
  index_face_i_sorted[1] = &index_face_i_flux_y_sorted;
  index_face_i_sorted[2] = &index_face_i_flux_z_sorted;

  FixedVector<std::vector<ArrOfInt>*,3> index_face_j_sorted;
  index_face_j_sorted[0] = &index_face_j_flux_x_sorted;
  index_face_j_sorted[1] = &index_face_j_flux_y_sorted;
  index_face_j_sorted[2] = &index_face_j_flux_z_sorted;

  const int nb_fluxes = (*(subgrid_fluxes[dir]))[k_layer].size_array();
  for (int i_flux=0; i_flux<nb_fluxes; i_flux++)
    {
      const int i=(*(index_face_i_sorted[dir]))[k_layer][i_flux];
      const int j=(*(index_face_j_sorted[dir]))[k_layer][i_flux];
      const double flux_ij = (*(subgrid_fluxes[dir]))[k_layer][i_flux];
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
