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
      const int nb_subproblems = thermal_subproblems_->size();
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
  const int flux_sign[6] = FLUX_SIGN;
  const int nb_faces_to_correct = intersection_ijk_cell_->get_nb_faces_to_correct();
  fluxes.reset();
  fluxes.resize(nb_faces_to_correct);
  dist_.reset();
  dist_.resize(nb_faces_to_correct);
  int counter_faces = 0;
  const DoubleTab dist_interf = intersection_ijk_cell_->dist_pure_faces_interf();
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
              surf_face *= flux_sign[l];
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
      return thermal_subproblems_->get_temperature_gradient_times_diffusivity_profile_at_point(index_subproblem, dist, dir);
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
  const int flux_sign[6] = FLUX_SIGN;
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
              flux_face *= flux_sign[l];
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
      return thermal_subproblems_->get_temperature_gradient_times_diffusivity_profile_discrete_integral_at_point(index_subproblem, dist, levels_, dir);
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

void Corrige_flux_FT_temperature_subresolution::compute_ijk_pure_faces_indices()
{
  /*
   * Be careful, the ijk_intersection class is not sorting the faces the same way
   */
//	FixedVector<DoubleVect, 3>& ijk_faces_to_correct = intersection_ijk_cell_->get_set_ijk_pure_face_to_correct();
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

void Corrige_flux_FT_temperature_subresolution::sort_ijk_intersections_subproblems_indices_by_k_layers()
{
  // TODO: Get the numbers of k layers from another field ?
  const int nb_i_layer = ref_ijk_ft_->itfce().I().ni();
  const int nb_j_layer = ref_ijk_ft_->itfce().I().nj();
  const int nb_k_layer = ref_ijk_ft_->itfce().I().nk();
  // const int nb_faces_per_dir[3] = {nb_k_layer, nb_k_layer, nb_k_layer + 1};
  index_face_i_flux_x_sorted_.resize(nb_k_layer);
  index_face_j_flux_x_sorted_.resize(nb_k_layer);
  index_face_i_flux_y_sorted_.resize(nb_k_layer);
  index_face_j_flux_y_sorted_.resize(nb_k_layer);
  index_face_i_flux_z_sorted_.resize(nb_k_layer + 1);
  index_face_j_flux_z_sorted_.resize(nb_k_layer + 1);

  convective_flux_x_sorted_.resize(nb_k_layer);
  convective_flux_y_sorted_.resize(nb_k_layer);
  convective_flux_z_sorted_.resize(nb_k_layer + 1);
  diffusive_flux_x_sorted_.resize(nb_k_layer);
  diffusive_flux_y_sorted_.resize(nb_k_layer);
  diffusive_flux_z_sorted_.resize(nb_k_layer + 1);

  FixedVector<std::vector<ArrOfInt>*,3> index_face_i_sorted;
  index_face_i_sorted[0] = &index_face_i_flux_x_sorted_;
  index_face_i_sorted[1] = &index_face_i_flux_y_sorted_;
  index_face_i_sorted[2] = &index_face_i_flux_z_sorted_;

  FixedVector<std::vector<ArrOfInt>*,3> index_face_j_sorted;
  index_face_j_sorted[0] = &index_face_j_flux_x_sorted_;
  index_face_j_sorted[1] = &index_face_j_flux_y_sorted_;
  index_face_j_sorted[2] = &index_face_j_flux_z_sorted_;

  FixedVector<std::vector<ArrOfDouble>*,3> convective_fluxes;
  convective_fluxes[0] = &convective_flux_x_sorted_;
  convective_fluxes[1] = &convective_flux_y_sorted_;
  convective_fluxes[2] = &convective_flux_z_sorted_;

  FixedVector<std::vector<ArrOfDouble>*,3> diffusive_fluxes;
  diffusive_fluxes[0] = &diffusive_flux_x_sorted_;
  diffusive_fluxes[1] = &diffusive_flux_y_sorted_;
  diffusive_fluxes[2] = &diffusive_flux_z_sorted_;

  for (int dir=0; dir<3; dir++)
    for (int k_layer=0; k_layer<nb_k_layer+1; k_layer++)
      {
        if ((dir==DIRECTION_I || dir==DIRECTION_J) && k_layer==nb_k_layer)
          break;
        (*(index_face_i_sorted[dir]))[k_layer].reset();
        (*(index_face_j_sorted[dir]))[k_layer].reset();
        (*(index_face_i_sorted[dir]))[k_layer].set_smart_resize(1);
        (*(index_face_j_sorted[dir]))[k_layer].set_smart_resize(1);

        (*(convective_fluxes[dir]))[k_layer].reset();
        (*(diffusive_fluxes[dir]))[k_layer].reset();
        (*(convective_fluxes[dir]))[k_layer].set_smart_resize(1);
        (*(diffusive_fluxes[dir]))[k_layer].set_smart_resize(1);
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
      if (!convection_negligible_)
        (*(convective_fluxes[dir]))[k].append_array(convective_fluxes_[i_flux]);
      if (!diffusion_negligible_)
        (*(diffusive_fluxes[dir]))[k].append_array(diffusive_fluxes_[i_flux]);
      (*(index_face_i_sorted[dir]))[k].append_array(i_pure_face_to_correct[i_flux]);
      (*(index_face_j_sorted[dir]))[k].append_array(j_pure_face_to_correct[i_flux]);

      /*
       * TODO: Rewrite for convection and diffusion separately
       * Do something close to boundaries
       */
      if (dir==0)
        {
          if (i_pure_face_to_correct[i_flux] == 0)
            {
              if (!convection_negligible_)
                (*(convective_fluxes[dir]))[k].append_array(convective_fluxes_[i_flux]);
              if (!diffusion_negligible_)
                (*(diffusive_fluxes[dir]))[k].append_array(diffusive_fluxes_[i_flux]);
              (*(index_face_i_sorted[dir]))[k].append_array(nb_i_layer);
              (*(index_face_j_sorted[dir]))[k].append_array(j_pure_face_to_correct[i_flux]);
            }
          if (i_pure_face_to_correct[i_flux] == nb_i_layer)
            {
              if (!convection_negligible_)
                (*(convective_fluxes[dir]))[k].append_array(convective_fluxes_[i_flux]);
              if (!diffusion_negligible_)
                (*(diffusive_fluxes[dir]))[k].append_array(diffusive_fluxes_[i_flux]);
              (*(index_face_i_sorted[dir]))[k].append_array(0);
              (*(index_face_j_sorted[dir]))[k].append_array(j_pure_face_to_correct[i_flux]);
            }
        }
      if (dir==1)
        {
          if (j_pure_face_to_correct[i_flux] == 0)
            {
              if (!convection_negligible_)
                (*(convective_fluxes[dir]))[k].append_array(convective_fluxes_[i_flux]);
              if (!diffusion_negligible_)
                (*(diffusive_fluxes[dir]))[k].append_array(diffusive_fluxes_[i_flux]);
              (*(index_face_i_sorted[dir]))[k].append_array(i_pure_face_to_correct[i_flux]);
              (*(index_face_j_sorted[dir]))[k].append_array(nb_j_layer);
            }
          if (j_pure_face_to_correct[i_flux] == nb_j_layer)
            {
              if (!convection_negligible_)
                (*(convective_fluxes[dir]))[k].append_array(convective_fluxes_[i_flux]);
              if (!diffusion_negligible_)
                (*(diffusive_fluxes[dir]))[k].append_array(diffusive_fluxes_[i_flux]);
              (*(index_face_i_sorted[dir]))[k].append_array(i_pure_face_to_correct[i_flux]);
              (*(index_face_j_sorted[dir]))[k].append_array(0);
            }
        }
      if (dir==2)
        {
          if (k == nb_k_layer)
            {
              if (!convection_negligible_)
                (*(convective_fluxes[dir]))[0].append_array(convective_fluxes_[i_flux]);
              if (!diffusion_negligible_)
                (*(diffusive_fluxes[dir]))[0].append_array(diffusive_fluxes_[i_flux]);
              (*(index_face_i_sorted[dir]))[0].append_array(i_pure_face_to_correct[i_flux]);
              (*(index_face_j_sorted[dir]))[0].append_array(j_pure_face_to_correct[i_flux]);
            }
          if (k == 0)
            {
              if (!convection_negligible_)
                (*(convective_fluxes[dir]))[nb_k_layer].append_array(convective_fluxes_[i_flux]);
              if (!diffusion_negligible_)
                (*(diffusive_fluxes[dir]))[nb_k_layer].append_array(diffusive_fluxes_[i_flux]);
              (*(index_face_i_sorted[dir]))[nb_k_layer].append_array(i_pure_face_to_correct[i_flux]);
              (*(index_face_j_sorted[dir]))[nb_k_layer].append_array(j_pure_face_to_correct[i_flux]);
            }
        }
    }
  Cerr << "Thermal Sub-resolutions fluxes are now sorted" << finl;
}

void Corrige_flux_FT_temperature_subresolution::corrige_flux_faceIJ(IJK_Field_local_double *const flux,
                                                                    const int k_layer,
                                                                    const int dir)
{
  corrige_flux_faceIJ_any_flux(flux, convective_flux_x_sorted_, convective_flux_y_sorted_, convective_flux_z_sorted_, k_layer, dir);
}

void Corrige_flux_FT_temperature_subresolution::corrige_flux_diff_faceIJ(IJK_Field_local_double *const flux,
                                                                         const int k_layer,
                                                                         const int dir)
{
  corrige_flux_faceIJ_any_flux(flux, diffusive_flux_x_sorted_, diffusive_flux_y_sorted_, diffusive_flux_z_sorted_, k_layer, dir);
}

void Corrige_flux_FT_temperature_subresolution::corrige_flux_faceIJ_any_flux(IJK_Field_local_double *const flux,
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
  index_face_i_sorted[0] = &index_face_i_flux_x_sorted_;
  index_face_i_sorted[1] = &index_face_i_flux_y_sorted_;
  index_face_i_sorted[2] = &index_face_i_flux_z_sorted_;

  FixedVector<std::vector<ArrOfInt>*,3> index_face_j_sorted;
  index_face_j_sorted[0] = &index_face_j_flux_x_sorted_;
  index_face_j_sorted[1] = &index_face_j_flux_y_sorted_;
  index_face_j_sorted[2] = &index_face_j_flux_z_sorted_;

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
