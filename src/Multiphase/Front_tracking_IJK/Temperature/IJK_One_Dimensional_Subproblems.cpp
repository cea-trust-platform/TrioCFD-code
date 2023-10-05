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

Implemente_instanciable(IJK_One_Dimensional_Subproblems, "IJK_One_Dimensional_Subproblems", LIST(IJK_One_Dimensional_Subproblem));

IJK_One_Dimensional_Subproblems::IJK_One_Dimensional_Subproblems(const IJK_FT_double& ijk_ft) : IJK_One_Dimensional_Subproblems()
{
  ref_ijk_ft_ = ijk_ft;
}

void IJK_One_Dimensional_Subproblems::clean()
{
  max_subproblems_ = 0;
  subproblems_counter_ = 0;
  vide();
}

void IJK_One_Dimensional_Subproblems::clean_add()
{

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
                                                                      const IJK_Field_double& eulerian_compo_connex,
                                                                      const IJK_Field_double& eulerian_distance,
                                                                      const IJK_Field_double& eulerian_curvature,
                                                                      const IJK_Field_double& eulerian_interfacial_area,
                                                                      FixedVector<IJK_Field_double, 3> eulerian_facets_barycentre,
                                                                      FixedVector<IJK_Field_double, 3> eulerian_normal_vectors,
                                                                      ArrOfDouble rising_velocities,
                                                                      DoubleTab rising_vectors,
                                                                      DoubleTab bubbles_barycentre,
                                                                      int advected_frame_of_reference,
                                                                      int neglect_frame_of_reference_radial_advection,
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
  bool create_subproblems_iteratively = true;
  debug_ = debug;
  if (subproblems_counter_ < max_subproblems_ || create_subproblems_iteratively)
    {
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
      //
      IJK_Splitting splitting = eulerian_compo_connex.get_splitting();
      const double bubble_rising_velocity = rising_velocities(compo_connex);
      for (int dir=0; dir < 3; dir++)
        {
          facet_barycentre(dir) = eulerian_facets_barycentre[dir](i, j, k);
          normal_vector(dir) = eulerian_normal_vectors[dir](i, j, k);
          bubble_rising_vector(dir) = rising_vectors(compo_connex, dir);
          bubble_barycentre(dir) = bubbles_barycentre(compo_connex, dir);
        }
      IJK_One_Dimensional_Subproblem subproblem(ref_ijk_ft_);
      (*this).add(subproblem);
      (*this)[subproblems_counter_].associate_sub_problem_to_inputs(debug,
                                                                    subproblems_counter_,
                                                                    i, j, k, compo_connex,
                                                                    distance, curvature, interfacial_area,
                                                                    facet_barycentre, normal_vector,
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
                                                                    radial_first_order_operator_raw,
                                                                    radial_second_order_operator_raw,
                                                                    radial_first_order_operator,
                                                                    radial_second_order_operator,
                                                                    radial_diffusion_matrix,
                                                                    radial_convection_matrix,
                                                                    interfaces,
                                                                    eulerian_distance,
                                                                    eulerian_curvature,
                                                                    eulerian_interfacial_area,
                                                                    eulerian_normal_vectors,
                                                                    eulerian_facets_barycentre,
                                                                    temperature,
                                                                    temperature_ft,
                                                                    velocity,
                                                                    velocity_ft,
                                                                    grad_T_elem,
                                                                    hess_diag_T_elem,
                                                                    hess_cross_T_elem,
                                                                    finite_difference_assembler,
                                                                    thermal_subproblems_matrix_assembly,
                                                                    thermal_subproblems_rhs_assembly,
                                                                    thermal_subproblems_temperature_solution,
                                                                    source_terms_type,
                                                                    source_terms_correction);
      subproblems_counter_++;
    }
  else
    {
      Cerr << max_subproblems_ << "subproblems were expected but" << subproblems_counter_ << "subproblem try to be associated with the list of subproblems";
    }
}

void IJK_One_Dimensional_Subproblems::compute_radial_convection_diffusion_operators()
{
  for (auto& itr : *this)
    itr.compute_radial_convection_diffusion_operators();
}

void IJK_One_Dimensional_Subproblems::impose_boundary_conditions(DoubleVect& thermal_subproblems_rhs_assembly,
                                                                 const int& boundary_condition_interface,
                                                                 const double& interfacial_boundary_condition_value,
                                                                 const int& impose_boundary_condition_interface_from_simulation,
                                                                 const int& boundary_condition_end,
                                                                 const double& end_boundary_condition_value,
                                                                 const int& impose_user_boundary_condition_end_value)
{
  for (auto& itr : *this)
    itr.impose_boundary_conditions(thermal_subproblems_rhs_assembly,
                                   boundary_condition_interface,
                                   interfacial_boundary_condition_value,
                                   impose_boundary_condition_interface_from_simulation,
                                   boundary_condition_end,
                                   end_boundary_condition_value,
                                   impose_user_boundary_condition_end_value);
}

void IJK_One_Dimensional_Subproblems::compute_add_source_terms()
{
  for (auto& itr : *this)
    itr.compute_add_source_terms();
}

void IJK_One_Dimensional_Subproblems::retrieve_radial_quantities()
{
  for (auto& itr : *this)
    itr.retrieve_radial_quantities();
}

void IJK_One_Dimensional_Subproblems::retrieve_temperature_solutions()
{
  for (auto& itr : *this)
    itr.retrieve_temperature_solution();
}

void IJK_One_Dimensional_Subproblems::compute_local_temperature_gradient_solutions()
{
  for (auto& itr : *this)
    itr.compute_local_temperature_gradient_solution();
}

void IJK_One_Dimensional_Subproblems::compute_local_velocity_gradient()
{
  for (auto& itr : *this)
    itr.compute_local_velocity_gradient();
}

void IJK_One_Dimensional_Subproblems::get_subproblem_ijk_indices(int& i, int& j, int& k, int& subproblem_index) const
{
  (*this)[subproblem_index].get_ijk_indices(i,j,k);
}

double IJK_One_Dimensional_Subproblems::get_interfacial_gradient_corrected(int i)
{
  return (*this)[i].get_interfacial_gradient_corrected();
}

double IJK_One_Dimensional_Subproblems::get_temperature_profile_at_point(const int& i, const double& dist) const
{
  return (*this)[i].get_temperature_profile_at_point(dist);
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

double IJK_One_Dimensional_Subproblems::get_temperature_gradient_times_diffusivity_profile_at_point(const int& i, const double& dist, const int& dir) const
{
  return (*this)[i].get_temperature_gradient_times_diffusivity_profile_at_point(dist, dir);
}

DoubleVect IJK_One_Dimensional_Subproblems::get_temperature_gradient_profile_discrete_integral_at_point(const int& i, const double& dist, const int& level, const int& dir) const
{
  return (*this)[i].get_temperature_gradient_profile_discrete_integral_at_point(dist, level, dir);
}

DoubleVect IJK_One_Dimensional_Subproblems::get_temperature_gradient_times_diffusivity_profile_discrete_integral_at_point(const int& i, const double& dist, const int& level, const int& dir) const
{
  return (*this)[i].get_temperature_gradient_times_diffusivity_profile_discrete_integral_at_point(dist, level, dir);
}

void IJK_One_Dimensional_Subproblems::thermal_subresolution_outputs(SFichier& fic, const int rank)
{
  for (auto& itr : *this)
    itr.thermal_subresolution_outputs(fic, rank);
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

