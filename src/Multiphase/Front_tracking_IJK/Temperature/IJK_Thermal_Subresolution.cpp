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
// File      : IJK_Thermal_Subresolution.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Thermal_Subresolution.h>
#include <Param.h>
#include <IJK_Navier_Stokes_tools.h>
#include <DebogIJK.h>
#include <stat_counters.h>
#include <IJK_FT.h>
#include <Corrige_flux_FT.h>
#include <OpConvDiscQuickIJKScalar.h>
#include <IJK_Ghost_Fluid_tools.h>

Implemente_instanciable_sans_constructeur( IJK_Thermal_Subresolution, "IJK_Thermal_Subresolution", IJK_Thermal_base ) ;

IJK_Thermal_Subresolution::IJK_Thermal_Subresolution()
{
  convective_flux_correction_ = 0;
  diffusive_flux_correction_ = 0;
  override_vapour_mixed_values_ = 0;
  compute_eulerian_compo_ = 1;
  points_per_thermal_subproblem_ = 32;
}

Sortie& IJK_Thermal_Subresolution::printOn( Sortie& os ) const
{
  IJK_Thermal_base::printOn( os );
  os << "  {\n";
  os << "    convective_flux_correction" <<  " " << convective_flux_correction_ << "\n";
  os << "    diffusive_flux_correction" <<  " " << diffusive_flux_correction_ << "\n";
  os << "    override_vapour_mixed_values" <<  " " << override_vapour_mixed_values_ << "\n";
  os << "  \n}";
  return os;
}

Entree& IJK_Thermal_Subresolution::readOn( Entree& is )
{
  IJK_Thermal_base::readOn( is );
  if (ghost_fluid_)
    override_vapour_mixed_values_ = 1;
  if (convective_flux_correction_ || diffusive_flux_correction_)
    compute_grad_T_elem_ = 1;
  return is;
}

void IJK_Thermal_Subresolution::set_param( Param& param )
{
  IJK_Thermal_base::set_param(param);
  param.ajouter_flag("convective_flux_correction", &convective_flux_correction_);
  param.ajouter_flag("diffusive_flux_correction", &diffusive_flux_correction_);
  param.ajouter_flag("override_vapour_mixed_values", &override_vapour_mixed_values_);
  param.ajouter("points_per_thermal_subproblem", &points_per_thermal_subproblem_);
  param.ajouter("coeff_distance_diagonal", &coeff_distance_diagonal_);
}

int IJK_Thermal_Subresolution::initialize(const IJK_Splitting& splitting, const int idx)
{
  Cout << que_suis_je() << "::initialize()" << finl;

  uniform_lambda_ = lambda_liquid_;
  uniform_alpha_ =	lambda_liquid_ / (ref_ijk_ft_->get_rho_l() * cp_liquid_);
  //  calulate_grad_T_ = 1;
  // TODO: Reused grad T calculation
  calulate_grad_T_=0;
  /*
   * Be careful, it plays a role for allocating the fields in
   * IJK_Thermal_base::initialize
   */
  if (ghost_fluid_)
    {
      compute_grad_T_interface_ = 1;
      compute_curvature_ = 1;
      compute_distance_= 1;
    }

  int nalloc = 0;
  nalloc = IJK_Thermal_base::initialize(splitting, idx);
  //  temperature_diffusion_op_.set_uniform_lambda(uniform_lambda_);
  temperature_diffusion_op_.set_conductivity_coefficient(uniform_lambda_, temperature_, temperature_, temperature_, temperature_);

  /*TODO:
   * Change the operators to add fluxes corrections (Maybe not)
   */
  if (diffusive_flux_correction_)
    {
      temperature_diffusion_op_.typer("OpDiffUniformIJKScalarCorrection_double");
      temperature_convection_op_.initialize(splitting);
    }
  if (convective_flux_correction_)
    {
      temperature_convection_op_.typer("OpConvIJKQuickScalar_double");
      temperature_diffusion_op_.initialize(splitting);
    }

  thermal_local_subproblems_.associer(ref_ijk_ft_);

  Cout << "End of " << que_suis_je() << "::initialize()" << finl;
  return nalloc;
}

void IJK_Thermal_Subresolution::update_thermal_properties()
{
  IJK_Thermal_base::update_thermal_properties();
}

void IJK_Thermal_Subresolution::compute_diffusion_increment()
{
  /*
   * Update d_temperature
   * d_temperature_ += div_lambda_grad_T_volume_;
   * It depends on the nature of the properties (one-fluid or single-fluid)
   */
  const int ni = d_temperature_.ni();
  const int nj = d_temperature_.nj();
  const int nk = d_temperature_.nk();
  const double rhocp_l = ref_ijk_ft_->get_rho_l() * cp_liquid_;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          double rhocpVol;
          rhocpVol = rhocp_l * vol_;
          const double ope = div_coeff_grad_T_volume_(i,j,k);
          const double resu = ope/rhocpVol;
          d_temperature_(i,j,k) += resu;
        }
}

void IJK_Thermal_Subresolution::correct_temperature_for_eulerian_fluxes()
{
  if (override_vapour_mixed_values_)
    {
      const int ni = temperature_.ni();
      const int nj = temperature_.nj();
      const int nk = temperature_.nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double indic = ref_ijk_ft_->itfce().I(i,j,k);
              if (fabs(1.-indic) > 1.e-8) // Mixed cells and pure vapour cells
                { temperature_(i,j,k) = 0.; }
            }
    }
// TODO: Temperature sub-resolution (Weak coupling)
//			else
//				{
//
//				}
}

void IJK_Thermal_Subresolution::correct_temperature_for_visu()
{
  /*
   * Correct only if the temperature gradient is post-processed
   * If not we may want to reconstruct the gradient field in Python
   * using the ghost temperature !
   */
  if (liste_post_instantanes_.contient_("GRAD_T_ELEM"))
    {
      const int ni = temperature_.ni();
      const int nj = temperature_.nj();
      const int nk = temperature_.nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double temperature = temperature_(i,j,k);
              if (temperature > 0)
                temperature_(i,j,k) = 0;
            }
      temperature_.echange_espace_virtuel(temperature_.ghost());
    }
}

void IJK_Thermal_Subresolution::compute_ghost_cell_numbers_for_subproblems(const IJK_Splitting& splitting, int ghost_init)
{
  int ghost_cells;
  compute_cell_diagonal(splitting);
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  double maximum_distance = cell_diagonal_ * coeff_distance_diagonal_;
  const int max_dx = (int) trunc(maximum_distance / dx + 1);
  const int max_dy = (int) trunc(maximum_distance / dy + 1);
  const int max_dz = (int) trunc(maximum_distance / dz + 1);
  ghost_cells = std::max(max_dx, std::max(max_dy, max_dz));
  if (ghost_cells < ghost_init)
    ghost_cells = ghost_init;
  ghost_cells_ = ghost_cells;
}

void IJK_Thermal_Subresolution::compute_overall_probes_parameters()
{
  dr_ = probe_length_ / (points_per_thermal_subproblem_ - 1);
  for (int i=0; i < points_per_thermal_subproblem_; i++)
    radial_coordinates_(i) = i * dr_;

  /*
   * Compute the matrices for Finite-Differences
   */
  compute_radial_convection_diffusion_operators(radial_first_order_operator_raw_, radial_second_order_operator_raw_,
                                                radial_first_order_operator_, radial_second_order_operator_,
                                                radial_convection_matrix_, radial_diffusion_matrix_);
}

void IJK_Thermal_Subresolution::compute_radial_convection_diffusion_operators(DoubleTab& radial_first_order_operator_raw,
                                                                              DoubleTab& radial_second_order_operator_raw,
                                                                              DoubleTab& radial_first_order_operator,
                                                                              DoubleTab& radial_second_order_operator,
                                                                              DoubleTab& radial_diffusion_matrix,
                                                                              DoubleTab& radial_convection_matrix)
{
  /*
   * Compute the matrices for Finite-Differences
   */
  compute_first_order_operator_raw(radial_first_order_operator_raw);
  compute_second_order_operator_raw(radial_second_order_operator_raw);
  compute_first_order_operator(radial_first_order_operator, dr_);
  compute_second_order_operator(radial_second_order_operator, dr_);
  compute_radial_convection_operator(radial_first_order_operator, radial_convection_matrix);
  compute_radial_diffusion_operator(radial_second_order_operator, radial_convection_matrix);
}

/*
 * FIXME : Should be moved to a class of tools
 */

void IJK_Thermal_Subresolution::compute_first_order_operator_raw(DoubleTab& radial_first_order_operator)
{
  /*
   * Compute the first-order matrix for Finite-Differences
   */
  // TODO:
}

void IJK_Thermal_Subresolution::compute_first_order_operator(DoubleTab& radial_first_order_operator, double dr)
{
  radial_first_order_operator /= dr;
}

void IJK_Thermal_Subresolution::compute_second_order_operator_raw(DoubleTab& radial_second_order_operator)
{
  /*
   * Compute the second-order matrix for Finite-Differences
   */
  // TODO:
}

void IJK_Thermal_Subresolution::compute_second_order_operator(DoubleTab& radial_second_order_operator, double dr)
{
  radial_second_order_operator /= pow(dr, 2);
}

void IJK_Thermal_Subresolution::compute_radial_convection_operator(const DoubleTab& radial_first_order_operator,
                                                                   DoubleTab& radial_convection_matrix)
{
  /*
   * Compute the matrices for Finite-Differences
   */

}

void IJK_Thermal_Subresolution::compute_radial_diffusion_operator(const DoubleTab& radial_second_order_operator,
                                                                  DoubleTab& radial_diffusion_matrix)
{
  /*
   * Compute the matrices for Finite-Differences
   */

}

void IJK_Thermal_Subresolution::initialise_thermal_subproblems()
{
  // FIXME : Should I use IJK_Field_local_double
  //	IJK_Field_local_double eulerian_normal_vectors_ns;
  //	IJK_Field_local_double eulerian_facets_barycentre_ns;
  //	const int ni = eulerian_normal_vectors_[0].ni;
  //	const int nj = eulerian_normal_vectors_[0].nj;
  //	const int nk = eulerian_normal_vectors_[0].nk;
  //	const int ng= eulerian_normal_vectors_[0].ghost();
  //	eulerian_normal_vectors_ns.allocate(ni, nj, nk, ng);
  //	eulerian_facets_barycentre_ns.allocate(ni, nj, nk, ng);
  // FIXME : Work with FT fields
  //	compute_mixed_cells_number(ref_ijk_ft_->itfce().I());
  //  thermal_local_subproblems_.add_subproblems(mixed_cells_number_);
  //  FixedVector<IJK_Field_double, 3> eulerian_facets_barycentre_ns;
  //  FixedVector<IJK_Field_double, 3> eulerian_normal_vectors_ns;
  //  IJK_Splitting splitting = temperature_.get_splitting();
  //  allocate_cell_vector(eulerian_facets_barycentre_ns, splitting, 0);
  //  allocate_cell_vector(eulerian_normal_vectors_ns, splitting, 0);
  //  eulerian_facets_barycentre_ns.echange_espace_virtuel();
  //  eulerian_normal_vectors_ns.echange_espace_virtuel();
  //  for (int dir=0; dir < 3; dir++)
  //    {
  //      ref_ijk_ft_->redistribute_from_splitting_ft_elem(eulerian_facets_barycentre_[dir], eulerian_facets_barycentre_ns[dir]);
  //      ref_ijk_ft_->redistribute_from_splitting_ft_elem(eulerian_normal_vectors_[dir], eulerian_normal_vectors_ns[dir]);
  //    }
  // FIXME: should I perform this on the extended fields ? or original ?
  // If so, I need to convert, distance, curvature, interfacial_area earlier !
  // Or I need to calculate gradient and hessian fields using the temperature_ft_ mesh !
  const IJK_Field_double& indicator = ref_ijk_ft_->itfce().I_ft();
  const int ni = temperature_ft_.ni();
  const int nj = temperature_ft_.nj();
  const int nk = temperature_ft_.nk();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        if (fabs(indicator(i,j,k)) > VAPOUR_INDICATOR_TEST && fabs(indicator(i,j,k)) < LIQUID_INDICATOR_TEST)
          {
            thermal_local_subproblems_.associate_sub_problem_to_inputs(i, j, k,
                                                                       eulerian_compo_connex_,
                                                                       eulerian_distance_,
                                                                       eulerian_curvature_,
                                                                       eulerian_interfacial_area_,
                                                                       eulerian_facets_barycentre_,
                                                                       eulerian_normal_vectors_,
                                                                       rising_velocities_,
                                                                       rising_vectors_,
                                                                       points_per_thermal_subproblem_,
                                                                       uniform_alpha_,
                                                                       coeff_distance_diagonal_,
                                                                       cell_diagonal_,
                                                                       dr_,
                                                                       radial_coordinates_,
                                                                       radial_first_order_operator_raw_,
                                                                       radial_second_order_operator_raw_,
                                                                       radial_first_order_operator_,
                                                                       radial_second_order_operator_,
                                                                       radial_diffusion_matrix_,
                                                                       radial_convection_matrix_,
                                                                       ref_ijk_ft_->itfce(),
                                                                       temperature_,
                                                                       temperature_ft_,
                                                                       ref_ijk_ft_->get_velocity(),
                                                                       ref_ijk_ft_->get_velocity_ft(),
                                                                       grad_T_elem_,
                                                                       hess_diag_T_elem_,
                                                                       hess_cross_T_elem_);

          }
}

/*
 * TODO:
 */
void IJK_Thermal_Subresolution::solve_thermal_subproblems()
{

}

/*
 * TODO:
 */
void IJK_Thermal_Subresolution::apply_thermal_flux_correction()
{

}

void IJK_Thermal_Subresolution::clean_thermal_subproblems()
{
  thermal_local_subproblems_.clean();
}

