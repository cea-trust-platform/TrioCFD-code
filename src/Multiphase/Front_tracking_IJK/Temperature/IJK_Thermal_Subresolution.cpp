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
  disable_mixed_cells_increment_=1;
  allow_temperature_correction_for_visu_=0;
  subproblem_temperature_extension_=0;
  disable_subresolution_=0;
  convective_flux_correction_ = 0;
  diffusive_flux_correction_ = 0;
  override_vapour_mixed_values_ = 0;
  compute_eulerian_compo_ = 1;

  points_per_thermal_subproblem_ = 32;
  dr_ = 0.;
  probe_length_=0.;

  boundary_condition_interface_ = -1;
  interfacial_boundary_condition_value_ = 0.;
  impose_boundary_condition_interface_from_simulation_=0;
  boundary_condition_end_ = 0;
  end_boundary_condition_value_ = -1;
  impose_user_boundary_condition_end_value_=0;

  source_terms_type_=2;
  source_terms_correction_=0;

  fd_solver_type_ = "Solv_Gen";
  fd_solvers_ = Motcles(4);
  {
    fd_solvers_[0] = "Solv_Gen";
    fd_solvers_[1] = "Solv_GCP";
    fd_solvers_[2] = "Solv_Cholesky";
    fd_solvers_[3] = "Solv_Cholesky";
  }
  fd_solvers_jdd_ = Motcles(4);
  {
    fd_solvers_jdd_[0] = "thermal_fd_solver_1";
    fd_solvers_jdd_[1] = "thermal_fd_solver_2";
    fd_solvers_jdd_[2] = "thermal_fd_solver_3";
    fd_solvers_jdd_[3] = "thermal_fd_solver_4";
  }
  fd_solver_rank_ = 0;
  // one_dimensional_advection_diffusion_thermal_solver_.nommer("finite_difference_solver");
  // one_dimensional_advection_diffusion_thermal_solver_.typer("Solv_GCP");
  discrete_integral_ = 0;
  quadtree_levels_ = 1;
  advected_frame_of_reference_=0;
  neglect_frame_of_reference_radial_advection_=1;
  approximate_temperature_increment_=0;
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
  if (!ghost_fluid_)
    approximate_temperature_increment_ = 0;
  if (convective_flux_correction_ || diffusive_flux_correction_)
    compute_grad_T_elem_ = 1;
  if (boundary_condition_interface_ == -1 && (!impose_boundary_condition_interface_from_simulation_ && interfacial_boundary_condition_value_!=0.))
    {
      boundary_condition_interface_ = 0; // Dirichlet
      impose_boundary_condition_interface_from_simulation_ = 0;
    }
  if (boundary_condition_end_ == -1 && (!impose_user_boundary_condition_end_value_ && end_boundary_condition_value_!=-1.))
    {
      boundary_condition_end_ = 0; // Dirichlet
      impose_user_boundary_condition_end_value_ = 1;
    }

  return is;
}

void IJK_Thermal_Subresolution::set_param( Param& param )
{
  IJK_Thermal_base::set_param(param);
  param.ajouter_flag("disable_subresolution", &disable_subresolution_);
  param.ajouter_flag("convective_flux_correction", &convective_flux_correction_);
  param.ajouter_flag("diffusive_flux_correction", &diffusive_flux_correction_);
  param.ajouter_flag("override_vapour_mixed_values", &override_vapour_mixed_values_);
  param.ajouter("points_per_thermal_subproblem", &points_per_thermal_subproblem_);
  param.ajouter("coeff_distance_diagonal", &coeff_distance_diagonal_);
  param.ajouter("finite_difference_assembler", &finite_difference_assembler_);
  param.ajouter_flag("disable_mixed_cells_increment", &disable_mixed_cells_increment_);
  param.ajouter_flag("allow_temperature_correction_for_visu", &allow_temperature_correction_for_visu_);

  param.ajouter("boundary_condition_interface", &boundary_condition_interface_);
  param.dictionnaire("dirichlet", 0);
  param.dictionnaire("neumann", 1);
  param.dictionnaire("flux_jump", 2);
  param.ajouter("interfacial_boundary_condition_value", &interfacial_boundary_condition_value_);
  param.ajouter_flag("impose_boundary_condition_interface_from_simulation", &impose_boundary_condition_interface_from_simulation_);
  param.ajouter("boundary_condition_end", &boundary_condition_end_);
  param.dictionnaire("dirichlet", 0);
  param.dictionnaire("neumann",1);
  param.ajouter("end_boundary_condition_value", &end_boundary_condition_value_);
  param.ajouter_flag("impose_user_boundary_condition_end_value", &impose_user_boundary_condition_end_value_);

  param.ajouter("source_terms_type", &source_terms_type_);
  param.dictionnaire("linear_diffusion", 0);
  param.dictionnaire("spherical_diffusion",1);
  param.dictionnaire("tangential_conv_2D", 2);
  param.dictionnaire("tangential_conv_3D", 3);
  param.dictionnaire("tangential_conv_2D_tangential_diffusion_2D", 4);
  param.dictionnaire("tangential_conv_3D_tangentual_diffusion_3D", 5);
  param.ajouter_flag("source_terms_correction", &source_terms_correction_);

  // param.ajouter("thermal_fd_solver", &one_dimensional_advection_diffusion_thermal_solver_);
  param.ajouter("thermal_fd_solver", &one_dimensional_advection_diffusion_thermal_solver_);

  param.ajouter_flag("discrete_integral", &discrete_integral_);
  param.ajouter("quadtree_levels", &quadtree_levels_);
  param.ajouter_flag("advected_frame_of_reference", &advected_frame_of_reference_);
  param.ajouter_flag("neglect_frame_of_reference_radial_advection", &neglect_frame_of_reference_radial_advection_);
  param.ajouter_flag("approximate_temperature_increment", &approximate_temperature_increment_);

  //  for (int i=0; i<fd_solvers_jdd_.size(); i++)
  //    param.ajouter_non_std(fd_solvers_jdd_[i], this);
}

int IJK_Thermal_Subresolution::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  read_fd_solver(mot, is);
  return 1;
}

void IJK_Thermal_Subresolution::read_fd_solver(const Motcle& mot, Entree& is)
{
  Nom type = "";
  fd_solver_rank_ = fd_solvers_jdd_.search(mot);
  type += fd_solvers_[fd_solver_rank_];
  Cerr << "Finite difference solver: " << type << finl;
  fd_solver_type_ = fd_solvers_[fd_solver_rank_];
  one_dimensional_advection_diffusion_thermal_solver_.typer(type);
  // is >> one_dimensional_advection_diffusion_thermal_solver_;
  // Cerr << "Finite difference solver: " << fd_solver_type_ << finl;
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

  corrige_flux_.typer("Corrige_flux_FT_temperature_subresolution");
  corrige_flux_.set_fluxes_feedback_params(discrete_integral_, quadtree_levels_);
  //  temperature_diffusion_op_.set_uniform_lambda(uniform_lambda_);
  temperature_diffusion_op_.set_conductivity_coefficient(uniform_lambda_, temperature_, temperature_, temperature_, temperature_);

  /*TODO:
   * Change the operators to add fluxes corrections (Maybe not)
   */
  if (diffusive_flux_correction_)
    {
      // Use an operator that override compute_set() with corrige_flux;
      temperature_diffusion_op_.typer("OpDiffUniformIJKScalarCorrection_double");
      temperature_diffusion_op_.initialize(splitting);
      temperature_diffusion_op_.set_corrige_flux(corrige_flux_);
    }
  if (convective_flux_correction_)
    {
      // Use an operator that override compute_set() with corrige_flux;
      temperature_convection_op_.typer("OpConvQuickInterfaceIJKScalar_double");
      temperature_convection_op_.initialize(splitting);
      temperature_convection_op_.set_corrige_flux(corrige_flux_);
    }

  if (approximate_temperature_increment_)
    {
      d_temperature_uncorrected_.allocate(splitting, IJK_Splitting::ELEM, 2);
      div_coeff_grad_T_volume_uncorrected_.allocate(splitting, IJK_Splitting::ELEM, 0);
      nalloc += 2;
      // Centre2
      temperature_diffusion_op_.typer("standard");
      temperature_diffusion_op_.initialize(splitting);
      // Quick
      temperature_convection_op_.typer("standard");
      temperature_convection_op_.initialize(splitting);
    }

  thermal_local_subproblems_.associer(ref_ijk_ft_);
  compute_overall_probes_parameters();

  if (one_dimensional_advection_diffusion_thermal_solver_.est_nul())
    one_dimensional_advection_diffusion_thermal_solver_.cast_direct_solver_by_default();

  /*
   * Considered constant values
   */
  corrige_flux_.set_physical_parameters(cp_liquid_ * ref_ijk_ft_->get_rho_l(),
                                        cp_vapour_ * ref_ijk_ft_->get_rho_v(),
                                        lambda_liquid_,
                                        lambda_vapour_);
  corrige_flux_.initialize_with_subproblems(
    ref_ijk_ft_->get_splitting_ns(),
    temperature_,
    ref_ijk_ft_->itfce(),
    ref_ijk_ft_,
    ref_ijk_ft_->get_set_interface().get_set_intersection_ijk_face(),
    ref_ijk_ft_->get_set_interface().get_set_intersection_ijk_cell(),
    thermal_local_subproblems_);
  corrige_flux_->set_convection_negligible(!convective_flux_correction_);
  corrige_flux_->set_diffusion_negligible(!diffusive_flux_correction_);

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
              if (fabs(indic) < LIQUID_INDICATOR_TEST) // Mixed cells and pure vapour cells
                { temperature_(i,j,k) = 0.; }
            }
    }
}

void IJK_Thermal_Subresolution::correct_temperature_increment_for_interface_leaving_cell()
{
  /*
   * Correct only if we have not extended the temperature field across the interface (no fluxes calculation)
   */
  const int ni = d_temperature_.ni();
  const int nj = d_temperature_.nj();
  const int nk = d_temperature_.nk();
  if (disable_mixed_cells_increment_ && (!subproblem_temperature_extension_ || ghost_fluid_))
    {
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double indic = ref_ijk_ft_->itfce().I(i,j,k);
              if (fabs(indic)<LIQUID_INDICATOR_TEST) // Mixed cells and pure vapour cells
                { d_temperature_(i,j,k) = 0; }
            }
    }
}

void IJK_Thermal_Subresolution::correct_temperature_for_visu()
{
  /*
   * Correct only if the temperature gradient is post-processed
   * If not we may want to reconstruct the gradient field in Python
   * using the ghost temperature !
   */
  if (liste_post_instantanes_.contient_("GRAD_T_ELEM") && allow_temperature_correction_for_visu_)
    {
      const int ni = temperature_.ni();
      const int nj = temperature_.nj();
      const int nk = temperature_.nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              // const double temperature = temperature_(i,j,k);
              const double indic = ref_ijk_ft_->itfce().I(i,j,k);
              // if (temperature > 0)
              if (indic < VAPOUR_INDICATOR_TEST)
                temperature_(i,j,k) = 0;
            }
      temperature_.echange_espace_virtuel(temperature_.ghost());
    }
}

void IJK_Thermal_Subresolution::compute_thermal_subproblems()
{
  /*
   * FIXME: Do the first step only at the first iteration
   */
  if (debug_)
    Cerr << "Initialise thermal subproblems" << finl;
  initialise_thermal_subproblems();

  if (debug_)
    Cerr << "Compute radial subresolution convection and diffusion operators" << finl;
  compute_radial_subresolution_convection_diffusion_operators();

  if (debug_)
    Cerr << "Impose boundary conditions" << finl;
  temperature_.echange_espace_virtuel(temperature_.ghost());
  impose_subresolution_boundary_conditions();

  if (debug_)
    Cerr << "Add source terms" << finl;
  compute_add_subresolution_source_terms();

  if (debug_)
    Cerr << "Compute material derivative (modelling)" << finl;
  approximate_temperature_increment_material_derivative();

  if (debug_)
    Cerr << "Solve thermal subproblems" << finl;
  solve_thermal_subproblems();

  if (debug_)
    Cerr << "Prepare thermal flux correction" << finl;
  update_intersections();
  prepare_thermal_flux_correction();

  // apply_thermal_flux_correction();
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
  // Add one or more cell ??
  if (ghost_cells < ghost_init)
    ghost_cells = ghost_init;
  ghost_cells_ = ghost_cells;
}

void IJK_Thermal_Subresolution::compute_overall_probes_parameters()
{
  if (!disable_subresolution_)
    {
      probe_length_ = coeff_distance_diagonal_ * cell_diagonal_;
      dr_ = probe_length_ / (points_per_thermal_subproblem_ - 1);
      radial_coordinates_ = DoubleVect(points_per_thermal_subproblem_);
      for (int i=0; i < points_per_thermal_subproblem_; i++)
        radial_coordinates_(i) = i * dr_;

      /*
       * Compute the matrices for Finite-Differences (first iteration only)
       */
      compute_radial_first_second_order_operators(radial_first_order_operator_raw_, radial_second_order_operator_raw_,
                                                  radial_first_order_operator_, radial_second_order_operator_);
    }
}

void IJK_Thermal_Subresolution::compute_radial_first_second_order_operators(Matrice& radial_first_order_operator_raw,
                                                                            Matrice& radial_second_order_operator_raw,
                                                                            Matrice& radial_first_order_operator,
                                                                            Matrice& radial_second_order_operator)
{
  /*
   * Compute the matrices for Finite-Differences
   */
  compute_first_order_operator_raw(radial_first_order_operator_raw);
  compute_second_order_operator_raw(radial_second_order_operator_raw);
  radial_first_order_operator = Matrice(radial_first_order_operator_raw);
  radial_second_order_operator = Matrice(radial_second_order_operator_raw);
  compute_first_order_operator(radial_first_order_operator, dr_);
  compute_second_order_operator(radial_second_order_operator, dr_);
}

/*
 * FIXME : Should be moved to a class of tools
 */

void IJK_Thermal_Subresolution::compute_first_order_operator_raw(Matrice& radial_first_order_operator)
{
  /*
   * Compute the first-order matrix for Finite-Differences
   */
  // TODO:
  int check_nb_elem;
  check_nb_elem = finite_difference_assembler_.build(radial_first_order_operator, points_per_thermal_subproblem_, 0);
  Cerr << "Check_nb_elem" << check_nb_elem << finl;

}

void IJK_Thermal_Subresolution::compute_first_order_operator(Matrice& radial_first_order_operator, double dr)
{
  const double dr_inv = 1 / dr;
  radial_first_order_operator *= dr_inv;
}

void IJK_Thermal_Subresolution::compute_second_order_operator_raw(Matrice& radial_second_order_operator)
{
  /*
   * Compute the second-order matrix for Finite-Differences
   */
  // TODO:
  int check_nb_elem;
  check_nb_elem = finite_difference_assembler_.build(radial_second_order_operator, points_per_thermal_subproblem_, 1);
  Cerr << "Check_nb_elem" << check_nb_elem << finl;
}

void IJK_Thermal_Subresolution::compute_second_order_operator(Matrice& radial_second_order_operator, double dr)
{
  const double dr_squared_inv = 1 / pow(dr, 2);
  radial_second_order_operator *= dr_squared_inv;
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
  if (!disable_subresolution_)
    {
      const IJK_Field_double& indicator = ref_ijk_ft_->itfce().I();
      const int ni = temperature_.ni();
      const int nj = temperature_.nj();
      const int nk = temperature_.nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            if (fabs(indicator(i,j,k)) > VAPOUR_INDICATOR_TEST && fabs(indicator(i,j,k)) < LIQUID_INDICATOR_TEST)
              {
                if (debug_)
                  {
                    Cerr << "Liquid Indicator (Old): " << indicator(i,j,k) << finl;
                    Cerr << "Liquid Indicator (Next): " << ref_ijk_ft_->itfce().In()(i,j,k) << finl;
                    Cerr << "Indicator test: " << indicator_test_(i,j,k) << finl;
                    for (int i_bulle=0; i_bulle < bounding_box_.dimension(0); i_bulle++)
                      {
                        Cerr << "Bounding box : " << i_bulle << "; " << bounding_box_(i_bulle,0,0) << finl;
                        Cerr << "Bounding box : " << i_bulle << "; " << bounding_box_(i_bulle,0,1) << finl;
                        Cerr << "Bounding box : " << i_bulle << "; " << bounding_box_(i_bulle,1,0) << finl;
                        Cerr << "Bounding box : " << i_bulle << "; " << bounding_box_(i_bulle,1,1) << finl;
                        Cerr << "Bounding box : " << i_bulle << "; " << bounding_box_(i_bulle,2,0) << finl;
                        Cerr << "Bounding box : " << i_bulle << "; " << bounding_box_(i_bulle,2,1) << finl;

                        Cerr << "Min-Max Bounding box : " << i_bulle << "; " << min_max_larger_box_(i_bulle,0,0) << finl;
                        Cerr << "Min-Max Bounding box : " << i_bulle << "; " << min_max_larger_box_(i_bulle,0,1) << finl;
                        Cerr << "Min-Max Bounding box : " << i_bulle << "; " << min_max_larger_box_(i_bulle,1,0) << finl;
                        Cerr << "Min-Max Bounding box : " << i_bulle << "; " << min_max_larger_box_(i_bulle,1,1) << finl;
                        Cerr << "Min-Max Bounding box : " << i_bulle << "; " << min_max_larger_box_(i_bulle,2,0) << finl;
                        Cerr << "Min-Max Bounding box : " << i_bulle << "; " << min_max_larger_box_(i_bulle,2,1) << finl;
                      }

                  }
                thermal_local_subproblems_.associate_sub_problem_to_inputs(debug_,
                                                                           i, j, k,
                                                                           *eulerian_compo_connex_ns_,
                                                                           eulerian_distance_ns_,
                                                                           eulerian_curvature_ns_,
                                                                           eulerian_interfacial_area_ns_,
                                                                           eulerian_facets_barycentre_ns_,
                                                                           eulerian_normal_vectors_ns_,
                                                                           rising_velocities_,
                                                                           rising_vectors_,
                                                                           bubbles_barycentre_,
                                                                           advected_frame_of_reference_,
                                                                           neglect_frame_of_reference_radial_advection_,
                                                                           points_per_thermal_subproblem_,
                                                                           uniform_alpha_,
                                                                           uniform_lambda_,
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
                                                                           hess_cross_T_elem_,
                                                                           finite_difference_assembler_,
                                                                           thermal_subproblems_matrix_assembly_,
                                                                           thermal_subproblems_rhs_assembly_,
                                                                           thermal_subproblems_temperature_solution_,
                                                                           source_terms_type_,
                                                                           source_terms_correction_);

              }
    }
}

void IJK_Thermal_Subresolution::compute_radial_subresolution_convection_diffusion_operators()
{
  if (!disable_subresolution_)
    {
      initialise_radial_convection_operator(radial_first_order_operator_, radial_convection_matrix_);
      initialise_radial_diffusion_operator(radial_second_order_operator_, radial_diffusion_matrix_);
      thermal_local_subproblems_.compute_radial_convection_diffusion_operators();

      thermal_subproblems_matrix_assembly_ = radial_convection_matrix_;
      finite_difference_assembler_.sum_matrices_subproblems(thermal_subproblems_matrix_assembly_, radial_diffusion_matrix_);
    }
}

void IJK_Thermal_Subresolution::initialise_radial_convection_operator(const Matrice& radial_first_order_operator,
                                                                      Matrice& radial_convection_matrix)
{
  /*
   * Compute the convection matrices for Finite-Differences
   */
  if (debug_)
    Cerr << "Initialise the Radial convection operator" << finl;
  const int nb_subproblems = thermal_local_subproblems_.get_subproblems_counter();
  finite_difference_assembler_.initialise_matrix_subproblems(radial_convection_matrix, radial_first_order_operator_, nb_subproblems);
  if (debug_)
    Cerr << "Radial convection operator has been initialised." << finl;
}

void IJK_Thermal_Subresolution::initialise_radial_diffusion_operator(const Matrice& radial_second_order_operator,
                                                                     Matrice& radial_diffusion_matrix)
{
  /*
   * Compute the diffusion matrices for Finite-Differences
   */
  if (debug_)
    Cerr << "Initialise the Radial diffusion operator" << finl;
  const int nb_subproblems = thermal_local_subproblems_.get_subproblems_counter();
  finite_difference_assembler_.initialise_matrix_subproblems(radial_diffusion_matrix, radial_second_order_operator_, nb_subproblems);
  if (debug_)
    Cerr << "Radial diffusion operator has been initialised." << finl;
}

void IJK_Thermal_Subresolution::impose_subresolution_boundary_conditions()
{
  if (!disable_subresolution_)
    {
      thermal_subproblems_rhs_assembly_.reset();
      thermal_subproblems_rhs_assembly_.set_smart_resize(1);
      thermal_local_subproblems_.impose_boundary_conditions(thermal_subproblems_rhs_assembly_,
                                                            boundary_condition_interface_,
                                                            interfacial_boundary_condition_value_,
                                                            impose_boundary_condition_interface_from_simulation_,
                                                            boundary_condition_end_,
                                                            end_boundary_condition_value_,
                                                            impose_user_boundary_condition_end_value_);
    }
}

void IJK_Thermal_Subresolution::compute_add_subresolution_source_terms()
{
  if (!disable_subresolution_)
    {
      thermal_local_subproblems_.compute_add_source_terms();
    }
}

void IJK_Thermal_Subresolution::approximate_temperature_increment_material_derivative()
{
  if (!disable_subresolution_)
    {
      thermal_local_subproblems_.approximate_temperature_increment_material_derivative();
    }
}

/*
 * TODO: DEBUG
 */
void IJK_Thermal_Subresolution::solve_thermal_subproblems()
{
  if (!disable_subresolution_)
    {
      thermal_subproblems_temperature_solution_.reset();
      thermal_subproblems_temperature_solution_.resize(thermal_subproblems_rhs_assembly_.size());
      check_wrong_values_rhs();
      /*
       * Convert into a huge sparse matrix
       */
      thermal_subproblems_matrix_assembly_for_solver_.typer("Matrice_Morse");
      Matrice_Morse& sparse_matrix_solver  = ref_cast(Matrice_Morse, thermal_subproblems_matrix_assembly_for_solver_.valeur());
      Matrice_Bloc& bloc_matrix_solver  = ref_cast(Matrice_Bloc, thermal_subproblems_matrix_assembly_.valeur());
      bloc_matrix_solver.block_to_morse(sparse_matrix_solver);

      // const MD_Vector& md_vect = temperature_.get_md_vector();
      ArrOfInt pe_voisins_dummy;
      ArrsOfInt items_to_send_dummy;
      ArrsOfInt items_to_recv_dummy;
      ArrsOfInt blocs_to_recv_dummy;
      pe_voisins_dummy.resize(0);
      items_to_send_dummy.resize(0);
      items_to_recv_dummy.resize(0);
      blocs_to_recv_dummy.resize(0);
      MD_Vector_std md_std = MD_Vector_std(thermal_subproblems_rhs_assembly_.size(), thermal_subproblems_rhs_assembly_.size(),
                                           pe_voisins_dummy, items_to_send_dummy, items_to_recv_dummy, blocs_to_recv_dummy);
      md_.copy(md_std);
      MD_Vector_tools::creer_tableau_distribue(md_, thermal_subproblems_rhs_assembly_);
      MD_Vector_tools::creer_tableau_distribue(md_, thermal_subproblems_temperature_solution_);
      Cerr << "Finite-difference thermal sub-resolution has started !" << finl;
      one_dimensional_advection_diffusion_thermal_solver_.resoudre_systeme(thermal_subproblems_matrix_assembly_for_solver_.valeur(),
                                                                           thermal_subproblems_rhs_assembly_,
                                                                           thermal_subproblems_temperature_solution_);
      Cerr << "Finite-difference thermal sub-resolution has finished !" << finl;
      thermal_local_subproblems_.retrieve_radial_quantities();
    }
}

void IJK_Thermal_Subresolution::check_wrong_values_rhs()
{
  Cerr << "Check the modified RHS: INI" << finl;
  const double fd_second_order_magnitude = 1 / pow(dr_,2) * 1e3;
  for (int i=0; i<thermal_subproblems_rhs_assembly_.size(); i++)
    if (fabs(thermal_subproblems_rhs_assembly_[i]) > fd_second_order_magnitude)
      Cerr << "Check the modified RHS values: " << thermal_subproblems_rhs_assembly_[i] << finl;
  Cerr << "Check the modified RHS: END" << finl;
}

/*
 * TODO:
 */

void IJK_Thermal_Subresolution::update_intersections()
{
  if (!disable_subresolution_)
    corrige_flux_->update_intersections();
}

void IJK_Thermal_Subresolution::prepare_thermal_flux_correction()
{
  if (!disable_subresolution_)
    corrige_flux_->update();
}

void IJK_Thermal_Subresolution::compute_convective_fluxes_face_centre()
{
  if (!disable_subresolution_ && convective_flux_correction_)
    corrige_flux_->compute_thermal_convective_fluxes();
}

void IJK_Thermal_Subresolution::compute_diffusive_fluxes_face_centre()
{
  if (!disable_subresolution_ && diffusive_flux_correction_)
    corrige_flux_->compute_thermal_diffusive_fluxes();
}

void IJK_Thermal_Subresolution::compute_temperature_cell_centres()
{
  if (!disable_subresolution_)
    corrige_flux_->compute_temperature_cell_centre(temperature_);
}

void IJK_Thermal_Subresolution::prepare_ij_fluxes_k_layers()
{
  if ((!disable_subresolution_) && (convective_flux_correction_ || diffusive_flux_correction_))
    {
      corrige_flux_->compute_ijk_pure_faces_indices();
      corrige_flux_->sort_ijk_intersections_subproblems_indices_by_k_layers();
    }
}

void IJK_Thermal_Subresolution::set_zero_temperature_increment()
{
  if (!disable_subresolution_)
    corrige_flux_->set_zero_temperature_increment(d_temperature_);
}

void IJK_Thermal_Subresolution::clean_thermal_subproblems()
{
  if (!disable_subresolution_)
    thermal_local_subproblems_.clean();
}

void IJK_Thermal_Subresolution::clean_ijk_intersections()
{
  if (!disable_subresolution_)
    corrige_flux_->clean();
}

void IJK_Thermal_Subresolution::clean_add_thermal_subproblems()
{
  if (!disable_subresolution_)
    thermal_local_subproblems_.clean_add();
}

void IJK_Thermal_Subresolution::set_thermal_subresolution_outputs(SFichier& fic)
{
  if (!disable_subresolution_)
    thermal_local_subproblems_.thermal_subresolution_outputs(fic, rang_);
}



