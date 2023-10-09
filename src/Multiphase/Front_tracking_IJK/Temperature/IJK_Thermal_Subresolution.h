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
// File      : IJK_Thermal_Subresolution.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Thermal_Subresolution_included
#define IJK_Thermal_Subresolution_included

#include <IJK_Thermal_base.h>
#include <IJK_Field.h>
#include <Boundary_Conditions_Thermique.h>
#include <IJK_Splitting.h>
#include <Parser.h>
#include <IJK_Lata_writer.h>
#include <OpConvQuickIJKScalar.h>
#include <OpConvCentre2IJKScalar.h>
#include <Ouvrir_fichier.h>
#include <Corrige_flux_FT.h>
#include <TRUST_Ref.h>
#include <Operateur_IJK_elem_diff_base.h>
#include <OpConvAmontIJK.h>
#include <OpConvDiscQuickIJKScalar.h>
#include <OpConvCentre4IJK.h>
#include <IJK_One_Dimensional_Subproblems.h>
#include <IJK_Finite_Difference_One_Dimensional_Matrix_Assembler.h>
#include <IJK_SolveSys_FD_thermal.h>
#include <MD_Vector_tools.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_Thermal_Subresolution
//
// <Description of class IJK_Thermal_Subresolution>
//
/////////////////////////////////////////////////////////////////////////////


class IJK_Thermal_Subresolution : public IJK_Thermal_base
{

  Declare_instanciable( IJK_Thermal_Subresolution ) ;

public :

  int initialize(const IJK_Splitting& splitting, const int idx) override;
  void update_thermal_properties() override;
  void set_param(Param& param) override;
  void compute_ghost_cell_numbers_for_subproblems(const IJK_Splitting& splitting, int ghost_init) override;
  // Entree& read_fd_solver(Entree& is);
  void read_fd_solver(const Motcle& mot, Entree& is);
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void set_thermal_subresolution_outputs(SFichier& fic) override;

protected :
  void compute_thermal_subproblems() override;
  void compute_diffusion_increment() override;
  void correct_temperature_increment_for_interface_leaving_cell() override;
  void correct_temperature_for_eulerian_fluxes() override;
  void correct_temperature_for_visu() override;
  void compute_overall_probes_parameters();
  void compute_radial_subresolution_convection_diffusion_operators();
  void compute_source_terms_impose_subresolution_boundary_conditions();
  void compute_add_subresolution_source_terms();
  void compute_subresolution_temporal_explicit_implicit_matrices();
  void approximate_temperature_increment_material_derivative();
  void compute_radial_first_second_order_operators(Matrice& radial_first_order_operator_raw,
                                                   Matrice& radial_second_order_operator_raw,
                                                   Matrice& radial_first_order_operator,
                                                   Matrice& radial_second_order_operator);
  void compute_first_order_operator_raw(Matrice& radial_first_order_operator);
  void compute_first_order_operator(Matrice& radial_first_order_operator, double dr);
  void compute_second_order_operator(Matrice& radial_second_order_operator, double dr);
  void compute_second_order_operator_raw(Matrice& radial_second_order_operator);
  void initialise_identity_matrices(Matrice& identity_matrix,
                                    Matrice& identity_matrix_subproblems);
  void initialise_radial_convection_operator(Matrice& radial_first_order_operator,
                                             Matrice& radial_convection_matrix);
  void initialise_radial_diffusion_operator(Matrice& radial_second_order_operator,
                                            Matrice& radial_diffusion_matrix);
  void check_wrong_values_rhs();
  void initialise_thermal_subproblems();
  void solve_thermal_subproblems();
  void prepare_thermal_flux_correction();
  void update_intersections() override;
  void compute_convective_fluxes_face_centre() override;
  void compute_diffusive_fluxes_face_centre() override;
  void compute_temperature_cell_centres() override;
  void prepare_ij_fluxes_k_layers() override;
  void set_zero_temperature_increment() override;
  void clean_thermal_subproblems() override;
  void clean_ijk_intersections() override;
  void clean_add_thermal_subproblems();

  /* compute_rho_cp_u_mean() May be clearly overridden later */
  double compute_rho_cp_u_mean(const IJK_Field_double& vx) override { return IJK_Thermal_base::compute_rho_cp_u_mean(vx); };

  int disable_mixed_cells_increment_;
  int allow_temperature_correction_for_visu_;
  int disable_subresolution_;
  int diffusive_flux_correction_;
  int convective_flux_correction_;
  int subproblem_temperature_extension_; // ghost fluid extension based on the interfacial gradient computed with the subproblem

  int override_vapour_mixed_values_; // For debug purposes

  IJK_One_Dimensional_Subproblems thermal_local_subproblems_;
  int points_per_thermal_subproblem_;
  double coeff_distance_diagonal_ = 3.;
  double probe_length_;
  double dr_;
  DoubleVect radial_coordinates_;
  Matrice identity_matrix_explicit_implicit_;
  Matrice radial_first_order_operator_raw_;
  Matrice radial_second_order_operator_raw_;
  Matrice radial_first_order_operator_;
  Matrice radial_second_order_operator_;
  Matrice identity_matrix_subproblems_;
  Matrice radial_diffusion_matrix_;
  Matrice radial_convection_matrix_;
  /*
   * Thermal subproblems are regrouped in a single linear system AX=b
   * on each processor !
   */
  IJK_Finite_Difference_One_Dimensional_Matrix_Assembler finite_difference_assembler_;
  Matrice thermal_subproblems_matrix_assembly_;
  Matrice thermal_subproblems_matrix_assembly_for_solver_;
  DoubleVect thermal_subproblems_temperature_solution_ini_;
  DoubleVect thermal_subproblems_rhs_assembly_;
  DoubleVect thermal_subproblems_temperature_solution_;
  // SolveurSys one_dimensional_advection_diffusion_thermal_solver_;
  IJK_SolveSys_FD_thermal one_dimensional_advection_diffusion_thermal_solver_;
  IJK_SolveSys_FD_thermal one_dimensional_advection_diffusion_thermal_solver_implicit_;
  MD_Vector md_;
  Motcles fd_solvers_;
  Motcles fd_solvers_jdd_;
  int fd_solver_rank_;
  Nom fd_solver_type_;
  int discrete_integral_;
  int quadtree_levels_;
  // DoubleVect radial_convective_vector_prefactor_;
  // DoubleVect diffusive_vector_prefactor_;

  int boundary_condition_interface_;
  double interfacial_boundary_condition_value_;
  int impose_boundary_condition_interface_from_simulation_;
  int boundary_condition_end_;
  double end_boundary_condition_value_;
  int impose_user_boundary_condition_end_value_;
  int source_terms_type_;
  int source_terms_correction_;
  int advected_frame_of_reference_;
  int neglect_frame_of_reference_radial_advection_;
  int approximate_temperature_increment_;

  bool is_first_time_step_;
  int first_time_step_temporal_;
  int first_time_step_explicit_;
};

#endif /* IJK_Thermal_Subresolution_included */
