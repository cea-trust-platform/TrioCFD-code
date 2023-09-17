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
// File      : IJK_One_Dimensional_Subproblem.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_One_Dimensional_Subproblem_included
#define IJK_One_Dimensional_Subproblem_included

#include <Objet_U.h>
#include <IJK_Field.h>
#include <IJK_Interfaces.h>
#include <Linear_algebra_tools.h>
#include <FixedVector.h>
#include <TRUSTArrays.h>
#include <TRUSTTab.h>
#include <Vecteur3.h>
#include <Matrice33.h>
#include <Matrice.h>
#include <IJK_Finite_Difference_One_Dimensional_Matrix_Assembler.h>

// #include <TRUSTArray.h>
#define INVALID_TEMPERATURE 1e10
#define INVALID_FIELD 1e10
#define INVALID_VELOCITY 1e-12
#define INVALID_INTERP 1.e20
#define INVALID_INTERP_TEST 1.e19

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_One_Dimensional_Subproblem
//
// <Description of class IJK_One_Dimensional_Subproblem>
//
/////////////////////////////////////////////////////////////////////////////
class IJK_FT_double;

class IJK_One_Dimensional_Subproblem : public Objet_U
{

  Declare_instanciable( IJK_One_Dimensional_Subproblem ) ;

public :
  IJK_One_Dimensional_Subproblem(const IJK_FT_double& ijk_ft);
  void associer(const IJK_FT_double& ijk_ft) { ref_ijk_ft_ = ijk_ft; };
  void associate_sub_problem_to_inputs(int debug,
                                       int sub_problem_index,
                                       int i, int j, int k, int compo_connex,
                                       double distance,
                                       double curvature,
                                       double interfacial_area,
                                       ArrOfDouble facet_barycentre,
                                       ArrOfDouble normal_vector,
                                       double bubble_rising_velocity,
                                       ArrOfDouble bubble_rising_vector,
                                       ArrOfDouble bubbles_barycentre,
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
                                       const IJK_Field_double& eulerian_distance,
                                       const IJK_Field_double& eulerian_curvature,
                                       const IJK_Field_double& eulerian_interfacial_area,
                                       const FixedVector<IJK_Field_double, 3>& eulerian_normal_vect,
                                       const FixedVector<IJK_Field_double, 3>& eulerian_facets_barycentre,
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
                                       const int& source_terms_correction);

  void compute_radial_convection_diffusion_operators();
  void impose_boundary_conditions(DoubleVect& thermal_subproblems_rhs_assembly,
                                  const int& boundary_condition_interface,
                                  const double& interfacial_boundary_condition_value,
                                  const int& impose_boundary_condition_interface_from_simulation,
                                  const int& boundary_condition_end,
                                  const double& end_boundary_condition_value,
                                  const int& impose_user_boundary_condition_end_value);
  void compute_add_source_terms();
  void retrieve_temperature_solution();
  void compute_local_temperature_gradient_solution();
  double get_interfacial_gradient_corrected() const;

  void compute_local_velocity_gradient();
  double get_normal_velocity_normal_gradient() const;
  double get_tangential_velocity_normal_gradient() const;
  double get_second_tangential_velocity_normal_gradient() const;
  double get_azymuthal_velocity_normal_gradient() const;

  void get_ijk_indices(int& i, int& j, int& k) const;
  double get_field_profile_at_point(const double& dist, const DoubleVect& field) const;
  double get_temperature_profile_at_point(const double& dist) const;
  double get_temperature_gradient_profile_at_point(const double& dist, const int& dir) const;
  void thermal_subresolution_outputs();

  double get_min_temperature() const;
  double get_max_temperature() const;
  double get_min_temperature_domain_ends() const;
  double get_max_temperature_domain_ends() const;
protected :
  void associate_cell_ijk(int i, int j, int k) { index_i_ = i; index_j_=j; index_k_=k; };
  void associate_compos(int compo_connex) { compo_connex_ = compo_connex; };
  void associate_compos(int compo_connex, int compo_group) { compo_connex_ = compo_connex; compo_group_ = compo_group; };
  void associate_interface_related_parameters(double distance, double curvature, double interfacial_area, ArrOfDouble facet_barycentre, ArrOfDouble normal_vector)
  {
    distance_ = distance;
    curvature_ = curvature;
    interfacial_area_ = interfacial_area;
    facet_barycentre_ = Vecteur3(facet_barycentre);
    normal_vector_compo_ = Vecteur3(normal_vector);
  };
  void associate_rising_velocity(double bubble_rising_velocity, ArrOfDouble bubble_rising_vector, ArrOfDouble bubble_barycentre)
  {
    bubble_rising_velocity_ = bubble_rising_velocity;
    bubble_rising_vector_ = Vecteur3(bubble_rising_vector);
    bubble_barycentre_ = Vecteur3(bubble_barycentre);
  };

  void associate_eulerian_fields_references(const IJK_Interfaces& interfaces,
                                            const IJK_Field_double& eulerian_distance,
                                            const IJK_Field_double& eulerian_curvature,
                                            const IJK_Field_double& eulerian_interfacial_area,
                                            const FixedVector<IJK_Field_double, 3>& eulerian_normal_vect,
                                            const FixedVector<IJK_Field_double, 3>& eulerian_facets_barycentre,
                                            const IJK_Field_double& temperature,
                                            const IJK_Field_double& temperature_ft,
                                            const FixedVector<IJK_Field_double, 3>& velocity,
                                            const FixedVector<IJK_Field_double, 3>& velocity_ft,
                                            const FixedVector<IJK_Field_double, 3>& grad_T_elem,
                                            const FixedVector<IJK_Field_double, 3>& hess_diag_T_elem,
                                            const FixedVector<IJK_Field_double, 3>& hess_cross_T_elem);
  void associate_probe_parameters(const int& points_per_thermal_subproblem,
                                  const double& alpha,
                                  const double& lambda,
                                  const double& coeff_distance_diagonal,
                                  const double& cell_diagonal,
                                  const double& dr_base,
                                  const DoubleVect& radial_coordinates);
  void associate_finite_difference_operators(const Matrice& radial_first_order_operator_raw,
                                             const Matrice& radial_second_order_operator_raw,
                                             const Matrice& radial_first_order_operator,
                                             const Matrice& radial_second_order_operator,
                                             Matrice& radial_diffusion_matrix,
                                             Matrice& radial_convection_matrix);
  void initialise_thermal_probe();
  void compute_interface_basis_vectors();
  void compute_pure_spherical_basis_vectors();
  const int *  increase_number_of_points();
  void compute_first_order_operator_local(Matrice& radial_first_order_operator);
  void compute_second_order_operator_local(Matrice& second_first_order_operator);
  void recompute_finite_difference_matrices();
  void initialise_radial_convection_operator_local();
  void initialise_radial_diffusion_operator_local();
  void interpolate_project_velocities_on_probes();
  void interpolate_cartesian_velocities_on_probes();
  void project_velocities_on_probes();
  void correct_radial_velocity();
  void project_cartesian_onto_basis_vector(const DoubleVect& compo_x, const DoubleVect& compo_y, const DoubleVect& compo_z, const Vecteur3& basis, DoubleVect& projection);
  void project_basis_vector_onto_cartesian_dir(const int& dir, const DoubleVect& compo_u, const DoubleVect& compo_v, const DoubleVect& compo_w,
                                               const Vecteur3& basis_u, const Vecteur3& basis_v, const Vecteur3& basis_w,
                                               DoubleVect& projection);
  void interpolate_temperature_on_probe();
  void interpolate_temperature_gradient_on_probe();
  void project_temperature_gradient_on_probes();
  void interpolate_temperature_hessian_on_probe();
  void project_temperature_hessian_on_probes();
  void correct_tangential_temperature_gradient(DoubleVect& tangential_convection_source_terms);
  void correct_tangential_temperature_hessian(DoubleVect& tangential_diffusion_source_terms);
  void find_interval(const double& dist, int& left_interval, int& right_interval) const;

  void post_process_radial_quantities();
  void post_process_interfacial_quantities();

  int debug_;
  /*
   * FIXME: Should I use only references or just for IJK_Field_double ?
   * Should I use IJK_Field_local_double or IJK_Field_double as pointers ?
   */
  int sub_problem_index_ = 0;
  int index_i_ = 0, index_j_ = 0, index_k_ = 0;
  int compo_connex_ = -1;
  int compo_group_ = -1;
  double distance_ = 0.;
  double curvature_ = 0.;
  double interfacial_area_ = 0.;
  double osculating_radius_ = 0.;
  Vecteur3 facet_barycentre_;
  Vecteur3 normal_vector_compo_;

  double bubble_rising_velocity_ = 0.;
  Vecteur3 bubble_rising_vector_;
  Vecteur3 bubble_barycentre_;
  Vecteur3 facet_barycentre_relative_;

  Vecteur3 osculating_sphere_centre_;

  /*
   * Several ways to calculate the tangential vector !
   * Either by considering a unique tangential vector (pure spherical)
   * Or two tangential vectors (osculating sphere)
   */
  Vecteur3 interfacial_velocity_compo_;
  Vecteur3 interfacial_tangential_velocity_compo_;

  Vecteur3 first_tangential_vector_compo_;
  Vecteur3 second_tangential_vector_compo_;

  Vecteur3 first_tangential_vector_compo_from_rising_dir_;
  Vecteur3 azymuthal_vector_compo_raw_;
  Vecteur3 azymuthal_vector_compo_;

  bool tangential_from_rising_vel = false;
  Vecteur3 * first_tangential_vector_compo_solver_;
  Vecteur3 * second_tangential_vector_compo_solver_;

  Vecteur3 interfacial_temperature_gradient_compo_;
  Matrice33 interfacial_temperature_hessian_compo_;

  double normal_interfacial_gradient_ = 0;
  Vecteur3 normal_interfacial_gradient_compo_;

  // FIXME: Should each probes have their own number of points ?
  bool global_probes_characteristics = true;

  const int * points_per_thermal_subproblem_base_;
  const int * points_per_thermal_subproblem_;
  int increased_point_numbers_ = 32;
  // FIXME: Should alpha_liq be constant, or a reference ?
  const double * alpha_;
  const double * lambda_;
  double Pr_l_ = 0.;
  const double * coeff_distance_diagonal_;
  const double * cell_diagonal_;
  double probe_length_ = 0.;
  double surface_ = 0.;

  double r_sph_ = 0.;
  double theta_sph_ = 0.;
  double phi_sph_ = 0.;
  Vecteur3 er_sph_;
  Vecteur3 etheta_sph_;
  Vecteur3 ephi_sph_;

  /*
   * References to IJK_Field_double to avoid copy of large fields
   * Similar to operators !
   */
  const IJK_Interfaces * interfaces_;
  const IJK_Field_double * eulerian_distance_;
  const IJK_Field_double * eulerian_curvature_;
  const IJK_Field_double * eulerian_interfacial_area_;
  const FixedVector<IJK_Field_double, 3> * eulerian_normal_vect_;
  const FixedVector<IJK_Field_double, 3> * eulerian_facets_barycentre_;

  const IJK_Field_double * temperature_;
  const IJK_Field_double * temperature_ft_;
  const FixedVector<IJK_Field_double, 3> * velocity_;
  const FixedVector<IJK_Field_double, 3> * velocity_ft_;

  const FixedVector<IJK_Field_double, 3> * grad_T_elem_;
  const FixedVector<IJK_Field_double, 3> * hess_diag_T_elem_;
  const FixedVector<IJK_Field_double, 3> * hess_cross_T_elem_;

  const double * dr_base_ = 0;
  // FIXME: Should I use DoubleTab instead ?
  const DoubleVect* radial_coordinates_base_;

  const Matrice *radial_first_order_operator_raw_base_;
  const Matrice *radial_second_order_operator_raw_base_;
  const Matrice *radial_first_order_operator_base_;
  const Matrice *radial_second_order_operator_base_;
  const Matrice *radial_first_order_operator_;
  const Matrice *radial_second_order_operator_;
  Matrice radial_first_order_operator_local_;
  Matrice radial_second_order_operator_local_;

  /*
   * Pointers to non-constant matrice
   * FIXME: Should I declare constant pointers ?
   */
  Matrice *radial_diffusion_matrix_base_;
  Matrice *radial_convection_matrix_base_;
  const Matrice *radial_velocity_convection_matrix_base_;
//  const Matrice* tangential_velocity_convection_matrix_base_;
//  const Matrice* azymuthal_velocity_convection_matrix_base_;

  double dr_=0.;
  double dr_inv_=0.;
  const DoubleVect * radial_coordinates_;
  DoubleVect radial_coordinates_modified_;
  DoubleVect osculating_radial_coordinates_;
  DoubleVect osculating_radial_coordinates_inv_;
  DoubleTab radial_coordinates_cartesian_compo_;
  DoubleTab osculating_radial_coordinates_cartesian_compo_;
  DoubleTab coordinates_cartesian_compo_;

  DoubleVect x_velocity_;
  DoubleVect y_velocity_;
  DoubleVect z_velocity_;
  DoubleVect radial_velocity_;
  DoubleVect radial_velocity_corrected_;
  DoubleVect first_tangential_velocity_;
  DoubleVect azymuthal_velocity_;
  DoubleVect second_tangential_velocity_;
  DoubleVect radial_convection_prefactor_;
  DoubleVect temperature_interp_;
  FixedVector<DoubleVect, 3> grad_T_elem_interp_;
  FixedVector<DoubleVect, 3> hess_diag_T_elem_interp_;
  FixedVector<DoubleVect, 3> hess_cross_T_elem_interp_;

  int source_terms_type_=0;
  enum Source_terms { linear_diffusion, spherical_diffusion, tangential_conv_2D, tangential_conv_3D, tangential_conv_2D_tangential_diffusion_2D, tangential_conv_3D_tangentual_diffusion_3D};
  DoubleVect normal_temperature_gradient_;
  DoubleVect tangential_temperature_gradient_first_;
  DoubleVect tangential_temperature_gradient_second_;
  DoubleVect azymuthal_temperature_gradient_;
  DoubleVect tangential_hessian_contribution_;
  DoubleVect tangential_convection_source_terms_first_;
  DoubleVect tangential_convection_source_terms_second_;
  DoubleVect tangential_convection_source_terms_;
  DoubleVect tangential_diffusion_source_terms_;
  DoubleVect source_terms_;

  int correct_radial_velocity_;
  int correct_tangential_temperature_gradient_;
  int correct_tangential_temperature_hessian_;

  IJK_Finite_Difference_One_Dimensional_Matrix_Assembler * finite_difference_assembler_;
  Matrice * thermal_subproblems_matrix_assembly_;
  DoubleVect * thermal_subproblems_rhs_assembly_;
  DoubleVect * thermal_subproblems_temperature_solution_;
  DoubleVect rhs_assembly_;
  double interfacial_boundary_condition_value_;
  double end_boundary_condition_value_;
  int start_index_;
  int end_index_;
  DoubleVect temperature_solution_;
  FixedVector<DoubleVect, 3> temperature_gradient_solution_;
  DoubleVect normal_temperature_gradient_solution_;
  DoubleVect normal_temperature_double_derivative_solution_;
  DoubleVect temperature_x_gradient_solution_;
  DoubleVect temperature_y_gradient_solution_;
  DoubleVect temperature_z_gradient_solution_;
  DoubleVect thermal_flux_;

  DoubleVect normal_velocity_normal_gradient_;
  DoubleVect tangential_velocity_normal_gradient_;
  DoubleVect second_tangential_velocity_normal_gradient_;
  DoubleVect azymuthal_velocity_normal_gradient_;

  REF(IJK_FT_double) ref_ijk_ft_;

};

#endif /* IJK_One_Dimensional_Subproblem_included */
