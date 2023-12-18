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
// File      : IJK_Finite_Difference_One_Dimensional_Matrix_Assembler.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Finite_Difference_One_Dimensional_Matrix_Assembler_included
#define IJK_Finite_Difference_One_Dimensional_Matrix_Assembler_included

#include <Objet_U.h>
#include <Matrice_Bloc.h>
#include <Matrice_Morse_Sym.h>
#include <Matrice.h>
#include <FixedVector.h>
#include <Param.h>

#define MAX_ORDER_DERIVATIVE 4
#define FORTRAN_INDEX_INI 1
/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_Finite_Difference_One_Dimensional_Matrix_Assembler
//
// <Description of class IJK_Finite_Difference_One_Dimensional_Matrix_Assembler>
//
/////////////////////////////////////////////////////////////////////////////

class IJK_Finite_Difference_One_Dimensional_Matrix_Assembler : public Objet_U
{

  Declare_instanciable( IJK_Finite_Difference_One_Dimensional_Matrix_Assembler ) ;

public :
  void set_param(Param& param);
  int build(Matrice& matrix, const int& nb_elem, const int& derivative_order);
  int build_with_unknown_pattern(Matrice& matrix, const int& nb_elem, const int& derivative_order);
  int build_with_known_pattern(Matrice& matrix, const int& nb_elem, const int& derivative_order);
  void modify_rhs_for_bc(Matrice& modified_matrix,
                         DoubleVect& modified_rhs,
                         const int& ini_boundary_conditions,
                         const int& end_boundary_conditions);
  void scale_matrix_by_vector(Matrice& matrix,
                              const DoubleVect& vector,
                              const int& boundary_conditions);
  void scale_matrix_subproblem_by_vector(Matrice * matrix,
                                         const DoubleVect& vector,
                                         const int& subproblem_index,
                                         const int& boundary_conditions,
                                         const int& use_sparse_matrix,
                                         const FixedVector<ArrOfInt,6> * first_indices_sparse_matrix);
  void make_operation_on_sub_matrix_sparse(Matrice& local_sub_matrix,
                                           Matrice * matrix,
                                           const int& subproblem_index,
                                           const FixedVector<ArrOfInt,6> * first_indices_sparse_matrix);
  void recombined_local_sub_matrix_with_matrix(Matrice& local_sub_matrix,
                                               Matrice * matrix,
                                               const int& subproblem_index,
                                               const FixedVector<ArrOfInt,6> * first_indices_sparse_matrix);
  void impose_boundary_conditions(Matrice& modified_matrix,
                                  DoubleVect& mdified_rhs,
                                  const int& ini_boundary_conditions,
                                  const double& interfacial_value,
                                  const int& end_boundary_conditions,
                                  const double& end_values,
                                  const double& dr_inv,
                                  const int& first_time_step_temporal,
                                  const int& first_time_step_explicit,
                                  const DoubleVect& temperature_ini_temporal_schemes);
  void impose_boundary_conditions_subproblem(Matrice * matrix,
                                             DoubleVect * global_rhs,
                                             DoubleVect& local_rhs,
                                             const int& ini_boundary_conditions,
                                             const double& interfacial_value,
                                             const int& end_boundary_conditions,
                                             const double& end_value,
                                             const int& subproblem_index,
                                             const double& dr_inv,
                                             const int& first_time_step_temporal,
                                             const int& first_time_step_explicit,
                                             const DoubleVect& temperature_ini_temporal_schemes,
                                             const int& start_index,
                                             const FixedVector<ArrOfInt, 6> * first_indices_sparse_matrix,
                                             const int& use_sparse_matrix);
  void sum_any_matrices_subproblems(Matrice& matrix_A, Matrice& matrix_B, const int& use_sparse_matrix);
  void sum_sparse_matrices_subproblems(Matrice& matrix_A, Matrice& matrix_B);
  void sum_matrices_subproblems(Matrice& matrix_A, Matrice& matrix_B);
  void sum_matrices(Matrice& matrix_A, Matrice& matrix_B);
  void initialise_sparse_matrix_subproblems(Matrice& matrix_subproblems,
                                            Matrice& fd_operator,
                                            const int& nb_subproblems,
                                            const int& first_time_step_varying_probes,
                                            FixedVector<ArrOfInt, 6>& first_indices_sparse_matrix,
                                            const int& first_initialisation,
                                            int& initialise_sparse_indices);
  void initialise_matrix_subproblems(Matrice& matrix_subproblems,
                                     Matrice& fd_operator,
                                     const int& subproblems,
                                     const int& first_time_step_varying_probes);
  void reinitialise_any_matrix_subproblem(Matrice * matrix_subproblems,
                                          const Matrice * fd_operator,
                                          const int& nb_subproblems,
                                          const int& use_sparse_matrix,
                                          FixedVector<ArrOfInt,6> * first_indices_sparse_matrix,
                                          const int& first_initialisation);
  void reinitialise_sparse_matrix_subproblem(Matrice * matrix_subproblems,
                                             const Matrice * fd_operator,
                                             const int& nb_subproblems,
                                             FixedVector<ArrOfInt,6> * first_indices_sparse_matrix,
                                             const int& first_initialisation);
  void reinitialise_matrix_subproblem(Matrice * matrix_subproblems,
                                      const Matrice * fd_operator,
                                      const int& nb_subproblems);
  void add_source_terms(DoubleVect * thermal_subproblems_rhs_assembly,
                        const DoubleVect& rhs_assembly);
  void compute_operator(const Matrice * fd_operator, const DoubleVect& solution, DoubleVect& res);
  void apply_euler_time_step(Matrice * convection_matrix,
                             Matrice * diffusion_matrix,
                             const int& subproblem_index,
                             const double& local_time_step,
                             const double& alpha);
  void correct_sign_temporal_schemes_subproblems(Matrice * convection_matrix,
                                                 Matrice * diffusion_matrix,
                                                 const int& subproblem_index,
                                                 const double& local_time_step,
                                                 const double& alpha);

  void pre_initialise_matrix_subproblems(Matrice& matrice, Matrice& fd_operator, const int& max_subproblems_predicted);
  void pre_initialise_sparse_matrix_subproblems(Matrice& matrix_subproblems,
                                                Matrice& fd_operator,
                                                const int& max_subproblems_predicted);
  void complete_empty_matrices_initialisation(Matrice& matrix_subproblems,
                                              Matrice& fd_operator,
                                              const int& empty_problem_start_index,
                                              const int& empty_problem_end_index);

  void reduce_solver_matrix(const Matrice_Morse& thermal_subproblems_matrix_assembly_for_solver,
                            Matrice_Morse& thermal_subproblems_matrix_assembly_for_solver_reduced,
                            const int& nb_points,
                            const int& pre_initialise_thermal_subproblems_list);

protected :
  enum Fd_coefficient_type_ { identity=-1, forward, centred, backward };
  enum Precision_Order_ { first_order, second_order };
  enum Derivative_Order_ { first, second };
  // enum Boundary_conditions { dirichlet, neumann, flux_jump };
  enum Boundary_conditions { dirichlet, neumann, flux_jump };
  /*
   *  -elim_coeff_nul=0, on ne supprime pas les coefficients nuls de la matrice
   *  -elim_coeff_nul=1, on supprime les coefficients nuls de la matrice
   *  -elim_coeff_nul=2, on supprime les coefficients nuls et quasi-nuls de la matrice
   */
  enum Remove_zero_coeff { keep, remove, remove_close_zero };
  enum Equation_type { advection_diffusion, linear_diffusion };

  /*
   * Identity matrix
   */
  FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE> identity_coefficient_;
  /*
   * First order derivatives
   */
  // Forward
  const std::vector<std::vector<double>>  first_order_derivative_forward_vector_ = {{-1., 1.},
    {-3./2., 2., -1./2.},
    {-11./6., 3., -3./2., 1./3.},
    {-25./12., 4., -3., 4./3., -1./4.}
  };
  FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE> first_order_derivative_forward_;

  // Centred
  const std::vector<std::vector<double>> first_order_derivative_centred_vector_ = {{-1./2., 0., 1./2.},
    {-1./2., 0., 1./2.},
    {1./12., -2./3., 0., 2./3., -1./12.},
    {1./12., -2./3., 0., 2./3., -1./12.}
  };
  FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE> first_order_derivative_centred_;

  // Backward
  const std::vector<std::vector<double>> first_order_derivative_backward_vector_ = {{1., -1.},
    {3./2., -2., 1./2.},
    {11./6., -3., 3./2., -1./3.},
    {25./12., -4., 3., -4./3., 1./4.}
  };
  FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE> first_order_derivative_backward_;

  /*
   * Second order derivatives
   */
  // Forward
  const std::vector<std::vector<double>> second_order_derivative_forward_vector_ = {{1., -2., 1.},
    {2., -5., 4., -1.},
    {35./12., -26./3., 19./2., -14./3., 11./12.},
    {15./4., -77./6., 107./6., -13., 61./12., -5./6.}
  };
  FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE> second_order_derivative_forward_;

  // Centred
  const std::vector<std::vector<double>> second_order_derivative_centred_vector_ = {{1., -2., 1.},
    {1., -2., 1.},
    {-2., -1., 0., 1., 2.},
    {-2., -1., 0., 1., 2.}
  };
  FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE> second_order_derivative_centred_;

  // Backward
  const std::vector<std::vector<double>> second_order_derivative_backward_vector_ = {{1., -2., 1.},
    {-2., 5., -4., 1.},
    {35./12., -26./3., 19./2., -14./3., 11./12.},
    {15./4., -77./6., 107./6., -13., 61./12., -5./6.}
  };
  FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE> second_order_derivative_backward_;

  int precision_order_;
  int reduce_side_precision_;
  int core_matrix_type_;
  int equation_type_;
  int nb_elem_;
  int non_zero_elem_;
  int stencil_forward_max_;
  int stencil_centred_max_;
  int stencil_backward_max_;
  int forward_left_offset_[3] = {0, 0, 0};
  int forward_right_offset_[3] = {0, 0, 0};
  int centred_left_offset_[3] = {0, 0, 0};
  int centred_right_offset_[3] = {0, 0, 0};
  int backward_left_offset_[3] = {0, 0, 0};
  int backward_right_offset_[3] = {0, 0, 0};
  bool computed_stencil_;
  int known_pattern_=1;

  void set_operators_size(const std::vector<std::vector<double>>& fd_operator_vector, FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE>& fd_operator);
  void set_operators_indices(const std::vector<std::vector<double>>& fd_operator_vector,
                             FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE>& fd_operator,
                             const int& fd_coefficient_type);
  int non_zero_stencil_values(const FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE>& fd_operator);
  int get_max_operators(const FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE>& fd_operator_conv,
                        const FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE>& fd_operator_diff,
                        int (&stencil_left_offset) [3], int (&stencil_right_offset) [3]);
  void get_max_stencil(const int& nb_elem,
                       int& non_zero_elem_max,
                       int& stencil_forward_max,
                       int& stencil_centred_max,
                       int& stencil_backward_max);
  void compute_stencil_properties(const int nb_elem);
  void fill_lines(Matrice& matrix);

};

#endif /* IJK_Finite_Difference_One_Dimensional_Matrix_Assembler_included */
