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

// #include <TRUSTArray.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_One_Dimensional_Subproblem
//
// <Description of class IJK_One_Dimensional_Subproblem>
//
/////////////////////////////////////////////////////////////////////////////

class IJK_One_Dimensional_Subproblem : public Objet_U
{

  Declare_instanciable( IJK_One_Dimensional_Subproblem ) ;

public :
  void associate_sub_problem_to_inputs(int i, int j, int k, int compo_connex,
                                       double distance,
                                       double curvature,
                                       double interfacial_area,
                                       ArrOfDouble facet_barycentre,
                                       ArrOfDouble normal_vector,
                                       double bubble_rising_velocity,
                                       ArrOfDouble bubble_rising_vector,
                                       const int& points_per_thermal_subproblem,
                                       const double& alpha,
                                       const double& coeff_distance_diagonal,
                                       const double& cell_diagonal,
                                       const double& dr_base,
                                       const DoubleVect& radial_coordinates,
                                       const DoubleTab& radial_first_order_operator_raw,
                                       const DoubleTab& radial_second_order_operator_raw,
                                       const DoubleTab& radial_first_order_operator,
                                       const DoubleTab& radial_second_order_operator,
                                       const DoubleTab& radial_diffusion_matrix,
                                       const DoubleTab& radial_convection_matrix,
                                       const IJK_Interfaces& interfaces,
                                       const IJK_Field_local_double& eulerian_distance,
                                       const IJK_Field_local_double& eulerian_curvature,
                                       const IJK_Field_local_double& eulerian_interfacial_area,
                                       const FixedVector<IJK_Field_double, 3>& eulerian_normal_vect,
                                       const FixedVector<IJK_Field_double, 3>& eulerian_facets_barycentre,
                                       const IJK_Field_local_double& temperature,
                                       const IJK_Field_local_double& temperature_ft,
                                       const FixedVector<IJK_Field_double, 3>& velocity,
                                       const FixedVector<IJK_Field_double, 3>& velocity_ft,
                                       const FixedVector<IJK_Field_double, 3>& grad_T_elem,
                                       const FixedVector<IJK_Field_double, 3>& hess_diag_T_elem,
                                       const FixedVector<IJK_Field_double, 3>& hess_cross_T_elem);

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
  void associate_rising_velocity(double bubble_rising_velocity, ArrOfDouble bubble_rising_vector)
  {
    bubble_rising_velocity_ = bubble_rising_velocity;
    bubble_rising_vector_ = Vecteur3(bubble_rising_vector);
  };

  void associate_eulerian_fields_references(const IJK_Interfaces& interfaces,
                                            const IJK_Field_local_double& eulerian_distance,
                                            const IJK_Field_local_double& eulerian_curvature,
                                            const IJK_Field_local_double& eulerian_interfacial_area,
                                            const FixedVector<IJK_Field_double, 3>& eulerian_normal_vect,
                                            const FixedVector<IJK_Field_double, 3>& eulerian_facets_barycentre,
                                            const IJK_Field_local_double& temperature,
                                            const IJK_Field_local_double& temperature_ft,
                                            const FixedVector<IJK_Field_double, 3>& velocity,
                                            const FixedVector<IJK_Field_double, 3>& velocity_ft,
                                            const FixedVector<IJK_Field_double, 3>& grad_T_elem,
                                            const FixedVector<IJK_Field_double, 3>& hess_diag_T_elem,
                                            const FixedVector<IJK_Field_double, 3>& hess_cross_T_elem);
  void associate_probe_parameters(const int& points_per_thermal_subproblem,
                                  const double& alpha,
                                  const double& coeff_distance_diagonal,
                                  const double& cell_diagonal,
                                  const double& dr_base,
                                  const DoubleVect& radial_coordinates);
  void associate_finite_difference_operators(const DoubleTab& radial_first_order_operator_raw,
                                             const DoubleTab& radial_second_order_operator_raw,
                                             const DoubleTab& radial_first_order_operator,
                                             const DoubleTab& radial_second_order_operator,
                                             const DoubleTab& radial_diffusion_matrix,
                                             const DoubleTab& radial_convection_matrix);
  void initialise_thermal_probe();
  const int *  increase_number_of_points();
  void initialise_finite_difference_operators();

  /*
   * FIXME: Should I use only references or just for IJK_Field_local_double ?
   * Should I use IJK_Field_local_double or IJK_Field_double as pointers ?
   */
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
  Vecteur3 bubble_centre_;

  Vecteur3 osculating_sphere_centre_;

  /*
   * Several ways to calculate the tangential vector !
   * Either by considering a unique tangential vector (pure spherical)
   * Or two tangential vectors (osculating sphere)
   */
  Vecteur3 first_tangential_vector_compo_;
  Vecteur3 azymuthal_vector_compo_;
  Vecteur3 second_tangential_vector_compo_;

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
  const double * coeff_distance_diagonal_;
  const double * cell_diagonal_;
  double probe_length_ = 0.;

  /*
   * References to IJK_Field_double to avoid copy of large fields
   * Similar to operators !
   */
  const IJK_Interfaces * interfaces_;
  const IJK_Field_local_double * eulerian_distance_;
  const IJK_Field_local_double * eulerian_curvature_;
  const IJK_Field_local_double * eulerian_interfacial_area_;
  const FixedVector<IJK_Field_double, 3> * eulerian_normal_vect_;
  const FixedVector<IJK_Field_double, 3> * eulerian_facets_barycentre_;

  const IJK_Field_local_double * temperature_;
  const IJK_Field_local_double * temperature_ft_;
  const FixedVector<IJK_Field_double, 3> * velocity_;
  const FixedVector<IJK_Field_double, 3> * velocity_ft_;

  const FixedVector<IJK_Field_double, 3> * grad_T_elem_;
  const FixedVector<IJK_Field_double, 3> * hess_diag_T_elem_;
  const FixedVector<IJK_Field_double, 3> * hess_cross_T_elem_;

  // FIXME: Should I use DoubleTab instead ?
  const double * dr_base_ = 0;
  const DoubleVect* radial_coordinates_base_;

  const DoubleTab* radial_first_order_operator_raw_base_;
  const DoubleTab* radial_second_order_operator_raw_base_;
  const DoubleTab* radial_first_order_operator_base_;
  const DoubleTab* radial_second_order_operator_base_;
  const DoubleTab* radial_diffusion_matrix_base_;
  const DoubleTab* radial_convection_matrix_base_;

  double dr_=0.;
  const DoubleVect * radial_coordinates_;
  DoubleVect radial_coordinates_modified_;
  DoubleVect osculating_radial_coordinates_;
  DoubleTab radial_coordinates_cartesian_compo_;
  DoubleTab osculating_radial_coordinates_cartesian_compo_;

};

#endif /* IJK_One_Dimensional_Subproblem_included */
