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
  void associate_cell_ijk(int i, int j, int k) { index_i_ = i; index_j_=j; index_k_=k; };
  void associate_compos(int compo_connex) { compo_connex_ = compo_connex; };
  void associate_compos(int compo_connex, int compo_group) { compo_connex_ = compo_connex; compo_group_ = compo_group; };
  void associate_interface_related_parameters(double distance, double curvature, ArrOfDouble facet_barycentre, ArrOfDouble normal_vector)
  {
    distance_ = distance;
    curvature_ = curvature;
    facet_barycentre_ = facet_barycentre;
    normal_vector_compo_ = normal_vector;
  };
  void associate_rising_velocity(double bubble_rising_velocity, ArrOfDouble bubble_rising_vector)
  {
    bubble_rising_velocity_ = bubble_rising_velocity;
    bubble_rising_vector_ = bubble_rising_vector;
  };
  void associate_sub_problem_to_inputs(int i, int j, int k, int compo_connex,
                                       double distance, double curvature,
                                       ArrOfDouble facet_barycentre,
                                       ArrOfDouble normal_vector,
                                       double bubble_rising_velocity, ArrOfDouble bubble_rising_vector);
protected :
  int index_i_ = 0, index_j_ = 0, index_k_ = 0;
  int compo_connex_ = -1;
  int compo_group_ = -1;
  double distance_ = 0.;
  double curvature_ = 0.;
  double osculating_radius_ = 0.;
  ArrOfDouble normal_vector_compo_ = ArrOfDouble(3);

  double bubble_rising_velocity_ = 0.;
  ArrOfDouble bubble_rising_vector_ = ArrOfDouble(3);

  ArrOfDouble osculating_sphere_centre_ = ArrOfDouble(3);
  /*
   * Several ways to calculate the tangential vector !
   * Either by considering a unique tangential vector (pure spherical)
   * Or two tangential vectors (osculating sphere)
   */
  ArrOfDouble first_tangential_vector_compo_ = ArrOfDouble(3);
  ArrOfDouble azymuthal_vector_compo = ArrOfDouble(3);
  ArrOfDouble second_tangential_vector_compo_ = ArrOfDouble(3);
  //
  ArrOfDouble facet_barycentre_ = ArrOfDouble(3);
  ArrOfDouble temperature_gradient_compo_ = ArrOfDouble(3);
  DoubleTab temperature_hessian_compo_ = DoubleTab(3,3);
  double normal_interfacial_gradient_ = 0;
  ArrOfDouble normal_interfacial_gradient_compo_ = ArrOfDouble(3);
};

#endif /* IJK_One_Dimensional_Subproblem_included */
