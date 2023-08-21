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

protected :
  int index_i_ = 0, index_j_ = 0, index_k_ = 0;
  int compo_connex_ = -1;
  int compo_group_ = -1;
  double distance_ = 0.;
  double curvature_ = 0.;
  double osculating_radius_ = 0.;
  double bubble_relative_velocity_ = 0.;
  double osculating_sphere_centre_[3] = {0., 0., 0.};
  double normal_vector_compo_[3] = {0., 0., 0.};
  /*
   * Several ways to calculate the tangential vector !
   * Either by considering a unique tangential vector (pure spherical)
   * Or two tangential vector (osculating sphere)
   */
  double first_tangential_vector_compo_[3] = {0., 0., 0.};
  double azymuthal_vector_compo[3] = {0., 0., 0.};
  double second_tangential_vector_compo_[3] = {0., 0., 0.};
  //
  double bary_compo_[3] = {0., 0., 0.};
  double temperature_gradient_compo_[3] = {0., 0., 0.};
  double temperature_hessian_compo_[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
  double normal_interfacial_gradient_ = 0;
  double normal_interfacial_gradient_compo_[3] = {0., 0., 0.};
};

#endif /* IJK_One_Dimensional_Subproblem_included */
