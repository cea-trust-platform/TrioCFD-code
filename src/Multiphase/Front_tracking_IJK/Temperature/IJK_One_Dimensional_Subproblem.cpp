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
// File      : IJK_One_Dimensional_Subproblem.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_One_Dimensional_Subproblem.h>

Implemente_instanciable( IJK_One_Dimensional_Subproblem, "IJK_One_Dimensional_Subproblem", Objet_U ) ;

// IJK_One_Dimensional_Subproblem::IJK_One_Dimensional_Subproblem() {}

Sortie& IJK_One_Dimensional_Subproblem::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& IJK_One_Dimensional_Subproblem::readOn( Entree& is )
{
  Objet_U::readOn( is );
  return is;
}

void IJK_One_Dimensional_Subproblem::associate_sub_problem_to_inputs(int i, int j, int k, int compo_connex,
                                                                     double distance, double curvature,
                                                                     ArrOfDouble facet_barycentre,
                                                                     ArrOfDouble normal_vector,
                                                                     double bubble_rising_velocity, ArrOfDouble bubble_rising_vector)
{
  associate_cell_ijk(i, j, k);
  associate_compos(compo_connex);
  associate_interface_related_parameters(distance, curvature, facet_barycentre, normal_vector);
  associate_rising_velocity(bubble_rising_velocity, bubble_rising_vector);
  /*
   *  Curvature is negative for a convex bubble
   *  but R should be positive in that case
   */
  if (fabs(curvature_) > DMINFLOAT)
    osculating_radius_ = fabs(2/curvature_);
  // FIXME: What happen with highly deformable bubbles (concave interface portions) ?
}

/*
 * TODO: Associate a basis to each subproblem
 * Use Rodrigues' rotation formula to determine ephi
 * Needs an axis of (rotation gravity_dir x relative_vectors)
 * and an angle (gravity_dir dot relative_vectors) / (norm(gravity_dir)*norm(relative_vectors))
 * ephi is determined in the gravity_align rising direction
 * 		 | gravity_dir
 * 		 |
 *   *****
 * ***   ***
 * **     **
 * ***   ***
 *   *****
 *     |
 *     |
 */


