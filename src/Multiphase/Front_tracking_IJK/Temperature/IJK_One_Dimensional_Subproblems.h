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
// File      : IJK_One_Dimensional_Subproblems.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_One_Dimensional_Subproblems_included
#define IJK_One_Dimensional_Subproblems_included

#include <IJK_One_Dimensional_Subproblem.h>
#include <TRUSTList.h>
#include <TRUST_List.h>


/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_One_Dimensional_Subproblems
//
// <Description of class IJK_One_Dimensional_Subproblems>
//
/////////////////////////////////////////////////////////////////////////////

class IJK_FT_double;
class Switch_FT_double;

class IJK_One_Dimensional_Subproblems : public LIST(IJK_One_Dimensional_Subproblem)
{

  Declare_instanciable_sans_constructeur(IJK_One_Dimensional_Subproblems);

public :
  IJK_One_Dimensional_Subproblems() { ; };
  IJK_One_Dimensional_Subproblems(const IJK_FT_double& ijk_ft);
  void associer(const IJK_FT_double& ijk_ft) { ref_ijk_ft_ = ijk_ft; };
  void add_subproblems(int n);
  void associate_sub_problem_to_inputs(int i, int j, int k,
                                       const IJK_Field_double& eulerian_compo_connex,
                                       const IJK_Field_double& eulerian_distance,
                                       const IJK_Field_double& eulerian_curvature,
                                       FixedVector<IJK_Field_double, 3> eulerian_normal_vectors,
                                       FixedVector<IJK_Field_double, 3> eulerian_facets_barycentre,
                                       ArrOfDouble rising_velocities,
                                       DoubleTab rising_vectors);
  const int& get_subproblems_counter() const
  {
    return subproblems_counter_;
  }

protected :
  int max_subproblems_ = 0;
  int subproblems_counter_ = 0;
  REF(IJK_FT_double) ref_ijk_ft_;

};

#endif /* IJK_One_Dimensional_Subproblems_included */
