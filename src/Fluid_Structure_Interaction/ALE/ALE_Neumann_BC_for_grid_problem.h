/****************************************************************************
* Copyright (c) 2021, CEA
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
// File      : ALE_Neumann_BC_for_grid_problem.h
// Directory : $ALE_MODEL_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef ALE_Neumann_BC_for_grid_problem_included
#define ALE_Neumann_BC_for_grid_problem_included

#include <Interprete_geometrique_base.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Projection_ALE_boundary
//
// <Description of class Projection_ALE_boundary>
// This class represents the projection of a function defined on a boundary in case of ALE calculation
// In general it is the modal deformation
// The fluid force is project on this function. For the modal deformation it results in the modal fluid force acting on a ALE boundary
/////////////////////////////////////////////////////////////////////////////

class ALE_Neumann_BC_for_grid_problem : Interprete_geometrique_base
{

  Declare_instanciable(ALE_Neumann_BC_for_grid_problem) ;

public :

protected :
  Entree& interpreter_(Entree&) override;
};


#endif /* ALE_Neumann_BC_for_grid_problem_included */
