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
// File      : ALE_Neumann_BC_for_grid_problem.cpp
// Directory : $BEAM_MODEL_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Domaine_ALE.h>
#include "ALE_Neumann_BC_for_grid_problem.h"

Implemente_instanciable(ALE_Neumann_BC_for_grid_problem, "ALE_Neumann_BC_for_grid_problem", Interprete_geometrique_base ) ;
//XD  ALE_Neumann_BC_for_grid_problem interprete ALE_Neumann_BC_for_grid_problem 0 block to indicates the names of the boundary with Neumann BC for the grid problem. By default, in the ALE grid problem, we impose a homogeneous Dirichelt-type BC on the fix boundary. This option allows you to impose also Neumann-type BCs on certain boundary.
//XD  attr dom ref_domaine dom 0 Name of domain.
//XD attr bloc bloc_lecture bloc 0 between the braces, you must specify the numbers of the mobile borders then list these mobile borders.  NL2 Example:  ALE_Neumann_BC_for_grid_problem dom_name  { 1  boundary_name }

Sortie& ALE_Neumann_BC_for_grid_problem::printOn( Sortie& os ) const
{
  return Interprete::printOn( os );
}

Entree& ALE_Neumann_BC_for_grid_problem::readOn( Entree& is )
{
  return Interprete::readOn( is );
}

Entree& ALE_Neumann_BC_for_grid_problem::interpreter_(Entree& is)
{
  associer_domaine(is);
  Domaine_ALE& dom=ref_cast(Domaine_ALE, domaine());
  dom.reading_ALE_Neumann_BC_for_grid_problem(is);
  return is;
}
