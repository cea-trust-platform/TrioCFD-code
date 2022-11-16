/****************************************************************************
* Copyright (c) 2019, CEA
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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Solver_moving_mesh_ALE.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/ALE/src
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////

#include <Domaine_ALE.h>
#include <Solver_moving_mesh_ALE.h>


Implemente_instanciable(Solver_moving_mesh_ALE,"Solver_moving_mesh_ALE",Interprete_geometrique_base);
//XD Solver_moving_mesh_ALE  interprete  Solver_moving_mesh_ALE 0 Solver used to solve the system giving the  mesh velocity for the ALE (Arbitrary Lagrangian-Eulerian) framework.
//XD attr  dom ref_domaine dom 0 Name of domain.
//XD  attr bloc bloc_lecture bloc 0   Example:  { PETSC GCP { precond ssor { omega 1.5 } seuil 1e-7 impr }  }


/*! @brief Simple appel a: Interprete::printOn(Sortie&)
 *
 *     Imprime l'interprete sur un flot de sortie
 *
 * @param (Sortie& os) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Solver_moving_mesh_ALE::printOn(Sortie& os) const
{
  return Interprete::printOn(os);
}


/*! @brief Simple appel a: Interprete::readOn(Entree&)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Solver_moving_mesh_ALE::readOn(Entree& is)
{
  return Interprete::readOn(is);
}


/*! @brief Fonction principale de l'interprete Solver_moving_mesh_ALE: Read the solver used to solve the system giving the moving mesh velocity
 *
 */
Entree& Solver_moving_mesh_ALE::interpreter_(Entree& is)
{
  associer_domaine(is);
  Domaine_ALE& dom=ref_cast(Domaine_ALE, domaine());
  dom.reading_solver_moving_mesh_ALE(is);
  return is;
}
