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
// File      : IJK_SolveSys_FD_thermal.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_SolveSys_FD_thermal.h>
#include <Param.h>
#include <Solv_Gmres.h>
#include <Solv_Petsc.h>

Implemente_instanciable( IJK_SolveSys_FD_thermal, "IJK_SolveSys_FD_thermal", SolveurSys ) ;

Sortie& IJK_SolveSys_FD_thermal::printOn( Sortie& os ) const
{
  SolveurSys::printOn( os );
  return os;
}

Entree& IJK_SolveSys_FD_thermal::readOn( Entree& is )
{
  Param param(que_suis_je());
  Nom solver_name;
  // Put any string to read the following word
  param.ajouter("thermal_fd_solver",&solver_name, Param::REQUIRED);
  param.lire_sans_accolade(is);

  if (solver_name == "direct")
    cast_direct_solver_by_default();
  else if (solver_name == "iterative")
    cast_iterative_solver_by_default();
  else
    {
      Nom type_solv_sys("Solv_");
      type_solv_sys+=solver_name;
      typer(type_solv_sys);
      nommer(Nom("thermal_fd_solver_") + solver_name);
      is >> valeur();
    }
  return is;
}

void IJK_SolveSys_FD_thermal::cast_iterative_solver_by_default()
{
  typer(iterative_solver_by_default_);
  nommer(Nom("thermal_fd_solver_") + iterative_solver_by_default_);
  /*
   * Fill the required params (example)
   */
  Solv_Gmres& gmres_solver = ref_cast(Solv_Gmres, valeur());
  gmres_solver.set_seuil(default_seuil);
}

void IJK_SolveSys_FD_thermal::cast_direct_solver_by_default()
{
  typer(direct_solver_by_default_);
  nommer(Nom("thermal_fd_solver_") + direct_solver_by_default_ + Nom("_") + petsc_solver_by_default_);
  Solv_Petsc& petsc_solver = ref_cast(Solv_Petsc, valeur());
  Nom solver_params("");
  Nom left_bracket("{");
  Nom right_bracket("}");
  Nom spacing(" ");
  solver_params += petsc_solver_by_default_ + spacing + left_bracket + spacing;
  /*
   * Add params by default
   */
  solver_params += Nom("nb_it_max") + spacing + Nom(default_nb_iter_max) + spacing;
  solver_params += right_bracket;
  Cerr << "Thermal F-D Solver:" << solver_params << finl;
  istringstream solver_params_istringstream(solver_params.getChar());
  istream& solver_params_istream = solver_params_istringstream;
  Entree solver_params_entry(solver_params_istream);
  /*
   * Fake a readOn
   */
  //solver_params_entry >> solver_params;
  petsc_solver.nommer(le_nom());
  solver_params_entry >> petsc_solver;
}
