/****************************************************************************
* Copyright (c) 2020, CEA
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
// File      : Pb_Thermohydraulique_sensibility.h
// Directory : $Sensitivity_analysis/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Pb_Thermohydraulique_sensibility_included
#define Pb_Thermohydraulique_sensibility_included

#include <Pb_Fluide_base.h>
#include <Navier_Stokes_std_sensibility.h>
#include <Convection_Diffusion_Temperature_sensibility.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Pb_Thermohydraulique_sensibility
//
// <Description of class Pb_Thermohydraulique_sensibility>
//  Cette classe represente un probleme de sensibilite thermohydraulique :
//      - Equations de Navier_Stokes_sensibility en regime laminaire
//        pour un fluide incompressible
//      - Equation sensibilite d'energie en regime laminaire
// .SECTION voir aussi
//     class Pb_Thermohydraulique, class Navier_Stokes_sensibility, class Convection_Diffusion_Temperature_sensibility
//
/////////////////////////////////////////////////////////////////////////////

class Pb_Thermohydraulique_sensibility : public Pb_Fluide_base
{

  Declare_instanciable( Pb_Thermohydraulique_sensibility ) ;

public:

  int nombre_d_equations() const;
  const Equation_base& equation(int) const ;
  Equation_base& equation(int);
  void associer_milieu_base(const Milieu_base& );
  int verifier();

protected:

  Navier_Stokes_std_sensibility eq_hydraulique;
  Convection_Diffusion_Temperature_sensibility eq_thermique;


};

#endif /* Pb_Thermohydraulique_sensibility_included */
