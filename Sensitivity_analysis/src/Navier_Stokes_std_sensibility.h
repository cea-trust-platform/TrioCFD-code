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
// File      : Navier_Stokes_std_sensibility.h
// Directory : $Sensitivity_analysis/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Navier_Stokes_std_sensibility_included
#define Navier_Stokes_std_sensibility_included

#include <Navier_Stokes_std.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Navier_Stokes_std_sensibility
//
// <Description of class Navier_Stokes_std_sensibility>
//
/////////////////////////////////////////////////////////////////////////////

class Navier_Stokes_std_sensibility : public Navier_Stokes_std
{

  Declare_instanciable( Navier_Stokes_std_sensibility ) ;

public :
  void set_param(Param& param);
  int lire_motcle_non_standard(const Motcle& mot, Entree& is);
  void associate_evaluator_field(const Nom& one_name_state_pb,const Motcle& one_name_state_field);
  void update_evaluator_field(const Nom& one_name_state_pb,const Motcle& one_name_state_field);
  virtual  void mettre_a_jour(double temps);
  const DoubleTab& get_state_field() const;
  const Champ_Inc_base& get_state() const;
  const Motcle& get_uncertain_variable_name() const;

protected :
  REF(Champ_Inc_base) state_field;  //Reference to the unknown field of the state problem
  Nom name_state_pb;                      //name of the problem state
  Motcle name_state_field;                 //name of the unknown field of the state problem
  Motcle uncertain_var;     				//name of the unknown field of the uncertain variable

};

#endif /* Navier_Stokes_std_sensibility_included */
