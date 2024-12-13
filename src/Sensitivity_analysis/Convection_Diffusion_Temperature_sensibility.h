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
// File      : Convection_Diffusion_Temperature_sensibility.h
// Directory : $Sensitivity_analysis/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Convection_Diffusion_Temperature_sensibility_included
#define Convection_Diffusion_Temperature_sensibility_included

#include <Convection_Diffusion_Temperature.h>

/*! @brief : class Convection_Diffusion_Temperature_sensibility
 *
 *  <Description of class Convection_Diffusion_Temperature_sensibility>
 *
 *
 *
 */
class Convection_Diffusion_Temperature_sensibility : public Convection_Diffusion_Temperature
{

  Declare_instanciable( Convection_Diffusion_Temperature_sensibility ) ;

public :
  void set_param(Param& param) override;
  int lire_motcle_non_standard(const Motcle& mot, Entree& is) override;
  void associate_evaluator_field(const Nom& one_name_state_pb,const Motcle& one_name_state_field);
  void update_evaluator_field(const Nom& one_name_state_pb,const Motcle& one_name_state_field);
  void mettre_a_jour(double temps) override;
  const DoubleTab& get_velocity_state_field() const;
  const DoubleTab& get_temperature_state_field() const;
  const Champ_Inc_base& get_velocity_state() const;
  const Champ_Inc_base& get_temperature_state() const;
  const Motcle& get_uncertain_variable_name() const;
  const double& get_poly_chaos_value() const;


protected :
  OBS_PTR(Champ_Inc_base) velocity_state_field;  //Reference to the unknown velocity field of the state problem
  OBS_PTR(Champ_Inc_base) temperature_state_field;  //Reference to the unknown temperature field of the state problem
  Nom name_state_pb;                      //name of the problem state
  Motcle name_velocity_state_field;         //name of the unknown velocity field of the state problem
  Motcle name_temperature_state_field;      //name of the unknown temperature field of the state problem
  Motcle uncertain_var;     				//name of the unknown field of the uncertain variable
  double poly_chaos;               //It is the method that we will use to study the sensitivity of the Convection_Diffusion_Temperature_sensibility

};

#endif /* Convection_Diffusion_Temperature_sensibility_included */
