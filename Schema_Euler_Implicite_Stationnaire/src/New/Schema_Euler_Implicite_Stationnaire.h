/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        Schema_Euler_Implicite_Stationnaire.h
// Directory:   $TRUST_ROOT/src/Kernel/Schemas_Temps
// Version:     /main/23
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Schema_Euler_Implicite_Stationnaire_included
#define Schema_Euler_Implicite_Stationnaire_included

//////////////////////////////////////////////////////////
// .DESCRIPTION class Schema_Euler_Implicite_Stationnaire
//  Il herite de schema implicite base et porte le solveur Implicit_steady
// pour effectuer les Faire_un_pas_de_temps..

#include <Schema_Euler_Implicite.h>

class Probleme_Couple;

class Schema_Euler_Implicite_Stationnaire : public Schema_Euler_Implicite
{

  Declare_instanciable(Schema_Euler_Implicite_Stationnaire);

public :
  void set_param(Param& );
  virtual int mettre_a_jour();
  virtual int reprendre(Entree& );
  void ajouter_inertie(Matrice_Base& mat_morse,DoubleTab& secmem,const Equation_base& eqn) const;
  void initialize();

  void mettre_a_jour_dt_stab();
  bool corriger_dt_calcule(double& dt_calc) const;
  void calculer_pas_de_temps_local_pb();

  // pour implicite ajoute l'inertie a la matrice et au scd membre avec pas de temps local
  // virtual void ajouter_inertie_pas_temps_locaux(Matrice_Base& mat_morse,DoubleTab& secmem,const Equation_base& eqn) const;
  /*
  inline int steady() const
  {
    return steady_;
  };
  inline double steady_security_facteur() const
  {
    return steady_security_facteur_;
  };
  inline double steady_global_dt() const
  {
    return steady_global_dt_;
  };
  inline void fix_dt_locaux(DoubleVect& dt_locaux)
  {
    dt_locaux_= dt_locaux;
    max_dt_locaux_= dt_locaux.mp_max_abs_vect();
    dt_locaux_.echange_espace_virtuel();
  }
  inline int nb_steps_div0_imposed() const
  {
    return nb_steps_div0_imposed_;
  };
  inline int div0_imposed() const
  {
    return div0_imposed_;
  };
  inline void set_div0_imposed(int value)
  {
    div0_imposed_ = value;
  };
  */
protected:
  double max_dt_locaux_;                     // max  of dt_locaux_
  //  int steady_;                          // 1 if one wants to use a local time step procedure
  double steady_security_facteur_;  // parameter used in the local time step calculation procedure
  double steady_global_dt_;         // in the local time step procedure, an overall time step is used
  int nb_steps_div0_imposed_;         //impose div u = 0 each nb_steps_div0_imposed
  int div0_imposed_;                  // 1 if one wants to impose div u = 0 only every nb_steps_div0_imposed time steps
};

#endif

