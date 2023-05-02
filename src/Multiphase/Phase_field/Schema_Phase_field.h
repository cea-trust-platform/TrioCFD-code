/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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

#ifndef Schema_Phase_field_included
#define Schema_Phase_field_included

#include <Convection_Diffusion_Phase_field.h>

/*! @brief classe Schema_Phase_field
 *
 * @sa Schema_Temps_base
 */
class Schema_Phase_field: public Schema_Temps_base
{

  Declare_instanciable(Schema_Phase_field);

public :

  ////////////////////////////////
  //                            //
  // Caracteristiques du schema //
  //                            //
  ////////////////////////////////

  int nb_valeurs_temporelles() const override;
  int nb_valeurs_futures() const override;
  double temps_futur(int i) const override;
  double temps_defaut() const override;

  /////////////////////////////////////////
  //                                     //
  // Fin des caracteristiques du schema  //
  //                                     //
  /////////////////////////////////////////

  void initialize() override;
  bool initTimeStep(double dt) override;
  int faire_un_pas_de_temps_eqn_base(Equation_base&) override;
  bool iterateTimeStep(bool& converged) override;
  virtual int faire_un_pas_de_temps_C_D_Phase_field(Convection_Diffusion_Phase_field&);
  virtual int premier_dt(Convection_Diffusion_Phase_field& eq_c);
  virtual int deuxieme_dt(Convection_Diffusion_Phase_field& eq_c);
  inline const DoubleTab& valeur_temps_intermediaire() const;
  //virtual void lire(const Motcle&, Entree&);
  void set_param(Param& titi) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  int mettre_a_jour() override;
  bool corriger_dt_calcule(double&) const override;
  void completer() override;
  void changer_temps_courant(const double) override;
  int stop() const override;
  void imprimer(Sortie&) const override;

  void associer_pb(const Probleme_base&) override;

protected:
  int nb_iterations_relax_;

  Schema_Temps sch2;
  Schema_Temps sch3;
};

#endif
