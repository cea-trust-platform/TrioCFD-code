//TRUST_NO_INDENT
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

#ifndef IJK_SWITCH_FT_H
#define IJK_SWITCH_FT_H

#include <IJK_switch.h>
#include <IJK_Interfaces.h>
#include <IJK_Thermique.h>
#include <TRUST_List.h>
#include <init_forcage_THI.h>
#include <Force_sp.h>
#include <corrections_qdm.h>
#include <IJK_Thermal.h>
#include <IJK_Thermals.h>

class Switch_FT_double : public Switch_double
{
  Declare_instanciable_sans_constructeur(Switch_FT_double);
  friend class IJK_Thermals;
public:
  Switch_FT_double();

  // Copy ctor is forbidden:
  Switch_FT_double(const Switch_FT_double&x) {
    Cerr << "Erreur Switch_double(Switch_double&) interdit a cause du membre IJK_Interfaces" << finl;
    Process::exit();
  }
  const Switch_FT_double& operator=(const Switch_FT_double&) {
    Cerr << "Erreur Switch_double=() interdit a cause du membre IJK_Interfaces" << finl;
    Process::exit();
    return *this;
  }

protected:
  // Override from base class:
  void initialise() override;
  int init_thermique() override;
  int init_thermals() override;
  void prepare_run() override;
  void set_param(Param& param) override;
  void ecrire_fichier_reprise(const char *fichier_sauvegarde, const bool and_lata=true) override;
  void set_param_reprise(Param& param) override;
  void prepare_thermique(const Nom lata_name) override;
  void prepare_thermals(const Nom lata_name) override;
  void ecrire_header_lata(const Nom lata_name) override; // const;
  void compute_and_write_extra_fields(const Nom& lata_name, DoubleTab& coeff_i, IntTab Indice_i, 
                                              DoubleTab& coeff_j, IntTab Indice_j,
                                              DoubleTab& coeff_k, IntTab Indice_k) override;
  void compute_and_write_extra_fields_direct(SFichier& file,
																						 DoubleTab& coeff_i, IntTab Indice_i,
																						 DoubleTab& coeff_j, IntTab Indice_j,
																						 DoubleTab& coeff_k, IntTab Indice_k) override;

  // Le maillage des interfaces:
  IJK_Interfaces interfaces_;
  int old_ijk_splitting_ft_extension_;
  
  // GAB : gabriel.ramirez@cea.fr
  init_forcage_THI old_forcage_;
  init_forcage_THI new_forcage_;
  corrections_qdm old_qdm_corrections_;
  corrections_qdm new_qdm_corrections_;
  double vap_velocity_tmoy_ = 0.;
  double liq_velocity_tmoy_ = 0.;
  double qdm_source_ = 0.;
  double last_source_qdm_update_time_ = 0.;
  int offset_list_index_ = 0.;
  int size_listes_source_;
  ArrOfDouble liste_instants_;
  ArrOfDouble liste_vap_dl_;  // liste des v_v * dt
  ArrOfDouble liste_liq_dl_;  // liste des v_l * dt
  double reprise_v_target_ = 0.;
  int list_index_ = 0.;
  // GAB : gabriel.ramirez@cea.fr

  // Dealing with thermal aspects:
  LIST(IJK_Thermique) thermique_;
  IJK_Thermals thermals_;
};

#endif
