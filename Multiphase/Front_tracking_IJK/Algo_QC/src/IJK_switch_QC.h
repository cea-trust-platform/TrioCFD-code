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

#ifndef IJK_SWITCH_QC_H
#define IJK_SWITCH_QC_H

#include <IJK_switch.h>


//
// WARNING: DO NOT EDIT THIS FILE! Only edit the template file IJK_switch_QC.h.P
//
class Switch_QC_double : public Switch_double
{
  Declare_instanciable_sans_constructeur(Switch_QC_double);

public:
  Switch_QC_double();

protected:
  // Overrides from base class
  virtual void initialise();
  virtual int  allocate_fields(double& sz_arr);
  virtual void prepare_run();
  virtual void set_param(Param& param);
  virtual void set_param_reprise(Param& param);
  virtual void extra_write_rho_for_qc(SFichier& file,
                                      DoubleTab& coeff_i, DoubleTab& coeff_j, DoubleTab& coeff_k,
                                      IntTab Indice_i, IntTab Indice_j, IntTab Indice_k);
  virtual void ecrire_fichier_reprise(const char *fichier_sauvegarde, const bool and_lata);
  virtual void compute_and_write_extra_fields(const Nom& lata_name, DoubleTab& coeff_i, IntTab Indice_i, 
                                              DoubleTab& coeff_j, IntTab Indice_j,
                                              DoubleTab& coeff_k, IntTab Indice_k);
  virtual void compute_and_write_extra_fields_direct(SFichier& file,
                                                     DoubleTab& coeff_i, IntTab Indice_i, 
                                                     DoubleTab& coeff_j, IntTab Indice_j,
                                                     DoubleTab& coeff_k, IntTab Indice_k);

  void remplir_ghost();  // UNUSED?

  Nom fichier_old_rho_;
  Nom fichier_new_rho_;

  int timestep_reprise_rho_;

  /// CL.
  double rho_kmin_;
  double rho_kmax_;

  double Cp_gaz_;
  double gamma_;
  double temperature_paroi_min_;
  double temperature_paroi_max_;
  double P_thermodynamique_;
  double constante_specifique_gaz_;
};

#endif
