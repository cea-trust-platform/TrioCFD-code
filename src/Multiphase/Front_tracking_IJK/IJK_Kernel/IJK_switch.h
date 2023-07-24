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

#ifndef IJK_SWITCH_H
#define IJK_SWITCH_H

#include <IJK_Splitting.h>
#include <FixedVector.h>
#include <Interprete.h>
#include <Linear_algebra_tools.h>
#include <Param.h>
#include <Interprete_bloc.h>
#include <SFichier.h>
#include <stat_counters.h>
#include <IJK_Lata_writer.h>
#include <IJK_Navier_Stokes_tools.h>
#include <communications.h>
#include <LecFicDiffuse_JDD.h>
#include <IJK_Field.h>
#include <TRUSTList.h>


/**
 * Interprete permettant d'interpoler vitesse et rho d'un maillage sur un autre
 */
class Switch_double : public Interprete
{
  Declare_base(Switch_double);

public:
  Entree& interpreter(Entree&) override;
  void run();

protected:
  virtual void initialise();
  virtual int allocate_fields(double& sz_arr);
  virtual void prepare_run();
  virtual void set_param(Param& param);
  virtual void ecrire_fichier_reprise(const char *fichier_sauvegarde, const bool and_lata=true) = 0;
  virtual void set_param_reprise(Param& param);
  virtual void lire_fichier_reprise(const char *fichier_reprise);
  virtual void ecrire_header_lata(const Nom lata_name); // const;
  virtual void compute_and_write_extra_fields(const Nom& lata_name, DoubleTab& coeff_i, IntTab Indice_i,
                                              DoubleTab& coeff_j, IntTab Indice_j,
                                              DoubleTab& coeff_k, IntTab Indice_k) = 0;
  virtual void compute_and_write_extra_fields_direct(SFichier& file,
                                                     DoubleTab& coeff_i, IntTab Indice_i,
                                                     DoubleTab& coeff_j, IntTab Indice_j,
                                                     DoubleTab& coeff_k, IntTab Indice_k) = 0;

  // overriden in FT - to be brought here at some point once merge with thermic BC is done ...
  virtual int init_thermique()  { return 0; };
  virtual void prepare_thermique(const Nom lata_name) {};
  virtual int init_thermals()  { return 0; };
  virtual void prepare_thermals(const Nom lata_name) {};

  void write_velocity(const Nom lata_name) const;
  void calculer_coords(const IJK_Splitting::Localisation loc);
  void calculer_coords_elem();
  void calculer_coords_Vi(const int dir);
  void remplir_gost();
  void calculer_coeff(DoubleTab& coeff_i,IntTab& Indice_i,
                      DoubleTab& coeff_j,IntTab& Indice_j ,
                      DoubleTab& coeff_k,IntTab& Indice_k ) ;
  void switch_vit(DoubleTab coeff_i, IntTab Indice_i,
                  DoubleTab coeff_j ,IntTab Indice_j,
                  DoubleTab coeff_k ,IntTab Indice_k,
                  const int dir);//faces

  void switch_scalar_field(const IJK_Field_double& oldf, IJK_Field_double& newf,
                           DoubleTab coeff_i, IntTab Indice_i,
                           DoubleTab coeff_j ,IntTab Indice_j,
                           DoubleTab coeff_k ,IntTab Indice_k) const;

  void switch_scalar_field_direct(SFichier& binary_file,
                                  const IJK_Field_double& fld,
                                  DoubleTab coeff_i, IntTab Indice_i,
                                  DoubleTab coeff_j ,IntTab Indice_j,
                                  DoubleTab coeff_k ,IntTab Indice_k);


  void switch_vit_direct(SFichier& binary_file);

  // maillages
  Nom old_mesh_name_, new_mesh_name_;
  IJK_Splitting new_mesh_;
  IJK_Splitting old_mesh_;

  // direct_write_ = 1 => ecrit les fichiers a la volee sans alouer la memoire
  // direct_write_ = 0 => ancien algo avec allocation de new_velocity_
  int direct_write_;

  int perio_k_;
  // vitesses
  FixedVector<IJK_Field_double, 3> new_velocity_;
  FixedVector<IJK_Field_double, 3> old_velocity_;

  // masse volumiques
  IJK_Field_double new_rho_;
  IJK_Field_double old_rho_;

  // coord en x,y et z;
  ArrOfDouble_with_ghost old_x_;
  ArrOfDouble_with_ghost new_x_;
  ArrOfDouble_with_ghost old_y_;
  ArrOfDouble_with_ghost new_y_;
  ArrOfDouble_with_ghost old_z_;
  ArrOfDouble_with_ghost new_z_;

  // Pour alleger le code.
  int old_ni_, old_nj_, old_nk_;

  int new_ni_, new_nj_, new_nk_;

  // pour compatibiliter des fonctions.
  double current_time_;
  double terme_source_acceleration_;

  // nom des fichiers.
  Nom nom_sauvegarde_;
  Nom nom_reprise_;

  /*
  // GAB : gabriel.ramirez@cea.fr
  init_forcage_THI forcage_;
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
  */

  // Le jeu de donnees doit fournir soit des fichiers de reprise: ..
  Nom fichier_old_vitesse_;
  Nom fichier_new_vitesse_;

  int timestep_reprise_vitesse_;
};

#endif
