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
/////////////////////////////////////////////////////////////////////////////
//
// File      : init_forcage.h
// Directory : $IJK_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef init_forcage_THI_included
#define init_forcage_THI_included

#include <FixedVector.h>
#include <IJK_Field.h>
#include <Objet_U.h>
#include <init_forcage_THI.h>
#include <string>
#include <iostream>
#include <math.h>
#include <fftw3.h>

#include <Force_ph.h>
#include <Force_sp.h>
#include <Random_process.h>


class init_forcage_THI : public Objet_U
{

  Declare_instanciable_sans_constructeur( init_forcage_THI ) ;

public :

  init_forcage_THI();
  void compute_initial_chouippe(int nproc_tot,
                                const IJK_Grid_Geometry& my_geom,
                                int my_ni, int my_nj, int my_nk,
                                const IJK_Splitting& my_splitting,
                                Nom nom_sauvegarde);
  // FixedVector<IJK_Field_double, 3> v_);

  void compute_THI_force(const int time_iteration,
                         const double tstep,
                         const double current_time,
                         const IJK_Splitting& my_splitting
                        );

  int get_type_forcage();
  int get_facteur_forcage();
  int get_forced_advection();
  int activate_forcage(const int current_time_step, const double current_time);
  int get_semi_gen();
  ArrOfDouble get_b_flt();
  FixedVector<IJK_Field_double, 3> get_force_ph();
  FixedVector<IJK_Field_double, 3>& get_force_ph2();

  void update_advection_velocity(ArrOfDouble& value);
  void update_advection_length(double dt);

protected :
  int type_forcage;
  int facteur_forcage_;
  int forced_advection_;
  // GAB : sarebbe più logico in vector3 ma non so come dare un vector3 in un jdd...
  ArrOfDouble advection_velocity_;
  ArrOfDouble advection_length_;
  double time_to_be_del_;
  int forcage_ts_stop;
  double forcage_t_stop;
  int mode_min;
  int mode_max;
  double amplitude;
  double eps_etoile;
  double tL;

  int i_offset;
  int j_offset;
  int k_offset;

//  int random_fixed_;
  Force_sp f_sp_THI;
  Force_ph f_ph_THI;
  Random_process random_;

};



#endif /* Init_spectral_included */
