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
// File      : IJK_Thermals.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Thermals_included
#define IJK_Thermals_included

#include <IJK_Thermal.h>
#include <TRUSTList.h>
#include <TRUST_List.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_Thermals
//
// <Description of class IJK_Thermals>
//
/////////////////////////////////////////////////////////////////////////////

class IJK_FT_double;
class Switch_FT_double;

class IJK_Thermals : public LIST(IJK_Thermal)
{

  Declare_instanciable_sans_constructeur( IJK_Thermals ) ;

public :
  IJK_Thermals() { ; };
  IJK_Thermals(const IJK_FT_double& ijk_ft);
  void associer(const IJK_FT_double& ijk_ft);
  void associer_post(const IJK_FT_Post& ijk_ft_post);
  void associer_switch(const Switch_FT_double& ijk_ft_switch);
  void sauvegarder_temperature(Nom& lata_name);
  void sauvegarder_thermals(SFichier& fichier);
  void compute_timestep(double& dt_thermals, const double dxmin);
  void initialize(const IJK_Splitting& splitting, int& nalloc);
  int size_thermal_problem(Nom thermal_problem);
  void update_thermal_properties();
  void euler_time_step(const double timestep);
  void euler_rustine_step(const double timestep);
  void rk3_sub_step(const int rk_step, const double total_timestep, const double time);
  void rk3_rustine_sub_step(const int rk_step, const double total_timestep,
                            const double fractionnal_timestep, const double time);
  void posttraiter_tous_champs_thermal(Motcles& liste_post_instantanes_);
  void posttraiter_champs_instantanes_thermal(const Motcles& liste_post_instantanes,
                                              const char *lata_name,
                                              const int latastep,
                                              const double current_time,
                                              int& n);
  int init_thermals(const IJK_Splitting& splitting);
  void prepare_thermals(const char *lataname);
  int ghost_fluid_flag();
  void ecrire_fichier_reprise(SFichier& fichier, const char *lata_name);
  int get_probes_ghost_cells(int ghost_init);
protected :
  REF(IJK_FT_double) ref_ijk_ft_;
  REF(IJK_FT_Post) ref_ijk_ft_post_;
  REF(Switch_FT_double) ref_ijk_ft_switch_;

};

#endif /* IJK_Thermals_included */
