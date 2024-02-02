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
#include <System.h>

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

  Declare_instanciable( IJK_Thermals ) ;

public :
  IJK_Thermals(const IJK_FT_double& ijk_ft);
  void set_fichier_reprise(const char *lataname);
  const Nom& get_fichier_reprise();
  void associer(const IJK_FT_double& ijk_ft);
  void associer_post(const IJK_FT_Post& ijk_ft_post);
  void associer_switch(const Switch_FT_double& ijk_ft_switch);
  void associer_interface_intersections(const Intersection_Interface_ijk_cell& intersection_ijk_cell_,
                                        const Intersection_Interface_ijk_face& intersection_ijk_face_);
  void retrieve_ghost_fluid_params();
  void sauvegarder_temperature(Nom& lata_name, const int& stop);
  void sauvegarder_thermals(SFichier& fichier);
  void compute_timestep(double& dt_thermals, const double dxmin);
  void initialize(const IJK_Splitting& splitting, int& nalloc);
  void recompute_temperature_init();
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
  int init_switch_thermals(const IJK_Splitting& splitting);
  void prepare_thermals(const char *lataname);
  int ghost_fluid_flag();
  void ecrire_fichier_reprise(SFichier& fichier, const char *lata_name);
  void compute_ghost_cell_numbers_for_subproblems(const IJK_Splitting& splitting, int ghost_init);
  int get_probes_ghost_cells(int ghost_init);

  void update_intersections();
  void clean_ijk_intersections();

  void compute_eulerian_distance();
  void compute_eulerian_curvature();
  void compute_eulerian_curvature_from_interface();
  void compute_eulerian_distance_curvature();

  void set_latastep_reprise(const bool stop);
  void thermal_subresolution_outputs(const int& dt_post_thermals_probes=0);
  int get_disable_post_processing_probes_out_files() const;
  double get_modified_time();
  void get_rising_velocities_parameters(int& compute_rising_velocities,
                                        int& fill_rising_velocities);
  void create_folders_for_probes();
  void create_folders(Nom folder_name_base);
  void set_first_step_thermals_post(int& first_step_thermals_post);
  void set_post_pro_first_call() { post_pro_first_call_ = 1; } ;
  void set_temperature_ini();
  void recompute_interface_smoothing();
  void compute_new_thermal_field(Switch_FT_double& switch_double_ft,
                                 const IJK_Splitting& new_mesh,
                                 const Nom& lata_name,
                                 DoubleTab& coeff_i,
                                 IntTab Indice_i,
                                 DoubleTab& coeff_j,
                                 IntTab Indice_j,
                                 DoubleTab& coeff_k,
                                 IntTab Indice_k);


protected :
  REF(IJK_FT_double) ref_ijk_ft_;
  REF(IJK_FT_Post) ref_ijk_ft_post_;
  REF(Switch_FT_double) ref_ijk_ft_switch_;
  REF(Intersection_Interface_ijk_cell) ref_intersection_ijk_cell_;
  REF(Intersection_Interface_ijk_face) ref_intersection_ijk_face_;

  IJK_Ghost_Fluid_Fields ghost_fluid_fields_;

  int post_pro_first_call_ = 0;

  System make_dir_for_out_files_;
  LIST(Nom) thermal_rank_folder_;
  Nom overall_bubbles_quantities_folder_;
  Nom interfacial_quantities_thermal_probes_folder_;
  Nom local_quantities_thermal_probes_folder_;
  Nom local_quantities_thermal_probes_time_index_folder_;
  int ini_folder_out_files_ = 0;

  bool is_diphasique_=false;
  std::vector<int> lata_step_reprise_ini_;
  std::vector<int> lata_step_reprise_;
};

#endif /* IJK_Thermals_included */
