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
// File      : IJK_Thermal.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Thermal_included
#define IJK_Thermal_included

#include <IJK_Thermal_base.h>
#include <TRUST_Deriv.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_Thermal
//
// <Description of class IJK_Thermal>
//
/////////////////////////////////////////////////////////////////////////////

class IJK_FT_double;
class Switch_FT_double;

class IJK_Thermal : public DERIV(IJK_Thermal_base)
{
  Declare_instanciable( IJK_Thermal );

public :
  /*
   * Getters
   */
  inline Nom& get_thermal_problem_type() { return thermal_problem_type_; };
  inline int& get_thermal_rank() { return thermal_rank_; };
  inline Motcles& get_thermal_words() { return thermal_words_; };
  inline Motcles& get_thermal_suffix() { return lata_suffix_; };
  inline const IJK_Field_double& get_temperature() const { return valeur().get_temperature(); };
  inline const IJK_Field_double& get_temperature_before_extrapolation() const { return valeur().get_temperature_before_extrapolation(); }
  inline const IJK_Field_double& get_temperature_cell_neighbours() const { return valeur().get_temperature_cell_neighbours(); }
  inline const IJK_Field_double& get_temperature_cell_neighbours_debug() const { return valeur().get_temperature_cell_neighbours_debug(); }
  inline const IJK_Field_int& get_cell_neighbours_corrected() const { return valeur().get_cell_neighbours_corrected(); }
  inline const IJK_Field_double& get_neighbours_temperature_colinearity_weighting() const { return valeur().get_neighbours_temperature_colinearity_weighting(); }
  inline IJK_Field_double& get_temperature_ft() { return valeur().get_temperature_ft(); }
  inline const IJK_Field_double& get_temperature_ana() const { return valeur().get_temperature_ana(); };
  inline const IJK_Field_double& get_ecart_t_ana() const { return valeur().get_ecart_t_ana(); }
  inline const FixedVector<IJK_Field_double, 3>& get_grad_T() const { return valeur().get_grad_T(); }
  inline const IJK_Field_double& get_div_lambda_grad_T() const { return valeur().get_div_lambda_grad_T(); }
  inline const IJK_Field_double& get_u_T_convective() const { return valeur().get_u_T_convective(); }
  inline const IJK_Field_double& get_eulerian_distance_ft() const { return valeur().get_eulerian_distance_ft(); }
  inline const IJK_Field_double& get_eulerian_curvature_ft() const { return valeur().get_eulerian_curvature_ft(); }
  inline const IJK_Field_double& get_interfacial_area_ft() const { return valeur().get_interfacial_area_ft(); }
  inline const IJK_Field_double& get_grad_T_interface_ft() const { return valeur().get_grad_T_interface_ft(); }
  inline const IJK_Field_double& get_eulerian_compo_connex_ft() const { return valeur().get_eulerian_compo_connex_ft(); }
  inline const IJK_Field_double& get_eulerian_compo_connex_ghost_ft() const { return valeur().get_eulerian_compo_connex_ghost_ft(); }
  inline const IJK_Field_double& get_eulerian_compo_connex_from_interface_ft() const { return valeur().get_eulerian_compo_connex_from_interface_ft(); }
  inline const IJK_Field_double& get_eulerian_compo_connex_ns() const { return valeur().get_eulerian_compo_connex_ns(); }
  inline const IJK_Field_double& get_eulerian_compo_connex_ghost_ns() const { return valeur().get_eulerian_compo_connex_ghost_ns(); }
  inline const IJK_Field_double& get_eulerian_compo_connex_from_interface_ns() const { return valeur().get_eulerian_compo_connex_from_interface_ns(); }
  inline const IJK_Field_int& get_eulerian_compo_connex_int_from_interface_ns() const { return valeur().get_eulerian_compo_connex_int_from_interface_ns(); }
  inline const IJK_Field_double& get_eulerian_distance_ns() const { return valeur().get_eulerian_distance_ns(); }
  inline const IJK_Field_double& get_eulerian_curvature_ns() const { return valeur().get_eulerian_curvature_ns(); }
  inline const IJK_Field_double& get_interfacial_area_ns() const { return valeur().get_interfacial_area_ns(); }
  inline const IJK_Field_double& get_grad_T_interface_ns() const { return valeur().get_grad_T_interface_ns(); }
  inline const IJK_Field_double& get_eulerian_rising_velocities() const {return valeur().get_eulerian_rising_velocities(); }
  inline const FixedVector<IJK_Field_double, 3>& get_bary() const { return valeur().get_bary(); }
  inline const FixedVector<IJK_Field_double, 3>& get_gradient_temperature_elem() { return valeur().get_gradient_temperature_elem(); }
  inline const FixedVector<IJK_Field_double, 3>& get_normal_vector_ns() const { return valeur().get_normal_vector_ns(); }
  inline const FixedVector<IJK_Field_double, 3>& get_normal_vector_ns_normed() const { return valeur().get_normal_vector_ns_normed();}
  inline const FixedVector<IJK_Field_double, 3>& get_normal_vector_ft() const { return valeur().get_normal_vector_ft(); }
  inline const FixedVector<IJK_Field_double, 3>& get_hessian_diag_temperature_elem() const { return valeur().get_hessian_diag_temperature_elem(); }
  inline const FixedVector<IJK_Field_double, 3>& get_hessian_cross_temperature_elem() const { return valeur().get_hessian_cross_temperature_elem(); }
  inline const FixedVector<IJK_Field_double,3>& get_cell_faces_corrected_diffusive() const { return valeur().get_cell_faces_corrected_diffusive(); }
  inline const FixedVector<IJK_Field_double,3>& get_cell_faces_corrected_convective() const { return valeur().get_cell_faces_corrected_convective(); }
  inline const FixedVector<IJK_Field_int,3>& get_cell_faces_corrected_bool() const { return valeur().get_cell_faces_corrected_bool(); }
  inline const double& get_E0() const { return valeur().get_E0(); };
  inline int& get_conserv_energy_global() { return valeur().get_conserv_energy_global(); };
  inline const char * get_fichier_sauvegarde() const { return valeur().get_fichier_sauvegarde(); };
  inline const int& get_ghost_fluid_flag() const { return valeur().get_ghost_fluid_flag(); };
  inline const int& get_probes_ghost_cells() { return valeur().get_ghost_cells();};
  inline void compute_ghost_cell_numbers_for_subproblems(const IJK_Splitting& splitting, int ghost_init) { return valeur().compute_ghost_cell_numbers_for_subproblems(splitting, ghost_init); };

  inline const int& get_debug() { return valeur().get_debug(); };
  inline const IJK_Field_double& get_debug_lrs_cells() { return valeur().get_debug_lrs_cells(); };

  inline int get_disable_post_processing_probes_out_files() const { return valeur().get_disable_post_processing_probes_out_files(); };
  /*
   * Setters
   */
  inline void set_field_T_ana() { return valeur().set_field_T_ana(); };
  void set_fichier_sauvegarde(const char *lata_name) { valeur().set_fichier_sauvegarde(lata_name); };

  inline int initialize(const IJK_Splitting& splitting, const int idx);
  inline void recompute_temperature_init();
  inline void update_thermal_properties();
  inline void euler_time_step(const double timestep);
  inline void rk3_sub_step(const int rk_step,
                           const double total_timestep,
                           const double time);
  inline void associer(const IJK_FT_double& ijk_ft);
  inline void associer_post(const IJK_FT_Post& ijk_ft_post);
  inline void associer_switch(const Switch_FT_double& ijk_ft_switch);
  inline void associer_interface_intersections(const Intersection_Interface_ijk_cell& intersection_ijk_cell,
                                               const Intersection_Interface_ijk_face& intersection_ijk_face);
  inline void sauvegarder_temperature(Nom& lata_name, int idx);
  inline double compute_timestep(const double timestep,
                                 const double dxmin);
  inline double compute_global_energy();
  inline void calculer_ecart_T_ana() { valeur().calculer_ecart_T_ana(); };
  inline void euler_rustine_step(const double timestep, const double dE);
  inline void rk3_rustine_sub_step(const int rk_step, const double total_timestep,
                                   const double fractionnal_timestep, const double time, const double dE);
  inline void compute_interfacial_temperature2(ArrOfDouble& interfacial_temperature,
                                               ArrOfDouble& flux_normal_interp);

  inline void update_intersections() { valeur().update_intersections(); };
  inline void clean_ijk_intersections() { valeur().clean_ijk_intersections(); };

  inline void compute_eulerian_curvature_from_interface();
  inline void compute_eulerian_distance();

  void posttraiter_tous_champs_thermal(Motcles& liste, const int idx) const;
  int posttraiter_champs_instantanes_thermal(const Motcles& liste_post_instantanes,
                                             const char * lata_name,
                                             const int latastep,
                                             const double current_time,
                                             const int idx);
  int posttraiter_champs_instantanes_thermal_interface(const Motcles& liste_post_instantanes,
                                                       const char *lata_name,
                                                       const int latastep,
                                                       const double current_time,
                                                       const int idx);
  int posttraiter_champs_instantanes_thermal_interface_ref(const Motcles& liste_post_instantanes,
                                                           const char *lata_name,
                                                           const int latastep,
                                                           const double current_time,
                                                           const int idx);
  Entree& typer_thermal( Entree& is );
  void thermal_subresolution_outputs(SFichier& fic);

protected:
  int thermal_rank_;
  Nom thermal_problem_type_;
  Nom prefix_;
  Motcles thermal_words_;
  Motcles lata_suffix_;

  REF(IJK_FT_double) ref_ijk_ft_;
  REF(IJK_FT_Post) ref_ijk_ft_post_;
  REF(Switch_FT_double) ref_ijk_ft_switch_;
  REF(Intersection_Interface_ijk_cell) ref_intersection_ijk_cell_;
  REF(Intersection_Interface_ijk_face) ref_intersection_ijk_face_;
};

inline int IJK_Thermal::initialize(const IJK_Splitting& splitting, const int idx)
{
  return valeur().initialize(splitting, idx);
}

inline void IJK_Thermal::recompute_temperature_init()
{
  valeur().recompute_temperature_init();
}

inline void IJK_Thermal::update_thermal_properties()
{
  return valeur().update_thermal_properties();
}

inline void IJK_Thermal::euler_time_step(const double timestep)
{
  valeur().euler_time_step(timestep);
}

inline void IJK_Thermal::associer(const IJK_FT_double& ijk_ft)
{
  ref_ijk_ft_ = ijk_ft;
  valeur().associer(ijk_ft);
}

inline void IJK_Thermal::associer_post(const IJK_FT_Post& ijk_ft_post)
{
  ref_ijk_ft_post_ = ijk_ft_post;
  valeur().associer_post(ijk_ft_post);
}

inline void IJK_Thermal::associer_switch(const Switch_FT_double& ijk_ft_switch)
{
  ref_ijk_ft_switch_ = ijk_ft_switch;
  valeur().associer_switch(ref_ijk_ft_switch_);
}

inline void IJK_Thermal::associer_interface_intersections(const Intersection_Interface_ijk_cell& intersection_ijk_cell,
                                                          const Intersection_Interface_ijk_face& intersection_ijk_face)
{
  ref_intersection_ijk_cell_ = intersection_ijk_cell;
  ref_intersection_ijk_face_ = intersection_ijk_face;
  valeur().associer_interface_intersections(intersection_ijk_cell, intersection_ijk_face);
}

inline void IJK_Thermal::sauvegarder_temperature(Nom& lata_name, int idx)
{
  valeur().sauvegarder_temperature(lata_name, idx);
}

inline double IJK_Thermal::compute_timestep(const double timestep, const double dxmin)
{
  return valeur().compute_timestep(timestep, dxmin);
}

inline double IJK_Thermal::compute_global_energy()
{
  return valeur().compute_global_energy();
}

inline void IJK_Thermal::euler_rustine_step(const double timestep, const double dE)
{
  valeur().euler_rustine_step(timestep, dE);
}

inline void IJK_Thermal::rk3_rustine_sub_step(const int rk_step, const double total_timestep,
                                              const double fractionnal_timestep, const double time, const double dE)
{
  valeur().rk3_rustine_sub_step(rk_step, total_timestep, fractionnal_timestep, time, dE);
}

inline void IJK_Thermal::rk3_sub_step(const int rk_step,
                                      const double total_timestep,
                                      const double time)
{
  valeur().rk3_sub_step(rk_step, total_timestep, time);
}

inline void IJK_Thermal::compute_interfacial_temperature2(ArrOfDouble& interfacial_temperature,
                                                          ArrOfDouble& flux_normal_interp)
{
  return valeur().compute_interfacial_temperature2(interfacial_temperature, flux_normal_interp);
}

inline void IJK_Thermal::compute_eulerian_distance()
{
  valeur().compute_eulerian_distance();
}

inline void IJK_Thermal::compute_eulerian_curvature_from_interface()
{
  valeur().compute_eulerian_curvature_from_interface();
}

#endif /* IJK_Thermal_included */
