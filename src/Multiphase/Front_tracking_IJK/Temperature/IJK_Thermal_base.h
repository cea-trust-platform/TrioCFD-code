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
// File      : IJK_Thermal_base.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Thermal_base_included
#define IJK_Thermal_base_included

#include <IJK_Field.h>
#include <Objet_U.h>
#include <Boundary_Conditions_Thermique.h>
#include <IJK_Splitting.h>
#include <Parser.h>
#include <IJK_Lata_writer.h>
#include <Operateur_IJK_elem_conv.h>
#include <Operateur_IJK_elem_diff.h>
#include <OpGradCentre2IJKScalar.h>
#include <OpHessCentre2IJKScalar.h>
#include <Ouvrir_fichier.h>
#include <Corrige_flux_FT.h>
#include <TRUST_Ref.h>
#include <IJK_FT_Post.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_Thermal_base
//
// <Description of class IJK_Thermal_base>
//
/////////////////////////////////////////////////////////////////////////////

class IJK_FT_double;
class Switch_FT_double;
class IJK_Interfaces;

class IJK_Thermal_base : public Objet_U
{
  Declare_base( IJK_Thermal_base ) ;
  friend class IJK_One_Dimensional_Subproblems;
  friend class IJK_One_Dimensional_Subproblem;
public:
  /*
   * Initialisation
   */
  virtual void set_param(Param& param);
  virtual int initialize(const IJK_Splitting& splitting, const int idx);
  virtual void update_thermal_properties();
  double compute_timestep(const double timestep,
                          const double dxmin);

  void associer(const IJK_FT_double& ijk_ft);
  void associer_post(const IJK_FT_Post& ijk_ft_post);
  void associer_switch(const Switch_FT_double& ijk_ft_switch);
  void associer_interface_intersections(const Intersection_Interface_ijk_cell& intersection_ijk_cell,
                                        const Intersection_Interface_ijk_face& intersection_ijk_face);
  void euler_time_step(const double timestep);
  void rk3_sub_step(const int rk_step,
                    const double total_timestep,
                    const double time);
  void sauvegarder_temperature(Nom& lata_name, int idx);  // const

  double compute_global_energy(const IJK_Field_double& temperature);
  double compute_global_energy()
  {
    return compute_global_energy(temperature_); // changes the attribute global_energy [J/m3]
  }
  int calculer_k_pour_bord(const IJK_Field_double& temperature, const bool bord_kmax);
  int calculer_flux_thermique_bord(const IJK_Field_double& temperature,
                                   const double lambda_de_t_paroi,
                                   const double T_paroi_impose,
                                   IJK_Field_local_double& flux_bord,
                                   const bool bord_kmax);
  int imposer_flux_thermique_bord(const IJK_Field_double& temperature,
                                  const double flux_paroi_impose,
                                  IJK_Field_local_double& flux_bord,
                                  const bool bord_kmax);
  virtual int get_first_step_thermals_post() { return 0; };
  /*
   * Getters and setters
   */
  double get_rhocp_l() const;
  double get_rhocp_v() const;

  const IJK_Field_double& get_temperature() const
  {
    return temperature_ ;
  }
  const IJK_Field_double& get_temperature_before_extrapolation() const
  {
    return temperature_before_extrapolation_ ;
  }
  IJK_Field_double& get_temperature_ft()
  {
    return temperature_ft_ ;
  }
  const FixedVector<IJK_Field_double, 3>& get_grad_T() const
  {
    return grad_T_ ;
  }
  IJK_Field_double& set_temperature()
  {
    return temperature_ ;
  }
  const IJK_Field_double& get_temperature_ana() const
  {
    return temperature_ana_ ;
  }
  const IJK_Field_double& get_ecart_t_ana() const
  {
    return ecart_t_ana_ ;
  }
  const IJK_Field_double& get_ecart_t_ana_rel() const
  {
    return ecart_t_ana_rel_;
  }
  const IJK_Field_double& get_div_lambda_grad_T() const
  {
    return div_coeff_grad_T_volume_ ;
  }
  const IJK_Field_double& get_u_T_convective() const
  {
    return u_T_convective_volume_;
  }
  const IJK_Field_double& get_u_T_convective_volume() const
  {
    return u_T_convective_volume_;
  }
  const IJK_Field_double& get_eulerian_distance_ft() const
  {
    return eulerian_distance_ft_;
  }
  const IJK_Field_double& get_eulerian_curvature_ft() const
  {
    return eulerian_curvature_ft_ ;
  }
  const IJK_Field_double& get_interfacial_area_ft() const
  {
    return eulerian_interfacial_area_ft_;
  }
  const IJK_Field_double& get_grad_T_interface_ft() const
  {
    return eulerian_grad_T_interface_ft_;
  }
  const IJK_Field_double& get_eulerian_compo_connex_ft() const
  {
    return *eulerian_compo_connex_ft_;
  }
  const IJK_Field_double& get_eulerian_compo_connex_ghost_ft() const
  {
    return *eulerian_compo_connex_ghost_ft_;
  }
  const IJK_Field_double& get_eulerian_compo_connex_from_interface_ft() const
  {
    return *eulerian_compo_connex_from_interface_ft_;
  }
  const IJK_Field_double& get_eulerian_compo_connex_from_interface_ghost_ft() const
  {
    return *eulerian_compo_connex_from_interface_ghost_ft_;
  }
  const IJK_Field_double& get_eulerian_compo_connex_ns() const
  {
    return *eulerian_compo_connex_ns_;
  }
  const IJK_Field_double& get_eulerian_compo_connex_ghost_ns() const
  {
    return *eulerian_compo_connex_ghost_ns_;
  }
  const IJK_Field_double& get_eulerian_compo_connex_from_interface_ns() const
  {
    return *eulerian_compo_connex_from_interface_ns_;
  }
  const IJK_Field_double& get_eulerian_compo_connex_from_interface_ghost_ns() const
  {
    return *eulerian_compo_connex_from_interface_ghost_ns_;
  }
  const IJK_Field_int& get_eulerian_compo_connex_int_from_interface_ns() const
  {
    return *eulerian_compo_connex_from_interface_int_ns_;
  }
  const IJK_Field_int& get_eulerian_compo_connex_int_from_interface_ghost_ns() const
  {
    return *eulerian_compo_connex_from_interface_ghost_int_ns_;
  }
  const IJK_Field_double& get_eulerian_distance_ns() const
  {
    return eulerian_distance_ns_;
  }
  const IJK_Field_double& get_eulerian_curvature_ns() const
  {
    return eulerian_curvature_ns_ ;
  }
  const IJK_Field_double& get_interfacial_area_ns() const
  {
    return eulerian_interfacial_area_ns_;
  }
  const IJK_Field_double& get_grad_T_interface_ns() const
  {
    return eulerian_grad_T_interface_ns_;
  }
  const IJK_Field_double& get_eulerian_rising_velocities() const
  {
    return eulerian_rising_velocities_;
  }
  const IJK_Field_double& get_temperature_adim_bulles() const
  {
    return temperature_adim_bulles_;
  }
  const FixedVector<IJK_Field_double, 3>& get_gradient_temperature() const
  {
    return grad_T_ ;
  }
  const FixedVector<IJK_Field_double, 3>& get_gradient_temperature_elem() const
  {
    return grad_T_elem_ ;
  }
  const FixedVector<IJK_Field_double, 3>& get_normal_vector_ns() const
  {
    return eulerian_normal_vectors_ns_;
  }
  const FixedVector<IJK_Field_double, 3>& get_normal_vector_ns_normed() const
  {
    return eulerian_normal_vectors_ns_normed_;
  }
  const FixedVector<IJK_Field_double, 3>& get_normal_vector_ft() const
  {
    return eulerian_normal_vectors_ft_;
  }
  const FixedVector<IJK_Field_double, 3>& get_hessian_diag_temperature_elem() const
  {
    return hess_diag_T_elem_ ;
  }
  const FixedVector<IJK_Field_double, 3>& get_hessian_cross_temperature_elem() const
  {
    return hess_cross_T_elem_ ;
  }
  const FixedVector<IJK_Field_double, 3>& get_bary() const
  {
    return eulerian_facets_barycentre_ft_;
  }
  const int& get_ghost_fluid_flag() const
  {
    return ghost_fluid_;
  };
  const int& get_ghost_cells() const
  {
    return ghost_cells_;
  };
  const int& get_debug() const
  {
    return debug_;
  };
  virtual const IJK_Field_double& get_temperature_cell_neighbours_debug() const
  {
    return dummy_double_field_; //dummy
  }
  virtual const IJK_Field_double& get_temperature_cell_neighbours() const
  {
    return dummy_double_field_; //dummy
  }
  virtual const IJK_Field_int& get_cell_neighbours_corrected() const
  {
    return dummy_int_field_; //dummy
  }
  virtual const IJK_Field_double& get_neighbours_temperature_colinearity_weighting() const
  {
    return dummy_double_field_; //dummy
  }
  virtual const IJK_Field_double& get_debug_lrs_cells() const
  {
    return dummy_double_field_; //dummy
  };
  virtual int get_disable_post_processing_probes_out_files() const
  {
    return 1;
  };
  virtual const FixedVector<IJK_Field_double,3>& get_cell_faces_corrected_diffusive() const
  {
    return dummy_double_vect_; //dummy
  }
  virtual const FixedVector<IJK_Field_double,3>& get_cell_faces_corrected_convective() const
  {
    return dummy_double_vect_; //dummy
  }
  virtual const FixedVector<IJK_Field_int,3>& get_cell_faces_corrected_bool() const
  {
    return dummy_int_vect_;
  }
  virtual const FixedVector<IJK_Field_int,3>& get_cell_faces_neighbours_corrected_diag_bool() const
  {
    return dummy_int_vect_;
  }
  virtual const FixedVector<IJK_Field_int,3>& get_cell_faces_neighbours_corrected_all_bool() const
  {
    return dummy_int_vect_;
  }
  virtual const FixedVector<IJK_Field_int,3>& get_cell_faces_neighbours_corrected_min_max_bool() const
  {
    return dummy_int_vect_;
  }
  virtual const FixedVector<IJK_Field_double,3>& get_cell_faces_neighbours_corrected_convective() const
  {
    return dummy_double_vect_;
  }
  virtual const FixedVector<IJK_Field_double,3>& get_cell_faces_neighbours_corrected_diffusive() const
  {
    return dummy_double_vect_;
  }
  virtual const FixedVector<IJK_Field_double,3>& get_neighbours_faces_weighting_colinearity() const
  {
    return dummy_double_vect_;
  }
  virtual const IJK_Field_int& get_cell_neighbours_corrected_trimmed() const
  {
    return dummy_int_field_;
  }
  virtual double get_modified_time();

  virtual double get_rho_cp_u_ijk(const IJK_Field_double& vx, int i, int j, int k) const;
  virtual double get_div_lambda_ijk(int i, int j, int k) const { return 0; };
  virtual double compute_temperature_dimensionless_theta_mean(const IJK_Field_double& vx);

  const char * get_fichier_sauvegarde() const
  {
    return fichier_reprise_temperature_;
  }
  void set_fichier_sauvegarde(const char *lataname)
  {
    fichier_reprise_temperature_ = lataname;
  }
  virtual void set_field_T_ana();

  /*
   * Patch from IJK_Thermique
   */
  void euler_rustine_step(const double timestep, const double dE);
  void rk3_rustine_sub_step(const int rk_step, const double total_timestep,
                            const double fractionnal_timestep, const double time, const double dE);
  virtual int& get_conserv_energy_global() { return conserv_energy_global_; };
  const double& get_E0() const { return E0_; };

  void compute_dT_rustine(const double dE);
  void compute_T_rust(const FixedVector<IJK_Field_double, 3>& velocity);

  virtual void calculer_ecart_T_ana();
  void compute_interfacial_temperature2(
    ArrOfDouble& interfacial_temperature,
    ArrOfDouble& flux_normal_interp); //const ;
#if 0
  void ecrire_reprise_thermique(SFichier& fichier);
#endif
  virtual void compute_ghost_cell_numbers_for_subproblems(const IJK_Splitting& splitting, int ghost_init) { ghost_cells_ = ghost_init; };
  void compute_eulerian_distance();
  void compute_eulerian_curvature_from_interface();
  virtual void update_intersections() { ; };
  virtual void clean_ijk_intersections() { ; };
  virtual void set_thermal_subresolution_outputs(const Nom& interfacial_quantities_thermal_probes,
                                                 const Nom& overall_bubbles_quantities,
                                                 const Nom& local_quantities_thermal_probes_time_index_folder) { ; };
  virtual void compute_temperature_init();
  virtual void recompute_temperature_init();

protected:

  void compute_cell_volume();
  void compute_min_cell_delta();
  void compute_cell_diagonal(const IJK_Splitting& splitting);
  void calculer_dT(const FixedVector<IJK_Field_double, 3>& velocity);
  virtual void post_process_after_temperature_increment();
  void compute_temperature_convection(const FixedVector<IJK_Field_double, 3>& velocity);
  virtual void add_temperature_diffusion();
  virtual void compute_diffusion_increment()=0;
  virtual void correct_temperature_for_eulerian_fluxes()=0;
  virtual void store_temperature_before_extrapolation() { ; };
  virtual void correct_temperature_increment_for_interface_leaving_cell() { ; };
  void enforce_zero_value_eulerian_distance();
  void compute_eulerian_curvature();
  void enforce_zero_value_eulerian_curvature();
  void enforce_max_value_eulerian_curvature();
  void compute_eulerian_grad_T_interface();
  void propagate_eulerian_grad_T_interface();
  void compute_eulerian_temperature_ghost();
  void compute_eulerian_bounding_box_fill_compo();
  void compute_rising_velocities();
  void enforce_zero_value_eulerian_field(IJK_Field_double& eulerian_field);
  void enforce_max_value_eulerian_field(IJK_Field_double& eulerian_field);
  void enforce_min_value_eulerian_field(IJK_Field_double& eulerian_field);
  void compute_temperature_gradient_elem();
  void compute_temperature_hessian_diag_elem();
  void compute_temperature_hessian_cross_elem();
  virtual void correct_temperature_for_visu() { ; };
  virtual void correct_operators_for_visu() { ; };
  virtual void clip_temperature_values() { ; };
  virtual void compute_thermal_subproblems() { ; };
  virtual void compute_convective_diffusive_fluxes_face_centre() { ; };
  virtual void compute_convective_fluxes_face_centre() { ; };
  virtual void compute_diffusive_fluxes_face_centre() { ; };
  virtual void prepare_ij_fluxes_k_layers() { ; };
  virtual void compute_temperature_cell_centres(const int first_corr) { ; };
  virtual void set_zero_temperature_increment() { ; };
  virtual void clean_thermal_subproblems() { ; };

  void calculer_gradient_temperature(const IJK_Field_double& temperature,
                                     FixedVector<IJK_Field_double, 3>& grad_T);
  void calculer_energies(double& E_liq_pure,
                         double& E_liq,
                         double& E_vap_pure,
                         double& E_vap,
                         double& E_mixt,
                         double& E_tot);

  void source_callback();
  void calculer_temperature_physique_T(const IJK_Field_double&  vx, const double dTm);
  void calculer_temperature_adim_bulles();
  void add_temperature_source();
  void calculer_Nusselt(const IJK_Field_double& vx);
  void calculer_temperature_adimensionnelle_theta(const IJK_Field_double&  vx, const double qw);
  void calculer_source_temperature_ana();
  virtual double compute_rho_cp_u_mean(const IJK_Field_double& vx);
  double compute_variable_wall_temperature(const int kmin, const int kmax);

  void force_upstream_temperature(IJK_Field_double& temperature, double T_imposed,
                                  const IJK_Interfaces& interfaces, double nb_diam, int upstream_dir,
                                  int gravity_dir, int upstream_stencil);
  virtual void enforce_periodic_temperature_boundary_value() { ; } ;

  int debug_;
  /*
   * Patch to conserve energy
   */
  double E0_; //volumique
  IJK_Field_double T_rust_;
  IJK_Field_double d_T_rustine_; // Temperature increment to conserve the energy.
  IJK_Field_double RK3_F_rustine_; // Temporary storage for substeps in the RK3 algorithm for the rustine calculation.

  REF(IJK_FT_double) ref_ijk_ft_;
  REF(IJK_FT_Post) ref_ijk_ft_post_;
  REF(Switch_FT_double) ref_ijk_ft_switch_;
  REF(Intersection_Interface_ijk_cell) ref_intersection_ijk_cell_;
  REF(Intersection_Interface_ijk_face) ref_intersection_ijk_face_;
  Corrige_flux_FT corrige_flux_;
  const IJK_Field_double& get_IJK_field(const Nom& nom) const;
  int rang_;

  /*
   * Physical parameters and inputs
   */
  double dt_fo_;
  double fo_;
  double cp_liquid_, cp_vapour_;
  double lambda_liquid_, lambda_vapour_;
  int single_phase_;
  double uniform_lambda_;
  double uniform_alpha_;
  /*
   * Initialisation (B.Cs, expression)
   */
  Boundary_Conditions_Thermique boundary_conditions_;
  IJK_Field_local_double boundary_flux_kmin_;
  IJK_Field_local_double boundary_flux_kmax_;
  Nom expression_T_init_;
  double upstream_temperature_;
  double nb_diam_upstream_;
  int side_temperature_;
  int stencil_side_;

  /*
   * Settings to resume a calculation
   */
  Nom fichier_reprise_temperature_;
  int timestep_reprise_temperature_;

  /*
   * Source of temperature, wall heating,
   * integral per boundary (Tryggvason)
   */
  Nom expression_source_temperature_;
  Nom type_T_source_;
  int lambda_variable_;
  int wall_flux_;
  IJK_Field_double source_temperature_;
  IJK_Field_double source_temperature_v_;
  IJK_Field_double source_temperature_l_;
  IJK_Field_double d_source_Tl_;
  IJK_Field_double d_source_Tv_;
  double dTl_;
  double dTv_;
  double Tl_;
  double Tv_;
  double Tv0_; // Serait-ce plutot Tref (une temperature de reference pour reconstruire le champ dim??)
  double kv_;
  double kl_;
  double T0v_;
  double T0l_;
  IJK_Field_double source_temperature_ana_;
  IJK_Field_double ecart_source_t_ana_;

  /*
   * Dimensionless temperature
   */
  IJK_Field_double temperature_physique_T_;
  IJK_Field_double temperature_adimensionnelle_theta_;
  IJK_Field_double temperature_adim_bulles_;

  /*
   * Storage for operators & time scheme
   */
  int diff_temperature_negligible_;
  int conv_temperature_negligible_;

  /*
   * type_temperature_convection_op_:
   * 1 : Quick
   * 2 : Centre2
   */
  Operateur_IJK_elem_conv temperature_convection_op_;
  Operateur_IJK_elem_diff temperature_diffusion_op_;
  IJK_Field_double div_coeff_grad_T_volume_;
  IJK_Field_double u_T_convective_volume_;
  OpGradCentre2IJKScalar_double temperature_grad_op_centre_;
  OpHessCentre2IJKScalar_double temperature_hess_op_centre_;


  /*
   * Fields
   */
  double vol_;
  double min_delta_xyz_;
  double cell_diagonal_;
  int ghost_cells_;
  IJK_Field_double rho_cp_;
  IJK_Field_double rho_cp_T_;
  IJK_Field_double temperature_;
  IJK_Field_double temperature_before_extrapolation_;
  IJK_Field_double d_temperature_; // Temperature increment.
  IJK_Field_double RK3_F_temperature_; // Temporary storage for substeps in the RK3 algorithm.
  FixedVector<IJK_Field_double, 3> storage_; // Temporary storage for fluxes calculation.
  int calculate_local_energy_;
  int conserv_energy_global_;

  /*
   * Fields FT
   * TODO: Clean FT_fields and redistribute curvature, interfacial_area if necessary
   */
  IJK_Field_double temperature_ft_;

  /*
   * Post-processing
   */
  Nom expression_T_ana_;
  Motcles liste_post_instantanes_; // liste des champs instantanes a postraiter
  IJK_Field_double temperature_ana_, ecart_t_ana_, ecart_t_ana_rel_;
  FixedVector<IJK_Field_double, 3> grad_T_;
  double global_energy_;
  int calulate_grad_T_;
  int rho_cp_post_;

  /*
   * For Ghost fluid method & Subresolution or Post-processing
   */
  int ghost_fluid_;
  int n_iter_distance_;
  int compute_distance_;
  int compute_curvature_;
  int compute_grad_T_interface_;
  /*
   * TODO: Move fields and avoid redundancies in IJK_Interfaces
   * Clean FT_fields
   */
  DoubleTab bounding_box_;
  DoubleTab min_max_larger_box_;
  IJK_Field_double eulerian_distance_ft_;
  IJK_Field_double eulerian_distance_ns_;
  FixedVector<IJK_Field_double, 3> eulerian_normal_vectors_ft_;
  FixedVector<IJK_Field_double, 3> eulerian_facets_barycentre_ft_;
  FixedVector<IJK_Field_double, 3> eulerian_normal_vectors_ns_;
  FixedVector<IJK_Field_double, 3> eulerian_normal_vectors_ns_normed_;
  FixedVector<IJK_Field_double, 3> eulerian_facets_barycentre_ns_;
  IJK_Field_double eulerian_curvature_ft_;
  IJK_Field_double eulerian_curvature_ns_;
  IJK_Field_double eulerian_interfacial_area_ft_;
  IJK_Field_double eulerian_interfacial_area_ns_;
  IJK_Field_double eulerian_grad_T_interface_ft_;
  IJK_Field_double eulerian_grad_T_interface_ns_;
  int compute_grad_T_elem_;
  FixedVector<IJK_Field_double, 3> grad_T_elem_;
  /*
   * hess(T) = grad(grad(T))
   * Only 6 coefficients in
   * cartesian coordinate system
   * | * * * |
   * | - * * |
   * | - - * |
   *   FixedVector<IJK_Field_double, 6> grad_grad_T_elem_;
   */
  FixedVector<IJK_Field_double, 3> hess_diag_T_elem_;
  FixedVector<IJK_Field_double, 3> hess_cross_T_elem_;
  FixedVector<IJK_Field_double, 3> facets_barycentre;
  IJK_Field_double d_temperature_uncorrected_;
  IJK_Field_double div_coeff_grad_T_volume_uncorrected_;
  Operateur_IJK_elem_conv temperature_convection_op_uncorrected_;
  Operateur_IJK_elem_diff temperature_diffusion_op_uncorrected_;
  int compute_hess_T_elem_;
  int compute_hess_diag_T_elem_;
  int compute_hess_cross_T_elem_;

  int mixed_cells_number_ = 0;
  void compute_mixed_cells_number(const IJK_Field_double& indicator);
  int compute_eulerian_compo_;

//  IJK_Field_double eulerian_compo_connex_ft_;
//  IJK_Field_double eulerian_compo_connex_ns_;
//  IJK_Field_double eulerian_compo_connex_ghost_ft_;
//  IJK_Field_double eulerian_compo_connex_ghost_ns_;
  int spherical_approx_;
  const IJK_Field_double * eulerian_compo_connex_ft_;
  const IJK_Field_double * eulerian_compo_connex_ns_;
  const IJK_Field_double * eulerian_compo_connex_ghost_ft_;
  const IJK_Field_double * eulerian_compo_connex_ghost_ns_;
  const IJK_Field_double * eulerian_compo_connex_from_interface_ft_;
  const IJK_Field_double * eulerian_compo_connex_from_interface_ns_;
  const IJK_Field_double * eulerian_compo_connex_from_interface_ghost_ft_;
  const IJK_Field_double * eulerian_compo_connex_from_interface_ghost_ns_;
  const IJK_Field_int * eulerian_compo_connex_from_interface_int_ns_;
  const IJK_Field_int * eulerian_compo_connex_from_interface_ghost_int_ns_;

  int compute_rising_velocities_;
  int fill_rising_velocities_;
  ArrOfDouble rising_velocities_;
  DoubleTab rising_vectors_;
  IJK_Field_double eulerian_rising_velocities_;
  ArrOfDouble bubbles_volume_;
  DoubleTab bubbles_barycentre_;

  FixedVector<IJK_Field_int,3> dummy_int_vect_;
  FixedVector<IJK_Field_double,3> dummy_double_vect_;
  IJK_Field_int dummy_int_field_;
  IJK_Field_double dummy_double_field_;
};

#endif /* IJK_Thermal_base_included */
