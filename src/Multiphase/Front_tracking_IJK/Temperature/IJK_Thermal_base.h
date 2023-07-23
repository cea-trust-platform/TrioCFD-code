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
#include <IJK_Field.h>
#include <Parser.h>
#include <IJK_Lata_writer.h>
#include <Operateur_IJK_elem_conv.h>
#include <Operateur_IJK_elem_diff.h>
#include <Ouvrir_fichier.h>
#include <Corrige_flux_FT.h>
#include <TRUST_Ref.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_Thermal_base
//
// <Description of class IJK_Thermal_base>
//
/////////////////////////////////////////////////////////////////////////////

class IJK_FT_double;

class IJK_Thermal_base : public Objet_U
{

//  friend class IJK_FT_Post;
//  friend class IJK_FT_double;
  Declare_base( IJK_Thermal_base ) ;

public:

  /*
   * Initialisation
   */
  virtual void set_param(Param& param);
  virtual int initialize(const IJK_Splitting& splitting, const int idx);
  virtual void update_thermal_properties();
  double compute_timestep(const double timestep,
                          const double dxmin) const;

  void associer(const IJK_FT_double& ijk_ft);
  void euler_time_step(const double timestep);
  void rk3_sub_step(const int rk_step,
                    const double total_timestep,
                    const double time);
  void sauvegarder_temperature(Nom& lata_name, int idx);// const

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
  /*
   * Getters and setters
   */
  double get_rhocp_l() const;
  double get_rhocp_v() const;
  //  int get_type_thermal_problem() const;
  //  std::string set_type_thermal_problem();

  const IJK_Field_double& get_temperature() const
  {
    return temperature_ ;
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
  const IJK_Field_double& get_div_lambda_grad_T() const
  {
    return div_coeff_grad_T_volume_ ;
  }
  const IJK_Field_double& get_temperature_adim_bulles() const
  {
    return temperature_adim_bulles_;
  }
  const FixedVector<IJK_Field_double, 3>& get_gradient_temperature() const
  {
    return grad_T_ ;
  }
  const char * get_fichier_sauvegarde() const
  {
    return fichier_reprise_temperature_;
  }
  void set_fichier_sauvegarde(const char *lataname)
  {
    fichier_reprise_temperature_ = lataname;
  }
  void set_field_T_ana();

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

  void calculer_ecart_T_ana();
  void compute_interfacial_temperature2(
    ArrOfDouble& interfacial_temperature,
    ArrOfDouble& flux_normal_interp); //const ;
#if 0
  void ecrire_reprise_thermique(SFichier& fichier);
#endif


protected:

  void calculer_dT(const FixedVector<IJK_Field_double, 3>& velocity);
  void compute_temperature_convection(const FixedVector<IJK_Field_double, 3>& velocity);
  virtual void add_temperature_diffusion();
  virtual void compute_diffusion_increment()=0;
  virtual void correct_temperature_for_eulerian_fluxes()=0;
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
  //  void calculer_Nusselt(const IJK_Field_double& vx);
  //  void calculer_temperature_adimensionnelle_theta(const IJK_Field_double&  vx, const double qw);
  //  void calculer_temperature_physique_T_dummy();


  /*
   * Patch to conserve energy
   */
  double E0_; //volumique
  IJK_Field_double T_rust_;
  IJK_Field_double d_T_rustine_; // Temperature increment to conserve the energy.
  IJK_Field_double RK3_F_rustine_; // Temporary storage for substeps in the RK3 algorithm for the rustine calculation.

  REF(IJK_FT_double) ref_ijk_ft_;
  Corrige_flux_FT corrige_flux_;
  const IJK_Field_double& get_IJK_field(const Nom& nom) const;
  int rang_;

  /*
   * Physical parameters and inputs
   */
  double fo_;
  double cp_liquid_, cp_vapour_;
  double lambda_liquid_, lambda_vapour_;
  int single_phase_;
  double uniform_lambda_;
  double uniform_alpha_;
//  int type_thermal_problem_;
//  std::string name_type_thermal_problem_;
  /*
   * Initialisation (B.Cs, expression)
   */
  Boundary_Conditions_Thermique boundary_conditions_;
  IJK_Field_local_double boundary_flux_kmin_;
  IJK_Field_local_double boundary_flux_kmax_;
  Nom expression_T_init_;

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

  /*
   * Dimensionless temperature
   */
  IJK_Field_double temperature_physique_T_;
  IJK_Field_double temperature_adimensionnelle_theta_;
  IJK_Field_double temperature_adim_bulles_;

  /*
   * Storage for operators & time scheme
   */
  int diff_temp_negligible_;
  int conv_temperature_negligible_;

  /*
   * type_temperature_convection_op_:
   * 1 : Quick
   * 2 : Centre2
   */
  int type_temperature_convection_op_;
  Operateur_IJK_elem_conv temperature_convection_op_;
  Operateur_IJK_elem_diff temperature_diffusion_op_;
  IJK_Field_double div_coeff_grad_T_volume_;

  /*
   * Fields
   */
  IJK_Field_double rho_cp_;
  IJK_Field_double rho_cp_T_;
  IJK_Field_double temperature_;
  IJK_Field_double d_temperature_; // Temperature increment.
  IJK_Field_double RK3_F_temperature_; // Temporary storage for substeps in the RK3 algorithm.
  FixedVector<IJK_Field_double, 3> storage_; // Temporary storage for fluxes calculation.
  int calculate_local_energy_;
  int conserv_energy_global_;

  /*
   * Fields FT
   */
  IJK_Field_double temperature_ft_;

  /*
   * Post-processing
   */
  Nom expression_T_ana_;
  Motcles liste_post_instantanes_; // liste des champs instantanes a postraiter
  IJK_Field_double temperature_ana_, ecart_t_ana_;
  FixedVector<IJK_Field_double, 3> grad_T_;
  double global_energy_;
  int calulate_grad_T_;

};

#endif /* IJK_Thermal_base_included */
