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
// File      : IJK_Thermique.h
// Directory : $IJK_ROOT/src/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Thermique_included
#define IJK_Thermique_included

#include <IJK_Field.h>
#include <Objet_U.h>
#include <Boundary_Conditions_Thermique.h>
#include <IJK_Splitting.h>
#include <IJK_Field.h>
#include <Parser.h>
//#include <Interprete.h>
#include <IJK_Lata_writer.h>
#include <OpConvIJKQuickScalar.h>
#include <OpConvIJKAmont.h>
#include <OpCentre4IJK.h>
#include <OpDiffTurbIJKScalar.h>
#include <OpConvDiscIJKQuickScalar.h>
#include <OpConvCentre2IJKScalar.h>
#include <Ref_IJK_FT_double.h>
#include <Ouvrir_fichier.h>
#include <Corrige_flux_FT.h>


/*
#include <OpConvIJKQuickScalar.h>
#include <OpDiffTurbIJK.h>
#include <OpCentre4IJK.h>
#include <OpConvIJKQuickScalar.h>
#include <OpDiffTurbIJKScalar.h>
*/
/*! @brief : class IJK_Thermique
 *
 *  <Description of class IJK_Thermique>
 *
 *
 *
 */
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

class IJK_FT_double;

class IJK_Thermique;
Declare_ref(IJK_Thermique);

class IJK_Thermique : public Objet_U
{

  friend class IJK_FT_Post;
  friend class IJK_FT_double;
  Declare_instanciable( IJK_Thermique ) ;

public :

  int initialize(const IJK_Splitting& splitting, const int idx);
  void update_thermal_properties();
  double compute_timestep(const double timestep,
                          const double rho_l, const double rho_v,
                          const double dxmin) const;
  void associer(const IJK_FT_double& ijk_ft);
  void euler_time_step(const double timestep);
  void euler_rustine_step(const double timestep, const double dE);
  void rk3_sub_step(const int rk_step, const double total_timestep,
                    const double time  );
  void rk3_rustine_sub_step(const int rk_step, const double total_timestep,
                            const double fractionnal_timestep, const double time, const double dE);
  const IJK_Field_double& get_temperature() const
  {
    // if (liste_post_instantanes_.contient_("TEMPERATURE_ADIM_BULLES"))
    //   {
    //     return temperature_adim_bulles_;
    //   }
    // else
    //   {
    return temperature_ ;
    //  }
  }
  const IJK_Field_double& get_temperature_adim_bulles() const
  {
    return temperature_adim_bulles_;
  }

  IJK_Field_double& set_temperature()
  {
    return temperature_ ;
  }
  FixedVector<IJK_Field_double, 3>& get_gradient_temperature()
  {
    return grad_T_ ;
  }

  void sauvegarder_temperature(Nom& lata_name, int idx);// const
  void set_fichier_sauvegarde(const char *lataname)
  {
    fichier_reprise_temperature_ = lataname;
  };
  const char * get_fichier_sauvegarde() const
  {
    return fichier_reprise_temperature_;
  };
  /*  void set_reprise(const int i)
    {
      reprise_ = i;
      return ;
    }; */
#if 0
  void ecrire_reprise_thermique(SFichier& fichier);
#endif
  void set_field_T_ana();

  double compute_global_energy(const IJK_Field_double& temperature);
  double compute_global_energy()
  {
    return compute_global_energy(temperature_); // changes the attribute global_energy [J/m3]
  }
  double& get_global_energy()
  {
    return global_energy_ ;
  };
  FixedVector<IJK_Field_double, 3>& get_storage()
  {
    return storage_;
  };
  IJK_Field_double& get_temperature_ft()
  {
    return temperature_ft_;
  };
  double get_rhocp_l() const;
  double get_rhocp_v() const;

protected :
  void calculer_dT(const FixedVector<IJK_Field_double, 3>& velocity);
  void compute_dT_rustine(const double dE);
  void add_temperature_diffusion();
  void compute_temperature_convection(const FixedVector<IJK_Field_double, 3>& velocity);
  void compute_temperature_convection_conservative(const FixedVector<IJK_Field_double, 3>& velocity);
  void add_temperature_source();
  void source_callback();
  void calculer_temperature_physique_T(const IJK_Field_double&  vx, const double dTm);
  void maj_Tl_Tv();
  void calculer_temperature_adim_bulles();
  void calculer_energies(double& E_liq_pure, double& E_lta, double& E_lth,
                         double& E_vap_pure, double& E_vta, double& E_vth,
                         double& E_mixt_arithm, double& E_mixt_harmo, double& E_tot, double& E_tot_h);
  void calculer_temperature_physique_T_dummy();
  void calculer_temperature_adimensionnelle_theta(const IJK_Field_double&  vx, const double qw);
  void calculer_Nusselt(const IJK_Field_double& vx);
// void calculer_terme_source();

  void calculer_ecart_T_ana();
  void calculer_source_temperature_ana();
  void calculer_gradient_temperature(const IJK_Field_double& temperature, FixedVector<IJK_Field_double, 3>& grad_T);
  void compute_interfacial_temperature(ArrOfDouble& interfacial_temperature, ArrOfDouble& interfacial_phin_ai,
                                       FixedVector<IJK_Field_double, 3> storage_) const ;
  // This method calls to the Correction_flux_FT static method to build and interfacial temperature
  // and heat flux field at the interface.
  void compute_interfacial_temperature2(
    ArrOfDouble& interfacial_temperature,
    ArrOfDouble& flux_normal_interp) const ;

  REF(IJK_FT_double) ref_ijk_ft_;
  int rang_;

  Boundary_Conditions_Thermique boundary_conditions_;
  Nom expression_T_init_;
  Nom fichier_reprise_temperature_;
  int timestep_reprise_temperature_;
  Nom expression_T_ana_;
  Nom expression_source_temperature_;
  Nom type_T_source_;

  //MR: get field from IJK_FT_Post
  const IJK_Field_double& get_IJK_field(const Nom& nom) const;

//  Nom fichier_sauvegarde_temperature_;
// int timestep_sauvegarde_temperature_;

  int type_temperature_convection_op_; // 1 : Amont / 2 : Centre2 / 3 : Quick  / 4 : Centre
  OpConvAmontIJK_double temperature_convection_op_amont_;
  OpConvCentre2IJKScalar_double temperature_convection_op_centre2_;
  OpConvIJKQuickScalar_double temperature_convection_op_quick_;
  OpConvCentre4IJK_double temperature_convection_op_centre4_;

  OpDiffIJKScalar_double diffusion_temperature_op_;

  OpConvDiscIJKQuickScalar_double rho_cp_convection_op_quick_;

  // Storage for operators & time scheme:
  int diff_temp_negligible_;
  int conv_temperature_negligible_;



  //MR: lambda variable dans le terme source de Tryggvason
  int lambda_variable_;
  int wall_flux_;
  int depracated_rho_cp_;
  int rho_cp_inv_;
  int conserv_energy_global_;
  double global_energy_;
  IJK_Field_double div_lambda_grad_T_volume_;

  // Thermal heat flux on the boundary (intergal per boundary face)
  IJK_Field_local_double boundary_flux_kmin_;
  IJK_Field_local_double boundary_flux_kmax_;

  IJK_Field_double temperature_;
  IJK_Field_double d_temperature_; // Temperature increment.
  IJK_Field_double d_T_rustine_; // Temperature increment to conserve the energy.
  IJK_Field_double RK3_F_temperature_; // Temporary storage for substeps in the RK3 algorithm.
  IJK_Field_double RK3_F_rustine_; // Temporary storage for substeps in the RK3 algorithm for the rustine calculation.
  IJK_Field_double cp_;
  IJK_Field_double lambda_;
  IJK_Field_double rho_cp_;
  IJK_Field_double rho_cp_T_;
  IJK_Field_double div_rho_cp_T_;
  FixedVector<IJK_Field_double, 3> storage_; // // Temporary storage for fluxes calculations for instance.
  IJK_Field_double temperature_ft_;

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

  //Rustine
  double E0_;//volumique
  IJK_Field_double T_rust_;
  void compute_T_rust(const FixedVector<IJK_Field_double, 3>& velocity);

  //forme de l'equation de la tempetature
  int type_temperature_convection_form_; // 1 : non conservative, 2 : conservative

  IJK_Field_double temperature_physique_T_;
  IJK_Field_double temperature_adimensionnelle_theta_;
  IJK_Field_double temperature_adim_bulles_;

  double fo_;
  double cp_liquid_, cp_vapor_;
  double lambda_liquid_, lambda_vapor_;
  int lambda_moy_arith_;

//  double wall_flux_;
  Motcles liste_post_instantanes_; // liste des champs instantanes a postraiter

  IJK_Field_double temperature_ana_, ecart_t_ana_;

  FixedVector<IJK_Field_double, 3> grad_T_;

//  IJK_Field_double source_temperature_ana_, ecart_source_t_ana_;

// int reprise_; // Flag indiquant si on fait une reprise


};

#endif /* IJK_Thermique_included */
