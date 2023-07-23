/****************************************************************************
 * Copyright (c) 2019, CEA
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 *modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *this list of conditions and the following disclaimer in the documentation
 *and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *may be used to endorse or promote products derived from this software without
 *specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *****************************************************************************/
/////////////////////////////////////////////////////////////////////////////
//
// File      : IJK_Energie.h
// Directory : $IJK_ROOT/src/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Energie_included
#define IJK_Energie_included

#include <Boundary_Conditions_Thermique.h>
#include <IJK_FT_Post.h>
#include <IJK_Field.h>
#include <IJK_Lata_writer.h>
#include <IJK_Splitting.h>
#include <MonofluidVar.h>
#include <Objet_U.h>
#include <OpConvQuickInterfaceIJKScalar.h>
#include <Ouvrir_fichier.h>
#include <Parser.h>
#include <TRUST_Ref.h>
#include <Corrige_flux_FT.h>
#include <Corrige_flux_FT_temperature_conv.h>
#include <Operateur_IJK_elem_diff_base.h>
// #include <Corrige_flux_FT_temperature_conv.h>

class IJK_FT_double;

/*! @brief : class IJK_Energie
 *
 *  <Description of class IJK_Energie>
 *
 *
 *
 */
class IJK_FT_double;


class IJK_Energie : public Objet_U
{
  friend class IJK_FT_Post;
  friend class IJK_FT_double;
  Declare_instanciable(IJK_Energie);

public:
  int initialize(const IJK_Splitting& splitting, const int idx);
  void update_thermal_properties();
  double compute_timestep(const double timestep, const double dxmin) const;
  void associer(const IJK_FT_double& ijk_ft);
  void euler_time_step(const FixedVector<IJK_Field_double, 3>& velocity);
  const IJK_Field_double& get_temperature() const { return temperature_; }
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
  IJK_Field_double& set_temperature() { return temperature_; }
  FixedVector<IJK_Field_double, 3>& get_gradient_temperature()
  {
    return grad_T_;
  }

  void sauvegarder_temperature(Nom& lata_name, int idx); // const
  void set_fichier_sauvegarde(const char *lataname)
  {
    fichier_reprise_temperature_ = lataname;
  };
  const char *get_fichier_sauvegarde() const
  {
    return fichier_reprise_temperature_;
  };
  void set_field_T_ana();

  double compute_global_energy(const IJK_Field_double& temperature);
  double compute_global_energy()
  {
    return compute_global_energy(
             temperature_); // changes the attribute global_energy [J/m3]
  }
  double& get_global_energy() { return global_energy_; };
  IJK_Field_double& get_temperature_ft() { return temperature_ft_; };
  double get_rhocp_l() const;
  double get_rhocp_v() const;
  double get_lda_l() const;
  double get_lda_v() const;

protected:
  void calculer_dT(const FixedVector<IJK_Field_double, 3>& velocity);
  void add_temperature_diffusion();
  void add_temporal_rho_cp_term();
  void divide_by_rho_cp_np1();
  void compute_energy_convection(
    const FixedVector<IJK_Field_double, 3>& velocity);
  void calculer_energies(double& E_liq_pure, double& E_lta, double& E_lth,
                         double& E_vap_pure, double& E_vta, double& E_vth,
                         double& E_mixt_arithm, double& E_mixt_harmo,
                         double& E_tot, double& E_tot_h);
  void calculer_Nusselt(const IJK_Field_double& vx);

  void calculer_ecart_T_ana();
  void calculer_gradient_temperature(const IJK_Field_double& temperature,
                                     FixedVector<IJK_Field_double, 3>& grad_T);
  void compute_interfacial_temperature2(ArrOfDouble& interfacial_temperature,
                                        ArrOfDouble& interfacial_phin_ai) const;

  REF(IJK_FT_double) ref_ijk_ft_;
  int rang_;

  Boundary_Conditions_Thermique boundary_conditions_;
  Nom expression_T_init_;
  Nom fichier_reprise_temperature_;
  int timestep_reprise_temperature_;
  Nom expression_T_ana_;

  // MR: get field from IJK_FT_Post
  const IJK_Field_double& get_IJK_field(const Nom& nom) const;

  //  Nom fichier_sauvegarde_temperature_;
  // int timestep_sauvegarde_temperature_;

  Corrige_flux_FT corrige_flux_;
//  Corrige_flux_FT_temperature_conv corrige_flux_temp_conv_;

  OpConvQuickInterfaceIJKScalar_double energy_convection_op_quick_interface_;
  OpDiffIJKScalar_double diffusion_temperature_op_;

  int diff_temp_negligible_;
  int conv_temperature_negligible_;

  double global_energy_;
  IJK_Field_double div_lambda_grad_T_volume_;

  // Thermal heat flux on the boundary (intergal per boundary face)
  IJK_Field_local_double boundary_flux_kmin_;
  IJK_Field_local_double boundary_flux_kmax_;

  double cp_liquid_, cp_vapor_;
  double lambda_liquid_, lambda_vapor_;

  IJK_MonofluidVar cp_;
  IJK_MonofluidVar lda_;
  IJK_MonofluidVar rho_cp_;
  IJK_MonofluidVar rho_;

  IJK_Field_double temperature_;
  IJK_Field_double d_temperature_; // Temperature increment.

  IJK_Field_double lambda_;
  IJK_Field_double rho_cp_T_;
  IJK_Field_double D_rhocp_T_;
  IJK_Field_double div_rho_cp_T_;
  IJK_Field_double temperature_ft_;

  double fo_;

  Motcles liste_post_instantanes_; // liste des champs instantanes a postraiter

  IJK_Field_double temperature_ana_, ecart_t_ana_;

  FixedVector<IJK_Field_double, 3> grad_T_;
};

#endif /* IJK_Energie_included */
