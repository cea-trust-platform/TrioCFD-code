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
//  inline const char* get_thermal_problem_type() { return thermal_problem_type_.c_str(); };
  inline Nom& get_thermal_problem_type() { return thermal_problem_type_; };
  inline int& get_thermal_rank() { return thermal_rank_; };
  inline Motcles& get_thermal_words() { return thermal_words_; };
  inline Motcles& get_thermal_suffix() { return lata_suffix_; };
  inline const IJK_Field_double& get_temperature() const { return valeur().get_temperature(); };
  inline IJK_Field_double& get_temperature_ft() { return valeur().get_temperature_ft(); }
  inline const IJK_Field_double& get_temperature_ana() const { return valeur().get_temperature_ana(); };
  inline const IJK_Field_double& get_ecart_t_ana() const { return valeur().get_ecart_t_ana(); }
  inline const FixedVector<IJK_Field_double, 3>& get_grad_T() const { return valeur().get_grad_T(); }
  inline const IJK_Field_double& get_div_lambda_grad_T() const { return valeur().get_div_lambda_grad_T(); }
  inline const double& get_E0() const { return valeur().get_E0(); };
  inline int& get_conserv_energy_global() { return valeur().get_conserv_energy_global(); };
  inline const char * get_fichier_sauvegarde() const { return valeur().get_fichier_sauvegarde(); };
  /*
   * Setters
   */
  inline void set_field_T_ana() { return valeur().set_field_T_ana(); };
  void set_fichier_sauvegarde(const char *lata_name) { valeur().set_fichier_sauvegarde(lata_name); };

  inline int initialize(const IJK_Splitting& splitting, const int idx);
  inline void update_thermal_properties();
  inline void euler_time_step(const double timestep);
  inline void rk3_sub_step(const int rk_step,
                           const double total_timestep,
                           const double time);
  inline void associer(const IJK_FT_double& ijk_ft);
  inline void associer_post(const IJK_FT_Post& ijk_ft_post);
  inline void associer_switch(const Switch_FT_double& ijk_ft_switch);
  inline void sauvegarder_temperature(Nom& lata_name, int idx);
  inline double compute_timestep(const double timestep,
                                 const double dxmin) const;
  inline double compute_global_energy();
  inline void calculer_ecart_T_ana() { valeur().calculer_ecart_T_ana(); };
  inline void euler_rustine_step(const double timestep, const double dE);
  inline void rk3_rustine_sub_step(const int rk_step, const double total_timestep,
                                   const double fractionnal_timestep, const double time, const double dE);
  inline void compute_interfacial_temperature2(ArrOfDouble& interfacial_temperature,
                                               ArrOfDouble& flux_normal_interp);
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

protected:
  int thermal_rank_;
  Nom thermal_problem_type_;
  Nom prefix_;
  Motcles thermal_words_;
  Motcles lata_suffix_;

  REF(IJK_FT_double) ref_ijk_ft_;
  REF(IJK_FT_Post) ref_ijk_ft_post_;
  REF(Switch_FT_double) ref_ijk_ft_switch_;
};

inline int IJK_Thermal::initialize(const IJK_Splitting& splitting, const int idx)
{
  return valeur().initialize(splitting, idx);
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

inline void IJK_Thermal::sauvegarder_temperature(Nom& lata_name, int idx)
{
  valeur().sauvegarder_temperature(lata_name, idx);
}

inline double IJK_Thermal::compute_timestep(const double timestep, const double dxmin) const
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

#endif /* IJK_Thermal_included */
