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

#ifndef OpCentre4IJK_H
#define OpCentre4IJK_H

#include <OpConvIJKFacesCommon.h>
#include <Boundary_Conditions.h>
#include <Boundary_Conditions_Thermique.h>

class OpConvCentre4IJK_double : public OpConvIJKFacesCommon_double
{
public:
  void initialize(const IJK_Splitting& splitting, const Boundary_Conditions& bc);
  void initialize(const IJK_Splitting& splitting, const Boundary_Conditions_Thermique& bc);
  void calculer(const IJK_Field_double& inputx, const IJK_Field_double& inputy, const IJK_Field_double& inputz,
                const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  void ajouter(const IJK_Field_double& inputx, const IJK_Field_double& inputy, const IJK_Field_double& inputz,
               const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
               IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);

  void calculer_avec_u_div_rhou(const IJK_Field_double& rhovx, const IJK_Field_double& rhovy, const IJK_Field_double& rhovz,
                                const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz,
                                IJK_Field_double& div_rho_u);
  void ajouter_avec_u_div_rhou(const IJK_Field_double& rhovx, const IJK_Field_double& rhovy, const IJK_Field_double& rhovz,
                               const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                               IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz,
                               IJK_Field_double& div_rho_u);
protected:
  inline void compute_flux_x_vx(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::X,DIRECTION::X>(resu,k_layer);
  }
  inline void compute_flux_x_vy(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::X,DIRECTION::Y>(resu,k_layer);
  }
  inline void compute_flux_x_vz(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::X,DIRECTION::Z>(resu,k_layer);
  }
  inline void compute_flux_y_vx(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Y,DIRECTION::X>(resu,k_layer);
  }
  inline void compute_flux_y_vy(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Y,DIRECTION::Y>(resu,k_layer);
  }
  inline void compute_flux_y_vz(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Y,DIRECTION::Z>(resu,k_layer);
  }
  inline void compute_flux_z_vx(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Z,DIRECTION::X>(resu,k_layer);
  }
  inline void compute_flux_z_vy(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Z,DIRECTION::Y>(resu,k_layer);
  }
  inline void compute_flux_z_vz(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Z,DIRECTION::Z>(resu,k_layer);
  }
  inline void exec_after_divergence_flux_x(IJK_Field_double& resu, const int k_layer) override
  {
    exec_after_divergence_flux_<DIRECTION::X>(resu, k_layer);
  }
  inline void exec_after_divergence_flux_y(IJK_Field_double& resu, const int k_layer) override
  {
    exec_after_divergence_flux_<DIRECTION::Y>(resu, k_layer);
  }
  inline void exec_after_divergence_flux_z(IJK_Field_double& resu, const int k_layer) override
  {
    exec_after_divergence_flux_<DIRECTION::Z>(resu, k_layer);
  }
  // First layer of non zero fluxes (eg, not in the walls)
  // (specific for wall boundary conditions at zmin and zmax)

  double get_g(int k_layer, int compo, int dir, int g_index) const
  {
    if (dir != 2)
      {
        Process::exit(); // error: non uniform mesh in i or j not coded
      }
    if (compo != 2)
      return g_compo_xy_dir_z_(k_layer, g_index);
    else
      return g_compo_z_dir_z_(k_layer, g_index);
  }

  // order 4 filtering coefficients for direction z fluxes
  //  (when non uniform mesh in z).
  DoubleTab g_compo_xy_dir_z_;
  DoubleTab g_compo_z_dir_z_;

  // Pointer to div_rho_u, filled by operator if ajouter_avec_u_div_rhou() is called
  IJK_Field_double *div_rho_u_;
  int last_computed_klayer_for_div_rhou_;
  REF(Boundary_Conditions) ref_bc_;
  REF(Boundary_Conditions_Thermique) ref_bc_Thermique_;
private:
  void calculer_div_rhou(const IJK_Field_double& rhovx, const IJK_Field_double& rhovy, const IJK_Field_double& rhovz,
                         IJK_Field_double& resu, int k_layer, const Operateur_IJK_data_channel& channel);
  template <DIRECTION _DIR_, DIRECTION _VCOMPO_>
  void compute_flux_(IJK_Field_local_double& resu, const int k_layer);
  template <DIRECTION _DIR_>
  void exec_after_divergence_flux_(IJK_Field_double& resu, const int k_layer);

};

#include <OpCentre4IJK.tpp>
#endif
