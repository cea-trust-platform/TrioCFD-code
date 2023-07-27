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
// File      : Deriv_Operateur_IJK_faces_diff.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Operateur_IJK_faces_diff_included
#define Operateur_IJK_faces_diff_included

#include <TRUST_Deriv.h>
#include <Operateur_IJK_faces_diff_base.h>

class Operateur_IJK_faces_diff : public DERIV( Operateur_IJK_faces_diff_base_double )
{
  Declare_instanciable( Operateur_IJK_faces_diff ) ;
public:
  inline double compute_dtstab_convection_local(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  inline void compute_set(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  inline void compute_add(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  inline void initialize(const IJK_Splitting& splitting);
  inline void compute_flux_x_vx(IJK_Field_local_double& resu, const int k_layer);
  inline void compute_flux_x_vy(IJK_Field_local_double& resu, const int k_layer);
  inline void compute_flux_x_vz(IJK_Field_local_double& resu, const int k_layer);
  inline void compute_flux_y_vx(IJK_Field_local_double& resu, const int k_layer);
  inline void compute_flux_y_vy(IJK_Field_local_double& resu, const int k_layer);
  inline void compute_flux_y_vz(IJK_Field_local_double& resu, const int k_layer);
  inline void compute_flux_z_vx(IJK_Field_local_double& resu, const int k_layer);
  inline void compute_flux_z_vy(IJK_Field_local_double& resu, const int k_layer);
  inline void compute_flux_z_vz(IJK_Field_local_double& resu, const int k_layer);
  inline void set_bc(const Boundary_Conditions& bc);
  inline void ajouter(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                      IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  inline void calculer(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                       IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  inline void set_coeff_x_y_z(const IJK_Field_double& coeff_tensor_xx,
                              const IJK_Field_double& coeff_tensor_xy,
                              const IJK_Field_double& coeff_tensor_xz,
                              const IJK_Field_double& coeff_tensor_yx,
                              const IJK_Field_double& coeff_tensor_yy,
                              const IJK_Field_double& coeff_tensor_yz,
                              const IJK_Field_double& coeff_tensor_zx,
                              const IJK_Field_double& coeff_tensor_zy,
                              const IJK_Field_double& coeff_tensor_zz);

  inline void set_nu(const IJK_Field_local_double& nu);
  inline void set_uniform_nu(const double& nu);
  inline void set_divergence(const IJK_Field_local_double& divergence);
  inline const IJK_Field_local_double& get_v(DIRECTION _DIR_);
  inline const IJK_Field_local_double& get_nu();
  inline const double& get_uniform_nu();
  inline const IJK_Field_local_double& get_divergence();
  inline const IJK_Field_local_double& get_coeff_tensor(DIRECTION _COMPO1_, DIRECTION _COMPO2_);

  /*
   * ReadOn
   */
  Entree& typer_diffusion_op( Entree& is );
  void typer_diffusion_op( const char * diffusion_op );
  int lire_motcle_non_standard(const Motcle& mot, Entree& is) override;
  void set_param(Param& param);
  Nom get_diffusion_op_type( Motcle word );
  /*
   * Getters
   */
  Nom get_diffusion_op_option() { return diffusion_option_; };
  Nom get_diffusion_op() { return diffusion_op_; };
  /*
   * Setters
   */

protected:
  Motcles diffusion_op_words_;
  Motcles diffusion_op_options_;
  Nom prefix_;
  Nom suffix_;
  int diffusion_rank_;
  Nom diffusion_op_;
  Nom diffusion_option_;
  bool is_cast_;
};

inline double Operateur_IJK_faces_diff::compute_dtstab_convection_local(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  return valeur().compute_dtstab_convection_local(dvx, dvy, dvz);
}

inline void Operateur_IJK_faces_diff::compute_set(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  valeur().compute_set(dvx, dvy, dvz);
}

inline void Operateur_IJK_faces_diff::compute_add(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  valeur().compute_add(dvx, dvy, dvz);
}

inline void Operateur_IJK_faces_diff::initialize(const IJK_Splitting& splitting)
{
  if (!is_cast_)
    typer_diffusion_op("standard");
  diffusion_option_ = diffusion_op_options_[0];
  valeur().initialize(splitting);
}

inline void Operateur_IJK_faces_diff::compute_flux_x_vx(IJK_Field_local_double& resu, const int k_layer)
{
  valeur().compute_flux_x_vx(resu, k_layer);
}

inline void Operateur_IJK_faces_diff::compute_flux_x_vy(IJK_Field_local_double& resu, const int k_layer)
{
  valeur().compute_flux_x_vy(resu, k_layer);
}

inline void Operateur_IJK_faces_diff::compute_flux_x_vz(IJK_Field_local_double& resu, const int k_layer)
{
  valeur().compute_flux_x_vz(resu, k_layer);
}

inline void Operateur_IJK_faces_diff::compute_flux_y_vx(IJK_Field_local_double& resu, const int k_layer)
{
  valeur().compute_flux_y_vx(resu, k_layer);
}

inline void Operateur_IJK_faces_diff::compute_flux_y_vy(IJK_Field_local_double& resu, const int k_layer)
{
  valeur().compute_flux_y_vy(resu, k_layer);
}

inline void Operateur_IJK_faces_diff::compute_flux_y_vz(IJK_Field_local_double& resu, const int k_layer)
{
  valeur().compute_flux_y_vz(resu, k_layer);
}

inline void Operateur_IJK_faces_diff::compute_flux_z_vx(IJK_Field_local_double& resu, const int k_layer)
{
  valeur().compute_flux_z_vx(resu, k_layer);
}

inline void Operateur_IJK_faces_diff::compute_flux_z_vy(IJK_Field_local_double& resu, const int k_layer)
{
  valeur().compute_flux_z_vy(resu, k_layer);
}

inline void Operateur_IJK_faces_diff::compute_flux_z_vz(IJK_Field_local_double& resu, const int k_layer)
{
  valeur().compute_flux_z_vz(resu, k_layer);
}

inline void Operateur_IJK_faces_diff::set_bc(const Boundary_Conditions& bc)
{
  valeur().set_bc(bc);
}

inline void Operateur_IJK_faces_diff::ajouter(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                              IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  valeur().ajouter(vx, vy, vz, dvx, dvy, dvz);
}

inline void Operateur_IJK_faces_diff::set_coeff_x_y_z(const IJK_Field_double& coeff_tensor_xx,
                                                      const IJK_Field_double& coeff_tensor_xy,
                                                      const IJK_Field_double& coeff_tensor_xz,
                                                      const IJK_Field_double& coeff_tensor_yx,
                                                      const IJK_Field_double& coeff_tensor_yy,
                                                      const IJK_Field_double& coeff_tensor_yz,
                                                      const IJK_Field_double& coeff_tensor_zx,
                                                      const IJK_Field_double& coeff_tensor_zy,
                                                      const IJK_Field_double& coeff_tensor_zz)
{
  valeur().set_coeff_x_y_z(coeff_tensor_xx, coeff_tensor_xy, coeff_tensor_xz,
                           coeff_tensor_yx, coeff_tensor_yy, coeff_tensor_yz,
                           coeff_tensor_zx, coeff_tensor_zy, coeff_tensor_zz);
}

inline void Operateur_IJK_faces_diff::set_divergence(const IJK_Field_local_double& divergence)
{
  valeur().set_divergence(divergence);
}

inline void Operateur_IJK_faces_diff::set_nu(const IJK_Field_local_double& nu)
{
  valeur().set_nu(nu);
}

inline void Operateur_IJK_faces_diff::set_uniform_nu(const double& uniform_nu)
{
  valeur().set_uniform_nu(uniform_nu);
}

inline const IJK_Field_local_double& Operateur_IJK_faces_diff::get_v(DIRECTION _DIR_)
{
  return valeur().get_v(_DIR_);
}

inline const IJK_Field_local_double& Operateur_IJK_faces_diff::get_divergence()
{
  return valeur().get_divergence();
}

inline const IJK_Field_local_double& Operateur_IJK_faces_diff::get_coeff_tensor(DIRECTION _COMPO1_, DIRECTION _COMPO2_)
{
  return valeur().get_coeff_tensor(_COMPO1_, _COMPO2_);
}

inline const IJK_Field_local_double& Operateur_IJK_faces_diff::get_nu()
{
  return valeur().get_nu();
}

inline const double& Operateur_IJK_faces_diff::get_uniform_nu()
{
  return valeur().get_uniform_nu();
}

#endif /* Operateur_IJK_faces_diff_included */
