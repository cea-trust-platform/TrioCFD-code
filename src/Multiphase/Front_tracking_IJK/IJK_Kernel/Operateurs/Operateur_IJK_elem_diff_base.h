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

#ifndef Operateur_IJK_elem_diff_base_H
#define Operateur_IJK_elem_diff_base_H

#include <Operateur_IJK_base.h>

class Operateur_IJK_elem_diff_base_double : public Operateur_IJK_elem_base_double
{
  Declare_base_sans_constructeur(Operateur_IJK_elem_diff_base_double);
public:
  Operateur_IJK_elem_diff_base_double();
  virtual void initialize(const IJK_Splitting& splitting) override;
  virtual void set_indicatrice(const IJK_Field_double& indicatrice) { indicatrice_= &indicatrice; };
  virtual void set_corrige_flux(Corrige_flux_FT& corrige_flux) { corrige_flux_ = &corrige_flux; };
  virtual void calculer(const IJK_Field_double& field,
                        IJK_Field_double& result,
                        const IJK_Field_local_double& boundary_flux_kmin,
                        const IJK_Field_local_double& boundary_flux_kmax);
  virtual void ajouter(const IJK_Field_double& field,
                       IJK_Field_double& result,
                       const IJK_Field_local_double& boundary_flux_kmin,
                       const IJK_Field_local_double& boundary_flux_kmax);

  inline void set_uniform_lambda(const double& uniform_lambda) { uniform_lambda_ = &uniform_lambda; };

  inline void set_lambda(const IJK_Field_local_double& lambda) { lambda_ = &lambda; };

  inline void set_coeff_x_y_z(const IJK_Field_local_double& coeff_field_x,
                              const IJK_Field_local_double& coeff_field_y,
                              const IJK_Field_local_double& coeff_field_z)
  {
    coeff_field_x_ = &coeff_field_x;
    coeff_field_y_ = &coeff_field_y;
    coeff_field_z_ = &coeff_field_z;
  }

  inline void compute_flux_x(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::X>(resu,k_layer);
  }
  inline void compute_flux_y(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Y>(resu,k_layer);
  }
  inline void compute_flux_z(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Z>(resu,k_layer);
  }

protected:
  template <DIRECTION _DIR_>
  inline void compute_flux_(IJK_Field_local_double& resu, const int k_layer);

  const IJK_Field_local_double& get_model(DIRECTION _DIR_);

  Operateur_IJK_data_channel channel_data_;
  bool perio_k_;

  // Pointers to input data (set by calculer, used by compute_flux_...)
  const IJK_Field_local_double *input_field_;

  const IJK_Field_local_double *lambda_;
  const double *uniform_lambda_;

  const IJK_Field_local_double *coeff_field_x_;
  const IJK_Field_local_double *coeff_field_y_;
  const IJK_Field_local_double *coeff_field_z_;

  const IJK_Field_local_double *boundary_flux_kmin_;
  const IJK_Field_local_double *boundary_flux_kmax_;

  Corrige_flux_FT *corrige_flux_;
  const IJK_Field_local_double *indicatrice_;

  bool is_anisotropic_;
  bool is_vectorial_;
  bool is_structural_;
  bool is_uniform_;
  bool is_corrected_;
  bool is_hess_;

};

class OpDiffUniformIJKScalar_double : public Operateur_IJK_elem_diff_base_double
{
  Declare_instanciable_sans_constructeur(OpDiffUniformIJKScalar_double);
public:
  OpDiffUniformIJKScalar_double() : Operateur_IJK_elem_diff_base_double() { is_uniform_ = true; }
};

class OpDiffUniformIJKScalarCorrection_double : public Operateur_IJK_elem_diff_base_double
{
  Declare_instanciable_sans_constructeur(OpDiffUniformIJKScalarCorrection_double);
public:
  OpDiffUniformIJKScalarCorrection_double() : Operateur_IJK_elem_diff_base_double() { is_uniform_ = true, is_corrected_ = true; }
private:
  void correct_flux(IJK_Field_local_double *const flux, const int k_layer, const int dir) override;
};

class OpDiffIJKScalar_double : public Operateur_IJK_elem_diff_base_double
{
  Declare_instanciable_sans_constructeur(OpDiffIJKScalar_double);
public:
  OpDiffIJKScalar_double() : Operateur_IJK_elem_diff_base_double() {}
};

class OpDiffAnisotropicIJKScalar_double : public Operateur_IJK_elem_diff_base_double
{
  Declare_instanciable_sans_constructeur(OpDiffAnisotropicIJKScalar_double);
public:
  OpDiffAnisotropicIJKScalar_double() : Operateur_IJK_elem_diff_base_double() { is_anisotropic_ = true; }
};

class OpDiffVectorialIJKScalar_double : public Operateur_IJK_elem_diff_base_double
{
  Declare_instanciable_sans_constructeur(OpDiffVectorialIJKScalar_double);
public:
  OpDiffVectorialIJKScalar_double() : Operateur_IJK_elem_diff_base_double() { is_vectorial_ = true; }
};

class OpDiffVectorialAnisotropicIJKScalar_double : public Operateur_IJK_elem_diff_base_double
{
  Declare_instanciable_sans_constructeur(OpDiffVectorialAnisotropicIJKScalar_double);
public:
  OpDiffVectorialAnisotropicIJKScalar_double() : Operateur_IJK_elem_diff_base_double() { is_vectorial_ = true, is_anisotropic_ = true; }
};

class OpDiffStructuralOnlyIJKScalar_double : public Operateur_IJK_elem_diff_base_double
{
  Declare_instanciable_sans_constructeur(OpDiffStructuralOnlyIJKScalar_double);
public:
  OpDiffStructuralOnlyIJKScalar_double() : Operateur_IJK_elem_diff_base_double() { is_structural_ = true; }
};

#include <Operateur_IJK_elem_diff_base.tpp>

#endif
