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

#include <Operateur_IJK_elem_diff_base.h>

Implemente_base_sans_constructeur(Operateur_IJK_elem_diff_base_double, "Operateur_IJK_elem_diff_base_double", Operateur_IJK_elem_base_double);

Operateur_IJK_elem_diff_base_double::Operateur_IJK_elem_diff_base_double()
{
  input_field_ = nullptr;
  uniform_lambda_ = nullptr;
  lambda_ = nullptr;

  boundary_flux_kmin_ = boundary_flux_kmax_ = nullptr;

  coeff_field_x_ = nullptr;
  coeff_field_y_ = nullptr;
  coeff_field_z_ = nullptr;

  is_anisotropic_ = false;
  is_vectorial_ = false;
  is_structural_ = false;
  is_uniform_ = false;
  is_corrected_ = false;

  perio_k_ = false;
  is_hess_ = false;

  corrige_flux_ = nullptr;
  indicatrice_ = nullptr;
}

Sortie& Operateur_IJK_elem_diff_base_double::printOn(Sortie& os) const
{
  return os;
}

Entree& Operateur_IJK_elem_diff_base_double::readOn(Entree& is)
{
  return is;
}

const IJK_Field_local_double& Operateur_IJK_elem_diff_base_double::get_model(DIRECTION _DIR_)
{
  assert(is_vectorial_);
  switch(_DIR_)
    {
    case DIRECTION::X:
      return *coeff_field_x_;
    case DIRECTION::Y:
      return *coeff_field_y_;
    case DIRECTION::Z:
      return *coeff_field_z_;
    default:
      Cerr << "Error in Operateur_IJK_elem_diff_base_double::get_model: wrong direction..." << finl;
      Process::exit();
    }
  return *coeff_field_x_;
}

void Operateur_IJK_elem_diff_base_double::initialize(const IJK_Splitting& splitting)
{
  channel_data_.initialize(splitting);
  perio_k_= splitting.get_grid_geometry().get_periodic_flag(DIRECTION_K);
}

void Operateur_IJK_elem_diff_base_double::calculer(const IJK_Field_double& field,
                                                   IJK_Field_double& result,
                                                   const IJK_Field_local_double& boundary_flux_kmin,
                                                   const IJK_Field_local_double& boundary_flux_kmax)
{
  input_field_ = &field;
  boundary_flux_kmin_ = &boundary_flux_kmin;
  boundary_flux_kmax_ = &boundary_flux_kmax;
  compute_set(result);
  input_field_ = 0;
  lambda_ = 0;
  coeff_field_x_ = 0;
  coeff_field_y_ = 0;
  coeff_field_z_ = 0;
  boundary_flux_kmin_ = boundary_flux_kmax_ = 0;
}

void Operateur_IJK_elem_diff_base_double::ajouter(const IJK_Field_double& field,
                                                  IJK_Field_double& result,
                                                  const IJK_Field_local_double& boundary_flux_kmin,
                                                  const IJK_Field_local_double& boundary_flux_kmax)
{
  input_field_ = &field;
  boundary_flux_kmin_ = &boundary_flux_kmin;
  boundary_flux_kmax_ = &boundary_flux_kmax;
  compute_add(result);
  input_field_ = 0;
  lambda_ = 0;
  coeff_field_x_ = 0;
  coeff_field_y_ = 0;
  coeff_field_z_ = 0;
  boundary_flux_kmin_ = boundary_flux_kmax_ = 0;
}

Implemente_instanciable_sans_constructeur(OpDiffUniformIJKScalar_double, "OpDiffUniformIJKScalar_double", Operateur_IJK_elem_diff_base_double);

Sortie& OpDiffUniformIJKScalar_double::printOn(Sortie& os) const
{
  //  Operateur_IJK_elem_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffUniformIJKScalar_double::readOn(Entree& is)
{
  //  Operateur_IJK_elem_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffUniformIJKScalarCorrection_double, "OpDiffUniformIJKScalarCorrection_double", Operateur_IJK_elem_diff_base_double);

Sortie& OpDiffUniformIJKScalarCorrection_double::printOn(Sortie& os) const
{
  //  Operateur_IJK_elem_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffUniformIJKScalarCorrection_double::readOn(Entree& is)
{
  //  Operateur_IJK_elem_diff_base_double::readOn(is);
  return is;
}

void OpDiffUniformIJKScalarCorrection_double::correct_flux(IJK_Field_local_double *const flux, const int k_layer, const int dir)
{
  corrige_flux_->corrige_flux_diff_faceIJ(flux, k_layer, dir);
}

Implemente_instanciable_sans_constructeur(OpDiffIJKScalar_double, "OpDiffIJKScalar_double", Operateur_IJK_elem_diff_base_double);

Sortie& OpDiffIJKScalar_double::printOn(Sortie& os) const
{
  //  Operateur_IJK_elem_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffIJKScalar_double::readOn(Entree& is)
{
  //  Operateur_IJK_elem_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffAnisotropicIJKScalar_double, "OpDiffAnisotropicIJKScalar_double", Operateur_IJK_elem_diff_base_double);

Sortie& OpDiffAnisotropicIJKScalar_double::printOn(Sortie& os) const
{
  //  Operateur_IJK_elem_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffAnisotropicIJKScalar_double::readOn(Entree& is)
{
  //  Operateur_IJK_elem_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffVectorialIJKScalar_double, "OpDiffVectorialIJKScalar_double", Operateur_IJK_elem_diff_base_double);

Sortie& OpDiffVectorialIJKScalar_double::printOn(Sortie& os) const
{
  //  Operateur_IJK_elem_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffVectorialIJKScalar_double::readOn(Entree& is)
{
  //  Operateur_IJK_elem_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffVectorialAnisotropicIJKScalar_double, "OpDiffVectorialAnisotropicIJKScalar_double", Operateur_IJK_elem_diff_base_double);

Sortie& OpDiffVectorialAnisotropicIJKScalar_double::printOn(Sortie& os) const
{
  //  Operateur_IJK_elem_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffVectorialAnisotropicIJKScalar_double::readOn(Entree& is)
{
  //  Operateur_IJK_elem_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffStructuralOnlyIJKScalar_double, "OpDiffStructuralOnlyIJKScalar_double", Operateur_IJK_elem_diff_base_double);

Sortie& OpDiffStructuralOnlyIJKScalar_double::printOn(Sortie& os) const
{
  //  Operateur_IJK_elem_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffStructuralOnlyIJKScalar_double::readOn(Entree& is)
{
  //  Operateur_IJK_elem_diff_base_double::readOn(is);
  return is;
}
