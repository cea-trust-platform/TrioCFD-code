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

#include <OpDiffTurbIJKScalar.h>

Implemente_base_sans_constructeur(OpDiffIJKScalarGeneric_double, "OpDiffIJKScalarGeneric_double", Operateur_IJK_elem_base_double);

Sortie& OpDiffIJKScalarGeneric_double::printOn(Sortie& os) const
{
  //  Operateur_IJK_elem_base_double::printOn(os);
  return os;
}

Entree& OpDiffIJKScalarGeneric_double::readOn(Entree& is)
{
  //  Operateur_IJK_elem_base_double::readOn(is);
  return is;
}

OpDiffIJKScalarGeneric_double::OpDiffIJKScalarGeneric_double()
{
  input_field_ = 0;
  uniform_lambda_ = 0;
  lambda_ = 0;

  boundary_flux_kmin_ = boundary_flux_kmax_ = 0;

  coeff_field_x_ =0;
  coeff_field_y_ =0;
  coeff_field_z_ =0;

  is_anisotropic_ = false;
  is_vectorial_ = false;
  is_structural_ = false;
  is_uniform_ = false;
  is_corrected_ = false;

  //  lambda_vector_x_ = 0;
  //  lambda_vector_y_ = 0;
  //  lambda_vector_z_ = 0;
  //
  //  structural_model_x_ = 0;
  //  structural_model_y_ = 0;
  //  structural_model_z_ = 0;

}

const IJK_Field_local_double& OpDiffIJKScalarGeneric_double::get_model(DIRECTION _DIR_)
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
      Cerr << "Error in OpDiffIJKScalarGeneric_double::get_model: wrong direction..." << finl;
      Process::exit();
    }
  return *coeff_field_x_;
}

//const IJK_Field_local_double& OpDiffIJKScalarGeneric_double::get_lambda_vectorial(DIRECTION _DIR_)
//{
//  assert(is_vectorial_);
//  switch(_DIR_)
//    {
//    case DIRECTION::X:
//      return *lambda_vector_x_;
//    case DIRECTION::Y:
//      return *lambda_vector_y_;
//    case DIRECTION::Z:
//      return *lambda_vector_z_;
//    default:
//      Cerr << "Error in OpDiffIJKScalarGeneric_double::get_lambda_vectorial: wrong direction..." << finl;
//      Process::exit();
//    }
//  return *lambda_vector_x_;
//}
//
//const IJK_Field_local_double& OpDiffIJKScalarGeneric_double::get_structural_model(DIRECTION _DIR_)
//{
//  assert(is_structural_);
//  switch(_DIR_)
//    {
//    case DIRECTION::X:
//      return *structural_model_x_;
//    case DIRECTION::Y:
//      return *structural_model_y_;
//    case DIRECTION::Z:
//      return *structural_model_z_;
//    default:
//      Cerr << "Error in OpDiffIJKScalarGeneric_double::get_strucutral_model: wrong direction..." << finl;
//      Process::exit();
//    }
//  return *structural_model_x_;
//}

void OpDiffIJKScalarGeneric_double::initialize(const IJK_Splitting& splitting)
{
  channel_data_.initialize(splitting);
  perio_k_= splitting.get_grid_geometry().get_periodic_flag(DIRECTION_K);
}

void OpDiffIJKScalarGeneric_double::calculer(const IJK_Field_double& field,
                                             IJK_Field_double& result,
                                             const IJK_Field_local_double& boundary_flux_kmin,
                                             const IJK_Field_local_double& boundary_flux_kmax)
{
  input_field_ = &field;
  boundary_flux_kmin_ = &boundary_flux_kmin;
  boundary_flux_kmax_ = &boundary_flux_kmax;
  compute_set(result);
  input_field_ = 0;
//  uniform_lambda_ = 0;
  lambda_ = 0;
  coeff_field_x_ = 0;
  coeff_field_y_ = 0;
  coeff_field_z_ = 0;
  boundary_flux_kmin_ = boundary_flux_kmax_ = 0;
}

void OpDiffIJKScalarGeneric_double::ajouter(const IJK_Field_double& field,
                                            IJK_Field_double& result,
                                            const IJK_Field_local_double& boundary_flux_kmin,
                                            const IJK_Field_local_double& boundary_flux_kmax)
{
  input_field_ = &field;
  boundary_flux_kmin_ = &boundary_flux_kmin;
  boundary_flux_kmax_ = &boundary_flux_kmax;
  compute_add(result);
  input_field_ = 0;
//  uniform_lambda_ = 0;
  lambda_ = 0;
  coeff_field_x_ = 0;
  coeff_field_y_ = 0;
  coeff_field_z_ = 0;
  boundary_flux_kmin_ = boundary_flux_kmax_ = 0;
}

Implemente_instanciable_sans_constructeur(OpDiffUniformIJKScalar_double, "OpDiffUniformIJKScalar_double", OpDiffIJKScalarGeneric_double);

Sortie& OpDiffUniformIJKScalar_double::printOn(Sortie& os) const
{
  //  OpDiffIJKScalarGeneric_double::printOn(os);
  return os;
}

Entree& OpDiffUniformIJKScalar_double::readOn(Entree& is)
{
  //  OpDiffIJKScalarGeneric_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffUniformIJKScalarCorrection_double, "OpDiffUniformIJKScalarCorrection_double", OpDiffIJKScalarGeneric_double);

Sortie& OpDiffUniformIJKScalarCorrection_double::printOn(Sortie& os) const
{
  //  OpDiffIJKScalarGeneric_double::printOn(os);
  return os;
}

Entree& OpDiffUniformIJKScalarCorrection_double::readOn(Entree& is)
{
  //  OpDiffIJKScalarGeneric_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffIJKScalar_double, "OpDiffIJKScalar_double", OpDiffIJKScalarGeneric_double);

Sortie& OpDiffIJKScalar_double::printOn(Sortie& os) const
{
  //  OpDiffIJKScalarGeneric_double::printOn(os);
  return os;
}

Entree& OpDiffIJKScalar_double::readOn(Entree& is)
{
  //  OpDiffIJKScalarGeneric_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffAnisotropicIJKScalar_double, "OpDiffAnisotropicIJKScalar_double", OpDiffIJKScalarGeneric_double);

Sortie& OpDiffAnisotropicIJKScalar_double::printOn(Sortie& os) const
{
  //  OpDiffIJKScalarGeneric_double::printOn(os);
  return os;
}

Entree& OpDiffAnisotropicIJKScalar_double::readOn(Entree& is)
{
  //  OpDiffIJKScalarGeneric_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffVectorialIJKScalar_double, "OpDiffVectorialIJKScalar_double", OpDiffIJKScalarGeneric_double);

Sortie& OpDiffVectorialIJKScalar_double::printOn(Sortie& os) const
{
  //  OpDiffIJKScalarGeneric_double::printOn(os);
  return os;
}

Entree& OpDiffVectorialIJKScalar_double::readOn(Entree& is)
{
  //  OpDiffIJKScalarGeneric_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffVectorialAnisotropicIJKScalar_double, "OpDiffVectorialAnisotropicIJKScalar_double", OpDiffIJKScalarGeneric_double);

Sortie& OpDiffVectorialAnisotropicIJKScalar_double::printOn(Sortie& os) const
{
  //  OpDiffIJKScalarGeneric_double::printOn(os);
  return os;
}

Entree& OpDiffVectorialAnisotropicIJKScalar_double::readOn(Entree& is)
{
  //  OpDiffIJKScalarGeneric_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffIJKScalarStructuralOnly_double, "OpDiffIJKScalarStructuralOnly_double", OpDiffIJKScalarGeneric_double);

Sortie& OpDiffIJKScalarStructuralOnly_double::printOn(Sortie& os) const
{
  //  OpDiffIJKScalarGeneric_double::printOn(os);
  return os;
}

Entree& OpDiffIJKScalarStructuralOnly_double::readOn(Entree& is)
{
  //  OpDiffIJKScalarGeneric_double::readOn(is);
  return is;
}

//void OpDiffUniformIJKScalar_double::calculer(const IJK_Field_double& field,
//																							const double& lambda,
//																							IJK_Field_double& result,
//																							const IJK_Field_local_double& boundary_flux_kmin,
//																							const IJK_Field_local_double& boundary_flux_kmax)
//{
//  input_field_ = &field;
//  uniform_lambda_ = lambda;
//  boundary_flux_kmin_ = &boundary_flux_kmin;
//  boundary_flux_kmax_ = &boundary_flux_kmax;
//  compute_set(result);
//  input_field_ = 0;
//  uniform_lambda_ = 0;
//  boundary_flux_kmin_ = boundary_flux_kmax_ = 0;
//}
//
//void OpDiffIJKScalar_double::calculer(const IJK_Field_double& field,
//                                      const IJK_Field_double& lambda,
//                                      IJK_Field_double& result,
//                                      const IJK_Field_local_double& boundary_flux_kmin,
//                                      const IJK_Field_local_double& boundary_flux_kmax)
//{
//  input_field_ = &field;
//  lambda_ = &lambda;
//  boundary_flux_kmin_ = &boundary_flux_kmin;
//  boundary_flux_kmax_ = &boundary_flux_kmax;
//  compute_set(result);
//  input_field_ = 0;
//  lambda_ = 0;
//  boundary_flux_kmin_ = boundary_flux_kmax_ = 0;
//}
//
//void OpDiffAnisotropicIJKScalar_double::calculer(const IJK_Field_double& field,
//                                                 const IJK_Field_double& lambda,
//                                                 IJK_Field_double& result,
//                                                 const IJK_Field_local_double& boundary_flux_kmin,
//                                                 const IJK_Field_local_double& boundary_flux_kmax)
//{
//  input_field_ = &field;
//  lambda_ = &lambda;
//  boundary_flux_kmin_ = &boundary_flux_kmin;
//  boundary_flux_kmax_ = &boundary_flux_kmax;
//  compute_set(result);
//  input_field_ = 0;
//  lambda_ = 0;
//  boundary_flux_kmin_ = boundary_flux_kmax_ = 0;
//}
//
//void OpDiffAnisotropicIJKScalar_double::ajouter(const IJK_Field_double& field,
//                                                const IJK_Field_double& lambda,
//                                                IJK_Field_double& result,
//                                                const IJK_Field_local_double& boundary_flux_kmin,
//                                                const IJK_Field_local_double& boundary_flux_kmax)
//{
//  input_field_ = &field;
//  lambda_ = &lambda;
//  boundary_flux_kmin_ = &boundary_flux_kmin;
//  boundary_flux_kmax_ = &boundary_flux_kmax;
//  compute_add(result);
//  input_field_ = 0;
//  lambda_ = 0;
//  boundary_flux_kmin_ = boundary_flux_kmax_ = 0;
//}
//
//
//void OpDiffVectorialIJKScalar_double::calculer(const IJK_Field_double& field,
//                                               const IJK_Field_double& lambda_vector_x,
//                                               const IJK_Field_double& lambda_vector_y,
//                                               const IJK_Field_double& lambda_vector_z,
//                                               IJK_Field_double& result,
//                                               const IJK_Field_local_double& boundary_flux_kmin,
//                                               const IJK_Field_local_double& boundary_flux_kmax)
//{
//  input_field_ = &field;
//  lambda_vector_x_ = &lambda_vector_x;
//  lambda_vector_y_ = &lambda_vector_y;
//  lambda_vector_z_ = &lambda_vector_z;
//  boundary_flux_kmin_ = &boundary_flux_kmin;
//  boundary_flux_kmax_ = &boundary_flux_kmax;
//  compute_set(result);
//  input_field_ = 0;
//  lambda_vector_x_ = 0;
//  lambda_vector_y_ = 0;
//  lambda_vector_z_ = 0;
//  boundary_flux_kmin_ = boundary_flux_kmax_ = 0;
//}
//
//
//
//void OpDiffVectorialAnisotropicIJKScalar_double::calculer(const IJK_Field_double& field,
//                                                          const IJK_Field_double& lambda_vector_x,
//                                                          const IJK_Field_double& lambda_vector_y,
//                                                          const IJK_Field_double& lambda_vector_z,
//                                                          IJK_Field_double& result,
//                                                          const IJK_Field_local_double& boundary_flux_kmin,
//                                                          const IJK_Field_local_double& boundary_flux_kmax)
//{
//  input_field_ = &field;
//  lambda_vector_x_ = &lambda_vector_x;
//  lambda_vector_y_ = &lambda_vector_y;
//  lambda_vector_z_ = &lambda_vector_z;
//  boundary_flux_kmin_ = &boundary_flux_kmin;
//  boundary_flux_kmax_ = &boundary_flux_kmax;
//  compute_set(result);
//  input_field_ = 0;
//  lambda_vector_x_ = 0;
//  lambda_vector_y_ = 0;
//  lambda_vector_z_ = 0;
//  boundary_flux_kmin_ = boundary_flux_kmax_ = 0;
//}
//
//void OpDiffVectorialAnisotropicIJKScalar_double::ajouter(const IJK_Field_double& field,
//                                                         const IJK_Field_double& lambda_vector_x,
//                                                         const IJK_Field_double& lambda_vector_y,
//                                                         const IJK_Field_double& lambda_vector_z,
//                                                         IJK_Field_double& result,
//                                                         const IJK_Field_local_double& boundary_flux_kmin,
//                                                         const IJK_Field_local_double& boundary_flux_kmax)
//{
//  input_field_ = &field;
//  lambda_vector_x_ = &lambda_vector_x;
//  lambda_vector_y_ = &lambda_vector_y;
//  lambda_vector_z_ = &lambda_vector_z;
//  boundary_flux_kmin_ = &boundary_flux_kmin;
//  boundary_flux_kmax_ = &boundary_flux_kmax;
//  compute_add(result);
//  input_field_ = 0;
//  lambda_vector_x_ = 0;
//  lambda_vector_y_ = 0;
//  lambda_vector_z_ = 0;
//  boundary_flux_kmin_ = boundary_flux_kmax_ = 0;
//}
//
//void OpDiffIJKScalarStructuralOnly_double::calculer(const IJK_Field_double& field,
//                                                    const IJK_Field_double& structural_model_x,
//                                                    const IJK_Field_double& structural_model_y,
//                                                    const IJK_Field_double& structural_model_z,
//                                                    IJK_Field_double& result,
//                                                    const IJK_Field_local_double& boundary_flux_kmin,
//                                                    const IJK_Field_local_double& boundary_flux_kmax)
//{
//  input_field_ = &field;
//  structural_model_x_ = &structural_model_x;
//  structural_model_y_ = &structural_model_y;
//  structural_model_z_ = &structural_model_z;
//  boundary_flux_kmin_ = &boundary_flux_kmin;
//  boundary_flux_kmax_ = &boundary_flux_kmax;
//  compute_set(result);
//  input_field_ = 0;
//  structural_model_x_ = 0;
//  structural_model_y_ = 0;
//  structural_model_z_ = 0;
//  boundary_flux_kmin_ = boundary_flux_kmax_ = 0;
//}
//
//void OpDiffIJKScalarStructuralOnly_double::ajouter(const IJK_Field_double& field,
//                                                   const IJK_Field_double& structural_model_x,
//                                                   const IJK_Field_double& structural_model_y,
//                                                   const IJK_Field_double& structural_model_z,
//                                                   IJK_Field_double& result,
//                                                   const IJK_Field_local_double& boundary_flux_kmin,
//                                                   const IJK_Field_local_double& boundary_flux_kmax)
//{
//  input_field_ = &field;
//  structural_model_x_ = &structural_model_x;
//  structural_model_y_ = &structural_model_y;
//  structural_model_z_ = &structural_model_z;
//  boundary_flux_kmin_ = &boundary_flux_kmin;
//  boundary_flux_kmax_ = &boundary_flux_kmax;
//  compute_add(result);
//  input_field_ = 0;
//  structural_model_x_ = 0;
//  structural_model_y_ = 0;
//  structural_model_z_ = 0;
//  boundary_flux_kmin_ = boundary_flux_kmax_ = 0;
//}
