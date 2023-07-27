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

#include <IJK_Splitting.h>
#include <stat_counters.h>
#include <Operateur_IJK_faces_diff_base.h>

/*
 * Options for CASE
 *
 *          Yes_Turb: the flux is turbulent_mu * (grad u + grad^T u)  -  2/3 * k  +  molecular_mu * grad u)'
 *          Yes_M_Grad: the flux is 'molecular_mu * grad u'
 *          Yes_M_Trans: the flux is 'molecular_mu * (grad u + grad^T u)'
 *          Yes_M_Div: the flux is 'molecular_mu * (grad u + grad^T u - 2/3 * div u * Id)'
 *          Yes_M_GradAnisotropic: the flux is 'molecular_mu^a * grad^a u' where (grad^a)_i = Delta_i (grad)_i
 *          Yes_M_TransAnisotropic: the flux is 'molecular_mu^a * (grad^a u + grad^a^T u)' where (grad^a)_i = Delta_i (grad)_i
 *          Yes_M_DivAnisotropic: the flux is 'molecular_mu^a * (grad^a u + grad^a^T u - 2/3 * div^a u * Id)' where (grad^a)_i = Delta_i (grad)_i
 *          Yes_M_GradTensorial: the flux is 'molecular_mu_tensor * grad u'
 *          Yes_M_TransTensorial: the flux is 'molecular_mu_tensor * (grad u + grad^T u)'
 *          Yes_M_DivTensorial: the flux is 'molecular_mu_tensor * (grad u + grad^T u - 2/3 * div u * Id)'
 *          Yes_M_GradTensorialAnisotropic: the flux is 'molecular_mu_tensor^a * grad^a u' where (grad^a)_i = Delta_i (grad)_i
 *          Yes_M_TransTensorialAnisotropic: the flux is 'molecular_mu_tensor^a * (grad^a u + grad^a^T u)' where (grad^a)_i = Delta_i (grad)_i
 *          Yes_M_DivTensorialAnisotropic: the flux is 'molecular_mu_tensor^a * (grad^a u + grad^a^T u - 2/3 * div^a u * Id)' where (grad^a)_i = Delta_i (grad)_i
 *          Yes_M_Struct: the flux is 'structural_model'
 */

Implemente_base_sans_constructeur(Operateur_IJK_faces_diff_base_double, "Operateur_IJK_faces_diff_base_double", Operateur_IJK_faces_base_double);

Sortie& Operateur_IJK_faces_diff_base_double::printOn(Sortie& os) const
{
  return os;
}

Entree& Operateur_IJK_faces_diff_base_double::readOn(Entree& is)
{
  return is;
}

Operateur_IJK_faces_diff_base_double::Operateur_IJK_faces_diff_base_double()
{
  vx_ = 0;
  vy_ = 0;
  vz_ = 0;

  uniform_nu_=0;
  nu_ = 0;
  divergence_ = 0;

  coeff_tensor_xx_ = 0;
  coeff_tensor_xy_ = 0;
  coeff_tensor_xz_ = 0;
  coeff_tensor_yx_ = 0;
  coeff_tensor_yy_ = 0;
  coeff_tensor_yz_ = 0;
  coeff_tensor_zx_ = 0;
  coeff_tensor_zy_ = 0;
  coeff_tensor_zz_ = 0;

  is_uniform_ = false;
  is_turb_= false;
  is_anisotropic_ = false;
  with_divergence_= false;
  with_transpose_= false;
  is_tensorial_= false;
  is_structural_= false;

  perio_k_=false;
}

void Operateur_IJK_faces_diff_base_double::set_coeff_x_y_z(const IJK_Field_double& coeff_tensor_xx,
                                                           const IJK_Field_double& coeff_tensor_xy,
                                                           const IJK_Field_double& coeff_tensor_xz,
                                                           const IJK_Field_double& coeff_tensor_yx,
                                                           const IJK_Field_double& coeff_tensor_yy,
                                                           const IJK_Field_double& coeff_tensor_yz,
                                                           const IJK_Field_double& coeff_tensor_zx,
                                                           const IJK_Field_double& coeff_tensor_zy,
                                                           const IJK_Field_double& coeff_tensor_zz)
{
  coeff_tensor_xx_ = &coeff_tensor_xx;
  coeff_tensor_xy_ = &coeff_tensor_xy;
  coeff_tensor_xz_ = &coeff_tensor_xz;
  coeff_tensor_yx_ = &coeff_tensor_yx;
  coeff_tensor_yy_ = &coeff_tensor_yy;
  coeff_tensor_yz_ = &coeff_tensor_yz;
  coeff_tensor_zx_ = &coeff_tensor_zx;
  coeff_tensor_zy_ = &coeff_tensor_zy;
  coeff_tensor_zz_ = &coeff_tensor_zz;
}

const IJK_Field_local_double& Operateur_IJK_faces_diff_base_double::get_v(DIRECTION _DIR_)
{
  switch(_DIR_)
    {
    case DIRECTION::X:
      return *vx_;
    case DIRECTION::Y:
      return *vy_;
    case DIRECTION::Z:
      return *vz_;
    default:
      Cerr << "Error in Operateur_IJK_faces_diff_base_double::get_v: wrong direction..." << finl;
      Process::exit();
    }
  return *vx_;
}

const IJK_Field_local_double& Operateur_IJK_faces_diff_base_double::get_coeff_tensor(DIRECTION _COMPO1_, DIRECTION _COMPO2_)
{
  assert(is_structural_ || is_tensorial_);

  switch(_COMPO1_)
    {
    case DIRECTION::X:
      {
        if(_COMPO2_ == DIRECTION::X)
          return *coeff_tensor_xx_;
        if(_COMPO2_ == DIRECTION::Y)
          return *coeff_tensor_xy_;
        if(_COMPO2_ == DIRECTION::Z)
          return *coeff_tensor_xz_;
        break;
      }
    case DIRECTION::Y:
      {
        if(_COMPO2_ == DIRECTION::X)
          return *coeff_tensor_yx_;
        if(_COMPO2_ == DIRECTION::Y)
          return *coeff_tensor_yy_;
        if(_COMPO2_ == DIRECTION::Z)
          return *coeff_tensor_yz_;
        break;
      }
    case DIRECTION::Z:
      {
        if(_COMPO2_ == DIRECTION::X)
          return *coeff_tensor_zx_;
        if(_COMPO2_ == DIRECTION::Y)
          return *coeff_tensor_zy_;
        if(_COMPO2_ == DIRECTION::Z)
          return *coeff_tensor_zz_;
        break;
      }
    default:
      Cerr << "Error in OpDiffStructuralOnlyZeroatwallIJK_double::get_structural_model: wrong direction..." << finl;
      Process::exit();
    }

  // for compilation only...
  return *coeff_tensor_xx_;
}

const IJK_Field_local_double& Operateur_IJK_faces_diff_base_double::get_nu()
{
  assert(!is_tensorial_);
  return *nu_;
}

const double& Operateur_IJK_faces_diff_base_double::get_uniform_nu()
{
  assert(is_uniform_);
  return *uniform_nu_;
}

const IJK_Field_local_double& Operateur_IJK_faces_diff_base_double::get_divergence()
{
  assert(with_divergence_);
  return *divergence_;
}

void Operateur_IJK_faces_diff_base_double::ajouter(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                   IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);
  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;
  compute_add(dvx, dvy, dvz);
  nu_=0;
  divergence_ = 0;
  coeff_tensor_xx_ = 0;
  coeff_tensor_xy_ = 0;
  coeff_tensor_xz_ = 0;
  coeff_tensor_yx_ = 0;
  coeff_tensor_yy_ = 0;
  coeff_tensor_yz_ = 0;
  coeff_tensor_zx_ = 0;
  coeff_tensor_zy_ = 0;
  coeff_tensor_zz_ = 0;
  statistiques().end_count(diffusion_counter_);
}

void Operateur_IJK_faces_diff_base_double::calculer(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                    IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);
  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;
  nu_=0;
  divergence_ = 0;
  coeff_tensor_xx_ = 0;
  coeff_tensor_xy_ = 0;
  coeff_tensor_xz_ = 0;
  coeff_tensor_yx_ = 0;
  coeff_tensor_yy_ = 0;
  coeff_tensor_yz_ = 0;
  coeff_tensor_zx_ = 0;
  coeff_tensor_zy_ = 0;
  coeff_tensor_zz_ = 0;
  compute_set(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);
}

/*
 * Definitions of the subclasses
 */

Implemente_instanciable_sans_constructeur(OpDiffIJK_double, "OpDiffIJK_double", Operateur_IJK_faces_diff_base_double);

Sortie& OpDiffIJK_double::printOn(Sortie& os) const
{
  Operateur_IJK_faces_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffIJK_double::readOn(Entree& is)
{
  Operateur_IJK_faces_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffTurbIJK_double, "OpDiffTurbIJK_double", Operateur_IJK_faces_diff_base_double);

Sortie& OpDiffTurbIJK_double::printOn(Sortie& os) const
{
  Operateur_IJK_faces_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffTurbIJK_double::readOn(Entree& is)
{
  Operateur_IJK_faces_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffStdWithLaminarTransposeIJK_double, "OpDiffStdWithLaminarTransposeIJK_double", Operateur_IJK_faces_diff_base_double);

Sortie& OpDiffStdWithLaminarTransposeIJK_double::printOn(Sortie& os) const
{
  Operateur_IJK_faces_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffStdWithLaminarTransposeIJK_double::readOn(Entree& is)
{
  Operateur_IJK_faces_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffStdWithLaminarTransposeAndDivergenceIJK_double, "OpDiffStdWithLaminarTransposeAndDivergenceIJK_double", Operateur_IJK_faces_diff_base_double);

Sortie& OpDiffStdWithLaminarTransposeAndDivergenceIJK_double::printOn(Sortie& os) const
{
  Operateur_IJK_faces_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffStdWithLaminarTransposeAndDivergenceIJK_double::readOn(Entree& is)
{
  Operateur_IJK_faces_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffAnisotropicIJK_double, "OpDiffAnisotropicIJK_double", Operateur_IJK_faces_diff_base_double);

Sortie& OpDiffAnisotropicIJK_double::printOn(Sortie& os) const
{
  Operateur_IJK_faces_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffAnisotropicIJK_double::readOn(Entree& is)
{
  Operateur_IJK_faces_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffStdWithLaminarTransposeAnisotropicIJK_double, "OpDiffStdWithLaminarTransposeAnisotropicIJK_double", Operateur_IJK_faces_diff_base_double);

Sortie& OpDiffStdWithLaminarTransposeAnisotropicIJK_double::printOn(Sortie& os) const
{
  Operateur_IJK_faces_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffStdWithLaminarTransposeAnisotropicIJK_double::readOn(Entree& is)
{
  Operateur_IJK_faces_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffStdWithLaminarTransposeAndDivergenceTensorialAnisotropicZeroatwallIJK_double, "OpDiffStdWithLaminarTransposeAndDivergenceTensorialAnisotropicZeroatwallIJK_double", Operateur_IJK_faces_diff_base_double);

Sortie& OpDiffStdWithLaminarTransposeAndDivergenceTensorialAnisotropicZeroatwallIJK_double::printOn(Sortie& os) const
{
  Operateur_IJK_faces_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffStdWithLaminarTransposeAndDivergenceTensorialAnisotropicZeroatwallIJK_double::readOn(Entree& is)
{
  Operateur_IJK_faces_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffStructuralOnlyZeroatwallIJK_double, "OpDiffStructuralOnlyZeroatwallIJK_double", Operateur_IJK_faces_diff_base_double);

Sortie& OpDiffStructuralOnlyZeroatwallIJK_double::printOn(Sortie& os) const
{
  Operateur_IJK_faces_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffStructuralOnlyZeroatwallIJK_double::readOn(Entree& is)
{
  Operateur_IJK_faces_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffStdWithLaminarTransposeTensorialAnisotropicZeroatwallIJK_double, "OpDiffStdWithLaminarTransposeTensorialAnisotropicZeroatwallIJK_double", Operateur_IJK_faces_diff_base_double);

Sortie& OpDiffStdWithLaminarTransposeTensorialAnisotropicZeroatwallIJK_double::printOn(Sortie& os) const
{
  Operateur_IJK_faces_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffStdWithLaminarTransposeTensorialAnisotropicZeroatwallIJK_double::readOn(Entree& is)
{
  Operateur_IJK_faces_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffStdWithLaminarTransposeAndDivergenceAnisotropicIJK_double, "OpDiffStdWithLaminarTransposeAndDivergenceAnisotropicIJK_double", Operateur_IJK_faces_diff_base_double);

Sortie& OpDiffStdWithLaminarTransposeAndDivergenceAnisotropicIJK_double::printOn(Sortie& os) const
{
  Operateur_IJK_faces_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffStdWithLaminarTransposeAndDivergenceAnisotropicIJK_double::readOn(Entree& is)
{
  Operateur_IJK_faces_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffTensorialZeroatwallIJK_double, "OpDiffTensorialZeroatwallIJK_double", Operateur_IJK_faces_diff_base_double);

Sortie& OpDiffTensorialZeroatwallIJK_double::printOn(Sortie& os) const
{
  Operateur_IJK_faces_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffTensorialZeroatwallIJK_double::readOn(Entree& is)
{
  Operateur_IJK_faces_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffStdWithLaminarTransposeTensorialZeroatwallIJK_double, "OpDiffStdWithLaminarTransposeTensorialZeroatwallIJK_double", Operateur_IJK_faces_diff_base_double);

Sortie& OpDiffStdWithLaminarTransposeTensorialZeroatwallIJK_double::printOn(Sortie& os) const
{
  Operateur_IJK_faces_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffStdWithLaminarTransposeTensorialZeroatwallIJK_double::readOn(Entree& is)
{
  Operateur_IJK_faces_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffTensorialAnisotropicZeroatwallIJK_double, "OpDiffTensorialAnisotropicZeroatwallIJK_double", Operateur_IJK_faces_diff_base_double);

Sortie& OpDiffTensorialAnisotropicZeroatwallIJK_double::printOn(Sortie& os) const
{
  Operateur_IJK_faces_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffTensorialAnisotropicZeroatwallIJK_double::readOn(Entree& is)
{
  Operateur_IJK_faces_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffStdWithLaminarTransposeAndDivergenceTensorialZeroatwallIJK_double, "OpDiffStdWithLaminarTransposeAndDivergenceTensorialZeroatwallIJK_double", Operateur_IJK_faces_diff_base_double);

Sortie& OpDiffStdWithLaminarTransposeAndDivergenceTensorialZeroatwallIJK_double::printOn(Sortie& os) const
{
  Operateur_IJK_faces_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffStdWithLaminarTransposeAndDivergenceTensorialZeroatwallIJK_double::readOn(Entree& is)
{
  Operateur_IJK_faces_diff_base_double::readOn(is);
  return is;
}
