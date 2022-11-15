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

#include <OpDiffTurbIJK.h>
#include <IJK_Splitting.h>
#include <stat_counters.h>


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

OpDiffIJKFacesGeneric_double::OpDiffIJKFacesGeneric_double()
{
  vx_ = 0;
  vy_ = 0;
  vz_ = 0;
  molecular_nu_ = 0;

  turbulent_nu_ = 0;
  turbulent_k_energy_ = 0;
  divergence_ = 0;

  molecular_nu_tensor_xx_ = 0;
  molecular_nu_tensor_xy_ = 0;
  molecular_nu_tensor_xz_ = 0;
  molecular_nu_tensor_yx_ = 0;
  molecular_nu_tensor_yy_ = 0;
  molecular_nu_tensor_yz_ = 0;
  molecular_nu_tensor_zx_ = 0;
  molecular_nu_tensor_zy_ = 0;
  molecular_nu_tensor_zz_ = 0;

  structural_model_xx_ = 0;
  structural_model_xy_ = 0;
  structural_model_xz_ = 0;
  structural_model_yx_ = 0;
  structural_model_yy_ = 0;
  structural_model_yz_ = 0;
  structural_model_zx_ = 0;
  structural_model_zy_ = 0;
  structural_model_zz_ = 0;

  is_turb_= false;
  is_anisotropic_ = false;
  with_divergence_= false;
  with_transpose_= false;
  is_tensorial_= false;
  is_structural_= false;


}


const IJK_Field_local_double& OpDiffIJKFacesGeneric_double::get_v(DIRECTION _DIR_)
{
  switch(_DIR_)
    {
    case DIRECTION::X:
      return *vx_;
      break;
    case DIRECTION::Y:
      return *vy_;
      break;
    case DIRECTION::Z:
      return *vz_;
      break;
    default:
      Cerr << "Error in OpDiffIJKFacesGeneric_double::get_v: wrong direction..." << finl;
      Process::exit();
      // for compilation only...
      return *vx_;
    }
}
const IJK_Field_local_double& OpDiffIJKFacesGeneric_double::get_molecular_nu_tensor(DIRECTION _COMPO1_, DIRECTION _COMPO2_)
{
  assert(is_tensorial_);
  switch(_COMPO1_)
    {
    case DIRECTION::X:
      {
        if(_COMPO2_ == DIRECTION::X)
          return *molecular_nu_tensor_xx_;
        if(_COMPO2_ == DIRECTION::Y)
          return *molecular_nu_tensor_xy_;
        if(_COMPO2_ == DIRECTION::Z)
          return *molecular_nu_tensor_xz_;
        break;
      }
    case DIRECTION::Y:
      {
        if(_COMPO2_ == DIRECTION::X)
          return *molecular_nu_tensor_yx_;
        if(_COMPO2_ == DIRECTION::Y)
          return *molecular_nu_tensor_yy_;
        if(_COMPO2_ == DIRECTION::Z)
          return *molecular_nu_tensor_yz_;
        break;
      }
    case DIRECTION::Z:
      {
        if(_COMPO2_ == DIRECTION::X)
          return *molecular_nu_tensor_zx_;
        if(_COMPO2_ == DIRECTION::Y)
          return *molecular_nu_tensor_zy_;
        if(_COMPO2_ == DIRECTION::Z)
          return *molecular_nu_tensor_zz_;
        break;
      }
    default:
      Cerr << "Error in OpDiffTensorial_base_double::get_molecular_nu: wrong direction..." << finl;
      Process::exit();
    }

  // for compilation only...
  return *molecular_nu_tensor_xx_;
}

const IJK_Field_local_double& OpDiffIJKFacesGeneric_double::get_molecular_nu()
{
  assert(!is_tensorial_);
  return *molecular_nu_;
}

const IJK_Field_local_double& OpDiffIJKFacesGeneric_double::get_structural_model(DIRECTION _COMPO1_, DIRECTION _COMPO2_)
{
  assert(is_structural_);

  switch(_COMPO1_)
    {
    case DIRECTION::X:
      {
        if(_COMPO2_ == DIRECTION::X)
          return *structural_model_xx_;
        if(_COMPO2_ == DIRECTION::Y)
          return *structural_model_xy_;
        if(_COMPO2_ == DIRECTION::Z)
          return *structural_model_xz_;
        break;
      }
    case DIRECTION::Y:
      {
        if(_COMPO2_ == DIRECTION::X)
          return *structural_model_yx_;
        if(_COMPO2_ == DIRECTION::Y)
          return *structural_model_yy_;
        if(_COMPO2_ == DIRECTION::Z)
          return *structural_model_yz_;
        break;
      }
    case DIRECTION::Z:
      {
        if(_COMPO2_ == DIRECTION::X)
          return *structural_model_zx_;
        if(_COMPO2_ == DIRECTION::Y)
          return *structural_model_zy_;
        if(_COMPO2_ == DIRECTION::Z)
          return *structural_model_zz_;
        break;
      }
    default:
      Cerr << "Error in OpDiffStructuralOnlyZeroatwallIJK_double::get_structural_model: wrong direction..." << finl;
      Process::exit();
    }

  // for compilation only...
  return *structural_model_zx_;
}

const IJK_Field_local_double& OpDiffIJKFacesGeneric_double::get_turbulent_nu()
{
  assert(is_turb_);
  return *turbulent_nu_;
}
const IJK_Field_local_double& OpDiffIJKFacesGeneric_double::get_turbulent_k_energy()
{
  assert(is_turb_);
  return *turbulent_k_energy_;
}

const IJK_Field_local_double& OpDiffIJKFacesGeneric_double::get_divergence()
{
  assert(with_divergence_);
  return *divergence_;
}

void OpDiffTurbIJK_double::ajouter(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                   const IJK_Field_double& molecular_nu,
                                   const IJK_Field_double& turbulent_nu,
                                   const IJK_Field_double& turbulent_k_energy,
                                   IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_ = &molecular_nu;

  turbulent_nu_ = &turbulent_nu;
  turbulent_k_energy_ = &turbulent_k_energy;
  compute_add(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}
void OpDiffTurbIJK_double::calculer(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                    const IJK_Field_double& molecular_nu,
                                    const IJK_Field_double& turbulent_nu,
                                    const IJK_Field_double& turbulent_k_energy,
                                    IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_ = &molecular_nu;

  turbulent_nu_ = &turbulent_nu;
  turbulent_k_energy_ = &turbulent_k_energy;
  compute_set(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}


void OpDiffIJK_double::ajouter(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                               const IJK_Field_double& molecular_nu,
                               IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_ = &molecular_nu;

  compute_add(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}
void OpDiffIJK_double::calculer(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                const IJK_Field_double& molecular_nu,
                                IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_ = &molecular_nu;

  compute_set(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}

void OpDiffStdWithLaminarTransposeIJK_double::ajouter(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                      const IJK_Field_double& molecular_nu,
                                                      IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_ = &molecular_nu;

  compute_add(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}
void OpDiffStdWithLaminarTransposeIJK_double::calculer(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                       const IJK_Field_double& molecular_nu,
                                                       IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_ = &molecular_nu;

  compute_set(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}


void OpDiffStdWithLaminarTransposeAndDivergenceIJK_double::ajouter(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                                   const IJK_Field_double& molecular_nu,
                                                                   const IJK_Field_double& divergence,
                                                                   IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_ = &molecular_nu;

  divergence_ = &divergence;
  compute_add(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}
void OpDiffStdWithLaminarTransposeAndDivergenceIJK_double::calculer(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                                    const IJK_Field_double& molecular_nu,
                                                                    const IJK_Field_double& divergence,
                                                                    IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_ = &molecular_nu;

  divergence_ = &divergence;
  compute_set(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}


void OpDiffAnisotropicIJK_double::ajouter(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                          const IJK_Field_double& molecular_nu,
                                          IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_ = &molecular_nu;

  compute_add(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}
void OpDiffAnisotropicIJK_double::calculer(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                           const IJK_Field_double& molecular_nu,
                                           IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_ = &molecular_nu;

  compute_set(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}


void OpDiffStdWithLaminarTransposeAnisotropicIJK_double::ajouter(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                                 const IJK_Field_double& molecular_nu,
                                                                 IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_ = &molecular_nu;

  compute_add(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}

void OpDiffStdWithLaminarTransposeAnisotropicIJK_double::calculer(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                                  const IJK_Field_double& molecular_nu,
                                                                  IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_ = &molecular_nu;

  compute_set(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}

void OpDiffStdWithLaminarTransposeAndDivergenceAnisotropicIJK_double::ajouter(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                                              const IJK_Field_double& molecular_nu,
                                                                              const IJK_Field_double& divergence,
                                                                              IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_ = &molecular_nu;

  divergence_ = &divergence;
  compute_add(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}
void OpDiffStdWithLaminarTransposeAndDivergenceAnisotropicIJK_double::calculer(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                                               const IJK_Field_double& molecular_nu,
                                                                               const IJK_Field_double& divergence,
                                                                               IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_ = &molecular_nu;

  divergence_ = &divergence;
  compute_set(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}

void OpDiffTensorialZeroatwallIJK_double::ajouter(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                  const IJK_Field_double& molecular_nu_tensor_xx,
                                                  const IJK_Field_double& molecular_nu_tensor_xy,
                                                  const IJK_Field_double& molecular_nu_tensor_xz,
                                                  const IJK_Field_double& molecular_nu_tensor_yx,
                                                  const IJK_Field_double& molecular_nu_tensor_yy,
                                                  const IJK_Field_double& molecular_nu_tensor_yz,
                                                  const IJK_Field_double& molecular_nu_tensor_zx,
                                                  const IJK_Field_double& molecular_nu_tensor_zy,
                                                  const IJK_Field_double& molecular_nu_tensor_zz,
                                                  IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_tensor_xx_ = &molecular_nu_tensor_xx;
  molecular_nu_tensor_xy_ = &molecular_nu_tensor_xy;
  molecular_nu_tensor_xz_ = &molecular_nu_tensor_xz;
  molecular_nu_tensor_yx_ = &molecular_nu_tensor_yx;
  molecular_nu_tensor_yy_ = &molecular_nu_tensor_yy;
  molecular_nu_tensor_yz_ = &molecular_nu_tensor_yz;
  molecular_nu_tensor_zx_ = &molecular_nu_tensor_zx;
  molecular_nu_tensor_zy_ = &molecular_nu_tensor_zy;
  molecular_nu_tensor_zz_ = &molecular_nu_tensor_zz;

  compute_add(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}
void OpDiffTensorialZeroatwallIJK_double::calculer(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                   const IJK_Field_double& molecular_nu_tensor_xx,
                                                   const IJK_Field_double& molecular_nu_tensor_xy,
                                                   const IJK_Field_double& molecular_nu_tensor_xz,
                                                   const IJK_Field_double& molecular_nu_tensor_yx,
                                                   const IJK_Field_double& molecular_nu_tensor_yy,
                                                   const IJK_Field_double& molecular_nu_tensor_yz,
                                                   const IJK_Field_double& molecular_nu_tensor_zx,
                                                   const IJK_Field_double& molecular_nu_tensor_zy,
                                                   const IJK_Field_double& molecular_nu_tensor_zz,
                                                   IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_tensor_xx_ = &molecular_nu_tensor_xx;
  molecular_nu_tensor_xy_ = &molecular_nu_tensor_xy;
  molecular_nu_tensor_xz_ = &molecular_nu_tensor_xz;
  molecular_nu_tensor_yx_ = &molecular_nu_tensor_yx;
  molecular_nu_tensor_yy_ = &molecular_nu_tensor_yy;
  molecular_nu_tensor_yz_ = &molecular_nu_tensor_yz;
  molecular_nu_tensor_zx_ = &molecular_nu_tensor_zx;
  molecular_nu_tensor_zy_ = &molecular_nu_tensor_zy;
  molecular_nu_tensor_zz_ = &molecular_nu_tensor_zz;

  compute_set(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}

void OpDiffStdWithLaminarTransposeTensorialZeroatwallIJK_double::ajouter(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                                         const IJK_Field_double& molecular_nu_tensor_xx,
                                                                         const IJK_Field_double& molecular_nu_tensor_xy,
                                                                         const IJK_Field_double& molecular_nu_tensor_xz,
                                                                         const IJK_Field_double& molecular_nu_tensor_yx,
                                                                         const IJK_Field_double& molecular_nu_tensor_yy,
                                                                         const IJK_Field_double& molecular_nu_tensor_yz,
                                                                         const IJK_Field_double& molecular_nu_tensor_zx,
                                                                         const IJK_Field_double& molecular_nu_tensor_zy,
                                                                         const IJK_Field_double& molecular_nu_tensor_zz,
                                                                         IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_tensor_xx_ = &molecular_nu_tensor_xx;
  molecular_nu_tensor_xy_ = &molecular_nu_tensor_xy;
  molecular_nu_tensor_xz_ = &molecular_nu_tensor_xz;
  molecular_nu_tensor_yx_ = &molecular_nu_tensor_yx;
  molecular_nu_tensor_yy_ = &molecular_nu_tensor_yy;
  molecular_nu_tensor_yz_ = &molecular_nu_tensor_yz;
  molecular_nu_tensor_zx_ = &molecular_nu_tensor_zx;
  molecular_nu_tensor_zy_ = &molecular_nu_tensor_zy;
  molecular_nu_tensor_zz_ = &molecular_nu_tensor_zz;

  compute_add(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}
void OpDiffStdWithLaminarTransposeTensorialZeroatwallIJK_double::calculer(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                                          const IJK_Field_double& molecular_nu_tensor_xx,
                                                                          const IJK_Field_double& molecular_nu_tensor_xy,
                                                                          const IJK_Field_double& molecular_nu_tensor_xz,
                                                                          const IJK_Field_double& molecular_nu_tensor_yx,
                                                                          const IJK_Field_double& molecular_nu_tensor_yy,
                                                                          const IJK_Field_double& molecular_nu_tensor_yz,
                                                                          const IJK_Field_double& molecular_nu_tensor_zx,
                                                                          const IJK_Field_double& molecular_nu_tensor_zy,
                                                                          const IJK_Field_double& molecular_nu_tensor_zz,
                                                                          IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_tensor_xx_ = &molecular_nu_tensor_xx;
  molecular_nu_tensor_xy_ = &molecular_nu_tensor_xy;
  molecular_nu_tensor_xz_ = &molecular_nu_tensor_xz;
  molecular_nu_tensor_yx_ = &molecular_nu_tensor_yx;
  molecular_nu_tensor_yy_ = &molecular_nu_tensor_yy;
  molecular_nu_tensor_yz_ = &molecular_nu_tensor_yz;
  molecular_nu_tensor_zx_ = &molecular_nu_tensor_zx;
  molecular_nu_tensor_zy_ = &molecular_nu_tensor_zy;
  molecular_nu_tensor_zz_ = &molecular_nu_tensor_zz;

  compute_set(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}


void OpDiffStdWithLaminarTransposeAndDivergenceTensorialZeroatwallIJK_double::ajouter(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                                                      const IJK_Field_double& molecular_nu_tensor_xx,
                                                                                      const IJK_Field_double& molecular_nu_tensor_xy,
                                                                                      const IJK_Field_double& molecular_nu_tensor_xz,
                                                                                      const IJK_Field_double& molecular_nu_tensor_yx,
                                                                                      const IJK_Field_double& molecular_nu_tensor_yy,
                                                                                      const IJK_Field_double& molecular_nu_tensor_yz,
                                                                                      const IJK_Field_double& molecular_nu_tensor_zx,
                                                                                      const IJK_Field_double& molecular_nu_tensor_zy,
                                                                                      const IJK_Field_double& molecular_nu_tensor_zz,
                                                                                      const IJK_Field_double& divergence,
                                                                                      IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_tensor_xx_ = &molecular_nu_tensor_xx;
  molecular_nu_tensor_xy_ = &molecular_nu_tensor_xy;
  molecular_nu_tensor_xz_ = &molecular_nu_tensor_xz;
  molecular_nu_tensor_yx_ = &molecular_nu_tensor_yx;
  molecular_nu_tensor_yy_ = &molecular_nu_tensor_yy;
  molecular_nu_tensor_yz_ = &molecular_nu_tensor_yz;
  molecular_nu_tensor_zx_ = &molecular_nu_tensor_zx;
  molecular_nu_tensor_zy_ = &molecular_nu_tensor_zy;
  molecular_nu_tensor_zz_ = &molecular_nu_tensor_zz;

  divergence_ = &divergence;
  compute_add(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}
void OpDiffStdWithLaminarTransposeAndDivergenceTensorialZeroatwallIJK_double::calculer(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                                                       const IJK_Field_double& molecular_nu_tensor_xx,
                                                                                       const IJK_Field_double& molecular_nu_tensor_xy,
                                                                                       const IJK_Field_double& molecular_nu_tensor_xz,
                                                                                       const IJK_Field_double& molecular_nu_tensor_yx,
                                                                                       const IJK_Field_double& molecular_nu_tensor_yy,
                                                                                       const IJK_Field_double& molecular_nu_tensor_yz,
                                                                                       const IJK_Field_double& molecular_nu_tensor_zx,
                                                                                       const IJK_Field_double& molecular_nu_tensor_zy,
                                                                                       const IJK_Field_double& molecular_nu_tensor_zz,
                                                                                       const IJK_Field_double& divergence,
                                                                                       IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_tensor_xx_ = &molecular_nu_tensor_xx;
  molecular_nu_tensor_xy_ = &molecular_nu_tensor_xy;
  molecular_nu_tensor_xz_ = &molecular_nu_tensor_xz;
  molecular_nu_tensor_yx_ = &molecular_nu_tensor_yx;
  molecular_nu_tensor_yy_ = &molecular_nu_tensor_yy;
  molecular_nu_tensor_yz_ = &molecular_nu_tensor_yz;
  molecular_nu_tensor_zx_ = &molecular_nu_tensor_zx;
  molecular_nu_tensor_zy_ = &molecular_nu_tensor_zy;
  molecular_nu_tensor_zz_ = &molecular_nu_tensor_zz;

  divergence_ = &divergence;
  compute_set(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}


void OpDiffTensorialAnisotropicZeroatwallIJK_double::ajouter(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                             const IJK_Field_double& molecular_nu_tensor_xx,
                                                             const IJK_Field_double& molecular_nu_tensor_xy,
                                                             const IJK_Field_double& molecular_nu_tensor_xz,
                                                             const IJK_Field_double& molecular_nu_tensor_yx,
                                                             const IJK_Field_double& molecular_nu_tensor_yy,
                                                             const IJK_Field_double& molecular_nu_tensor_yz,
                                                             const IJK_Field_double& molecular_nu_tensor_zx,
                                                             const IJK_Field_double& molecular_nu_tensor_zy,
                                                             const IJK_Field_double& molecular_nu_tensor_zz,
                                                             IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_tensor_xx_ = &molecular_nu_tensor_xx;
  molecular_nu_tensor_xy_ = &molecular_nu_tensor_xy;
  molecular_nu_tensor_xz_ = &molecular_nu_tensor_xz;
  molecular_nu_tensor_yx_ = &molecular_nu_tensor_yx;
  molecular_nu_tensor_yy_ = &molecular_nu_tensor_yy;
  molecular_nu_tensor_yz_ = &molecular_nu_tensor_yz;
  molecular_nu_tensor_zx_ = &molecular_nu_tensor_zx;
  molecular_nu_tensor_zy_ = &molecular_nu_tensor_zy;
  molecular_nu_tensor_zz_ = &molecular_nu_tensor_zz;

  compute_add(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}
void OpDiffTensorialAnisotropicZeroatwallIJK_double::calculer(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                              const IJK_Field_double& molecular_nu_tensor_xx,
                                                              const IJK_Field_double& molecular_nu_tensor_xy,
                                                              const IJK_Field_double& molecular_nu_tensor_xz,
                                                              const IJK_Field_double& molecular_nu_tensor_yx,
                                                              const IJK_Field_double& molecular_nu_tensor_yy,
                                                              const IJK_Field_double& molecular_nu_tensor_yz,
                                                              const IJK_Field_double& molecular_nu_tensor_zx,
                                                              const IJK_Field_double& molecular_nu_tensor_zy,
                                                              const IJK_Field_double& molecular_nu_tensor_zz,
                                                              IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_tensor_xx_ = &molecular_nu_tensor_xx;
  molecular_nu_tensor_xy_ = &molecular_nu_tensor_xy;
  molecular_nu_tensor_xz_ = &molecular_nu_tensor_xz;
  molecular_nu_tensor_yx_ = &molecular_nu_tensor_yx;
  molecular_nu_tensor_yy_ = &molecular_nu_tensor_yy;
  molecular_nu_tensor_yz_ = &molecular_nu_tensor_yz;
  molecular_nu_tensor_zx_ = &molecular_nu_tensor_zx;
  molecular_nu_tensor_zy_ = &molecular_nu_tensor_zy;
  molecular_nu_tensor_zz_ = &molecular_nu_tensor_zz;

  compute_set(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}

void OpDiffStdWithLaminarTransposeTensorialAnisotropicZeroatwallIJK_double::ajouter(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                                                    const IJK_Field_double& molecular_nu_tensor_xx,
                                                                                    const IJK_Field_double& molecular_nu_tensor_xy,
                                                                                    const IJK_Field_double& molecular_nu_tensor_xz,
                                                                                    const IJK_Field_double& molecular_nu_tensor_yx,
                                                                                    const IJK_Field_double& molecular_nu_tensor_yy,
                                                                                    const IJK_Field_double& molecular_nu_tensor_yz,
                                                                                    const IJK_Field_double& molecular_nu_tensor_zx,
                                                                                    const IJK_Field_double& molecular_nu_tensor_zy,
                                                                                    const IJK_Field_double& molecular_nu_tensor_zz,
                                                                                    IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_tensor_xx_ = &molecular_nu_tensor_xx;
  molecular_nu_tensor_xy_ = &molecular_nu_tensor_xy;
  molecular_nu_tensor_xz_ = &molecular_nu_tensor_xz;
  molecular_nu_tensor_yx_ = &molecular_nu_tensor_yx;
  molecular_nu_tensor_yy_ = &molecular_nu_tensor_yy;
  molecular_nu_tensor_yz_ = &molecular_nu_tensor_yz;
  molecular_nu_tensor_zx_ = &molecular_nu_tensor_zx;
  molecular_nu_tensor_zy_ = &molecular_nu_tensor_zy;
  molecular_nu_tensor_zz_ = &molecular_nu_tensor_zz;

  compute_add(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}
void OpDiffStdWithLaminarTransposeTensorialAnisotropicZeroatwallIJK_double::calculer(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                                                     const IJK_Field_double& molecular_nu_tensor_xx,
                                                                                     const IJK_Field_double& molecular_nu_tensor_xy,
                                                                                     const IJK_Field_double& molecular_nu_tensor_xz,
                                                                                     const IJK_Field_double& molecular_nu_tensor_yx,
                                                                                     const IJK_Field_double& molecular_nu_tensor_yy,
                                                                                     const IJK_Field_double& molecular_nu_tensor_yz,
                                                                                     const IJK_Field_double& molecular_nu_tensor_zx,
                                                                                     const IJK_Field_double& molecular_nu_tensor_zy,
                                                                                     const IJK_Field_double& molecular_nu_tensor_zz,
                                                                                     IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_tensor_xx_ = &molecular_nu_tensor_xx;
  molecular_nu_tensor_xy_ = &molecular_nu_tensor_xy;
  molecular_nu_tensor_xz_ = &molecular_nu_tensor_xz;
  molecular_nu_tensor_yx_ = &molecular_nu_tensor_yx;
  molecular_nu_tensor_yy_ = &molecular_nu_tensor_yy;
  molecular_nu_tensor_yz_ = &molecular_nu_tensor_yz;
  molecular_nu_tensor_zx_ = &molecular_nu_tensor_zx;
  molecular_nu_tensor_zy_ = &molecular_nu_tensor_zy;
  molecular_nu_tensor_zz_ = &molecular_nu_tensor_zz;

  compute_set(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}

void OpDiffStdWithLaminarTransposeAndDivergenceTensorialAnisotropicZeroatwallIJK_double::ajouter(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                                                                 const IJK_Field_double& molecular_nu_tensor_xx,
                                                                                                 const IJK_Field_double& molecular_nu_tensor_xy,
                                                                                                 const IJK_Field_double& molecular_nu_tensor_xz,
                                                                                                 const IJK_Field_double& molecular_nu_tensor_yx,
                                                                                                 const IJK_Field_double& molecular_nu_tensor_yy,
                                                                                                 const IJK_Field_double& molecular_nu_tensor_yz,
                                                                                                 const IJK_Field_double& molecular_nu_tensor_zx,
                                                                                                 const IJK_Field_double& molecular_nu_tensor_zy,
                                                                                                 const IJK_Field_double& molecular_nu_tensor_zz,
                                                                                                 const IJK_Field_double& divergence,
                                                                                                 IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_tensor_xx_ = &molecular_nu_tensor_xx;
  molecular_nu_tensor_xy_ = &molecular_nu_tensor_xy;
  molecular_nu_tensor_xz_ = &molecular_nu_tensor_xz;
  molecular_nu_tensor_yx_ = &molecular_nu_tensor_yx;
  molecular_nu_tensor_yy_ = &molecular_nu_tensor_yy;
  molecular_nu_tensor_yz_ = &molecular_nu_tensor_yz;
  molecular_nu_tensor_zx_ = &molecular_nu_tensor_zx;
  molecular_nu_tensor_zy_ = &molecular_nu_tensor_zy;
  molecular_nu_tensor_zz_ = &molecular_nu_tensor_zz;

  divergence_ = &divergence;
  compute_add(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}
void OpDiffStdWithLaminarTransposeAndDivergenceTensorialAnisotropicZeroatwallIJK_double::calculer(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                                                                  const IJK_Field_double& molecular_nu_tensor_xx,
                                                                                                  const IJK_Field_double& molecular_nu_tensor_xy,
                                                                                                  const IJK_Field_double& molecular_nu_tensor_xz,
                                                                                                  const IJK_Field_double& molecular_nu_tensor_yx,
                                                                                                  const IJK_Field_double& molecular_nu_tensor_yy,
                                                                                                  const IJK_Field_double& molecular_nu_tensor_yz,
                                                                                                  const IJK_Field_double& molecular_nu_tensor_zx,
                                                                                                  const IJK_Field_double& molecular_nu_tensor_zy,
                                                                                                  const IJK_Field_double& molecular_nu_tensor_zz,
                                                                                                  const IJK_Field_double& divergence,
                                                                                                  IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;

  molecular_nu_tensor_xx_ = &molecular_nu_tensor_xx;
  molecular_nu_tensor_xy_ = &molecular_nu_tensor_xy;
  molecular_nu_tensor_xz_ = &molecular_nu_tensor_xz;
  molecular_nu_tensor_yx_ = &molecular_nu_tensor_yx;
  molecular_nu_tensor_yy_ = &molecular_nu_tensor_yy;
  molecular_nu_tensor_yz_ = &molecular_nu_tensor_yz;
  molecular_nu_tensor_zx_ = &molecular_nu_tensor_zx;
  molecular_nu_tensor_zy_ = &molecular_nu_tensor_zy;
  molecular_nu_tensor_zz_ = &molecular_nu_tensor_zz;

  divergence_ = &divergence;
  compute_set(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}

void OpDiffStructuralOnlyZeroatwallIJK_double::ajouter(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                       const IJK_Field_double& structural_model_xx,
                                                       const IJK_Field_double& structural_model_xy,
                                                       const IJK_Field_double& structural_model_xz,
                                                       const IJK_Field_double& structural_model_yx,
                                                       const IJK_Field_double& structural_model_yy,
                                                       const IJK_Field_double& structural_model_yz,
                                                       const IJK_Field_double& structural_model_zx,
                                                       const IJK_Field_double& structural_model_zy,
                                                       const IJK_Field_double& structural_model_zz,
                                                       IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;


  structural_model_xx_ = &structural_model_xx;
  structural_model_xy_ = &structural_model_xy;
  structural_model_xz_ = &structural_model_xz;
  structural_model_yx_ = &structural_model_yx;
  structural_model_yy_ = &structural_model_yy;
  structural_model_yz_ = &structural_model_yz;
  structural_model_zx_ = &structural_model_zx;
  structural_model_zy_ = &structural_model_zy;
  structural_model_zz_ = &structural_model_zz;
  compute_add(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}
void OpDiffStructuralOnlyZeroatwallIJK_double::calculer(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                        const IJK_Field_double& structural_model_xx,
                                                        const IJK_Field_double& structural_model_xy,
                                                        const IJK_Field_double& structural_model_xz,
                                                        const IJK_Field_double& structural_model_yx,
                                                        const IJK_Field_double& structural_model_yy,
                                                        const IJK_Field_double& structural_model_yz,
                                                        const IJK_Field_double& structural_model_zx,
                                                        const IJK_Field_double& structural_model_zy,
                                                        const IJK_Field_double& structural_model_zz,
                                                        IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  statistiques().begin_count(diffusion_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;


  structural_model_xx_ = &structural_model_xx;
  structural_model_xy_ = &structural_model_xy;
  structural_model_xz_ = &structural_model_xz;
  structural_model_yx_ = &structural_model_yx;
  structural_model_yy_ = &structural_model_yy;
  structural_model_yz_ = &structural_model_yz;
  structural_model_zx_ = &structural_model_zx;
  structural_model_zy_ = &structural_model_zy;
  structural_model_zz_ = &structural_model_zz;
  compute_set(dvx, dvy, dvz);
  statistiques().end_count(diffusion_counter_);

}

