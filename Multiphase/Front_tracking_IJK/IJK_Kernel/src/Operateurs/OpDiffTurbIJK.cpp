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


/*

void OpDiffIJK_double::flux_loop_(DIRECTION _DIR_, DIRECTION _VCOMPO_, IJK_Field_local_double& resu, int k_layer, int top_wall, int bottom_wall)
{
  Cerr << "IN OPDIFFIJK FLUX LOOP!" << finl;
  const int idir = (int)_DIR_;
  const int nx = _DIR_ == DIRECTION::X ? get_v(_VCOMPO_).ni() + 1 : get_v(_VCOMPO_).ni() ;
  const int ny = _DIR_ == DIRECTION::Y ? get_v(_VCOMPO_).nj() + 1 : get_v(_VCOMPO_).nj();
  const int icompo = (int)_VCOMPO_;

  ConstIJK_double_ptr molecular_nu(*molecular_nu_, 0, 0, k_layer);
  ConstIJK_double_ptr vCOMPO_ptr(get_v(_VCOMPO_), 0, 0, k_layer);
  ConstIJK_double_ptr vDIR_ptr(get_v(_DIR_), 0, 0, k_layer);

  // Result (fluxes in direction DIR for component COMPO of the convected field)
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);

  const double surface = - channel_data_.get_surface(k_layer, icompo, idir);
  const double inv_distance_COMPO = channel_data_.inv_distance_for_gradient(k_layer, idir, icompo);
  const double inv_distance_DIR = channel_data_.inv_distance_for_gradient(k_layer, icompo, idir);

  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double flux = 0.;
          if(_DIR_ == _VCOMPO_)
            {
              // recoding of Eval_Dift_VDF_var_Face::flux_fa7_elem
              Simd_double m_nu, m_nu_dummy;
              // flux is between left face and center face, hence, it is centered on the "left" cell
              molecular_nu.get_left_center(_DIR_,i, m_nu, m_nu_dummy);
              Simd_double v1, v2;
              vCOMPO_ptr.get_left_center(_VCOMPO_,i, v1, v2);
              Simd_double tau = (v2 - v1) * inv_distance_COMPO;
              flux = m_nu * tau * surface;
            }
          else
            {
              // Interpolate diffusion coefficient from values at elements:
              Simd_double m_nu1, m_nu2, m_nu3, m_nu4;
              molecular_nu.get_left_center_c1c2(_DIR_, _VCOMPO_, i, m_nu1, m_nu2, m_nu3, m_nu4);
              // recoding of Eval_Dift_VDF_var_Face::flux_arete_paroi

              double mult_coeff = 0.25;
              if(bottom_wall)
                {
                  // bottom wall (z=0)
                  // nu1 and nu2 are "left" in direction z, hence in the wall:
                  m_nu1 = 0., m_nu2 = 0.;
                  mult_coeff = 0.5;
                }
              if(top_wall)
                {
                  // top wall (z=zmax)
                  // nu3 and nu4 are "center" in direction z, hence in the wall:
                  m_nu3 = 0., m_nu4 = 0.;
                  mult_coeff = 0.5;
                }
              Simd_double m_nu = (m_nu1 + m_nu2 + m_nu3 + m_nu4) * mult_coeff;

              // gradient in direction DIR of component VCOMPO
              Simd_double v3, v4;
              vCOMPO_ptr.get_left_center(_DIR_,i, v3, v4);
              if(top_wall)
                v4 = 0.;
              if(bottom_wall)  // bottom wall (z=0), v3 is left, hence in the wall)
                v3 = 0.;
              Simd_double tau = (v4 - v3) * inv_distance_DIR; //inv_distance_COMPO : distance to wall already contains 1/2
              flux = tau * surface * m_nu;
            }
          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      // instructions to perform to jump to next row
      vCOMPO_ptr.next_j();
      if(_DIR_ != _VCOMPO_)
        vDIR_ptr.next_j();
      molecular_nu.next_j();
      resu_ptr.next_j();
    }
}

/*! @brief compute fluxes in direction x for velocity component x for the layer of fluxes k_layer
 *
 *  Diffusion
 *
 */
void OpDiffStdWithLaminarTransposeIJK_double::flux_loop_(DIRECTION _DIR_, DIRECTION _VCOMPO_, IJK_Field_local_double& resu, int k_layer, int top_wall, int bottom_wall)
{
  const int idir = (int)_DIR_;
  const int nx = _DIR_ == DIRECTION::X ? get_v(_VCOMPO_).ni() + 1 : get_v(_VCOMPO_).ni() ;
  const int ny = _DIR_ == DIRECTION::Y ? get_v(_VCOMPO_).nj() + 1 : get_v(_VCOMPO_).nj();
  const int icompo = (int)_VCOMPO_;

  const double surface = - channel_data_.get_surface(k_layer, icompo, idir);
  const double inv_distance_COMPO = channel_data_.inv_distance_for_gradient(k_layer, idir, icompo);
  const double inv_distance_DIR = channel_data_.inv_distance_for_gradient(k_layer, icompo, idir);

  // Velocity field
  ConstIJK_double_ptr vCOMPO_ptr(get_v(_VCOMPO_), 0, 0, k_layer);
  ConstIJK_double_ptr vDIR_ptr(get_v(_DIR_), 0, 0, k_layer);
  ConstIJK_double_ptr molecular_nu(*molecular_nu_, 0, 0, k_layer);

  // Result (fluxes in direction DIR for component VCOMPO of the convected field)
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);

  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double flux = 0.;
          if(_DIR_ == _VCOMPO_)
            {
              // recoding of Eval_Dift_VDF_var_Face::flux_fa7_elem
              Simd_double m_nu, m_nu_dummy;
              // flux is between left face and center face, hence, it is centered on the "left" cell
              molecular_nu.get_left_center(_DIR_, i, m_nu, m_nu_dummy);
              Simd_double v1, v2;
              vCOMPO_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              Simd_double tau = (v2 - v1) * inv_distance_DIR;
              flux = m_nu * 2. * tau * surface; // The factor 2. is because tau_tr = tau
            }
          else
            {
              // Interpolate diffusion coefficient from values at elements:
              Simd_double m_nu1, m_nu2, m_nu3, m_nu4;
              molecular_nu.get_left_center_c1c2(_DIR_, _VCOMPO_, i, m_nu1, m_nu2, m_nu3, m_nu4);
              double mult_coeff = 0.25;
              if(bottom_wall)
                {
                  // bottom wall (z=0)
                  // nu1 and nu2 are "left" in direction z, hence in the wall:
                  m_nu1 = 0., m_nu2 = 0.;
                  mult_coeff = 0.5;
                }
              if(top_wall)
                {
                  // top wall (z=zmax)
                  // nu3 and nu4 are "center" in direction z, hence in the wall:
                  m_nu3 = 0., m_nu4 = 0.;
                  mult_coeff = 0.5;
                }

              // recoding of Eval_Dift_VDF_var_Face::flux_arete_interne
              Simd_double m_nu = (m_nu1 + m_nu2 + m_nu3 + m_nu4) * mult_coeff;
              // gradient in direction DIR of component COMPO
              Simd_double v3, v4;
              vCOMPO_ptr.get_left_center(_DIR_, i, v3, v4);
              if(top_wall)
                v4 = 0.;
              if(bottom_wall)  // bottom wall (z=0), v3 is left, hence in the wall)
                v3 = 0.;
              Simd_double tau = (v4 - v3) * inv_distance_DIR;

              // gradient in direction COMPO of component DIR
              Simd_double v1, v2;
              vDIR_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              Simd_double tau_tr = (v2 - v1) * inv_distance_COMPO;
              if(top_wall || bottom_wall)
                tau_tr = 0.;
              flux = m_nu * (tau + tau_tr) * surface;
            }
          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      // instructions to perform to jump to next row
      vCOMPO_ptr.next_j();
      if(_DIR_ != _VCOMPO_)
        vDIR_ptr.next_j();
      molecular_nu.next_j();
      resu_ptr.next_j();
    }

}


/*! @brief compute fluxes in direction DIR for velocity component VCOMPO for the layer of fluxes k_layer
 *
 *  Diffusion
 *
 */
void OpDiffStdWithLaminarTransposeAndDivergenceIJK_double::flux_loop_(DIRECTION _DIR_, DIRECTION _VCOMPO_, IJK_Field_local_double& resu, int k_layer, int top_wall, int bottom_wall)
{
  const int idir = (int)_DIR_;
  const int nx = _DIR_ == DIRECTION::X ? get_v(_VCOMPO_).ni() + 1 : get_v(_VCOMPO_).ni() ;
  const int ny = _DIR_ == DIRECTION::Y ? get_v(_VCOMPO_).nj() + 1 : get_v(_VCOMPO_).nj();
  const int icompo = (int)_VCOMPO_;

  const double surface = - channel_data_.get_surface(k_layer, icompo, idir);
  const double inv_distance_COMPO = channel_data_.inv_distance_for_gradient(k_layer, idir, icompo);
  const double inv_distance_DIR = channel_data_.inv_distance_for_gradient(k_layer, icompo, idir);

  // Velocity field
  ConstIJK_double_ptr vCOMPO_ptr(get_v(_VCOMPO_), 0, 0, k_layer);
  ConstIJK_double_ptr vDIR_ptr(get_v(_DIR_), 0, 0, k_layer);

  ConstIJK_double_ptr molecular_nu(*molecular_nu_, 0, 0, k_layer);
  ConstIJK_double_ptr div_ptr(*divergence_, 0, 0, k_layer);

  // Result (fluxes in direction x for component x of the convected field)
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);

  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double flux = 0.;

          if(_DIR_ == _VCOMPO_)
            {
              // recoding of Eval_Dift_VDF_var_Face::flux_fa7_elem
              Simd_double m_nu, m_nu_dummy;
              // flux is between left face and center face, hence, it is centered on the "left" cell
              molecular_nu.get_left_center(_DIR_, i, m_nu, m_nu_dummy);
              Simd_double v1, v2;
              vCOMPO_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              Simd_double tau = (v2 - v1) * inv_distance_DIR;
              Simd_double v5, v6_dummy;
              div_ptr.get_left_center(_DIR_,i, v5, v6_dummy);
              flux = m_nu * surface * (2. * tau - 0.66666666666666666 * v5); // The factor 2. is because tau_tr = tau
            }
          else
            {
              // Interpolate diffusion coefficient from values at elements:
              Simd_double m_nu1, m_nu2, m_nu3, m_nu4;
              molecular_nu.get_left_center_c1c2(_DIR_, _VCOMPO_, i, m_nu1, m_nu2, m_nu3, m_nu4);

              double mult_coeff = 0.25;
                // bottom wall (z=0)
                // nu1 and nu2 are "left" in direction z, hence in the wall:
                if(bottom_wall)
                  {
                    m_nu1 = 0., m_nu2 = 0.;
                    mult_coeff = 0.5;
                  }
                if(top_wall)
                  {
                    m_nu3 = 0., m_nu4 = 0.;
                    mult_coeff = 0.5;
                  }
              // recoding of Eval_Dift_VDF_var_Face::flux_arete_interne
              Simd_double m_nu = (m_nu1 + m_nu2 + m_nu3 + m_nu4) * mult_coeff;
              // gradient in direction DIR of component COMPO
              Simd_double v3, v4;
              vCOMPO_ptr.get_left_center(_DIR_,i, v3, v4);
              if(top_wall)
                           v4 = 0.;
                         if(bottom_wall)  // bottom wall (z=0), v3 is left, hence in the wall)
                           v3 = 0.;
              Simd_double tau = (v4 - v3) * inv_distance_DIR;

              Simd_double v1, v2;
                      vDIR_ptr.get_left_center(_VCOMPO_, i, v1, v2);
                      Simd_double tau_tr = (v2 - v1) * inv_distance_COMPO;
                      if(top_wall || bottom_wall)
                        tau_tr = 0.;
                      flux = m_nu * (tau + tau_tr) * surface;

            }

          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      // instructions to perform to jump to next row
      vCOMPO_ptr.next_j();
      if(_DIR_ != _VCOMPO_)
        vDIR_ptr.next_j();
      molecular_nu.next_j();
      div_ptr.next_j();
      resu_ptr.next_j();
    }
}


/*! @brief compute fluxes in direction x for velocity component x for the layer of fluxes k_layer
 *
 *  Diffusion
 *
 */
void OpDiffAnisotropicIJK_double::flux_loop_(DIRECTION _DIR_, DIRECTION _VCOMPO_, IJK_Field_local_double& resu, int k_layer, int top_wall, int bottom_wall)
{
  const int idir = (int)_DIR_;
  const int nx = _DIR_ == DIRECTION::X ? get_v(_VCOMPO_).ni() + 1 : get_v(_VCOMPO_).ni() ;
  const int ny = _DIR_ == DIRECTION::Y ? get_v(_VCOMPO_).nj() + 1 : get_v(_VCOMPO_).nj();
  const int icompo = (int)_VCOMPO_;

  const double surface = - channel_data_.get_surface(k_layer, icompo, idir);
  const double inv_distance_COMPO = channel_data_.inv_distance_for_gradient(k_layer, idir, icompo);
  const double inv_distance_DIR = channel_data_.inv_distance_for_gradient(k_layer, icompo, idir);

  // Velocity field
  ConstIJK_double_ptr vCOMPO_ptr(get_v(_VCOMPO_), 0, 0, k_layer);
  ConstIJK_double_ptr vDIR_ptr(get_v(_DIR_), 0, 0, k_layer);
  ConstIJK_double_ptr molecular_nu(*molecular_nu_, 0, 0, k_layer);

  // Result (fluxes in direction x for component x of the convected field)
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);

  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double flux = 0.;

          if(_DIR_ == _VCOMPO_)
            {
              // recoding of Eval_Dift_VDF_var_Face::flux_fa7_elem
              Simd_double m_nu, m_nu_dummy;
              // flux is between left face and center face, hence, it is centered on the "left" cell
              molecular_nu.get_left_center(_DIR_, i, m_nu, m_nu_dummy);
              Simd_double v1, v2;
              vCOMPO_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              Simd_double tau = (v2 - v1);
              flux = m_nu * tau * surface;
            }
          else
            {
              // Interpolate diffusion coefficient from values at elements:
              Simd_double m_nu1, m_nu2, m_nu3, m_nu4;
              molecular_nu.get_left_center_c1c2(_DIR_, _VCOMPO_, i, m_nu1, m_nu2, m_nu3, m_nu4);
              // recoding of Eval_Dift_VDF_var_Face::flux_arete_interne
              Simd_double m_nu = (m_nu1 + m_nu2 + m_nu3 + m_nu4) * 0.25;
              // gradient in direction DIR of component COMPO
              Simd_double v3, v4;
              vCOMPO_ptr.get_left_center(_DIR_, i, v3, v4);
              Simd_double tau = (v4 - v3);
              flux = m_nu * tau * surface;
            }

          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      // instructions to perform to jump to next row
      vCOMPO_ptr.next_j();
      if(_DIR_ != _VCOMPO_)
        vDIR_ptr.next_j();
      molecular_nu.next_j();
      resu_ptr.next_j();
    }

}


/*! @brief compute fluxes in direction x for velocity component x for the layer of fluxes k_layer
 *
 *  Diffusion
 *
 */
void OpDiffStdWithLaminarTransposeAnisotropicIJK_double::flux_loop_(DIRECTION _DIR_, DIRECTION _VCOMPO_, IJK_Field_local_double& resu, int k_layer, int top_wall, int bottom_wall)
{
  // Velocity field
  ConstIJK_double_ptr vCOMPO_ptr(get_v(_VCOMPO_), 0, 0, k_layer);
  ConstIJK_double_ptr vDIR_ptr(get_v(_DIR_), 0, 0, k_layer);

  const int idir = (int)_DIR_;
  const int nx = _DIR_ == DIRECTION::X ? get_v(_VCOMPO_).ni() + 1 : get_v(_VCOMPO_).ni() ;
  const int ny = _DIR_ == DIRECTION::Y ? get_v(_VCOMPO_).nj() + 1 : get_v(_VCOMPO_).nj();
  const int icompo = (int)_VCOMPO_;

  // Result (fluxes in direction x for component x of the convected field)
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);

  const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
  // global index of the layer of flux of the wall
  //  (assume one walls at zmin and zmax)
  const int first_global_k_layer = channel_data_.first_global_k_layer_flux(icompo, idir);
  const int last_global_k_layer = channel_data_.last_global_k_layer_flux(icompo, idir);

  const double surface = - channel_data_.get_surface(k_layer, icompo, idir);
  const double inv_distance_COMPO = channel_data_.inv_distance_for_gradient(k_layer, idir, icompo);
  const double inv_distance_DIR = channel_data_.inv_distance_for_gradient(k_layer, icompo, idir);

  if(!perio_k_)
    {
      if (global_k_layer < first_global_k_layer || global_k_layer > last_global_k_layer)
        {
          // We are inside the wall
          putzero(resu);
          return;
        }
      // Variable diffusion coefficients: molecular and turbulent, at elements
    }
  ConstIJK_double_ptr molecular_nu(*molecular_nu_, 0, 0, k_layer);

  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double flux = 0.;
          if(_DIR_ == _VCOMPO_)
            {
              // recoding of Eval_Dift_VDF_var_Face::flux_fa7_elem
              Simd_double m_nu, m_nu_dummy;
              // flux is between left face and center face, hence, it is centered on the "left" cell
              molecular_nu.get_left_center(_DIR_, i, m_nu, m_nu_dummy);
              Simd_double v1, v2;
              vCOMPO_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              Simd_double tau = (v2 - v1);
              flux = m_nu * 2. * tau * surface; // The factor 2. is because tau_tr = tau
            }
          else
            {
              // Interpolate diffusion coefficient from values at elements:
              Simd_double m_nu1, m_nu2, m_nu3, m_nu4;
              molecular_nu.get_left_center_c1c2(_DIR_, _VCOMPO_, i, m_nu1, m_nu2, m_nu3, m_nu4);
              // recoding of Eval_Dift_VDF_var_Face::flux_arete_interne
              Simd_double m_nu = (m_nu1 + m_nu2 + m_nu3 + m_nu4) * 0.25;
              // gradient in direction DIR of component COMPO
              Simd_double v3, v4;
              vCOMPO_ptr.get_left_center(_DIR_, i, v3, v4);
              Simd_double tau = (v4 - v3);
              // gradient in direction COMPO of component DIR
              Simd_double v1, v2;
              vDIR_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              Simd_double tau_tr = (v2 - v1);
              flux = m_nu * (tau + tau_tr) * surface;
            }
          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      // instructions to perform to jump to next row
      vCOMPO_ptr.next_j();
      if(_DIR_ != _VCOMPO_)
        vDIR_ptr.next_j();
      molecular_nu.next_j();
      resu_ptr.next_j();
    }

}


/*! @brief compute fluxes in direction x for velocity component x for the layer of fluxes k_layer
 *
 *  Diffusion
 *
 */
void OpDiffStdWithLaminarTransposeAndDivergenceAnisotropicIJK_double::flux_loop_(DIRECTION _DIR_, DIRECTION _VCOMPO_, IJK_Field_local_double& resu, int k_layer, int top_wall, int bottom_wall)
{
  // Velocity field
  ConstIJK_double_ptr vCOMPO_ptr(get_v(_VCOMPO_), 0, 0, k_layer);
  ConstIJK_double_ptr vDIR_ptr(get_v(_DIR_), 0, 0, k_layer);

  const int idir = (int)_DIR_;
  const int nx = _DIR_ == DIRECTION::X ? get_v(_VCOMPO_).ni() + 1 : get_v(_VCOMPO_).ni() ;
  const int ny = _DIR_ == DIRECTION::Y ? get_v(_VCOMPO_).nj() + 1 : get_v(_VCOMPO_).nj();
  const int icompo = (int)_VCOMPO_;


  // Result (fluxes in direction x for component x of the convected field)
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);

  const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
  // global index of the layer of flux of the wall
  //  (assume one walls at zmin and zmax)
  const int first_global_k_layer = channel_data_.first_global_k_layer_flux(icompo, idir);
  const int last_global_k_layer = channel_data_.last_global_k_layer_flux(icompo, idir);

  const double surface = - channel_data_.get_surface(k_layer, icompo, idir);
  const double inv_distance_COMPO = channel_data_.inv_distance_for_gradient(k_layer, idir, icompo);
  const double inv_distance_DIR = channel_data_.inv_distance_for_gradient(k_layer, icompo, idir);

  if(!perio_k_)
    {
      if (global_k_layer < first_global_k_layer || global_k_layer > last_global_k_layer)
        {
          // We are inside the wall
          putzero(resu);
          return;
        }
      // Variable diffusion coefficients: molecular and turbulent, at elements
    }
  ConstIJK_double_ptr molecular_nu(*molecular_nu_, 0, 0, k_layer);
  ConstIJK_double_ptr div_ptr(*divergence_, 0, 0, k_layer);

  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {

          Simd_double flux = 0.;
          if(_DIR_ == _VCOMPO_)
            {
              // recoding of Eval_Dift_VDF_var_Face::flux_fa7_elem
              Simd_double m_nu, m_nu_dummy;
              // flux is between left face and center face, hence, it is centered on the "left" cell
              molecular_nu.get_left_center(_DIR_, i, m_nu, m_nu_dummy);
              Simd_double v1, v2;
              vCOMPO_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              Simd_double tau = (v2 - v1);
              Simd_double v5, v6_dummy;
              div_ptr.get_left_center(_DIR_, i, v5, v6_dummy);
              flux = m_nu * surface * (2. * tau - 0.66666666666666666 * v5); // The factor 2. is because tau_tr = tau
            }
          else
            {
              // Interpolate diffusion coefficient from values at elements:
              Simd_double m_nu1, m_nu2, m_nu3, m_nu4;
              molecular_nu.get_left_center_c1c2(_DIR_, _VCOMPO_, i, m_nu1, m_nu2, m_nu3, m_nu4);
              // recoding of Eval_Dift_VDF_var_Face::flux_arete_interne
              Simd_double m_nu = (m_nu1 + m_nu2 + m_nu3 + m_nu4) * 0.25;
              // gradient in direction DIR of component COMPO
              Simd_double v3, v4;
              vCOMPO_ptr.get_left_center(_DIR_, i, v3, v4);
              Simd_double tau = (v4 - v3);
              // gradient in direction COMPO of component DIR
              Simd_double v1, v2;
              vDIR_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              Simd_double tau_tr = (v2 - v1);
              flux = m_nu * (tau + tau_tr) * surface;

            }
          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      // instructions to perform to jump to next row
      vCOMPO_ptr.next_j();
      if(_DIR_ != _VCOMPO_)
        vDIR_ptr.next_j();
      molecular_nu.next_j();
      div_ptr.next_j();
      resu_ptr.next_j();
    }

}



/*! @brief compute fluxes in direction x for velocity component x for the layer of fluxes k_layer
 *
 *  Diffusion
 *
 */
void OpDiffTensorialZeroatwallIJK_double::flux_loop_(DIRECTION _DIR_, DIRECTION _VCOMPO_, IJK_Field_local_double& resu, int k_layer, int top_wall, int bottom_wall)
{
  // Velocity field
  ConstIJK_double_ptr vCOMPO_ptr(get_v(_VCOMPO_), 0, 0, k_layer);
  ConstIJK_double_ptr vDIR_ptr(get_v(_DIR_), 0, 0, k_layer);

  const int idir = (int)_DIR_;
  const int nx = _DIR_ == DIRECTION::X ? get_v(_VCOMPO_).ni() + 1 : get_v(_VCOMPO_).ni() ;
  const int ny = _DIR_ == DIRECTION::Y ? get_v(_VCOMPO_).nj() + 1 : get_v(_VCOMPO_).nj();
  const int icompo = (int)_VCOMPO_;

  // Result (fluxes in direction x for component x of the convected field)
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);

  const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
  // global index of the layer of flux of the wall
  //  (assume one walls at zmin and zmax)
  const int first_global_k_layer = channel_data_.first_global_k_layer_flux(icompo, idir);
  const int last_global_k_layer = channel_data_.last_global_k_layer_flux(icompo, idir);

  const double surface = - channel_data_.get_surface(k_layer, icompo, idir);
  const double inv_distance_COMPO = channel_data_.inv_distance_for_gradient(k_layer, idir, icompo);
  const double inv_distance_DIR = channel_data_.inv_distance_for_gradient(k_layer, icompo, idir);

  if(!perio_k_)
    {
      if (global_k_layer < first_global_k_layer || global_k_layer > last_global_k_layer)
        {
          // We are inside the wall
          putzero(resu);
          return;
        }
      // Variable diffusion coefficients: molecular and turbulent, at elements
    }
  ConstIJK_double_ptr molecular_nu(get_molecular_nu(_VCOMPO_, _DIR_), 0, 0, k_layer);

  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double flux = 0.;
          if(_DIR_ == _VCOMPO_)
            {
              // recoding of Eval_Dift_VDF_var_Face::flux_fa7_elem
              Simd_double m_nu, m_nu_dummy;
              // flux is between left face and center face, hence, it is centered on the "left" cell
              molecular_nu.get_left_center(_DIR_, i, m_nu, m_nu_dummy);
              Simd_double v1, v2;
              vCOMPO_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              Simd_double tau = (v2 - v1) * inv_distance_DIR;
              flux = m_nu * tau * surface;
            }
          else
            {
              // Interpolate diffusion coefficient from values at elements:
              Simd_double m_nu1, m_nu2, m_nu3, m_nu4;
              molecular_nu.get_left_center_c1c2(_DIR_, _VCOMPO_, i, m_nu1, m_nu2, m_nu3, m_nu4);
              // recoding of Eval_Dift_VDF_var_Face::flux_arete_interne
              Simd_double m_nu = (m_nu1 + m_nu2 + m_nu3 + m_nu4) * 0.25;
              // gradient in direction DIR of component COMPO
              Simd_double v3, v4;
              vCOMPO_ptr.get_left_center(_DIR_, i, v3, v4);
              Simd_double tau = (v4 - v3) * inv_distance_DIR;
              flux = m_nu * tau * surface;

            }
          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      // instructions to perform to jump to next row
      vCOMPO_ptr.next_j();
      if(_DIR_ != _VCOMPO_)
        vDIR_ptr.next_j();
      molecular_nu.next_j();
      resu_ptr.next_j();
    }

}


/*! @brief compute fluxes in direction x for velocity component x for the layer of fluxes k_layer
 *
 *  Diffusion
 *
 */
void OpDiffStdWithLaminarTransposeTensorialZeroatwallIJK_double::flux_loop_(DIRECTION _DIR_, DIRECTION _VCOMPO_, IJK_Field_local_double& resu, int k_layer, int top_wall, int bottom_wall)
{
  // Velocity field
  ConstIJK_double_ptr vCOMPO_ptr(get_v(_VCOMPO_), 0, 0, k_layer);
  ConstIJK_double_ptr vDIR_ptr(get_v(_DIR_), 0, 0, k_layer);

  const int idir = (int)_DIR_;
  const int nx = _DIR_ == DIRECTION::X ? get_v(_VCOMPO_).ni() + 1 : get_v(_VCOMPO_).ni() ;
  const int ny = _DIR_ == DIRECTION::Y ? get_v(_VCOMPO_).nj() + 1 : get_v(_VCOMPO_).nj();
  const int icompo = (int)_VCOMPO_;


  // Result (fluxes in direction x for component x of the convected field)
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);

  const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
  // global index of the layer of flux of the wall
  //  (assume one walls at zmin and zmax)
  const int first_global_k_layer = channel_data_.first_global_k_layer_flux(icompo, idir);
  const int last_global_k_layer = channel_data_.last_global_k_layer_flux(icompo, idir);

  const double surface = - channel_data_.get_surface(k_layer, icompo, idir);
  const double inv_distance_COMPO = channel_data_.inv_distance_for_gradient(k_layer, idir, icompo);
  const double inv_distance_DIR = channel_data_.inv_distance_for_gradient(k_layer, icompo, idir);

  if(!perio_k_)
    {
      if (global_k_layer < first_global_k_layer || global_k_layer > last_global_k_layer)
        {
          // We are inside the wall
          putzero(resu);
          return;
        }
      // Variable diffusion coefficients: molecular and turbulent, at elements
    }
  ConstIJK_double_ptr molecular_nu(get_molecular_nu(_VCOMPO_, _DIR_), 0, 0, k_layer);

  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double flux = 0.;
          if(_VCOMPO_ == _DIR_)
            {
              // recoding of Eval_Dift_VDF_var_Face::flux_fa7_elem
              Simd_double m_nu, m_nu_dummy;
              // flux is between left face and center face, hence, it is centered on the "left" cell
              molecular_nu.get_left_center(_DIR_, i, m_nu, m_nu_dummy);
              Simd_double v1, v2;
              vCOMPO_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              Simd_double tau = (v2 - v1) * inv_distance_DIR;
              flux = m_nu * 2. * tau * surface; // The factor 2. is because tau_tr = tau
            }
          else
            {
              // Interpolate diffusion coefficient from values at elements:
              Simd_double m_nu1, m_nu2, m_nu3, m_nu4;
              molecular_nu.get_left_center_c1c2(_DIR_, _VCOMPO_, i, m_nu1, m_nu2, m_nu3, m_nu4);
              // recoding of Eval_Dift_VDF_var_Face::flux_arete_interne
              Simd_double m_nu = (m_nu1 + m_nu2 + m_nu3 + m_nu4) * 0.25;
              // gradient in direction DIR of component COMPO
              Simd_double v3, v4;
              vCOMPO_ptr.get_left_center(_DIR_, i, v3, v4);
              Simd_double tau = (v4 - v3) * inv_distance_DIR;
              // gradient in direction COMPO of component DIR
              Simd_double v1, v2;
              vDIR_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              Simd_double tau_tr = (v2 - v1) * inv_distance_COMPO;
              flux = m_nu * (tau + tau_tr) * surface;

            }
          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      // instructions to perform to jump to next row
      vCOMPO_ptr.next_j();
      if(_DIR_ != _VCOMPO_)
        vDIR_ptr.next_j();
      molecular_nu.next_j();
      resu_ptr.next_j();
    }

}

/*! @brief compute fluxes in direction x for velocity component x for the layer of fluxes k_layer
 *
 *  Diffusion
 *
 */
void OpDiffStdWithLaminarTransposeAndDivergenceTensorialZeroatwallIJK_double::flux_loop_(DIRECTION _DIR_, DIRECTION _VCOMPO_, IJK_Field_local_double& resu, int k_layer, int top_wall, int bottom_wall)
{
  // Velocity field
  ConstIJK_double_ptr vCOMPO_ptr(get_v(_VCOMPO_), 0, 0, k_layer);
  ConstIJK_double_ptr vDIR_ptr(get_v(_DIR_), 0, 0, k_layer);

  const int idir = (int)_DIR_;
  const int nx = _DIR_ == DIRECTION::X ? get_v(_VCOMPO_).ni() + 1 : get_v(_VCOMPO_).ni() ;
  const int ny = _DIR_ == DIRECTION::Y ? get_v(_VCOMPO_).nj() + 1 : get_v(_VCOMPO_).nj();
  const int icompo = (int)_VCOMPO_;


  // Result (fluxes in direction x for component x of the convected field)
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);

  const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
  // global index of the layer of flux of the wall
  //  (assume one walls at zmin and zmax)
  const int first_global_k_layer = channel_data_.first_global_k_layer_flux(icompo, idir);
  const int last_global_k_layer = channel_data_.last_global_k_layer_flux(icompo, idir);

  const double surface = - channel_data_.get_surface(k_layer, icompo, idir);
  const double inv_distance_COMPO = channel_data_.inv_distance_for_gradient(k_layer, idir, icompo);
  const double inv_distance_DIR = channel_data_.inv_distance_for_gradient(k_layer, icompo, idir);

  if(!perio_k_)
    {
      if (global_k_layer < first_global_k_layer || global_k_layer > last_global_k_layer)
        {
          // We are inside the wall
          putzero(resu);
          return;
        }
      // Variable diffusion coefficients: molecular and turbulent, at elements
    }
  ConstIJK_double_ptr molecular_nu(get_molecular_nu(_VCOMPO_, _DIR_), 0, 0, k_layer);
  ConstIJK_double_ptr div_ptr(*divergence_, 0, 0, k_layer);

  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double flux = 0.;
          if(_DIR_ == _VCOMPO_)
            {
              // recoding of Eval_Dift_VDF_var_Face::flux_fa7_elem
              Simd_double m_nu, m_nu_dummy;
              // flux is between left face and center face, hence, it is centered on the "left" cell
              molecular_nu.get_left_center(_DIR_, i, m_nu, m_nu_dummy);
              Simd_double v1, v2;
              vCOMPO_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              Simd_double tau = (v2 - v1) * inv_distance_DIR;
              Simd_double v5, v6_dummy;
              div_ptr.get_left_center(_DIR_, i, v5, v6_dummy);
              flux = m_nu * surface * (2. * tau - 0.66666666666666666 * v5); // The factor 2. is because tau_tr = tau
            }
          else
            {
              // Interpolate diffusion coefficient from values at elements:
              Simd_double m_nu1, m_nu2, m_nu3, m_nu4;
              molecular_nu.get_left_center_c1c2(_DIR_, _VCOMPO_, i, m_nu1, m_nu2, m_nu3, m_nu4);
              // recoding of Eval_Dift_VDF_var_Face::flux_arete_interne
              Simd_double m_nu = (m_nu1 + m_nu2 + m_nu3 + m_nu4) * 0.25;
              // gradient in direction x of component y
              Simd_double v3, v4;
              vCOMPO_ptr.get_left_center(_DIR_, i, v3, v4);
              Simd_double tau = (v4 - v3) * inv_distance_DIR;
              // gradient in direction y of component x
              Simd_double v1, v2;
              vDIR_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              Simd_double tau_tr = (v2 - v1) * inv_distance_COMPO;
              flux = m_nu * (tau + tau_tr) * surface;
            }
          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      // instructions to perform to jump to next row
      vCOMPO_ptr.next_j();
      if(_DIR_ != _VCOMPO_)
        vDIR_ptr.next_j();
      molecular_nu.next_j();
      div_ptr.next_j();
      resu_ptr.next_j();
    }
}


/*! @brief compute fluxes in direction x for velocity component x for the layer of fluxes k_layer
 *
 *  Diffusion
 *
 */
void OpDiffTensorialAnisotropicZeroatwallIJK_double::flux_loop_(DIRECTION _DIR_, DIRECTION _VCOMPO_, IJK_Field_local_double& resu, int k_layer, int top_wall, int bottom_wall)
{
  // Velocity field
  ConstIJK_double_ptr vCOMPO_ptr(get_v(_VCOMPO_), 0, 0, k_layer);
  ConstIJK_double_ptr vDIR_ptr(get_v(_DIR_), 0, 0, k_layer);

  const int idir = (int)_DIR_;
  const int nx = _DIR_ == DIRECTION::X ? get_v(_VCOMPO_).ni() + 1 : get_v(_VCOMPO_).ni() ;
  const int ny = _DIR_ == DIRECTION::Y ? get_v(_VCOMPO_).nj() + 1 : get_v(_VCOMPO_).nj();
  const int icompo = (int)_VCOMPO_;


  // Result (fluxes in direction x for component x of the convected field)
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);

  const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
  // global index of the layer of flux of the wall
  //  (assume one walls at zmin and zmax)
  const int first_global_k_layer = channel_data_.first_global_k_layer_flux(icompo, idir);
  const int last_global_k_layer = channel_data_.last_global_k_layer_flux(icompo, idir);

  const double surface = - channel_data_.get_surface(k_layer, icompo, idir);
  const double inv_distance_COMPO = channel_data_.inv_distance_for_gradient(k_layer, idir, icompo);
  const double inv_distance_DIR = channel_data_.inv_distance_for_gradient(k_layer, icompo, idir);

  if(!perio_k_)
    {
      if (global_k_layer < first_global_k_layer || global_k_layer > last_global_k_layer)
        {
          // We are inside the wall
          putzero(resu);
          return;
        }
      // Variable diffusion coefficients: molecular and turbulent, at elements
    }
  ConstIJK_double_ptr molecular_nu(get_molecular_nu(_VCOMPO_, _DIR_), 0, 0, k_layer);

  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double flux = 0.;

          if(_DIR_ == _VCOMPO_)
            {
              // recoding of Eval_Dift_VDF_var_Face::flux_fa7_elem
              Simd_double m_nu, m_nu_dummy;
              // flux is between left face and center face, hence, it is centered on the "left" cell
              molecular_nu.get_left_center(_DIR_, i, m_nu, m_nu_dummy);
              Simd_double v1, v2;
              vCOMPO_ptr.get_left_center(_VCOMPO_,i, v1, v2);
              Simd_double tau = (v2 - v1);
              flux = m_nu * tau * surface;
            }
          else
            {
              // Interpolate diffusion coefficient from values at elements:
              Simd_double m_nu1, m_nu2, m_nu3, m_nu4;
              molecular_nu.get_left_center_c1c2(_DIR_, _VCOMPO_, i, m_nu1, m_nu2, m_nu3, m_nu4);
              // recoding of Eval_Dift_VDF_var_Face::flux_arete_interne
              Simd_double m_nu = (m_nu1 + m_nu2 + m_nu3 + m_nu4) * 0.25;
              // gradient in direction DIR of component COMPO
              Simd_double v3, v4;
              vCOMPO_ptr.get_left_center(_DIR_, i, v3, v4);
              Simd_double tau = (v4 - v3);
              flux = m_nu * tau * surface;
            }
          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      // instructions to perform to jump to next row
      vCOMPO_ptr.next_j();
      if(_DIR_ != _VCOMPO_)
        vDIR_ptr.next_j();
      molecular_nu.next_j();
      resu_ptr.next_j();
    }
}


/*! @brief compute fluxes in direction x for velocity component x for the layer of fluxes k_layer
 *
 *  Diffusion
 *
 */
void OpDiffStdWithLaminarTransposeTensorialAnisotropicZeroatwallIJK_double::flux_loop_(DIRECTION _DIR_, DIRECTION _VCOMPO_, IJK_Field_local_double& resu, int k_layer, int top_wall, int bottom_wall)
{
  // Velocity field
  ConstIJK_double_ptr vCOMPO_ptr(get_v(_VCOMPO_), 0, 0, k_layer);
  ConstIJK_double_ptr vDIR_ptr(get_v(_DIR_), 0, 0, k_layer);

  const int idir = (int)_DIR_;
  const int nx = _DIR_ == DIRECTION::X ? get_v(_VCOMPO_).ni() + 1 : get_v(_VCOMPO_).ni() ;
  const int ny = _DIR_ == DIRECTION::Y ? get_v(_VCOMPO_).nj() + 1 : get_v(_VCOMPO_).nj();
  const int icompo = (int)_VCOMPO_;


  // Result (fluxes in direction x for component x of the convected field)
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);

  const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
  // global index of the layer of flux of the wall
  //  (assume one walls at zmin and zmax)
  const int first_global_k_layer = channel_data_.first_global_k_layer_flux(icompo, idir);
  const int last_global_k_layer = channel_data_.last_global_k_layer_flux(icompo, idir);

  const double surface = - channel_data_.get_surface(k_layer, icompo, idir);
  const double inv_distance_COMPO = channel_data_.inv_distance_for_gradient(k_layer, idir, icompo);
  const double inv_distance_DIR = channel_data_.inv_distance_for_gradient(k_layer, icompo, idir);

  if(!perio_k_)
    {
      if (global_k_layer < first_global_k_layer || global_k_layer > last_global_k_layer)
        {
          // We are inside the wall
          putzero(resu);
          return;
        }
      // Variable diffusion coefficients: molecular and turbulent, at elements
    }
  ConstIJK_double_ptr molecular_nu(get_molecular_nu(_VCOMPO_, _DIR_), 0, 0, k_layer);

  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double flux = 0.;
          if(_DIR_ == _VCOMPO_)
            {
              // recoding of Eval_Dift_VDF_var_Face::flux_fa7_elem
              Simd_double m_nu, m_nu_dummy;
              // flux is between left face and center face, hence, it is centered on the "left" cell
              molecular_nu.get_left_center(_DIR_,i, m_nu, m_nu_dummy);
              Simd_double v1, v2;
              vCOMPO_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              Simd_double tau = (v2 - v1);
              flux = m_nu * 2. * tau * surface; // The factor 2. is because tau_tr = tau
            }
          else
            {
              // Interpolate diffusion coefficient from values at elements:
              Simd_double m_nu1, m_nu2, m_nu3, m_nu4;
              molecular_nu.get_left_center_c1c2(_DIR_, _VCOMPO_, i, m_nu1, m_nu2, m_nu3, m_nu4);
              // recoding of Eval_Dift_VDF_var_Face::flux_arete_interne
              Simd_double m_nu = (m_nu1 + m_nu2 + m_nu3 + m_nu4) * 0.25;
              // gradient in direction x of component y
              Simd_double v3, v4;
              vCOMPO_ptr.get_left_center(_DIR_, i, v3, v4);
              Simd_double tau = (v4 - v3);
              // gradient in direction y of component x
              Simd_double v1, v2;
              vDIR_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              Simd_double tau_tr = (v2 - v1);
              flux = m_nu * (tau + tau_tr) * surface;

            }
          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      // instructions to perform to jump to next row
      vCOMPO_ptr.next_j();
      if(_DIR_ != _VCOMPO_)
        vDIR_ptr.next_j();
      molecular_nu.next_j();
      resu_ptr.next_j();
    }

}

/*! @brief compute fluxes in direction x for velocity component x for the layer of fluxes k_layer
 *
 *  Diffusion
 *
 */
void OpDiffStdWithLaminarTransposeAndDivergenceTensorialAnisotropicZeroatwallIJK_double::flux_loop_(DIRECTION _DIR_, DIRECTION _VCOMPO_, IJK_Field_local_double& resu, int k_layer, int top_wall, int bottom_wall)
{
  // Velocity field
  ConstIJK_double_ptr vCOMPO_ptr(get_v(_VCOMPO_), 0, 0, k_layer);
  ConstIJK_double_ptr vDIR_ptr(get_v(_DIR_), 0, 0, k_layer);

  const int idir = (int)_DIR_;
  const int nx = _DIR_ == DIRECTION::X ? get_v(_VCOMPO_).ni() + 1 : get_v(_VCOMPO_).ni() ;
  const int ny = _DIR_ == DIRECTION::Y ? get_v(_VCOMPO_).nj() + 1 : get_v(_VCOMPO_).nj();
  const int icompo = (int)_VCOMPO_;


  // Result (fluxes in direction x for component x of the convected field)
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);

  const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
  // global index of the layer of flux of the wall
  //  (assume one walls at zmin and zmax)
  const int first_global_k_layer = channel_data_.first_global_k_layer_flux(icompo, idir);
  const int last_global_k_layer = channel_data_.last_global_k_layer_flux(icompo, idir);

  const double surface = - channel_data_.get_surface(k_layer, icompo, idir);
  const double inv_distance_COMPO = channel_data_.inv_distance_for_gradient(k_layer, idir, icompo);
  const double inv_distance_DIR = channel_data_.inv_distance_for_gradient(k_layer, icompo, idir);

  if(!perio_k_)
    {
      if (global_k_layer < first_global_k_layer || global_k_layer > last_global_k_layer)
        {
          // We are inside the wall
          putzero(resu);
          return;
        }
      // Variable diffusion coefficients: molecular and turbulent, at elements
    }
  ConstIJK_double_ptr molecular_nu(get_molecular_nu(_VCOMPO_, _DIR_), 0, 0, k_layer);
  ConstIJK_double_ptr div_ptr(*divergence_, 0, 0, k_layer);

  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double flux = 0.;
          if(_DIR_ == _VCOMPO_)
            {
              // recoding of Eval_Dift_VDF_var_Face::flux_fa7_elem
              Simd_double m_nu, m_nu_dummy;
              // flux is between left face and center face, hence, it is centered on the "left" cell
              molecular_nu.get_left_center(_DIR_, i, m_nu, m_nu_dummy);
              Simd_double v1, v2;
              vCOMPO_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              Simd_double tau = (v2 - v1);
              Simd_double v5, v6_dummy;
              div_ptr.get_left_center(_DIR_, i, v5, v6_dummy);
              flux = m_nu * surface * (2. * tau - 0.66666666666666666 * v5); // The factor 2. is because tau_tr = tau
            }
          else
            {
              // Interpolate diffusion coefficient from values at elements:
              Simd_double m_nu1, m_nu2, m_nu3, m_nu4;
              molecular_nu.get_left_center_c1c2(_DIR_, _VCOMPO_, i, m_nu1, m_nu2, m_nu3, m_nu4);
              // recoding of Eval_Dift_VDF_var_Face::flux_arete_interne
              Simd_double m_nu = (m_nu1 + m_nu2 + m_nu3 + m_nu4) * 0.25;
              // gradient in direction x of component y
              Simd_double v3, v4;
              vCOMPO_ptr.get_left_center(_DIR_, i, v3, v4);
              Simd_double tau = (v4 - v3);
              // gradient in direction y of component x
              Simd_double v1, v2;
              vDIR_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              Simd_double tau_tr = (v2 - v1);
              flux = m_nu * (tau + tau_tr) * surface;

            }
          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      // instructions to perform to jump to next row
      vCOMPO_ptr.next_j();
      if(_DIR_ != _VCOMPO_)
        vDIR_ptr.next_j();
      molecular_nu.next_j();
      div_ptr.next_j();
      resu_ptr.next_j();
    }

}

/*! @brief compute fluxes in direction x for velocity component x for the layer of fluxes k_layer
 *
 *  Diffusion
 *
 */
void OpDiffStructuralOnlyZeroatwallIJK_double::flux_loop_(DIRECTION _DIR_, DIRECTION _VCOMPO_, IJK_Field_local_double& resu, int k_layer, int top_wall, int bottom_wall)
{
  // Velocity field
  ConstIJK_double_ptr vCOMPO_ptr(get_v(_VCOMPO_), 0, 0, k_layer);
  ConstIJK_double_ptr vDIR_ptr(get_v(_DIR_), 0, 0, k_layer);

  const int idir = (int)_DIR_;
  const int nx = _DIR_ == DIRECTION::X ? get_v(_VCOMPO_).ni() + 1 : get_v(_VCOMPO_).ni() ;
  const int ny = _DIR_ == DIRECTION::Y ? get_v(_VCOMPO_).nj() + 1 : get_v(_VCOMPO_).nj();
  const int icompo = (int)_VCOMPO_;


  // Result (fluxes in direction x for component x of the convected field)
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);

  const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
  // global index of the layer of flux of the wall
  //  (assume one walls at zmin and zmax)
  const int first_global_k_layer = channel_data_.first_global_k_layer_flux(icompo, idir);
  const int last_global_k_layer = channel_data_.last_global_k_layer_flux(icompo, idir);

  const double surface = - channel_data_.get_surface(k_layer, icompo, idir);
  const double inv_distance_COMPO = channel_data_.inv_distance_for_gradient(k_layer, idir, icompo);
  const double inv_distance_DIR = channel_data_.inv_distance_for_gradient(k_layer, icompo, idir);

  if(!perio_k_)
    {
      if (global_k_layer < first_global_k_layer || global_k_layer > last_global_k_layer)
        {
          // We are inside the wall
          putzero(resu);
          return;
        }
      // Variable diffusion coefficients: molecular and turbulent, at elements
    }
  ConstIJK_double_ptr structural_model(get_structural_model(_DIR_, _VCOMPO_), 0, 0, k_layer);

  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double flux = 0.;
          if(_DIR_ == _VCOMPO_)
            {
              Simd_double s_mo, s_mo_dummy;
              structural_model.get_left_center(_DIR_, i, s_mo, s_mo_dummy);
              flux = s_mo * surface;
            }
          else
            {
              // Interpolate diffusion coefficient from values at elements:
              Simd_double s_mo, s_mo_dummy;
              structural_model.get_left_center(_DIR_, i, s_mo_dummy, s_mo);
              // recoding of Eval_Dift_VDF_var_Face::flux_arete_interne
              // gradient in direction x of component y
              Simd_double v3, v4;
              vCOMPO_ptr.get_left_center(_DIR_, i, v3, v4);
              flux = s_mo * surface;

            }
          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      // instructions to perform to jump to next row
      vCOMPO_ptr.next_j();
      if(_DIR_ != _VCOMPO_)
        vDIR_ptr.next_j();
      structural_model.next_j();
      resu_ptr.next_j();
    }
}
 */
