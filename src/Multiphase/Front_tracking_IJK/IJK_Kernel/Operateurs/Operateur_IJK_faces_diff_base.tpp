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

#ifndef Operateur_IJK_faces_diff_base_TPP_included
#define Operateur_IJK_faces_diff_base_TPP_included
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

/*! @brief compute fluxes in direction DIR for velocity component COMPO for the layer of fluxes k_layer
 *
 *  Diffusion
 *
 */
template<DIRECTION _DIR_, DIRECTION _VCOMPO_>
void Operateur_IJK_faces_diff_base_double::compute_flux_(IJK_Field_local_double& resu, const int k_layer)
{
  const int idir = (int)_DIR_;
  const int icompo = (int)_VCOMPO_;

  const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
  // global index of the layer of flux of the wall
  //  (assume one walls at zmin and zmax)
  const int first_global_k_layer = channel_data_.first_global_k_layer_flux(icompo, idir);
  const int last_global_k_layer = channel_data_.last_global_k_layer_flux(icompo, idir);

  Boundary_Conditions::BCType bc_type = ref_bc_.valeur().get_bctype_k_min();

  // for mixte_shear boundary condition >> need to localize to and bottom boundary even if perio_k
  if(!perio_k_ || bc_type == Boundary_Conditions::Mixte_shear)
    {
      // For this direction and this component, we possibly have a wall boundary condition to treat:
      if(_DIR_ == DIRECTION::Z && _VCOMPO_ != DIRECTION::Z)
        {
          int top_wall = 0, bottom_wall = 0;
          // We are at bottom wall: z=0
          if (global_k_layer == first_global_k_layer-1)
            {
              bc_type = ref_bc_.valeur().get_bctype_k_min();
              bottom_wall = 1;
            }
          if (global_k_layer == last_global_k_layer + 1)
            {
              bc_type = ref_bc_.valeur().get_bctype_k_max();
              top_wall = 1;
            }

          if(bottom_wall || top_wall)
            {
              switch(bc_type)
                {
                case Boundary_Conditions::Paroi:
                  {
                    flux_loop_<_DIR_, _VCOMPO_>(resu, k_layer, top_wall, bottom_wall);
                    break;
                  }
                case Boundary_Conditions::Mixte_shear:
                  {
                    flux_loop_<_DIR_, _VCOMPO_>(resu, k_layer, top_wall, bottom_wall);
                    break;
                  }
                case Boundary_Conditions::Symetrie:
                  {
                    // Symetry boundary condition: momentum flux is zero
                    putzero(resu);
                    break;
                  }
                default:
                  Cerr << "Operateur_IJK_faces_diff_base_double::compute_flux_ wrong boundary condition." << finl;
                  Process::exit();
                }
              return;
            }
        }

      if (global_k_layer < first_global_k_layer || global_k_layer > last_global_k_layer)
        {
          // We are inside the wall
          putzero(resu);
          return;
        }
    }

  flux_loop_<_DIR_, _VCOMPO_>(resu, k_layer);
}


template<DIRECTION _DIR_, DIRECTION _VCOMPO_>
void Operateur_IJK_faces_diff_base_double::flux_loop_(IJK_Field_local_double& resu, int k_layer, int top_wall, int bottom_wall )
{
  const IJK_Field_local_double& vCOMPO = get_v(_VCOMPO_);
  const IJK_Field_local_double& vDIR = get_v(_DIR_);

  const int idir = (int)_DIR_;
  const int nx = _DIR_ == DIRECTION::X ? vCOMPO.ni() + 1 : vCOMPO.ni() ;
  const int ny = _DIR_ == DIRECTION::Y ? vCOMPO.nj() + 1 : vCOMPO.nj();
  const int icompo = (int)_VCOMPO_;

  const double surface = - channel_data_.get_surface(k_layer, icompo, idir);
  const double inv_distance_COMPO = channel_data_.inv_distance_for_gradient(k_layer, idir, icompo);
  const double inv_distance_DIR = channel_data_.inv_distance_for_gradient(k_layer, icompo, idir);

  // Velocity field
  ConstIJK_double_ptr vCOMPO_ptr(vCOMPO, 0, 0, k_layer);
  ConstIJK_double_ptr vDIR_ptr(vDIR, 0, 0, k_layer);

  const IJK_Field_local_double& dummy_field = vCOMPO;

  ConstIJK_double_ptr molecular_nu(!is_structural_ ? (is_tensorial_? get_coeff_tensor(_VCOMPO_, _DIR_) : get_nu()) : dummy_field , 0, 0, k_layer);

  ConstIJK_double_ptr div_ptr(with_divergence_ ? get_divergence() : dummy_field, 0, 0, k_layer);

  // Turbulent kinetic energy from turbulence model
  //  ConstIJK_double_ptr turbulent_nu(is_turb_ ? get_turbulent_nu() : dummy_field, 0, 0, k_layer);
  //  ConstIJK_double_ptr turbulent_k_energy(is_turb_ ? get_turbulent_k_energy() : dummy_field, 0, 0, k_layer );
  //  ConstIJK_double_ptr structural_model(is_structural_ ? get_structural_model(_DIR_, _VCOMPO_) : dummy_field, 0, 0, k_layer);
  //  ConstIJK_double_ptr structural_model(is_structural_ ? get_coeff_tensor(_DIR_, _VCOMPO_) : dummy_field, 0, 0, k_layer);

  ConstIJK_double_ptr turbulent_nu(is_turb_ ? get_nu() : *nu_, 0, 0, k_layer);
  ConstIJK_double_ptr turbulent_k_energy(is_turb_ ? get_nu() : *nu_, 0, 0, k_layer );
  ConstIJK_double_ptr structural_model(is_structural_ ? get_coeff_tensor(_DIR_, _VCOMPO_) : *nu_, 0, 0, k_layer);

  // Result (fluxes in direction DIR for component COMPO of the convected field)
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
              flux_loop_same_dir_compo_<_DIR_, _VCOMPO_>(i, surface, inv_distance_DIR, inv_distance_COMPO,
                                                         vCOMPO_ptr, vDIR_ptr,
                                                         molecular_nu, div_ptr,
                                                         turbulent_nu, turbulent_k_energy,
                                                         structural_model, flux);
            }
          else
            {
              flux_loop_different_dir_compo_<_DIR_, _VCOMPO_>(i, surface, inv_distance_DIR, inv_distance_COMPO,
                                                              top_wall, bottom_wall,
                                                              vCOMPO_ptr, vDIR_ptr,
                                                              molecular_nu, div_ptr,
                                                              turbulent_nu, turbulent_k_energy,
                                                              structural_model, flux);

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
      if(!is_structural_)
        molecular_nu.next_j();
      if(is_turb_)
        {
          turbulent_nu.next_j();
          turbulent_k_energy.next_j();
        }
      if(with_divergence_)
        div_ptr.next_j();
      if(is_structural_)
        structural_model.next_j();

      resu_ptr.next_j();

    }
}

template<DIRECTION _DIR_, DIRECTION _VCOMPO_>
void Operateur_IJK_faces_diff_base_double::flux_loop_same_dir_compo_(int i, double surface, double inv_distance_DIR, double inv_distance_COMPO,
                                                                     const ConstIJK_double_ptr& vCOMPO_ptr, const ConstIJK_double_ptr& vDIR_ptr,
                                                                     const ConstIJK_double_ptr& molecular_nu, const ConstIJK_double_ptr& div_ptr,
                                                                     const ConstIJK_double_ptr& turbulent_nu, const ConstIJK_double_ptr& turbulent_k_energy,
                                                                     const ConstIJK_double_ptr& structural_model, Simd_double& flux )
{
  if(is_structural_)
    {
      Simd_double s_mo, s_mo_dummy;
      structural_model.get_left_center(_DIR_, i, s_mo, s_mo_dummy);
      flux = s_mo * surface;
    }
  else
    {
      // recoding of Eval_Dift_VDF_var_Face::flux_fa7_elem
      Simd_double m_nu, m_nu_dummy;
      // flux is between left face and center face, hence, it is centered on the "left" cell
      molecular_nu.get_left_center(_DIR_, i, m_nu, m_nu_dummy);
      Simd_double v1, v2;
      vCOMPO_ptr.get_left_center(_VCOMPO_, i, v1, v2);
      Simd_double tau = (v2 - v1);
      if(!is_anisotropic_) tau *= inv_distance_COMPO;

      Simd_double minus_reyn = 0.;
      if(is_turb_)
        {
          Simd_double t_nu, t_nu_dummy;
          turbulent_nu.get_left_center(_DIR_, i, t_nu, t_nu_dummy);
          minus_reyn = t_nu * tau * (2.);
          // Shall we include "k" or not ? If we do, it is equivalent to adding
          //  grad(k) to the rhs of Navier Stokes, il we be "eaten" by the pressure
          //  solver.
          Simd_double k, k_dummy;
          // flux is between left face and center face, hence, it is centered on the "left" cell
          turbulent_k_energy.get_left_center(_DIR_, i, k, k_dummy);
          minus_reyn += (-0.66666666666666666) * k;
          // Check that diagonal terms of the Reynolds tensor are positive
          //  (exactly copied from Eval_Dift_VDF_var_Face::flux_fa7_elem)
          Simd_double zeroVec = 0.;
          minus_reyn = SimdMin(minus_reyn, zeroVec);
        }

      flux = (minus_reyn + m_nu * tau) * surface;

      if(with_transpose_)
        flux *= 2.; // The factor 2. is because tau_tr = tau

      if(with_divergence_)
        {
          Simd_double v5, v6_dummy;
          div_ptr.get_left_center(_DIR_,i, v5, v6_dummy);
          flux -= m_nu * surface * 0.66666666666666666 * v5; // The factor 2. is because tau_tr = tau
        }

    }

}

template<DIRECTION _DIR_, DIRECTION _VCOMPO_>
void Operateur_IJK_faces_diff_base_double::flux_loop_different_dir_compo_(int i, double surface, double inv_distance_DIR, double inv_distance_COMPO, int top_wall, int bottom_wall,
                                                                          const ConstIJK_double_ptr& vCOMPO_ptr, const ConstIJK_double_ptr& vDIR_ptr,
                                                                          const ConstIJK_double_ptr& molecular_nu, const ConstIJK_double_ptr& div_ptr,
                                                                          const ConstIJK_double_ptr& turbulent_nu, const ConstIJK_double_ptr& turbulent_k_energy,
                                                                          const ConstIJK_double_ptr& structural_model, Simd_double& flux )
{

  if(is_structural_)
    {
      Simd_double s_mo, s_mo_dummy;
      structural_model.get_left_center(_DIR_, i, s_mo_dummy, s_mo);
      // recoding of Eval_Dift_VDF_var_Face::flux_arete_interne
      // gradient in direction x of component y
      Simd_double v3, v4;
      vCOMPO_ptr.get_left_center(_DIR_, i, v3, v4);
      flux = s_mo * surface;

    }
  else
    {

      Boundary_Conditions::BCType bc_type = ref_bc_.valeur().get_bctype_k_min();
      // Interpolate diffusion coefficient from values at elements:
      Simd_double m_nu1, m_nu2, m_nu3, m_nu4;
      molecular_nu.get_left_center_c1c2(_DIR_, _VCOMPO_, i, m_nu1, m_nu2, m_nu3, m_nu4);
      double mult_coeff = 0.25;

      // for wall boundary conditions
      if(bottom_wall && bc_type!=Boundary_Conditions::Mixte_shear)
        {
          // bottom wall (z=0)
          // nu1 and nu2 are "left" in direction z, hence in the wall:
          m_nu1 = 0., m_nu2 = 0.;
          mult_coeff = 0.5;
        }
      // for wall boundary conditions
      if(top_wall && bc_type!=Boundary_Conditions::Mixte_shear)
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

      // for wall(or mooving wall) boundary conditions
      if(top_wall && bc_type!=Boundary_Conditions::Mixte_shear)
        {
          if(_VCOMPO_ == DIRECTION::X)
            {
              v4 = ref_bc_.valeur().get_vx_kmax();
            }
          else
            {
              v4 = 0.;
            }
        }
      if(bottom_wall && bc_type!=Boundary_Conditions::Mixte_shear)
        {
          if(_VCOMPO_ == DIRECTION::X)
            {
              v3 = ref_bc_.valeur().get_vx_kmin();
            }
          else
            {
              v3 = 0.;
            }
        }


      Simd_double tau = (v4 - v3);
      if(!is_anisotropic_)
        tau *= inv_distance_DIR;

      Simd_double tau_tr = 0.;
      if(with_transpose_)
        {
          // gradient in direction COMPO of component DIR
          if(!bottom_wall && !top_wall)
            {
              Simd_double v1, v2;
              vDIR_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              tau_tr = (v2 - v1);
              if(!is_anisotropic_)
                tau_tr *= inv_distance_COMPO;
            }
        }

      Simd_double reyn = 0.;
      if(is_turb_)
        {
          Simd_double t_nu1, t_nu2, t_nu3, t_nu4;
          turbulent_nu.get_left_center_c1c2(_DIR_, _VCOMPO_, i, t_nu1, t_nu2, t_nu3, t_nu4);

          mult_coeff = 0.25;
          if(bottom_wall)
            {
              t_nu1 = 0., t_nu2 = 0.;
              mult_coeff = 0.5;
            }
          if(top_wall)
            {
              t_nu3 = 0., t_nu4 = 0.;
              mult_coeff = 0.5;
            }
          Simd_double t_nu = (t_nu1 + t_nu2 + t_nu3 + t_nu4) * mult_coeff;
          // gradient in direction COMPO of component DIR
          if(!bottom_wall && !top_wall)
            {
              Simd_double v1, v2;
              vDIR_ptr.get_left_center(_VCOMPO_, i, v1, v2);
              tau_tr = (v2 - v1) * inv_distance_COMPO;
            }
          reyn = (tau + tau_tr) * t_nu;
        }

      flux = (m_nu * (tau+tau_tr) + reyn) * surface;
    }
}

#endif
