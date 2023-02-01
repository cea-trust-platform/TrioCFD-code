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
#ifndef OpCentre4IJK_H_TPP
#define OpCentre4IJK_H_TPP
#include <iostream>

// Methode appelee a chaque couche de vitesses calculees, apres le calcul de la divergence du flux
//  div(u x rho_u) pour la composante _DIR_ de vitesse, par Operateur_IJK_faces_base_double::compute_
// On ajoute ici u * div(rho_u) si c'est necessaire
template <DIRECTION _DIR_>
void OpConvCentre4IJK_double::exec_after_divergence_flux_(IJK_Field_double& resu, const int k_layer)
{

  if (div_rho_u_ == 0)
    return;
  std::cout << " " << std::endl;
  std::cout << "exec_after_divergence_flux_" << std::endl;
  std::cout << " " << std::endl;
  if(_DIR_==DIRECTION::Z)
    {
      const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
      // global index of the layer of flux of the wall
      //  (assume one walls at zmin and zmax)
      const int first_global_k_layer = 0;
      const int last_global_k_layer = resu.get_splitting().get_nb_items_global(IJK_Splitting::FACES_K, DIRECTION_K) - 1;

      if (!perio_k_ && (global_k_layer <= first_global_k_layer || global_k_layer >= last_global_k_layer))
        {
          return; // pas de calcul pour les faces en paroi (de toutes facons, v=0)
        }
    }
  // Calcul de div(rho_u) sur la couche, si besoin
  if (last_computed_klayer_for_div_rhou_ < k_layer)
    {
      // Il faut calculer une couche de div_rhou:
      calculer_div_rhou(*inputx_, *inputy_, *inputz_, *div_rho_u_, k_layer, channel_data_);
      last_computed_klayer_for_div_rhou_ = k_layer;
    }
  ConstIJK_double_ptr vitesse(get_v(_DIR_), 0, 0, k_layer);
  ConstIJK_double_ptr div_rhou(*div_rho_u_, 0, 0, k_layer);
  IJK_double_ptr resu_ptr(resu, 0, 0, k_layer);

  const int nx = resu.ni();
  const int ny = resu.nj();

  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double v, x_left, x_right;
          vitesse.get_center(i, v);
          std::cout << v.data_  << std::endl;
          // on prend div(rho_u) dans les elements a gauche et a droite de la face
          // rappel: l'element a gauche de la face i est a l'indice i-1
          div_rhou.get_left_center(_DIR_, i, x_left, x_right);
          // calcul du produit vitesse * div(rho_u)
          // en prenant div(rho_u) sur le volume de controle de la face (c'est l'integrale de div(rho_u)
          // qui est stocke, donc ici moyenne au sens volume fini)
          Simd_double a;
          resu_ptr.get_center(i, a);
          a += (x_left + x_right) * 0.5 * v;
          resu_ptr.put_val(i, a);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      // instructions to perform to jump to next row
      vitesse.next_j();
      div_rhou.next_j();
      resu_ptr.next_j();
    }


}

/*! @brief compute fluxes in direction DIR for velocity component COMPO for the layer of fluxes k_layer
 *
 *  4-th order centered convection scheme
 *
 */
template <DIRECTION _DIR_, DIRECTION _VCOMPO_>
void OpConvCentre4IJK_double::compute_flux_(IJK_Field_local_double& resu, const int k_layer)
{
  // convected field
  const IJK_Field_local_double& src = get_input(_VCOMPO_);
  // Convected vector field:
  ConstIJK_double_ptr src_ptr(src, 0, 0, k_layer);
  // Velocity in direction _DIR_ (convecting velocity)
  ConstIJK_double_ptr vconv_ptr(get_v(_DIR_), 0, 0, k_layer);

  const int idir = (int)_DIR_;
  const int nx = _DIR_ == DIRECTION::X ? src.ni() + 1 : src.ni();
  const int ny = _DIR_ == DIRECTION::Y ? src.nj() + 1 : src.nj();
  const int icompo = (int)_VCOMPO_;

  // Result (fluxes in direction x for component x of the convected field)
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);

  const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
  // global index of the layer of flux of the wall
  //  (assume one walls at zmin and zmax)
  const int first_global_k_layer = channel_data_.first_global_k_layer_flux(icompo, idir);
  const int last_global_k_layer = channel_data_.last_global_k_layer_flux(icompo, idir);
  Boundary_Conditions::BCType bc_type = ref_bc_.valeur().get_bctype_k_min();


  if (!perio_k_ && (global_k_layer <= first_global_k_layer || global_k_layer >= last_global_k_layer))
    {
      // We are in the wall
      putzero(resu);
      return;
    }

  // constant or variable coefficients depending on mesh.
  // second order degeneration is handeled here by setting appropriate coefficients:
  double g1, g2, g3, g4;
  double constant_factor0, constant_factor1;
  // surface of left face will be multiplied by velocity of right face, and reversed...
  channel_data_.get_surface_leftright(k_layer, icompo, idir, constant_factor1, constant_factor0);
  constant_factor0 *= 0.5;
  constant_factor1 *= 0.5;

  if(_DIR_ != DIRECTION::Z)
    {
      // Specific coding for uniform mesh in i and j, variable mesh in k,
      //  periodic in i and j, walls at bottom and top of k
      // Uniform mesh periodic everywhere so no special case.
      // We always have enough data for order 4:
      g1 = g4 = -0.0625;
      g2 = g3 = 0.5625;
    }
  else
    {
      // Specific coding for uniform mesh in i and j, variable mesh in k,
      //  periodic in i and j, walls at bottom and top of k
      g1 = get_g(k_layer,icompo,idir,0);
      g2 = get_g(k_layer,icompo,idir,1);
      g3 = get_g(k_layer,icompo,idir,2);
      g4 = get_g(k_layer,icompo,idir,3);
    }

  {
    const int imax = nx;
    const int jmax = ny;
    const int vsize = Simd_double::size();
    for (int j = 0; ; j++)
      {
        for (int i = 0; i < imax; i += vsize)     // specific coding for uniform mesh in x and y: surfaces are constant on an xy plane
          {
            Simd_double vit_0_0,vit_0,vit_1,vit_1_1; // 4 adjacent velocity values
            src_ptr.get_leftleft_left_center_right(_DIR_,i,vit_0_0,vit_0,vit_1,vit_1_1);


            if(_DIR_ == DIRECTION::Z && _VCOMPO_ == DIRECTION::X && bc_type==Boundary_Conditions::Mixte_shear)
              {
                if(global_k_layer == first_global_k_layer-1)
                  {
                    // std::cout << "first global k layer" << std::endl;
                    // std::cout << vit_0_0.data_ << vit_0.data_ << vit_1.data_ << vit_1_1.data_ << std::endl;
                    vit_0_0 +=ref_bc_.valeur().get_sh_t() * 1.;
                    vit_0 +=ref_bc_.valeur().get_sh_t() * 1.;
                    // std::cout << vit_0_0.data_ << vit_0.data_ << vit_1.data_ << vit_1_1.data_ << std::endl;
                  }

                else if(global_k_layer == first_global_k_layer)
                  {
                    // std::cout << "Second global k layer" << std::endl;
                    // std::cout << vit_0_0.data_ << vit_0.data_ << vit_1.data_ << vit_1_1.data_ << std::endl;
                    vit_0_0 +=ref_bc_.valeur().get_sh_t() * 1.;
                    // std::cout << vit_0_0.data_ << vit_0.data_ << vit_1.data_ << vit_1_1.data_ << std::endl;
                  }

                else if(global_k_layer == last_global_k_layer + 1)
                  {
                    // std::cout << "last global k layer + 1" << std::endl;
                    // std::cout << vit_0_0.data_ << vit_0.data_ << vit_1.data_ << vit_1_1.data_ << std::endl;
                    vit_1_1 -=ref_bc_.valeur().get_sh_t() * 1.;
                    vit_1 -=ref_bc_.valeur().get_sh_t() * 1.;
                    // std::cout << vit_0_0.data_ << vit_0.data_ << vit_1.data_ << vit_1_1.data_ << std::endl;
                  }
                else if(global_k_layer == last_global_k_layer)
                  {
                    // std::cout << "last global k layer" << std::endl;
                    // std::cout << vit_0_0.data_ << vit_0.data_ << vit_1.data_ << vit_1_1.data_ << std::endl;
                    vit_1_1 -=ref_bc_.valeur().get_sh_t() * 1.;
                    // std::cout << vit_0_0.data_ << vit_0.data_ << vit_1.data_ << vit_1_1.data_ << std::endl;
                  }

              }

            // if (vit_0_0.data_<0. && vit_0.data_ >0. && !(global_k_layer==23 || global_k_layer==24 || global_k_layer==25))
            //  {
            //    std::cout << global_k_layer << vit_0_0.data_ << vit_0.data_ << std::endl;
            //   }
            //  f (vit_0_0.data_>0. && vit_0.data_ <0. && !(global_k_layer==23 || global_k_layer==24 || global_k_layer==25))
            //   {
            //     std::cout << global_k_layer << vit_0_0.data_ << vit_0.data_ << std::endl;
            //   }
            //  if (vit_1.data_<0. && vit_0.data_ >0. && !(global_k_layer==23 || global_k_layer==24 || global_k_layer==25))
            //    {
            //      std::cout << global_k_layer << vit_1.data_ << vit_0.data_ << std::endl;
            //    }
            //  if (vit_1.data_>0. && vit_0.data_ <0. && !(global_k_layer==23 || global_k_layer==24 || global_k_layer==25))
            //    {
            //       std::cout << global_k_layer << vit_1.data_ << vit_0.data_ << std::endl;
            //     }
            //   if (vit_1.data_<0. && vit_1_1.data_ >0. && !(global_k_layer==23 || global_k_layer==24 || global_k_layer==25))
            //     {
            //      std::cout << global_k_layer << vit_1.data_ << vit_1_1.data_ << std::endl;
            //     }
            //   if (vit_1.data_>0. && vit_1_1.data_ <0. && !(global_k_layer==23 || global_k_layer==24 || global_k_layer==25))
            //     {
            //       std::cout << global_k_layer << vit_1.data_ << vit_1_1.data_ << std::endl;
            //    }

            Simd_double order4_velocity = g1 * vit_0_0 + g2 * vit_0 + g3 * vit_1 + g4 * vit_1_1;

            // get convecting velocity
            Simd_double vconv0, vconv1;
            vconv_ptr.get_left_center(_VCOMPO_,i, vconv0, vconv1);
            if(_DIR_ == DIRECTION::X && _VCOMPO_ == DIRECTION::Z && bc_type==Boundary_Conditions::Mixte_shear)
              {
                if(global_k_layer == first_global_k_layer-1)
                  {
                    // std::cout << "first global k layer" << std::endl;
                    // std::cout << vconv0.data_  << vconv1.data_  << std::endl;
                    vconv0 +=ref_bc_.valeur().get_sh_t() * 1.;
                    //        std::cout << vconv0.data_  << vconv1.data_  << std::endl;
                    // std::cout << " "<< std::endl;
                  }
              }

            //   if (vconv1.data_<0. && vconv0.data_ >0. && global_k_layer!=24)
            //    {
            //      std::cout << global_k_layer << vconv1.data_  << vconv0.data_  << std::endl;
            //    }
            //  if (vconv1.data_>0. && vconv0.data_ <0. && global_k_layer!=24)
            //    {
            //      std::cout << global_k_layer << vconv1.data_  << vconv0.data_  << std::endl;
            //    }

            Simd_double psc = vconv0 * constant_factor0 + vconv1 * constant_factor1;
            // with porosity we would code this: vconv = (vconv0 * porosity0 + vconv1 * porosity1) * constant_factor;
            Simd_double flux_conv = order4_velocity * psc;

            resu_ptr.put_val(i, flux_conv);
          }
        // do not execute end_iloop at last iteration (because of assert on valid j+1)
        if (j+1==jmax)
          break;
        // instructions to perform to jump to next row
        src_ptr.next_j();
        resu_ptr.next_j();
        vconv_ptr.next_j();
      }
  }
}


#endif
