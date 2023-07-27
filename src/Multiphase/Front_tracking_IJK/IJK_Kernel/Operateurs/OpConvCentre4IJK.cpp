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

#include <OpConvCentre4IJK.h>

Implemente_instanciable_sans_constructeur(OpConvCentre4IJK_double, "OpConvCentre4IJK_double", Operateur_IJK_faces_conv_base_double);

OpConvCentre4IJK_double::OpConvCentre4IJK_double()
{
  div_rho_u_=0;
  last_computed_klayer_for_div_rhou_=0;
}

Sortie& OpConvCentre4IJK_double::printOn(Sortie& os) const
{
  return os;
}

Entree& OpConvCentre4IJK_double::readOn(Entree& is)
{
  return is;
}

inline void calcul_g(const double& dxam, const double& dx, const double& dxav, double& g1, double& g2, double& g3, double& g4)
{
  g1 = -dx*dx*(dx/2+dxav)/(4*(dx+dxam+dxav)*(dx+dxam)*dxam);
  g2 =  (dx+2*dxam)*(dx+2*dxav)/(8*dxam*(dx+dxav));
  g3 =  (dx+2*dxam)*(dx+2*dxav)/(8*dxav*(dx+dxam));
  g4 = -dx*dx*(dx/2+dxam)/(4*(dx+dxam+dxav)*dxav*(dx+dxav));
}

// g contains one line per flux computed in the z direction.
// The first flux at index 0 is juste before the first face owned by the processor.
// offset is the index in the global mesh of the first element on this processor in direction z
// istart, iend: index of the first and last fluxes computed with 4-th order, others are 2-nd order.
// is_z_component: shall we compute coefficients for interpolation for velocity_z or for velocity_x and velocity_y ?
// delta_z is the size of the cells on this processor, we need 2 ghost cells
static void fill_g_compo(DoubleTab& g, int nb_values, int offset,
                         int istart, int iend,
                         const ArrOfDouble_with_ghost& delta_z, bool is_z_component)
{
  g.resize(nb_values, 4);
  for (int i = 0; i < nb_values; i++)
    {
      if (i + offset < istart || i + offset > iend)
        {
          // We are in the wall or in the first layer: degenerate coefficients for 2nd order interpolation
          g(i,0) = 0.;
          g(i,1) = 0.5;
          g(i,2) = 0.5;
          g(i,3) = 0.;
        }
      else
        {
          double d1, d2, d3;
          if (is_z_component)
            {
              // Filtering coefficients for faces oriented in z: interval between faces is the cell size
              d1 = delta_z[i-2];
              d2 = delta_z[i-1];
              d3 = delta_z[i];
            }
          else
            {
              // Filtering coefficients for faces oriented in x or y: interval between centers of faces is this:
              d1 = (delta_z[i-2] + delta_z[i-1]) * 0.5;
              d2 = (delta_z[i-1] + delta_z[i]) * 0.5;
              d3 = (delta_z[i]   + delta_z[i+1]) * 0.5;
            }
          calcul_g(d1, d2, d3, g(i,0), g(i,1), g(i,2), g(i,3));
        }
    }
}

void OpConvCentre4IJK_double::initialize(const IJK_Splitting& splitting) //, const Boundary_Conditions& bc)
{
  Operateur_IJK_faces_conv_base_double::initialize(splitting);

  // Fill 4-th order filtering coefficients for z direction:
  const int nb_xfaces = splitting.get_nb_faces_local(0 /* for component x */, 2 /* in direction z */);
  const int nb_zfaces = splitting.get_nb_faces_local(2 /* for component z */, 2 /* in direction z */);
  // number of flux values computed on this processor in direction k
  // equals the number of faces owned by the processor, plus 1

  const int offset_to_global_k_layer = channel_data_.offset_to_global_k_layer();
  const ArrOfDouble_with_ghost& delta_z = channel_data_.get_delta_z();

  // first flux computed with 4th order is 1 layer after the first non zero flux, after the wall.
  // SPECIFIC FOR CHANNEL WITH WALLS IN K DIRECTION !
  int istart, iend;
  istart = channel_data_.first_global_k_layer_flux(0 /* compo */, 2 /* dir */) + 1;
  iend = channel_data_.last_global_k_layer_flux(0 /* compo */, 2 /* dir */) - 1;
  fill_g_compo(g_compo_xy_dir_z_, nb_xfaces + 1, offset_to_global_k_layer, istart, iend, delta_z, false);

  istart = channel_data_.first_global_k_layer_flux(2 /* compo */, 2 /* dir */) + 1;
  iend = channel_data_.last_global_k_layer_flux(2 /* compo */, 2 /* dir */) - 1;
  fill_g_compo(g_compo_z_dir_z_, nb_zfaces + 1, offset_to_global_k_layer, istart, iend, delta_z, true);

//  ref_bc_ = bc;
}

void OpConvCentre4IJK_double::calculer(const IJK_Field_double& inputx, const IJK_Field_double& inputy, const IJK_Field_double& inputz,
                                       const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                       IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  Operateur_IJK_faces_conv_base_double::calculer(inputx, inputy, inputz, vx, vy, vz, dvx, dvy, dvz);
  div_rho_u_ = 0;
}

void OpConvCentre4IJK_double::ajouter(const IJK_Field_double& inputx, const IJK_Field_double& inputy, const IJK_Field_double& inputz,
                                      const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                      IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  Operateur_IJK_faces_conv_base_double::ajouter(inputx, inputy, inputz, vx, vy, vz, dvx, dvy, dvz);
  div_rho_u_ = 0;
}

// Calcule, sur la couche k_layer de mailles, l'integrale sur la maille de div(rho_v)
// On calcule une epaisseur de mailles ghost a gauche dans les directions i et j
// (pour calcul de div(rhou) aux faces)
// On suppose que c'est periodique en i en j
void OpConvCentre4IJK_double::calculer_div_rhou(const IJK_Field_double& rhovx, const IJK_Field_double& rhovy, const IJK_Field_double& rhovz,
                                                IJK_Field_double& resu, int k_layer, const Operateur_IJK_data_channel& channel)
{
  statistiques().begin_count(convection_counter_);
  const double surface_x = channel.get_delta_y() * channel.get_delta_z()[k_layer];
  const double surface_y = channel.get_delta_x() * channel.get_delta_z()[k_layer];
  const double surface_z = channel.get_delta_x() * channel.get_delta_y();
  const int ni = resu.ni();
  const int nj = resu.nj();
  // codage simple, non vectorise:
  for (int j = -1; j < nj; j++)
    for (int i = -1; i < ni; i++)
      {
        double divergence =
          (rhovx(i+1,j,k_layer) - rhovx(i,j,k_layer)) * surface_x
          + (rhovy(i,j+1,k_layer) - rhovy(i,j,k_layer)) * surface_y
          + (rhovz(i,j,k_layer+1) - rhovz(i,j,k_layer)) * surface_z;
        resu(i,j,k_layer) = divergence;
      }
  statistiques().end_count(convection_counter_);

}

void OpConvCentre4IJK_double::calculer_avec_u_div_rhou(const IJK_Field_double& rhovx, const IJK_Field_double& rhovy, const IJK_Field_double& rhovz,
                                                       const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                       IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz,
                                                       IJK_Field_double& div_rho_u)
{
  statistiques().begin_count(convection_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;
  inputx_ = &rhovx;
  inputy_ = &rhovy;
  inputz_ = &rhovz;
  div_rho_u_ = &div_rho_u;
  last_computed_klayer_for_div_rhou_ = -1;
  calculer_div_rhou(*inputx_, *inputy_, *inputz_, *div_rho_u_, -1, channel_data_);

  compute_set(dvx, dvy, dvz);

  vx_ = vy_ = vz_ = inputx_ = inputy_ = inputz_ = 0;
  div_rho_u_ = 0;
  statistiques().end_count(convection_counter_);

}

void OpConvCentre4IJK_double::ajouter_avec_u_div_rhou(const IJK_Field_double& rhovx, const IJK_Field_double& rhovy, const IJK_Field_double& rhovz,
                                                      const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                      IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz,
                                                      IJK_Field_double& div_rho_u)
{
  statistiques().begin_count(convection_counter_);

  vx_ = &vx;
  vy_ = &vy;
  vz_ = &vz;
  inputx_ = &rhovx;
  inputy_ = &rhovy;
  inputz_ = &rhovz;
  div_rho_u_ = &div_rho_u;
  last_computed_klayer_for_div_rhou_ = -1;
  calculer_div_rhou(*inputx_, *inputy_, *inputz_, *div_rho_u_, -1, channel_data_);

  compute_add(dvx, dvy, dvz);

  vx_ = vy_ = vz_ = inputx_ = inputy_ = inputz_ = 0;
  div_rho_u_ = 0;
  statistiques().end_count(convection_counter_);

}
