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

#ifndef Operateur_IJK_data_channel_h
#define Operateur_IJK_data_channel_h
#include <IJK_Field.h>
#include <IJK_Splitting.h>

// Data used by IJK operators in the case of a
//  biperiodic channel uniform mesh in i and j,
//  wall in k direction.
class Operateur_IJK_data_channel
{
public:
  void initialize(const IJK_Splitting& splitting);

  int nb_elem_k_tot() const
  {
    return nb_elem_k_tot_;
  }

  // First layer of non zero fluxes (eg, not in or at the walls)
  // (specific for wall boundary conditions at zmin and zmax)
  int first_global_k_layer_flux(int icompo, int idir) const
  {
    if (idir != 2 && icompo != 2) // not a trivial thing to determine...
      return 0;
    else
      return 1;
  }
  int last_global_k_layer_flux(int icompo, int idir) const
  {
    if (idir == 2 && icompo == 2) // not a trivial thing to determine...
      return nb_elem_k_tot_;
    else
      return nb_elem_k_tot_ - 1;
  }
  int offset_to_global_k_layer() const
  {
    return offset_to_global_k_layer_;
  }
  // Returns the surface between the two control volumes for the requested flux
  // (flux in direction dir for velocity component compo, flux layer k_layer)
  double get_surface(int k_layer, int compo, int dir)
  {
    // specific coding for uniform mesh in i and j, variable mesh in k:
    switch(compo)
      {
      case 0:
      case 1:
        switch(dir)
          {
          case 0:
            return delta_y_ * delta_z_[k_layer];
          case 1:
            return delta_x_ * delta_z_[k_layer];
          case 2:
          default:
            return delta_x_ * delta_y_;
          }
      case 2:
      default:
        switch(dir)
          {
          case 0:
            return delta_y_ * (delta_z_[k_layer - 1] + delta_z_[k_layer]) * 0.5;
          case 1:
            return delta_x_ * (delta_z_[k_layer - 1] + delta_z_[k_layer]) * 0.5;
          case 2:
          default:
            return delta_x_ * delta_y_;
          }
      }
  }
  // Returns the surfaces of the two ajacent faces to the edge for the requested flux,
  // (faces oriented in direction "dir", adjacent in direction "compo")
  void get_surface_leftright(int k_layer, int compo, int dir, double& left, double& right)
  {
    // specific coding for uniform mesh in i and j, variable mesh in k:
    switch(compo)
      {
      case 0:
      case 1:
        switch(dir)
          {
          case 0:
            left = right = delta_y_ * delta_z_[k_layer];
            break;
          case 1:
            left = right = delta_x_ * delta_z_[k_layer];
            break;
          case 2:
          default:
            left = right = delta_x_ * delta_y_;
          }
        break;
      case 2:
      default:
        switch(dir)
          {
          case 0:
            left  = delta_y_ * delta_z_[k_layer - 1];
            right = delta_y_ * delta_z_[k_layer];
            break;
          case 1:
            left  = delta_x_ * delta_z_[k_layer - 1];
            right = delta_x_ * delta_z_[k_layer];
            break;
          case 2:
          default:
            left = right = delta_x_ * delta_y_;
          }
        break;
      }
  }
  // Returns the inverse of the distance between velocity dof (component "compo",
  //  for gradient in direction "dir")
  double inv_distance_for_gradient(int k_layer, int compo, int dir)
  {
    switch(dir)
      {
      case 0:
        return inv_delta_x_;
      case 1:
        return inv_delta_y_;
      default:
        if (compo == 2)
          // flux layer k is, by convention, between face k-1 and face k, which is within element k-1:
          return inv_elem_size_z_[k_layer-1];
        else
          return inv_dist_z_elemcenter_[k_layer];
      }
  }
  const ArrOfDouble_with_ghost& get_delta_z() const
  {
    return delta_z_;
  }
  double get_delta_x() const
  {
    return delta_x_;
  }
  double get_delta_y() const
  {
    return delta_y_;
  }
  double get_delta(int dir, int k_layer=-1) const
  {
    switch (dir)
      {
      case 0:
        return delta_x_;
      case 1:
        return delta_y_;
      case 2:
        {
          if(k_layer==-1)
            {
              Cerr << "Operateur_IJK_data_channel:get_delta:: invalid k_layer" << finl;
              Process::exit();
            }
          return delta_z_[k_layer];
        }
      default:
        {
          Cerr << "Operateur_IJK_data_channel:get_delta:: invalid direction" << finl;
          Process::exit();
          return 0.;
        }
      }
  }

  double get_delta_xyz(int k_layer, int dir)
  {
    switch (dir)
      {
      case 0:
        return delta_x_;
      case 1:
        return delta_y_;
      case 2:
        return delta_z_[k_layer];
      default:
        return delta_x_;
      }
  }

protected:
  // Total number of mesh cells in the k direction (in global mesh)
  int nb_elem_k_tot_;
  // adding this offset to a local k_layer we known where we are in the global mesh
  int offset_to_global_k_layer_;
  // Mesh sizes for surface and distance computations:
  double delta_x_; // uniform in i and j directions
  double delta_y_;
  ArrOfDouble_with_ghost delta_z_; // Size of mesh cells in z direction

  double inv_delta_x_; // shortcut for inverse of cell size...
  double inv_delta_y_;
  // inverse of size of elements (distance between k-oriented faces in k direction)
  ArrOfDouble_with_ghost inv_elem_size_z_;
  // inverse of distance between element centers (and distance to wall at boundaries)
  // at index k we have the distance between elements k-1 and k
  ArrOfDouble_with_ghost inv_dist_z_elemcenter_;
};
#endif
