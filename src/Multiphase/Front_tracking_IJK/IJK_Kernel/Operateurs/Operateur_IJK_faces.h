/****************************************************************************
* Copyright (c) 2023, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Operateur_IJK_faces.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Operateur_IJK_faces_included
#define Operateur_IJK_faces_included

#include <Operateur_IJK_base.h>
#include <TRUST_Deriv.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Operateur_IJK_faces
//
// <Description of class Operateur_IJK_faces>
//
/////////////////////////////////////////////////////////////////////////////

class Operateur_IJK_faces : public DERIV(Operateur_IJK_faces_base_double)
{

  Declare_instanciable( Operateur_IJK_faces ) ;

public:
  inline double compute_dtstab_convection_local(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  inline void compute_set(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  inline void compute_add(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);

//protected:
//  // The derived class must implement the computation of fluxes (9 methods)
//  // See for example classes OpDiffTurbIJK_double and OpConvCentre4IJK_double.
//  inline void compute_flux_x_vx(IJK_Field_local_double& resu, const int k_layer);
//  inline void compute_flux_x_vy(IJK_Field_local_double& resu, const int k_layer);
//  inline void compute_flux_x_vz(IJK_Field_local_double& resu, const int k_layer);
//  inline void compute_flux_y_vx(IJK_Field_local_double& resu, const int k_layer);
//  inline void compute_flux_y_vy(IJK_Field_local_double& resu, const int k_layer);
//  inline void compute_flux_y_vz(IJK_Field_local_double& resu, const int k_layer);
//  inline void compute_flux_z_vx(IJK_Field_local_double& resu, const int k_layer);
//  inline void compute_flux_z_vy(IJK_Field_local_double& resu, const int k_layer);
//  inline void compute_flux_z_vz(IJK_Field_local_double& resu, const int k_layer);
//
//  // The derived class might implement (optional) supplemental computations,
//  // the method is called after computing the divergence of the flux
//  // (example, use it to add velocity * divergence(velocity) for compressible flows)
//  // Default implementation does nothing.
//  inline void exec_after_divergence_flux_x(IJK_Field_double& resu, const int k_layer);
//  inline void exec_after_divergence_flux_y(IJK_Field_double& resu, const int k_layer);
//  inline void exec_after_divergence_flux_z(IJK_Field_double& resu, const int k_layer);
//
//private:
//  inline void compute_(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz, bool add);

};

inline double Operateur_IJK_faces::compute_dtstab_convection_local(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  return valeur().compute_dtstab_convection_local(dvx, dvy, dvz);
}

inline void Operateur_IJK_faces::compute_set(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  valeur().compute_set(dvx, dvy, dvz);
}

inline void Operateur_IJK_faces::compute_add(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz)
{
  valeur().compute_add(dvx, dvy, dvz);
}

//inline void Operateur_IJK_faces::compute_flux_x_vx(IJK_Field_local_double& resu, const int k_layer)
//{
//  valeur().compute_flux_x_vx(resu, k_layer);
//}
//
//inline void Operateur_IJK_faces::compute_flux_x_vy(IJK_Field_local_double& resu, const int k_layer)
//{
//  valeur().compute_flux_x_vy(resu, k_layer);
//}
//
//inline void Operateur_IJK_faces::compute_flux_x_vz(IJK_Field_local_double& resu, const int k_layer)
//{
//  valeur().compute_flux_x_vz(resu, k_layer);
//}
//
//inline void Operateur_IJK_faces::compute_flux_y_vx(IJK_Field_local_double& resu, const int k_layer)
//{
//  valeur().compute_flux_y_vx(resu, k_layer);
//}
//
//inline void Operateur_IJK_faces::compute_flux_y_vy(IJK_Field_local_double& resu, const int k_layer)
//{
//  valeur().compute_flux_y_vy(resu, k_layer);
//}
//
//inline void Operateur_IJK_faces::compute_flux_y_vz(IJK_Field_local_double& resu, const int k_layer)
//{
//  valeur().compute_flux_y_vz(resu, k_layer);
//}
//
//inline void Operateur_IJK_faces::compute_flux_z_vx(IJK_Field_local_double& resu, const int k_layer)
//{
//  valeur().compute_flux_z_vx(resu, k_layer);
//}
//
//inline void Operateur_IJK_faces::compute_flux_z_vy(IJK_Field_local_double& resu, const int k_layer)
//{
//  valeur().compute_flux_z_vy(resu, k_layer);
//}
//
//inline void Operateur_IJK_faces::compute_flux_z_vz(IJK_Field_local_double& resu, const int k_layer)
//{
//  valeur().compute_flux_z_vz(resu, k_layer);
//}
//
//inline void Operateur_IJK_faces::exec_after_divergence_flux_x(IJK_Field_double& resu, const int k_layer)
//{
//  valeur().exec_after_divergence_flux_x(resu, k_layer);
//}
//
//inline void Operateur_IJK_faces::exec_after_divergence_flux_y(IJK_Field_double& resu, const int k_layer)
//{
//  valeur().exec_after_divergence_flux_y(resu, k_layer);
//}
//
//inline void Operateur_IJK_faces::exec_after_divergence_flux_z(IJK_Field_double& resu, const int k_layer)
//{
//  valeur().exec_after_divergence_flux_y(resu, k_layer);
//}
//
//inline void Operateur_IJK_faces::compute_(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz, bool add)
//{
//  valeur().compute_(dvx, dvy, dvz, add);
//}

#endif /* Operateur_IJK_faces_included */
