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

#ifndef Operateur_IJK_base_included
#define Operateur_IJK_base_included

#include <IJK_Field.h>
#include <IJK_Field_simd_tools.h>
#include <IJK_Splitting.h>
#include <Operateur_IJK_data_channel.h>
#include <Corrige_flux_FT.h>

inline void putzero(IJK_Field_local_double& flux)
{
  const int vsize = Simd_double::size();
  const int ni = flux.ni();
  const int nj = flux.nj();
  double *ptr = &flux(0,0,0);
  double *end_ptr = &flux(ni-1,nj-1,0);
  const int n = (int)(end_ptr - ptr + 1);
  Simd_double zero = 0.;
  for (int i = 0; i < n; i += vsize)
    SimdPut(ptr+i, zero);
}

// Base class for operators on velocity field.
// It provides the structure for "conservative finite volume" operators:
//  The operator computes fluxes between the control volumes, then computes the
//  divergence of the fluxes (div_set() and div_add() method).
// The derived class must implement a public "compute" method that will call compute_set()
// or compute_add() (see for example OpDiffTurbIJK_double and OpConvCentre4IJK_double)
// It must also implement the compute_flux...() methods to compute fluxes between the control volumes.
class Operateur_IJK_faces_base_double : public Objet_U
{
  Declare_base( Operateur_IJK_faces_base_double );

public:
  virtual double compute_dtstab_convection_local(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  void compute_set(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
  void compute_add(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz);
protected:
  // The derived class must implement the computation of fluxes (9 methods)
  // See for example classes OpDiffTurbIJK_double and OpConvCentre4IJK_double.
  virtual void compute_flux_x_vx(IJK_Field_local_double& resu, const int k_layer) = 0;
  virtual void compute_flux_x_vy(IJK_Field_local_double& resu, const int k_layer) = 0;
  virtual void compute_flux_x_vz(IJK_Field_local_double& resu, const int k_layer) = 0;
  virtual void compute_flux_y_vx(IJK_Field_local_double& resu, const int k_layer) = 0;
  virtual void compute_flux_y_vy(IJK_Field_local_double& resu, const int k_layer) = 0;
  virtual void compute_flux_y_vz(IJK_Field_local_double& resu, const int k_layer) = 0;
  virtual void compute_flux_z_vx(IJK_Field_local_double& resu, const int k_layer) = 0;
  virtual void compute_flux_z_vy(IJK_Field_local_double& resu, const int k_layer) = 0;
  virtual void compute_flux_z_vz(IJK_Field_local_double& resu, const int k_layer) = 0;

  // The derived class might implement (optional) supplemental computations,
  // the method is called after computing the divergence of the flux
  // (example, use it to add velocity * divergence(velocity) for compressible flows)
  // Default implementation does nothing.
  virtual void exec_after_divergence_flux_x(IJK_Field_double& resu, const int k_layer) { };
  virtual void exec_after_divergence_flux_y(IJK_Field_double& resu, const int k_layer) { };
  virtual void exec_after_divergence_flux_z(IJK_Field_double& resu, const int k_layer) { };

private:
  void compute_(IJK_Field_double& dvx, IJK_Field_double& dvy, IJK_Field_double& dvz, bool add);

};

// Same, but base class for operator on cell centered scalar datapeeb
class Operateur_IJK_elem_base_double : public Objet_U
{
  Declare_base( Operateur_IJK_elem_base_double );

public:
  virtual void initialize(const IJK_Splitting& splitting)=0;
  virtual void compute_set(IJK_Field_double& dx);
  virtual void compute_add(IJK_Field_double& dx);
protected:
  // The derived class must implement the computation of fluxes (3 fluxes, one per direction)
  virtual void compute_flux_x(IJK_Field_local_double& resu, const int k_layer) = 0;
  virtual void compute_flux_y(IJK_Field_local_double& resu, const int k_layer) = 0;
  virtual void compute_flux_z(IJK_Field_local_double& resu, const int k_layer) = 0;

private:
  void compute_(IJK_Field_double& dx, bool add);

};

#endif
