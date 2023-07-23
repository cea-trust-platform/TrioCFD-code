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

#ifndef OpConvDiscQuickIJKScalar_included
#define OpConvDiscQuickIJKScalar_included

#include <Operateur_IJK_elem_conv_base.h>

class OpConvDiscQuickIJKScalar_double : public Operateur_IJK_elem_conv_base_double
{
  Declare_instanciable_sans_constructeur(OpConvDiscQuickIJKScalar_double);
public:
  OpConvDiscQuickIJKScalar_double() : Operateur_IJK_elem_conv_base_double() {  };
  void calculer(const IJK_Field_double& field,
                const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                IJK_Field_double& result) override;
  void ajouter(const IJK_Field_double& field,
               const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
               IJK_Field_double& result) override;
protected:

  inline void compute_flux_x(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::X>(resu,k_layer);
  }
  inline void compute_flux_y(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Y>(resu,k_layer);
  }
  inline void compute_flux_z(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Z>(resu,k_layer);
  }

  const IJK_Field_local_double *input_indicatrice_=0;

private:

  template <DIRECTION _DIR_>
  void compute_flux_(IJK_Field_local_double& resu, const int k_layer);

};

#include <OpConvDiscQuickIJKScalar.tpp>

#endif /* OpConvDiscQuickIJKScalar_included */
