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
// File      : OpHessCentre2IJKScalar.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#ifndef OpHessCentre2IJKScalar_included
#define OpHessCentre2IJKScalar_included

#include <Operateur_IJK_elem_diff_base.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class OpHessCentre2IJKScalar
//
// <Description of class OpHessCentre2IJKScalar>
//
/////////////////////////////////////////////////////////////////////////////

class OpHessCentre2IJKScalar_double : public OpDiffUniformIJKScalar_double
{

  Declare_instanciable( OpHessCentre2IJKScalar_double ) ;

public:
  void calculer_hess(const IJK_Field_double& field,
                     FixedVector<IJK_Field_double, 3>& result,
                     const IJK_Field_local_double& boundary_flux_kmin,
                     const IJK_Field_local_double& boundary_flux_kmax);

  void calculer_hess_xx(const IJK_Field_double& field,
                        IJK_Field_double& result,
                        const IJK_Field_local_double& boundary_flux_kmin,
                        const IJK_Field_local_double& boundary_flux_kmax);

  void calculer_hess_yy(const IJK_Field_double& field,
                        IJK_Field_double& result,
                        const IJK_Field_local_double& boundary_flux_kmin,
                        const IJK_Field_local_double& boundary_flux_kmax);

  void calculer_hess_zz(const IJK_Field_double& field,
                        IJK_Field_double& result,
                        const IJK_Field_local_double& boundary_flux_kmin,
                        const IJK_Field_local_double& boundary_flux_kmax);
protected:
  const double unit_lambda_ = 1.;
  void fill_grad_field_x_y_(IJK_Field_local_double& flux, IJK_Field_double& resu, int k, int dir) override;
  void fill_grad_field_z_(IJK_Field_local_double& flux_min, IJK_Field_local_double& flux_max, IJK_Field_double& resu, int k) override;
};

#endif /* OpHessCentre2IJKScalar_included */
