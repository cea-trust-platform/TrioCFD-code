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
// File      : IJK_Thermal_Multiple_Subresolutions.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Thermal_Multiple_Subresolutions_included
#define IJK_Thermal_Multiple_Subresolutions_included

#include <IJK_Thermal_Subresolution.h>
#include <IJK_Field.h>
#include <Boundary_Conditions_Thermique.h>
#include <IJK_Splitting.h>
#include <IJK_Field.h>
#include <Parser.h>
#include <IJK_Lata_writer.h>
#include <OpConvQuickIJKScalar.h>
#include <OpConvCentre2IJKScalar.h>
#include <Ouvrir_fichier.h>
#include <Corrige_flux_FT.h>
#include <TRUST_Ref.h>
#include <Operateur_IJK_elem_diff_base.h>
#include <OpConvAmontIJK.h>
#include <OpConvDiscQuickIJKScalar.h>
#include <OpConvCentre4IJK.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_Thermal_Multiple_Subresolutions
//
// <Description of class IJK_Thermal_Multiple_Subresolutions>
//
/////////////////////////////////////////////////////////////////////////////


class IJK_Thermal_Multiple_Subresolutions : public IJK_Thermal_Subresolution
{

  Declare_instanciable( IJK_Thermal_Multiple_Subresolutions ) ;

public :

  int initialize(const IJK_Splitting& splitting, const int idx) override;
  void update_thermal_properties() override;
  void set_param(Param& param) override;

protected :

  /*
   * Treat the second phase like it is vapour !
   */

  IJK_Field_double temperature_vapour_;
  IJK_Field_double div_coeff_grad_T_vapour_volume_;
  FixedVector<IJK_Field_double, 3> grad_T_vapour_;

  void correct_temperature_vapour_for_eulerian_fluxes();

  OpDiffUniformIJKScalar_double diffusion_temperature_vapour_op_;

  int main_phase_;
  int diffusion_flux_vapour_correction_;
  int convective_flux_vapour_correction_;

  double uniform_lambda_vap_;
  double uniform_alpha_vap_;

};

#endif /* IJK_Thermal_Multiple_Subresolutions_included */
