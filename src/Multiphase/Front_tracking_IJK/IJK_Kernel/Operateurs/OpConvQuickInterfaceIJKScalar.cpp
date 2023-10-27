/****************************************************************************
 * Copyright (c) 2015 - 2016, CEA
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 *modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *this list of conditions and the following disclaimer in the documentation
 *and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *may be used to endorse or promote products derived from this software without
 *specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *****************************************************************************/
/////////////////////////////////////////////////////////////////////////////
//
// File      : OpConvIJKQuickScalar_interface.cpp
// Directory : $IJK_ROOT/src/IJK/OpVDF
//
/////////////////////////////////////////////////////////////////////////////
#include <IJK_Field_simd_tools.h>
#include <OpConvQuickInterfaceIJKScalar.h>
#include <Operateur_IJK_base.h>

Implemente_instanciable_sans_constructeur(OpConvQuickInterfaceIJKScalar_double, "OpConvQuickInterfaceIJKScalar_double", OpConvQuickIJKScalar_double);

Sortie& OpConvQuickInterfaceIJKScalar_double::printOn(Sortie& os) const
{
  OpConvQuickIJKScalar_double::printOn(os);
  return os;
}

Entree& OpConvQuickInterfaceIJKScalar_double::readOn(Entree& is)
{
  OpConvQuickIJKScalar_double::readOn(is);
  return is;
}

void OpConvQuickInterfaceIJKScalar_double::correct_flux(IJK_Field_local_double *const flux,	const int k_layer, const int dir)
{
  corrige_flux_->corrige_flux_faceIJ(flux, k_layer, dir);
}

void OpConvQuickInterfaceIJKScalar_double::correct_flux_spherical(Simd_double& a, Simd_double& b, const int& i, const int& j, int k_layer, int dir)
{
  corrige_flux_->correct_flux_spherical(a, b, i, j, k_layer, dir);
}
