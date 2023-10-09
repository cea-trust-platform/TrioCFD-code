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
// File      : OpConvQuickInterfaceOnefluidIJKScalar.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#include <OpConvQuickInterfaceOnefluidIJKScalar.h>

Implemente_instanciable_sans_constructeur( OpConvQuickInterfaceOnefluidIJKScalar_double, "OpConvQuickInterfaceOnefluidIJKScalar_double", OpConvQuickIJKScalar_double ) ;

Sortie& OpConvQuickInterfaceOnefluidIJKScalar_double::printOn( Sortie& os ) const
{
  OpConvQuickIJKScalar_double::printOn( os );
  return os;
}

Entree& OpConvQuickInterfaceOnefluidIJKScalar_double::readOn( Entree& is )
{
  OpConvQuickIJKScalar_double::readOn( is );
  return is;
}

void OpConvQuickInterfaceOnefluidIJKScalar_double::correct_flux(IJK_Field_local_double *const flux,	const int k_layer, const int dir)
{
  /*
   * Aymeric indices to correct fluxes
   */
  if (dir==2)
    if (k_layer==0)
      {
        corrige_flux_->update();
        corrige_flux_->corrige_flux_faceIJ(flux, -1, dir);
      }
    else
      corrige_flux_->corrige_flux_faceIJ(flux, k_layer, dir);
  else
    corrige_flux_->corrige_flux_faceIJ(flux, k_layer, dir);
}
