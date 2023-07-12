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
// File      : Operateur_IJK_elem_conv.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Operateur_IJK_elem_conv_included
#define Operateur_IJK_elem_conv_included

#include <Operateur_IJK_elem.h>
#include <OpConvIJKElemCommon.h>
#include <TRUST_Deriv.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Operateur_IJK_elem_conv
//
// <Description of class Operateur_IJK_elem_conv>
//
/////////////////////////////////////////////////////////////////////////////

class Operateur_IJK_elem_conv : public Operateur_IJK_elem
{
  Declare_instanciable( Operateur_IJK_elem_conv ) ;
public :
  DERIV(OpConvIJKElemCommon_double) OpConvIJKElemDeriv;
  inline void typer(const char * type);
  inline void set_corrige_flux(Corrige_flux_FT& corrige_flux);
  inline void calculer(const IJK_Field_double& field,
                       const IJK_Field_double& vx,
                       const IJK_Field_double& vy,
                       const IJK_Field_double& vz,
                       IJK_Field_double& result);
  inline void ajouter(const IJK_Field_double& field,
                      const IJK_Field_double& vx,
                      const IJK_Field_double& vy,
                      const IJK_Field_double& vz,
                      IJK_Field_double& result);
};

inline void Operateur_IJK_elem_conv::typer(const char * type)
{
  Operateur_IJK_elem::typer(type);
  OpConvIJKElemDeriv.typer(type);
}

inline void Operateur_IJK_elem_conv::set_corrige_flux(Corrige_flux_FT& corrige_flux)
{
  OpConvIJKElemDeriv.valeur().set_corrige_flux(corrige_flux);
}

inline void Operateur_IJK_elem_conv::calculer(const IJK_Field_double& field,
                                              const IJK_Field_double& vx,
                                              const IJK_Field_double& vy,
                                              const IJK_Field_double& vz,
                                              IJK_Field_double& result)
{
  OpConvIJKElemDeriv.valeur().calculer(field, vx, vy, vz, result);
}


inline void Operateur_IJK_elem_conv::ajouter(const IJK_Field_double& field,
                                             const IJK_Field_double& vx,
                                             const IJK_Field_double& vy,
                                             const IJK_Field_double& vz,
                                             IJK_Field_double& result)
{
  OpConvIJKElemDeriv.valeur().ajouter(field, vx, vy, vz, result);
}

#endif /* Operateur_IJK_elem_conv_included */
