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
// File      : Operateur_IJK_elem.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Operateur_IJK_included
#define Operateur_IJK_included

#include <Operateur_IJK_base.h>
#include <TRUST_Deriv.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Operateur_IJK_elem
//
// <Description of class Operateur_IJK_elem>
//
/////////////////////////////////////////////////////////////////////////////

class Operateur_IJK_elem : public DERIV(Operateur_IJK_elem_base_double)
{

  Declare_instanciable( Operateur_IJK_elem );

public:
  inline void initialize(const IJK_Splitting& splitting);
  inline void compute_set(IJK_Field_double& dx);
  inline void compute_add(IJK_Field_double& dx);
//  inline void compute_flux_x(IJK_Field_local_double& resu, const int k_layer);
//  inline void compute_flux_y(IJK_Field_local_double& resu, const int k_layer);
//  inline void compute_flux_z(IJK_Field_local_double& resu, const int k_layer);
//  inline void compute_(IJK_Field_double& dx, bool add);

};

inline void Operateur_IJK_elem::initialize(const IJK_Splitting& splitting)
{
  valeur().initialize(splitting);
}

inline void Operateur_IJK_elem::compute_set(IJK_Field_double& dx)
{
  valeur().compute_set(dx);
}

inline void Operateur_IJK_elem::compute_add(IJK_Field_double& dx)
{
  valeur().compute_add(dx);
}

//inline void Operateur_IJK_elem::compute_flux_x(IJK_Field_local_double& resu, const int k_layer)
//{
//  valeur().compute_flux_x(resu, k_layer);
//}
//
//inline void Operateur_IJK_elem::compute_flux_y(IJK_Field_local_double& resu, const int k_layer)
//{
//  valeur().compute_flux_y(resu, k_layer);
//}
//
//inline void Operateur_IJK_elem::compute_flux_z(IJK_Field_local_double& resu, const int k_layer)
//{
//  valeur().compute_flux_z(resu, k_layer);
//}
//
//inline void Operateur_IJK_elem::compute_(IJK_Field_double& dx, bool add)
//{
//  valeur().compute_(dx, add);
//}

#endif /* Operateur_IJK_elem_included */
