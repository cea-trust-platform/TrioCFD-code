/****************************************************************************
* Copyright (c) 2019, CEA
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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Source_Transport_K_KEps_VDF_Elem.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Sources/UNUSED
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Source_Transport_K_KEps_VDF_Elem_included
#define Source_Transport_K_KEps_VDF_Elem_included

#include <Source_Transport_VDF_Elem_base.h>
#include <Ref_Transport_K_KEps.h>

//.DESCRIPTION class Source_Transport_K_KEps_VDF_Elem
// Cette classe represente le terme source qui figure dans l'equation de transport du couple (k,eps) avec le modele a deux couches et dans
// le cas ou les equations de Navier-Stokes ne sont pas couplees a la thermique ou a l'equation de convection-diffusion d'une concentration.
class Source_Transport_K_KEps_VDF_Elem : public Source_Transport_VDF_Elem_base
{
  Declare_instanciable_sans_constructeur(Source_Transport_K_KEps_VDF_Elem);
public:
  Source_Transport_K_KEps_VDF_Elem(double cte1 = C1__, double cte2 = C2__ ) : Source_Transport_VDF_Elem_base(cte1,cte2) { }
  DoubleTab& ajouter(DoubleTab& ) const override;

protected:
  REF(Transport_K_KEps)  mon_eq_transport_K_Eps;
  void associer_pb(const Probleme_base& pb) override;
};

#endif /* Source_Transport_K_KEps_VDF_Elem_included */
