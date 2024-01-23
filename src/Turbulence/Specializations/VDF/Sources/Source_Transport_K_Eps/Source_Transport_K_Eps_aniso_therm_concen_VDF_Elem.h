/****************************************************************************
* Copyright (c) 2022, CEA
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
// File      : Source_Transport_K_Eps_aniso_therm_concen_VDF_Elem.h
// Directory : $TURBULENCE_ROOT/src/Specializations/VDF/Sources
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Source_Transport_K_Eps_aniso_therm_concen_VDF_Elem_included
#define Source_Transport_K_Eps_aniso_therm_concen_VDF_Elem_included

#include <Source_Transport_K_Eps_VDF_Elem.h>

// CLASS: Source_Transport_K_Eps_aniso_therm_concen_VDF_Elem
// Cette classe represente le terme source qui figure dans l'equation de transport du couple (k,eps) dans le cas ou les equations de
// Navier_Stokes sont couplees a l'equation de convection diffusion d'une concentration et a l'equation de la thermique. Les champs beta_t et beta_c sont uniformes
class Source_Transport_K_Eps_aniso_therm_concen_VDF_Elem : public Source_Transport_K_Eps_VDF_Elem
{
  Declare_instanciable_sans_constructeur(Source_Transport_K_Eps_aniso_therm_concen_VDF_Elem);
public:
  Source_Transport_K_Eps_aniso_therm_concen_VDF_Elem(double cte1 = C1__, double cte2 = C2__, double cte3 = C3__) : Source_Transport_K_Eps_VDF_Elem(cte1,cte2) { C3 = cte3; }
  void associer_pb(const Probleme_base& ) override;
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const override;

private:
  void fill_resu_anisotherme_concen(const DoubleVect& , const DoubleVect& , const DoubleVect& , const DoubleVect& , DoubleTab& ) const override;
};

#endif /* Source_Transport_K_Eps_aniso_therm_concen_VDF_Elem_included */
