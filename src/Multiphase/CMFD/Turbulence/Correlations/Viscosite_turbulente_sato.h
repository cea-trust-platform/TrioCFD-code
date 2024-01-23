/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Viscosite_turbulente_sato.h
// Directory:   $TRUST_ROOT/src/Turbulence/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Viscosite_turbulente_sato_included
#define Viscosite_turbulente_sato_included
#include <TRUSTTab.h>
#include <Viscosite_turbulente_base.h>

/*! @brief classe Viscosite_turbulente_sato Viscosite turbulente induite par les bulles : nu_b = coef_sato*alpha*d_B*norm(u_R)
 *
 *     (Taux de vide*diametre des bulles*norme de la vitesse relative)
 *
 *
 */
class Viscosite_turbulente_sato : public Viscosite_turbulente_base
{
  Declare_instanciable(Viscosite_turbulente_sato);
public:
  void eddy_viscosity(DoubleTab& nu_t) const override;
  void reynolds_stress(DoubleTab& R_ij) const override;
  void k_over_eps(DoubleTab& k_sur_eps) const override;
  void eps(DoubleTab& eps) const override;

private:
  double coef_sato = 0.6; // coefficient du modèle ajusté par Sato
  int n_l = -1; // phase liquide
};

#endif
