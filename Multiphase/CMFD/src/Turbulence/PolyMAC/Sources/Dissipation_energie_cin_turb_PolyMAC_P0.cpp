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
// File:        Dissipation_energie_cin_turb_PolyMAC_P0.cpp
// Directory:   $TRUST_ROOT/src/Turbulence/PolyMAC_P0/Sources
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Dissipation_energie_cin_turb_PolyMAC_P0.h>

Implemente_instanciable(Dissipation_energie_cin_turb_PolyMAC_P0,"Dissipation_energie_cin_turb_Elem_PolyMAC_P0|Terme_dissipation_energie_cinetique_turbulente_Elem_PolyMAC_P0", Source_Dissipation_energie_cin_turb);
// XD Terme_dissipation_energie_cinetique_turbulente source_base Terme_dissipation_energie_cinetique_turbulente -1 Dissipation source term used in the TKE equation
// XD attr beta_k floattant beta_k 1 Constant for the used model

// XD Production_echelle_temp_taux_diss_turb source_base Production_echelle_temp_taux_diss_turb -1 Production source term used in the tau and omega equations
// XD attr alpha_omega floattant alpha_omega 1 Constant for the used model

// XD Dissipation_echelle_temp_taux_diss_turb source_base Dissipation_echelle_temp_taux_diss_turb -1 Dissipation source term used in the tau and omega equations
// XD attr beta_omega floattant beta_omega 1 Constant for the used model

// XD Diffusion_croisee_echelle_temp_taux_diss_turb source_base Diffusion_croisee_echelle_temp_taux_diss_turb -1 Cross-diffusion source term used in the tau and omega equations
// XD attr sigma_d floattant sigma_d 1 Constant for the used model


Sortie& Dissipation_energie_cin_turb_PolyMAC_P0::printOn(Sortie& os) const { return Source_Dissipation_energie_cin_turb::printOn(os); }
Entree& Dissipation_energie_cin_turb_PolyMAC_P0::readOn(Entree& is) { return Source_Dissipation_energie_cin_turb::readOn(is); }