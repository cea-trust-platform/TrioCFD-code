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
// File:        Dissipation_echelle_temp_taux_diss_turb_VDF.cpp
// Directory:   $TRUST_ROOT/src/Turbulence/PolyMAC_P0/Sources
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Dissipation_echelle_temp_taux_diss_turb_VDF.h>

Implemente_instanciable(Dissipation_echelle_temp_taux_diss_turb_VDF,"Dissipation_echelle_temp_taux_diss_turb_VDF_P0_VDF", Source_Dissipation_echelle_temp_taux_diss_turb);
// X_D Dissipation_echelle_temp_taux_diss_turb_VDF_P0_VDF source_base Dissipation_echelle_temp_taux_diss_turb_VDF_P0_VDF 0 Source term which corresponds to the dissipation source term that appears in the transport equation for tau (in the k-tau turbulence model)

Sortie& Dissipation_echelle_temp_taux_diss_turb_VDF::printOn(Sortie& os) const {  return Source_Dissipation_echelle_temp_taux_diss_turb::printOn(os); }
Entree& Dissipation_echelle_temp_taux_diss_turb_VDF::readOn(Entree& is) { return Source_Dissipation_echelle_temp_taux_diss_turb::readOn(is); }
