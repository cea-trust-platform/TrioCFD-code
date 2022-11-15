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

#include <Problemes_rayo.h>

// XD Pb_Rayo_Conduction Pb_Conduction Pb_Rayo_Conduction -1 Resolution of the heat equation with rayonnement.
// XD Pb_Rayo_Hydraulique pb_hydraulique Pb_Rayo_Hydraulique -1 Resolution of the Navier-Stokes equations with rayonnement.
// XD Pb_Rayo_Thermohydraulique pb_thermohydraulique Pb_Rayo_Thermohydraulique -1 Resolution of pb_thermohydraulique with rayonnement.
// XD Pb_Rayo_Hydraulique_Turbulent pb_hydraulique_turbulent Pb_Rayo_Hydraulique_Turbulent -1 Resolution of pb_hydraulique_turbulent with rayonnement.
// XD Pb_Rayo_Thermohydraulique_Turbulent pb_thermohydraulique_turbulent Pb_Rayo_Thermohydraulique_Turbulent -1 Resolution of pb_thermohydraulique_turbulent with rayonnement.

Implemente_instanciable(Pb_Rayo_Conduction,"Pb_Rayo_Conduction",Pb_Conduction);
Implemente_instanciable(Pb_Rayo_Hydraulique,"Pb_Rayo_Hydraulique",Pb_Hydraulique);
Implemente_instanciable(Pb_Rayo_Thermohydraulique,"Pb_Rayo_Thermohydraulique",Pb_Thermohydraulique);
Implemente_instanciable(Pb_Rayo_Hydraulique_Turbulent,"Pb_Rayo_Hydraulique_Turbulent",Pb_Hydraulique_Turbulent);
Implemente_instanciable(Pb_Rayo_Thermohydraulique_Turbulent,"Pb_Rayo_Thermohydraulique_Turbulent",Pb_Thermohydraulique_Turbulent);

Sortie& Pb_Rayo_Conduction::printOn(Sortie& os) const { return Pb_Conduction::printOn(os); }
Entree& Pb_Rayo_Conduction::readOn(Entree& is) { return Pb_Conduction::readOn(is); }

Sortie& Pb_Rayo_Hydraulique::printOn(Sortie& os) const { return Pb_Hydraulique::printOn(os); }
Entree& Pb_Rayo_Hydraulique::readOn(Entree& is) { return Pb_Hydraulique::readOn(is); }

Sortie& Pb_Rayo_Thermohydraulique::printOn(Sortie& os) const { return Pb_Thermohydraulique::printOn(os); }
Entree& Pb_Rayo_Thermohydraulique::readOn(Entree& is) { return Pb_Thermohydraulique::readOn(is); }

Sortie& Pb_Rayo_Hydraulique_Turbulent::printOn(Sortie& os) const { return Pb_Hydraulique_Turbulent::printOn(os); }
Entree& Pb_Rayo_Hydraulique_Turbulent::readOn(Entree& is) { return Pb_Hydraulique_Turbulent::readOn(is); }

Sortie& Pb_Rayo_Thermohydraulique_Turbulent::printOn(Sortie& os) const { return Pb_Thermohydraulique_Turbulent::printOn(os); }
Entree& Pb_Rayo_Thermohydraulique_Turbulent::readOn(Entree& is) { return Pb_Thermohydraulique_Turbulent::readOn(is); }
