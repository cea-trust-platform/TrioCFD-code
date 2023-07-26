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
// File:        Source_Diffusion_croisee_echelle_temp_taux_diss_turb.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Sources
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Diffusion_croisee_echelle_temp_taux_diss_turb.h>

#include <Domaine_Cl_PolyMAC.h>
#include <Pb_Multiphase.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <Domaine_VF.h>

Implemente_base(Source_Diffusion_croisee_echelle_temp_taux_diss_turb,"Source_Diffusion_croisee_echelle_temp_taux_diss_turb", Sources_Multiphase_base);

Sortie& Source_Diffusion_croisee_echelle_temp_taux_diss_turb::printOn(Sortie& os) const
{
  return os;
}

Entree& Source_Diffusion_croisee_echelle_temp_taux_diss_turb::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("sigma_d", &sigma_d);
  param.lire_avec_accolades_depuis(is);
  return is;
}