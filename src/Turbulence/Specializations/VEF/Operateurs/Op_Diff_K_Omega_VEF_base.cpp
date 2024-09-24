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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Op_Diff_K_Omega_VEF_base.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Operateurs
//
//////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_K_Omega_VEF_base.h>
#include <Modele_turbulence_hyd_K_Omega.h>
#include <Champ_P1NC.h>
#include <Paroi_hyd_base_VEF.h>

Implemente_base_sans_constructeur(Op_Diff_K_Omega_VEF_base,
                                  "Op_Diff_K_Omega_VEF_base",
                                  Op_Diff_K_Omega_base);

Sortie& Op_Diff_K_Omega_VEF_base::printOn(Sortie& s) const
{
  return s << que_suis_je();
}

Entree& Op_Diff_K_Omega_VEF_base::readOn(Entree& s)
{
  return s;
}

void Op_Diff_K_Omega_VEF_base::completer()
{
  Operateur_base::completer();

  if (sub_type(Transport_K_Omega, mon_equation.valeur()))
    {
      const Transport_K_Omega& eqn_transport = ref_cast(Transport_K_Omega, mon_equation.valeur());
      const Modele_turbulence_hyd_K_Omega& mod_turb = ref_cast(Modele_turbulence_hyd_K_Omega,
                                                               eqn_transport.modele_turbulence());
      const Champ_Fonc_base& visc_turb = mod_turb.viscosite_turbulente();
      associer_diffusivite_turbulente(visc_turb);
      Op_Diff_K_Omega_VEF_base::associer_Pr_K_Omega(mod_turb.get_Prandtl_K(),
                                                    mod_turb.get_Prandtl_Omega());
    }
}

void  Op_Diff_K_Omega_VEF_base::associer_diffusivite(const Champ_base& diffu)
{
  diffusivite_ = diffu;
}
