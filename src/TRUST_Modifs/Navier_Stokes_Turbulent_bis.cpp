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

#include <Navier_Stokes_Turbulent.h>
#include <Modele_turbulence_hyd_K_Eps.h>
#include <Modele_turbulence_hyd_K_Eps_Realisable.h>
#include <Modele_turbulence_hyd_K_Eps_Bicephale.h>
#include <Modele_turbulence_hyd_K_Omega.h>

void Navier_Stokes_Turbulent::creer_champ(const Motcle& motlu)
{
  Navier_Stokes_std::creer_champ(motlu);

  if (le_modele_turbulence.non_nul())
    le_modele_turbulence->creer_champ(motlu);

  // to create k_eps_residu field
  if(le_modele_turbulence.non_nul())
    {
      if (sub_type(Modele_turbulence_hyd_K_Eps, le_modele_turbulence.valeur()))
        ref_cast(Modele_turbulence_hyd_K_Eps,
                 le_modele_turbulence.valeur()).eqn_transp_K_Eps().creer_champ(motlu);
      else if (sub_type(Modele_turbulence_hyd_K_Eps_Realisable, le_modele_turbulence.valeur()))
        ref_cast(Modele_turbulence_hyd_K_Eps_Realisable,
                 le_modele_turbulence.valeur()).eqn_transp_K_Eps().creer_champ(motlu);
      else if (sub_type(Modele_turbulence_hyd_K_Omega, le_modele_turbulence.valeur()))
        ref_cast(Modele_turbulence_hyd_K_Omega,
                 le_modele_turbulence.valeur()).eqn_transp_K_Omega().creer_champ(motlu);
    }
}

// Impression du residu dans fic (generalement dt_ev)
// Cette methode peut etre surchargee par des equations
// imprimant des residus particuliers (K-Eps, Concentrations,...)
void Navier_Stokes_Turbulent::imprime_residu(SFichier& fic)
{
  Equation_base::imprime_residu(fic);
  // Si K-Eps, on imprime
  if (sub_type(Modele_turbulence_hyd_K_Eps,le_modele_turbulence.valeur()))
    ref_cast(Modele_turbulence_hyd_K_Eps,le_modele_turbulence.valeur()).eqn_transp_K_Eps().imprime_residu(fic);
  else if (sub_type(Modele_turbulence_hyd_K_Eps_Realisable,le_modele_turbulence.valeur()))
    ref_cast(Modele_turbulence_hyd_K_Eps_Realisable,le_modele_turbulence.valeur()).eqn_transp_K_Eps().imprime_residu(fic);
  else if (sub_type(Modele_turbulence_hyd_K_Eps_Bicephale,le_modele_turbulence.valeur()))
    {
      ref_cast(Modele_turbulence_hyd_K_Eps_Bicephale,le_modele_turbulence.valeur()).eqn_transp_K().imprime_residu(fic);
      ref_cast(Modele_turbulence_hyd_K_Eps_Bicephale,le_modele_turbulence.valeur()).eqn_transp_Eps().imprime_residu(fic);
    }
  else if (sub_type(Modele_turbulence_hyd_K_Omega,le_modele_turbulence.valeur()))
    ref_cast(Modele_turbulence_hyd_K_Omega,le_modele_turbulence.valeur()).eqn_transp_K_Omega().imprime_residu(fic);
}

// Retourne l'expression du residu (de meme peut etre surcharge)
Nom Navier_Stokes_Turbulent::expression_residu()
{
  Nom tmp=Equation_base::expression_residu();
  if (sub_type(Modele_turbulence_hyd_K_Eps,le_modele_turbulence.valeur()))
    tmp+=ref_cast(Modele_turbulence_hyd_K_Eps,le_modele_turbulence.valeur()).eqn_transp_K_Eps().expression_residu();
  else if (sub_type(Modele_turbulence_hyd_K_Eps_Realisable,le_modele_turbulence.valeur()))
    tmp+=ref_cast(Modele_turbulence_hyd_K_Eps_Realisable,le_modele_turbulence.valeur()).eqn_transp_K_Eps().expression_residu();
  else if (sub_type(Modele_turbulence_hyd_K_Omega,le_modele_turbulence.valeur()))
    tmp+=ref_cast(Modele_turbulence_hyd_K_Omega,le_modele_turbulence.valeur()).eqn_transp_K_Omega().expression_residu();
  else if (sub_type(Modele_turbulence_hyd_K_Eps_Bicephale,le_modele_turbulence.valeur()))
    {
      tmp+=ref_cast(Modele_turbulence_hyd_K_Eps_Bicephale,le_modele_turbulence.valeur()).eqn_transp_K().expression_residu();
      tmp+=ref_cast(Modele_turbulence_hyd_K_Eps_Bicephale,le_modele_turbulence.valeur()).eqn_transp_Eps().expression_residu();
    }
  return tmp;
}

