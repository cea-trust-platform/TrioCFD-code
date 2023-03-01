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
// File:        Convection_diffusion_turbulence_multiphase.h
// Directory:   $TRUST_ROOT/src/Turbulence/Equations
// Version:     /main/20
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Convection_diffusion_turbulence_multiphase_included
#define Convection_diffusion_turbulence_multiphase_included

#include <Convection_Diffusion_std.h>
#include <Operateur_Grad.h>
#include <Fluide_base.h>
#include <TRUST_Ref.h>

/*! @brief classe Convection_diffusion_turbulence_multiphase Equation de transport des quantites turbulentes (k, omega, epsilon, tau)
 *
 * @sa Conv_Diffusion_std Convection_Diffusion_Temperature
 */
class Convection_diffusion_turbulence_multiphase : public Convection_Diffusion_std
{
  Declare_base(Convection_diffusion_turbulence_multiphase);

public :

  void completer() override;
  inline const Champ_Inc& inconnue() const override;
  inline Champ_Inc& inconnue() override;
  void associer_milieu_base(const Milieu_base& ) override;
  const Milieu_base& milieu() const override;
  Milieu_base& milieu() override;
  int impr(Sortie& os) const override;
  const Motcle& domaine_application() const override;
  int positive_unkown() override {return 1;};
  inline const Operateur_Grad& operateur_gradient_inconnue() const { return Op_Grad_;}

protected :

  Champ_Inc l_inco_ch;
  REF(Fluide_base) le_fluide;
  Operateur_Grad Op_Grad_; // Pour calculer le gradient en VDF

};




/*! @brief Renvoie le champ inconnue representant l'inconnue (k, omega, epsilon, tau) (version const)
 *
 * @return (Champ_Inc&) le champ inconnue representant la temperature (GP) ou l'enthalpie (GR)
 */
inline const Champ_Inc& Convection_diffusion_turbulence_multiphase::inconnue() const
{
  return l_inco_ch;
}


/*! @brief Renvoie le champ inconnue representant l'inconnue (k, omega, epsilon, tau)
 *
 * @return (Champ_Inc&) le champ inconnue representant la temperature (GP) ou l'enthalpie (GR)
 */
inline Champ_Inc& Convection_diffusion_turbulence_multiphase::inconnue()
{
  return l_inco_ch;
}

#endif
