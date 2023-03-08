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
// File:        Taux_dissipation_turbulent.h
// Directory:   $TRUST_ROOT/src/Turbulence/Equations
// Version:     /main/20
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Taux_dissipation_turbulent_included
#define Taux_dissipation_turbulent_included

#include <Convection_diffusion_turbulence_multiphase.h>
#include <Operateur_Grad.h>
#include <Fluide_base.h>
#include <TRUST_Ref.h>

/*! @brief classe Taux_dissipation_turbulent Equation de transport du taux de dissipation turbulen (modele k-omega)
 *
 * @sa Conv_Diffusion_std Convection_Diffusion_Temperature Convection_diffusion_turbulence_multiphase
 */
class Taux_dissipation_turbulent : public Convection_diffusion_turbulence_multiphase
{
  Declare_instanciable_sans_constructeur(Taux_dissipation_turbulent);

public :

  Taux_dissipation_turbulent();

  void discretiser() override;

  const Champ_Don& diffusivite_pour_transport() const override;
  const Champ_base& diffusivite_pour_pas_de_temps() const override;

  /* champ convecte : alpha (si Pb_Multiphase) * rho * k */
  static void calculer_omega(const Objet_U& obj, DoubleTab& val, DoubleTab& bval, tabs_t& deriv);
  virtual std::pair<std::string, fonc_calc_t> get_fonc_champ_conserve() const override
  {
    return { "omega", calculer_omega };
  }
};

#endif
