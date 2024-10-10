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
// File:        Energie_cinetique_turbulente_WIT.h
// Directory:   $TRUST_ROOT/src/Turbulence/Equations
// Version:     /main/20
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Energie_cinetique_turbulente_WIT_included
#define Energie_cinetique_turbulente_WIT_included

#include <Convection_Diffusion_std.h>
#include <Fluide_base.h>
#include <TRUST_Ref.h>

/*! @brief classe Energie_cinetique_turbulente_WIT Equation de transport d'une energie cinetique turbulente WIT pour la turbulence avec bulles
 *
 * @sa Conv_Diffusion_std Convection_Diffusion_Temperature
 */
class Energie_cinetique_turbulente_WIT : public Convection_Diffusion_std
{
  Declare_instanciable_sans_constructeur(Energie_cinetique_turbulente_WIT);

public :

  Energie_cinetique_turbulente_WIT();

  void associer_fluide(const Fluide_base& );
  inline const Champ_Inc_base& inconnue() const override;
  inline Champ_Inc_base& inconnue() override;
  void discretiser() override;
  const Milieu_base& milieu() const override;
  Milieu_base& milieu() override;
  void associer_milieu_base(const Milieu_base& ) override;
  int impr(Sortie& os) const override;
  const Champ_Don_base& diffusivite_pour_transport() const override;
  const Champ_base& diffusivite_pour_pas_de_temps() const override;

  const Motcle& domaine_application() const override;

  /* champ convecte : alpha (si Pb_Multiphase) * rho * k */
  static void calculer_alpha_rho_k_WIT(const Objet_U& obj, DoubleTab& val, DoubleTab& bval, tabs_t& deriv);
  virtual std::pair<std::string, fonc_calc_t> get_fonc_champ_conserve() const override
  {
    return { "alpha_rho_k_WIT", calculer_alpha_rho_k_WIT };
  }

  virtual int positive_unkown() override {return 1;};

protected :

  OWN_PTR(Champ_Inc_base) l_inco_ch;
  OBS_PTR(Fluide_base) le_fluide;
};




/*! @brief Renvoie le champ inconnue representant l'inconnue (T ou H) (version const)
 *
 * @return (Champ_Inc_base&) le champ inconnue representant la temperature (GP) ou l'enthalpie (GR)
 */
inline const Champ_Inc_base& Energie_cinetique_turbulente_WIT::inconnue() const
{
  return l_inco_ch;
}


/*! @brief Renvoie le champ inconnue representant l'inconnue (T ou H)
 *
 * @return (Champ_Inc_base&) le champ inconnue representant la temperature (GP) ou l'enthalpie (GR)
 */
inline Champ_Inc_base& Energie_cinetique_turbulente_WIT::inconnue()
{
  return l_inco_ch;
}

#endif
