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
// File:        Echelle_temporelle_turbulente.h
// Directory:   $TRUST_ROOT/src/Turbulence/Equations
// Version:     /main/20
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Echelle_temporelle_turbulente_included
#define Echelle_temporelle_turbulente_included

#include <Convection_Diffusion_std.h>
#include <Fluide_base.h>
#include <Ref_Fluide_base.h>

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//     classe Echelle_temporelle_turbulente
//     Equation de transport de l'echelle temporelle turbulente (modele k-tau)
// .SECTION voir aussi
//     Conv_Diffusion_std Convection_Diffusion_Temperature
//////////////////////////////////////////////////////////////////////////////
class Echelle_temporelle_turbulente : public Convection_Diffusion_std
{
  Declare_instanciable_sans_constructeur(Echelle_temporelle_turbulente);

public :

  Echelle_temporelle_turbulente();

  void associer_fluide(const Fluide_base& );
  inline const Champ_Inc& inconnue() const;
  inline Champ_Inc& inconnue();
  void discretiser();
  const Milieu_base& milieu() const;
  const Fluide_base& fluide() const;
  Fluide_base& fluide();
  Milieu_base& milieu();
  void associer_milieu_base(const Milieu_base& );
  virtual int impr(Sortie& os) const;
  const Champ_Don& diffusivite_pour_transport();
  const Champ_base& diffusivite_pour_pas_de_temps();

  virtual const Motcle& domaine_application() const;

  /* champ convecte : alpha (si Pb_Multiphase) * rho * k */
  static void calculer_alpha_rho_tau(const Objet_U& obj, DoubleTab& val, DoubleTab& bval, tabs_t& deriv);
  virtual std::pair<std::string, fonc_calc_t> get_fonc_champ_conserve() const
  {
    return { "alpha_rho_tau", calculer_alpha_rho_tau };
  }


protected :

  Champ_Inc l_inco_ch;
  REF(Fluide_base) le_fluide;
};




// Description:
//    Renvoie le champ inconnue representant l'inconnue (T ou H)
//    (version const)
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Champ_Inc&
//    Signification: le champ inconnue representant la temperature (GP)
//      ou l'enthalpie (GR)
//    Contraintes: reference constante
// Exception:
// Effets de bord:
// Postcondition:
inline const Champ_Inc& Echelle_temporelle_turbulente::inconnue() const
{
  return l_inco_ch;
}


// Description:
//    Renvoie le champ inconnue representant l'inconnue (T ou H)
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Champ_Inc&
//    Signification: le champ inconnue representant la temperature (GP)
//      ou l'enthalpie (GR)
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
inline Champ_Inc& Echelle_temporelle_turbulente::inconnue()
{
  return l_inco_ch;
}

#endif
