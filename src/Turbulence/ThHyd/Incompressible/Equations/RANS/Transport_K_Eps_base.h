/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
// File:        Transport_K_Eps_base.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Equations/RANS
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Transport_K_Eps_base_included
#define Transport_K_Eps_base_included

#include <Mod_turb_hyd_RANS_keps.h>
#include <Transport_2eq_base.h>
#include <TRUST_Ref.h>

class Champ_Inc;
class Milieu_base;
class Champ_Inc_base;

/*! @brief Classe Transport_K_Eps_base Classe de base pour les equations
 *
 *     de transport des modeles k_Epsilon.
 *
 */
class Transport_K_Eps_base: public Transport_2eq_base
{

  Declare_base_sans_constructeur(Transport_K_Eps_base);

public:

  Transport_K_Eps_base();
  virtual void associer_modele_turbulence(const Mod_turb_hyd_RANS_keps& )=0;
  void discretiser() override;
  void discretiser_K_Eps(const Schema_Temps_base&, Domaine_dis&, Champ_Inc&) const;

  int controler_K_Eps();
  void valider_iteration() override;
  inline const Champ_Inc& inconnue() const override;
  inline Champ_Inc& inconnue() override;
  inline const Mod_turb_hyd_RANS_keps& modele_turbulence() const;
  inline Mod_turb_hyd_RANS_keps& modele_turbulence();

  void get_position_cells(Nom&, int&);
  void get_position_faces(Nom&, int&);


protected:

  Champ_Inc le_champ_K_Eps;
  REF(Mod_turb_hyd_RANS_keps) mon_modele;

};

/*! @brief Renvoie le champ inconnue de l'equation.
 *
 * Un champ vecteur contenant K et epsilon.
 *
 * @return (Champ_Inc&) le champ inconnue de l'equation
 */
inline Champ_Inc& Transport_K_Eps_base::inconnue() { return le_champ_K_Eps; }


/*! @brief Renvoie le champ inconnue de l'equation.
 *
 * Un champ vecteur contenant K et epsilon.
 *     (version const)
 *
 * @return (Champ_Inc&) le champ inconnue de l'equation
 */
inline const Champ_Inc& Transport_K_Eps_base::inconnue() const { return le_champ_K_Eps; }

/*! @brief Renvoie le modele de turbulence associe a l'equation.
 *
 * (version const)
 *
 * @return (Modele_turbulence_hyd_K_Eps&) le modele de turbulence associe a l'equation
 */
inline const Mod_turb_hyd_RANS_keps& Transport_K_Eps_base::modele_turbulence() const
{
  assert(mon_modele.non_nul());
  return mon_modele.valeur();
}

/*! @brief Renvoie le modele de turbulence associe a l'equation.
 *
 * @return (Modele_turbulence_hyd_K_Eps&) le modele de turbulence associe a l'equation
 */
inline Mod_turb_hyd_RANS_keps& Transport_K_Eps_base::modele_turbulence()
{
  assert(mon_modele.non_nul());
  return mon_modele.valeur();
}

#endif
