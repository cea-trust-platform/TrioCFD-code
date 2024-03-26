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
// File:        Transport_K_ou_Eps_base.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Equations/RANS
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Transport_K_ou_Eps_base_included
#define Transport_K_ou_Eps_base_included

#include <Equation_base.h>
#include <Mod_turb_hyd_RANS_Bicephale.h>
#include <TRUST_Ref.h>

class Champ_Inc;
class Milieu_base;
class Champ_Inc_base;

/*! @brief Classe Transport_K_ou_Eps_base Classe de base pour l'equation
 *
 *     de transport des modeles k_Epsilon dans une approche ou K et Epsilon sont traites par deux equations differentes.
 *
 */
class Transport_K_ou_Eps_base: public Equation_base
{

  Declare_base(Transport_K_ou_Eps_base);

public:

  void set_param(Param&) override;
  double calculer_pas_de_temps() const override;
  inline void associer_vitesse(const Champ_base& );
  void associer_milieu_base(const Milieu_base&) override;
  virtual void associer_modele_turbulence(const Mod_turb_hyd_RANS_Bicephale& )=0;
  const Milieu_base& milieu() const override ;
  Milieu_base& milieu() override ;
  void associer(const Equation_base&);
  void discretiser() override;
  virtual void discretiser_K_Eps(const Schema_Temps_base&, Domaine_dis&, Champ_Inc&) const;

  virtual int controler_variable();
  void valider_iteration() override;
  inline const Champ_Inc& inconnue() const override;
  inline Champ_Inc& inconnue() override;
  inline const Mod_turb_hyd_RANS_Bicephale& modele_turbulence() const;
  inline Mod_turb_hyd_RANS_Bicephale& modele_turbulence();

  const Champ_base& get_champ( const Motcle& nom ) const override;
//  void creer_champ( const Motcle& motlu );

  inline void transporte_K();
  inline void transporte_Eps();
  inline bool transporte_t_il_K() const;

protected:

  Champ_Inc le_champ_;
  Champ_Fonc residu_;

  REF(Milieu_base) le_fluide;
  REF(Champ_Inc_base) la_vitesse_transportante;
  REF(Mod_turb_hyd_RANS_Bicephale) mon_modele;

  bool transporte_K_;
};
/*! @brief Renvoie le champ inconnue de l'equation.
 *
 * Un champ vecteur contenant K ou epsilon.
 *
 * @return (Champ_Inc&) le champ inconnue de l'equation
 */
inline Champ_Inc& Transport_K_ou_Eps_base::inconnue()
{
  return le_champ_;
}


/*! @brief Renvoie le champ inconnue de l'equation.
 *
 * Un champ vecteur contenant K ou epsilon.
 *     (version const)
 *
 * @return (Champ_Inc&) le champ inconnue de l'equation
 */
inline const Champ_Inc& Transport_K_ou_Eps_base::inconnue() const
{
  return le_champ_;
}

/*! @brief Renvoie le modele de turbulence associe a l'equation.
 *
 * (version const)
 *
 * @return (Modele_turbulence_hyd_K_Eps&) le modele de turbulence associe a l'equation
 */
inline const Mod_turb_hyd_RANS_Bicephale& Transport_K_ou_Eps_base::modele_turbulence() const
{
  assert(mon_modele.non_nul());
  return mon_modele.valeur();
}


/*! @brief Renvoie le modele de turbulence associe a l'equation.
 *
 * @return (Modele_turbulence_hyd_K_Eps&) le modele de turbulence associe a l'equation
 */
inline Mod_turb_hyd_RANS_Bicephale& Transport_K_ou_Eps_base::modele_turbulence()
{
  assert(mon_modele.non_nul());
  return mon_modele.valeur();
}

inline void Transport_K_ou_Eps_base::associer_vitesse(const Champ_base& vit)
{
  la_vitesse_transportante = ref_cast(Champ_Inc_base,vit);
}

inline void Transport_K_ou_Eps_base::transporte_K()
{
  transporte_K_ = true;
}

inline void Transport_K_ou_Eps_base::transporte_Eps()
{
  transporte_K_ = false;
}

inline bool Transport_K_ou_Eps_base::transporte_t_il_K() const
{
  return transporte_K_;
}

#endif
