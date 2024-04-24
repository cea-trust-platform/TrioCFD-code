/****************************************************************************
* Copyright (c) 2017, CEA
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
// File:        Modele_turbulence_hyd_RANS_K_Eps_base.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Modele_turbulence_hyd_RANS_K_Eps_base_included
#define Modele_turbulence_hyd_RANS_K_Eps_base_included

#include <Modele_turbulence_hyd_2_eq_base.h>
#include <Modele_turbulence_hyd_RANS_Gen.h>
#include <Modele_Fonc_Bas_Reynolds.h>

class Transport_K_Eps_base;

/*! @brief Classe Modele_turbulence_hyd_RANS_K_Eps_base Classe de base des modeles de type RANS_keps
 *
 * @sa Modele_turbulence_hyd_base
 */
class Modele_turbulence_hyd_RANS_K_Eps_base : public Modele_turbulence_hyd_2_eq_base
{
  Declare_base(Modele_turbulence_hyd_RANS_K_Eps_base);
public:

  void set_param(Param& param) override;
  void completer() override;
  int sauvegarder(Sortie& os) const override;
  int reprendre(Entree& is) override;

  const Champ_base& get_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;

  virtual Transport_K_Eps_base& eqn_transp_K_Eps()=0;
  virtual const Transport_K_Eps_base& eqn_transp_K_Eps() const=0;
  virtual const Equation_base& equation_k_eps(int) const=0 ;

  inline Modele_Fonc_Bas_Reynolds& associe_modele_fonction() { return mon_modele_fonc_; }
  inline const Modele_Fonc_Bas_Reynolds& associe_modele_fonction() const { return mon_modele_fonc_; }

  bool calcul_tenseur_Re(const DoubleTab& nu_turb, const DoubleTab& grad, DoubleTab& Re) const override
  {
    if (associe_modele_fonction().non_nul() && associe_modele_fonction().Calcul_is_Reynolds_stress_isotrope()==0)
      return associe_modele_fonction().valeur().calcul_tenseur_Re(nu_turb, grad, Re);
    else
      return false;
  }

protected:
  Modele_Fonc_Bas_Reynolds mon_modele_fonc_;
};

#endif /* Modele_turbulence_hyd_RANS_K_Eps_base_included */
