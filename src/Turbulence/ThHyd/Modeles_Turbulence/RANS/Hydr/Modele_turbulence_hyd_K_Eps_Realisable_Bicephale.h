/****************************************************************************
* Copyright (c) 2019, CEA
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

#ifndef Modele_turbulence_hyd_K_Eps_Realisable_Bicephale_included
#define Modele_turbulence_hyd_K_Eps_Realisable_Bicephale_included

#include <Modele_turbulence_hyd_RANS_Bicephale_base.h>
#include <Transport_K_ou_Eps_Realisable.h>
#include <Modele_Fonc_Realisable_base.h>



/*! @brief class Modele_turbulence_hyd_K_Eps_Realisable_Bicephale
 *
 *  Decrire ici la classe Modele_turbulence_hyd_K_Eps_Realisable_Bicephale
 *  Cette classe represente le modele de turbulence (k,eps) realisable ou les 2 equations de k et eps sont gerees separement du point de vue informatique.
 *
 */
class Modele_turbulence_hyd_K_Eps_Realisable_Bicephale: public Modele_turbulence_hyd_RANS_Bicephale_base, public Modele_turbulence_hyd_RANS_Gen<Modele_turbulence_hyd_K_Eps_Realisable_Bicephale>
{
  Declare_instanciable(Modele_turbulence_hyd_K_Eps_Realisable_Bicephale);
public:
  void set_param(Param& param) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  int preparer_calcul() override;
  void verifie_loi_paroi() override;
  bool initTimeStep(double dt) override;
  void mettre_a_jour(double) override;
  const Equation_base& equation_k_eps(int) const override;

  virtual inline Champ_Inc& K() { return eqn_transport_K_Rea_.inconnue(); }
  virtual inline const Champ_Inc& K() const { return eqn_transport_K_Rea_.inconnue(); }
  virtual inline Champ_Inc& Eps() { return eqn_transport_Eps_Rea_.inconnue(); }
  virtual inline const Champ_Inc& Eps() const { return eqn_transport_Eps_Rea_.inconnue(); }

  inline Transport_K_ou_Eps_base& eqn_transp_K() override { return eqn_transport_K_Rea_; }
  inline const Transport_K_ou_Eps_base& eqn_transp_K() const override { return eqn_transport_K_Rea_; }
  inline Transport_K_ou_Eps_base& eqn_transp_Eps() override { return eqn_transport_Eps_Rea_; }
  inline const Transport_K_ou_Eps_base& eqn_transp_Eps() const override { return eqn_transport_Eps_Rea_; }

  inline OWN_PTR(Modele_Fonc_Realisable_base)& associe_modele_fonction() { return mon_modele_fonc_; }
  inline const OWN_PTR(Modele_Fonc_Realisable_base)& associe_modele_fonction() const { return mon_modele_fonc_; }
  inline OWN_PTR(Modele_Fonc_Realisable_base)& get_modele_fonction() { return mon_modele_fonc_; }
  inline const OWN_PTR(Modele_Fonc_Realisable_base)& get_modele_fonction() const { return mon_modele_fonc_; }

  const Champ_base& get_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom, Option opt = NONE) const override;

  virtual Champ_Fonc& calculer_viscosite_turbulente(double temps);
  void controler()
  {
    eqn_transport_K_Rea_.controler_variable();
    eqn_transport_Eps_Rea_.controler_variable();
  }

private:
  OWN_PTR(Modele_Fonc_Realisable_base) mon_modele_fonc_;
  Transport_K_ou_Eps_Realisable eqn_transport_K_Rea_, eqn_transport_Eps_Rea_;
  void fill_turbulent_viscosity_tab(const int , const DoubleTab&, const DoubleTab&, DoubleTab& );
};

#endif /* Modele_turbulence_hyd_K_Eps_Realisable_Bicephale_included */
