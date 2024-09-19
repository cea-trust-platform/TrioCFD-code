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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Modele_turbulence_scal_Fluctuation_Temperature_W.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/Common/Scal
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_turbulence_scal_Fluctuation_Temperature_W_included
#define Modele_turbulence_scal_Fluctuation_Temperature_W_included

#include <Equation_base.h>
#include <Transport_K_Eps_Bas_Reynolds.h>
#include <Modele_turbulence_scal_base.h>
#include <Transport_Fluctuation_Temperature_W.h>
#include <TRUST_Ref.h>
#include <Champ_Fonc.h>

class Transport_Fluctuation_Temperature_W;


class Modele_turbulence_scal_Fluctuation_Temperature_W :  public Modele_turbulence_scal_base
{
  Declare_instanciable(Modele_turbulence_scal_Fluctuation_Temperature_W);

public:
  void completer() override;
  int preparer_calcul() override;
  bool initTimeStep(double dt) override;
  void mettre_a_jour(double ) override;
  //virtual void associer_eqn(const Equation_base&);
  void associer_viscosite_turbulente(const Champ_Fonc& );
  inline virtual Champ_Inc_base& Fluctu_Temperature();
  inline virtual const Champ_Inc_base& Fluctu_Temperature() const;
  inline virtual Transport_Fluctuation_Temperature_W& equation_Fluctu();
  inline virtual const Transport_Fluctuation_Temperature_W& equation_Fluctu() const;

  int sauvegarder(Sortie& os) const override;
  int reprendre(Entree& is) override;
  virtual Champ_Fonc& calculer_diffusivite_turbulente();
  void imprimer(Sortie&) const override;

  void set_param(Param&) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;

  //////////////////////////////////////////////////////
  //Methode creer_champ pas codee a surcharger si necessaire
  //virtual void creer_champ(const Motcle& motlu);
  const Champ_base& get_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;
  /////////////////////////////////////////////////////

private :

//Entree& lire(const Motcle&, Entree&);
  REF(Transport_Fluctuation_Temperature_W) eqn_transport_Fluctu_Temp;



protected :
  OWN_PTR(Equation_base) eqn;
  REF(Champ_Fonc) la_viscosite_turbulente;
  // nous n'avons plus alpha_turb = visco_turb/Prdt_turb


};


//
//  Fonctions inline de la classe Modele_turbulence_hyd_K_Eps
//

inline Transport_Fluctuation_Temperature_W& Modele_turbulence_scal_Fluctuation_Temperature_W::equation_Fluctu()
{
  Cerr << "on est dans W11 " << finl;
  Cerr << " eqn =  " << eqn_transport_Fluctu_Temp.valeur();
  return eqn_transport_Fluctu_Temp.valeur();
}

inline const Transport_Fluctuation_Temperature_W& Modele_turbulence_scal_Fluctuation_Temperature_W::equation_Fluctu() const
{
  Cerr << "on est dans W11 " << finl;
  Cerr << " eqn =  " << eqn_transport_Fluctu_Temp.valeur();

  return eqn_transport_Fluctu_Temp.valeur();
}


inline const Champ_Inc_base& Modele_turbulence_scal_Fluctuation_Temperature_W::Fluctu_Temperature() const
{
  return eqn_transport_Fluctu_Temp->inconnue();
}

inline Champ_Inc_base& Modele_turbulence_scal_Fluctuation_Temperature_W::Fluctu_Temperature()
{
  return eqn_transport_Fluctu_Temp->inconnue();
}


#endif
