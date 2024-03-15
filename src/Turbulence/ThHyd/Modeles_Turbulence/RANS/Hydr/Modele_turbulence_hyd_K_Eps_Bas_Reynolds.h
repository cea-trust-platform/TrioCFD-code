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
/*! @brief class Modele_turbulence_hyd_K_Eps_Bas_Reynolds
 *
 *  Decrire ici la classe Modele_turbulence_hyd_K_Eps_Bas_Reynolds
 *  Cette classe represente le modele de turbulence (k,eps) pour de
 *  faible nombre de Reynolds.
 *
 *
 *
 */
#ifndef Modele_turbulence_hyd_K_Eps_Bas_Reynolds_included
#define Modele_turbulence_hyd_K_Eps_Bas_Reynolds_included

#include <Modele_turbulence_hyd_K_Eps.h>
#include <Transport_K_Eps_Bas_Reynolds.h>
class Domaine_dis;
class Domaine_Cl_dis;


class Modele_turbulence_hyd_K_Eps_Bas_Reynolds : public  Modele_turbulence_hyd_RANS_keps_base
{
  Declare_instanciable(Modele_turbulence_hyd_K_Eps_Bas_Reynolds);

public:

  void set_param(Param& param) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  int preparer_calcul() override;
  bool initTimeStep(double dt) override;
  void mettre_a_jour(double ) override;
  virtual inline Champ_Inc& K_Eps();
  virtual inline const Champ_Inc& K_Eps() const;
  void completer() override;

  inline int nombre_d_equations() const override;
  inline Transport_K_Eps_base& eqn_transp_K_Eps() override;
  inline const Transport_K_Eps_base& eqn_transp_K_Eps() const override;
  const Equation_base& equation_k_eps(int) const override ;

  //Methodes de l interface des champs postraitables
  /////////////////////////////////////////////////////
  //Methode creer_champ pas codde a surcharger
  //virtual void creer_champ(const Motcle& motlu);
  const Champ_base& get_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;
  /////////////////////////////////////////////////////

  void imprimer(Sortie&) const override;
private :
  Transport_K_Eps_Bas_Reynolds  eqn_transport_K_Eps_Bas_Re;
  virtual Champ_Fonc& calculer_viscosite_turbulente(double temps);

};

//
//  Fonctions inline de la classe Modele_turbulence_hyd_K_Eps
//

inline const Transport_K_Eps_base& Modele_turbulence_hyd_K_Eps_Bas_Reynolds::eqn_transp_K_Eps() const
{
  return eqn_transport_K_Eps_Bas_Re;
}

inline  Transport_K_Eps_base& Modele_turbulence_hyd_K_Eps_Bas_Reynolds::eqn_transp_K_Eps()
{
  return eqn_transport_K_Eps_Bas_Re;
}

inline const Champ_Inc& Modele_turbulence_hyd_K_Eps_Bas_Reynolds::K_Eps() const
{
  return eqn_transport_K_Eps_Bas_Re.inconnue();
}

inline Champ_Inc& Modele_turbulence_hyd_K_Eps_Bas_Reynolds::K_Eps()
{
  return eqn_transport_K_Eps_Bas_Re.inconnue();
}

inline int Modele_turbulence_hyd_K_Eps_Bas_Reynolds::nombre_d_equations() const
{
  return 1;
}

#endif
