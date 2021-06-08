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
// File:        Fluide_Diphasique.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////
#ifndef Fluide_Diphasique_included
#define Fluide_Diphasique_included
#include <Objet_U.h>
#include <Fluide_Incompressible.h>
#include <Milieu_base.h>

class Entree;
class Motcle;

class Fluide_Diphasique : public Milieu_base
{
  Declare_instanciable_sans_constructeur(Fluide_Diphasique);
public:
  Fluide_Diphasique();
  // Renvoie le fluide de la phase 0 ou 1.
  const Fluide_Incompressible& fluide_phase(int la_phase) const;
  double sigma() const;
  double chaleur_latente() const;
  int formule_mu() const;
  // Methode utilisee pour le calcul du mu du domaine

  // Surcharge des methodes standard du milieu_base :
  void set_param(Param& param);
  void verifier_coherence_champs(int& err,Nom& message);
  int lire_motcle_non_standard(const Motcle&, Entree&);
  int initialiser(const double& temps);
  void mettre_a_jour(double temps);
  virtual void discretiser(const Probleme_base& pb, const  Discretisation_base& dis);

  // L'appel a ces methodes est invalide et genere une erreur
  const Champ_base& masse_volumique() const;
  Champ_base&       masse_volumique();
  const Champ_Don& diffusivite() const;
  Champ_Don&       diffusivite();
  const Champ_Don& conductivite() const;
  Champ_Don&       conductivite();
  const Champ_Don& capacite_calorifique() const;
  Champ_Don&       capacite_calorifique();
  const Champ_Don& beta_t() const;
  Champ_Don&       beta_t();


protected:

private:

  Fluide_Incompressible phase0_;
  Fluide_Incompressible phase1_;
  // Tension de surface (J/m^2)
  Champ_Don sigma_;
  // Enthalpie de changement de phase h(phase1_) - h(phase0_) (J/kg/K)
  Champ_Don chaleur_latente_;
  // Formule utilisee pour le calcul de la moyenne de mu
  Motcle formule_mu_;
};
#endif
