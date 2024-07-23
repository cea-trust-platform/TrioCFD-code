/****************************************************************************
* Copyright (c) 2024, CEA
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

#include <Fluide_Incompressible.h>
#include <Milieu_base.h>
#include <Objet_U.h>
#include <Milieu.h>

class Fluide_Diphasique: public Milieu_base
{
  Declare_instanciable_sans_constructeur(Fluide_Diphasique);
public:

  Fluide_Diphasique()
  {
    indic_rayo_ = NONRAYO;
    formule_mu_ = "standard";
  }

  const Fluide_Incompressible& fluide_phase(int la_phase) const;
  double sigma() const;
  double chaleur_latente() const;
  int formule_mu() const;

  // Surcharge des methodes standard du milieu_base :
  void set_param(Param& param) override;
  void verifier_coherence_champs(int& err, Nom& message) override;
  int initialiser(const double temps) override;
  void mettre_a_jour(double temps) override;
  void discretiser(const Probleme_base& pb, const Discretisation_base& dis) override;

  // L'appel a ces methodes est invalide et genere une erreur
  const Champ& masse_volumique() const override { return invalid_<const Champ&>(__func__); }
  Champ& masse_volumique() override { return invalid_<Champ&>(__func__); }
  const Champ_Don& diffusivite() const override { return invalid_<const Champ_Don&>(__func__); }
  Champ_Don& diffusivite() override { return invalid_<Champ_Don&>(__func__); }
  const Champ_Don& conductivite() const override  { return invalid_<const Champ_Don&>(__func__); }
  Champ_Don& conductivite() override { return invalid_<Champ_Don&>(__func__); }
  const Champ_Don& capacite_calorifique() const override  { return invalid_<const Champ_Don&>(__func__); }
  Champ_Don& capacite_calorifique() override { return invalid_<Champ_Don&>(__func__); }
  const Champ_Don& beta_t() const override  { return invalid_<const Champ_Don&>(__func__); }
  Champ_Don& beta_t() override { return invalid_<Champ_Don&>(__func__); }

private:
  Milieu phase0_, phase1_;
  Champ_Don sigma_; // Tension de surface (J/m^2)
  Champ_Don chaleur_latente_; // Enthalpie de changement de phase h(phase1_) - h(phase0_) (J/kg/K)
  Motcle formule_mu_; // Formule utilisee pour le calcul de la moyenne de mu

  template <typename RETURN_TYPE>
  RETURN_TYPE invalid_(const char * nom_funct) const
  {
    Cerr << "Invalid call to the method Fluide_Diphasique::" << nom_funct << " !!!" << finl;
    throw;
  }
};

#endif /* Fluide_Diphasique_included */
