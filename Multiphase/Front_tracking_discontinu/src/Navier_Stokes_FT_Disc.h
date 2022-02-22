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
// File:        Navier_Stokes_FT_Disc.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/22
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Navier_Stokes_FT_Disc_included
#define Navier_Stokes_FT_Disc_included

#include <Navier_Stokes_Turbulent.h>
#include <Ref_Probleme_FT_Disc_gen.h>
#include <Champ_Don.h>

class Navier_Stokes_FT_Disc_interne;
class Maillage_FT_Disc;
class Fluide_Diphasique;

class Navier_Stokes_FT_Disc : public Navier_Stokes_Turbulent
{
  Declare_instanciable_sans_constructeur(Navier_Stokes_FT_Disc);
public:

  Navier_Stokes_FT_Disc();
  //
  // Methodes surchargees
  //
  void set_param(Param& titi) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  const Milieu_base& milieu() const override;
  Milieu_base&        milieu() override;
  void                associer_pb_base(const Probleme_base& probleme) override;
  void                discretiser() override;
  int              preparer_calcul() override;
  void                preparer_pas_de_temps();
  void mettre_a_jour(double temps) override;
  void                calculer_la_pression_en_pa() override;
  DoubleTab&          derivee_en_temps_inco(DoubleTab& vpoint) override;
  void                projeter() override;
  virtual const Champ_base& calculer_div_normale_interface();
  void calculer_delta_u_interface(Champ_base& u0, int phase_pilote, int ordre);
  const Champ_Don& diffusivite_pour_transport() const override;

  virtual const Champ_base * get_delta_vitesse_interface() const;
  virtual const Fluide_Diphasique&     fluide_diphasique() const;

  int is_terme_gravite_rhog() const;
  const Champ_Fonc& champ_rho_faces() const;

  virtual void calculer_dI_dt(DoubleVect& dI_dt) const;
  const int& get_is_penalized() const;

protected:
  // Methode surchargee de Navier_Stokes_std :
  void discretiser_assembleur_pression() override;
  void associer_milieu_base(const Milieu_base& fluide) override;

  // Nouvelles methodes
  virtual const Probleme_FT_Disc_gen& probleme_ft() const;
  virtual Probleme_FT_Disc_gen&        probleme_ft() ;
  virtual void calculer_champ_forces_superficielles(const Maillage_FT_Disc& maillage,
                                                    const Champ_base& gradient_indicatrice,
                                                    Champ_base& potentiel_elements,
                                                    Champ_base& potentiel_faces,
                                                    Champ_base& champ) const;
  virtual void calculer_gradient_indicatrice(const Champ_base& indicatrice,
                                             const DoubleTab& distance_interface_sommets,
                                             Champ_base& gradient_i);

  REF(Probleme_FT_Disc_gen)  probleme_ft_;

  // Masse volumique calculee aux elements
  Champ_Fonc champ_rho_elem_;
  // Masse volumique calculee pour les volumes de controle de la vitesse
  // (pour division   v = (rho.v) / rho et pour matrice de pression)
  Champ_Fonc champ_rho_faces_;
  // Viscosite dynamique (calcul dans preparer_pas_de_temps)
  // champ du type requis pour l'operateur diffusion.
  Champ_Don champ_mu_;
  // Viscosite cinematique pour le calcul du pas de temps de diffusion
  Champ_Don champ_nu_;

protected:


private:
  const Navier_Stokes_FT_Disc_interne& variables_internes() const;
  Navier_Stokes_FT_Disc_interne& variables_internes();

  // Ne pas utiliser ce pointeur : utiliser variables_internes() a la place !
  Navier_Stokes_FT_Disc_interne *variables_internes_;

  double minx,maxx,pente;
  int is_repulsion;
};


#endif
