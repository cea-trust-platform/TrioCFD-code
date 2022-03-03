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
// File:        Navier_Stokes_phase_field.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Multiphase/Phase_field/src
// Version:     /main/15
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Navier_Stokes_phase_field_included
#define Navier_Stokes_phase_field_included

#include <Navier_Stokes_std.h>
#include <Champ_Don.h>
#include <time.h>
#include <Ref_Champ_Don.h>
#include <Constituant.h>


/*! @brief classe Navier_Stokes_phase_field Cette classe porte les termes de l'equation de la dynamique
 *
 *     pour un fluide sans modelisation de la turbulence.
 *     On suppose l'hypothese de fluide incompressible: div U = 0
 *     On peut 1) soit utiliser un modele a rho variable
 *             2) soit utiliser l'hypothese de Boussinesq
 *     Dans ce cas, on considere la masse volumique constante (egale a rho_0) sauf dans le
 *     terme des forces de gravite.
 *     Sous ces hypotheses, on utilise la forme suivante des equations de
 *     Navier_Stokes:
 *        DU/dt = div(terme visqueux) - gradP/rho_0 + Bt(T-T0)g + autres sources/rho_0
 *        div U = 0
 *     avec DU/dt : derivee particulaire de la vitesse
 *          rho_0 : masse volumique de reference
 *          T0    : temperature de reference
 *          Bt    : coefficient de dilatabilite du fluide
 *          g     : vecteur gravite
 *     Rq : l'implementation de la classe permet bien sur de negliger
 *          certains termes de l'equation (le terme visqueux, le terme
 *          convectif, tel ou tel terme source).
 *     L'inconnue est le champ de vitesse.
 *
 *     Pour le traitement des cas un peu particulier : ajout de Traitement_particulier
 *     exemple : THI, canal (CA)
 *
 * @sa Equation_base Pb_Hydraulique Pb_Thermohydraulique, Post-traitement du potentiel chimique generalise
 */
class Navier_Stokes_phase_field : public Navier_Stokes_std
{
  Declare_instanciable_sans_constructeur_ni_destructeur(Navier_Stokes_phase_field);

public :

  Navier_Stokes_phase_field();
  ~Navier_Stokes_phase_field() override;
  void set_param(Param& titi) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void discretiser() override;
  void completer() override;
  int preparer_calcul() override;
  virtual void rho_aux_faces(const DoubleTab&, Champ_Don&);
  inline const Champ_Don& rho() const;
  inline const Champ_Don& drhodc() const;
  inline const double& rho0() const;
  inline Champ_Don& mu();
  inline const Champ_Don& mu() const;
  inline int& get_boussi_();
  inline DoubleVect& get_g_();
  inline const DoubleVect& get_g_() const;
  DoubleTab& derivee_en_temps_inco(DoubleTab& ) override;
  void mettre_a_jour(double ) override;
  const Champ_Don& diffusivite_pour_transport() const override;
  void creer_champ(const Motcle& motlu) override;

  /////////////////////////////////////////////////////

  virtual void calculer_rho(const bool init=false);
  virtual void calculer_mu(const bool init=false);

  double calculer_pas_de_temps() const override;

  inline int& getset_terme_source_();
  inline int& getset_compressible_();
  inline int& get_diff_boussi_();

protected :

  Champ_Don rho_; // soit lu dans le jdd (boussi_==0)
  // soit construit comme rho0_*(1+betac_*c) (boussi_==1)
  Champ_Don drhodc_; // soit lu dans le jdd (boussi_==0)
  // soit construit comme rho0_*betac_ (boussi_==1 mais attention strictement valide seulement si betac_= Champs_Uniforme)
  double rho0_;
  REF(Champ_Don) betac_;
  Champ_Don mu_;// lu dans le jdd (boussi_==0) seulement pour boussi_==0 && diff_boussi_==0

  int boussi_;
  int diff_boussi_;
  DoubleVect g_;
  int terme_source_;
  int compressible_;

  void _aff_donnee_P0(const DoubleTab& , const Motcle& ) const;
  void _aff_donnee_P1(const DoubleTab& , const Motcle& ) const;
  void _aff_donnee_P1Bulle(const DoubleTab& , const Motcle& ) const;


};

// Methode de calcul de la valeur sur un champ aux elements d'un champ uniforme ou non a plusieurs composantes
inline double& valeur(DoubleTab& valeurs, const int elem, const int dim)
{
  return valeurs(elem,dim);
}

inline const Champ_Don& Navier_Stokes_phase_field::rho() const
{
  return rho_;
}

inline const Champ_Don& Navier_Stokes_phase_field::drhodc() const
{
  return drhodc_;
}

inline const double& Navier_Stokes_phase_field::rho0() const
{
  return rho0_;
}

inline Champ_Don& Navier_Stokes_phase_field::mu()
{
  return mu_;
}

inline const Champ_Don& Navier_Stokes_phase_field::mu() const
{
  return mu_;
}

inline int& Navier_Stokes_phase_field::get_boussi_()
{
  return boussi_;
}

inline DoubleVect& Navier_Stokes_phase_field::get_g_()
{
  return g_;
}

inline const DoubleVect& Navier_Stokes_phase_field::get_g_() const
{
  return g_;
}

inline int& Navier_Stokes_phase_field::get_diff_boussi_()
{
  return diff_boussi_;
}

inline int& Navier_Stokes_phase_field::getset_terme_source_()
{
  return terme_source_;
}

inline int& Navier_Stokes_phase_field::getset_compressible_()
{
  return compressible_;
}

#endif
