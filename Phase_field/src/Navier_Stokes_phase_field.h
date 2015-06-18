/****************************************************************************
* Copyright (c) 2015, CEA
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
// Directory:   $TRUST_ROOT/src/Phase_field
// Version:     /main/15
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Navier_Stokes_phase_field_included
#define Navier_Stokes_phase_field_included

#include <Navier_Stokes_std.h>
#include <Champ_Don.h>
#include <time.h>

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    classe Navier_Stokes_phase_field
//    Cette classe porte les termes de l'equation de la dynamique
//    pour un fluide sans modelisation de la turbulence.
//    On suppose l'hypothese de fluide incompressible: div U = 0
//    On peut 1) soit utiliser un modele a rho variable
//            2) soit utiliser l'hypothese de Boussinesq
//    Dans ce cas, on considere la masse volumique constante (egale a rho_0) sauf dans le
//    terme des forces de gravite.
//    Sous ces hypotheses, on utilise la forme suivante des equations de
//    Navier_Stokes:
//       DU/dt = div(terme visqueux) - gradP/rho_0 + Bt(T-T0)g + autres sources/rho_0
//       div U = 0
//    avec DU/dt : derivee particulaire de la vitesse
//         rho_0 : masse volumique de reference
//         T0    : temperature de reference
//         Bt    : coefficient de dilatabilite du fluide
//         g     : vecteur gravite
//    Rq : l'implementation de la classe permet bien sur de negliger
//         certains termes de l'equation (le terme visqueux, le terme
//         convectif, tel ou tel terme source).
//    L'inconnue est le champ de vitesse.
//
//    Pour le traitement des cas un peu particulier : ajout de Traitement_particulier
//    exemple : THI, canal (CA)
// .SECTION voir aussi
//      Equation_base Pb_Hydraulique Pb_Thermohydraulique
//
// Post-traitement du potentiel chimique generalise
//////////////////////////////////////////////////////////////////////////////
class Navier_Stokes_phase_field : public Navier_Stokes_std
{
  Declare_instanciable_sans_constructeur(Navier_Stokes_phase_field);

public :

  Navier_Stokes_phase_field();
  void set_param(Param& titi);
  int lire_motcle_non_standard(const Motcle&, Entree&);
  virtual void discretiser();
  virtual int preparer_calcul();
  virtual void rho_aux_faces(const DoubleTab&, Champ_Don&);
  inline Champ_Don& rho();
  inline const Champ_Don& rho() const;
  inline Champ_Don& mu();
  inline const Champ_Don& mu() const;
  inline int& get_boussi_();
  inline DoubleVect& get_g_();
  inline const DoubleVect& get_g_() const;
  DoubleTab& derivee_en_temps_inco(DoubleTab& );
  virtual void mettre_a_jour(double );
  const Champ_Don& diffusivite_pour_transport();


  /////////////////////////////////////////////////////

  double calculer_pas_de_temps() const;

  inline int& getset_terme_source_();
  inline int& getset_compressible_();
  inline int& get_diff_boussi_();

protected :

  Champ_Don rho_;
  Champ_Don mu_;
  int boussi_;
  int diff_boussi_;
  DoubleVect g_;
  double rho0_;
  int terme_source_;
  int compressible_;

  void _aff_donnee_P0(const DoubleTab& , const Motcle& ) const;
  void _aff_donnee_P1(const DoubleTab& , const Motcle& ) const;
  void _aff_donnee_P1Bulle(const DoubleTab& , const Motcle& ) const;


};

inline Champ_Don& Navier_Stokes_phase_field::rho()
{
  return rho_;
}

inline const Champ_Don& Navier_Stokes_phase_field::rho() const
{
  return rho_;
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
