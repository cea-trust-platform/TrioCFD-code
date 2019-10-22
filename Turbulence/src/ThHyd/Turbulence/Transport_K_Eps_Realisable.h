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
// File:        Transport_K_Eps_Realisable.h
// Directory:   $TRUST_ROOT/src/ThHyd/Turbulence
// Version:     /main/15
//
//////////////////////////////////////////////////////////////////////////////

// .SECTION voir aussi
//  Transport_K_Eps
//////////////////////////////////////////////////////////////////////////////

#ifndef Transport_K_Eps_Realisable_included
#define Transport_K_Eps_Realisable_included

#include <Transport_K_Eps_base.h>
#include <Op_Diff_K_Eps_Bas_Re_base.h>
#include <Op_Diff_K_Eps_base.h>
#include <Operateur_Conv.h>
#include <Modele_Fonc_Realisable.h>
#include <Ref_Modele_Fonc_Realisable.h>

class Motcle;

class Transport_K_Eps_Realisable : public Transport_K_Eps_base
{

  Declare_instanciable(Transport_K_Eps_Realisable);

public :

  void set_param(Param& titi);
  int lire_motcle_non_standard(const Motcle&, Entree&);
  virtual const Champ_Don& diffusivite_pour_transport();
  virtual const Champ_base& vitesse_pour_transport();
  int nombre_d_operateurs() const;
  const Operateur& operateur(int) const;
  Operateur& operateur(int);

  void associer_milieu_base(const Milieu_base&);
  void associer_modele_turbulence(const Mod_turb_hyd_RANS&);
  inline const Modele_Fonc_Realisable& modele_fonc() const;
  inline  Modele_Fonc_Realisable& modele_fonc();
//   inline const Champ_Inc& vitesse_transportante();
  const Motcle& domaine_application() const;
  void completer();

protected:
  int with_nu_;
  Op_Diff_K_Eps terme_diffusif;
  Operateur_Conv terme_convectif;

  REF(Champ_Inc) inco_eqn_associee;
  Champ_Don Champ_don_nul_;  // on y met 0 si on ne veut pas de nu

private :

  REF( Modele_Fonc_Realisable ) mon_modele_fonc;

};


// Fonctions inline:

inline const Modele_Fonc_Realisable& Transport_K_Eps_Realisable::modele_fonc() const
{
  return mon_modele_fonc.valeur();
}

inline Modele_Fonc_Realisable& Transport_K_Eps_Realisable::modele_fonc()
{
  return mon_modele_fonc.valeur();
}


#endif


