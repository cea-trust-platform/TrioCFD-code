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
// File:        Transport_K_Eps_Bas_Reynolds.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Equations/RANS
//
//////////////////////////////////////////////////////////////////////////////

// .SECTION voir aussi
//  Transport_K_Eps
//////////////////////////////////////////////////////////////////////////////

#ifndef Transport_K_Eps_Bas_Reynolds_included
#define Transport_K_Eps_Bas_Reynolds_included

#include <Transport_K_Eps_non_std.h>
#include <Modele_Fonc_Bas_Reynolds_Base.h>
#include <TRUST_Ref.h>

class Motcle;

class Transport_K_Eps_Bas_Reynolds : public Transport_K_Eps_non_std
{

  Declare_instanciable(Transport_K_Eps_Bas_Reynolds);

public :

  void associer_milieu_base(const Milieu_base&) override;
  void associer_modele_turbulence(const Modele_turbulence_hyd_RANS_K_Eps_base&) override;
  inline const OWN_PTR(Modele_Fonc_Bas_Reynolds_Base)& modele_fonc() const;
  inline  OWN_PTR(Modele_Fonc_Bas_Reynolds_Base)& modele_fonc();
  inline const Champ_Inc_base& vitesse_transportante();
  const Motcle& domaine_application() const override;
  void completer() override;

private :

  REF(OWN_PTR(Modele_Fonc_Bas_Reynolds_Base)) mon_modele_fonc;

};


// Fonctions inline:

inline const OWN_PTR(Modele_Fonc_Bas_Reynolds_Base)& Transport_K_Eps_Bas_Reynolds::modele_fonc() const
{
  return mon_modele_fonc.valeur();
}

inline OWN_PTR(Modele_Fonc_Bas_Reynolds_Base)& Transport_K_Eps_Bas_Reynolds::modele_fonc()
{
  return mon_modele_fonc.valeur();
}


#endif


