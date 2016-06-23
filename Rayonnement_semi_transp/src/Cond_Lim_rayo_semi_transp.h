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
// File:        Cond_Lim_rayo_semi_transp.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Cond_Lim_rayo_semi_transp_included
#define Cond_Lim_rayo_semi_transp_included


#include <Ref_Modele_rayo_semi_transp.h>
#include <Neumann_sortie_libre.h>
#include <Champ_front.h>

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//     classe Cond_Lim_rayo_semi_transp
//
// .SECTION voir aussi
//
//     Ce n'est pas une classe de l'arbre TRUST a elle seule.
//     Cette classe est faite etre une classe mere d'une classe
//     qui heritera par ailleurs d'Objet_U
//////////////////////////////////////////////////////////////////////////////
class Cond_Lim_rayo_semi_transp
{

public:
  virtual void associer_modele(const Modele_rayo_semi_transp&);
  inline const Modele_rayo_semi_transp& modele() const;
  inline Modele_rayo_semi_transp& modele();

  virtual void recherche_emissivite_et_A();

  inline Champ_front& emissivite();
  inline const Champ_front& emissivite() const;
  inline double& A();
  inline const double& A() const;

  virtual const Cond_lim_base& la_cl() const=0;

protected :
  REF(Modele_rayo_semi_transp) mon_modele;
  Champ_front emissivite_;
  double A_;

  inline virtual ~Cond_Lim_rayo_semi_transp();
};

Cond_Lim_rayo_semi_transp::~Cond_Lim_rayo_semi_transp()
{}

// Description:
//    Renvoie la reference sur le modele pointe par
//    Cond_Lim_rayo_semi_transp::mon_modele.
//    (version const)
// Precondition: on doit avoir associe un modele
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Equation_base&
//    Signification: le modele associe a l'objet
//    Contraintes: reference constante
// Exception: pas de modele associe
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
inline const Modele_rayo_semi_transp& Cond_Lim_rayo_semi_transp::modele() const
{
  assert (mon_modele.non_nul());
  return mon_modele.valeur();
}


// Description:
//    Renvoie la reference sur le modele pointe par
//    Cond_Lim_rayo_semi_transp::mon_modele.
// Precondition: on doit avoir associe un modele
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Equation_base&
//    Signification: le modele associe a l'objet
//    Contraintes:
// Exception: pas de modele associe
// Effets de bord:
// Postcondition:
inline Modele_rayo_semi_transp& Cond_Lim_rayo_semi_transp::modele()
{
  assert (mon_modele.non_nul());
  return mon_modele.valeur();
}

inline Champ_front& Cond_Lim_rayo_semi_transp::emissivite()
{
  return emissivite_;
}

inline const Champ_front& Cond_Lim_rayo_semi_transp::emissivite() const
{
  return emissivite_;
}

inline double& Cond_Lim_rayo_semi_transp::A()
{
  return A_;
}

inline const double& Cond_Lim_rayo_semi_transp::A() const
{
  return A_;
}

#endif
