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
// File:        Champ_front_contact_rayo_semi_transp_VEF.h
// Directory:   $TRUST_ROOT/src/Rayonnement_semi_transp/VEF
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Champ_front_contact_rayo_semi_transp_VEF_included
#define Champ_front_contact_rayo_semi_transp_VEF_included

#include <Champ_front_contact_VEF.h>
#include <Ref_Modele_rayo_semi_transp.h>


//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    classe Champ_front_contact_rayo_semi_transp_VEF
//
// .SECTION
//////////////////////////////////////////////////////////////////////////////
class Champ_front_contact_rayo_semi_transp_VEF: public Champ_front_contact_VEF
{
  Declare_instanciable(Champ_front_contact_rayo_semi_transp_VEF);

public :
  virtual int initialiser(double temps, const Champ_Inc_base& inco);
  void calculer_temperature_bord(double temps);
  Champ_front_base& affecter_(const Champ_front_base& ch) ;
  void mettre_a_jour(double temps);
  inline Champ_Inc_base& inconnue1();
  inline const Champ_Inc_base& inconnue1() const;

  inline Champ_Inc_base& inconnue2();
  inline const Champ_Inc_base& inconnue2() const;

  inline Nom& nom_prob1();
  inline const Nom& nom_prob1() const;
  inline Nom& nom_prob2();
  inline const Nom& nom_prob2() const;

  inline void associer_modele_rayo(const Modele_rayo_semi_transp& modele);
  inline const Modele_rayo_semi_transp& modele_rayo() const;
  void mettre_a_jour_flux_radiatif(); // le copie du modele
  void calcul_grads_locaux(double temps);
  void modifie_gradients_pour_rayonnement(DoubleVect& gradient_num_transf, DoubleVect& gradient_num_transf_autre_pb);


protected :

  DoubleVect flux_radiatif;
  REF(Modele_rayo_semi_transp) le_modele_rayo;

};


inline void Champ_front_contact_rayo_semi_transp_VEF::associer_modele_rayo(const Modele_rayo_semi_transp& modele)
{
  le_modele_rayo = modele;
}

inline const Modele_rayo_semi_transp& Champ_front_contact_rayo_semi_transp_VEF::modele_rayo() const
{
  return le_modele_rayo.valeur();
}

inline Nom& Champ_front_contact_rayo_semi_transp_VEF::nom_prob1()
{
  return nom_pb1;
}

inline const Nom& Champ_front_contact_rayo_semi_transp_VEF::nom_prob1() const
{
  return nom_pb1;
}

inline Nom& Champ_front_contact_rayo_semi_transp_VEF::nom_prob2()
{
  return nom_pb2;
}

inline const Nom& Champ_front_contact_rayo_semi_transp_VEF::nom_prob2() const
{
  return nom_pb2;
}

inline Champ_Inc_base& Champ_front_contact_rayo_semi_transp_VEF::inconnue1()
{
  return l_inconnue1.valeur();
}

inline const Champ_Inc_base& Champ_front_contact_rayo_semi_transp_VEF::inconnue1() const
{
  return l_inconnue1.valeur();
}

inline Champ_Inc_base& Champ_front_contact_rayo_semi_transp_VEF::inconnue2()
{
  return l_inconnue2.valeur();
}

inline const Champ_Inc_base& Champ_front_contact_rayo_semi_transp_VEF::inconnue2() const
{
  return l_inconnue2.valeur();
}

#endif
