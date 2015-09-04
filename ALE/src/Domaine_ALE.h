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
// File:        Domaine_ALE.h
// Directory:   $TRUST_ROOT/src/ALE
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Domaine_ALE_included
#define Domaine_ALE_included

#include <Domaine.h>
#include <DoubleLists.h>
#include <IntLists.h>
#include <Champs_front.h>
#include <Champ_P1NC.h>
class Domaine_dis;
//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    Classe Domaine_ALE
//    Cette classe est un interprete qui sert a lire l'attribut axi.
//    Directive:
//        Domaine_ALE
//    Cette directive optionelle permets de faire les calculs en
//    coordonnees cylindriques. En l'absence de cette directive les calculs
//    se font en coordonnees cartesiennes.
// .SECTION voir aussi
//    Interprete Objet_U
//////////////////////////////////////////////////////////////////////////////
class Domaine_ALE : public Domaine
{
  Declare_instanciable(Domaine_ALE);

public :

  inline const double& get_dt() const;
  virtual void set_dt(double& dt_);
  inline const DoubleTab& vitesse() const;
  inline DoubleTab& vitesse_faces();
  inline const DoubleTab& vitesse_faces() const;
  virtual void mettre_a_jour (double temps, Domaine_dis&, Probleme_base&);
  virtual void initialiser (double temps, Domaine_dis&, Probleme_base&);
  DoubleTab calculer_vitesse(double temps,Domaine_dis&, Probleme_base&);
  DoubleTab& calculer_vitesse_faces(DoubleTab&, int, int, IntTab&);
  virtual void lecture_vit_bords_ALE(Entree& is);
  virtual DoubleTab& laplacien(Domaine_dis&, Probleme_base&, const DoubleTab&, DoubleTab&);


protected:

  double dt;
  DoubleTab la_vitesse;
  DoubleTab vf;
  IntTab som_faces_bords;
  SolveurSys solv;
  Matrice_Morse_Sym mat;
  Champs_front les_champs_front;
  int nb_bords_ALE;
  Bords les_bords_ALE;
};


inline const DoubleTab& Domaine_ALE::vitesse() const
{
  return la_vitesse;
}

inline const double& Domaine_ALE::get_dt() const
{
  return dt;
}

inline DoubleTab& Domaine_ALE::vitesse_faces()
{
  return vf;
}
inline const DoubleTab& Domaine_ALE::vitesse_faces() const
{
  return vf;
}
#endif
