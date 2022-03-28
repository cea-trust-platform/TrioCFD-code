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
// File:        Proprietes_part_vol.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Proprietes_part_vol_included
#define Proprietes_part_vol_included

#include <Objet_U.h>
#include <TRUSTTabFT.h>
class Param;

//Description
//Classe qui porte les proprietes de particules
//Actuellement seules les particules materielles sont dotees de proprietes
//Proprietes dynamiques   : vitesse_p_, delat_v_
//Proprietes energetiques : temperature_p_
//Proprietes physiques    : masse_volumique_p_
//Proprietes geometriques : diametre_p_, volume_p_,

//Toute operation modifiant le nombre de particules dans le domaine doit s accompagner
//d une actualisation de l attribut nb_particules_ de cette classe
//Les operations susceptibles de modifier le nombre de particules sont actuellement :
//-lecture de la condition initiale
//-injection-transformation
//-nettoyage
//-passage d un sous domaine a un autre dans le cas d un calcul parallele


class Proprietes_part_vol : public Objet_U
{
  Declare_instanciable_sans_constructeur(Proprietes_part_vol);
public:

  Proprietes_part_vol();

  void set_param(Param& param);
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void completer();
  void fixer_nb_particules(const int& nb_part);
  void lire_distribution(Entree& is);
  void nettoyer(const ArrOfInt& som_utilises);
  void ajouter_proprietes(const Proprietes_part_vol& proprietes_tmp);

  inline DoubleTab& vitesse_particules();
  inline const DoubleTab& vitesse_particules() const;
  inline DoubleTab& delta_v();
  inline const DoubleTab& delta_v() const;
  inline DoubleTab& temperature_particules();
  inline const DoubleTab& temperature_particules() const;
  inline DoubleTab& masse_vol_particules();
  inline const DoubleTab& masse_vol_particules() const;
  inline DoubleTab& diametre_particules();
  inline const DoubleTab& diametre_particules() const;
  inline DoubleTab& volume_particules();
  inline const DoubleTab& volume_particules() const;
  inline int nb_particules() const;

protected:

  DoubleTabFT vitesse_p_;           //vitesse des particule
  DoubleTabFT deltat_v_;           //difference entre vitesse vitesse du fluide
  //(a la position de la particule) et vitesse de la particule
  DoubleTabFT temperature_p_;     //temperature
  DoubleTabFT masse_volumique_p_; //densite
  DoubleTabFT diametre_p_;           //masse_volumique
  DoubleTabFT volume_p_;           //volume

  int nb_particules_;           //nombre de particules portees par l objet

};

inline DoubleTab& Proprietes_part_vol::vitesse_particules()
{
  return vitesse_p_;
}

inline const DoubleTab& Proprietes_part_vol::vitesse_particules() const
{
  return vitesse_p_;
}

inline DoubleTab& Proprietes_part_vol::delta_v()
{
  return deltat_v_;
}

inline const DoubleTab& Proprietes_part_vol::delta_v() const
{
  return deltat_v_;
}

inline DoubleTab& Proprietes_part_vol::temperature_particules()
{
  return temperature_p_;
}

inline const DoubleTab& Proprietes_part_vol::temperature_particules() const
{
  return temperature_p_;
}

inline DoubleTab& Proprietes_part_vol::masse_vol_particules()
{
  return masse_volumique_p_;
}

inline const DoubleTab& Proprietes_part_vol::masse_vol_particules() const
{
  return masse_volumique_p_;
}

inline DoubleTab& Proprietes_part_vol::diametre_particules()
{
  return diametre_p_;
}

inline const DoubleTab& Proprietes_part_vol::diametre_particules() const
{
  return diametre_p_;
}

inline DoubleTab& Proprietes_part_vol::volume_particules()
{
  return volume_p_;
}

inline const DoubleTab& Proprietes_part_vol::volume_particules() const
{
  return volume_p_;
}

inline int Proprietes_part_vol::nb_particules() const
{
  return nb_particules_;
}

#endif
