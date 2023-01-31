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
// File:        Remailleur_Collision_FT_Thomas.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Remailleur_Collision_FT_Thomas_included
#define Remailleur_Collision_FT_Thomas_included

#include <Remailleur_Collision_FT_Juric.h>
#include <Transport_Interfaces_FT_Disc.h>

#include <TRUSTTabs_forward.h>
#include <TRUSTTabs_forward.h>
#include <TRUSTTabs_forward.h>
class Domaine_dis;
class Domaine_dis_base;

/*! @brief : class Remailleur_Collision_FT_Thomas Cette classe implemente les procedures de remaillage des interfaces pour le Front-Tracking :
 *
 *
 *
 * @sa Transport_Interfaces_FT_Disc Maillage_FT_Disc
 */
class Remailleur_Collision_FT_Thomas : public Remailleur_Collision_FT_Juric
{
  Declare_instanciable_sans_constructeur(Remailleur_Collision_FT_Thomas);

public:
  Remailleur_Collision_FT_Thomas();
  int traite_RuptureCoalescenceInterfaces_Conservatif(Maillage_FT_Disc&, Champ_base&) override;

  //Fonctions d'acces aux attributs
  inline const IntTab& voisinage_sommet() const;
  inline IntTab& voisinage_sommet();
  inline const IntTab& next_elem() const;
  inline IntTab& next_elem();
  inline const IntTab& distance_interface_element_eulerien() const;
  inline IntTab& distance_interface_element_eulerien();
  inline const IntTab& nombre_de_voisins_plus_proches() const;
  inline IntTab& nombre_de_voisins_plus_proches();
  inline const DoubleTab& volume_perdu() const;
  inline DoubleTab& volume_perdu();
  inline const DoubleTab& surface_interface_elements_voisins() const;
  inline DoubleTab& surface_interface_elements_voisins();

  inline int plus_grande_distance_interface_element_eulerien() const;
  inline int distance_utilisateur() const;
  inline int est_dimensionne() const;

protected:

  // Fonction qui dimensionne les attributs de type IntTab et DoubleTab
  int initialiser_data(const Maillage_FT_Disc&);

  // Fonction qui construit le voisinage de tous les sommets "s" du maillage eulerien
  // Voisinage d'un sommet "s" := ensemble des elements du maillage eulerien possedant ce sommet "s"
  int construire_voisinage_sommet(const Maillage_FT_Disc&);

  // Fonction qui calcule la distance separant les elements euleriens de la
  // facette de l'interface qui leur est la plus proche, le nombre des voisins
  // d'un element et la surface d'interface coupant certains elements
  int mettre_a_jour_data(const Maillage_FT_Disc&);

  // Fonction qui distribue la perte de volume due au remaillage de l'interface,
  // a un element "elem" donne qui est traverse par la nouvelle interface remaillee.
  int transport_volume_perdu_sur_element(const int,const Maillage_FT_Disc&) ;

  // Fonction qui distribue la perte de volume due au remaillage de l'interface,
  double transport_volume_perdu_sur_element(const Maillage_FT_Disc&) ;

  //Fonction qui a chaque sommet du maillage de l'interface associe une partie du
  //volume perdu dans l'element "elem" lors du remaillage
  int transport_volume_perdu_sur_sommet(const int, ArrOfDouble&, const Maillage_FT_Disc&) const;

  //Fonction qui a chaque sommet du maillage de l'interface associe une partie du
  //volume perdu lors du remaillage
  double transport_volume_perdu_sur_sommet(ArrOfDouble&, const Maillage_FT_Disc&) const;

  //Fonction qui donne la liste des elements voisins a un element "elem" (numero global) donne
  int elements_voisins(const int, ArrOfInt&, const Domaine_dis_base&) const;

  //Fonction qui donne la liste des elements voisins a un element "elem" (numero global) donne et
  //situes a une distance de l'interface strictement plus petite que la distance separant "elem"
  //de l'interface. (normalement la difference des deux distances doit etre de 1)
  int elements_voisins_a_distance_plus_petite2(const int, ArrOfInt&) const;

  // Fonction qui donne le nombre d'elements voisins a un element "elem" donne
  int nb_elements_voisins(const int, const Domaine_dis_base&) const;

  //Soit "elem" (numero global) un element du maillage eulerien a distance "n" de l'interface
  //La fonction renvoie le nombre d'elements voisins de "elem" a distance "n-1"
  int nb_elements_voisins_a_distance_plus_petite(const int, const Domaine_dis_base&) const;
  int nb_elements_voisins_a_distance_plus_petite(const int) const;


private :

  //Attributs de classe
  IntTab voisinage_sommet_;
  IntTab next_elem_;
  IntTab distance_interface_element_eulerien_;
  IntTab nombre_de_voisins_plus_proches_;//pour le parallele

  DoubleTab volume_perdu_;
  DoubleTab surface_interface_elements_voisins_;//pour le parallele

  int plus_grande_distance_interface_element_eulerien_;
  int distance_utilisateur_;
  //int conservation_volume_autorisee_;
  int tester_;
  int est_dimensionne_;

  ArrOfBit tmp_flag_elements_; // tableau temporaire utilise dans elements_voisins()
  //Fonction qui renvoie la surface occupee par l'interface
  //dans un element donne
  double surface_intersection(const int, const Maillage_FT_Disc&) const;

  //Fonction qui renvoie la domaine discrete dans laquelle evolue l'interface
  inline const Domaine_dis& domaine_dis(const Maillage_FT_Disc&) const;

  //Fonction qui renvoie le Remaillage conservatif en volume
  inline const Remaillage_FT& remaillage_FT(const Maillage_FT_Disc&) const;
  inline Remaillage_FT& remaillage_FT(Maillage_FT_Disc&);

  //Fonctions destinees a tester les autres fonctions membre de la classe
  void tester(const Maillage_FT_Disc&) const;
  void tester_interface(const Maillage_FT_Disc&) const;
  void tester_voisinage(const Maillage_FT_Disc&) const;
  void tester_distance(const Maillage_FT_Disc&) const;
  void tester_liste(const Maillage_FT_Disc&) const;
  void tester_transport(const Maillage_FT_Disc&) const;
  void tester_transport_complet(const Maillage_FT_Disc&, DoubleTab&) const;
  void tester_volume_par_sommet(const Maillage_FT_Disc&, const DoubleTab&) const;
};


inline const IntTab& Remailleur_Collision_FT_Thomas::voisinage_sommet() const
{
  return voisinage_sommet_;
}

inline IntTab& Remailleur_Collision_FT_Thomas::voisinage_sommet()
{
  return voisinage_sommet_;
}

inline const IntTab& Remailleur_Collision_FT_Thomas::next_elem() const
{
  return next_elem_;
}

inline IntTab& Remailleur_Collision_FT_Thomas::next_elem()
{
  return next_elem_;
}

inline const IntTab& Remailleur_Collision_FT_Thomas::distance_interface_element_eulerien() const
{
  return distance_interface_element_eulerien_;
}

inline IntTab& Remailleur_Collision_FT_Thomas::distance_interface_element_eulerien()
{
  return distance_interface_element_eulerien_;
}

inline const DoubleTab& Remailleur_Collision_FT_Thomas::volume_perdu() const
{
  return volume_perdu_;
}

inline DoubleTab& Remailleur_Collision_FT_Thomas::volume_perdu()
{
  return volume_perdu_;
}

inline int Remailleur_Collision_FT_Thomas::plus_grande_distance_interface_element_eulerien() const
{
  return plus_grande_distance_interface_element_eulerien_;
}

inline const Domaine_dis& Remailleur_Collision_FT_Thomas::domaine_dis(const Maillage_FT_Disc& maillage) const
{
  return maillage.equation_transport().domaine_dis();
}

inline const Remaillage_FT& Remailleur_Collision_FT_Thomas::remaillage_FT(const Maillage_FT_Disc& maillage) const
{
  return maillage.equation_transport().remaillage_interface();
}

inline Remaillage_FT& Remailleur_Collision_FT_Thomas::remaillage_FT(Maillage_FT_Disc& maillage)
{
  return maillage.equation_transport().remaillage_interface();
}

inline const DoubleTab& Remailleur_Collision_FT_Thomas::surface_interface_elements_voisins() const
{
  return surface_interface_elements_voisins_;
}

inline DoubleTab& Remailleur_Collision_FT_Thomas::surface_interface_elements_voisins()
{
  return surface_interface_elements_voisins_;
}

inline const IntTab& Remailleur_Collision_FT_Thomas::nombre_de_voisins_plus_proches() const
{
  return nombre_de_voisins_plus_proches_;
}

inline IntTab& Remailleur_Collision_FT_Thomas::nombre_de_voisins_plus_proches()
{
  return nombre_de_voisins_plus_proches_;
}

inline int Remailleur_Collision_FT_Thomas::distance_utilisateur() const
{
  return distance_utilisateur_;
}

inline int Remailleur_Collision_FT_Thomas::est_dimensionne() const
{
  return est_dimensionne_;
}

#endif
