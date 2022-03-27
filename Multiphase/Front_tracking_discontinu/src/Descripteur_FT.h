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
// File:        Descripteur_FT.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////

// C'est un premier jet.
// On peut envisager de stocker autrement :
// pour chaque PE_voisin : un espace distant + un espace virtuel.
// Dans ce cas les communications sont toujours symetriques.

#ifndef Descripteur_FT_included
#define Descripteur_FT_included

#include <Objet_U.h>
#include <ArrOfIntFT.h>
#include <Schema_Comm_FT.h>
#include <Vect_ArrOfInt.h>
#include <MD_Vector.h>
#include <MD_Vector_tools.h>

class ArrOfDoubleFT;
class IntTabFT;
class DoubleTabFT;
#include <TRUSTTabs_forward.h>
class Comm_Group;

// .DESCRIPTION        : class Descripteur_FT
//   Descripteur_FT stocke pour chaque PE une liste de numeros d'elements.
// .SECTION voir aussi
//   Desc_Structure_FT

class Descripteur_FT : public Objet_U
{
  Declare_instanciable_sans_constructeur(Descripteur_FT);
public:
  // Constructeur
  Descripteur_FT();
  Descripteur_FT(const Descripteur_FT&) = default;
  Descripteur_FT& operator=(const Descripteur_FT&);
  void reset();

  // Ajoute l'element au PE_voisin, renvoie le rang de l'element ajoute
  // dans le tableau d'elements du PE.
  // Il faut appeler calcul_liste_pe_voisins()
  int ajoute_element(int PE_voisin, int element);
  int ajoute_elements(int PE_voisin, const ArrOfInt& elements);

  // renvoie les PEs pour lesquels le tableau d'elements n'est pas vide,
  // tries dans l'ordre croissant de numero de PE.
  inline const ArrOfInt& pe_voisins() const;

  // Renvoie la liste des elements distants/virtuels du pe en parametre.
  inline const ArrOfInt& elements(int pe_voisin) const;

  inline const VECT(ArrOfInt)& all_elements() const
  {
    return elements_;
  }

  // Remplace la liste des elements par celle en parametre.
  // Si on vide un tableau ou si on en remplit un qui etait vide, il
  // faut recalculer la liste de pe voisins.
  void set_elements(int PE_voisin, const ArrOfInt& elements);

  // Renvoie "pas zero" si l'element est deja dans le descripteur pour le pe donne,
  // 0 sinon.
  int contient_element(int pe, int element) const;

  // Recalcul de la liste des pe_voisins (triee par ordre croissant)
  // A faire apres avoir modifie les voisinages.
  void calcul_liste_pe_voisins();

  void retirer_elements(const ArrOfInt& nouveau_pe,
                        Descripteur_FT&    elements_retires);

private:

  // Tableau des numeros des PEs voisins (contient exactement la
  // liste des PEs pour lesquels le tableau elements_ est non vide,
  // dans l'ordre croissant).
  ArrOfIntFT pe_voisins_;
  // Tableau des numeros des elements distants ou virtuels
  // (sommets ou facettes) (le vecteur a autant d'elements que le nb de procs,
  // la plupart des tableaux sont generalement vides).
  VECT(ArrOfInt) elements_;

  // Si on ajoute des elements ou si on remplace la liste d'elements,
  // le statut devient BAD et on n'a plus le droit d'acceder a la liste
  // des PEs. Celle-ci doit etre recalculee avec calcul_liste_pe_voisins();
  enum Status { BAD=1, OK=2 } status_;
};

// ******************************************************************************

// .DESCRIPTION        : class Desc_Structure_FT
//
// Desc_Structure_FT est un descripteur adapte aux maillages lagrangiens des
// interfaces. Il contient un tableau d'elements et non une serie d'intervalles
// (le plus souvent, les elements distants et virtuels ne sont pas contigus dans
// les tableaux). De plus, l'allocation memoire est specialisee (les tableaux
// croissent geometriquement mais ne decroissent jamais, leur taille est
// generalement superieure au nombre d'elements qu'ils contiennent).
//
// Definition de la correspondance:
// Pour deux processeurs numerotes A et B quelconques et un indice i,
// on considere les deux numeros suivants
// nA = (Sur processeur A, espace_distant.elements(B) [i])
// nB = (Sur processeur B, espace_virtuel.elements(A) [i])
// On dit que l'espace distant et l'espace virtuel sont en correspondance
// si l'element nA sur le processeur A et l'element nB sur le processeur B
// representent le meme element (sommet ou facette ...)

class Desc_Structure_FT : public Objet_U
{
  Declare_instanciable_sans_constructeur(Desc_Structure_FT);
public:
  // Constructeur par defaut (initialisation des structures)
  Desc_Structure_FT();

  // Vide la structure (la memoire n'est pas liberee)
  void reset(void);

  // Echange des espaces virtuels pour differents types de tableaux
  void echange_espace_virtuel(ArrOfDouble& tab) const;
  void echange_espace_virtuel(ArrOfInt& tab) const;

  // Collection des donnees virtuelles (sommation des valeurs virtuelles
  // sur les differents processeurs et collection sur le noeud reel,
  // puis echange_espace_virtuel)
  void collecter_espace_virtuel(ArrOfDouble& tab, MD_Vector_tools::Operations_echange op) const;
  void collecter_espace_virtuel(ArrOfInt& tab, MD_Vector_tools::Operations_echange op) const;

  Descripteur_FT& espace_distant();               // Status -> BAD
  const Descripteur_FT& espace_distant() const ;
  Descripteur_FT& espace_virtuel();               // Status -> BAD
  const Descripteur_FT& espace_virtuel() const ;

  // Cette fonction doit etre appelee chaque fois qu'on modifie les
  // espaces distants et virtuels
  void calcul_schema_comm(const int nb_items_tot);

  void echanger_elements(const ArrOfInt& nouveau_pe);
  void remplir_element_pe(ArrOfInt& element_pe) const;

  int check() const;

  const Schema_Comm_FT& schema_comm() const;
  const Schema_Comm_FT& schema_comm_inverse() const;

  const MD_Vector& get_md_vector() const;

protected:

  Descripteur_FT espace_distant_;
  Descripteur_FT espace_virtuel_;

  // Si on modifie le descripteur (utilisation de espace_distant()
  // ou espace_virtuel() en version non-const), alors le
  // statut passe a BAD. On recalcule le descripteur au besoin
  enum Status { BAD, OK } status_md_;

  ArrOfInt pe_voisins_;
  VECT(ArrOfInt) blocs_to_recv_;
  MD_Vector md_vector_;
  mutable DoubleVect tmp_doublevect_;
  mutable IntVect tmp_intvect_;

  Schema_Comm_FT schema_comm_;
  Schema_Comm_FT schema_comm_inverse_;
};

// Description:
// Renvoie la liste des PE pour lesquels la liste d'elements est non vide,
// dans l'ordre croissant des numeros de PE.
inline const ArrOfInt& Descripteur_FT::pe_voisins() const
{
  assert(status_ == OK);
  return pe_voisins_;
}

// Description:
// Renvoie la liste des elements distants/virtuels du pe en parametre.
inline const ArrOfInt& Descripteur_FT::elements(int pe_voisin) const
{
  return elements_[pe_voisin];
}

#endif
