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
// File:        Intersections_Elem_Facettes_Data.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/2
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Intersections_Elem_Facettes_Data_included
#define Intersections_Elem_Facettes_Data_included
#include <TRUSTTabFT.h>
// ====================================================================
class Intersections_Elem_Facettes_Data
{
public:
#ifdef AVEC_BUG_SURFACES
  // Surface de l'intersection facette/element
  // !!!!! gros bug : en fait, quand on remplit la structure, on y met la fraction
  //  de surface.
  double surface_intersection_;
#else
  // Fraction surface_intersection_facette_element / surface_totale_facette:
  double fraction_surface_intersection_;
#endif
  double contrib_volume_phase1_;
  // Coordonnees barycentriques du centre de gravite de l'intersection
  // par rapport aux trois sommets de la facette. On a toujours
  // barycentre[0] + barycentre[1] + barycentre[2] = 1.
  // En 2D, on a barycentre_[2] = 0.;
  double barycentre_[3];

  int index_element_suivant_;  // -1 si dernier element de la liste
  int index_facette_suivante_; // idem.
  int numero_facette_;
  int numero_element_;
};

// ====================================================================
// .DESCRIPTION        : class Intersections_Elem_Facettes
//
// Cette classe contient les donnees des intersections entre les facettes
// de l'interface et les elements euleriens sous la forme d'une liste
// doublement chainee.
// Pour parcourir les facettes qui coupent un element "elem", on fait:
//
//   int index=index_elem()[elem];
//   while (index >= 0) {
//     const Intersections_Elem_Facettes_Data & data = data_intersection(index);
//     ... // faire quelque chose avec data
//     index = data.index_facette_suivante_;
//   }
//
// Pour parcourir les elements qui sont coupes par une facette "facette":
//
//   int index=index_facette()[facette];
//   while (index >= 0) {
//     const Intersections_Elem_Facettes_Data & data = data_intersection(index);
//     ... // faire quelque chose avec data
//     index = data.index_element_suivant_;
//   }

class Intersections_Elem_Facettes
{
public:
  Intersections_Elem_Facettes();
  ~Intersections_Elem_Facettes();
  void get_liste_elements_traverses(int num_facette,
                                    ArrOfInt& liste_elements) const;
  void get_liste_facettes_traversantes(int num_element,
                                       ArrOfInt& liste_facettes) const;
  void ajoute_intersection(int num_facette,
                           int num_element,
                           double surface_intersection,
                           double contrib_volume_phase1,
                           double barycentre_u,
                           double barycentre_v,
                           double barycentre_w);

  // A faire avant la premiere utilisation et pour remettre a zero
  // nb_elements_euleriens doit etre correct.
  // nb_facettes peut etre une estimation seulement (en pratique lors du
  // parcours on ne peut pas prevoir le nombre final de facettes)
  void reset(int nb_elements_euleriens=0, int nb_facettes=0);

  const ArrOfInt& index_elem() const;
  const ArrOfInt& index_facette() const;
  inline const Intersections_Elem_Facettes_Data& data_intersection(int index) const;
  inline Intersections_Elem_Facettes_Data& get_set_data_intersection(int index);

  //operateur de copie
  const Intersections_Elem_Facettes& operator=(const Intersections_Elem_Facettes& ief);

private:
  // Tableau de taille "nombre d'elements euleriens"
  // Pour chaque element, index de l'entree correspondant a la premiere facette
  // dans data.
  //  -1 si l'element ne contient pas de facette
  ArrOfIntFT index_elem_facette_;
  // pour chaque facette, index de l'entree correspondant au premier element
  //  -1 si la facette ne traverse aucun element (facette qui devrait etre effacee)
  ArrOfIntFT index_facette_element_;

  //
  // Les donnees : un element par couple facette/element traverse.
  // (politique d'allocation : on agrandit d'un facteur 2, on ne diminue jamais).
  int data_allocated_size; // Le nombre d'elements alloues
  int data_real_size;      // Le nombre d'elements reellement utilises
  Intersections_Elem_Facettes_Data * data;
};

// Description:
//  Renvoie les donnees de l'intersection stockee a l'indice "index"
//  dans le tableau "data" ( 0 <= index < data_real_size_ )
inline const Intersections_Elem_Facettes_Data& Intersections_Elem_Facettes::data_intersection(int index) const
{
  assert(index >= 0 && index < data_real_size);
  return data[index];
}
// Description:
//  Renvoie les donnees de l'intersection stockee a l'indice "index"
//  dans le tableau "data" ( 0 <= index < data_real_size_ )
// ATTENTION A SON UTILISATION !!!
inline Intersections_Elem_Facettes_Data& Intersections_Elem_Facettes::get_set_data_intersection(int index)
{
  assert(index >= 0 && index < data_real_size);
  return data[index];
}
#endif
