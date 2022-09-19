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
// File:        Intersections_Elem_Facettes_Data.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/1
//
//////////////////////////////////////////////////////////////////////////////

#include <Intersections_Elem_Facettes_Data.h>

// ===========================================================================
// Description:
// Ajoute une entree a la liste doublement chainee d'intersections entre
// la facette d'interface num_facette et l'element eulerien num_element.
// Le numero d'element doit verifier 0 <= num_element < zone.nb_elem()
// Si le numero de facette est superieur a la taille de l'index des facettes,
// on agrandit l'index.
void Intersections_Elem_Facettes::ajoute_intersection(int num_facette,
                                                      int num_element,
                                                      double fraction_surface_intersection,
                                                      double contrib_volume_phase1,
                                                      double barycentre_phase1[3],
                                                      double barycentre_u,
                                                      double barycentre_v,
                                                      double barycentre_w)
{
  // Verification de taille de l'index des facettes
  int nb_facettes = index_facette_element_.size_array();
  if (num_facette >= nb_facettes)
    {
      index_facette_element_.resize_array(num_facette+1);
      for (int i = nb_facettes; i <= num_facette; i++)
        index_facette_element_[i] = -1;
    }

  // Verification de taille du tableau de donnees des intersections
  if (data_allocated_size < data_real_size+1)
    {
      data_allocated_size = (data_real_size+1)*2;
      Intersections_Elem_Facettes_Data *new_data =
        new Intersections_Elem_Facettes_Data[data_allocated_size];
      if (data_real_size > 0)
        for (int i=0; i<data_real_size; i++)
          new_data[i] = data[i];
      if (data)
        delete[] data;
      data = new_data;
    }
  // Ajout de la nouvelle entree en debut de liste (pour eviter d'avoir a
  // chercher la fin...)
  {
    // Indice de l'intersection entre l'element et la premiere facette
    int& lindex_elem = index_elem_facette_[num_element];
    // Indice de l'intersection entre la facette et le premier element
    int& lindex_facette = index_facette_element_[num_facette];
    Intersections_Elem_Facettes_Data& new_entry = data[data_real_size];
    // Indice de l'intersection de la meme facette avec l'element suivant
    new_entry.index_element_suivant_ = lindex_facette;
    // Indice de l'intersection du meme element avec la facette suivante
    new_entry.index_facette_suivante_ = lindex_elem;
    new_entry.numero_facette_ = num_facette;
    new_entry.numero_element_ = num_element;
#ifdef AVEC_BUG_SURFACES
    new_entry.surface_intersection_ = fraction_surface_intersection;
#else
    new_entry.fraction_surface_intersection_ = fraction_surface_intersection;
#endif
    new_entry.contrib_volume_phase1_ = contrib_volume_phase1;
    new_entry.barycentre_phase1_[0] = barycentre_phase1[0];
    new_entry.barycentre_phase1_[1] = barycentre_phase1[1];
    new_entry.barycentre_phase1_[2] = barycentre_phase1[2];
    new_entry.barycentre_[0] = barycentre_u;
    new_entry.barycentre_[1] = barycentre_v;
    new_entry.barycentre_[2] = barycentre_w;
    lindex_elem = data_real_size;
    lindex_facette = data_real_size;
  }
  data_real_size++;
}

Intersections_Elem_Facettes::Intersections_Elem_Facettes()
  : data_allocated_size(0), data_real_size(0), data(0)
{
}

Intersections_Elem_Facettes::~Intersections_Elem_Facettes()
{
  if (data)
    delete[] data;
  data = 0;
}

void Intersections_Elem_Facettes::reset(int nb_elements_euleriens,
                                        int nb_facettes)
{
  data_real_size = 0;
  index_elem_facette_.resize_array(nb_elements_euleriens);
  index_elem_facette_ = -1;
  index_facette_element_.resize_array(nb_facettes);
  index_facette_element_ = -1;
}

void Intersections_Elem_Facettes::get_liste_elements_traverses(int num_facette,
                                                               ArrOfInt& liste_elements) const
{
  liste_elements.resize_array(0);
  if (num_facette >= index_facette_element_.size_array())
    // Aucune intersection n'a ete ajoutee pour cette facette,
    // elle ne figure meme pas dans l'index
    return;

  int index = index_facette_element_[num_facette];
  for (; index >= 0; index = data[index].index_element_suivant_)
    {
      const int element = data[index].numero_element_;
      liste_elements.append_array(element);
    }
}

void Intersections_Elem_Facettes::get_liste_facettes_traversantes(int num_element,
                                                                  ArrOfInt& liste_facettes) const
{
  liste_facettes.resize_array(0);
  int index = index_elem_facette_[num_element];
  for (; index >= 0; index = data[index].index_facette_suivante_)
    {
      const int facette = data[index].numero_facette_;
      liste_facettes.append_array(facette);
    }
}

// Description:
// Renvoie un tableau de taille zone.nb_elem():
//  pour un element 0 <= elem < zone.nb_elem(),
//  index_elem()[elem] est l'indice de la premiere intersection entre l'element
//  et les facettes du maillage lagrangien (voir description de la classe)
const ArrOfInt& Intersections_Elem_Facettes::index_elem() const
{
  return index_elem_facette_;
}

// Description:
// Renvoie un tableau de taille "nombre de facettes de l'interface"
//  pour un element 0 <= facette < nb_facettes,
//  index_facette()[facette] est l'indice de la premiere intersection entre
//  la facette et les elements du maillage lagrangien
//  (voir description de la classe)
const ArrOfInt& Intersections_Elem_Facettes::index_facette() const
{
  return index_facette_element_;
}

// Description:
//operateur de copie
const Intersections_Elem_Facettes& Intersections_Elem_Facettes::operator=(const Intersections_Elem_Facettes& ief)
{
  if (&ief != this)
    {
      index_elem_facette_ = ief.index_elem();
      assert(index_elem_facette_.size_array()>0);
      index_facette_element_ = ief.index_facette();
      data_allocated_size = ief.data_allocated_size;
      data_real_size = ief.data_real_size;

      if (data)
        delete data;
      data = new Intersections_Elem_Facettes_Data[data_allocated_size];
      int i;
      for (i=0 ; i<data_real_size; i++)
        {
          data[i] = ief.data[i];
        }
    }
  return *this;
}
