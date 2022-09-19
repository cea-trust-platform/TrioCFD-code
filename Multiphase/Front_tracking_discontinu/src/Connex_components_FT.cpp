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
// File:        Connex_components_FT.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/5
//
//////////////////////////////////////////////////////////////////////////////
#include <Connex_components_FT.h>
#include <communications.h>
#include <ArrOfBit.h>
#include <Connectivite_som_elem.h>
#include <Maillage_FT_Disc.h>
#include <Connex_components.h>
#include <Static_Int_Lists.h>
#include <TRUSTArray.h>

/*! @brief Calcul des ensembles connexes par sommets d'elements non "marques" (les elements sont relies entre eux par un graphe
 *
 *   symetrique passant par les faces).
 *   Une portion de domaine connexe porte un numero 0 <= i < N
 *   unique et est delimite soit par un bord, soit par un element voisin "marque"
 *   par num_compo[elem] = -1.
 *   Cette methode est sequentielle (peut etre appelee sur un seul processeur)
 *
 * @param (elem_faces)
 * @param (faces_elem)
 * @param (num_compo)
 */
int search_connex_components_local_FT(const Maillage_FT_Disc& mesh, ArrOfInt& num_compo)
{

  const int dimension=mesh.dimension;
  const int nbelem = mesh.nb_facettes();
  // Construction de la connectivite element-element (pour chaque facette, numeros des 2 ou 3 facettes voisines
  //  -1 si pas de voisin)

  Static_Int_Lists som_elem;
  const IntTab& facettes = mesh.facettes();
  IntTab elems(nbelem, dimension);
  int i, j, k;
  for (i=0; i < nbelem; i++)
    for (j=0; j < dimension; j++)
      elems(i,j) = facettes(i,j);

  construire_connectivite_som_elem(mesh.nb_sommets(), elems, som_elem, 0 /* include_virtual=0 */);
  assert(som_elem.get_nb_lists ()==mesh.nb_sommets());

  for (i = 0; i < nbelem; i++)
    if (num_compo[i] != -1)
      num_compo[i] = -2;

  int start_element = 0;
  int num_compo_courant = 0;
  ArrOfInt liste_elems;
  liste_elems.set_smart_resize(1);
  ArrOfInt tmp_liste;
  tmp_liste.set_smart_resize(1);
  //const int nb_som_elem=dimension-1;
  const int nb_som_elem=dimension;

  do
    {
      // Cherche le prochain element non attribue a une composante connexe
      while (start_element < nbelem && num_compo[start_element] >= -1)
        start_element++;
      if (start_element == nbelem)
        break;
      // Recherche des elements de la composante connexe a partir de cet element
      liste_elems.resize_array(1);
      liste_elems[0] = start_element;
      num_compo[start_element] = num_compo_courant;
      while (liste_elems.size_array() > 0)
        {
          tmp_liste.resize_array(0);
          const int liste_elems_size = liste_elems.size_array();
          for (int i_elem = 0; i_elem < liste_elems_size; i_elem++)
            {
              const int elem = liste_elems[i_elem];
              // Ajout des voisins non attribues de cet element dans la liste a
              // traiter a l'etape suivante
              for (int j2 = 0; j2 <  nb_som_elem; j2++)
                {
                  const int mon_som= facettes(elem, j2);
                  const int nb_elem_som= som_elem.get_list_size(mon_som);
                  for(k=0; k< nb_elem_som; k++)
                    {
                      const int voisin= som_elem(mon_som,k);
                      if (voisin >= 0)
                        {
                          const int num = num_compo[voisin];
                          if (num == -2)
                            {
                              num_compo[voisin] = num_compo_courant;
                              tmp_liste.append_array(voisin);
                            }
                        }
                    }
                }
            }
          liste_elems = tmp_liste;
        }
      num_compo_courant++;
    }
  while (1);
  // Renvoie le nombre de composantes connexes locales trouvees
  return num_compo_courant;
}

/*! @brief Recherche les composantes connexes d'un ensemble d'elements distribue sur tous les processeurs.
 *
 * Cette methode est parallele et doit etre appelee en
 *   meme temps sur tous les processeurs.
 *
 * @param (num_compo)
 * @param (nb_local_components)
 */
int compute_global_connex_components_FT(const Maillage_FT_Disc& mesh, ArrOfInt& num_compo, int nb_local_components)
{
  const int nbelem_tot = num_compo.size_array();
  int i;

  // Transformation des indices locaux de composantes connexes en un indice global
  // (on ajoute un decalage aux indices globaux avec mppartial_sum())
  const int decalage = mppartial_sum(nb_local_components);
  const int nb_total_components = Process::mp_sum(nb_local_components);
  for (i = 0; i < nbelem_tot; i++)
    if (num_compo[i] >= 0)
      num_compo[i] += decalage;

  // Pour trouver les correspondances entre un numero de composante locale et un
  // numero de la meme composante sur le processeur voisin, on cree une copie du
  // tableau num_compo sur laquelle on fait un echange_espace_virtuel(). Ainsi,
  // sur les cases virtuelles du tableau, on a dans num_compo le numero de la
  // composante locale et dans copie_compo le numero de cette meme composante sur
  // le processeur proprietaire de l'element. Donc ces deux numeros designent
  // la meme composante connexe.
  ArrOfInt copie_compo(num_compo);
  mesh.desc_facettes().echange_espace_virtuel(copie_compo);

  // Recherche des equivalences entre les numeros des composantes locales et
  // les numeros des composantes voisines. On construit un graphe dont les
  // liens relient les composantes equivalentes.
  // Tableau de marqeurs pour les equivalences deja trouvees.
  // Dimensions = nb composantes locales * nb_composantes total
  //  (pour ne pas prendre en compte la meme composante plusieurs fois).
  ArrOfBit markers(nb_local_components * nb_total_components);
  markers = 0;
  // Tableau de correspondances entre composantes connexes locales et distantes
  IntTab graph;
  graph.set_smart_resize(1);
  int graph_size = 0;
  // Parcours des elements virtuels uniquement
  for (i = 0; i < nbelem_tot; i++)
    {
      if (!mesh.facette_virtuelle(i))
        // Ne traiter que les facettes virtuelles:
        continue;

      int compo = num_compo[i];
      if (compo < 0)
        continue;
      int compo2 = copie_compo[i];
      // Index du couple compo2/compo dans le tableau markers
      // Le tableau num_compo ne doit contenir que des composantes locales:
      assert(compo >= decalage && compo - decalage < nb_local_components);
      // compo2 est forcement une composante distante.
      assert(compo2 < decalage || compo2 - decalage >= nb_local_components);
      const int index = (compo - decalage) * nb_total_components + compo2;
      if (!markers.testsetbit(index))
        {
          graph.resize(graph_size+1, 2);
          // On met le plus petit numero de composante en colonne 0:
          if (compo2 < compo)
            {
              int tmp = compo;
              compo = compo2;
              compo2 = tmp;
            }
          graph(graph_size, 0) = compo;
          graph(graph_size, 1) = compo2;
          graph_size++;
        }
    }

  ArrOfInt renum;
  if (Process::je_suis_maitre())
    {
      // Reception des portions de graphe des autres processeurs
      IntTab tmp;
      const int nproc = Process::nproc();
      int pe;
      for (pe = 1; pe < nproc; pe++)
        {
          recevoir(tmp, pe, 54 /* tag */);
          const int n2 = tmp.dimension(0);
          graph.resize(graph_size + n2, 2);
          for (int i2 = 0; i2 < n2; i2++)
            {
              graph(graph_size, 0) = tmp(i2, 0);
              graph(graph_size, 1) = tmp(i2, 1);
              graph_size++;
            }
        }
      // Calcul des composantes connexes du graphe
      renum.resize_array(nb_total_components);
      const int n = compute_graph_connex_components(graph, renum);
      Process::Journal() << "compute_global_connex_components: nb_components=" << n << finl;
      // Envoi des composantes connexes aux autres processeurs
    }
  else
    {
      // Envoi du graphe local au processeur 0
      envoyer(graph, 0, 54 /* tag */);
    }

  envoyer_broadcast(renum, 0 /* processeur source */);

  // Renumerotation des composantes dans num_compo
  for (i = 0; i < nbelem_tot; i++)
    {
      const int x = num_compo[i];
      if (x >= 0)
        {
          const int new_x = renum[x];
          num_compo[i] = new_x;
        }
    }
  // Verification: si on fait un echange espace virtuel,
  //  cela ne doit par changer le numero des composantes
  //  connexes !

  int nb_components = 0;
  if (renum.size_array() > 0)
    nb_components = max_array(renum) + 1;
  return nb_components;
}

