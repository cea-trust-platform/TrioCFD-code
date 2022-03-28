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
// File:        Descripteur_FT.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/16
//
//////////////////////////////////////////////////////////////////////////////

#include <Descripteur_FT.h>
#include <Comm_Group.h>
#include <TRUSTTabs.h>
#include <communications.h>
#include <Array_tools.h>
#include <MD_Vector_std.h>

Implemente_instanciable_sans_constructeur(Desc_Structure_FT,"Desc_Structure_FT",Objet_U);

Implemente_instanciable_sans_constructeur(Descripteur_FT,"Descripteur_FT",Objet_U);

Descripteur_FT::Descripteur_FT() :
  status_(OK)
{
  // Les tableaux sont initialises avec la taille maximale : nombre
  // total de processeurs
  const int n_proc = nproc();
  // On reserve la memoire pour le maximum...
  pe_voisins_.resize_array(n_proc);
  pe_voisins_.resize_array(0);
  elements_.dimensionner(n_proc);
  for (int i = 0; i < n_proc; i++)
    elements_[i].set_smart_resize(1);
}

Descripteur_FT& Descripteur_FT::operator=(const Descripteur_FT& src)
{
  const int n = nproc();
  // Pour ne pas casser l'attribut smart resize:
  for (int i = 0; i < n; i++)
    elements_[i] = src.elements_[i];
  pe_voisins_ = src.pe_voisins_;
  return *this;
}

Entree& Descripteur_FT::readOn(Entree& is)
{
  reset();
  Nom motlu;
  is >> motlu;
  if (motlu != que_suis_je())
    {
      Cerr << "Erreur dans Descripteur_FT::readOn\n";
      Cerr << " On attendait le motcle " << que_suis_je();
      Cerr << "\n On a trouve " << motlu << finl;
      exit();
    }
  is >> pe_voisins_;
  int i;
  const int nb_pe_voisins = pe_voisins_.size_array();
  for (i = 0; i < nb_pe_voisins; i++)
    {
      int pe = pe_voisins_[i];
      is >> elements_[pe];
    }
  return is;
}

Sortie& Descripteur_FT::printOn(Sortie& os) const
{
  const int nb_pe_voisins = pe_voisins_.size_array();
  os << que_suis_je() << finl;
  os << pe_voisins_;
  int i;
  for (i = 0; i < nb_pe_voisins; i++)
    {
      int pe = pe_voisins_[i];
      os << elements_[pe];
    }
  return os;
}

void Descripteur_FT::reset()
{
  int i;
  const int n = elements_.size();
  for (i = 0; i < n; i++)
    {
      elements_[i].resize_array(0);
    }
  pe_voisins_.resize_array(0);
  status_ = OK;
}

// Description: Ajoute l'element au tableau du PE_voisin.
// renvoie le rang de l'element ajoute dans le tableau d'elements du PE.
int Descripteur_FT::ajoute_element(int PE_voisin, int element)
{
  ArrOfInt& elements_voisins = elements_[PE_voisin];
  int n = elements_voisins.size_array();
  elements_voisins.append_array(element);
  status_ = BAD;
  return n;
}
// Description: Ajoute l'element au tableau du PE_voisin.
// renvoie la taiile du tableau d'elements du PE.
int Descripteur_FT::ajoute_elements(int PE_voisin, const ArrOfInt& NVelements)
{
  ArrOfInt& elements_voisins = elements_[PE_voisin];
  append_array_to_array(elements_voisins, NVelements);
  status_ = BAD;
  return elements_voisins.size_array();
}

// Description:
// Remplace la liste des elements par celle en parametre.
void Descripteur_FT::set_elements(int PE_voisin,
                                  const ArrOfInt& new_elements)
{
  ArrOfInt& elements_voisins = elements_[PE_voisin];
  elements_voisins = new_elements;
  status_ = BAD;
}

// Description: Renvoie "pas zero" si l'element est deja dans le
// descripteur pour le pe donne, 0 sinon.
int Descripteur_FT::contient_element(int pe, int element) const
{
  const ArrOfInt& elems = elements_[pe];
  const int n = elems.size_array();
  int i;
  for (i = 0; i < n ; i++)
    if (elems[i] == element) return 1;
  return 0;
}

// Description: Calcule la liste des PEs dont la liste d'elements est
// non vide, tries dans l'ordre croissant de numero de PE.
// Le statut passe a OK.
void Descripteur_FT::calcul_liste_pe_voisins()
{
  pe_voisins_.resize_array(0);
  const int n = elements_.size();
  // Boucle sur tous les PEs. Si la liste d'elements est non vide, on
  // l'ajoute a la fin du tableau pe_voisins_.
  for (int i = 0; i < n; i++)
    if (elements_[i].size_array() > 0)
      pe_voisins_.append_array(i);
  status_ = OK;
}

// Description:
// Pour chaque PE du descripteur, et pour 0 <= i < elements_[pe].size_array()
// on considere l'element suivante:
//  elem = elements_[pe][i]
// Si nouveau_pe[elem] >= 0, on retire elem du tableau elements_[pe]
// et on l'ajoute au tableau elements_retires.elements_[pe].
//
// Parametre: nouveau_pe
//  Signification: pour chaque element i (sommet ou facette), nouveau_pe[i] indique si
//                 l'element doit etre retire (nouveau_pe[i] >= 0) ou non
//                 (nouveau_pe[i] < 0) du descripteur
//  Contrainte:    la taille du tableau doit etre superieure au plus grand
//                 numero d'element stocke dans le descripteur.
// Parametre: elements_retires
//  Signification: descripteur contenant les elements retires.
void Descripteur_FT::retirer_elements(const ArrOfInt& nouveau_pe,
                                      Descripteur_FT&    elements_retires)
{
  static ArrOfIntFT new_elements;

  elements_retires.reset();

  int index_pe, i;
  int nb_voisins = pe_voisins_.size_array();
  // Boucle sur les PEs du descripteur.
  for (index_pe = 0; index_pe < nb_voisins; index_pe++)
    {

      int pe = pe_voisins_[index_pe];
      const ArrOfInt& elements_voisins = elements_[pe];
      ArrOfInt& elem_retires = elements_retires.elements_[pe];
      new_elements.resize_array(0);
      const int n = elements_voisins.size_array();
      // Boucle sur les elements du descripteur
      for (i = 0; i < n; i++)
        {
          const int element = elements_voisins[i];
          const int new_pe = nouveau_pe[element];
          if (new_pe < 0)
            new_elements.append_array(element);
          else
            elem_retires.append_array(element);
        }
      // Remplace la liste d'elements du descripteur par la nouvelle liste
      set_elements(pe, new_elements);
    }
  // Comme on a change les voisinages, on met a jour
  calcul_liste_pe_voisins();
  elements_retires.calcul_liste_pe_voisins();
}

//******************************************************************************
//******************************************************************************

// Constructeur par defaut

Desc_Structure_FT::Desc_Structure_FT() :
  status_md_(BAD)
{
}

Entree& Desc_Structure_FT::readOn(Entree& is)
{
  reset();
  Nom motlu;
  is >> motlu;
  if (motlu != que_suis_je())
    {
      Cerr << "Erreur dans Desc_Structure_FT::readOn\n";
      Cerr << " On attendait le motcle " << que_suis_je();
      Cerr << "\n On a trouve " << motlu << finl;
      exit();
    }
  is >> espace_distant_;
  is >> espace_virtuel_;

  check();
  return is;
}

Sortie& Desc_Structure_FT::printOn(Sortie& os) const
{
  os << que_suis_je() << finl;
  os << espace_distant_;
  os << espace_virtuel_;
  return os;
}

// Description : vide toutes les listes d'elements des descripteurs.
void Desc_Structure_FT::reset()
{
  espace_distant_.reset();
  espace_virtuel_.reset();
  status_md_ = BAD;
}

// Description: Renvoie une reference non-const a l'espace distant.
// Le statut du descripteur passe a BAD => il faudra recalculer le schema de com.
Descripteur_FT& Desc_Structure_FT::espace_distant()
{
  status_md_ = BAD;
  return espace_distant_;
}

// Description: Renvoie une reference const a l'espace distant.
// Le statut reste OK.
const Descripteur_FT& Desc_Structure_FT::espace_distant() const
{
  return espace_distant_;
}

// Description: Renvoie une reference non-const a l'espace virtuel.
// Le statut du descripteur passe a BAD => il faudra recalculer le schema de com.
Descripteur_FT& Desc_Structure_FT::espace_virtuel()
{
  status_md_ = BAD;
  return espace_virtuel_;
}

// Description: Renvoie une reference const a l'espace virtuel.
// Le statut reste OK.
const Descripteur_FT& Desc_Structure_FT::espace_virtuel() const
{
  return espace_virtuel_;
}


// Description:
// Verification de la coherence de la structure (graphe des voisins
// et taille des espaces virtuels et distants)
// On envoie tout au processeur 0.
int Desc_Structure_FT::check() const
{
  const int nb_proc = nproc();
  int i;

  // Pour etre sur d'arriver en meme temps ici:
  barrier();

  // Premiere etape du test: un element est
  // * soit dans le tableau d'elements_ virtuels d'un unique pe,
  // * soit dans un nombre quelconque de tableaux d'elements distants.
  {
    // Calcul du plus grand numero d'element mentionne dans les descripteurs
    int max_num_element = 0;
    {
      const ArrOfInt& pe_list = espace_distant_.pe_voisins();
      const int npe = pe_list.size_array();
      for (i = 0; i < npe; i++)
        {
          const int pe_voisin = pe_list[i];
          const int valeur_max = max_array(espace_distant_.elements(pe_voisin));
          if (valeur_max > max_num_element)
            max_num_element = valeur_max;
        }
    }
    {
      const ArrOfInt& pe_list = espace_virtuel_.pe_voisins();
      const int npe = pe_list.size_array();
      for (i = 0; i < npe; i++)
        {
          const int pe_voisin = pe_list[i];
          const int valeur_max = max_array(espace_virtuel_.elements(pe_voisin));
          if (valeur_max > max_num_element)
            max_num_element = valeur_max;
        }
    }
    ArrOfInt nb_tableaux_distants(max_num_element+1);
    nb_tableaux_distants = 0;
    // On compte le nombre de fois ou les elements apparaissent dans l'espace distant
    {
      const ArrOfInt& pe_list = espace_distant_.pe_voisins();
      const int npe = pe_list.size_array();
      for (i = 0; i < npe; i++)
        {
          const int pe_voisin = pe_list[i];
          const ArrOfInt& elements = espace_distant_.elements(pe_voisin);
          const int n = elements.size_array();
          for (int j = 0; j < n; j++)
            {
              const int elem = elements[j];
              ++nb_tableaux_distants[elem];
            }
        }
    }
    // On compte le nombre de fois ou les elements apparaissent dans l'espace virtuel
    {
      const ArrOfInt& pe_list = espace_virtuel_.pe_voisins();
      const int npe = pe_list.size_array();
      for (i = 0; i < npe; i++)
        {
          const int pe_voisin = pe_list[i];
          const ArrOfInt& elements = espace_virtuel_.elements(pe_voisin);
          const int n = elements.size_array();
          for (int j = 0; j < n; j++)
            {
              const int elem = elements[j];
              const int x = nb_tableaux_distants[elem];
              if (x != 0)
                {
                  // Soit l'element est dans un espace distant, soit il est dans un autre
                  // espace virtuel. Ca va pas.
                  Cerr << "(PE " << Process::me();
                  Cerr << ") l'element " << elem;
                  Cerr << " apparait dans ";
                  if (x > 0)
                    Cerr << " des espaces distants et un espace virtuel" << finl;
                  else
                    Cerr << " plusieurs espaces virtuels" << finl;
                  exit();
                }
              else
                {
                  nb_tableaux_distants[elem] = -1;
                }
            }
        }
    }
  }

  VECT(IntTab) procs_distants(nb_proc);
  VECT(IntTab) procs_virtuels(nb_proc);

  // Tout le monde envoie ses descripteurs au processeur 0
  {
    // Construction d'un tableau : pour chaque pe voisin, numero du pe (colonne 0)
    // et nombre d'elements a envoyer/recevoir (colonne 1)
    int moi = me();
    IntTab& procs_dist = procs_distants[moi];
    const ArrOfInt& pe_distants = espace_distant_.pe_voisins();
    int nb_procs_distants = pe_distants.size_array();
    procs_dist.resize(nb_procs_distants, 2);
    for (i = 0; i < nb_procs_distants; i++)
      {
        int pe_voisin = pe_distants[i];
        procs_dist(i, 0) = pe_voisin;
        procs_dist(i, 1) = espace_distant_.elements(pe_voisin).size_array();
      }

    IntTab& procs_virt = procs_virtuels[moi];
    const ArrOfInt& pe_virtuels = espace_virtuel_.pe_voisins();
    int nb_procs_virtuels = pe_virtuels.size_array();
    procs_virt.resize(nb_procs_virtuels, 2);
    for (i = 0; i < nb_procs_virtuels; i++)
      {
        int pe_voisin = pe_virtuels[i];
        procs_virt(i, 0) = pe_voisin;
        procs_virt(i, 1) = espace_virtuel_.elements(pe_voisin).size_array();
      }
  }

  if (!je_suis_maitre())
    {
      int moi = me();
      envoyer(procs_distants[moi], 0, 0);
      envoyer(procs_virtuels[moi], 0, 1);
    }
  else
    {
      for (i = 1; i < nb_proc; i++)
        {
          recevoir(procs_distants[i], i, 0);
          recevoir(procs_virtuels[i], i, 1);
        }
    }
  // Le maitre reste tard le soir pour finir le travail
  if (me() != 0)
    return 1;

  // On verifie que le nombre d'elements distants correspond bien
  // au nombre d'elements virtuels du processeur voisin...

  // Boucle sur les PEs,
  //   et pour chaque PE, boucle sur les voisins de l'espace distant.
  //   On cherche ce voisin dans l'espace virtuel du PE...
  {
    // Nombre de voisins deja parcourus dans les espaces virtuels...
    ArrOfInt curseurs(nb_proc);
    curseurs = 0;

    for (i = 0; i < nb_proc; i++)
      {

        // La liste des voisins du proc i par l'espace distant
        const IntTab& procs_dist = procs_distants[i];
        int nb_pe_voisins = procs_dist.dimension(0);

        for (int j = 0; j < nb_pe_voisins; j++)
          {

            int pe_voisin = procs_dist(j, 0);
            int nb_elements = procs_dist(j, 1);
            // On verifie que l'espace virtuel du pe_voisin correspond
            const int index = curseurs[pe_voisin];
            curseurs[pe_voisin]++;
            const IntTab& procs_virt = procs_virtuels[pe_voisin];
            if (index >= procs_virt.dimension(0) || procs_virt(index, 0) != i)
              {
                Cerr << "Descripteur_FT::check : erreur sur le graphe de voisinage.\n";
                Cerr << " Les espaces distants et virtuels ne se correspondent pas.\n";
                Cerr << pe_voisin << finl;
                exit();
              }
            if (procs_virt(index, 1) != nb_elements)
              {
                Cerr << "Descripteur_FT::check : erreur sur une liste d'elements.\n";
                Cerr << " Espace distant du PE " << i << " : PE voisin " << pe_voisin;
                Cerr << ", nombre d'elements " << nb_elements;
                Cerr << "\n Espace virtuel du PE " << pe_voisin << " : PE voisin " << i;
                Cerr << ", nombre d'elements " << procs_virt(index, 1) << finl;
                exit();
              }
          }
      }
    // On verifie que tous les espaces virtuels ont ete balayes
    for (i = 0; i < nb_proc; i++)
      {
        const int n1 = curseurs[i];
        const int n2 = procs_virtuels[i].dimension(0);
        if (n1 < n2)
          {
            Cerr << "Descripteur_FT::check : erreur sur le graphe de voisinage.\n";
            Cerr << " Les espaces distants et virtuels ne se correspondent pas." << finl;
            exit();
          }
      }
  }
  return 1;
}

// Description:
// Correction des espaces distants et virtuels lors d'un changement de proprietaire
// (noeud ou facette). Les operations sont les suivantes :
// - sur l'ancien proprietaire de l'element, on retire l'element de l'espace distant.
// - sur le nouveau proprietaire, on ajoute l'element a l'espace distant, pour tous
//   les PEs ou cet element est virtuel.
// - sur les PEs ou cet element est virtuel, on retire l'element de la liste
//   d'elements virtuels de l'ancien PE et on l'ajoute a la fin de la liste d'elements
//   virtuels du nouveau proprietaire.
// L'ordre d'ajout des elements dans les espaces distants et virtuels est identique
// pour que les espaces soient en correspondance.
//
// Parametre: nouveau_pe
// Signification: tableau de taille nb_sommets ou nb_facettes,
//            qui contient le numero du nouveau pe proprietaire,
//            ou -1 si l'element ne change pas de main.
//            L'espace virtuel de ce tableau doit etre a jour
//            (echange_espace_virtuel realise).

void Desc_Structure_FT::echanger_elements(const ArrOfInt& nouveau_pe)
{
  static Descripteur_FT elements_distants_retires;
  static Descripteur_FT elements_virtuels_retires;
  static Descripteur_FT elements_recus;
  static Descripteur_FT elements_envoyes;
  static ArrOfInt cles;

  const int moi = me();

  cles.resize_array(nouveau_pe.size_array());
  cles = -1;

  // On retire des espaces distants et virtuels les elements qui ont change de main.
  // Chaque element est retire a la fois de l'espace distant et de l'espace virtuel,
  // ces deux espaces restent donc en correspondance.
  // De plus elements_distants_retires et elements_virtuels_retires sont aussi
  // en correspondance.

  espace_distant_.retirer_elements(nouveau_pe, elements_distants_retires);
  espace_virtuel_.retirer_elements(nouveau_pe, elements_virtuels_retires);

  // Traitement des elements_virtuels_retires. Pour chaque element:
  // * soit le nouveau pe de l'element n'est pas moi : on cree un element virtuel
  // * soit c'est moi qui le recoit : on ajoute l'element a la liste d'elements recus
  // A la fin de ce traitement, elements_recus contient la liste de tous les elements
  // qu'il faut ajouter a des espaces distants.
  // On a ajoute a l'espace virtuel les elements virtuels->virtuels
  // (qui etaient virtuels avant et qui sont virtuels apres).
  {
    const ArrOfInt& pe_voisins = elements_virtuels_retires.pe_voisins();
    const int nb_voisins = pe_voisins.size_array();
    for (int index_pe = 0; index_pe < nb_voisins; index_pe++)
      {
        const int pe_origine = pe_voisins[index_pe];
        const ArrOfInt& elements = elements_virtuels_retires.elements(pe_origine);
        const int nb_elements = elements.size_array();
        for (int i = 0; i < nb_elements; i++)
          {
            const int element = elements[i];
            int pe_destination = nouveau_pe[element];
            if (pe_destination == moi)
              {
                elements_recus.ajoute_element(pe_origine, element);
              }
            else
              {
                espace_virtuel_.ajoute_element(pe_destination, element);
              }
          }
      }
  }
  elements_recus.calcul_liste_pe_voisins();

  // Remplissage d'un descripteur contenant, pour chaque PE, la liste des
  // elements qu'on lui envoie (elements_envoyes). Ces elements sont les
  // elements_distants_retires dont le nouveau_pe est le pe_distant.
  // On remplit en meme temps la "cle" qui represente l'element transmis
  // au pe destination. Cette cle est le rang de l'element dans
  // elements_envoyes.elements(pe_destination).
  // On ajoute a l'espace virtuel les elements   reel->virtuel
  // (qui etaient reels avant et qui sont virtuels maintenant).
  {
    const ArrOfInt& pe_voisins = elements_distants_retires.pe_voisins();
    const int nb_voisins = pe_voisins.size_array();
    for (int index_pe = 0; index_pe < nb_voisins; index_pe++)
      {
        const int pe_distant = pe_voisins[index_pe];
        const ArrOfInt& elements = elements_distants_retires.elements(pe_distant);
        const int nb_elements = elements.size_array();
        for (int i = 0; i < nb_elements; i++)
          {
            const int element = elements[i];
            const int pe_destination = nouveau_pe[element];
            if (pe_destination == pe_distant)
              {
                const int cle = elements_envoyes.ajoute_element(pe_distant, element);
                assert(cles[element] < 0);
                cles[element] = cle;
                // Cet ajout est symetrique de l'ajout (*) a la fin de la fonction:
                espace_virtuel_.ajoute_element(pe_distant, element);
              }
          }
      }
  }
  elements_envoyes.calcul_liste_pe_voisins();

  // Traitement des espaces distants: il faut envoyer au nouveau proprietaire
  // les numeros des PEs voisins pour chaque noeud qu'on lui envoie.
  // Le noeud est repere par la "cle", car les descripteurs elements_envoyes
  // et elements_recus sont en correspondance.
  schema_comm_.set_send_recv_pe_list(elements_envoyes.pe_voisins(),
                                     elements_recus.pe_voisins());

  schema_comm_.begin_comm();
  {
    const ArrOfInt& pe_voisins = elements_distants_retires.pe_voisins();
    const int nb_voisins = pe_voisins.size_array();
    for (int index_pe = 0; index_pe < nb_voisins; index_pe++)
      {
        const int pe_distant = pe_voisins[index_pe];
        const ArrOfInt& elements = elements_distants_retires.elements(pe_distant);
        const int nb_elements = elements.size_array();
        for (int i = 0; i < nb_elements; i++)
          {
            const int element = elements[i];
            const int pe_destination = nouveau_pe[element];
            if (pe_destination != pe_distant)
              {
                const int cle = cles[element];
                assert(cle >= 0);
                schema_comm_.send_buffer(pe_destination) << cle << pe_distant;
              }
          }
      }
  }
  schema_comm_.echange_taille_et_messages();
  // Recuperation des espaces distants :
  // Si le processeur distant est moi, alors je mets le noeud dans l'espace
  // distant du PE d'origine.
  // On ajoute a l'espace distant les elements virtuel->virtuel
  {
    const ArrOfInt& pe_voisins = elements_recus.pe_voisins();
    const int nb_voisins = pe_voisins.size_array();
    for (int index_pe = 0; index_pe < nb_voisins; index_pe++)
      {
        const int pe_origine = pe_voisins[index_pe];
        Entree& recv_buffer = schema_comm_.recv_buffer(pe_origine);
        const ArrOfInt& elements = elements_recus.elements(pe_origine);
        while (1)
          {
            int cle, pe_distant;
            recv_buffer >> cle >> pe_distant;
            if (recv_buffer.eof())
              break;
            const int element = elements[cle];
            espace_distant_.ajoute_element(pe_distant, element);
          }
      }
  }
  schema_comm_.end_comm();

  // Ajout (*) des elements virtuels->distants a l'espace distant
  // (qui etaient virtuels et qui sont reels maintenant)
  {
    const ArrOfInt& pe_voisins = elements_virtuels_retires.pe_voisins();
    const int nb_voisins = pe_voisins.size_array();
    for (int index_pe = 0; index_pe < nb_voisins; index_pe++)
      {
        const int pe_distant = pe_voisins[index_pe];
        const ArrOfInt& elements = elements_virtuels_retires.elements(pe_distant);
        const int nb_elements = elements.size_array();
        for (int i = 0; i < nb_elements; i++)
          {
            const int element = elements[i];
            int pe_destination = nouveau_pe[element];
            if (pe_destination == moi)
              espace_distant_.ajoute_element(pe_distant, element);
          }
      }
  }

  espace_distant_.calcul_liste_pe_voisins();
  espace_virtuel_.calcul_liste_pe_voisins();
  if (Comm_Group::check_enabled()) check();

  status_md_ = BAD;
  calcul_schema_comm(nouveau_pe.size_array());
  cles.resize_array(0); // Libere la memoire du tableau
}

// Description: remplit le tableau qui donne pour chaque element le
// numero du pe proprietaire. Le tableau fourni doit avoir la bonne
// taille (nombre d'elements = sommets ou facettes par ex.)
// Ce tableau est rempli en utilisant l'espace virtuel (tous les elements
// sont a moi, sauf ceux qui sont dans l'espace virtuel).

void Desc_Structure_FT::remplir_element_pe(ArrOfInt& element_pe) const
{
  // Tout est a moi ...
  element_pe = me();

  // ... sauf ce qui est dans les espaces virtuels
  const ArrOfInt& pe_voisins = espace_virtuel_.pe_voisins();
  const int nb_voisins = pe_voisins.size_array();
  for (int i = 0; i < nb_voisins; i++)
    {
      const int pe = pe_voisins[i];
      const ArrOfInt& elements = espace_virtuel_.elements(pe);
      const int nb_elements = elements.size_array();
      for (int j = 0; j < nb_elements; j++)
        {
          int element = elements[j];
          element_pe[element] = pe;
        }
    }
}

void Desc_Structure_FT::calcul_schema_comm(const int nb_items_tot)
{
  const int n = nproc();
  if (pe_voisins_.size_array() != n)
    {
      pe_voisins_.resize_array(n);
      blocs_to_recv_.dimensionner(n);
      for (int i = 0; i < n; i++)
        pe_voisins_[i] = i;
    }

  const VECT(ArrOfInt) & items_to_send_ = espace_distant().all_elements();
  const VECT(ArrOfInt) & items_to_recv_ = espace_virtuel().all_elements();

  MD_Vector_std md(nb_items_tot,
                   -1, /* nb_items reels */
                   pe_voisins_,
                   items_to_send_,
                   items_to_recv_,
                   blocs_to_recv_);

  md_vector_.copy(md);

  const ArrOfInt& send_procs = espace_distant().pe_voisins();
  const ArrOfInt& recv_procs = espace_virtuel().pe_voisins();
  schema_comm_.set_send_recv_pe_list(send_procs, recv_procs);
  schema_comm_inverse_.set_send_recv_pe_list(recv_procs, send_procs);

  status_md_ = OK;
}

const Schema_Comm_FT& Desc_Structure_FT::schema_comm() const
{
  assert(status_md_ == OK);
  return schema_comm_;
}

const Schema_Comm_FT& Desc_Structure_FT::schema_comm_inverse() const
{
  assert(status_md_ == OK);
  return schema_comm_inverse_;
}

void Desc_Structure_FT::echange_espace_virtuel(ArrOfInt& x) const
{
  collecter_espace_virtuel(x, MD_Vector_tools::ECHANGE_EV);
}

void Desc_Structure_FT::echange_espace_virtuel(ArrOfDouble& x) const
{
  collecter_espace_virtuel(x, MD_Vector_tools::ECHANGE_EV);
}

void Desc_Structure_FT::collecter_espace_virtuel(ArrOfInt& x, MD_Vector_tools::Operations_echange op) const
{
  assert(status_md_ == OK);
  IntVect& y = tmp_intvect_;
  IntVect* intV = dynamic_cast<IntVect*>(&x);

  if (intV) y.ref(*intV);
  else y.ref_array(x);

  y.set_md_vector(md_vector_);
  MD_Vector_tools::echange_espace_virtuel(y, op);
  y.reset();
}

void Desc_Structure_FT::collecter_espace_virtuel(ArrOfDouble& x, MD_Vector_tools::Operations_echange op) const
{
  assert(status_md_ == OK);
  DoubleVect& y = tmp_doublevect_;
  DoubleVect* doubleV = dynamic_cast<DoubleVect*>(&x);

  if (doubleV) y.ref(*doubleV);
  else y.ref_array(x);

  y.set_md_vector(md_vector_);
  MD_Vector_tools::echange_espace_virtuel(y, op);
  y.reset();
}

const MD_Vector& Desc_Structure_FT::get_md_vector() const
{
  assert(status_md_ == OK);
  return md_vector_;
}
