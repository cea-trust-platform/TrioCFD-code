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
// File:        Sauvegarde_Reprise_Maillage_FT.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////
#include <Sauvegarde_Reprise_Maillage_FT.h>
#include <Maillage_FT_Disc.h>
#include <Zone_VF.h>
#include <communications.h>
#include <Array_tools.h>
#include <Domaine.h>

void ecrire_tableau(Sortie& os, const DoubleTab& tab)
{
  const int dim0 = tab.dimension(0);
  const int dimtot = Process::mp_sum(dim0);
  if (Process::je_suis_maitre())
    os << dimtot << space << tab.dimension(1) << finl;
  os.put(tab.addr(), tab.size_array());
  os.syncfile();
}

void ecrire_tableau(Sortie& os, const IntTab& tab)
{
  const int dim0 = tab.dimension(0);
  const int dimtot = Process::mp_sum(dim0);
  if (Process::je_suis_maitre())
    os << dimtot << space << tab.dimension(1) << finl;
  os.put(tab.addr(), tab.size_array());
  os.syncfile();
}

void ecrire_tableau(Sortie& os, const ArrOfInt& tab)
{
  const int dim0 = tab.size_array();
  const int dimtot = Process::mp_sum(dim0);
  if (Process::je_suis_maitre())
    os << dimtot << finl;
  os.put(tab.addr(), tab.size_array());
  os.syncfile();
}

void Sauvegarde_Reprise_Maillage_FT::ecrire_xyz(const Maillage_FT_Disc& mesh, const Zone_VF& zone_vf, Sortie& fichier)
{
  if (Process::je_suis_maitre())
    fichier << mesh.temps_physique_ << finl;
  // ========== ETAPE 1 : SOMMETS ===========
  // Sauvegarde des sommets et valeurs associees.
  // Pour chaque sommet, ecriture de :
  //  - coordonnees des coordonnees du centre de gravite des elements contenant le sommet
  //  - coordonnees des sommets
  //  - si le sommet est sur une face de bord, indice de la face dans l'element, sinon -1
  const int nb_sommets = mesh.nb_sommets();
  DoubleTab coord_elem;
  DoubleTab coord_som;
  ArrOfInt  num_face_bord;
  ArrOfIntFT indice_global_sommet(nb_sommets);
  coord_elem.set_smart_resize(1);
  coord_som.set_smart_resize(1);
  num_face_bord.set_smart_resize(1);
  // On surdimensionne les tableaux, on fera resize a la fin
  const int dim = mesh.sommets_.dimension(1);
  coord_elem.resize(nb_sommets, dim);
  coord_som.resize(nb_sommets, dim);
  num_face_bord.resize_array(nb_sommets);
  const DoubleTab& xp = zone_vf.xp();
  // Compteur de sommets reels
  int nb_sommets_reels = 0;
  int i, j;
  for (i = 0; i < nb_sommets; i++)
    {
      if (mesh.sommet_virtuel(i))
        {
          indice_global_sommet[i] = -1;
          continue;
        }
      const int indice_element = mesh.sommet_elem_[i];
      for (j = 0; j < dim; j++)
        coord_elem(nb_sommets_reels, j) = xp(indice_element, j);
      for (j = 0; j < dim; j++)
        coord_som(nb_sommets_reels, j) = mesh.sommets_(i, j);
      const int face_bord = mesh.sommet_face_bord_[i];
      int i_face_bord = -1;
      if (face_bord >= 0)
        {
          i_face_bord = zone_vf.numero_face_local(face_bord, indice_element);
          assert(face_bord == zone_vf.elem_faces(indice_element, i_face_bord));
        }
      num_face_bord[nb_sommets_reels] = i_face_bord;
      indice_global_sommet[i] = nb_sommets_reels;
      nb_sommets_reels++;
    }
  coord_elem.resize(nb_sommets_reels, dim);
  coord_som.resize(nb_sommets_reels, dim);
  num_face_bord.resize_array(nb_sommets_reels);
  // Calcul de l'offset a ajouter pour avoir l'indice global du sommet:
  const int offset = mppartial_sum(nb_sommets_reels);
  for (i = 0; i < nb_sommets; i++)
    if (indice_global_sommet[i] >= 0)
      indice_global_sommet[i] += offset;
  mesh.desc_sommets().echange_espace_virtuel(indice_global_sommet);

  ecrire_tableau(fichier, coord_elem);
  ecrire_tableau(fichier, coord_som);
  ecrire_tableau(fichier, num_face_bord);

  // ============== ETAPE 2 : FACETTES =============
  // Sauvegarde des facettes : conversion des indices de sommets des facettes en indices globaux
  // et sauvegarde des facettes reeles uniquement.
  IntTab facettes_to_save;
  const int nb_facettes = mesh.nb_facettes();
  facettes_to_save.set_smart_resize(1);
  facettes_to_save.resize(nb_facettes, dim);
  int nb_facettes_reeles = 0;
  for (i = 0; i < nb_facettes; i++)
    {
      if (mesh.facette_virtuelle(i))
        continue;
      for (j = 0; j < dim; j++)
        {
          const int i_som = mesh.facettes_(i, j);
          const int i_som_global = indice_global_sommet[i_som];
          assert(i_som_global >= 0);
          facettes_to_save(nb_facettes_reeles, j) = i_som_global;
        }
      nb_facettes_reeles++;
    }
  facettes_to_save.resize(nb_facettes_reeles, dim);
  ecrire_tableau(fichier, facettes_to_save);
}

// Description: Remplissage du maillage mesh a partir d'un fichier de sauvegarde xyz (fichier)
//  ou d'un domaine source (domaine_src pour importation d'une cao)
//  zone_vf est la zone support du maillage (pour calcul de l'indice de l'element contenant chaque
//  sommet du mesh. zone_vf peut etre nul dans le cas d'une reprise xyz au cours de l'etape 'avancer'.
//  Dans ce cas, on lit les coordonnees et les elements mais on ne remplit pas mesh.
void Sauvegarde_Reprise_Maillage_FT::lire_xyz(Maillage_FT_Disc& mesh,
                                              const Zone_VF * zone_vf,
                                              Entree * fichier,
                                              const Domaine * domaine_src)
{
  // On doit donner un fichier OU un domaine_src, mais pas les deux !
  assert((fichier && (!domaine_src)) || ((!fichier) && domaine_src));
  const int moi = Process::me();
  const int nproc = Process::nproc();
  mesh.statut_ = Maillage_FT_Disc::MINIMAL;
  if (fichier)
    (*fichier) >> mesh.temps_physique_;
  // ============== ETAPE 1 : lecture des sommets reels ===============
  // On lit les sommets, et on ne conserve que les sommets dont on est proprietaire
  const int dim = Objet_U::dimension;
  mesh.sommets_.resize(0, dim);

  int nb_som_tot, dummy_int;
  if (fichier)
    {
      // Reprise xyz : on lit le nombre de sommets
      (*fichier) >> nb_som_tot >> dummy_int;
      assert(dummy_int == dim);
    }
  else
    {
      // lecture cao : on utilise le domaine
      nb_som_tot = domaine_src->nb_som();
      if (domaine_src->coord_sommets().dimension(1) != dim)
        {
          Cerr << "Error in Sauvegarde_Reprise_Maillage_FT::lire_xyz : source domain has bad dimension" << finl;
          Process::exit();
        }
    }

  const int chunk_size = 65536; // Lire les sommets par blocs
  // Pour chaque sommet conserve, quel est son indice global
  ArrOfInt indice_global_sommet;
  indice_global_sommet.set_smart_resize(1);
  const Zone * const zone = (zone_vf ? &zone_vf->zone() : 0);
  // Lecture des indices des elements contenant le sommet
  int erreur_sommets_exterieurs = 0;
  {
    DoubleTab tmp;
    // tmp2 : pour chaque sommet, -1 si je ne le garde pas, sinon numero de l'element qui contient le sommet
    ArrOfInt tmp2;
    // tmp3 : -1 ou numero du processeur qui garde le sommet.
    ArrOfInt tmp3;
    int i = 0;
    const int nb_elements_reels = (zone_vf ? zone->nb_elem() : 0);
    while (i < nb_som_tot)
      {
        const int n_to_read = std::min(chunk_size, nb_som_tot - i);
        tmp.resize(n_to_read, dim);
        tmp2.resize_array(n_to_read);
        tmp3.resize_array(n_to_read);
        if (fichier)
          {
            // Lecture des coordonnees des centres des elements dans le fichier
            fichier->get(tmp.addr(), tmp.size_array());
          }
        else
          {
            // Copie des coordonnees des sommets dans le domaine
            const DoubleTab& coord = domaine_src->coord_sommets();
            for (int j = 0; j < n_to_read; j++)
              for (int k = 0; k < dim; k++)
                tmp(j, k) = coord(i+j, k);
          }
        if (zone_vf)
          {
            zone->chercher_elements(tmp, tmp2);
            // Ne pas tenir compte des sommets dans les elements virtuels
            for (int j = 0; j < n_to_read; j++)
              if (tmp2[j] >= nb_elements_reels)
                tmp2[j] = -1;
            // Verifier qu'un processeur et un seul prend la main.
            // Le premier processeur initialise tmp3 en mettant son numero
            // pour les sommets qu'il possede et -1 pour les autres.
            // Puis la liste est passee au processeur suivant qui la complete.
            if (moi == 0)
              {
                // Preparer la premiere liste
                for (int j = 0; j < n_to_read; j++)
                  if (tmp2[j] >= 0)
                    tmp3[j] = moi;
                  else
                    tmp3[j] = -1;
              }
            else
              {
                // Recevoir une liste d'attribution du processeur precedent
                // et la completer
                recevoir(tmp3, moi - 1, moi, 0 /* canal */);
                for (int j = 0; j < n_to_read; j++)
                  if (tmp3[j] >= 0) // Le sommet a ete pris par un autre processeur
                    tmp2[j] = -1;
                  else if (tmp2[j] >= 0) // Le sommet n'a pas ete pris et est chez moi
                    tmp3[j] = moi;
              }
            if (moi < nproc - 1)
              {
                envoyer(tmp3, moi, moi + 1, 0 /* canal */);
              }
            else
              {
                // Le dernier processeur prend les sommets qui restent
                for (int j = 0; j < n_to_read; j++)
                  if (tmp3[j] < 0)
                    {
                      tmp3[j] = moi;
                      // Pour l'instant on ne sait pas traiter des maillages qui depassent du
                      // domaine:
                      erreur_sommets_exterieurs=1;
                    }
              }

            for (int j = 0; j < n_to_read; j++)
              {
                // i+j est le rang de ce sommet dans le fichier lu
                if (tmp3[j] == moi)
                  {
                    indice_global_sommet.append_array(i+j);
                    mesh.sommet_elem_.append_array(tmp2[j]);
                  }
              }
          }
        i += n_to_read;
      }
  }

// if (Process::mp_sum(erreur_sommets_exterieurs))
// GF bloque sinon dans avancer en // de plis c'est inutil
  if (erreur_sommets_exterieurs)
    {
      Cerr << "Erreur a la lecture d'un maillage front-tracking :\n"
           << " Certains sommets sont a l'exterieur du maillage eulerien.\n"
           << " On ne sait pas encore traiter ce cas." << finl;
      Process::barrier();
      Process::exit();
    }

  // Lecture des coordonnees des sommets
  if (fichier)
    {
      int n;
      (*fichier) >> n >> dummy_int;
    }
  const int nb_sommets_locaux = mesh.sommet_elem_.size_array();
  {
    int i = 0;
    int i_sommet = 0;
    DoubleTab tmp;
    mesh.sommets_.resize(nb_sommets_locaux, dim);
    while (i < nb_som_tot)
      {
        const int n_to_read = std::min(chunk_size, nb_som_tot - i);
        tmp.resize(n_to_read, dim);
        if (fichier)
          {
            fichier->get(tmp.addr(), tmp.size_array());
          }
        else
          {
            // Copie des coordonnees des sommets dans le domaine
            const DoubleTab& coord = domaine_src->coord_sommets();
            for (int j = 0; j < n_to_read; j++)
              for (int k = 0; k < dim; k++)
                tmp(j, k) = coord(i+j, k);
          }
        if (zone_vf)
          {
            for (int j = 0; j < n_to_read; j++)
              {
                const int i_glob = i+j;
                // Ce processeur possede-t-il ce sommet ?
                if (i_sommet >= nb_sommets_locaux || i_glob != indice_global_sommet[i_sommet])
                  continue;
                for (int k = 0; k < dim; k++)
                  mesh.sommets_(i_sommet, k) = tmp(j, k);
                i_sommet++;
              }
          }
        i += n_to_read;
      }
    assert(zone_vf == 0 || i_sommet == mesh.nb_sommets());
  }
  // Lecture des faces de bord
  if (fichier)
    {
      int dummy;
      (*fichier) >> dummy;
      int i = 0;
      int i_sommet = 0;
      ArrOfInt tmp;
      mesh.sommet_face_bord_.resize_array(nb_sommets_locaux);
      const IntTab * elem_faces = (zone_vf ? &zone_vf->elem_faces() : 0);
      const int nb_faces_bord = (zone ? zone->nb_faces_frontiere() : 0);
      while (i < nb_som_tot)
        {
          const int n_to_read = std::min(chunk_size, nb_som_tot - i);
          tmp.resize_array(n_to_read);
          fichier->get(tmp.addr(), tmp.size_array());
          if (zone_vf)
            {
              for (int j = 0; j < n_to_read; j++)
                {
                  const int i_glob = i+j;
                  // Ce processeur possede-t-il ce sommet ?
                  if (i_sommet >= nb_sommets_locaux || i_glob != indice_global_sommet[i_sommet])
                    continue;
                  int i_face = -1;
                  const int i_face_local = tmp[j];
                  if (i_face_local >= 0)
                    {
                      const int elem = mesh.sommet_elem_[i_sommet];
                      i_face = (*elem_faces)(elem, i_face_local);
                      if (i_face >= nb_faces_bord)
                        {
                          // La face trouvee n'est pas une face de bord
                          // Seule explication possible: les faces des elements
                          // ont ete renumerotees.
                          Cerr << "Erreur a la relecture des faces de bord" << finl;
                          Process::exit();
                        }
                    }
                  mesh.sommet_face_bord_[i_sommet] = i_face;
                  i_sommet++;
                }
            }
          i += n_to_read;
        }
      assert(zone_vf == 0 || i_sommet == mesh.nb_sommets());
    }
  else
    {
      // Pas de sommet au bord pour l'instant
      mesh.sommet_face_bord_.resize_array(nb_sommets_locaux);
      mesh.sommet_face_bord_ = -1;
    }
  // Lecture des facettes
  // On lit dans itmp et on copie dans facettes les facettes
  // qui appartiennent au processeur local.
  IntTab facettes(0, dim);
  int nb_faces_tot;
  if (fichier)
    {
      (*fichier) >> nb_faces_tot >> dummy_int;
      assert(dummy_int == dim);
    }
  else
    {
      // lecture cao : on utilise le domaine
      nb_faces_tot = domaine_src->zone(0).nb_elem();
      if (domaine_src->zone(0).les_elems().dimension(1) != dim)
        {
          Cerr << "Error in Sauvegarde_Reprise_Maillage_FT::lire_xyz : source domain elements have bad dimension" << finl;
          Process::exit();
        }
    }
  {
    IntTab tmp;
    int i = 0;
    facettes.set_smart_resize(1);
    while (i < nb_faces_tot)
      {
        const int n_to_read = std::min(chunk_size, nb_faces_tot - i);
        tmp.resize(n_to_read, dim);
        if (fichier)
          fichier->get(tmp.addr(), tmp.size_array());
        else
          {
            // Copie des elements dans le domaine
            const IntTab& elems = domaine_src->zone(0).les_elems();
            for (int j = 0; j < n_to_read; j++)
              for (int k = 0; k < dim; k++)
                tmp(j, k) = elems(i+j, k);
          }
        if (zone_vf)
          {
            for (int j = 0; j < n_to_read; j++)
              {
                // Indice global du premier sommet de la facette.
                const int i_som = tmp(j, 0);
                // Le processeur local possede-t-il ce sommet ?
                const int i_som_loc = array_bsearch(indice_global_sommet, i_som);
                if (i_som_loc >= 0)
                  {
                    // oui
                    const int n = facettes.dimension(0);
                    facettes.resize(n+1, dim);
                    for (int k = 0; k < dim; k++)
                      facettes(n, k) = tmp(j, k);
                  }
              }
          }
        i += n_to_read;
      }
  }
  if (!zone_vf)
    return;

  {
    const int n = mesh.nb_sommets();
    mesh.sommet_PE_owner_.resize_array(n);
    mesh.sommet_num_owner_.resize_array(n);
    mesh.drapeaux_sommets_.resize_array(n);
    for (int i = 0; i < n; i++)
      {
        mesh.sommet_PE_owner_[i] = moi;
        mesh.sommet_num_owner_[i] = i;
        mesh.drapeaux_sommets_[i] = 0;
      }
    mesh.desc_sommets_.calcul_schema_comm(mesh.nb_sommets());
    mesh.desc_facettes_.calcul_schema_comm(mesh.nb_facettes());
  }
  // Recuperer les sommets manquants chez les autres processeurs:
  // Chaque processeur envoie au processeur suivant les indices globaux des
  // sommets qui lui manquent. Le processeur traite la liste qu'il
  // a recu et envoie les sommets restants au processeur suivant.
  {
    // Liste d'indices de sommets locaux a envoyer
    ArrOfIntFT liste_sommets_to_send;
    // Pour chaque sommet de la liste precedente, numero du pe destination
    ArrOfIntFT liste_pe_to_send;
    // Liste d'indices des indices de sommets globaux a transmettre
    // (au depart, la liste dont le processeur local a besoin, puis on
    // fait passer la liste au processeur suivant en chaine).
    ArrOfInt liste_sommets_globaux;
    liste_sommets_globaux.set_smart_resize(1);
    // Recherche des sommets inconnus parmi les sommets des facettes
    {
      // Tableau facettes, vu comme un tableau monodimensionnel:
      const ArrOfInt& liste_tous_sommets = facettes;
      mesh.facettes_.resize(facettes.dimension(0), dim);
      ArrOfInt& liste_sommets_locaux = mesh.facettes_;
      const int n = liste_tous_sommets.size_array();
      for (int i = 0; i < n; i++)
        {
          const int som = liste_tous_sommets[i];
          const int i_som_loc = array_bsearch(indice_global_sommet, som);
          liste_sommets_locaux[i] = i_som_loc;
          if (i_som_loc < 0)
            liste_sommets_globaux.append_array(som);
        }
      array_trier_retirer_doublons(liste_sommets_globaux);
    }
    // Numero du processeur qui a genere la liste en cours de traitement
    int pe_sender = moi;
    int next_pe = (pe_sender + 1) % nproc;
    int prec_pe = (pe_sender + nproc - 1) % nproc;
    // Liste des processeurs a qui on envoie : un seul, le suivant
    // Processeurs de qui on recoit : le precedent
    const int send_list_size = (pe_sender == next_pe ? 0 : 1);
    ArrOfInt send_list(send_list_size);
    send_list.set_smart_resize(1);
    ArrOfInt recv_list(send_list_size);
    if (send_list_size)
      {
        send_list[0] = next_pe;
        recv_list[0] = prec_pe;
      }
    Schema_Comm schema_comm;
    schema_comm.set_send_recv_pe_list(send_list, recv_list);
    //send_list.set_smart_resize(1);
    send_list.resize_array(0);
    recv_list.resize_array(0);
    ArrOfInt new_list;
    new_list.set_smart_resize(1);
    for (int num_proc = 0; num_proc < nproc - 1; num_proc++)
      {
        // A chaque iteration, on transmet la liste au processeur suivant
        // Le processeur qui a envoye la liste a l'origine est donc le
        // precedent.
        pe_sender = (pe_sender + nproc - 1) % nproc;
        // Envoi de la liste de sommets au processeur suivant
        schema_comm.begin_comm();
        schema_comm.send_buffer(next_pe) << liste_sommets_globaux;
        schema_comm.echange_taille_et_messages();
        schema_comm.recv_buffer(prec_pe) >> liste_sommets_globaux;
        schema_comm.end_comm();
        // Pour chaque sommet de la liste, si je le possede,
        // je l'envoie a pe_sender et je le retire de la liste
        const int n = liste_sommets_globaux.size_array();
        new_list.resize_array(0);
        int flag = 0; // Va-on envoyer un sommet a pe_sender ?
        for (int i = 0; i < n; i++)
          {
            const int som = liste_sommets_globaux[i];
            const int i_som_loc = array_bsearch(indice_global_sommet, som);
            if (i_som_loc >= 0)
              {
                liste_sommets_to_send.append_array(i_som_loc);
                liste_pe_to_send.append_array(pe_sender);
                flag = 1;
              }
            else
              {
                new_list.append_array(som);
              }
          }
        liste_sommets_globaux = new_list;
        if (flag)
          send_list.append_array(pe_sender);
      }
    // A la fin, on doit avoir trouve tous les sommets
    if (liste_sommets_globaux.size_array() > 0)
      {
        Cerr << "Erreur sur PE " << moi << " : les sommets d'indices suivants"
             << " ne sont sur aucun processeur\n" << liste_sommets_globaux << finl;
        Process::exit();
      }
    // Creer les sommets virtuels
    reverse_send_recv_pe_list(send_list, recv_list);
    // TEMPORAIRE:
    Schema_Comm_FT sch;
    sch.set_send_recv_pe_list(send_list, recv_list);
    mesh.creer_sommets_virtuels(liste_sommets_to_send, liste_pe_to_send, sch);
  }

  // Echanger les indices globaux
  const int nb_sommets_reels = indice_global_sommet.size_array();
  indice_global_sommet.resize_array(mesh.nb_sommets());
  mesh.desc_sommets().echange_espace_virtuel(indice_global_sommet);

  // Traduire les indice globaux des sommets des facettes en indices locaux
  {
    const int nb_sommets_tot = indice_global_sommet.size_array();
    // Tableau facettes, vu comme un tableau monodimensionnel:
    const ArrOfInt& indices_globaux = facettes;
    ArrOfInt& indices_locaux = mesh.facettes_;
    const int n = indices_globaux.size_array();
    for (int i = 0; i < n; i++)
      {
        if (indices_locaux[i] < 0)
          {
            const int indice_global = indices_globaux[i];
            // Cherche l'indice local du sommet parmi les sommets virtuels
            int indice_local;
            for (indice_local = nb_sommets_reels; indice_local < nb_sommets_tot; indice_local++)
              if (indice_global_sommet[indice_local] == indice_global)
                break;
            assert(indice_local < nb_sommets_tot);
            indices_locaux[i] = indice_local;
          }
      }
  }
  mesh.desc_facettes_.calcul_schema_comm(mesh.nb_facettes());
  {
    const int n = mesh.facettes_.dimension(0);
    mesh.facette_num_owner_.resize_array(n);
    for (int i = 0; i < n; i++)
      mesh.facette_num_owner_[i] = i;
  }
}
