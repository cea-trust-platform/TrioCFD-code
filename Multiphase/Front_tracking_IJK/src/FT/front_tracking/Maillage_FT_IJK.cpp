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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Maillage_FT_IJK.cpp
// Directory : $IJK_ROOT/src/FT/front_tracking
//
/////////////////////////////////////////////////////////////////////////////

#include <Maillage_FT_IJK.h>
#include <IJK_Field.h>
#include <IJK_Lata_writer.h>
#include <LataDB.h>
#include <Process.h>
#include <communications.h>
#include <TRUSTTabFT.h>
#include <Comm_Group.h>
#include <Statistiques.h>
#include <TRUSTTabs.h>
#include <TRUSTArrays.h>

Implemente_instanciable(Maillage_FT_IJK,"Maillage_FT_IJK",Objet_U) ;

Sortie& Maillage_FT_IJK::printOn(Sortie& os) const
{
  Objet_U::printOn(os);
  return os;
}

Entree& Maillage_FT_IJK::readOn(Entree& is)
{
  Objet_U::readOn(is);
  return is;
}

void Maillage_FT_IJK::initialize(const IJK_Splitting& s, const Zone_dis& zone_dis, const Parcours_interface& parcours)
{
  ref_splitting_ = s;
  // Mise a jour des tableaux de processeurs voisins :
  initialize_processor_neighbourhood();
  nbmailles_euler_i_ = s.get_nb_elem_local(DIRECTION_I);
  nbmailles_euler_j_ = s.get_nb_elem_local(DIRECTION_J);
  nbmailles_euler_k_ = s.get_nb_elem_local(DIRECTION_K);

  associer_domaine_dis_parcours(zone_dis, parcours);
}

// Deplace tous les sommets lagrangiens dont les indices locaux sont mentionnes dans liste_sommets_initiale
//  d'un vecteur deplacement_initial.
// Methode parallele qui doit etre appelee simultanement par tous les processeurs.
// Elle gere les sommets qui changent de processeurs et met a jour le maillage en consequence.
// Numero_face_sortie servait en vdf et vef a savoir par quelle face on sortait d'un processeur
//  il fallait ensuite relancer un deplacement pour continuer jusqu'a atteindre le point d'arrivee
void Maillage_FT_IJK::deplacer_sommets(const ArrOfInt& liste_sommets,
                                       const DoubleTab& deplacement,
                                       ArrOfInt& liste_sommets_sortis,
                                       ArrOfInt& numero_face_sortie, int skip_facettes)
{
  if (Comm_Group::check_enabled()) check_mesh();
  const IJK_Splitting& splitting = ref_splitting_.valeur();
  assert(deplacement.dimension(0) == liste_sommets.size_array());
  assert(deplacement.dimension(1) == 3);

  // Donnees de decoupage du maillage dans les trois directions, pour trouver
  // sur quel processeur se trouve chaque point:
  VECT(ArrOfInt) slice_offsets(3);
  Int3 my_slice, my_new_slice;
  //int i;
  for (int i = 0; i < 3; i++)
    {
      splitting.get_slice_offsets(i, slice_offsets[i]);
      // On ajoute un dernier element au tableau, une tranche virtuelle tout a la fin
      // pour que la taille d'une tranche soit calculable par offset[n+1] - offset[n]
      const int n = slice_offsets[i].size_array();
      slice_offsets[i].resize_array(n+1);
      slice_offsets[i][n] = splitting.get_nb_items_global(IJK_Splitting::ELEM, i);
      my_slice[i] = splitting.get_local_slice_index(i);
    }
  // Pour chaque sommet, changer sa coordonnee et trouver l'indice de maille ou il se trouve,
  // puis faire les echanges de donnees paralleles pour attribuer le sommet au bon processeur s'il change
  // de processeur:
  const int my_processor = Process::me();
  const int nb_som = liste_sommets.size_array();
  ArrOfInt liste_processeurs_destination; // Pour chaque sommet envoye, a qui on l'envoie
  ArrOfInt send_pe_list; // Liste des processeurs a qui on envoie des sommets
  ArrOfInt drapeau_processeur_destination(Process::nproc());
  int max_voisinage = 0;
  send_pe_list.set_smart_resize(1);
  liste_sommets_sortis.set_smart_resize(1);
  liste_processeurs_destination.set_smart_resize(1);// algorithme d'allocation "predictif" pour append_array()
  for (int i_sommet = 0; i_sommet < nb_som; i_sommet++)
    {
      // Indice dans le maillage du sommet a deplacer:
      const int num_sommet = liste_sommets[i_sommet];
      if (sommet_virtuel(num_sommet))
        {
          Cerr << "Erreur dans  Maillage_FT_IJK::deplacer_sommets: \n"
               << " le sommet "<<num_sommet<< "est dans la liste a deplacer et est virtuel" << finl;
          Process::exit();
        }

      // Coordonnee courante du sommet (on la deplace d'element en element)
      Vecteur3 pos(sommets_, num_sommet);
      // Position du point d'arrivee
      Vecteur3 d(deplacement, i_sommet);
      pos += d;
      // Indice ijk de l'element contenant le sommet au depart
      Int3 ijk_pos = get_ijk_cell_index(num_sommet);
      // Trouver l'indice de la maille contenant le sommet:
      // si on passe une paroi, bloquer le sommet sur la paroi, si on passe une frontiere periodique, passer de l'autre cote:
      for (int direction = 0; direction < 3; direction++)
        {
          double coord = pos[direction];
          int i_maille = ijk_pos[direction] + splitting.get_offset_local(direction);
          const ArrOfDouble& coord_noeuds = splitting.get_grid_geometry().get_node_coordinates(direction);
          const int derniere_maille = coord_noeuds.size_array() - 2;
          // Avancer vers la droite si besoin
          while (coord_noeuds[i_maille+1] < coord)
            {
              // Le numero de la derniere maille est nbmailles-1 ou encore nbsommets-2
              if (i_maille == derniere_maille)
                {
                  //On est au bout du domaine... on colle le noeud sur la frontiere
                  coord = coord_noeuds[derniere_maille + 1];
                  break;
                }
              i_maille++;
            }
          // Avancer vers la gauche si besoin
          while (coord_noeuds[i_maille] > coord)
            {
              // Le numero de la derniere maille est nbmailles-1 ou encore nbsommets-2
              if (i_maille == 0)
                {
                  //On est au bout du domaine:
                  coord = coord_noeuds[0];
                  break;
                }
              i_maille--;
            }
          //      Cerr << "deplacement " << num_sommet << " " << direction << " de " << sommets_(num_sommet,direction) << " vers " << coord << finl;
          sommets_(num_sommet,direction) = coord;

          int current_slice	= my_slice[direction];
          while (i_maille >= slice_offsets[direction][current_slice+1])   // avancer vers la droite si besoin
            {
              current_slice++;
            }
          while (i_maille < slice_offsets[direction][current_slice])   // avancer vers la gauche si besoin
            {
              current_slice--;
            }
          my_new_slice[direction] = current_slice;
          ijk_pos[direction] = i_maille - slice_offsets[direction][current_slice];
        }

      // Pour le assert :
      Int3 nbmailles_euler_in_new_slice;
      for (int direction = 0; direction < 3; direction++)
        {
          int new_slice = my_new_slice[direction];
          nbmailles_euler_in_new_slice[direction] = slice_offsets[direction][new_slice+1] - slice_offsets[direction][new_slice] ;
          // Si on n'a pas change de slice dans cette direction, on doit retrouver nbmailles_euler_ijk_
          assert(!( (my_slice[direction] == my_new_slice[direction]) &&
                    (nbmailles_euler_in_new_slice[direction] != splitting.get_nb_elem_local(direction))));
        }
#if 0
      assert((ijk_pos[0] < 0 && ijk_pos[1] < 0 && ijk_pos[2] < 0)
             || (ijk_pos[0] >= 0 && ijk_pos[0] < nbmailles_euler_i_
                 && ijk_pos[1] >= 0 && ijk_pos[1] < nbmailles_euler_j_
                 && ijk_pos[2] >= 0 && ijk_pos[2] < nbmailles_euler_k_));
#else
      // Dans le cas ou le sommet a change de tranche, et que celle-ci n'a pas le meme nombre de mailles (plus de mailles)
      // et que le sommet est dans la derniere maille, il ne faut pas se comparer a nbmailles_euler_ijk_
      // mais plutot a nbmailles_euler_in_new_slice[direction] :
      assert((ijk_pos[0] < 0 && ijk_pos[1] < 0 && ijk_pos[2] < 0)
             || (ijk_pos[0] >= 0 && ijk_pos[0] < nbmailles_euler_in_new_slice[0]
                 && ijk_pos[1] >= 0 && ijk_pos[1] < nbmailles_euler_in_new_slice[1]
                 && ijk_pos[2] >= 0 && ijk_pos[2] < nbmailles_euler_in_new_slice[2]));

#endif
      // Traitement des changements de processeurs :
      int new_processor = splitting.get_processor_by_ijk(my_new_slice[0],my_new_slice[1],my_new_slice[2]);
      if (new_processor != my_processor)
        {
          // Le sommet change de processeur
          // Il faut l'envoyer au voisin:
          liste_sommets_sortis.append_array(num_sommet);
          liste_processeurs_destination.append_array(new_processor);
          if (drapeau_processeur_destination[new_processor] == 0)
            {
              // On n'a pas encore trouve de sommet a envoyer a ce processeur
              drapeau_processeur_destination[new_processor] = 1;
              send_pe_list.append_array(new_processor);
              max_voisinage = std::max(max_voisinage, voisinage_processeur_[new_processor]);
            }
        }
      // On stocke le numero local de l'element eulerien contenant le sommet
      //  sur le processeur qui contient ce sommet en reel.
      // (la structure de donnee est donc pour l'instant invalide car ce n'est pas
      //  l'indice local sur mon processeur mais sur le processeur ou le sommet sera reel
      //  apres avoir cree les sommets virtuels et echange les proprietaires des sommets)
#if 0
      set_ijk_cell_index(num_sommet, ijk_pos);
#else
      // Si on change de proc et que le nouveau n'a pas le meme nombre d'elem dans ses slices
      if ((new_processor != my_processor)
          && ( nbmailles_euler_i_ != nbmailles_euler_in_new_slice[0]
               || nbmailles_euler_j_ != nbmailles_euler_in_new_slice[1]
               || nbmailles_euler_k_ != nbmailles_euler_in_new_slice[2]))
        {
          // On fait comme set_ijk_cell_index mais pour le proc de destination, cad avec le nouveau nbmailles_euler...
          const int i = ijk_pos[0];
          const int j = ijk_pos[1];
          const int k = ijk_pos[2];
          assert((i < 0 && j < 0 && k < 0)
                 || (i >= 0 && i < nbmailles_euler_in_new_slice[0]
                     && j >= 0 && j < nbmailles_euler_in_new_slice[1]
                     && k >= 0 && k < nbmailles_euler_in_new_slice[2] ));
          if (i < 0)
            sommet_elem_[num_sommet]= -1;
          else
            sommet_elem_[num_sommet]= (k * nbmailles_euler_in_new_slice[1] + j) * nbmailles_euler_in_new_slice[0] + i;
        }
      else
        {
          set_ijk_cell_index(num_sommet, ijk_pos);
        }
#endif
    }
  max_voisinage = ::mp_max(max_voisinage);
  ArrOfInt recv_pe_list; // Liste des processeurs qui recoivent des sommets
  if (max_voisinage == 0)
    {
      // Il n'y a aucune communication de sommets
    }
  else if (max_voisinage == 1)
    {
      // On a trouve des communications par les faces... ajout des 6 voisins par faces
      send_pe_list = liste_processeurs_voisins_faces_;
      recv_pe_list = liste_processeurs_voisins_faces_;
    }
  else if (max_voisinage == 2)
    {
      // On a trouve des communications par les faces et/ou les aretes
      send_pe_list = liste_processeurs_voisins_aretes_;
      recv_pe_list = liste_processeurs_voisins_aretes_;
    }
  else if (max_voisinage == 3)
    {
      // On a trouve des communications par les faces/aretes/coins
      send_pe_list = liste_processeurs_voisins_coins_;
      recv_pe_list = liste_processeurs_voisins_coins_;
    }
  else
    {
      // On a trouve des communications autres qu'avec les voisins directs
      // On calcule la liste des communications a faire:
      reverse_send_recv_pe_list(send_pe_list, recv_pe_list);
    }
  // Mise a jour de l'espace virtuel des coodonnees des sommets
  desc_sommets_.echange_espace_virtuel(sommets_);

  if (max_voisinage > 0)
    {
      Schema_Comm schema_comm;
      schema_comm.set_send_recv_pe_list(send_pe_list, recv_pe_list);
      // Cree les sommets sur les processeurs destination, et initialise les coordonnees
      // le drapeau, le numero du processeur proprietaire (pour l'instant moi), etc...
      creer_sommets_virtuels(liste_sommets_sortis, liste_processeurs_destination, schema_comm);

      // Change le proprietaire des sommets qui ont ete envoyes:
      // On s'inspire de la methode echanger_sommets_PE()
      //  * Remplissage du tableau sommet_PE_owner avec le nouveau proprietaire
      //    et -1 si le sommet ne change pas de main,
      //  * Appel a la methode echanger_elements qui transforme
      //    les espaces distants et virtuels en fonction du nouveau proprietaire.
      sommet_PE_owner_ = -1;
      const int nechange = liste_sommets_sortis.size_array();
      for (int i = 0; i < nechange; i++)
        {
          const int sommet = liste_sommets_sortis[i];
          const int pe_destination = liste_processeurs_destination[i];
          sommet_PE_owner_[sommet] = pe_destination;
        }
      // L'echange espace virtuel va envoyer a tout le monde les nouveaux proprietaires de chaque sommet:
      desc_sommets_.echange_espace_virtuel(sommet_PE_owner_);
      // Envoi du numero d'element contenant les sommets a tous les processeurs
      desc_sommets_.echange_espace_virtuel(sommet_elem_);
      // Mise a jour des espaces distants et virtuels (c'est la qu'on change
      // reelement de proprietaire, on deplace des numeros de sommets d'une liste d'elements virtuels
      // a une autre (autre processeur) ou a une liste d'elements distants (si je deviens proprietaire)
      // ou inversement..., mais il n'y a pas d'envoi de message proprement dit dans cette etape)
      desc_sommets_.echanger_elements(sommet_PE_owner_);
      // Recalcul de sommet_PE_owner a partir des descripteurs (on l'avait ecrase
      // avec -1 partout sauf pour les sommets qui changent de procs).
      desc_sommets_.remplir_element_pe(sommet_PE_owner_);
      // On remet a -1 le numero d'element si le sommet n'est pas a moi
      // Mise a jour de sommet_num_owner_: on fait comme si les sommets etaient tous a moi (l'indice sur le
      // proprietaire est mon indice a moi), puis on echange pour ecraser les valeurs dont je ne suis pas proprietaire
      const int nsommets = sommets_.dimension(0);
      const int moi = me();
      for (int i = 0; i < nsommets; i++)
        {
          if (sommet_PE_owner_[i] != moi)
            sommet_elem_[i] = -1;
          sommet_num_owner_[i] = i;
        }
      // Envoi du nouveau sommet_num_owner_ a tous les autres procs qui possedent le sommet
      desc_sommets_.echange_espace_virtuel(sommet_num_owner_);
    }
  if (Comm_Group::check_enabled()) check_mesh(1, 1, 0); // ne pas tester que le proprietaire des facettes est conforme
  // a la convention definie dans maillage_ft_disc, ce n'est pas cense etre fait par cette fonction
  // On teste le compo_connexe : il n'a pas ete mis a jour, les facettes non plus... Ils sont donc en coherence.
}

void Maillage_FT_IJK::lire_maillage_ft_dans_lata(const char *filename_with_path, int tstep,
                                                 const char *geometryname)
{
  Cerr << "Maillage_FT_IJK::lire_maillage_ft_dans_lata fichier " << filename_with_path << " tstep=" << tstep
       << " geom=" << geometryname << finl;
  const IJK_Splitting& splitting = ref_splitting_.valeur();
  reset();

  const int master = Process::je_suis_maitre();
  Nom path, dbname;
  split_path_filename(filename_with_path, path, dbname);
  LataDB db;
  FloatTab coord;

  Vecteur3 origine;
  for (int j = 0; j < 3; j++)
    origine[j] = splitting.get_grid_geometry().get_origin(j);

  int nbsom = 0;

  if (master)
    {
      db.read_master_file(path, filename_with_path);

      const LataDBField& db_som = db.get_field(tstep, geometryname, "SOMMETS", "");
      db.read_data(db_som, coord);
      const LataDBField& db_elem = db.get_field(tstep, geometryname, "ELEMENTS", "");

      // Lecture de la composante connexe:
      {
        if (!db.field_exists(tstep, geometryname, "COMPO_CONNEXE"))
          {
            Cerr << "Erreur a la reprise d'une interface ft dans le fichier lata " << filename_with_path
                 << ":\n Le champ COMPO_CONNEXE est absent." << finl;
            Process::exit();
          }
        const LataDBField& db_field = db.get_field(tstep, geometryname, "COMPO_CONNEXE", "ELEM");
        IntTab tmp;
        db.read_data(db_field, tmp);
        compo_connexe_facettes_ = tmp;
      }

      db.read_data(db_elem, facettes_); // on lit directement les indices de sommets dans le tableau destination
      // Il faut traduire les indices de sommets virtuels en indices reels avec le tableau INDEX_SOMMET_REEL
      // En meme temps on supprime les facettes virtuelles
      // Il faut a ce moment faire suivre le tableau compo_connex pour qu'il soit coherent...
      // Si le champ INDEX_SOMMET_REEL n'est pas postraite, on suppose qu'il n'y a ni sommets virtuels
      // ni facettes virtuelles
      if (db.field_exists(tstep, geometryname, "INDEX_SOMMET_REEL"))
        {
          const LataDBField& db_index_sommet_reel = db.get_field(tstep, geometryname, "INDEX_SOMMET_REEL","");
          IntTab index_sommet_reel_tab;
          db.read_data(db_index_sommet_reel, index_sommet_reel_tab);
          // le tableau lu est dimensionne (nsom,1), on cast en arrofint pour ne pas avoir a donner 2 indices
          const ArrOfInt& index_sommet_reel = index_sommet_reel_tab;
          // Traduction. On cast le IntTab en ArrOfInt pour adresser lineairement tous les sommets
          // de tous les triangles. On a facettes_(i,j) = all_sommets[i*3+j]
          const int nb_facettes_tot = facettes_.dimension(0);
          // On va supprimer des facettes: index ou on ecrit la facette i:
          int idest = 0;
          for (int i = 0; i < nb_facettes_tot; i++)
            {
              // La facette est-elle reelle ? Oui si l'indice du premier sommet est un indice de sommet reel:
              const int i_som0 = facettes_(i,0);
              const int i_som0_renum = index_sommet_reel[i_som0];
              if (i_som0 == i_som0_renum)
                {
                  // Le sommet 0 de la facette est reel, on garde la facette
                  for (int j = 0; j < 3; j++)
                    {
                      const int i_sommet = facettes_(i,j);
                      const int i_som_renum = index_sommet_reel[i_sommet];
                      facettes_(idest,j) = i_som_renum;
                    }
                  // On met le tableau de compo connex en coherence :
                  compo_connexe_facettes_[idest] = compo_connexe_facettes_[i];
                  idest++;
                }
            }
          Journal() << " Renumerotation des facettes ft relues: " << nb_facettes_tot
                    << " facettes lues, " << idest << " facettes reelles conservees" << finl;
          facettes_.resize(idest, 3);
          // On resize aussi le tableau de compo connexe au cas ou il soit plus petit :
          compo_connexe_facettes_.resize_array(idest);
        }
      else
        {
          Journal() << " Pas de champ INDEX_SOMMET_REEL, on suppose qu'il n'y a pas de sommet virtuel" << finl;
        }
      nbsom = coord.dimension(0);
      Cerr << "Le maillage a nbsom=" << nbsom << " (avant suppression des sommets virtuels lus)" << finl;
      // Creation d'un maillage sur le proc0, tous les sommets a l'origine du maillage
      sommets_.resize(nbsom, 3);
      sommet_elem_.resize_array(nbsom);
      sommet_face_bord_.resize_array(nbsom);
      sommet_PE_owner_.resize_array(nbsom);
      sommet_num_owner_.resize_array(nbsom);
      drapeaux_sommets_.resize_array(nbsom);
      Int3 ijk; // indice de l'element d'origine
      ijk[0] = ijk[1] = ijk[2] = 0;
      int i;
      for (i = 0; i < nbsom; i++)
        {
          for (int j = 0; j < 3; j++)
            sommets_(i,j) = origine[j];
          set_ijk_cell_index(i, ijk);
          sommet_face_bord_[i] = -1; // le sommet n'est pas sur un bord
          sommet_PE_owner_[i] = 0; // Les sommets sont initialement sur le processeur maitre
          sommet_num_owner_[i] = i;
        }
      const int nbelem= facettes_.dimension(0);
      Cerr << "Le maillage a nbelem=" << nbelem << finl;
      facette_num_owner_.resize_array(nbelem);
      for (i = 0; i < nbelem; i++)
        facette_num_owner_[i] = i;

      // Avant, on lisait la compo_connexe ici,  apres avoir potentiellement reduit le nombre de facettes...
      // c'etait un bug. Il faut bien le faire avant comme c'est le cas maintenant
    }
  else
    {
      // On n'est pas sur le proc maitre. Pour l'instant, tous les elements lus sont sur le
      // proc maitre. On force donc la dimension des tableaux a zero :
      const int nbelem= 0;
      sommets_.resize(nbsom, 3);
      sommet_elem_.resize_array(nbsom);
      sommet_face_bord_.resize_array(nbsom);
      sommet_PE_owner_.resize_array(nbsom);
      sommet_num_owner_.resize_array(nbsom);
      drapeaux_sommets_.resize_array(nbsom);
      facette_num_owner_.resize_array(nbelem);
      compo_connexe_facettes_.resize_array(nbelem);
    }
  desc_sommets_.reset();
  desc_facettes_.reset();
  statut_ = MINIMAL;

  const int nbproc = Process::nproc();
  DoubleTab liste_sommets_recus;
  DoubleTab liste_coord_recues;
  IntTab  liste_facettes_recues;
  ArrOfInt  liste_compo_connexe_facettes_recues;
  if (master)
    {
      int nb_bubbles =0;
      if (compo_connexe_facettes_.size_array())
        // On ne peut pas faire un max_array sur un tableau vide
        nb_bubbles= max_array(compo_connexe_facettes_)+1;
      Cerr << "We read " << nb_bubbles  << " (real) bubbles on master proc and will redistribute them on procs." << finl;
      VECT(DoubleTab) liste_sommets(nbproc);
      VECT(DoubleTab) liste_coord(nbproc);
      VECT(IntTab) liste_facettes(nbproc);
      VECT(ArrOfInt) liste_compo_connexe_facettes(nbproc);
      for(int p=0; p<nbproc; p++)
        {
          liste_sommets[p].set_smart_resize(1);
          liste_sommets[p].resize(0,3);
          liste_coord[p].set_smart_resize(1);
          liste_coord[p].resize(0,3);
          liste_facettes[p].set_smart_resize(1);
          liste_facettes[p].resize(0,3);
          liste_compo_connexe_facettes[p].set_smart_resize(1);
          liste_compo_connexe_facettes[p].resize_array(0);
        }

      const int nbelem= facettes_.dimension(0);
      ArrOfInt dest_proc_sommets(nbsom);
      ArrOfInt sommets_renum(nbsom);
      for (int i = 0; i < nbelem; i++)
        {
          // On envoie les bulles sur d'autres proc par un modulo :
          const int icompo = compo_connexe_facettes_[i];
          int new_owner = icompo % nbproc;
          for (int j = 0; j < 3; j++)
            {
              const int i_som = facettes_(i,j);
              dest_proc_sommets[i_som] = new_owner;
            }
          liste_facettes[new_owner].append_line(facettes_(i,0),facettes_(i,1),facettes_(i,2));
          liste_compo_connexe_facettes[new_owner].append_array(icompo);
        }

      for (int i = 0; i < nbsom; i++)
        {
          const int new_owner = dest_proc_sommets[i];
          sommets_renum[i] = liste_sommets[new_owner].dimension(0); // il est mis a la fin de la liste du proc p
          liste_sommets[new_owner].append_line(sommets_(i,0),sommets_(i,1),sommets_(i,2));
          liste_coord[new_owner].append_line(coord(i,0),coord(i,1),coord(i,2));
        }
      // Correction des numeros des sommets formant les nouvelles facettes...
      for (int p=0; p< nbproc; p++)
        {
          const int nb_facettes_loc = liste_facettes[p].dimension(0);
          for (int i = 0; i < nb_facettes_loc; i++)
            for (int j=0; j<3; j++)
              {
                const int old_isom = liste_facettes[p](i,j); // le vieux numero du sommet
                const int new_isom = sommets_renum[old_isom]; // le nouveau numero auquel il correspond
                liste_facettes[p](i,j) = new_isom;// la correction de la facette.
              }
        }

      // Cheking for lost elements :
      {
        int sum_nbsom = 0;
        int sum_nbelem = 0;
        for(int p=0; p<nbproc; p++)
          {
            sum_nbsom += liste_sommets[p].dimension(0);
            sum_nbelem +=liste_facettes[p].dimension(0);
          }
        if (sum_nbsom != nbsom)
          {
            Cerr << "Some Front sommets have been lost!"  << finl;
            Process::exit();
          }
        if (sum_nbelem != nbelem)
          {
            Cerr << "Some Front facettes have been lost!"  << finl;
            Process::exit();
          }
      }
      // On a rempli les listes a envoyer, y compris pour le proc maitre;
      // Sauf que pour le mettre, on ne les envoie pas...
      Cerr << "Master starts sending packages..." << finl;
      // The loops starts at 1 because we dont want the master sending packages to himself!
      for(int p=1; p<nbproc; p++)
        {
          // On envoi chaque tableau sur un canal different :
          envoyer(liste_sommets[p],0,p,2000+4*p);
          envoyer(liste_facettes[p],0,p,2000+4*p+1);
          envoyer(liste_compo_connexe_facettes[p],0,p,2000+4*p+2);
          envoyer(liste_coord[p],0,p,2000+4*p+3);
        }
      Cerr << "Master is done sending packages..." << finl;

// Master fills lists for itself :
      liste_sommets_recus = liste_sommets[0];
      liste_facettes_recues = liste_facettes[0];
      liste_compo_connexe_facettes_recues = liste_compo_connexe_facettes[0];
      liste_coord_recues = liste_coord[0];
    }
  else
    {
      Journal() << "Receiving from master..." << finl;
      // On recoit les envois du maitre :
      const int p = Process::me();
      recevoir(liste_sommets_recus,0,p,2000+4*p);
      recevoir(liste_facettes_recues,0,p,2000+4*p+1 /* canal */);
      recevoir(liste_compo_connexe_facettes_recues,0,p,2000+4*p+2);
      recevoir(liste_coord_recues,0,p,2000+4*p+3);
      Journal() << "All data recieved..." << finl;
    }

  // All processes rebuild the listes sommets and facettes (master has not recieved it, but it filled it itself.
  const int p = Process::me();
  sommets_ = liste_sommets_recus;
  nbsom = sommets_.dimension(0);
// On ne peut pas faire (mismatched types) :  facettes_ = liste_facettes_recues;
  const int nbelem = liste_facettes_recues.dimension(0);
  assert(liste_compo_connexe_facettes_recues.size_array() == nbelem);
  facettes_.resize(nbelem,3);
  facette_num_owner_.resize_array(nbelem);
  compo_connexe_facettes_.resize_array(nbelem);
  for (int i=0; i< nbelem; i++)
    {
      for (int j = 0; j < 3; j++)
        facettes_(i,j) = liste_facettes_recues(i,j);
      facette_num_owner_[i] = i;
      compo_connexe_facettes_[i] = liste_compo_connexe_facettes_recues[i];
    }
// On remplit d'autres structures :
  // Les sommets et facettes sont tous reels et ils ont ete redistribues sur chaque processeurs.
  sommet_elem_.resize_array(nbsom);
  sommet_face_bord_.resize_array(nbsom);
  sommet_PE_owner_.resize_array(nbsom);
  sommet_num_owner_.resize_array(nbsom);
  drapeaux_sommets_.resize_array(nbsom); // Je ne sais pas a quoi il sert, donc je ne le remplis pas... Mais en tout cas, il faut qu'il ait la bonne taille!
  for (int i = 0; i < nbsom; i++)
    {
      sommet_elem_[i] = 0; // On le met en 0,0,0 provisoirement. Cela changera au deplacement.
      sommet_face_bord_[i] = -1; // le sommet n'est pas sur un bord
      sommet_PE_owner_[i] = p; // Les sommets sont redistribues
      sommet_num_owner_[i] = i;
      drapeaux_sommets_[i] = 0; // Je ne sais pas quoi mettre...
    }

  desc_sommets_.calcul_schema_comm(nb_sommets());
  desc_facettes_.calcul_schema_comm(nb_facettes());

  // GB 2020/04/20. Je ne pense pas sur que cet appel soit necessaire si ce qui precede est bien fait.
  //                De plus, la creation possible de sommets virtuels est un probleme pour le deplacement ci-dessous,
  //                Dont la taille doit etre la meme que liste_coord_recues, c'est a dire nbsom (avant creation de potentiels som virt)
  // corriger_proprietaires_facettes();
  DoubleTab deplacement(nbsom,3);
  for (int i = 0; i < nbsom; i++)
    {
      for (int j = 0; j < 3; j++)
        deplacement(i,j) = liste_coord_recues(i,j) - origine[j];
    }
  transporter(deplacement);
// Suppression des sommets inutilises et des facette virtuelles:
  nettoyer_maillage();

}
void Maillage_FT_IJK::transporter(const DoubleTab& deplacement)
{
  // La mise a jour du tableau compo_connexe_facette est realisee par une
  // surcharge de creer_facettes_virtuelles, c'est une methode de plus bas niveau
  // qui est utilisee aussi lors du parcours de l'interface, on centralise
  // la gestion des compo connexes a moins d'endroits du code.
  // Inconvenient: la mise a jour des compo connexes sera faite plus souvent,
  // pas optimal en nombre d'operations.

  // preparer_tableau_avant_transport(compo_connexe_facettes_, desc_facettes_);
  Maillage_FT_Disc::transporter(deplacement);
  // update_tableau_apres_transport(compo_connexe_facettes_, nb_facettes() , desc_facettes_);
}

void Maillage_FT_IJK::nettoyer_maillage()
{
  // Nettoyage des composantes connexes:
  {
    const int nb_fa7 = facettes_.dimension(0);
    int n = 0;
    for (int i = 0; i < nb_fa7; i++)
      {
        const int invalide = (facettes_(i,0) == facettes_(i,1));
        const int virtuelle = facette_virtuelle(i);
        if (!invalide && !virtuelle)
          {
            compo_connexe_facettes_[n] = compo_connexe_facettes_[i];
            n++;
          }
      }
    compo_connexe_facettes_.resize_array(n);
  }
  // Appel a l'ancienne methode
  Maillage_FT_Disc::nettoyer_maillage();
}

void Maillage_FT_IJK::creer_facettes_virtuelles(const ArrOfInt& liste_facettes,
                                                const ArrOfInt& liste_facettes_pe,
                                                const ArrOfInt& facettes_send_pe_list,
                                                const ArrOfInt& facettes_recv_pe_list)
{
  // Appel a l'ancienne methode
  Maillage_FT_Disc::creer_facettes_virtuelles(liste_facettes,
                                              liste_facettes_pe,
                                              facettes_send_pe_list,
                                              facettes_recv_pe_list);

  // Mise a jour des composantes connexes
// static Stat_Counter_Id cmpt = statistiques().new_counter(0, "Mise a jour compo connex dans creer facette virt");
  //statistiques().begin_count(cmpt);
  compo_connexe_facettes_.resize_array(nb_facettes());
  desc_facettes().echange_espace_virtuel(compo_connexe_facettes_);

  // La verification du tableau des compo_connexe n'a pas ete faite par
  // Maillage_FT_Disc::creer_facettes_virtuelles. On la fait ici :
  if (Comm_Group::check_enabled())
    check_mesh(1, 1 /* ne pas tester le proprietaire de la facette */);
  //statistiques().end_count(cmpt);

}

void Maillage_FT_IJK::initialize_processor_neighbourhood()
{
  const IJK_Splitting& splitting = ref_splitting_.valeur();
  int np = Process::nproc();
  int npx = splitting.get_nprocessor_per_direction(0);
  int npy = splitting.get_nprocessor_per_direction(1);
  int npz = splitting.get_nprocessor_per_direction(2);
  voisinage_processeur_.resize_array(np);
  liste_processeurs_voisins_faces_.set_smart_resize(1);
  liste_processeurs_voisins_coins_.set_smart_resize(1);
  liste_processeurs_voisins_aretes_.set_smart_resize(1);

  voisinage_processeur_ = 4; // On initialize le tableau comme non-voisins.
  Int3 my_ijk, voisin_ijk;
  for (int i = 0; i < 3; i++)
    my_ijk[i] = splitting.get_local_slice_index(i);

  for (int i = -1; i< 2; i++)
    {
      for (int j = -1; j< 2; j++)
        {
          for (int k = -1; k< 2; k++)
            {
              voisin_ijk = my_ijk;
              voisin_ijk[0] += i;
              voisin_ijk[1] += j;
              voisin_ijk[2] += k;
              if (voisin_ijk[0] >=0 && voisin_ijk[0]<npx &&
                  voisin_ijk[1] >=0 && voisin_ijk[1]<npy &&
                  voisin_ijk[2] >=0 && voisin_ijk[2]<npz )
                {
                  int rang_voisin = splitting.get_processor_by_ijk(voisin_ijk[0],voisin_ijk[1],voisin_ijk[2]);
                  int max_voisinage = abs(i)+abs(j)+abs(k);
                  voisinage_processeur_[rang_voisin] = max_voisinage;
                  if (max_voisinage == 0)
                    {
                      // Je suis sur me.
                    }
                  else if (max_voisinage == 1)
                    {
                      liste_processeurs_voisins_coins_.append_array(rang_voisin);
                      liste_processeurs_voisins_aretes_.append_array(rang_voisin);// On a un voisin par les aretes.
                      liste_processeurs_voisins_faces_.append_array(rang_voisin); // On a un voisin par les faces.
                    }
                  else if (max_voisinage == 2)
                    {
                      liste_processeurs_voisins_coins_.append_array(rang_voisin);
                      liste_processeurs_voisins_aretes_.append_array(rang_voisin);// On a un voisin par les aretes.
                    }
                  else if (max_voisinage == 3)
                    {
                      liste_processeurs_voisins_coins_.append_array(rang_voisin);
                    }
                  else
                    {
                      // Je ne suis pas voisin.
                    }
                }
            }
        }
    }
}

// Surcharge
int Maillage_FT_IJK::check_mesh(int error_is_fatal, int skip_facette_pe, int skip_facettes) const
{
  if (!Comm_Group::check_enabled())
    //Optimisation : Ne pas faire le check_mesh sauf si on est en debug...
    // L'applique partout dans le code...
    return -1;
  Maillage_FT_Disc::check_mesh(error_is_fatal, skip_facette_pe, skip_facettes);
  const int nf = nb_facettes();
  if (nf != compo_connexe_facettes_.size_array())
    {
      Journal() << "Erreur Maillage_FT_IJK::check_mesh : taille de compo_connexe_facettes_ invalide" << finl;
      Journal() << "nb_facettes = " << nf << " et compo_connexe_facettes_ = " <<  compo_connexe_facettes_.size_array() << finl ;
      assert(0);
      exit();
    }
  return -1;
}

// Surcharge de Maillage_FT_Disc::check_sommets, specialise pour le IJK:
int Maillage_FT_IJK::check_sommets(int error_is_fatal) const
{
  const double invalid_value = DMAXFLOAT*0.9;

  if (statut_ == RESET)
    {
      int ok = (nb_sommets() == 0);
      ok = ok && (sommet_elem_.size_array() == 0);
      ok = ok && (sommet_face_bord_.size_array() == 0);
      ok = ok && (sommet_PE_owner_.size_array() == 0);
      ok = ok && (sommet_num_owner_.size_array() == 0);
      ok = ok && (desc_sommets_.espace_virtuel().pe_voisins().size_array() == 0);
      ok = ok && (desc_sommets_.espace_distant().pe_voisins().size_array() == 0);
      if (!ok && error_is_fatal)
        {
          Journal() << "Erreur Maillage_FT_Disc::check_sommets : maillage RESET invalide" << finl;
          assert(0);
          exit();
        }
      return !ok;
    }
  // Verification que les espaces distants et virtuels sont coherents:
  desc_sommets_.check();
  const int nsom = sommets_.dimension(0);
  // Verification de sommet_PE_owner_ : on le recalcule et on compare
  {
    if (sommet_PE_owner_.size_array() != nsom)
      {
        if (error_is_fatal)
          {
            assert(0);
            exit();
          }
        Cerr << "" << finl;
      }
    ArrOfIntFT pe_owner(nsom);
    desc_sommets_.remplir_element_pe(pe_owner);
    int som ;
    for (som = 0; som < nsom; som++)
      {
        if (sommet_PE_owner_[som] != pe_owner[som])
          {
            if (error_is_fatal)
              {
                assert(0);
                exit();
              }
            break;
          }
      }
    if (som < nsom)
      Cerr << "Erreur Verification de sommet_PE_owner_" << finl;
  }
  // Verification de sommet_num_owner_ : on le recalcule et on compare
  {
    if (sommet_num_owner_.size_array() != nsom)
      {
        Journal() << "Erreur sommet_num_owner_.size_array() != nb_sommets" << finl;
        if (error_is_fatal)
          {
            assert(0);
            exit();
          }
      }
    ArrOfIntFT num_owner(nsom);
    for (int i = 0; i < nsom; i++)
      num_owner[i] = i;
    desc_sommets_.echange_espace_virtuel(num_owner);
    for (int i = 0; i < nsom; i++)
      {
        if (num_owner[i] != sommet_num_owner_[i])
          {
            Journal() << "Erreur num_owner[" << i << "] = " << num_owner[i];
            Journal() << "  sommet_num_owner_ = " << sommet_num_owner_[i] << finl;
            if (error_is_fatal)
              {
                assert(0);
                exit();
              }
          }
      }
  }
  // Verification des coordonnees des sommets virtuels : comparaison avec une copie echangee
  {
    DoubleTabFT copie_sommets = sommets_;
    // On invalide les elements virtuels. Cela permet de detecter le cas
    // ou un element se trouverait a la fois dans l'espace distant et dans
    // l'espace virtuel.
    const Descripteur_FT& espace_virtuel = desc_sommets_.espace_virtuel();
    const ArrOfInt& pe_voisins = espace_virtuel.pe_voisins();
    const int nb_pe_voisins = pe_voisins.size_array();
    for (int indice_pe = 0; indice_pe < nb_pe_voisins; indice_pe++)
      {
        const int pe = pe_voisins[indice_pe];
        const ArrOfInt& elements = espace_virtuel.elements(pe);
        const int n = elements.size_array();
        for (int i = 0; i < n; i++)
          {
            const int num_sommet = elements[i];
            copie_sommets(num_sommet, 0) = invalid_value;
            copie_sommets(num_sommet, 1) = invalid_value;
            copie_sommets(num_sommet, 2) = invalid_value;
          }
      }
    // Echange espace virtuel sur la copie
    desc_sommets_.echange_espace_virtuel(copie_sommets);
    // Comparaison de la copie et du tableau des sommets.
    for (int i = 0; i < nsom; i++)
      {
        for (int j = 0; j < Objet_U::dimension; j++)
          {
            if (copie_sommets(i,j) == invalid_value)
              {
                Journal() << "Erreur copie_sommets(" << i << ",";
                Journal() << j << ") == DMAX_FLOAT" << finl;
                if (error_is_fatal)
                  {
                    assert(0);
                    exit();
                  }
              }
            if (copie_sommets(i,j) != sommets_(i,j))
              {
                Journal() << "Erreur copie_sommets(" << i << ",";
                Journal() << j << ") = " << copie_sommets(i,j);
                Journal() << " sommets_ = " << sommets_(i,j) << finl;
                if (error_is_fatal)
                  {
                    assert(0);
                    exit();
                  }
              }
          }
      }
  }

  // Verification de sommet_elem_ :
  {
    if (sommet_elem_.size_array() != nsom)
      {
        Journal() << "Erreur sommet_elem_.size_array() != nb_sommets" << finl;
        if (error_is_fatal)
          {
            assert(0);
            exit();
          }
      }
    // On verifie que les sommets virtuels ont bien sommet_elem_ < 0:
    const Descripteur_FT& espace_virtuel = desc_sommets_.espace_virtuel();
    const ArrOfInt& pe_voisins = espace_virtuel.pe_voisins();
    const int nb_pe_voisins = pe_voisins.size_array();
    for (int indice_pe = 0; indice_pe < nb_pe_voisins; indice_pe++)
      {
        const int pe = pe_voisins[indice_pe];
        const ArrOfInt& elements = espace_virtuel.elements(pe);
        const int n = elements.size_array();
        for (int i = 0; i < n; i++)
          {
            const int num_sommet = elements[i];
            if (sommet_elem_[num_sommet] >= 0)
              {
                Journal() << "Erreur sommet_elem_[" << num_sommet << "] = " << sommet_elem_[num_sommet]
                          << " Il devrait etre negatif car le sommet est virtuel" << finl;
                if (error_is_fatal)
                  {
                    assert(0);
                    exit();
                  }
              }
          }
      }

  }
  return 0;

}

void Maillage_FT_IJK::calculer_compo_connexe_sommets(ArrOfIntFT& compo_connexe_sommets) const
{

  const int nb_fa7 = facettes_.dimension(0);
  const int nbsom = sommets().dimension(0);
  compo_connexe_sommets.resize_array(nbsom, Array_base::NOCOPY_NOINIT); // tous les sommets, y compris virtuels.
  compo_connexe_sommets = -10000000; // Force une initialisation bidon.
  // On parcours toutes les facettes.
  for (int i = 0; i < nb_fa7; i++)
    {
      const int icompo = compo_connexe_facettes_[i];
      for (int j = 0; j < 3; j++)
        {
          const int i_som = facettes_(i,j);
          compo_connexe_sommets[i_som] = icompo;
        }
    }


  // Echange entre les processeurs pour mettre a jour la valeur aux sommets isoles sur un proc
  // qui n'appartiendrait qu'a des facettes d'autres processeurs.
  // On prend le max sur tous les processeurs qui partagent le sommet pour les sommets isoles
  const Desc_Structure_FT& desc_som = desc_sommets();
  desc_som.collecter_espace_virtuel(compo_connexe_sommets, MD_Vector_tools::EV_MAX);
  desc_som.echange_espace_virtuel(compo_connexe_sommets);

}

// Fait appel a la methode Maillage_FT_Disc::recopie pour copier le maillage.
// Puis initialise le tableau de composantes connexes avec la valeur imposee icompo.
void Maillage_FT_IJK::recopie_force_compo(const Maillage_FT_IJK& source_mesh, const int icompo)
{
  Cerr << "Methode Maillage_FT_IJK::recopie_force_compo inusite jusqu'a present"
       << " donc a tester! " << finl;
  Process::exit();
  recopie(source_mesh, MINIMAL);
  const int nf = nb_facettes();
  compo_connexe_facettes_.resize_array(nf, Array_base::NOCOPY_NOINIT); // tous les sommets, y compris virtuels.
  for (int i_facette = 0; i_facette < nf; i_facette++)
    {
      compo_connexe_facettes_[i_facette] = icompo;
    }
}

// Surcharge la methode Maillage_FT_Disc::recopie pour copier le maillage.
// Puis initialise le tableau de composantes connexes avec la valeur dans source_mesh.
void Maillage_FT_IJK::recopie(const Maillage_FT_Disc& source_mesh, Statut_Maillage niveau_copie)
{
  if (sub_type(Maillage_FT_IJK, source_mesh))
    {
      const Maillage_FT_IJK& source_mesh_ijk = ref_cast(Maillage_FT_IJK, source_mesh);

      nbmailles_euler_i_ = source_mesh_ijk.nbmailles_euler_i_;
      nbmailles_euler_j_ = source_mesh_ijk.nbmailles_euler_j_;
      nbmailles_euler_k_ = source_mesh_ijk.nbmailles_euler_k_;
      voisinage_processeur_ = source_mesh_ijk.voisinage_processeur_;
      liste_processeurs_voisins_faces_=source_mesh_ijk.liste_processeurs_voisins_faces_;
      liste_processeurs_voisins_aretes_=source_mesh_ijk.liste_processeurs_voisins_aretes_;
      liste_processeurs_voisins_coins_=source_mesh_ijk.liste_processeurs_voisins_coins_;

      Maillage_FT_Disc::recopie(source_mesh, niveau_copie);
      compo_connexe_facettes_.copy_array(source_mesh_ijk.compo_connexe_facettes());
    }
  else
    {
      // soit erreur, soit initialise avec une valeur par defaut.
      Cerr << "Maillage_FT_IJK::recopie_maillage : cette surcharge s'attend a  "
           << "recevoir un Maillage_FT_IJK" << finl;
      Process::exit();
    }

}

// Surcharge de la methode Maillage_FT_Disc::ajouter_maillage
void Maillage_FT_IJK::ajouter_maillage(const Maillage_FT_Disc& maillage_tmp,int skip_facettes)
{

  if (sub_type(Maillage_FT_IJK, maillage_tmp))
    {
      ajouter_maillage_IJK(ref_cast(Maillage_FT_IJK, maillage_tmp)) ; // copie les compo connexes existantes
    }
  else
    {
      // soit erreur, soit initialise avec une valeur par defaut.
      Cerr << "Maillage_FT_IJK::ajouter_maillage : cette surcharge s'attend a  "
           << "recevoir un Maillage_FT_IJK" << finl;
      Process::exit();
    }
}


// Pour remplir correctement le tableau des composantes connexes
// a partir du tableau compo_connex du maillage source.
void Maillage_FT_IJK::ajouter_maillage_IJK(const Maillage_FT_IJK& added_mesh)
{

  const int nb_facettes_ini = nb_facettes(); // du maillage initial : *this.
  const int nb_facettes_add = added_mesh.nb_facettes(); // du maillage a ajouter.
  const ArrOfInt& compo_connexe_add = added_mesh.compo_connexe_facettes();  // Compo connexes a ajouter

  // On etends le tableau de composantes connexes en ajoutant celui de added_mesh a la fin:
  compo_connexe_facettes_.resize_array(nb_facettes_ini+nb_facettes_add);
  compo_connexe_facettes_.inject_array(compo_connexe_add, -1, /* tous les elements de la source */
                                       nb_facettes_ini, /* indice destination */
                                       0); /* indice source */

  // La methode Maillage_FT_Disc::ajouter_maillage va ajouter les facettes, les sommets...
  // et elle va mettre a jour les desc et les num_owner... en ajoutant des elements dans les decripteurs
  // Mais elle ne fait pas d'echange. Pas de renumerotation (seulement la gestion du decalage nb_facettes_ini
  // pour les elements ajoutes aux desc).
  // La compo connexe est donc valide apres :
  Maillage_FT_Disc::ajouter_maillage(added_mesh);

}


void Maillage_FT_IJK::calculer_costheta_minmax(DoubleTab& costheta) const
{
  const int nb_som = nb_sommets();
  costheta.resize(nb_som, 2);
  costheta = 0. ;
}

static double distance_sommets(const DoubleTab& sommets, const int s0, const int s1)
{
  const double d = (sommets(s0,0) - sommets(s1,0)) * (sommets(s0,0) - sommets(s1,0))
                   + (sommets(s0,1) - sommets(s1,1)) * (sommets(s0,1) - sommets(s1,1))
                   + (sommets(s0,2) - sommets(s1,2)) * (sommets(s0,2) - sommets(s1,2));
  return sqrt(d);
}

double Maillage_FT_IJK::minimum_longueur_arrete() const
{
  // nettoyer_maillage();
  const int nb_fa7 = facettes_.dimension(0);
  const DoubleTab& tab_sommets = sommets();
  double lg = 1.e10;
  double petit  = 1.e-7;

  for (int i = 0; i < nb_fa7; i++)
    {
      const int s0 = facettes_(i,0);
      const int s1 = facettes_(i,1);
      const int s2 = facettes_(i,2);
      double d1 = distance_sommets(tab_sommets, s0, s1);
      double d2 = distance_sommets(tab_sommets, s0, s2);
      double d3 = distance_sommets(tab_sommets, s1, s2);
      //Cerr << d1 << " " << d2 << " " << d3 << finl;

      if ((d1> petit) && (d1 < lg))
        {
          lg = d1;
        }
      if ((d2> petit) && (d2 < lg))
        {
          lg = d2;
        }
      if ((d3> petit) && (d3 < lg))
        {
          lg = d3;
        }
    }

  // Syncronisation :
  lg = Process::mp_min(lg);
  return lg;

}
