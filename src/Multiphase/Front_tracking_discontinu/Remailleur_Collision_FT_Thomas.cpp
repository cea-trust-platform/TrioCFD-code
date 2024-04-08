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
// File:        Remailleur_Collision_FT_Thomas.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Remailleur_Collision_FT_Thomas.h>
#include <Domaine_VF.h>
#include <Domaine.h>
#include <Param.h>
#include <communications.h>
#include <SFichier.h>

Implemente_instanciable_sans_constructeur(Remailleur_Collision_FT_Thomas,"Remailleur_Collision_FT_Thomas",Remailleur_Collision_FT_Juric);

Remailleur_Collision_FT_Thomas::Remailleur_Collision_FT_Thomas()
{
  tester_ = 0;
  est_dimensionne_ = 0;
  distance_utilisateur_ = 1;
}

/*! @brief Lecture des parametres dans le fichier .
 *
 * data
 *
 */
Entree& Remailleur_Collision_FT_Thomas::readOn(Entree& is)
{
  Cerr << "Remailleur_Collision_FT_Thomas::readOn()" << finl;
  Param param(que_suis_je());
  param.ajouter_flag("tester", & tester_);
  param.ajouter("distance_interface_element_max", & distance_utilisateur_);
  param.lire_avec_accolades_depuis(is);

  if (tester_)
    Cerr << " Test des fonctions de Remailleur_Collision_FT_Thomas." << finl;
  if (distance_utilisateur_ <= 0.)
    {
      Cerr << "Error: distance_interface_element_max must be positive" << finl;
      barrier();
      exit();
    }
  return is;
}

/*! @brief Erreur.
 *
 * Pas code.
 *
 */
Sortie& Remailleur_Collision_FT_Thomas::printOn(Sortie& os) const
{
  Cerr << "Erreur  Remailleur_Collision_FT_Thomas::printOn() n'est pas code." << finl;
  exit();
  return os;
}

//Initialise les attributs "voisinage_sommet_" et "next_elem_" qui sont des IntTab
//et qui forment a eux deux une liste chainee
//Principe : a un numero d'element eulerien "k" et a un numero local de sommet "sloc"
//correspondent deux numeros uniques ;
//-le numero global "sglob" associe a "sloc"
//-l'entier "p" definit par p=k*nb_som_elem+sloc.
//Si l'on realise la division euclidienne de "p" par "nb_som_elem", nous retrouvons
//facilement le numero globlal du sommet associe a "sloc" et l'element eulerien "k" qui le contient.
//Par consequent, dans "voisinage_sommet_[si]" on stocke le "p" qui correspond au premier
//element "k" contenant "si", dans "next_elem[p]" l'entier "p1" qui contient le second
//element "k1" contenant "si", dans "next_elem[p1]" l'entier "p2" qui contient le
//troisieme element "k2" contenant "si" ... jusqu'a la fin de liste repere par -1.
//Soit "s" un sommet du maillage eulerien.
//REMARQUE : il n'y a pas redondance d' elements dans les differentes
//listes "voisinage_sommet_[s]"
//REMARQUE : pour le parallele, on stocke egalement les elements voisins qui sont virtuels
int Remailleur_Collision_FT_Thomas::construire_voisinage_sommet(const Maillage_FT_Disc& maillage)
{
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,domaine_dis(maillage).valeur());
  const Domaine& domaine = domaine_VF.domaine();

  const int nb_elem_tot = domaine.nb_elem_tot();
  const int nb_som_elem = domaine.nb_som_elem();
  const int end_liste = -1;

  //Tableau des sommets par element
  const IntTab& elem_sommets = domaine.les_elems();

  assert(voisinage_sommet_.dimension_tot(0) == domaine.nb_som_tot());
  assert(next_elem_.dimension_tot(0) == domaine.nb_elem_tot()*domaine.nb_som_elem());

  voisinage_sommet_ = end_liste;

  //Construction du voisinage d'un sommet "s" :
  //On boucle sur un element "elem" dont on stocke le numero global
  //On boucle sur les sommets "s_i" de "elem" dont on recupere les numeros globaux
  //Si "voisinage_sommet_[s_i] = -1, on initialise "voisinage_sommet_[s_i]" a "elem*nb_som_elem+numero_local_sommet"
  //Si "voisinage_sommet_[s_i] = l != -1 :
  //on calcule "elem*nb_som_elem+numero_local_sommet"
  //on transfere la valeur "l"  dans "next_elem["elem*nb_som_elem+numero_local_sommet"]
  //on modifie "voisinage_sommet_[s_i]" que l'on fixe a "elem*nb_som_elem+numero_local_sommet"
  for (int elem=0; elem<nb_elem_tot; elem++)
    for (int som=0; som<nb_som_elem; som++)
      {
        const int som_global =
          domaine.get_renum_som_perio(elem_sommets(elem,som));

        const int p = nb_som_elem*elem+som;

        next_elem_[p] = voisinage_sommet_[som_global];
        voisinage_sommet_[som_global] = p;
      }//fin des for

  tmp_flag_elements_.resize_array(nb_elem_tot);
  tmp_flag_elements_ = 0;
  return 1;
}


//Mets a jour l'attribut "distance_interface_element_eulerien" qui est un tableau d'entiers
//Un element eulerien a la distance 0 de l'interface est un element traverse par l'interface
//et tel que la surface de cette portion d'interface soit non nulle
//Un element eulerien a la distance 1 de l'interface est un element voisin d'un element a distance 0
//qui n'est pas lui-meme a distance 0
//Et ainsi de suite ...
//Voisinage d'un sommet "s" := union des elements possedant ce sommet
//Voisinage d'un element "elem" := union des voisinages des sommets "s_i" de "elem"
//Mets a jour l'attribut "nombre_de_voisins_plus_proches_" qui est tableau d'entiers
//Est associe a un element eulerien, le nombre des ses elements voisins qui sont
//plus proches de l'interface qu'il ne l'est lui-meme
//Mets a jour l'attribut "surface_interface_elements_voisins_" qui est un tableau de double
//A un element eulerien a distance 0 de l'interface est associe la surface d'interface le traversant
//A un element eulerien a distance 1 de l'interface est associe la somme des surfaces d'interface
//traversant ses voisins a distance 0 de l'interface
int  Remailleur_Collision_FT_Thomas::mettre_a_jour_data(const Maillage_FT_Disc& maillage)
{
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,domaine_dis(maillage).valeur());
  const Domaine& domaine = domaine_VF.domaine();

  const Intersections_Elem_Facettes&  intersection_elem_facettes =
    maillage.intersections_elem_facettes();

  const ArrOfInt& index_facette = intersection_elem_facettes.index_facette();

  const int nb_elem_tot = domaine.nb_elem_tot();
  const int nb_elem = domaine.nb_elem();
  const int nb_facettes_interface = maillage.nb_facettes();
  const int seen = 1;

  int nb_elem_parcourus = 0;
  int numero_voisinage = 0;

  //Indispensable pour le parallele : il faut s'assurer
  //que les procs effectuent bien le meme nombre de boucle.
  //Ce nombre est cree pour assurer cette contrainte.
  int ready_to_quit = 0;

  //Indispensable pour le parallele mais LOCAL a chaque proc
  //Ne stocke des valeurs utiles que pour les elements VIRTUELS
  IntTab vu_localement(nb_elem_tot-nb_elem);

  //Par convention :
  //"elem_proches" pointe toujours sur la liste des elements les plus
  //proches de l'interface
  //"elem_eloignes" pointe toujours sur la liste des elements les plus
  //eloignes de l'interface
  ArrOfInt elements_proches(nb_facettes_interface),elements_eloignes(3*nb_facettes_interface);
  ArrOfInt *elem_proches, *elem_eloignes; //eloignes = distance_(x+1)

  elements_proches.resize_array(0);
  elements_eloignes.resize_array(0);

  elem_proches = &(elements_proches);
  elem_eloignes = &(elements_eloignes);


  //
  //Initialisation des attributs "distance_interface_element_eulerien_"
  //"nombre_de_voisins_plus_proches_" et "surface_interface_elements_voisins_"
  //
  distance_interface_element_eulerien_=-1;
  nombre_de_voisins_plus_proches_=0;
  surface_interface_elements_voisins_=0.;
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  //
  //Initialisation des elements a distance 0 de l'interface
  //Initialisation du nombre de voisins plus proches de ces memes elements
  //Initialisation de la surface totale d'interface traversant ces memes elements
  //
  for (int facette=0; facette<nb_facettes_interface; facette++)
    {
      const double surface_facette = surface_facettes[facette];
      //
      //Parcours des elements  REELS traverses par "facette"
      //et initialisation de la liste des elements  REELS a
      //distance 0 de l'interface, du nombre de leurs
      //voisins les plus proches, de la surface totale
      //d'interface les traversant
      //
      int index = index_facette[facette];

      while (index >= 0)
        {
          const Intersections_Elem_Facettes_Data& data =
            intersection_elem_facettes.data_intersection(index);

          const int elem = data.numero_element_;
          assert(elem<nb_elem);

          //Condition determinant les elements a distance 0 de l'interface
          if (distance_interface_element_eulerien_[elem]==-1 && data.fraction_surface_intersection_!=0.)
            {
              distance_interface_element_eulerien_[elem] = 0;
              elem_proches->append_array(elem);
            }

          //Calcul de la surface d'interface coupant les elements a distance
          //0 de cette meme interface
          surface_interface_elements_voisins_[elem] += data.fraction_surface_intersection_ * surface_facette;


          index = data.index_element_suivant_;
        }//fin du while

    }//fin du for sur facette


  //Pour le parallele : il faut s'assurer que les processus en cours
  //ont bien finis leurs boucles respectives
  distance_interface_element_eulerien_.echange_espace_virtuel();
  surface_interface_elements_voisins_.echange_espace_virtuel();

  //Nous comptons les elements deja parcourus :
  //nous ne pouvons pas utiliser seulement la fonction size_array()
  //des ArrOfInt appliquee a "elem_proches" parce que nous parcourons
  //des elements REELS par definition de la structure de donnee "data".
  //Or nous voulons parcourir tous les elements y compris les elements
  //VIRTUELS. D'ou l'astuce suivante qui doit etre realisee APRES
  //la fonction echange_espace_virtuel().
  //De meme, il faut rajouter dans la liste "elem_proches" les elements
  //VIRTUELS qui sont a distance 0 de l'interface
  nb_elem_parcourus+=elem_proches->size_array();

  for (int elem=nb_elem; elem<nb_elem_tot; elem++)
    if (surface_interface_elements_voisins_[elem]>0.)
      {
        elem_proches->append_array(elem);
        nb_elem_parcourus++;
      }

  //
  //Initialisation des autres elements a distance >= 1 de l'interface
  //Principe : supposant connu les element a distance n de l'interface,
  //nous recherchons les voisins de ces elements dont la distance n'a
  //pas encore ete calculee (distance = -1)
  //Nous initialisons alors l'attribut "distance_interface_element_eulerien_"
  //avec la valeur appropriee
  //Initialisation de "nombre_de_voisins_plus_proches"
  //Initialisation de "surface_interface_elements_voisins_"
  //REMARQUE : il faut s'assurer que chaque proc effectue le meme nombre de
  //boucles d'ou le test initial de la boucle while
  //
  ready_to_quit = (int) Process::mp_min(ready_to_quit);
  ArrOfInt elem_voisins;

  while (!ready_to_quit)
    {
      const int taille_liste = elem_proches->size_array();
      numero_voisinage++;

      //Pour le parallele : il faut reinitialiser la variable
      //a chaque nouveau front que l'on souhaite parcourir
      vu_localement = !seen;

      //Pour le parallele : il faut reinitialiser la variable
      //a chaque nouveau front que l'on souhaite parcourir
      for (int elem_loc=0; elem_loc<taille_liste; elem_loc++)
        {
          //Obtention des elements voisins de "elem"
          const int elem_glob = (*elem_proches)[elem_loc];
          elements_voisins(elem_glob,elem_voisins,domaine_VF);
          const int nb_elem_voisins = elem_voisins.size_array();
          const int distance_elem = distance_interface_element_eulerien_[elem_glob];

          //Boucle sur les elements voisins qui peuvent etre VIRTUELS
          for (int elem_voisin_loc=0; elem_voisin_loc<nb_elem_voisins; elem_voisin_loc++)
            {
              const int elem_voisin_glob = elem_voisins[elem_voisin_loc];

              //Initialisation locale de l'attribut "distance_interface_element_eulerien_"
              if (distance_interface_element_eulerien_[elem_voisin_glob]==-1)
                {
                  distance_interface_element_eulerien_[elem_voisin_glob] = numero_voisinage;
                  elem_eloignes->append_array(elem_voisin_glob);
                  if (elem_voisin_glob>=nb_elem) vu_localement[elem_voisin_glob-nb_elem]=seen;
                }

              //Initialisation des attributs "nombre_de_voisins_plus_proches"
              //et "surface_interface_elements_voisins_"
              //ATTENTION : cette petite boucle doit imperativement venir apres
              //l'initialisation locale de "distance_interface_element_eulerien_"
              const int distance_elem_voisin = distance_interface_element_eulerien_[elem_voisin_glob];

              if (distance_elem_voisin>distance_elem)
                {
                  //Il suffit d'ajouter 1 au nombre des voisins les plus proches
                  //de "elem_voisin_glob"
                  nombre_de_voisins_plus_proches_[elem_voisin_glob]++;

                  //Nous ne voulons tenir compte des surfaces que si "distance_elem_voisin==1"
                  if (distance_elem==0)
                    {
                      //pour s'assurer que l'on a pas fait d'erreur
                      assert(distance_elem_voisin==1);
                      surface_interface_elements_voisins_[elem_voisin_glob] += surface_interface_elements_voisins_[elem_glob];
                    }//fin du if sur "distance_elem==0"

                }//fin du if sur "distance_elem_voisin>distance_elem"

            }//fin du for sur "elem_voisin_loc"

        }//fin du for sur "elem_loc"

      //Pour le parallele : il faut s'assurer que les processus en cours
      //ont bien finis leurs boucles respectives
      distance_interface_element_eulerien_.echange_espace_virtuel();


      //Nous incrementons le compteur
      //Il faut faire les memes manip que pour les elements
      //a distance 0 de l'interface en prenant garde, en plus, de ne
      //pas rajouter des elements dans "elements_eloignes"
      //qui y ont deja ete place lors du traitement LOCAL, i.e.
      //avant l'usage de la fonction echange_espace_virtuel()
      nb_elem_parcourus+=elem_eloignes->size_array();

      for (int elem=nb_elem; elem<nb_elem_tot; elem++)
        {
          const int test_distance = (distance_interface_element_eulerien_[elem]==numero_voisinage);
          const int test_vu = (vu_localement[elem-nb_elem]!=seen);

          if ( test_distance && test_vu)
            {
              elem_eloignes->append_array(elem);
              nb_elem_parcourus++;
            }

        }//fin du for sur "elem"


      //
      //Swap des pointeurs pour ne pas couter trop cher en memoire
      //
      ArrOfInt *temporaire;
      elem_proches->resize_array(0);
      temporaire = &(*elem_proches);
      elem_proches = &(*elem_eloignes);
      elem_eloignes = &(*temporaire);

      //
      //Pour mettre a jour le test sur la boucle while
      //Un processeur est pret a quitter la boucle while() :
      //s' il a parcouru tous les elements du maillage local
      //OU s'il a atteint la distance imposee par l'utilisateur
      //
      ready_to_quit = (nb_elem_parcourus == nb_elem_tot) || (numero_voisinage == distance_utilisateur_);
      ready_to_quit = (int) Process::mp_min(ready_to_quit);

    }//fin du while

  assert(nb_elem_parcourus==nb_elem_tot || numero_voisinage==distance_utilisateur_);

  //Pour le parallele : il faut s'assurer que les processus en cours
  //ont bien finis leurs boucles respectives
  nombre_de_voisins_plus_proches_.echange_espace_virtuel();
  surface_interface_elements_voisins_.echange_espace_virtuel();

  //Inititalisation de l'attribut "plus_grande_distance_interface_element_eulerien_"
  //REMARQUE IMPORTANTE : cet attribut sera identique sur tous les procs
  plus_grande_distance_interface_element_eulerien_ = numero_voisinage ;
  assert_parallel(plus_grande_distance_interface_element_eulerien_);

  return 1;
}


//Soit un element "elem" (numero global) passe en argument.
//PRINCIPE
//Si le remaillage de l'interface ne conduit pas a une perte de volume dans "elem" : sortie
//Si le remaillage conduit a une perte de volume dans "elem " alors :
//   Si "elem" est a distance n > 1 de l'interface, on distribue equitablement le volume
//   perdu sur l'ensemble des elements voisins de "elem" a distance n-1
//   Si "elem" est a distance n = 1, on distribue le volume perdu sur l'ensemble
//   des voisins de "elem" a distance 0 de l'interface, proportionnellement a la surface
//   d'interface contenue dans ces memes voisins.
//   Si "elem" est a distance n = 0, on ne transporte pas le volume perdu : sortie
//Voisinage d'un sommet "s" := ensemble des elements contenant "s"
//Voisinage d'un element "K" := union des voisinages des sommets "si" de "K"
//Valeurs de sortie : 1 => pas de volume perdu
//                    2 => pas de transport du volume perdu pour une maille a distance 0 de l'interface
//                    3 => transport du volume perdu pour une maille a distance de l'interface >1
//                    4 => transport du volume perdu pour une maille a distance de l'interface =1
//                   -1 => situation problematique (une erreur s'est produite)
int  Remailleur_Collision_FT_Thomas::transport_volume_perdu_sur_element(const int elem,const Maillage_FT_Disc& maillage)
{
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,domaine_dis(maillage).valeur());
  const Domaine& domaine = domaine_VF.domaine();

  const int nb_elem = domaine.nb_elem();
  const int distance_interface_elem = distance_interface_element_eulerien_[elem];

  //int distance_interface_elem_voisin = -2;

  ArrOfInt voisins_a_distance_plus_petite;

  if (volume_perdu_[elem]==0) return 1; //pas de volume perdu => sortie

  //
  //Si nous arrivons ici, c'est que du volume a ete perdu.
  //Nous mettons en place l'algorithme de transport
  //

  //Construction de la liste des elements voisins de "elem"
  //a distance de l'interface strictement plus petite que la distance
  //separant "elem" de l'interface.
  elements_voisins(elem, voisins_a_distance_plus_petite, domaine_VF);
  elements_voisins_a_distance_plus_petite2(elem, voisins_a_distance_plus_petite);

  const int nb_voisins_a_distance_plus_petite = voisins_a_distance_plus_petite.size_array();


  switch (distance_interface_elem)
    {

    case -1 :

      Cerr << "Erreur Remailleur_Collision_FT_Thomas::transport_volume_perdu_sur_element()" << finl;
      Cerr << "La distance separant l'interface de l'element eulerien " << elem
           << " est strictement superieur a " << distance_utilisateur_ << finl;
      Cerr << "L'utilisateur a demande de ne considerer que les elements situes a une distance d'au plus "
           << distance_utilisateur_ << " de l'interface" << finl;
      Cerr << "Sortie du programme." << finl;
      Process::exit();

      // bfa : necessary to avoid a compilation error "[-Werror=implicit-fallthrough=]"
      return -1;

    case 0 :
      return 2;//pas de transport => sortie

    case 1 :

      assert(nb_voisins_a_distance_plus_petite!=0 || elem>=nb_elem);
      assert(nb_voisins_a_distance_plus_petite<=nombre_de_voisins_plus_proches_[elem]);
      if (elem<nb_elem) assert(nb_voisins_a_distance_plus_petite==nombre_de_voisins_plus_proches_[elem]);

      //Repartition du volume selon les principes enonces.
      //Utilisation de l'attribut "surface_interface_elements_voisins_" dont nous rappelons
      //la definition ici :
      //Aux elements "elem" a distance 0 correspond la surface d'interface qui les coupe
      //Aux elements "elem" a distance 1 correspond la somme des surfaces d'interface qui coupe
      //les voisins des "elem" a distance 0 de l'interface
      for (int elem_loc=0; elem_loc<nb_voisins_a_distance_plus_petite; elem_loc++)
        {
          const int elem_voisin = voisins_a_distance_plus_petite[elem_loc];


          assert(distance_interface_elem-distance_interface_element_eulerien_[elem_voisin]==1);
          volume_perdu_[elem_voisin]+=volume_perdu_[elem]*surface_interface_elements_voisins_[elem_voisin]
                                      /surface_interface_elements_voisins_[elem];

        }//fin du for sur "elem_loc"

      //remise a zero pour plus de surete
      volume_perdu_[elem]=0.;

      return 4;

    default :

      assert(nb_voisins_a_distance_plus_petite!=0 || elem>=nb_elem);
      assert(nb_voisins_a_distance_plus_petite<=nombre_de_voisins_plus_proches_[elem]);
      if (elem<nb_elem) assert(nb_voisins_a_distance_plus_petite==nombre_de_voisins_plus_proches_[elem]);

      //Repartition du volume selon les principes enonces.
      for (int elem_loc=0; elem_loc<nb_voisins_a_distance_plus_petite; elem_loc++)
        {
          const int elem_voisin = voisins_a_distance_plus_petite[elem_loc];

          assert(distance_interface_elem- distance_interface_element_eulerien_[elem_voisin]==1);
          volume_perdu_[elem_voisin]+=volume_perdu_[elem]/nombre_de_voisins_plus_proches_[elem];

        }//fin du for sur "elem_loc"

      //remise a zero pour plus de surete
      volume_perdu_[elem]=0.;

      return 3;

    }//fin du "switch"
}


//Fonction qui transporte le volume perdu lors du remaillage sur les elements
//euleriens a distance 0 de l'interface
//PRINCIPE
//Nous creons une liste chainee qui a partir d'une distance d'interface donnee "n"
//nous donnes l'ensemble des elements euleriens "elem" a distance "n"
//Ensuite, nous utilisons la fonction "transport_volume_perdu_sur_element(const int&, const Maillage_FT_Disc&)
//recursivement, en partant des elements a la plus grande distance de l'interface pour aller
//vers les elements a la plus petite distance de l'interface i.e. ceux a distance 0
//RETOUR : la valeur du volume perdu qui a ete transporte
//REMARQUE : l'attribut "distance_interface_element_eulerien_" doit etre initialise
double  Remailleur_Collision_FT_Thomas::transport_volume_perdu_sur_element(const Maillage_FT_Disc& maillage)
{
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,domaine_dis(maillage).valeur());
  const Domaine& domaine = domaine_VF.domaine();

  const int nb_elem = domaine.nb_elem();
  const int nb_elem_tot = domaine.nb_elem_tot();
  const int end_liste = -1;

  //Construction de la liste chainee distance->{ensemble des elements a cette distance}
  //Principe : Creation d'une premiere liste de taille "plus_grande_distance_interface_element_eulerien_"
  //           et initialisee a -1
  //           A la distance "n" sera stocke un element "elem" a distance "n" de l'interface
  //           Creation d'une deuxieme liste de taille "nb_elem_tot" initialisee a -1
  //           A la place "elem" est stocke le numero "elem2" d'un deuxieme element a distance "n"
  //           A la place "elem2" est stocke le numero "elem3" d'un deuxieme element a distance "n" etc
  //           La fin de la liste chainee est reperee par un -1

  IntTab distance_des_elements(plus_grande_distance_interface_element_eulerien_+1);
  IntTab suite_des_elements(nb_elem_tot);
  distance_des_elements = end_liste;
  suite_des_elements = end_liste;

  assert(nb_elem_tot==distance_interface_element_eulerien_.dimension_tot(0));

  //Boucle sur "distance_interface_element_eulerien_" afin de modifier nos listes
  for (int elem=0; elem<nb_elem_tot; elem++)
    {
      const int distance = distance_interface_element_eulerien_[elem];

      if (distance==-1) continue; //on passe tout de suite a l'element suivant
      assert(distance<=plus_grande_distance_interface_element_eulerien_);

      suite_des_elements(elem)=distance_des_elements(distance);
      distance_des_elements(distance)=elem;

    }//fin du for sur "elem"


  //
  //Transport du volume perdu proprement dit
  //L'attribut "plus_grande_distance_interface_element_eulerien_" etant identique
  //sur tous les procs, il n'y a pas de precautions particulieres a prendre
  //
  for (int distance=plus_grande_distance_interface_element_eulerien_; distance!=-1; distance--)
    {
      //Parcours de tous les elements "elem" situes a la distance "distance" de l'interface
      //et transport du volume perdu
      for (int elem=distance_des_elements(distance); elem!=-1; elem=suite_des_elements(elem))
        transport_volume_perdu_sur_element(elem,maillage);

      //Pour le parallele : il faut s'assurer que les processus en cours
      //ont bien finis leurs boucles respectives
      volume_perdu_.echange_espace_virtuel();

    }//fin du for sur "distance"

  //
  //Calcul du volume perdu transporte sur les elements
  //a distance 0 de l'interface
  //
  double volume_perdu_transporte = 0.;

  for (int elem=0; elem<nb_elem; elem++)
    if (!distance_interface_element_eulerien_[elem])
      volume_perdu_transporte += volume_perdu_[elem];

  return volume_perdu_transporte;
}



//Fonction qui donne la liste des elements voisins a un element "elem" (numero global) donne
//Voisinage d'un sommet "s" := ensemble des elements du maillage contenant "s"
//Voisinage d'un element "elem" := union des voisinages des sommets "si" de "elem"
//REMARQUE : un element de la liste "liste_voisins" n'est present qu'une seule fois
//REMARQUE : les attributs "voisinage_sommet_" et "next_elem_" doivent etre initialises
//REMARQUE : le parametre "liste_voisins" doit avoir une taille initiale differente de "0"
//(astuce de gestion memoire avec l'objet ArrOfInt)
//REMARQUE : pour le parallele, "liste_voisins" contient egalement des elements voisins virtuels
//VOIR : Remailleur_Collision_FT_Thomas::voisinage_sommets()
int Remailleur_Collision_FT_Thomas::elements_voisins(const int elem,
                                                     ArrOfInt& liste_voisins,
                                                     const Domaine_dis_base& un_domaine_dis) const
{
  //  const Maillage_FT_Disc& maillage = ;
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,un_domaine_dis);
  const Domaine& domaine = domaine_VF.domaine();

  const int nb_som_elem = domaine.nb_som_elem();
  const int end_liste = -1;

  //Tableau des sommets par element
  const IntTab& elem_sommets = domaine.les_elems();


  liste_voisins.resize_array(0);

  //Calcul du voisinage
  for (int som=0; som<nb_som_elem; som++)
    {
      const int som_global = domaine.get_renum_som_perio(elem_sommets(elem,som));

      //On parcourt la liste des elements "elem_voisin" qui contiennent "som_global"
      for (int p=voisinage_sommet_[som_global]; p!=end_liste; p=next_elem_[p])
        {
          const int elem_voisin = p/nb_som_elem;
          assert(som_global==domaine.get_renum_som_perio(elem_sommets(elem_voisin,p%nb_som_elem)));
          if (!tmp_flag_elements_.testsetbit(elem_voisin))
            liste_voisins.append_array(elem_voisin);
        }
    }

  const int n = liste_voisins.size_array();
  for (int i = 0; i < n; i++)
    tmp_flag_elements_.clearbit(liste_voisins[i]);

  return 1;
}


//Fonction qui donne la liste des elements voisins a un element "elem" (numero global) donne et
//situes a une distance de l'interface strictement plus petite que la distance separant "elem"
//de l'interface. (normalement la difference des deux distances doit etre de 1)
// Changement (BM): on suppose que liste_voisins contient le resultat de la fonction elements_voisins()
//Voisinage d'un sommet "s" := ensemble des elements du maillage contenant "s"
//Voisinage d'un element "elem" := union des voisinages des sommets "si" de "elem"
//Entree : l'element "elem" dont on veut la liste des voisins
//Entree ET Sortie : la liste des voisins a distance de l'interface plus petite que la distance
//separant "elem" de l'interface
//REMARQUE : un element de la liste "liste_voisins" n'est present qu'une seule fois
//REMARQUE : l'attribut "distance_interface_element_eulerien_" doit etre initialise
//REMARQUE : pour le parallele, "liste_voisins" contient egalement des elements voisins virtuels
//VOIR : Remailleur_Collision_FT_Thomas::voisinage_sommets() et Remailleur_Collision_FT_Thomas::elements_voisins()
int  Remailleur_Collision_FT_Thomas::elements_voisins_a_distance_plus_petite2(const int elem,
                                                                              ArrOfInt& liste_voisins) const
{
  const int distance_interface_elem = distance_interface_element_eulerien_[elem];
  const int n = liste_voisins.size_array();
  int j = 0;

  for (int i = 0; i < n; i++)
    {
      const int elem_voisin = liste_voisins[i];
      const int distance_interface_elem_voisin = distance_interface_element_eulerien_[elem_voisin];

      const bool is_plus_petit = distance_interface_elem_voisin<distance_interface_elem;
      const bool is_positif = distance_interface_elem_voisin>-1;

      if (is_positif && is_plus_petit)
        {
          assert(distance_interface_elem-distance_interface_elem_voisin==1);
          liste_voisins[j] = elem_voisin;
          j++;
        }
    }
  liste_voisins.resize_array(j);

  return 1;
}




//Fonction qui donne le nombre d'elements voisins a un element "elem" (numero global) donne
//Voisinage d'un sommet "s" := ensemble des elements du maillage contenant "s"
//Voisinage d'un element "elem" := union des voisinages des sommets "si" de "elem"
//REMARQUE : pour le parallele, cette fonction comptabilise egalement les elements voisins virtuels
//VOIR : Remailleur_Collision_FT_Thomas::elements_voisins()
int  Remailleur_Collision_FT_Thomas::nb_elements_voisins(const int elem, const Domaine_dis_base& un_domaine_dis) const
{
  ArrOfInt liste_voisins;
  elements_voisins(elem,liste_voisins,un_domaine_dis);
  return liste_voisins.size_array();
}

//Fonction qui donne le nombre d'elements voisins a un element "elem" (numero global) donne
//situes a une distance de l'interface plus petite que la distance separant "elem" de
//l'interface.
//Voisinage d'un sommet "s" := ensemble des elements du maillage contenant "s"
//Voisinage d'un element "elem" := union des voisinages des sommets "si" de "elem"
//REMARQUE : pour le parallele, cette fonction comptabilise egalement les elements voisins virtuels
//VOIR : Remailleur_Collision_FT_Thomas::elements_voisins_a_distance_plus_petite()
int Remailleur_Collision_FT_Thomas::nb_elements_voisins_a_distance_plus_petite(const int elem, const Domaine_dis_base& un_domaine_dis) const
{
  ArrOfInt liste_voisins;
  elements_voisins(elem,liste_voisins,un_domaine_dis);
  elements_voisins_a_distance_plus_petite2(elem,liste_voisins);
  return liste_voisins.size_array();
}

int Remailleur_Collision_FT_Thomas::nb_elements_voisins_a_distance_plus_petite(const int elem) const
{
  assert(nombre_de_voisins_plus_proches_.dimension_tot(0)>elem);
  return nombre_de_voisins_plus_proches_[elem];
}


//Fonction qui assure le remaillage de l'interface.
//Cette fonction calcule, en outre, le volume perdu dans l'operation de remaillage,
//et le redistribue de facon a ce que le volume de occupe par l'indicatrice de phase
//soit identique avant et apres l'operation de remaillage topologique.
//Remarque : on suppose que les valeurs de l'interface sont connues dans
//chaque element du maillage eulerien.
int Remailleur_Collision_FT_Thomas::traite_RuptureCoalescenceInterfaces_Conservatif(Maillage_FT_Disc& maillage,
                                                                                    Champ_base& indicatrice)
{
  Journal() << "Remailleur_Collision_FT_Thomas::traite_RuptureCoalescenceInterfaces_Conservatif" << finl;
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,domaine_dis(maillage).valeur());
  const Domaine& domaine = domaine_VF.domaine();

  const DoubleVect& volumes_elements_euleriens = domaine_VF.volumes();

  const int nb_elements_euleriens = domaine.nb_elem();
  const int nb_elements_euleriens_tot = domaine.nb_elem_tot();

  int ok_mise_a_jour = 0, ok_remailleur = 0;

  const double temps = maillage.temps();
  double volume_initial = 0.;
  double volume_final = 0.;
  double le_volume_perdu = 0.;
  double volume_perdu_elements = 0.;
  double volume_perdu_sommets = 0.;

  ArrOfDoubleFT variation_volume;

  Champ_base& indicatrice_finale = indicatrice; //pour l'initialisation slt

  assert(indicatrice.valeurs().dimension_tot(0)==nb_elements_euleriens_tot);
  assert(volumes_elements_euleriens.size_totale()==nb_elements_euleriens_tot);

  //On dimensionne les donnees qui nous interesse
  if (!est_dimensionne_) initialiser_data(maillage);
  assert(est_dimensionne_);

  //On initialise l'attribut "volume_perdu_" et
  //mise a jour entre procs
  //On calcule le volume initial avant remaillage
  volume_perdu_ = indicatrice.valeurs();
  volume_perdu_.echange_espace_virtuel();

  for (int elem=0; elem<nb_elements_euleriens; elem++)
    volume_initial += indicatrice.valeurs()[elem]*volumes_elements_euleriens[elem];

  volume_initial = Process::mp_sum(volume_initial);

  //On remaille topologiquement puis conservativement
  //pour obtenir une surface lisse, i.e. sans angle
  //trop prononce
  ok_remailleur =
    Remailleur_Collision_FT_Juric::traite_RuptureCoalescenceInterfaces(maillage, indicatrice);

  remaillage_FT(maillage).remaillage_local_interface(temps, maillage);

  //On modifie l'attribut "volume_perdu_" et on calcule
  //le volume total perdu dans l'operation de remaillage
  //REMARQUE : le maillage a change apres le remaillage
  //           et donc l'indicatrice egalement
  volume_perdu_ -= maillage.equation_transport().get_update_indicatrice().valeurs();
  volume_perdu_.echange_espace_virtuel();

  for (int elem=0; elem<nb_elements_euleriens_tot; elem++)
    volume_perdu_[elem] *= volumes_elements_euleriens[elem];

  volume_perdu_.echange_espace_virtuel();

  for (int elem=0; elem<nb_elements_euleriens; elem++)
    le_volume_perdu+=volume_perdu_[elem];

  le_volume_perdu = Process::mp_sum(le_volume_perdu);


  //Nous executons la suite de l'algo a savoir :
  //calcul des voisinages
  //calcul des distances a l'interface
  //transport des volumes perdus
  //calcul des volumes perdus par sommet de l'interface
  //L'etat du maillage est suppose etre au moins : PARCOURU
  variation_volume.resize_array(maillage.nb_sommets());
  variation_volume = 0.;

  maillage.parcourir_maillage();
  ok_mise_a_jour = mettre_a_jour_data(maillage);

  //Verification du bon deroulement des evenements
  int ok = ok_remailleur*ok_mise_a_jour;
  assert(ok);

  //Pour le debug
  if (tester_) tester(maillage);

  //On transporte le volume perdu jusqu'aux elements a distance 0 de l'interface
  //On en profite pour calculer le volume total transporte
  volume_perdu_elements = transport_volume_perdu_sur_element(maillage);
  volume_perdu_elements = Process::mp_sum(volume_perdu_elements);

  //On transporte le volume perdu aux sommets du maillage de l'interface
  //On en profite pour calculer le volume total redistribue sur les sommets
  volume_perdu_sommets = transport_volume_perdu_sur_sommet(variation_volume, maillage);
  volume_perdu_sommets = Process::mp_sum(volume_perdu_sommets);

  //On deplace les sommets de l'interface pour corriger le volume
  //puis on calcule le volume finalement obtenu
  //PARAMETRE : le tableau "variation_volume" calcule ci-dessus
  //PRINCIPE : lissage puis deplacement des sommets
  remaillage_FT(maillage).lisser_dvolume(maillage,variation_volume,5);
  variation_volume*=-1.;
  remaillage_FT(maillage).corriger_volume(maillage,variation_volume);

  indicatrice_finale = maillage.equation_transport().get_update_indicatrice();
  for (int elem=0; elem<nb_elements_euleriens; elem++)
    volume_final += indicatrice_finale.valeurs()[elem]*volumes_elements_euleriens[elem];

  volume_final = Process::mp_sum(volume_final);

  //On affiche la perte de volume calculee
  Journal() << "Perte de volume due au remaillage  : " << le_volume_perdu << " "
            << volume_perdu_elements << " " << volume_perdu_sommets
            << "\nVolume initial : " << volume_initial
            << "  volume final : " << volume_final << finl;

  return ok;
}

//Fonction qui a chaque sommet du maillage de l'interface associe une partie du
//volume perdu dans l'element "elem" lors du remaillage
//Principe : si du volume a ete perdu sur la maille eulerienne "elem" traversee par l'interface
//           alors une fraction de ce volume est distribue sur chaque sommet "s" de l'interface coupant "elem".
//           Pour cela, pour chaque facette "fa7" coupant "elem", nous calculons un volume perdu surfacique
//           "volume_perdu_surfacique" qui est le volume perdu dans "elem" divise par la surface de l'intersection de
//           "fa7" dans "elem", nous divisons "volume_perdu_surfacique" par la surface totale d'intersection entre
//           l'interface et "elem", nous recuperons les sommets "si" de "fa7", recuperons le barycentre de l'intersection
//           "G" de la "fa7" dans "elem", et distribuons aux sommets "si" de fa7 (meme si l'un des "si" n'appartient pas a "elem")
//           "volume_perdu_surfacique" proportionnellement aux coordonnees barycentriques de "G", definies
//           par rapport aux "si" de "fa7".
//REMARQUE : le volume perdu total est suppose avoir ete transfere (d'une maniere ou d'une autre) sur les mailles
//           euleriennes traversees par l'interface
//REMARQUE : l'attribut "volume_perdu_" est suppose avoir ete initialise
//REMARQUE : l'attribut "distance_interface_element_eulerien_" est suppose avoir ete initialise
int  Remailleur_Collision_FT_Thomas::transport_volume_perdu_sur_sommet(const int elem, ArrOfDouble& volume_perdu_sommet, const Maillage_FT_Disc& maillage) const
{
  const Intersections_Elem_Facettes& intersections =
    maillage.intersections_elem_facettes();

  const ArrOfInt& index_elem = intersections.index_elem();
  const IntTab& facettes = maillage.facettes();

  const int distance_interface_elem = distance_interface_element_eulerien_[elem];

  assert(volume_perdu_.dimension_tot(0)==domaine_dis(maillage).valeur().domaine().nb_elem_tot());

  const double volume_perdu_elem = volume_perdu_[elem];
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();

  //Boucle sur les faces qui traversent l'element:
  //ATTENTION : en parallele, le tableau "index_elem" est dimensionne
  //a nb_elem() et pas a nb_elem_tot()
  assert(elem<domaine_dis(maillage).valeur().domaine().nb_elem());
  int index = index_elem[elem];

  const double somme_coefs = surface_interface_elements_voisins_[elem];
  const double inv_somme_coefs = (somme_coefs > 0.) ? (1. / somme_coefs) : 0.;

  while (index >= 0 && distance_interface_elem==0)
    {
      const Intersections_Elem_Facettes_Data& data =
        intersections.data_intersection(index);

      //Calcul du volume perdu associe a la facette "facette"
      // i.e. volume perdu*(surface de "facette" dans "elem")/(surface totale de l'interface dans "elem")
      const int facette = data.numero_facette_;
      assert(facette<=maillage.nb_facettes());
      const double coeff_surface = data.fraction_surface_intersection_ * surface_facettes[facette] * inv_somme_coefs;

      const double volume_perdu_surfacique = volume_perdu_elem*coeff_surface;

      for (int som_loc=0; som_loc<Objet_U::dimension; som_loc++)
        {
          //Recuperation d'un sommet de "facette" et des coodonnees barycentrique
          //du barycentre du centre de gravite de l'intersection entre l'interface
          //et "elem"
          const int sommet = facettes(facette,som_loc);
          assert(sommet<maillage.nb_sommets());

          //Calcul du volume perdu associe au sommet
          volume_perdu_sommet[sommet]+=volume_perdu_surfacique*data.barycentre_[som_loc];

        }//fin du for sur "som_loc"

      index = data.index_facette_suivante_;

    }//fin du while

  return 1;
}



//Fonction qui a chaque sommet du maillage de l'interface associe une partie du
//volume perdu lors du remaillage
//Principe : Appel la fonction  Remailleur_Collision_FT_Thomas::transport_volume_perdu_sur_sommet(const int&, ArrOfDouble &, const Maillage_FT_Disc&)
//           pour chaque element "elem" du maillage eulerien.
//RETOUR : renvoie le volume perdu total qui a ete transporte jusqu'aux sommets du maillage de l'interface
double Remailleur_Collision_FT_Thomas::transport_volume_perdu_sur_sommet(ArrOfDouble& volume_perdu_sommet, const Maillage_FT_Disc& maillage) const
{
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,domaine_dis(maillage).valeur());
  const Domaine& domaine = domaine_VF.domaine();

  const int nb_elem = domaine.nb_elem();
  const int nb_sommets_interface = maillage.nb_sommets();
  // int ok = -1; //pour les assert()

  //Initialisation de "volume_perdu_sommet"
  if (volume_perdu_sommet.size_array()!=nb_sommets_interface)
    volume_perdu_sommet.resize_array(nb_sommets_interface);

  volume_perdu_sommet = 0.;


  //Appel de la fonction  Remailleur_Collision_FT_Thomas::transport_volume_perdu_sur_sommet(const int&, ArrOfDouble &, const Maillage_FT_Disc&)
  for (int elem=0; elem<nb_elem; elem++)
    {
#ifndef NDEBUG
      int ok=
#endif
        transport_volume_perdu_sur_sommet(elem,volume_perdu_sommet,maillage);
      assert(ok==1);
    }

  //Pour le parallele
  maillage.desc_sommets().collecter_espace_virtuel(volume_perdu_sommet, MD_Vector_tools::EV_SOMME);

  //
  //Calcul du volume perdu transporte sur les
  //sommets de l'interface
  //
  double volume_perdu_transporte = 0.;

  for (int som=0; som<nb_sommets_interface; som++)
    if (!maillage.sommet_virtuel(som)) volume_perdu_transporte+=volume_perdu_sommet[som];

  return volume_perdu_transporte;
}


//Fonction qui renvoie la surface occupee par l'interface
//dans un element "elem" REEL donne
//REMARQUE : renvoie 0 si l'element n'est pas traverse par l'interface
double  Remailleur_Collision_FT_Thomas::surface_intersection(const int elem, const Maillage_FT_Disc& maillage) const
{
  const Intersections_Elem_Facettes& intersections =
    maillage.intersections_elem_facettes();

  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const ArrOfInt& index_elem = intersections.index_elem();

  double surface_dans_elem = 0.;

  // Boucle sur les faces qui traversent l'element:
  //ATTENTION : en parallele, le tableau "index_elem" est dimensionne
  //a nb_elem() et pas a nb_elem_tot()
  assert(elem<domaine_dis(maillage).valeur().domaine().nb_elem());
  int index = index_elem[elem];

  while (index >= 0)
    {
      const Intersections_Elem_Facettes_Data& data =
        intersections.data_intersection(index);

      surface_dans_elem += data.fraction_surface_intersection_ * surface_facettes[data.numero_facette_];
      index = data.index_facette_suivante_;

    }//fin du while

  return surface_dans_elem;
}


//Fonctions destinees a tester les autres fonctions membre de la classe
void  Remailleur_Collision_FT_Thomas::tester(const Maillage_FT_Disc& maillage) const
{
  tester_interface(maillage);
  tester_distance(maillage);
  tester_voisinage(maillage);
  tester_liste(maillage);
  tester_transport(maillage);

  DoubleTab copie_volume_perdu(volume_perdu_);
  copie_volume_perdu.echange_espace_virtuel();

  tester_transport_complet(maillage, copie_volume_perdu);
  tester_volume_par_sommet(maillage, copie_volume_perdu);

  if (Process::je_suis_maitre())
    {
      Cerr << "Fin des tests." << finl;
      Cerr << "Sortie du programme." << finl;
    }

  barrier();
  Process::exit();
}

void  Remailleur_Collision_FT_Thomas::tester_interface(const Maillage_FT_Disc& maillage) const
{
  Nom nom("Elements_euleriens_proc_");
  nom+=Process::me();
  nom+=".txt";
  SFichier fichier_elements_euleriens(nom.getChar());

  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,domaine_dis(maillage).valeur());
  const Domaine& domaine = domaine_VF.domaine();

  const IntTab& facettes = maillage.facettes();
  const ArrOfInt& sommet_elem = maillage.sommet_elem();

  const DoubleVect& volumes_elements_euleriens = domaine_VF.volumes();
  const DoubleTab& indicatrice =
    maillage.equation_transport().inconnue().valeurs();

  const int nb_elements_euleriens = domaine.nb_elem();

  //On va afficher les elements traverses par l'interface
  //Les facettes coupant ces elements
  //Les sommets de ces facettes
  //Les sommets de ces facettes appartenant aux elements coupes par l'interface
  const Intersections_Elem_Facettes& intersections =
    maillage.intersections_elem_facettes();

  const ArrOfInt& index_element = intersections.index_elem();
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();

  for (int elem=0; elem<nb_elements_euleriens; elem++)
    {
      fichier_elements_euleriens << "Element eulerien reel numero : " << elem << finl;
      fichier_elements_euleriens << "Taux d'occupation de l'interface : " << indicatrice[elem] << finl;
      fichier_elements_euleriens << "Volume de l'element " << elem << " : " << volumes_elements_euleriens[elem] << finl;
      fichier_elements_euleriens << "Volume occupee par l'interface : " << indicatrice[elem]*volumes_elements_euleriens[elem] << finl;
      fichier_elements_euleriens << "L'element " << elem << " est coupe par les facettes : " << finl;

      double surface_intersection_totale = 0.;

      // Boucle sur les faces qui traversent l'element:
      int index = index_element[elem];

      while (index >= 0)
        {
          const Intersections_Elem_Facettes_Data& data =
            intersections.data_intersection(index);

          const int facette = data.numero_facette_;

          double surface_dans_elem = data.fraction_surface_intersection_ * surface_facettes[facette];
          surface_intersection_totale += surface_dans_elem;

          if (surface_dans_elem!=0.)
            {
              fichier_elements_euleriens << facette << " de surface d'intersection "
                                         << surface_dans_elem << finl;

              fichier_elements_euleriens << "De sommets : ( ";

              for (int som_loc=0; som_loc<Objet_U::dimension; som_loc++)
                {
                  //Recuperation d'un sommet de "facette" et des coodonnees barycentrique
                  //du barycentre du centre de gravite de l'intersection entre l'interface
                  //et "elem"
                  const int sommet = facettes(facette,som_loc);
                  assert(sommet<maillage.nb_sommets());

                  if (sommet_elem[sommet]==elem)
                    fichier_elements_euleriens << "{" << (int) sommet << "," << (int) 1 << "} ";
                  else
                    fichier_elements_euleriens << "{" << (int) sommet << "," << (int) 0 << "} ";

                }//fin du for sur "som_loc"

              fichier_elements_euleriens << ")" << finl;

            }//fin du if

          index = data.index_facette_suivante_;

        }//fin du while

      fichier_elements_euleriens << "Surface d'intersection totale : " << surface_intersection_totale << finl;
      fichier_elements_euleriens << "Erreur fonction  Remailleur_Collision_FT_Thomas::surface_intersection()"
                                 << surface_intersection_totale-surface_intersection(elem,maillage) << finl;

      fichier_elements_euleriens << "**************************************************************" << finl;
      fichier_elements_euleriens << "**************************************************************" << finl;

    }//fin du for sur "elem"

}


void  Remailleur_Collision_FT_Thomas::tester_voisinage(const Maillage_FT_Disc& maillage) const
{
  Nom nom_sommets("Voisinage_sommets_proc_");
  nom_sommets+=Process::me();
  nom_sommets+=".txt";

  Nom nom_elements("Voisinage_elements_proc_");
  nom_elements+=Process::me();
  nom_elements+=".txt";

  SFichier fichier_voisinage_sommet(nom_sommets.getChar());
  SFichier fichier_voisinage_element(nom_elements.getChar());

  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,domaine_dis(maillage).valeur());
  const Domaine& domaine = domaine_VF.domaine();

  //Affichage des voisinages des sommets
  const int nb_som_tot = domaine.nb_som_tot();
  const int nb_elem_tot = domaine.nb_elem_tot();
  const int nb_som_elem = domaine.nb_som_elem();
  const int end_liste = -1;

  fichier_voisinage_sommet << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  fichier_voisinage_sommet << "Nombre d'elements reels : " << domaine.nb_elem() << finl;
  fichier_voisinage_sommet << "Nombre d'elements virtuels : " << nb_elem_tot << finl;
  fichier_voisinage_sommet << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  fichier_voisinage_sommet << finl;
  fichier_voisinage_sommet << finl;

  for (int som=domaine.get_renum_som_perio(0); som<nb_som_tot; som++) // som=domaine.get_renum_som_perio(som+1))
    {
      fichier_voisinage_sommet << "Affichage des elements voisins du sommet " << som << " : {{{ ";

      for (int p=voisinage_sommet_[som]; p!=end_liste; p=next_elem_[p])
        {
          const int elem = p/nb_som_elem;

          assert(elem<nb_elem_tot);
          assert(som==domaine.get_renum_som_perio(domaine.les_elems()(elem,p%nb_som_elem)));
          fichier_voisinage_sommet << elem << " ";
        }

      fichier_voisinage_sommet << "}}}" << finl;
      fichier_voisinage_sommet << "*************************************************************" << finl;
      fichier_voisinage_sommet << "*************************************************************" << finl;
      fichier_voisinage_sommet << finl;

    }//fin du for sur "som"



  int taille_liste_voisins_plus_proches = -1;
  ArrOfInt elem_voisins;

  for (int elem=0; elem<nb_elem_tot; elem++)
    {
      elements_voisins(elem, elem_voisins, domaine_VF);
      const int taille_liste_voisins = elem_voisins.size_array();

      //
      //Affichage des voisins d'un element "elem"
      //
      fichier_voisinage_element << "Affichage des elements voisins de l'element " << elem << " : {{{ ";

      for (int elem_loc=0; elem_loc<taille_liste_voisins; elem_loc++)
        fichier_voisinage_element << elem_voisins[elem_loc] << " ";

      fichier_voisinage_element << "}}}" << finl;

      elements_voisins_a_distance_plus_petite2(elem, elem_voisins);

      taille_liste_voisins_plus_proches = elem_voisins.size_array();
      assert(taille_liste_voisins_plus_proches<=nombre_de_voisins_plus_proches_[elem]);
      if (elem<domaine.nb_elem()) assert(taille_liste_voisins_plus_proches==nombre_de_voisins_plus_proches_[elem]);

      //
      //Affichage des voisins + proche de l'interface que "elem
      //
      fichier_voisinage_element << "Affichage des elements voisins de l'element " << elem
                                << " a distance plus petite que : "
                                << distance_interface_element_eulerien_[elem] << " : {{{ ";

      for (int elem_loc=0; elem_loc<taille_liste_voisins_plus_proches; elem_loc++)
        {
          const int elem_voisin_plus_proche = elem_voisins[elem_loc];
          fichier_voisinage_element << "(" << elem_voisin_plus_proche << ","
                                    << distance_interface_element_eulerien_[elem_voisin_plus_proche]
                                    << ") ";
        }//fin du for sur "elem_loc"

      fichier_voisinage_element << "}}}" << finl;

      fichier_voisinage_element << "Nombre de voisins plus proches : " << nombre_de_voisins_plus_proches_[elem] << finl;

      fichier_voisinage_element << "*************************************************************" << finl;
      fichier_voisinage_element << "*************************************************************" << finl;
      fichier_voisinage_element << finl;

    }//fin du for sur "elem"
}



void  Remailleur_Collision_FT_Thomas::tester_distance(const Maillage_FT_Disc& maillage) const
{
  Nom nom("Distance_interface_proc_");
  nom+=Process::me();
  nom+=".txt";

  SFichier fichier_distance(nom.getChar());

  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,domaine_dis(maillage).valeur());
  const Domaine& domaine = domaine_VF.domaine();

  const int nb_elem_tot = domaine.nb_elem_tot();

  fichier_distance << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  fichier_distance << "Nombre d'elements reels : " << domaine.nb_elem() << finl;
  fichier_distance << "Nombre d'elements virtuels : " << nb_elem_tot << finl;
  fichier_distance << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  fichier_distance << finl;
  fichier_distance << finl;


  fichier_distance << "Affichage de la plus grande distance a l'interface." << plus_grande_distance_interface_element_eulerien_;
  fichier_distance << finl;
  fichier_distance << "*************************************************************" << finl;
  fichier_distance << "*************************************************************" << finl;
  fichier_distance << finl;

  for (int elem=0; elem<nb_elem_tot; elem++)
    {
      const int distance = distance_interface_element_eulerien_[elem];
      assert(distance<=plus_grande_distance_interface_element_eulerien_);
      fichier_distance << "Affichage de l'element " << elem << " ";
      fichier_distance << "et de sa distance a l'interface : " << distance << finl;
    }
}


void  Remailleur_Collision_FT_Thomas::tester_liste(const Maillage_FT_Disc& maillage) const
{
  Nom nom("Liste_chainee_proc_");
  nom+=Process::me();
  nom+=".txt";

  SFichier fichier_liste(nom.getChar());

  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,domaine_dis(maillage).valeur());
  const Domaine& domaine = domaine_VF.domaine();

  const int nb_elem_tot = domaine.nb_elem_tot();

  //
  //Test du chainage dans la fonction  Remailleur_Collision_FT_Thomas::transport_volume_perdu_sur_element(const Maillage_FT_Disc&)
  //
  IntTab distance_des_elements(plus_grande_distance_interface_element_eulerien_+1);
  IntTab suite_des_elements(nb_elem_tot);
  distance_des_elements = -1;
  suite_des_elements = -1;

  assert(nb_elem_tot==distance_interface_element_eulerien_.dimension_tot(0));


  fichier_liste << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  fichier_liste << "Nombre d'elements reels : " << domaine.nb_elem() << finl;
  fichier_liste << "Nombre d'elements virtuels : " << nb_elem_tot << finl;
  fichier_liste << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  fichier_liste << finl;
  fichier_liste << finl;


  //Boucle sur "distance_interface_element_eulerien_" afin de modifier nos listes
  for (int elem=0; elem<nb_elem_tot; elem++)
    {
      const int distance = distance_interface_element_eulerien_[elem];
      if (distance==-1) continue; //on passe directement a l'element suivant
      assert(distance<=plus_grande_distance_interface_element_eulerien_);

      suite_des_elements(elem)=distance_des_elements(distance);
      distance_des_elements(distance)=elem;

    }//fin du for sur "elem"


  fichier_liste << "Affichage de la premiere liste chainee" << finl;
  fichier_liste << finl;
  for (int distance=0; distance<plus_grande_distance_interface_element_eulerien_+1; distance++)
    fichier_liste << "Distance : " << distance << " element : " << distance_des_elements(distance) << finl;

  fichier_liste << finl;
  fichier_liste << "*************************************************************" << finl;
  fichier_liste << "*************************************************************" << finl;
  fichier_liste << finl;


  fichier_liste << "Affichage de la deuxieme liste chainee" << finl;
  fichier_liste << finl;

  for (int elem=0; elem<nb_elem_tot; elem++)
    fichier_liste << "Precedent : " << elem << " suivant : " << suite_des_elements(elem) << finl;

  fichier_liste << finl;
  fichier_liste << "*************************************************************" << finl;
  fichier_liste << "*************************************************************" << finl;
  fichier_liste << finl;

  fichier_liste << "Test du parcours de la liste chainee" << finl;
  fichier_liste << finl;


  for (int distance=plus_grande_distance_interface_element_eulerien_; distance!=-1; distance--)
    {
      //Parcours de tous les elements "elem" situes a la distance "distance" de l'interface
      //et transport du volume perdu
      fichier_liste << "Affichage des elements a distance " << distance << " de l'interface." << finl;
      fichier_liste << "{ ";

      for (int elem=distance_des_elements(distance); elem!=-1; elem=suite_des_elements(elem))
        fichier_liste << elem << " ";

      fichier_liste << "}" << finl;
      fichier_liste << "-------------------------------------------------------------" << finl;
      fichier_liste << finl;

    }//fin du for sur distance
}


void  Remailleur_Collision_FT_Thomas::tester_transport(const Maillage_FT_Disc& maillage) const
{
  Nom nom("Transport_proc_");
  nom+=Process::me();
  nom+=".txt";

  SFichier fichier_transport(nom.getChar());

  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,domaine_dis(maillage).valeur());
  const Domaine& domaine = domaine_VF.domaine();

  const int nb_elem = domaine.nb_elem();
  const int nb_elem_tot = domaine.nb_elem_tot();

  if (volume_perdu_.dimension_tot(0) != nb_elem_tot)
    {
      Cerr << "Erreur  Remailleur_Collision_FT_Thomas::transport_volume_perdu_sur_element()" << finl;
      Cerr << "Le tableau volume_distribue_ n'est pas correctement initialise." << finl;
      Cerr << "Sortie du programme." << finl;
      Process::exit();
    }

  if (distance_interface_element_eulerien_.dimension_tot(0) != nb_elem_tot)
    {
      Cerr << "Erreur  Remailleur_Collision_FT_Thomas::transport_volume_perdu_sur_element()" << finl;
      Cerr << "Le tableau distance_interface_element_eulerien_ n'est pas correctement initialise." << finl;
      Cerr << "Sortie du programme." << finl;
      Process::exit();
    }



  ArrOfInt voisins_a_distance_plus_petite(100);
  voisins_a_distance_plus_petite.resize_array(0);


  fichier_transport << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  fichier_transport << "Nombre d'elements reels : " << nb_elem << finl;
  fichier_transport << "Nombre d'elements virtuels : " << nb_elem_tot << finl;
  fichier_transport << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  fichier_transport << finl;
  fichier_transport << finl;

  for (int elem=0; elem<nb_elem_tot; elem++)
    {
      const int distance_interface_elem = distance_interface_element_eulerien_[elem];

      if (volume_perdu_[elem]==0)
        {
          fichier_transport << "Pas de volume perdu pour l'element " << elem << finl;
          fichier_transport << "*******************************************************" << finl;
          fichier_transport << "*******************************************************" << finl;
          fichier_transport << finl;
          continue;
        }

      elements_voisins(elem, voisins_a_distance_plus_petite,domaine_VF);
      elements_voisins_a_distance_plus_petite2(elem, voisins_a_distance_plus_petite);
      const int nb_voisins_a_distance_plus_petite = voisins_a_distance_plus_petite.size_array();

      ArrOfDouble surface_par_element(nb_voisins_a_distance_plus_petite);

      //Recuperation de la surface totale des facettes traversant
      //l'ensemble des elements a distance 0 de l'interface, ainsi que la surface
      //des facettes traversant chaque element a distance 0 de l'inteface
      double surface_totale = 0.;
      //      int distance_interface_elem_voisin = -2;


      switch(distance_interface_elem)
        {

        case -1 :

          fichier_transport << "L'element " << elem
                            << " est a une distance de l'interface superieure a " << distance_interface_elem << finl;
          fichier_transport << "*******************************************************" << finl;
          fichier_transport << "*******************************************************" << finl;
          fichier_transport << finl;

          break;



        case 0 :

          fichier_transport << "Pas de transport pour l'element " << elem
                            << " a distance " << distance_interface_elem << " de l'interface" << finl;
          fichier_transport << "*******************************************************" << finl;
          fichier_transport << "*******************************************************" << finl;
          fichier_transport << finl;

          break;



        case 1 : //cas ou distance_interface_element==1

          assert(nb_voisins_a_distance_plus_petite!=0 || elem>=nb_elem);
          assert(nb_voisins_a_distance_plus_petite<=nombre_de_voisins_plus_proches_[elem]);
          if (elem<nb_elem) assert(nb_voisins_a_distance_plus_petite==nombre_de_voisins_plus_proches_[elem]);

          for (int elem_loc=0; elem_loc<nb_voisins_a_distance_plus_petite; elem_loc++)
            {
              const int elem_voisin = voisins_a_distance_plus_petite[elem_loc];
              const double surface_dans_elem_voisin = surface_interface_elements_voisins_[elem_voisin];

              surface_totale += surface_dans_elem_voisin;
              surface_par_element[elem_loc] += surface_dans_elem_voisin;

            }//fin du for sur "elem_loc"


          fichier_transport << "L'element " << elem << " est a distance "
                            << distance_interface_elem << " case 1" << finl;
          fichier_transport << "Cet element a perdu le volume " << volume_perdu_[elem] << finl;
          fichier_transport << "Cet element a " << nombre_de_voisins_plus_proches_[elem] << " voisins plus proches" << finl;
          fichier_transport << "La surface d'interface dans l'ensemble de ces voisins les plus proches est : "
                            << surface_totale << " ou " << surface_interface_elements_voisins_[elem]  << finl;
          fichier_transport << "-------------------------------------------------------" << finl;

          //Repartition du volume selon les principes enonces.
          for (int elem_loc=0; elem_loc<nb_voisins_a_distance_plus_petite; elem_loc++)
            {
              const int elem_voisin = voisins_a_distance_plus_petite[elem_loc];


              assert(distance_interface_elem-distance_interface_element_eulerien_[elem_voisin]==1);

              fichier_transport << "La surface d'interface dans le voisin " << elem_voisin
                                << " de " << elem << " est : "
                                << surface_par_element[elem_loc] << " ou "
                                << surface_interface_elements_voisins_[elem_voisin] << finl;

              fichier_transport << "L'element voisin " << elem_voisin << " recoit le volume "
                                << volume_perdu_[elem]*surface_interface_elements_voisins_[elem_voisin]/surface_interface_elements_voisins_[elem]
                                << finl;
              fichier_transport << finl;

            }//fin du for sur "elem_loc"

          fichier_transport << "-------------------------------------------------------" << finl;
          fichier_transport << "*******************************************************" << finl;
          fichier_transport << "*******************************************************" << finl;
          fichier_transport << finl;

          break;


        default :

          assert(nb_voisins_a_distance_plus_petite!=0 || elem>=nb_elem);
          assert(nb_voisins_a_distance_plus_petite==nombre_de_voisins_plus_proches_[elem]);
          if (elem<nb_elem) assert(nb_voisins_a_distance_plus_petite==nombre_de_voisins_plus_proches_[elem]);

          fichier_transport << "L'element " << elem << " est a distance "
                            << distance_interface_elem << " default" << finl;

          fichier_transport << "Cet element a perdu le volume " << volume_perdu_[elem] << finl;
          fichier_transport << "Cet element a " << nombre_de_voisins_plus_proches_[elem] << " voisins plus proches" << finl;
          fichier_transport << "Et distribue a ces voisins les plus proches le volume : "
                            << volume_perdu_[elem]/nombre_de_voisins_plus_proches_[elem] << finl;
          fichier_transport << "*******************************************************" << finl;
          fichier_transport << "*******************************************************" << finl;
          fichier_transport << finl;

          break;

        }//fin du switch

    }//fin du for sur "elem"
}



void  Remailleur_Collision_FT_Thomas::tester_transport_complet(const Maillage_FT_Disc& maillage, DoubleTab& copie_volume_perdu_) const
{
  Nom nom("Transport_complet_proc_");
  nom+=Process::me();
  nom+=".txt";

  SFichier fichier_transport(nom.getChar());

  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,domaine_dis(maillage).valeur());
  const Domaine& domaine = domaine_VF.domaine();

  const int nb_elem = domaine.nb_elem();
  const int nb_elem_tot = domaine.nb_elem_tot();
  const int end_liste = -1;

  //Construction de la liste chainee distance->{ensemble des elements a cette distance}
  //Principe : Creation d'une premiere liste de taille "plus_grande_distance_interface_element_eulerien_"
  //           et initialisee a -1
  //           A la distance "n" sera stocke un element "elem" a distance "n" de l'interface
  //           Creation d'une deuxieme liste de taille "nb_elem_tot" initialisee a -1
  //           A la place "elem" est stocke le numero "elem2" d'un deuxieme element a distance "n"
  //           A la place "elem2" est stocke le numero "elem3" d'un deuxieme element a distance "n" etc
  //           La fin de la liste chainee est reperee par un -1

  IntTab distance_des_elements(plus_grande_distance_interface_element_eulerien_+1);
  IntTab suite_des_elements(nb_elem_tot);
  distance_des_elements = end_liste;
  suite_des_elements = end_liste;

  assert(nb_elem_tot==distance_interface_element_eulerien_.dimension_tot(0));


  fichier_transport << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  fichier_transport << "Nombre d'elements reels : " << nb_elem << finl;
  fichier_transport << "Nombre d'elements virtuels : " << nb_elem_tot << finl;
  fichier_transport << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  fichier_transport << finl;
  fichier_transport << finl;


  fichier_transport << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  double volume_total_avant = 0.;
  for (int elem=0; elem<nb_elem; elem++)
    volume_total_avant += copie_volume_perdu_[elem];

  fichier_transport << "Affichage du volume perdu avant tester_transport_complet() : "
                    << Process::mp_sum(volume_total_avant) << finl;
  fichier_transport << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  fichier_transport << finl;
  fichier_transport << finl;

  //Boucle sur "distance_interface_element_eulerien_" afin de modifier nos listes
  for (int elem=0; elem<nb_elem_tot; elem++)
    {
      const int distance = distance_interface_element_eulerien_[elem];
      if (distance==-1) continue; //on passe directement a l'element suivant
      assert(distance<=plus_grande_distance_interface_element_eulerien_);

      suite_des_elements(elem)=distance_des_elements(distance);
      distance_des_elements(distance)=elem;

    }//fin du for sur "elem"


  //
  //Transport du volume perdu proprement dit
  //

  ArrOfInt voisins_a_distance_plus_petite(100);

  //  int distance_interface_elem_voisin = -2;

  for (int distance=plus_grande_distance_interface_element_eulerien_; distance!=-1; distance--)
    {

      //Parcours de tous les elements "elem" situes a la distance "distance" de l'interface
      //et transport du volume perdu
      for (int elem=distance_des_elements(distance); elem!=-1; elem=suite_des_elements(elem))
        {
          voisins_a_distance_plus_petite.resize_array(0);

          elements_voisins(elem, voisins_a_distance_plus_petite,domaine_VF);
          elements_voisins_a_distance_plus_petite2(elem, voisins_a_distance_plus_petite);
          const int nb_voisins_a_distance_plus_petite = voisins_a_distance_plus_petite.size_array();

          const int distance_interface_elem = distance_interface_element_eulerien_[elem];

          ArrOfDouble surface_par_element(nb_voisins_a_distance_plus_petite);
          double surface_totale = 0.;

          if (copie_volume_perdu_[elem]==0)
            {
              fichier_transport << "L'element " << elem << " n' a pas perdu de volume" << finl;
              fichier_transport << "*******************************************************" << finl;
              fichier_transport << "*******************************************************" << finl;

              continue;
            }; //pas de volume perdu => sortie


          switch(distance_interface_elem)
            {

            case -1 :

              fichier_transport << "L'element " << elem << " est a une distance de l'interface superieure "
                                << "a la distance prescrite par l'utilisateur." << finl;
              fichier_transport << "*******************************************************" << finl;
              fichier_transport << "*******************************************************" << finl;

              break;


            case 0 :

              fichier_transport << "L'element " << elem << " n' est pas traverse par l'interface" << finl;
              fichier_transport << "*******************************************************" << finl;
              fichier_transport << "*******************************************************" << finl;

              break;


            case 1 :

              assert(nb_voisins_a_distance_plus_petite!=0 || elem>=nb_elem);
              assert(nb_voisins_a_distance_plus_petite<=nombre_de_voisins_plus_proches_[elem]);
              if (elem<nb_elem) assert(nb_voisins_a_distance_plus_petite==nombre_de_voisins_plus_proches_[elem]);

              fichier_transport << "Case 1" << finl;
              fichier_transport << "L'element " << elem << " est a distance " << distance_interface_elem
                                << " de l'interface et a perdu le volume : " << copie_volume_perdu_[elem] << finl;
              fichier_transport << "Il a " << nombre_de_voisins_plus_proches_[elem]
                                << " voisins plus proches de l'interface." << finl;
              fichier_transport << "Ces voisins sont : " << finl;


              //Recuperation de la surface totale des facettes traversant
              //l'ensemble des elements a distance 0 de l'interface, ainsi que la surface
              //des facettes traversant chaque element a distance 0 de l'interface
              for (int elem_loc=0; elem_loc<nb_voisins_a_distance_plus_petite; elem_loc++)
                {
                  const int elem_voisin = voisins_a_distance_plus_petite[elem_loc];
                  const double surface_dans_elem_voisin = surface_interface_elements_voisins_[elem_voisin];

                  surface_totale += surface_dans_elem_voisin;
                  surface_par_element[elem_loc] += surface_dans_elem_voisin;
                  assert (std::fabs(surface_dans_elem_voisin-surface_interface_elements_voisins_[elem_voisin])<=1.e-15);

                  fichier_transport << "l'element " << elem_voisin
                                    << " traverse par l'interface de surface : "
                                    << surface_interface_elements_voisins_[elem_voisin] << finl;

                }//fin du for sur "elem_loc"


              fichier_transport << "La surface totale est donc de : " << surface_totale
                                << " ou " << surface_interface_elements_voisins_[elem]
                                << finl;

              assert(std::fabs(surface_totale-surface_interface_elements_voisins_[elem])<=1.e-14);

              //Repartition du volume selon les principes enonces.
              for (int elem_loc=0; elem_loc<nb_voisins_a_distance_plus_petite; elem_loc++)
                {
                  const int elem_voisin = voisins_a_distance_plus_petite[elem_loc];

                  //  assert(distance_interface_elem_voisin!=-1);
                  assert(distance_interface_elem-distance_interface_element_eulerien_[elem_voisin]==1);

                  copie_volume_perdu_[elem_voisin]+=copie_volume_perdu_[elem]
                                                    *surface_interface_elements_voisins_[elem_voisin]/surface_interface_elements_voisins_[elem];

                  fichier_transport << "L'element " << elem_voisin << " recoit de la part de l'element " << elem
                                    << " le volume : "
                                    << copie_volume_perdu_[elem]*surface_interface_elements_voisins_[elem_voisin]/surface_interface_elements_voisins_[elem]
                                    << finl;
                  fichier_transport << "Le nouveau volume perdu par l'element " << elem_voisin
                                    << " est donc : " << copie_volume_perdu_[elem_voisin] << finl;
                }//fin du for sur "elem_loc"

              //remise a zero pour plus de surete
              copie_volume_perdu_[elem]=0.;

              fichier_transport << "Le volume perdu par l'element " << elem
                                << " est maintenant nul : " << copie_volume_perdu_[elem] << finl;

              fichier_transport << "*******************************************************" << finl;
              fichier_transport << "*******************************************************" << finl;

              break;



            default :

              assert(nb_voisins_a_distance_plus_petite!=0 || elem>=nb_elem);
              assert(nb_voisins_a_distance_plus_petite<=nombre_de_voisins_plus_proches_[elem]);
              if (elem<nb_elem) assert(nb_voisins_a_distance_plus_petite==nombre_de_voisins_plus_proches_[elem]);

              fichier_transport << "Default" << finl;
              fichier_transport << "L'element " << elem << " est a distance " << distance_interface_elem
                                << " de l'interface et a perdu le volume : " << copie_volume_perdu_[elem] << finl;
              fichier_transport << "Il a " << nombre_de_voisins_plus_proches_[elem]
                                << " voisins plus proches de l'interface." << finl;
              fichier_transport << "Ces voisins sont : " << finl;

              //Repartition du volume selon les principes enonces.
              for (int elem_loc=0; elem_loc<nombre_de_voisins_plus_proches_[elem]; elem_loc++)
                {
                  const int elem_voisin = voisins_a_distance_plus_petite[elem_loc];

                  //                  assert(distance_interface_elem_voisin!=-1);
                  assert(distance_interface_elem-distance_interface_element_eulerien_[elem_voisin]==1);
                  copie_volume_perdu_[elem_voisin]+=copie_volume_perdu_[elem]/nombre_de_voisins_plus_proches_[elem];

                  fichier_transport << "l'element " << elem_voisin
                                    << " qui recoit le volume perdu : "
                                    <<  copie_volume_perdu_[elem]/nombre_de_voisins_plus_proches_[elem] << finl;
                  fichier_transport << "Le nouveau volume perdu par l'element " << elem_voisin
                                    << " est donc : " << copie_volume_perdu_[elem_voisin] << finl;
                }//fin du for sur "elem_loc"

              //remise a zero pour plus de surete
              copie_volume_perdu_[elem]=0.;

              fichier_transport << "Le volume perdu par l'element " << elem
                                << " est maintenant nul : " << copie_volume_perdu_[elem] << finl;

              fichier_transport << "*******************************************************" << finl;
              fichier_transport << "*******************************************************" << finl;

              break;

            }//fin du switch sur "distance_interface_elem"

        }//fin du for sur "elem"

      copie_volume_perdu_.echange_espace_virtuel();

    }//fin du for sur "distance

  fichier_transport << finl;
  fichier_transport << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  fichier_transport << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;

  double volume_total_perdu = 0.;
  for (int elem=0; elem<nb_elem; elem++)
    volume_total_perdu += copie_volume_perdu_[elem];

  fichier_transport << "Affichage du volume perdu total transporte : " << Process::mp_sum(volume_total_perdu) << finl;
}



void  Remailleur_Collision_FT_Thomas::tester_volume_par_sommet(const Maillage_FT_Disc& maillage, const DoubleTab& copie_volume_perdu_) const
{
  Nom nom("Volume_sommet_proc_");
  nom+=Process::me();
  nom+=".txt";

  SFichier fichier_volume(nom.getChar());

  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,domaine_dis(maillage).valeur());
  const Domaine& domaine = domaine_VF.domaine();

  const Intersections_Elem_Facettes& intersections =
    maillage.intersections_elem_facettes();

  const ArrOfInt& index_elem = intersections.index_elem();
  const IntTab& facettes = maillage.facettes();

  const int nb_elem = domaine.nb_elem();
  const int nb_elem_tot = domaine.nb_elem_tot();
  const int nb_sommets_interface = maillage.nb_sommets();

  ArrOfDoubleFT volume_perdu_par_sommet(nb_sommets_interface);
  volume_perdu_par_sommet = 0.;


  fichier_volume << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  fichier_volume << "Nombre d'elements reels : " << nb_elem << finl;
  fichier_volume << "Nombre d'elements virtuels : " << nb_elem_tot << finl;
  fichier_volume << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  fichier_volume << finl;
  fichier_volume << finl;


  fichier_volume << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  double volume_total_avant = 0.;
  for (int elem=0; elem<nb_elem; elem++)
    volume_total_avant += copie_volume_perdu_[elem];

  fichier_volume << "Affichage du volume perdu avant tester_volume_par_sommet() : "
                 << Process::mp_sum(volume_total_avant) << finl;
  fichier_volume << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  fichier_volume << finl;
  fichier_volume << finl;
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();

  for (int elem=0; elem<nb_elem; elem++)
    {
      assert(copie_volume_perdu_.dimension_tot(0)==nb_elem_tot);

      const double surface_interface_dans_elem = surface_interface_elements_voisins_[elem];

      const double volume_perdu_elem = copie_volume_perdu_[elem];

      fichier_volume << "L'element " << elem << " a perdu le volume " << copie_volume_perdu_[elem] << finl;
      fichier_volume << "La surface d'interface le coupant est de : " << surface_interface_dans_elem << finl;

      //Boucle sur les faces qui traversent l'element:
      int index = index_elem[elem];

      while (index >= 0)
        {
          const Intersections_Elem_Facettes_Data& data =
            intersections.data_intersection(index);

          //Calcul du volume perdu associe a "facette"
          // i.e. volume perdu*(surface de "facette" dans "elem")/(surface totale de l'interface dans "elem")
          if (data.fraction_surface_intersection_!=0)
            {
              const int facette = data.numero_facette_;
              const double coeff_surface = data.fraction_surface_intersection_*surface_facettes[facette]/surface_interface_dans_elem;
              const double volume_perdu_surfacique = volume_perdu_elem*coeff_surface;

              fichier_volume << "La facette " << facette << " coupant l'element " << elem
                             << " a une surface dans l'element de " << data.fraction_surface_intersection_*surface_facettes[facette] << finl;
              fichier_volume << "Cela represente " << coeff_surface << " de la surface totale" << finl;
              fichier_volume << "La facette " << facette << " recoit donc un volume surfacique de "
                             << volume_perdu_surfacique << finl;
              fichier_volume << "-------------------------------------------------------" << finl;

              for (int som_loc=0; som_loc<Objet_U::dimension; som_loc++)
                {
                  //Recuperation d'un sommet de "facette" et des coodonnees barycentrique
                  //du barycentre du centre de gravite de l'intersection entre l'interface
                  //et "elem"
                  const int sommet = facettes(facette,som_loc);
                  assert(sommet<nb_sommets_interface);

                  fichier_volume << "La coordonnee barycentrique du barycentre "
                                 << "de la surface de l'interface dans l'element " << elem
                                 << " par rapport au sommet " << sommet << " est : " <<  data.barycentre_[som_loc]
                                 << finl;
                  fichier_volume << "Le sommet " << sommet << " de la facette " << facette
                                 << " recoit donc le volume : " << volume_perdu_surfacique*data.barycentre_[som_loc];
                  fichier_volume << finl;

                  volume_perdu_par_sommet[sommet]+=volume_perdu_surfacique*data.barycentre_[som_loc];

                }//fin du for sur "som_loc"

              fichier_volume << "-------------------------------------------------------" << finl;

            }//fin du if

          index = data.index_facette_suivante_;

        }//fin du while

      fichier_volume << "*******************************************************" << finl;
      fichier_volume << "*******************************************************" << finl;
      fichier_volume << finl;

    }//fin du for sur "elem"

  fichier_volume << finl;
  fichier_volume << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
  fichier_volume << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;

  //Pour le parallele
  maillage.desc_sommets().collecter_espace_virtuel(volume_perdu_par_sommet, MD_Vector_tools::EV_SOMME);

  double volume_perdu_total = 0.;
  for (int som=0; som<nb_sommets_interface; som++)
    if (!maillage.sommet_virtuel(som)) volume_perdu_total+=volume_perdu_par_sommet[som];

  fichier_volume << "Le volume perdu par sommet total est : " << Process::mp_sum(volume_perdu_total) << finl;
}

// Fonction qui dimensionne les attributs de type IntTab et DoubleTab
int Remailleur_Collision_FT_Thomas::initialiser_data(const Maillage_FT_Disc& maillage)
{
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,domaine_dis(maillage).valeur());
  const Domaine& domaine = domaine_VF.domaine();

  const int nb_elem_tot = domaine.nb_elem_tot();
  const int nb_som_elem = domaine.nb_som_elem();
  const int nb_som_tot = domaine.nb_som_tot();

  //On dimensionne les donnees
  distance_interface_element_eulerien_.reset();
  nombre_de_voisins_plus_proches_.reset();
  surface_interface_elements_voisins_.reset();
  volume_perdu_.reset();
  domaine.creer_tableau_elements(distance_interface_element_eulerien_, RESIZE_OPTIONS::NOCOPY_NOINIT);
  distance_interface_element_eulerien_=-1;
  domaine.creer_tableau_elements(nombre_de_voisins_plus_proches_);
  domaine.creer_tableau_elements(surface_interface_elements_voisins_);
  domaine.creer_tableau_elements(volume_perdu_);

  //On dimensionne les donnees
  voisinage_sommet_.resize(nb_som_tot, RESIZE_OPTIONS::NOCOPY_NOINIT);
  next_elem_.resize(nb_elem_tot*nb_som_elem, RESIZE_OPTIONS::NOCOPY_NOINIT);

  //On leur attribue une valeur
  voisinage_sommet_ = -1;
  construire_voisinage_sommet(maillage) ;

  return (est_dimensionne_=1);
}
