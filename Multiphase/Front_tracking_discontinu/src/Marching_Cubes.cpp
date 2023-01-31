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
// File:        Marching_Cubes.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/19
//
//////////////////////////////////////////////////////////////////////////////

#include <Marching_Cubes.h>
#include <TRUST_Deriv.h>
#include <Marching_Cubes_data.h>
#include <ArrOfBit.h>
#include <TRUSTVect.h>
#include <Zone_VF.h>
#include <Zone.h>
#include <Rectangle.h>
#include <Rectangle_2D_axi.h>
#include <Triangle.h>
#include <Tetraedre.h>
#include <Hexaedre.h>
#include <Parser.h>
#include <Comm_Group.h>
#include <communications.h>

Implemente_instanciable_sans_constructeur(Marching_Cubes,"Marching_Cubes",Objet_U);


// **********************************************************************
//                       Methodes publiques

Marching_Cubes::Marching_Cubes() :
  last_def_noeud_size(0)
{
}

/*! @brief Constructeur par copie interdit !
 *
 */
Marching_Cubes::Marching_Cubes(const Marching_Cubes& a): Objet_U(a)
{
  Cerr << "Erreur : Marching_Cubes::Marching_Cubes  Constructeur par copie interdit !" << finl;
  exit();
}

/*! @brief Operateur copie interdit !
 *
 */
const Marching_Cubes& Marching_Cubes::operator=(const Marching_Cubes&)
{
  Cerr << "Erreur : Marching_Cubes::operator= Operateur copie interdit !" << finl;
  exit();
  return *this;
}

Entree& Marching_Cubes::readOn(Entree& is)
{
  Cerr << "Erreur : Marching_Cubes::readOn n'est pas code." << finl;
  exit();
  return is;
}

Sortie& Marching_Cubes::printOn(Sortie& os) const
{
  Cerr << "Erreur : Marching_Cubes::printOn n'est pas code." << finl;
  exit();
  return os;
}

void Marching_Cubes::associer_domaine_vf(const Zone_VF& zone_vf)
{
  const Zone& zone = zone_vf.zone();
  ref_zone_vf_ = zone_vf;
  remplir_data_marching_cubes(zone);
  remplir_renum_virt_loc(zone);
}

/*! @brief Construction d'un maillage en segments ou en triangles comme l'isovaleur d'une fonction discretisee aux sommets du maillage eulerien (associer_domaine_vf).
 *
 *   L'algorithme est un "marching cubes" generalise pour travailler avec un
 *   maillage vf en triangles, rectangles, cubes ou tetraedres.
 *   En cas d'erreur, le maillage est remis a zero (reset()).
 *
 * @param (valeurs_sommets) Un vecteur distribue de valeurs aux sommets du maillage vf dont on va construire l'isovaleur. L'espace virtuel doit etre a jour.
 * @param (isovaleur) L'isovaleur a construire
 * @param (maillage) L'objet maillage a construire. On efface completement son contenu et on le remplace par l'isovaleur.
 * @param (indicatrice_approchee) On y stocke 0 si valeurs_sommets[i] < isovaleur pour tous les sommets de la maille et si valeurs_sommets[i] > isovaleur Certaines mailles (traversees par l'interface) ont une valeur differente qu'on ne calcule pas (on met 0.5). Si phase != AJOUTE_TOUTES_PHASES, on suppose que indicatrice_approchee contient une approximation de l'indicatrice des interfaces existantes avant l'entree dans cette fonction (elle doit contenir 0 ou 1 dans les mailles non traversees par une interface, et une valeur intermediaire dans les mailles traversees). Ce tableau est ecrase quoi qu'il arrive (y compris en cas de valeur de retour=0) (faire une copie du tableau avant d'appeler la fonction si on veut revenir en arriere !)
 * @param (phase) Indique comment mettre a jour indicatrice_approchee. Ce parametre est utilise pour construire une interface en plusieurs etapes (reunion de plusieurs bulles, soustraction d'une phase, ...) phase == AJOUTE_TOUTES_PHASES : la valeur de l'indicatrice approchee est calculee partout phase == 0 : on met a 0 l'indicatrice si valeur<isovaleur, l'indicatrice est inchangee si valeur>isovaleur phase == 1 : on met a 1 l'indicatrice si valeur>isovaleur, l'indicatrice l'indicatrice est inchangee si valeur<isovaleur Valeur de retour: 1: construction reussie. 0: erreur: l'interface qu'on vient de construire entre en collision avec une interface existante (teste base sur la valeur de "indicatrice_approchee" a l'entree de la fonction)
 */
int Marching_Cubes::construire_iso(const DoubleVect& valeurs_sommets,
                                   double isovaleur,
                                   Maillage_FT_Disc& maillage,
                                   DoubleVect& indicatrice_approchee,
                                   const Maillage_FT_Disc::AjoutPhase phase,
                                   int ignorer_collision) const
{
  // Note sur le parallelisme:
  //  La difficulte de l'algorithme consiste a creer un sommet unique aux
  //  sommets des triangles (pour retrouver la connectivite du maillage
  //  qui est ignoree dans la methode de Juric d'origine)
  //  Les sommets du maillage cree se trouvent obligatoirement
  //  sur des aretes du maillage VF, et il y en a au plus un par arete (aretes
  //  sur lesquelles valeurs_sommets-isovaleur change de signe).
  //  Donc un sommet du maillage lagrangien est repere de facon unique par
  //  le couple des sommets de l'arete VF sur laquelle il se trouve
  //  (couples ranges dans le tableau def_noeud).
  //  Dans un premier temps, on cree toutes les facettes de l'isovaleur
  //  et on duplique tous les sommets (on cree 3 sommets par facette).
  //  ("construire_noeuds_et_facettes").
  //  Afin de savoir quels sommets sont sur un bord du domaine, on cree des
  //  sommets supplementaires sur les faces de bord, marques de facon
  //  particuliere ("construire_noeuds_joints").
  //  Puis on va trier tous ces sommets par ordre lexicographique et supprimer
  //  les doublons ("trier_les_noeuds", et "construire_noeuds_uniques")
  //  Or il faut faire le lien logique entre les sommets de joints situes
  //  aux frontieres entre deux processeurs (le noeud i du PE 1 est identique
  //  au noeud j du PE 2) : "correspondance_espaces_distant_virtuel".
  //
  //  On parcourt aussi les elements virtuels et on cree des sommets
  //  supplementaires marques par le numero du processeur. Apres avoir
  //  trie les sommets, on connait pour chaque sommet les processeurs voisins
  //  qui ont cree le meme sommet.
  //
  //  Donc pour que l'algorithme fonctionne, il faut que tous les elements
  //  virtuels voisins d'un element reel par une arete soient connus.
  //  Dans la version actuelle, c'est ok: on connait tous les elements voisins
  //  par au moins un sommet.

  if (! ref_zone_vf_.non_nul())
    {
      Cerr << "Marching_Cubes::construire_iso : Erreur :" << finl;
      Cerr << " Aucune zone n'a ete associee a Marching_Cubes" << finl;
      assert(0);
      exit();
    }

  ArrOfBit signe; // Taille : nombre de sommets euleriens
  // Estimation de la taille necessaire du tableau
  // Jusqu'a "construire_noeuds_uniques" la signification de def_noeud
  // est la suivante:
  //  def_noeud(i,0) = numero du premier sommet du segment qui porte le noeud
  //  def_noeud(i,1) = numero du deuxieme sommet
  //                    (les sommets sont tries: def_noeud(i,1) > def_noeud(i,0))
  //  Pour un sommet interne (construit a partir d'un element reel)
  //   def_noeud(i,2) = me()
  //   def_noeud(i,3) = i
  //   def_noeud(i,4) = numero de l'element qui a servi a construire le noeud
  //  Pour un sommet sur une face de joint avec PE_voisin
  //   def_noeud(i,2) = PE_voisin
  //   def_noeud(i,3) = i
  //   def_noeud(i,4) = -1
  //  Pour un sommet sur une face de bord
  //   def_noeud(i,2) = nproc()
  //   def_noeud(i,3) = i
  //   def_noeud(i,4) = numero de la face de bord dans la zone
  // Un meme sommet (ie, meme segment) peut figurer plusieurs fois comme sommet
  // interne, sommet de bord et/ou sommet de joint.
  IntTab def_noeud;
  def_noeud.set_smart_resize(1);
  def_noeud.resize(last_def_noeud_size * 2, 4, Array_base::NOCOPY_NOINIT);

  maillage.reset();

  IntTab& facettes = maillage.facettes_;
  DoubleTab& coord_noeuds = maillage.sommets_;

  calculer_signe(valeurs_sommets, isovaleur, signe);

  const int resultat_ok = construire_noeuds_et_facettes(signe, def_noeud, facettes,
                                                        indicatrice_approchee, phase);

  if (resultat_ok || ignorer_collision)
    {
      construire_noeuds_joints(signe, def_noeud);

      trier_les_noeuds(def_noeud);

      construire_noeuds_uniques(def_noeud, maillage);
      // Ici def_noeud a change de definition:
      //  def_noeud(i,0) = numero du premier sommet du segment qui porte le noeud
      //  def_noeud(i,1) = numero du deuxieme sommet
      //                    (les sommets sont tries: def_noeud(i,1) > def_noeud(i,0))
      //  def_noeud(i,2) = PE proprietaire du noeud pour les noeuds crees en parcourant
      //                   les elements du domaine.
      //                   Si le noeud est un duplicata cree sur les faces de bord,
      //                   alors def_noeud(i,2) == nb_procs.
      //  def_noeud(i,3) = numero de l'element (-1 si noeud virtuel)
      //  def_noeud(i,4) = numero de la face de bord (-1 si noeud virtuel)

      correspondance_espaces_distant_virtuel(def_noeud, maillage.desc_sommets_);

      calculer_coord_noeuds(valeurs_sommets, isovaleur,
                            def_noeud, maillage);

      {
        int nb_sommets = maillage.sommets_.dimension(0);
        maillage.sommet_PE_owner_. resize_array(nb_sommets);
        maillage.sommet_num_owner_.resize_array(nb_sommets);
        maillage.sommet_elem_.     resize_array(nb_sommets);
        maillage.sommet_face_bord_.resize_array(nb_sommets);
        maillage.drapeaux_sommets_.resize_array(nb_sommets);

        maillage.desc_sommets_.remplir_element_pe(maillage.sommet_PE_owner_);
        // Le descripteur des facettes est vide, mais on calcule quand meme
        // le schema de comm pour qu'il soit valide
        maillage.desc_facettes_.espace_distant().calcul_liste_pe_voisins();
        maillage.desc_facettes_.espace_virtuel().calcul_liste_pe_voisins();
        maillage.desc_facettes_.calcul_schema_comm(maillage.facettes_.dimension(0));

        const int moi = Process::me();

        for (int i = 0; i < nb_sommets; i++)
          {
            maillage.sommet_num_owner_[i] = i;
            const int pe_owner = maillage.sommet_PE_owner_[i];
            const int elem = (pe_owner == moi) ? def_noeud(i, 3) : -1;
            const int face = def_noeud(i, 4);
            maillage.sommet_elem_[i] = elem;
            maillage.sommet_face_bord_[i] = face;
            maillage.drapeaux_sommets_[i] = 0;
          }
        maillage.desc_sommets_.echange_espace_virtuel(maillage.sommet_num_owner_);
        maillage.desc_sommets_.echange_espace_virtuel(maillage.sommet_face_bord_);
      }

      // il faut mettre le statut minimal pour corriger_proprietaire_facettes
      maillage.statut_ = Maillage_FT_Disc::MINIMAL;
      maillage.corriger_proprietaires_facettes();

      // On conserve la taille du tableau temporaire pour la prochaine execution.
      last_def_noeud_size = def_noeud.dimension(0);

      Journal() << "Marching_Cubes::construire_iso" << finl;
      Journal() << " " << facettes.dimension(0) << " facettes, ";
      Journal() << coord_noeuds.dimension(0) << " noeuds, ";

      maillage.maillage_modifie(Maillage_FT_Disc::MINIMAL);
    }
  else
    {
      Journal() << "Marching_Cubes::construire_iso erreur: collision avec une interface existante" << finl;
      maillage.reset();
    }
  return resultat_ok;
}

// Construction d'une interface comme l'isovaleur zero d'une fonction
// dont l'expression est donnee en parametre. L'expression est evaluee
// aux sommets du maillage eulerien et l'interface est construite par
// l'algorithme des marching cubes.
// Le maillage est reinitialise au debut de l'operation.
//
// Parametre: expression
// Signification: une expression mathematique f(x,y,z) comprise par le parser
//                exemple : "x+2*y+z*z"
// Parametre: isovaleur
// Signification: on va construire le maillage de la surface definie par
//                f(x,y,z)=isovaleur
// Parametre: maillage
// Signification: l'objet Maillage dans lequel on va stocker le resultat
// Parametre: indicatrice_approchee
// Signification: voir Marching_Cubes::construire_iso
// Parametre: phase
// Signification: voir Marching_Cubes::construire_iso
// Parametre: eval_expression_sommets
// Signification: un tableau a une dimension de taille nb_sommets (euleriens),
//                et dont les items communs sont correctement initialises
//                pour pouvoir faire un echange_espace_virtuel de valeurs aux
//                sommets.
// Valeur de retour: voir  Marching_Cubes::construire_iso
int Marching_Cubes::construire_iso(const Nom& expression, double isovaleur,
                                   Maillage_FT_Disc& maillage,
                                   DoubleVect& indicatrice_approchee,
                                   const Maillage_FT_Disc::AjoutPhase phase,
                                   DoubleTab& eval_expression_sommets,
                                   int ignorer_collision) const
{
  const int dimension3 = (dimension == 3);

  if (! ref_zone_vf_.non_nul())
    {
      Cerr << "Marching_Cubes::construire_iso : Erreur :" << finl;
      Cerr << " Aucune zone n'a ete associee a Marching_Cubes" << finl;
      assert(0);
      exit();
    }

  std::string expr_chaine(expression);
  Parser parser(expr_chaine, dimension);
  parser.addVar("x");
  parser.addVar("y");
  if (dimension3)
    parser.addVar("z");
  parser.parseString();

  // Construction d'un tableau de valeurs aux sommets euleriens
  const Zone& domaine = ref_zone_vf_.valeur().zone();
  const int nb_sommets = domaine.nb_som();

  for (int i = 0; i < nb_sommets; i++)
    {
      double x, y, z = 0.;
      x = domaine.coord(i, 0);
      y = domaine.coord(i, 1);
      if (dimension3)
        z = domaine.coord(i, 2);
      parser.setVar("x", x);
      parser.setVar("y", y);
      if (dimension3)
        parser.setVar("z", z);
      double valeur = parser.eval();
      eval_expression_sommets(i) = valeur;
    }
  eval_expression_sommets.echange_espace_virtuel();

  // Construction de l'interface
  const int ok = construire_iso(eval_expression_sommets, isovaleur, maillage, indicatrice_approchee, phase,
                                ignorer_collision);
  // l'appel a maillage_modifie() est fait dans construire_iso(valeurs_sommets, ...)
  return ok;
}

// **********************************************************************
//                 Methodes protegees

// Structures de donnees de l'algorithme marching cubes
// (en gros, description des facettes a creer en fonction
//  du type d'element eulerien et du signe de la fonction aux
//  sommets de l'element, voir Marching_Cubes_data.h).
void Marching_Cubes::remplir_data_marching_cubes(const Zone& zone)
{
  // Detection du type d'element eulerien
  const Elem_geom_base& elem_geom = zone.type_elem().valeur();

  const int (*def_aretes)[2]=0;      // Pointeur sur un tableau de 2 entiers
  const int (*def_aretes_faces)[2]=0;
  const int *nb_facettes=0;
  const int *facettes=0;

  int nb_cas_marching_cubes=0;

  // Selection des tableaux a copier...
  // Pour ajouter de nouveaux elements euleriens, il "suffit" d'ajouter
  // un test if ci-dessous et de definir les tableaux mcubes_*
  // dans Marching_Cubes_data.h
  if (sub_type(Rectangle, elem_geom))
    {

      nb_sommets_element = 4;
      nb_aretes_element = 4;
      nb_aretes_faces = 1;
      nb_sommets_par_face = 2;
      nb_sommets_facette = 2;
      nb_cas_marching_cubes = 16;

      def_aretes       = mcubes_def_aretes_vdf_2d;
      def_aretes_faces = mcubes_def_aretes_faces_vdf_2d;
      nb_facettes      = mcubes_nb_facettes_vdf_2d;
      facettes         = mcubes_facettes_vdf_2d;

    }
  else if (sub_type(Triangle, elem_geom))
    {

      nb_sommets_element = 3;
      nb_aretes_element = 3;
      nb_aretes_faces = 1;
      nb_sommets_par_face = 2;
      nb_sommets_facette = 2;
      nb_cas_marching_cubes = 8;

      def_aretes       = mcubes_def_aretes_vef_2d;
      def_aretes_faces = mcubes_def_aretes_faces_vef_2d;
      nb_facettes      = mcubes_nb_facettes_vef_2d;
      facettes         = mcubes_facettes_vef_2d;

    }
  else if (sub_type(Tetraedre, elem_geom))
    {

      nb_sommets_element = 4;
      nb_aretes_element = 6;
      nb_aretes_faces = 3;
      nb_sommets_par_face = 3;
      nb_sommets_facette = 3;
      nb_cas_marching_cubes = 16;

      def_aretes       = mcubes_def_aretes_vef_3d;
      def_aretes_faces = mcubes_def_aretes_faces_vef_3d;
      nb_facettes      = mcubes_nb_facettes_vef_3d;
      facettes         = mcubes_facettes_vef_3d;

    }
  else if (sub_type(Hexaedre, elem_geom))
    {

      nb_sommets_element = 8;
      nb_aretes_element = 12;
      nb_aretes_faces = 4;
      nb_sommets_par_face = 4;
      nb_sommets_facette = 3;
      nb_cas_marching_cubes = 256;

      def_aretes       = mcubes_def_aretes_vdf_3d;
      def_aretes_faces = mcubes_def_aretes_faces_vdf_3d;
      nb_facettes      = mcubes_nb_facettes_vdf_3d;
      facettes         = mcubes_facettes_vdf_3d;

    }
  else
    {
      Cerr << "Erreur dans Marching_Cubes::remplir_data_marching_cubes :" << finl;
      Cerr << " l'element geometrique " << elem_geom.que_suis_je();
      Cerr << " n'est pas pris en charge." << finl;
      exit();
    }

  int i;

  // Creation des tables de donnees
  // Definition des aretes de l'element geometrique eulerien
  mcubes_def_aretes.resize(nb_aretes_element, 2);
  for (i = 0; i < nb_aretes_element; i++)
    {
      mcubes_def_aretes(i, 0) = def_aretes[i][0];
      mcubes_def_aretes(i, 1) = def_aretes[i][1];
    }
  // Definition des aretes d'une face de l'element geometrique
  mcubes_def_aretes_faces.resize(nb_aretes_faces, 2);
  for (i = 0; i < nb_aretes_faces; i++)
    {
      mcubes_def_aretes_faces(i, 0) = def_aretes_faces[i][0];
      mcubes_def_aretes_faces(i, 1) = def_aretes_faces[i][1];
    }

  // Remplissage de l'index des facettes
  mcubes_index_facettes.resize_array(nb_cas_marching_cubes + 1);
  mcubes_nb_facettes.resize_array(nb_cas_marching_cubes);
  int index = 0;
  for (i = 0; i < nb_cas_marching_cubes ; i++)
    {
      mcubes_index_facettes[i] = index;
      int n = nb_facettes[i];
      mcubes_nb_facettes[i] = n;
      index += n * nb_sommets_facette;
    }
  mcubes_index_facettes[i] = index;

  // Remplissage du tableau de description des facettes dans les
  // differents cas marching_cubes
  mcubes_facettes.resize_array(index);
  for (i = 0; i < index; i++)
    mcubes_facettes[i] = facettes[i];
  assert(facettes[i] == -1); // Signature de fin de tableau...
}

True_int fonction_tri_mcubes_renum_virt_loc(const void *pt1, const void *pt2)
{
  int x = *(const int *) pt1;
  int y = *(const int *) pt2;
  return x - y;
}

// Pour chaque joint, on construit un tableau de correspondance :
//  Pour i=0..nombre de sommets de joint,
//   renum_virt_loc(joint)(i,0) = numero distant d'un sommet eulerien de joint
//   renum_virt_loc(joint)(i,1) = numero local du meme sommet
// La difference avec le tableau dans Joint est la suivante:
// Propriete : le tableau est trie par ordre croissant des numeros distants
//             (pour recherche binaire rapide dans renum_sommets_dist_loc).
void Marching_Cubes::remplir_renum_virt_loc(const Zone& zone)
{
  const int nb_joints = zone.nb_joints();
  renum_virt_loc_.dimensionner(nb_joints);
  indice_joint_.resize_array(nproc());
  indice_joint_ = -1;

  for (int num_joint = 0; num_joint < nb_joints; num_joint++)
    {

      // renum_unsorted contient :
      //  renum_unsorted(i,0) = numero de sommet sur le PE voisin
      //  renum_unsorted(i,1) = numero local du sommet

      const Joint& joint = zone.joint(num_joint);
      const IntTab& renum_unsorted = joint.renum_virt_loc();
      const int nb_sommets_joint = renum_unsorted.dimension(0);
      IntTab& renum_sorted = renum_virt_loc_[num_joint];

      // On copie le tableau, puis on le trie par ordre croissant
      // de la premiere colonne.
      renum_sorted = renum_unsorted;
      assert(renum_sorted.dimension_tot(1) == 2);

      qsort(renum_sorted.addr(),
            nb_sommets_joint,
            2 * sizeof(int),
            fonction_tri_mcubes_renum_virt_loc);

      const int PE_voisin = joint.PEvoisin();
      indice_joint_[PE_voisin] = num_joint;
    }
}

// Cette fonction convertit les numeros de sommets euleriens distants stockes
// dans le tableau en numeros locaux. Les numeros doivent correspondre
// a des sommets du joint avec le pe_voisin.
// entree : num_sommet contient des numeros de sommets euleriens du joint
//          avec le pe_voisin. Ce sont les numeros des sommets sur le pe_voisin
// sortie : on remplace les numeros distants par les numeros locaux.
// Comportement indefini si le sommet n'est pas dans le joint.

void Marching_Cubes::renum_sommets_dist_loc(const int pe_voisin,
                                            ArrOfInt& num_sommets) const
{
  // Recherche de l'indice du joint avec pe_voisin
  const int indice_joint = indice_joint_[pe_voisin];
  assert(indice_joint >= 0);
  // Tableau de correspondance numero distant-numero local
  const IntTab& renum_virt_loc = renum_virt_loc_[indice_joint];
  const int dernier_sommet_joint = renum_virt_loc.dimension(0) - 1;
  const int nb_sommets = num_sommets.size_array();

  for(int i = 0; i < nb_sommets; i++)
    {
      const int numero_distant = num_sommets[i];
      // Recherche binaire du sommet dans renum_virt_loc
      int min = 0;
      int max = dernier_sommet_joint;
      int valeur = -1;
      int milieu;
      while (min < max)
        {
          milieu = (min + max) >> 1; // (min + max) / 2
          valeur = renum_virt_loc(milieu, 0);
          if (numero_distant > valeur)
            min = milieu + 1;
          else if (numero_distant < valeur)
            max = milieu - 1;
          else
            min = max = milieu;
        }
      // Teste si le sommet a ete trouve
      assert(renum_virt_loc(min, 0) == numero_distant);
      const int numero_local = renum_virt_loc(min, 1);
      num_sommets[i] = numero_local;
    }
}

// #########################################################################
// Premiere etape : determination du signe de la fonction (valeur-isovaleur)
// en chaque sommet (reel)
// Complexite : N (N=nombre de sommets du maillage eulerien)

void Marching_Cubes::calculer_signe(const DoubleVect& valeurs_sommets,
                                    const double isovaleur,
                                    ArrOfBit& signe) const
{
  int i;

  const int nb_sommets = ref_zone_vf_.valeur().nb_som();
  assert(valeurs_sommets.size() == nb_sommets);
  signe.resize_array(nb_sommets);
  signe = 0;
  // signe vaut 1 si valeur - isovaleur > 0, 0 sinon.
  for(i = 0; i < nb_sommets; i++)
    if (valeurs_sommets[i] - isovaleur > 0.)
      signe.setbit(i);
}

// #########################################################################
// Deuxieme etape : Parcours des elements reels, construction des facettes et
// des noeuds de l'isovaleur. Les noeuds sont situes sur des aretes du maillage
// eulerien et caracterises par les numeros des deux sommets de l'arete.
// Les noeuds sont dupliques (pour chaque arete,
// il y a autant de noeuds identiques que d'elements voisins de l'arete).
// Complexite : N (N=nombre d'elements du maillage eulerien)
// Valeur de retour: 1=>construction reussie
//  0=>Rate:l'interface qui vient d'etre creee entre en collision avec une
//   interface existante. def_noeud et facettes est tout de meme rempli
//   et indicatrice_approchee mis a jour (il faut donc en creer une copie
//   avant d'appeler cette methode pour restaurer la valeur initiale en cas
//   de probleme).
int Marching_Cubes::construire_noeuds_et_facettes(const ArrOfBit& signe,
                                                  IntTab& def_noeud,
                                                  IntTab& facettes,
                                                  DoubleVect& indicatrice_approchee,
                                                  const Maillage_FT_Disc::AjoutPhase phase) const
{
  int resultat_ok = 1; // Valeur de retour de la fonction
  int arete, elem, sommet;

  const Zone& zone = ref_zone_vf_.valeur().zone();
  // Pour chaque element virtuel, numero du PE proprietaire :
  const IntTab& elem_virt_pe_num = zone.elem_virt_pe_num();
  // Raccourci vers les numeros des sommets des elements
  const IntTab& elem_sommets = zone.les_elems();
  const int nb_elements_reels = zone.nb_elem();
  const int nb_elem_tot = zone.nb_elem_tot();
  const int nb_sommets_reels = zone.nb_som();
  // Mon numero de PE
  const int mon_PE = Process::me();

  // Tableaux de travail :
  // Numero attribue au noeud cree sur l'arete i de l'element
  // On initialise a -1 pour assert...
  ArrOfInt numero_noeud_arete(nb_aretes_element);
  numero_noeud_arete = -1;
  // Numero de chaque sommet de l'element en cours de traitement
  ArrOfInt numero_sommet(nb_sommets_element);
  // Signe de la fonction sur chaque sommet de l'element
  ArrOfInt signe_sommet(nb_sommets_element);

  def_noeud.resize(0, 4);
  facettes.resize(0, nb_sommets_facette);
  int nb_facettes = 0;
  int nb_noeuds = 0;
  const int numero_dernier_cas = (1 << nb_sommets_element) - 1;

  // Boucle sur les elements du maillage.
  //  Pour les elements reels, on cree des noeuds sur les aretes et on
  //  cree les faces.
  //  Pour les elements virtuels, on cree seulement les noeuds sur les
  //  aretes reelles, et on enregistre ces noeuds comme appartenant
  //  au processeur voisin.

  for (elem = 0; elem < nb_elem_tot; elem++)
    {

      int cas_marching_cubes = 0;
      int facteur = 1;
      // On recupere les numeros des sommets et le signe en local
      for (sommet = 0; sommet < nb_sommets_element; sommet++)
        {
          const int n = elem_sommets(elem, sommet);
          numero_sommet[sommet] = n;
          int s;
          if (n < nb_sommets_reels)   /* Est-ce un sommet reel ? */
            {
              s = signe[n];
              // Calcul du numero du "cas marching_cubes"
              // C'est une suite de bits 0 ou 1 selon le signe de la fonction
              cas_marching_cubes += s * facteur;
              facteur *= 2;
            }
          else
            {
              // facteur negatif => le cas marching_cubes sera negatif a la fin.
              cas_marching_cubes = -1;
              facteur = -1;
              s = -1;
            }
          signe_sommet[sommet] = s;
        }

      if (elem < nb_elements_reels)
        {
          assert(cas_marching_cubes >= 0);
          double indic;
          if (cas_marching_cubes == 0)
            indic = 0.;
          else if (cas_marching_cubes == numero_dernier_cas)
            indic = 1.;
          else
            indic = 0.5;

          if (phase == Maillage_FT_Disc::AJOUTE_TOUTES_PHASES)
            {
              // On ecrase toutes les interfaces existantes (ce n'est pas un ajout d'une interface)
              indicatrice_approchee[elem] = indic;
            }
          else
            {
              if (indic == (1. - phase))
                {
                  // Ne rien faire
                }
              else if (indicatrice_approchee[elem] != (1. - phase))
                {
                  // Collision !
                  resultat_ok = 0;
                }
              else
                {
                  // Je mets une nouvelle interface
                  indicatrice_approchee[elem] = indic;
                }
            }
        }
      else
        {
          cas_marching_cubes = -1;
        }

      // Early quit : si l'element ne contient pas de facette, c'est a dire
      // si tous les sommets ont le meme signe, on se casse
      if (cas_marching_cubes == 0 || cas_marching_cubes == numero_dernier_cas)
        continue;

      int numero_element_a_stocker = elem;
      int PE_element = mon_PE;
      if (elem >= nb_elements_reels)
        {
          // C'est un element virtuel
          numero_element_a_stocker = -1;
          PE_element = elem_virt_pe_num(elem - nb_elements_reels,
                                        0 /* colonne 0=numero du PE */);
        }

      // Creation des noeuds sur les segments de l'element dont les sommets
      // sont reels et de signe different.

      // (Idee pour ameliorer : precalculer la liste des noeuds a creer
      //  en fonction du cas marching cubes comme pour les facettes)
      // Boucle sur les aretes de l'element, creation des noeuds
      for (arete = 0; arete < nb_aretes_element; arete++)
        {
          numero_noeud_arete[arete] = -1;
          const int s1 = mcubes_def_aretes(arete, 0);
          const int s2 = mcubes_def_aretes(arete, 1);
          const int signe1 = signe_sommet[s1];
          const int signe2 = signe_sommet[s2];
          // Si on est sur une arete virtuelle (un des deux sommets est virtuel)
          // alors on n'ajoute pas de sommet
          // Un noeud est present si la fonction change de signe
          if (signe1 != signe2 && signe1 >= 0 && signe2 >= 0)
            {
              numero_noeud_arete[arete] = nb_noeuds;
              int n1 = numero_sommet[s1];
              int n2 = numero_sommet[s2];
              if (n1 > n2)
                {
                  int n = n1;
                  n1 = n2;
                  n2 = n;
                }
              // Le noeud est caracterise par les numeros des deux sommets
              // de l'arete qui le porte.
              def_noeud.resize(nb_noeuds+1, 5);
              def_noeud(nb_noeuds, 0) = n1;
              def_noeud(nb_noeuds, 1) = n2;
              def_noeud(nb_noeuds, 2) = PE_element;
              def_noeud(nb_noeuds, 3) = nb_noeuds;
              def_noeud(nb_noeuds, 4) = numero_element_a_stocker;
              nb_noeuds++;
            }
        }

      // Creation des facettes si l'element est reel :
      if (elem < nb_elements_reels)
        {
          // Index dans le tableau de description des facettes a creer
          int index = mcubes_index_facettes[cas_marching_cubes];
          const int index_fin = mcubes_index_facettes[cas_marching_cubes+1];
          const int newsize = nb_facettes + mcubes_nb_facettes[cas_marching_cubes];
          facettes.resize(newsize, nb_sommets_facette);

          const int dim = Objet_U::dimension;

          for (; index < index_fin; index += nb_sommets_facette)
            {
              int j;
              for (j = 0; j < dim; j++)
                {
                  const int larete = mcubes_facettes[index+j];
                  const int numero_noeud = numero_noeud_arete[larete];
                  assert(numero_noeud >= 0);
                  facettes(nb_facettes, j) = numero_noeud;
                }
              nb_facettes++;
            }
        }
    }
  // S'il y a une erreur sur un processeur, tout le monde renvoie 0
  resultat_ok = ::mp_min(resultat_ok);
  return resultat_ok;
}

/*! @brief Ajout des sommets situes sur des faces (bords ou joints) dans le tableau def_noeud.
 *
 *   Soit faces_sommets est une liste de faces d'un joint, dans ce cas, numero_PE
 *   est le PE_voisin du joint.
 *   Soit faces_sommets est la liste des faces de la zone et nb_faces_a_traiter
 *   est le nombre de faces de bord. Dans ce cas, numero_PE = nproc().
 *
 */
void Marching_Cubes::construire_noeuds_liste_faces(const ArrOfBit& signe,
                                                   const IntTab& faces_sommets,
                                                   const int nb_faces_a_traiter,
                                                   const int numero_PE,
                                                   IntTab& def_noeud) const
{
  // ATTENTION : static : il faut resizer explicitement au cas ou la
  // classe est utilisee avec des discretisations differentes.
  // Numero attribue au noeud cree sur l'arete i de la face
  // On initialise a -1 pour assert...
  static ArrOfInt numero_noeud_arete;
  numero_noeud_arete.resize_array(nb_aretes_element);
  numero_noeud_arete = -1;
  // Numero de chaque sommet de la face
  static ArrOfInt numero_sommet;
  numero_sommet.resize_array(nb_sommets_par_face);
  // Signe de la fonction sur chaque sommet de la face
  static ArrOfInt signe_sommet;
  signe_sommet.resize_array(nb_sommets_par_face);

  const int nb_procs = Process::nproc();

  for (int face = 0; face < nb_faces_a_traiter; face++)
    {

      int somme_signe = 0;
      for (int sommet = 0; sommet < nb_sommets_par_face; sommet++)
        {
          const int n = faces_sommets(face, sommet);
          const int s = signe[n];
          numero_sommet[sommet] = n;
          signe_sommet[sommet] = s;
          somme_signe += s;
        }

      // Si tous les sommets ont le meme signe, on arrete pour cette face
      if (somme_signe == 0 || somme_signe == nb_sommets_par_face)
        continue;

      // Boucle sur les aretes, creation des noeuds

      for (int arete = 0; arete < nb_aretes_faces; arete++)
        {
          const int s1 = mcubes_def_aretes_faces(arete, 0);
          const int s2 = mcubes_def_aretes_faces(arete, 1);
          const int signe1 = signe_sommet[s1];
          const int signe2 = signe_sommet[s2];
          if (signe1 != signe2)
            {
              const int nb_noeuds = def_noeud.dimension(0);
              numero_noeud_arete[arete] = nb_noeuds;
              int n1 = numero_sommet[s1];
              int n2 = numero_sommet[s2];
              if (n1 > n2)
                {
                  int n = n1;
                  n1 = n2;
                  n2 = n;
                }
              // Le noeud est caracterise par les numeros des deux sommets
              // de l'arete qui le porte, le premier numero etant tjs le + petit.
              def_noeud.resize(nb_noeuds+1, 5);
              def_noeud(nb_noeuds, 0) = n1;
              def_noeud(nb_noeuds, 1) = n2;
              def_noeud(nb_noeuds, 2) = numero_PE;
              def_noeud(nb_noeuds, 3) = nb_noeuds;
              if (numero_PE == nb_procs)
                // On balaie les faces de bord
                def_noeud(nb_noeuds, 4) = face;
              else
                // On balaie les faces de joint
                def_noeud(nb_noeuds, 4) = -1;
            }
        } // Fin de la boucle sur les aretes
    } // Fin de la boucle sur les faces
}

void Marching_Cubes::construire_noeuds_joints(const ArrOfBit& signe,
                                              IntTab& def_noeud) const
{
  const Zone_VF& zone_vf = ref_zone_vf_.valeur();
  // Creation des sommets situes sur les faces de bord
  {
    // Les faces de bord sont les premieres faces de la zone

    const IntTab& faces_sommets = zone_vf.face_sommets();
    const int numero_PE = Process::nproc();
    const int nb_faces_bord = zone_vf.nb_faces_bord();

    construire_noeuds_liste_faces(signe,
                                  faces_sommets,
                                  nb_faces_bord,
                                  numero_PE,
                                  def_noeud);
  }
}

// #########################################################################
// Quatrieme etape
// Tri lexicographique des noeuds crees. Le critere est le suivant:
// - ordre croissant du numero du premier sommet,
// - si 1er sommet identique, ordre croissant du deuxieme sommet,
// - si 2eme sommet identique, ordre croissant du PE associe
// Complexite : N * log(N) (N=nombre de noeuds dupliques du maillage lagrangien)

// Fonction de comparaison lexicographique de deux noeuds:
// Tri par ordre croissant de premiere colonne,
// puis deuxieme, puis troisieme colonne.

True_int fct_compar_mcubes_noeuds(const void *pt1, const void *pt2)
{
  const int * const p1 = (const int *) pt1;
  const int * const p2 = (const int *) pt2;
  int x;
  x = p1[0] - p2[0];
  if (x) return x;
  x = p1[1] - p2[1];
  if (x) return x;
  x = p1[2] - p2[2];
  return x;
}

void Marching_Cubes::trier_les_noeuds(IntTab& def_noeud) const
{
  int nb_lignes = def_noeud.dimension(0);
  int nb_colonnes = def_noeud.dimension(1);

  qsort(def_noeud.addr(),             // tableau a trier
        nb_lignes,                    // nombre d'elements du tableau
        nb_colonnes * sizeof(int), // taille de chaque element
        fct_compar_mcubes_noeuds);    // fonction de comparaison
}

// #########################################################################
// Cinquieme etape :
// Reduction de la liste des noeuds a une liste ou chaque noeud est unique,
// creation de la structure de descripteur noeuds distants/virtuels,
// et mise a jour des facettes avec les numeros definitifs des noeuds.
// A la fin, def_noeud contient une liste unique d'aretes.

// On suppose que la liste de noeuds est triee (voir ci-dessus).
// Ainsi, les noeuds identiques se suivent et le PE proprietaire apparait
// en premier.
// Les noeuds distants et virtuels sont crees dans l'ordre
// croissant du numero local du noeud.

void Marching_Cubes::construire_noeuds_uniques(IntTab& def_noeud,
                                               Maillage_FT_Disc& maillage) const
{
  const int mon_PE   = Process::me();
  const int nb_procs = Process::nproc();
  const Zone_VF& zone_vf = ref_zone_vf_.valeur();
  const IntTab& face_voisins = zone_vf.face_voisins();
  const int nb_noeuds = def_noeud.dimension(0);
  ArrOfInt renumerotation(nb_noeuds);
  IntTab& facettes = maillage.facettes_;
  Descripteur_FT& espace_distant = maillage.desc_sommets_.espace_distant();
  Descripteur_FT& espace_virtuel = maillage.desc_sommets_.espace_virtuel();

  // L'initialisation ne sert que pour le test assert car
  // normalement le tableau est entierement rempli.
  renumerotation = -1;

  // Indice du dernier noeud unique ajoute au tableau
  int dernier = -1;
  // Numeros des sommets du dernier noeud ajoute
  int dernier_som0 = -1;
  int dernier_som1 = -1;
  int dernier_PEvoisin_ajoute = -1;

  int noeud;
  enum {NOEUD_REEL, NOEUD_VIRTUEL} type_noeud = NOEUD_REEL;

  // Certains sommets sont crees sur les elements virtuels alors qu'ils n'existent
  // pas sur les elements reels (element virtuel qui possede une arete virtuelle
  // dont les deux sommets sont reels). Il faut supprimer ces sommets car le processeur
  // voisin ne trouvera pas de noeud correspondant dans son espace virtuel.
  // On met donc keep_last_node a 1 des qu'on trouve me() dans le numero du processeur.
  int keep_last_node = 1;

  for (noeud = 0; noeud < nb_noeuds; noeud++)
    {
      const int noeud_som0 = def_noeud(noeud, 0);
      const int noeud_som1 = def_noeud(noeud, 1);
      const int noeud_PE = def_noeud(noeud, 2);
      const int noeud_ancien_numero = def_noeud(noeud, 3);
      const int noeud_ligne_contact = (noeud_PE == nb_procs);
      const int noeud_elem_ou_face = def_noeud(noeud, 4);
      // Numero de la face ou se trouve le noeud
      const int noeud_face = noeud_ligne_contact ? noeud_elem_ou_face : -1;
      // Numero de l'element ou se trouve le noeud
      int noeud_element;
      if (noeud_ligne_contact)
        {
          // Le noeud est sur une face de bord => un des deux voisins est -1.
          assert(face_voisins(noeud_face, 0) == -1
                 || face_voisins(noeud_face, 1) == -1);
          noeud_element =
            face_voisins(noeud_face, 0) + face_voisins(noeud_face, 1) + 1;
        }
      else
        {
          noeud_element = noeud_elem_ou_face;
        }

      // Le noeud est-il different du precedent ?
      if (noeud_som0 != dernier_som0 || noeud_som1 != dernier_som1)
        {
          // Oui, on cree un nouveau noeud. Si le noeud precedent doit etre
          // ignore, on le supprime:
          if (keep_last_node)
            {
              if (dernier >= 0 && type_noeud == NOEUD_VIRTUEL)
                {
                  // Quel est le proprietaire du noeud ?
                  const int pe_owner = def_noeud(dernier, 2);
                  // On ajoute le noeud a l'espace virtuel
                  espace_virtuel.ajoute_element(pe_owner, dernier);
                }
              dernier++;
              keep_last_node = 0;
            }
          // On ajoute le nouveau noeud a la liste des noeuds uniques
          def_noeud(dernier, 0) = dernier_som0 = noeud_som0;
          def_noeud(dernier, 1) = dernier_som1 = noeud_som1;
          def_noeud(dernier, 2) = noeud_PE;
          def_noeud(dernier, 3) = noeud_element;
          def_noeud(dernier, 4) = noeud_face;

          if (noeud_PE != mon_PE)
            {
              // Le plus petit PE qui possede le noeud n'est pas moi.
              // Donc le noeud est virtuel.
              type_noeud = NOEUD_VIRTUEL;
            }
          else
            {
              // Le plus petit PE qui possede le noeud est moi.
              // Le noeud est reel (ou distant si on trouve d'autres PEs
              // par la suite).
              type_noeud = NOEUD_REEL;
              dernier_PEvoisin_ajoute = noeud_PE;
            }
        }
      else
        {
          // C'est le meme noeud que le precedent.
          // On cherche si le noeud est sur des faces de bord et/ou
          // sur des faces de joint.
          if (type_noeud == NOEUD_REEL)
            {
              if (noeud_PE == nb_procs)
                {
                  // Le noeud est sur une face de bord. On connait maintenant
                  // la face de bord et on met le noeud dansl'element adjacent
                  def_noeud(dernier, 3) = noeud_element;
                  def_noeud(dernier, 4) = noeud_face;
                }
              else if (noeud_PE != dernier_PEvoisin_ajoute)
                {
                  // Le noeud appartient a une face de joint avec noeud_PE :
                  // il y a un nouveau voisin a ajouter a l'espace distant
                  espace_distant.ajoute_element(noeud_PE, dernier);
                  dernier_PEvoisin_ajoute = noeud_PE;
                }
            }
        }
      // Si le noeud existe sur le processeur local, on peut le garder :
      if (noeud_PE == mon_PE)
        keep_last_node = 1;

      renumerotation[noeud_ancien_numero] = dernier;
    }
  if (keep_last_node)
    {
      if (dernier >= 0 && type_noeud == NOEUD_VIRTUEL)
        {
          // Quel est le proprietaire du noeud ?
          const int pe_owner = def_noeud(dernier, 2);
          // On ajoute le noeud a l'espace virtuel
          espace_virtuel.ajoute_element(pe_owner, dernier);
        }
      dernier++;
    }

  espace_distant.calcul_liste_pe_voisins();
  espace_virtuel.calcul_liste_pe_voisins();
  maillage.desc_sommets_.calcul_schema_comm(dernier);

  def_noeud.resize(dernier, 5);

  // #########################################################################
  // Sixieme etape :
  // On remplace les numeros des noeuds des facettes par les nouveaux
  // numeros uniques.
  // Complexite : N (N=nombre de facettes du maillage lagrangien)

  {
    ArrOfInt& tab = facettes;  // tableau vu comme unidimensionnel
    int i;
    // Nb total de numeros de sommets = nb de facettes * nb sommets par facette
    const int nb_sommets = tab.size_array();
    for (i = 0; i < nb_sommets; i++)
      {
        int num_sommet = tab[i];
        num_sommet = renumerotation[num_sommet];
        assert(num_sommet > -1);
        tab[i] = num_sommet;
      }
  }
}

// #########################################################################
// Septieme etape :
// Tri des listes d'elements virtuels de sorte que les sommets partages
// apparaissent dans le meme ordre dans la liste des elements distants d'un PE
// et dans la liste des elements virtuels du PE voisin.
// Les noeuds virtuels sont toujours sur des aretes de joint.
// A l'entree dans cette fonction, le descripteur contient les bons elements pour
// chaque PE voisins (distant et virtuel), mais ils sont a priori dans le desordre.
// Exemple:
// * le processeur 1 a les noeuds 0:A, 1:B, 2:C, 3:D, 4:E, 5:F, 6:G.
// * le processeur 2 a les noeuds 0:a, 1:b, 2:c, 3:d, 4:e, 5:f.
// * les noeuds communs sont A=d, C=b, F=c, G=f
// * au depart:
//     le tableau des elements distants de 1 est (0,2,5,6) (indices des noeuds communs)
//     le tableau des elements virtuels de 2 est (1,2,3,5) (indices des noeuds communs)
// * Quand le processeur 1 envoie son espace distant, il envoie (A,C,F,G) dans cet
//   ordre. Le but est que le processeur 2 associe les donnees recues a (d,b,c,f)
//   dans cet ordre. Le tableau des elements virtuels du processeur 2 doit etre
//   (3,1,2,5). On doit reordonner les elements virtuels.
//
// Pour faire ca, 2 recoit les noeuds de 1 (une cle qui decrit le noeud de
// facon unique). Il recoit (A,C,F,G). Il traduit cette cle sous la forme locale:
// (d,b,c,f). Puis il cherche les numeros des noeuds correspondant : (3,1,2,5),
// c'est son espace virtuel.

void Marching_Cubes::correspondance_espaces_distant_virtuel(const IntTab& def_noeud,
                                                            Desc_Structure_FT& desc) const
{
  // Pour chaque espace distant, on construit un tableau contenant une cle qui decrit
  // les noeuds d'interface de facon unique: la cle est constituee des numeros
  // des deux sommets euleriens de l'arete qui porte le noeud.
  // On aurait aussi pu utiliser les coordonnees x,y,z comme cle mais c'est moins sur...
  // En pratique, il s'agit simplement des deux premieres colonnes du tableau def_noeud.

  // Copie du tableau
  IntTab cles(def_noeud);

  // On echange ce tableau : pour chaque noeud virtuel, on recupere la cle distante.
  desc.echange_espace_virtuel(cles);

  // Boucle sur les PEs voisins de l'espace virtuel
  Descripteur_FT& esp_virtuel = desc.espace_virtuel();
  const ArrOfInt& pe_voisins = esp_virtuel.pe_voisins();
  const int nb_pe_voisins = pe_voisins.size_array();
  int index_pe;
  IntTabFT cles_pe;    // Les cles des noeuds virtuels d'un voisin particulier
  IntTabFT temp_renum; // Tableau temporaire pour l'envoi a la fonction de renumerotation
  // Elements de l'espace virtuel tries dans le bon ordre
  ArrOfIntFT nouvel_espace_virtuel;

  for(index_pe = 0; index_pe < nb_pe_voisins; index_pe++)
    {

      const int pe_voisin = pe_voisins[index_pe];
      const ArrOfInt& elements = esp_virtuel.elements(pe_voisin);
      const int nb_elements = elements.size_array();

      cles_pe.resize(nb_elements, 3);
      temp_renum.resize(nb_elements, 2);
      nouvel_espace_virtuel.resize_array(nb_elements);

      // On copie les cles, on ajoute un numero d'ordre et on remplace le numero
      // des sommets euleriens de l'arete (qui est actuellement un numero sur
      // le processeur distant) par le numero local.
      // Copie des cles de l'espace virtuel
      for (int i = 0; i < nb_elements; i++)
        {
          const int n = elements[i];
          temp_renum(i, 0) = cles(n, 0);
          temp_renum(i, 1) = cles(n, 1);
        }
      // Transformation des numeros de sommets euleriens distants de chaque cle
      // en numeros de sommets locaux.
      renum_sommets_dist_loc(pe_voisin, temp_renum);
      // Creation d'un tableau de cles avec une colonne supplementaire contenant
      // le numero d'ordre de la cle.
      for (int i = 0; i < nb_elements; i++)
        {
          // On corrige la cle : sommet de plus petit numero en premier
          int s0 = temp_renum(i, 0);
          int s1 = temp_renum(i, 1);
          if (s0 < s1)
            {
              cles_pe(i, 0) = s0;
              cles_pe(i, 1) = s1;
            }
          else
            {
              cles_pe(i, 0) = s1;
              cles_pe(i, 1) = s0;
            }
          cles_pe(i, 2) = i;
        }
      // Par exemple, on aurait : cles_pe = { (d,0), (b,1), (c,2), (f,3) }

      // Maintenant, on trie ces cles par ordre lexicographique (meme tri que
      // le tableau def_noeud)
      qsort(cles_pe.addr(), nb_elements, 3 * sizeof(int),
            fct_compar_mcubes_noeuds);
      // Les cles apparaissent maintenant dans le meme ordre que dans le
      // tableau def_noeud, donc dans le meme ordre que dans l'espace virtuel, soit:
      //   (exemple : cles_pe = { (b,1), (c,2), (d,0), (f,3) })
      assert(cles_pe.dimension(0) == elements.size_array());
      for (int i = 0; i < nb_elements; i++)
        {
          const int element = elements[i];
          // Si les espaces distants et virtuels sont coherents, les cles
          // sont identiques.
          assert(cles_pe(i, 0) == def_noeud(element, 0));
          assert(cles_pe(i, 1) == def_noeud(element, 1));
          // Mettre ce noeud au bon endroit dans l'espace virtuel :
          // le noeud b en 1, le noeud c en 2, le noeud d en 0, le noeud f en 3.
          int position = cles_pe(i, 2);
          nouvel_espace_virtuel[position] = element;
        }
      // On remplace l'espace virtuel par le nouvel espace trie.
      esp_virtuel.set_elements(pe_voisin, nouvel_espace_virtuel);
    }

  // Les espaces virtuels et distants sont maintenant coherents, on
  // peut calculer le schema de comm :
  desc.espace_virtuel().calcul_liste_pe_voisins();
  desc.calcul_schema_comm(def_noeud.dimension(0));
}

// #########################################################################
// Huitieme etape :
// Calcul des coordonnees des noeuds sur les aretes.
// Complexite : N (N=nombre de noeuds uniques du maillage lagrangien)

void Marching_Cubes::calculer_coord_noeuds(const DoubleVect& valeurs_sommets,
                                           const double isovaleur,
                                           const IntTab& def_noeud,
                                           Maillage_FT_Disc& maillage) const
{
  // Raccourci vers les coordonnees des sommets du maillage eulerien
  const DoubleTab& coord_ = ref_zone_vf_.valeur().zone().coord_sommets();
  const int nb_noeuds = def_noeud.dimension(0);
  DoubleTab& coord_noeuds = maillage.sommets_;
  coord_noeuds.resize(nb_noeuds, dimension);
  int noeud;
  int dimension3 = (dimension == 3);

  for (noeud = 0; noeud < nb_noeuds; noeud++)
    {

      int s0 = def_noeud(noeud, 0);
      int s1 = def_noeud(noeud, 1);
      double valeur0 = valeurs_sommets(s0);
      double valeur1 = valeurs_sommets(s1);
      // valeur1-valeur0 est forcement non nul, sinon les signes seraient les memes
      double alpha = (isovaleur - valeur0) / (valeur1 - valeur0);
      // Sanity checking
      if (alpha < 0.) alpha = 0.;
      if (alpha > 1.) alpha = 1.;
      // Coordonnees (on deroule, c'est plus efficace ...?)
      coord_noeuds(noeud, 0) = coord_(s0,0)*(1.-alpha)+coord_(s1,0)*alpha;
      coord_noeuds(noeud, 1) = coord_(s0,1)*(1.-alpha)+coord_(s1,1)*alpha;
      if (dimension3)
        coord_noeuds(noeud, 2) = coord_(s0,2)*(1.-alpha)+coord_(s1,2)*alpha;
    }

  const Desc_Structure_FT& desc = maillage.desc_sommets_;
  desc.echange_espace_virtuel(coord_noeuds);
}

