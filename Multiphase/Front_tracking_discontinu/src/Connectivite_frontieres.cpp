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
// File:        Connectivite_frontieres.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////
#include <Connectivite_frontieres.h>
#include <TRUST_Deriv.h>
#include <Ref_Connectivite_frontieres.h>
#include <Zone_VF.h>

Implemente_instanciable(Connectivite_frontieres,"Connectivite_frontieres",Objet_U);
Implemente_ref(Connectivite_frontieres);

Entree& Connectivite_frontieres::readOn(Entree& is)
{
  Cerr << "Connectivite_frontieres::readOn" << finl;
  assert(0);
  exit();
  return is;
}

Sortie& Connectivite_frontieres::printOn(Sortie& is) const
{
  Cerr << "Connectivite_frontieres::printOn" << finl;
  assert(0);
  exit();
  return is;
}

static True_int fonction_tri_lexicographique_3_colonnes(const void *ptr1,
                                                        const void *ptr2)
{
  const int * tab1 = (const int *) ptr1;
  const int * tab2 = (const int *) ptr2;
  int delta;
  delta = tab1[0] - tab2[0];
  if (delta) return delta;
  delta = tab1[1] - tab2[1];
  if (delta) return delta;
  delta = tab1[2] - tab2[2];
  return delta;
}

void Connectivite_frontieres::remplir_def_face_aretes(const Zone_VF& zone_vf)
{
  const Nom& nom_elem = zone_vf.zone().type_elem().valeur().que_suis_je();

  // Definition des aretes si la face est un segment
  static const int segment[2][2] = { { 0, -1 },
    { 1, -1 }
  };
  // ... face est un triangle
  static const int triangle[3][2] = { { 0, 1 },
    { 1, 2 },
    { 2, 0 }
  };
  // .. face est un quadrangle
  static const int quad[4][2] = { { 0, 1 },
    { 1, 3 },
    { 3, 2 },
    { 2, 0 }
  };
  int nb_aretes_par_face = -1;
  const int (*tableau)[2] = 0; // Pointeur sur un tableau de 2 entiers

  if (nom_elem == "Rectangle")
    {
      tableau = segment;
      nb_aretes_par_face = 2;
    }
  else if (nom_elem == "Rectangle_2D_axi")
    {
      tableau = segment;
      nb_aretes_par_face = 2;
    }
  else if (nom_elem == "Triangle")
    {
      tableau = segment;
      nb_aretes_par_face = 2;
    }
  else if (nom_elem == "Tetraedre")
    {
      tableau = triangle;
      nb_aretes_par_face = 3;
    }
  else if (nom_elem == "Hexaedre")
    {
      tableau = quad;
      nb_aretes_par_face = 4;
    }
  else
    {
      Cerr << "Connectivite_frontieres::remplir_def_face_aretes\n";
      Cerr << " Le type d'element " << nom_elem;
      Cerr << " n'est pas reconnu !\n";
      exit();
    }
  def_face_aretes_.resize(nb_aretes_par_face,2);
  for (int i = 0; i < nb_aretes_par_face; i++)
    {
      def_face_aretes_(i,0) = tableau[i][0];
      def_face_aretes_(i,1) = tableau[i][1];
    }
}

void Connectivite_frontieres::associer_zone_vf(const Zone_VF& zone_vf)
{
  refzone_vf_ = zone_vf;

  remplir_def_face_aretes(zone_vf);
  remplir_faces_voisins(zone_vf);
}

/*! @brief Remplissage de faces_voisins_
 *
 */
void Connectivite_frontieres::remplir_faces_voisins(const Zone_VF& zone_vf)
{
  int i_frontiere; // Un numero de frontiere (bord, raccord ou faces internes)
  int i;
  const Zone& zone = zone_vf.zone();
  const int nb_frontieres = zone.nb_front_Cl();
  assert(def_face_aretes_.dimension(0) > 0);
  const int nb_aretes_par_face = def_face_aretes_.dimension(0);
  // Nombre de faces de bord reelles
  const int nb_faces_front = zone.nb_faces_frontiere();

  // Comptage des faces de bord virtuelles
  int nb_faces_frontiere_virt = 0;
  for (i_frontiere = 0; i_frontiere < nb_frontieres; i_frontiere++)
    {
      const ArrOfInt& liste_faces_virt =
        zone.frontiere(i_frontiere).get_faces_virt();
      nb_faces_frontiere_virt += liste_faces_virt.size_array();
    }

  // Nombre total de faces
  const int nb_faces_front_tot = nb_faces_front + nb_faces_frontiere_virt;

  // Remplissage du tableau liste_faces contenant la liste de toutes les
  // faces frontiere (reelles et virtuelles)
  ArrOfInt liste_faces(nb_faces_front_tot);

  // Faces reelles = les nb_faces_frontiere premieres faces de la zone
  // n = nombre de faces remplies dans liste_faces
  int n;
  for (n = 0; n < nb_faces_front; n++)
    liste_faces[n] = n;
  // Maintenant, n==nb_faces_frontiere ...

  // Faces virtuelles :
  {
    const int nb_frontieres_cl = zone.nb_front_Cl();
    //int i_frontiere;
    for (i_frontiere = 0; i_frontiere < nb_frontieres_cl; i_frontiere++)
      {
        const Frontiere& frontiere = zone.frontiere(i_frontiere);
        const ArrOfInt& faces_virtuelles = frontiere.get_faces_virt();
        const int nb_faces_frontiere = faces_virtuelles.size_array();
        for (int ii = 0; ii < nb_faces_frontiere; ii++)
          {
            liste_faces[n] = faces_virtuelles[ii];
            n++;
          }
      }
  }

  // Construction d'une liste d'aretes : pour chaque face frontiere,
  // 2, 3 ou 4 aretes. Pour chaque arete, on construit une ligne du tableau les_aretes
  // colonne 1 et 2 : numeros des deux sommets extremites, le + petit numero en premier
  //                  et en 2D, colonne2 = -1
  // colonne 3 : indice de la face adjacente a l'arete dans zone_vf.face_voisins_
  // colonne 4 : numero de l'arete sur la face

  IntTab les_aretes(nb_faces_front_tot * nb_aretes_par_face, 4);
  int nb_aretes = 0; // Nombre d'aretes rangees dans le tableau

  {
    const IntTab& face_sommets = zone_vf.face_sommets();
    int i_face;        // indice de la face dans liste_faces

    for (i_face = 0; i_face < nb_faces_front_tot; i_face++)
      {

        // Indice de la face a traiter dans zone_vf.face_voisins_ :
        const int face = liste_faces[i_face];

        // Boucle sur les aretes de la face (les aretes de la face sont definie
        //  par le contenu du tableau def_face_aretes)
        for (int i_arete = 0; i_arete < nb_aretes_par_face; i_arete++)
          {
            // Numeros des deux sommets de l'arete sur la face
            int i_sommet0 = def_face_aretes_(i_arete, 0);
            int i_sommet1 = def_face_aretes_(i_arete, 1);
            // Numeros des deux sommets dans la zone
            int sommet0 = face_sommets(face, i_sommet0);
            int sommet1 = (i_sommet1 >= 0) ? face_sommets(face, i_sommet1) : -1;
            // On classe les deux sommets dans l'ordre croissant:
            if (sommet1 >= 0 && sommet1 < sommet0)
              {
                int tmp = sommet0;
                sommet0 = sommet1;
                sommet1 = tmp;
              }
            // Les deux sommets definissent un arete, on la range dans les_aretes :
            les_aretes(nb_aretes, 0) = sommet0;
            les_aretes(nb_aretes, 1) = sommet1;
            les_aretes(nb_aretes, 2) = face;
            les_aretes(nb_aretes, 3) = i_arete;
            nb_aretes++;
          }
      }
  }

  // On trie la liste d'aretes dans l'ordre lexicographique:
  //  par ordre croissant de sommet0, puis de sommet1, puis de i_face
  assert(les_aretes.dimension(1) == 4);
  if (nb_aretes > 0)
    qsort(les_aretes.addr(),
          nb_aretes,
          sizeof(int) * 4,
          fonction_tri_lexicographique_3_colonnes);

  // Maintenant, les aretes identiques se suivent dans le tableau.
  // On trouve donc des couples de numeros de faces adjacentes.
  // On remplit le tableau des faces voisines...

  faces_voisins_.resize(nb_faces_front, nb_aretes_par_face);
  faces_voisins_ = -1;
  for (i = 0; i < nb_aretes - 1; i++)
    {
      // Recherche d'un couple d'aretes identiques
      if (les_aretes(i, 0) == les_aretes(i+1, 0)
          && les_aretes(i, 1) == les_aretes(i+1, 1))
        {
          const int face0 = les_aretes(i, 2);
          const int arete0= les_aretes(i, 3);
          const int face1 = les_aretes(i+1, 2);
          const int arete1= les_aretes(i+1, 3);

          assert(face0 != face1);
          assert(face0 >= nb_faces_front || faces_voisins_(face0, arete0) < 0);
          assert(face1 >= nb_faces_front || faces_voisins_(face1, arete1) < 0);

          if (face0 < nb_faces_front)
            faces_voisins_(face0, arete0) = face1;
          if (face1 < nb_faces_front)
            faces_voisins_(face1, arete1) = face0;
        }
    }

  // Il ne doit pas rester de face de bord reelle sans voisin :
  for (i = 0; i < nb_faces_front; i++)
    {
      for (int j = 0; j < nb_aretes_par_face; j++)
        {
          if (faces_voisins_(i,j) < 0)
            {
              Cerr << "(PE" << me();
              Cerr << ") Erreur dans Connectivite_frontieres::associer_zone_vf\n";
              Cerr << " faces_voisins_(" << i << "," << j << ") < 0" << finl;
              assert(0);
              exit();
            }
        }
    }
}

/*! @brief Constructeur par copie : produit une erreur.
 *
 */
Connectivite_frontieres::Connectivite_frontieres(const Connectivite_frontieres& a): Objet_U(a)
{
  Cerr << "Erreur : Connectivite_frontieres::Connectivite_frontieres(const Connectivite_frontieres &)"
       << finl;
  assert(0);
  exit();
}

/*! @brief Operateur copie : produit une erreur.
 *
 */
const Connectivite_frontieres& Connectivite_frontieres::operator=(const Connectivite_frontieres&)
{
  Cerr << "Erreur : Connectivite_frontieres::operator=" << finl;
  assert(0);
  exit();
  return *this;
}
