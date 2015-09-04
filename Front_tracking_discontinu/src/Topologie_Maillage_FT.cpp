/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        Topologie_Maillage_FT.cpp
// Directory:   $TRUST_ROOT/src/Front_tracking_discontinu
// Version:     /main/16
//
//////////////////////////////////////////////////////////////////////////////

#include <Topologie_Maillage_FT.h>
#include <Deriv_Topologie_Maillage_FT.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Motcle.h>
#include <Remailleur_Collision_FT_Juric.h>
#include <time.h>
#include <SFichier.h>
#include <Zone_VF.h>
#include <Connex_components.h>
#include <Array_tools.h>
#include <Param.h>
#include <Check_espace_virtuel.h>
#include <communications.h>

Implemente_instanciable_sans_constructeur(Topologie_Maillage_FT,"Topologie_Maillage_FT",Objet_U);

Implemente_deriv(Topologie_Maillage_FT);

Topologie_Maillage_FT::Topologie_Maillage_FT() :
  active_(0),
  Erreur_max_coordonnees_(-1.),
  juric_local_(0),
  phase_continue_(0) // B.M. j'aurais bien mis -1 mais il faut rendre coherent avec phase_marquee de Transport_Marqueur_FT
{
}

Sortie& Topologie_Maillage_FT::printOn(Sortie& os) const
{
  Cerr << "Topologie_Maillage_FT::printOn" << finl;
  exit();
  return os;
}

Entree& Topologie_Maillage_FT::readOn(Entree& is)
{
  int juric_pour_tout=0;
  int phase_continue;
  Param param(que_suis_je());
  param.ajouter_flag("active", &active_);
  param.ajouter_deriv("type_remaillage", "Remailleur_Collision_FT_", &remailleur_Collision_, Param::REQUIRED);
  param.ajouter_flag("juric_pour_tout", &juric_pour_tout);
  param.ajouter_flag("juric_local", &juric_local_);
  param.ajouter("phase_continue", &phase_continue);
  param.dictionnaire("0", 0);
  param.dictionnaire("1", 1);
  param.lire_avec_accolades_depuis(is);
  phase_continue_ = phase_continue;

  if (juric_pour_tout)
    {
      Cerr<<"Warning : juric_pour_tout is obsolete "<<finl;
    }
  // Verifications qu'on a bien tout lu :
  if (!active_)
    {
      if (Process::je_suis_maitre())
        Cerr << "ATTENTION : vous n'avez pas active la gestion des ruptures-coalescences..." << finl;
    }
  else
    {
      if (! remailleur_Collision_.non_nul())
        {
          Cerr << "Erreur dans Topologie_Maillage_FT::readOn :\n"
               << " Il faut fournir le type du remailleur et ses parametres.\n"
               << " Valeur conseillee: Juric { ... }" << finl;
          exit();
        }
      if (!sub_type(Remailleur_Collision_FT_Juric, remailleur_Collision_.valeur()))
        {
          Cerr << "Erreur dans Topologie_Maillage_FT::readOn :\n"
               << " Il faut un remailleur de type Juric" << finl;
          exit();
        }
    }
  return is;
}

int Topologie_Maillage_FT::get_phase_continue() const
{
  if (phase_continue_ < 0.)
    {
      Cerr << "Error in Topologie_Maillage_FT::get_phase_continue(): phase_continue must be specified in the data file !" << finl;
      exit();
    }
  return phase_continue_;
}

// Description: tri du tableau tab dans l'ordre croissant.
static inline void trier_trois_entiers(int tab[3])
{
  int tmp;
  int a = tab[0];
  int b = tab[1];
  int c = tab[2];
  if (a > b)
    {
      tmp = a;
      a = b;
      b = tmp;
    }
  if (b > c)
    {
      tmp = b;
      b = c;
      c = tmp;
    }
  if (a > b)
    {
      tmp = a;
      a = b;
      b = tmp;
    }
  tab[0] = a;
  tab[1] = b;
  tab[2] = c;
}

// Description:
//  Copie les coordonnees des trois sommets d'indices facette[i] dans
//  la matrice coord.
static inline void get_coord_sommets_triangle(const int facette[3],
                                              const DoubleTab& sommets,
                                              double coord[3][3])
{
  int i;
  for (i = 0; i < 3; i++)
    {
      const int som = facette[i];
      coord[i][0] = sommets(som,0);
      coord[i][1] = sommets(som,1);
      coord[i][2] = sommets(som,2);
    }
}

static inline void calcul_equation_plan(const double *coord,
                                        const DoubleTab& normales,
                                        const int facette,
                                        double plan[4])
{
  double a = normales(facette, 0);
  double b = normales(facette, 1);
  double c = normales(facette, 2);
  plan[0] = a;
  plan[1] = b;
  plan[2] = c;
  plan[3] = - (a * coord[0] + b * coord[1] + c * coord[2]);
}

static inline void calcul_distance_plan(const double coord[3][3],
                                        const double plan[4],
                                        double distance[3])
{
  double a = plan[0];
  double b = plan[1];
  double c = plan[2];
  double d = plan[3];
  int i;
  for (i = 0; i < 3; i++)
    distance[i] = a * coord[i][0] + b * coord[i][1] + c * coord[i][2] + d;
}

// Description: Calcul de la distance s_min et s_max par rapport au plan
//  "2" des extremites du segment d'intersection entre un triangle est un plan "1".
// distance_1 est la distance entre les sommets du triangle et le plan "1".
//  On suppose que les trois valeurs ne sont pas toutes de meme signe
//  (donc l'intersection est non vide).
// distance_2 est la distance entre les sommets du triangle et le plan "2".
// On renvoie s_min et s_max avec s_min < s_max.
static inline void calcul_intersection_plan_triangle(const double distance_1[3],
                                                     const double distance_2[3],
                                                     double& s_min,
                                                     double& s_max)
{
  int n = 0;
  double s[2];
  int i;
  double z0 = distance_1[2];
  double d0 = distance_2[2];
  for (i = 0; i < 3; i++)
    {
      double z1 = distance_1[i];
      double d1 = distance_2[i];
      if (((z0 > 0.) + (z1 > 0.)) == 1)
        {
          // Changement de signe entre le sommet i et le sommet j
          // Le denominateur ne peut pas etre nul :
          const double inv_delta = 1. / (z1 - z0);
          // Calcul de la distance au plan orthogonal aux deux facettes
          s[n] = (d0 * z1 - d1 * z0) * inv_delta;
          n++;
        }
      z0 = z1;
      d0 = d1;
    }
  assert(n == 2);
  if (s[0] < s[1])
    {
      s_min = s[0];
      s_max = s[1];
    }
  else
    {
      s_min = s[1];
      s_max = s[0];
    }
}

// Recodage B. Mathieu


// Calcul du determinant
//  det( (p2 - p1) , (p3 - p1) )
// pour p1, p2 et p3 des vecteurs a deux composantes
static inline double prod_vect_triangle(double *p1, double *p2, double *p3)
{
  return (p2[0]-p1[0])*(p3[1]-p1[1]) - (p2[1]-p1[1])*(p3[0]-p1[0]);
}

// Description:
//  Voir test_intersection_facettes_3D. Algorithme B. Mathieu
//  Renvoie 0 si les facettes ne se coupent pas
//         -1 si 1 sommet commun
//         -2 si 2 sommets communs
//          1 si les facettes se coupent
int Topologie_Maillage_FT::test_intersection_facettes_2D(
  int fa70, int fa71,
  const Maillage_FT_Disc& maillage) const
{
  const IntTab& facettes = maillage.facettes();
  // Sommet commun ?
  int facette[2][2];
  facette[0][0] = facettes(fa70, 0);
  facette[0][1] = facettes(fa70, 1);
  facette[1][0] = facettes(fa71, 0);
  facette[1][1] = facettes(fa71, 1);
  int nb_sommets_communs = 0;
  if (facette[0][0] == facette[1][0] || facette[0][0] == facette[1][1])
    nb_sommets_communs++;
  if (facette[0][1] == facette[1][0] || facette[0][1] == facette[1][1])
    nb_sommets_communs++;

  if (nb_sommets_communs)
    return -nb_sommets_communs;

  const DoubleTab& sommets = maillage.sommets();
  double c[2][2][2];
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      {
        int s = facette[i][j];
        c[i][j][0] = sommets(s, 0);
        c[i][j][1] = sommets(s, 1);
      }

  const double v00 = prod_vect_triangle(c[0][1],c[0][0],c[1][0]);
  const double v01 = prod_vect_triangle(c[0][1],c[0][0],c[1][1]);
  if (v00 * v01 > 0.)
    return 0;
  const double v10 = prod_vect_triangle(c[1][1],c[1][0],c[0][0]);
  const double v11 = prod_vect_triangle(c[1][1],c[1][0],c[0][1]);
  if (v10 * v11 > 0.)
    return 0;

  return 1;
}

// Description:
//  Teste si les facettes fa70 et fa71 se coupent.
//  Renvoie 0 si elles ne se coupent pas,
//         -1 si elles ont deux sommets communs,
//         -2 si elles ont trois sommets communs,
//         -3 si le test preliminaire est positif et l'autre test non
//         1  si les facettes se coupent et n'ont aucun sommet commun
//         -4 si les facettes ont un sommet commun.
// L'algorithme est identique a celui decrit dans
//  http://www.cs.lth.se/home/Tomas_Akenine_Moller/pubs/tritri.pdf
//  "a fast triangle-triangle intersection test (tomas moller)"
int Topologie_Maillage_FT::test_intersection_facettes_3D(
  int fa70, int fa71,
  const Maillage_FT_Disc& maillage) const
{
  // LA PREMIERE PARTIE DU TEST EST SPECIFIQUE AU FRONT-TRACKING:
  //  ignorer les intersections entre facettes voisines.

  // Copie locale des indices des sommets des facettes
  int facette0[3];
  int facette1[3];
  {
    const IntTab& facettes = maillage.facettes();
    facette0[0] = facettes(fa70,0);
    facette0[1] = facettes(fa70,1);
    facette0[2] = facettes(fa70,2);
    facette1[0] = facettes(fa71,0);
    facette1[1] = facettes(fa71,1);
    facette1[2] = facettes(fa71,2);
  }
  // On teste si les facettes ont deux sommets en commun:
  //int liste_sommets_communs[3];
  int nb_sommets_communs = 0;
  trier_trois_entiers(facette0);
  trier_trois_entiers(facette1);
  {
    int i0 = 0;
    int i1 = 0;
    while (i0 < 3 && i1 < 3)
      {
        const int f0 = facette0[i0];
        const int f1 = facette1[i1];
        if (f0 <= f1)
          i0++;
        if (f0 >= f1)
          i1++;
        if (f0 == f1)  // On a trouve un sommet commun :
          {
            // liste_sommets_communs[nb_sommets_communs] = f0;
            nb_sommets_communs++;
          }
      }
  }
  // Trois sommets communs: c'est une facette dupliquee.
  // Pas d'intersection...
  if (nb_sommets_communs == 3)
    {
      return -3;
    }
  // Si les facettes sont voisines par une arete,
  // on dit qu'elles n'ont pas d'intersection.
  if (nb_sommets_communs == 2)
    {
      return -2;
    }
  // Essai: si les facettes ont un sommet commun,
  // on dit qu'elles n'ont pas d'intersection
  if (nb_sommets_communs == 1)
    {
      return -4;
    }

  // ICI COMMENCE L'ALGORITHME DE TOMAS MOLLER...

  // Un sommet ou moins, on cherche une intersection geometrique
  // On recupere les coordonnees des sommets
  // coord[i][j] = coordonnee j du sommet i
  double coord_0[3][3];
  double coord_1[3][3];
  {
    const DoubleTab& sommets = maillage.sommets();
    get_coord_sommets_triangle(facette0, sommets, coord_0);
    get_coord_sommets_triangle(facette1, sommets, coord_1);
  }
  // Equation du plan contenant la facette
  // Le plan s'ecrit a[0]*x+a[1]*y+a[2]*z+a[3]=0
  double equation_plan_0[4];
  double equation_plan_1[4];
  {
    const DoubleTab& facettes_normales = maillage.get_update_normale_facettes();
    calcul_equation_plan(coord_0[0], facettes_normales, fa70, equation_plan_0);
    calcul_equation_plan(coord_1[0], facettes_normales, fa71, equation_plan_1);
  }
  // Pour chaque sommet, calcul de la distance au plan
  // equation du plan 0 appliquee aux sommets de la facette 1
  double distance_plan0_som1[3];
  // equation du plan 1 appliquee aux sommets de la facette 0
  double distance_plan1_som0[3];
  calcul_distance_plan(coord_1, equation_plan_0, distance_plan0_som1);
  calcul_distance_plan(coord_0, equation_plan_1, distance_plan1_som0);
  // Y a-t-il changement de signe (sommets de part et d'autre du plan) ?
  int test_intersection = 0;
  {
    int n_positifs_0 = 0;
    int n_positifs_1 = 0;
    int i;
    for (i = 0; i < 3; i++)
      {
        n_positifs_0 += (distance_plan0_som1[i] > 0.);
        n_positifs_1 += (distance_plan1_som0[i] > 0.);
      }
    if ((n_positifs_0 * n_positifs_1) % 3 != 0)
      {
        test_intersection = -3;
      }
  }
  // Changement de signe => possibilite d'une intersection, on fait le test complet.
  if (test_intersection)
    {
      static const double TOLERANCE_FACETTES_PARALLELES = 1e-12;
      // Equation d'un plan orthogonal aux deux facettes
      //  normale = produit vectoriel des normales aux facettes
      double equation_plan_2[4];
      double a = equation_plan_0[1] * equation_plan_1[2] - equation_plan_0[2] * equation_plan_1[1];
      double b = equation_plan_0[2] * equation_plan_1[0] - equation_plan_0[0] * equation_plan_1[2];
      double c = equation_plan_0[0] * equation_plan_1[1] - equation_plan_0[1] * equation_plan_1[0];
      equation_plan_2[0] = a;
      equation_plan_2[1] = b;
      equation_plan_2[2] = c;
      equation_plan_2[3] = 0.; // Constante : indifferent
      double norme_carre = a * a + b * b + c * c;
      // Les normales aux facettes sont unitaires, norme_carre est donc adimensionnel,
      // c'est le carre du sinus de l'angle entre les facettes.
      if (norme_carre < TOLERANCE_FACETTES_PARALLELES)
        {
          // Les deux facettes sont paralleles et coplanaires
          test_intersection = 0;
        }
      else
        {
          double distance_plan2_som0[3];
          double distance_plan2_som1[3];
          calcul_distance_plan(coord_0, equation_plan_2, distance_plan2_som0);
          calcul_distance_plan(coord_1, equation_plan_2, distance_plan2_som1);
          // Calcul de l'intersection entre la facette_0 et le plan contenant facette_1:
          // c'est un segment dont les extremites sont a une 'distance' s0_min et s0_max
          // du plan_2
          double s0_min, s0_max, s1_min, s1_max;
          calcul_intersection_plan_triangle(distance_plan1_som0, distance_plan2_som0,
                                            s0_min, s0_max);
          // Idem avec la facette_1 et le plan contenant facette_0:
          calcul_intersection_plan_triangle(distance_plan0_som1, distance_plan2_som1,
                                            s1_min, s1_max);
          // Les deux segments ont-ils une intersection non vide ? Si oui alors les
          // triangles se coupent.
          {
            double s_min = (s0_min > s1_min) ? s0_min : s1_min;
            double s_max = (s0_max < s1_max) ? s0_max : s1_max;
            if (s_min < s_max)
              // Intersection non vide => les facettes se coupent
              test_intersection = 1 + nb_sommets_communs;
          }
        }
    }
  return test_intersection;
}

// Description:
//  Teste s'il existe deux facettes du maillage qui se coupent (avec
//   test_intersection_facettes_2D/3D). On ne teste pas tous les couples
//   possibles, seulement les couples de facettes qui coupent un meme
//   element eulerien.
// Precondition: le maillage doit etre parcouru.
int Topologie_Maillage_FT::test_collision_facettes(const Maillage_FT_Disc& maillage,
                                                   ArrOfInt& liste_elements_collision) const
{
  // Note BM 08/2009: on peut reecrire cet algorithme en utilisant un OctreeDouble des
  //  facettes du maillage. Cela permet de ne pas utiliser le maillage eulerien, temps cpu
  //  independant du rapport lagrangien/eulerien et on s'affranchit des erreurs potentielles sur le parcours.
  //  Attention quand-meme a avoir les mailles virtuelles pour ne pas oublier des collisions
  //  avec des mailles d'un autre proc.

  int i_facette;
  const int nb_facettes = maillage.nb_facettes();
  ArrOfBit facettes_testees(nb_facettes);
  facettes_testees = 0;
  ArrOfIntFT liste_facettes_testees;
  const Intersections_Elem_Facettes& intersection =
    maillage.intersections_elem_facettes();

  const ArrOfInt& index_elem = intersection.index_elem();
  const ArrOfInt& index_facettes = intersection.index_facette();
  int collision_trouvee = 0;

  liste_elements_collision.set_smart_resize(1);
  liste_elements_collision.resize_array(0);

  // Methode de recherche des intersections:
  //  pour chaque facette,
  //    quels elements euleriens sont traverses ?
  //    pour chaque element eulerien traverse,
  //       quelles sont les facettes_2 qui coupent cet element ?
  //       pour chaque facette_2,
  //          cette facette_2, coupe-t-elle la facette ?
  //  pour ne pas tester plusieurs fois le meme couple:
  //    on ne teste que si facette_2 > facette
  //    on maintient pour chaque facette un tableau de marqueurs
  //      des facettes_2 deja testees.

  for (i_facette = 0; i_facette < nb_facettes && !collision_trouvee; i_facette++)
    {

      // Boucle sur les elements traverses par la facette
      int index = index_facettes[i_facette];
      while (index >= 0 && !collision_trouvee)
        {
          const Intersections_Elem_Facettes_Data& data = intersection.data_intersection(index);
          // -----------------------------------------------------
          // Numero de l'element traverse
          const int i_elem = data.numero_element_;
          // Boucle sur les facettes traversees par l'element
          int index2 = index_elem[i_elem];
          while (index2 >= 0 && !collision_trouvee)
            {
              const Intersections_Elem_Facettes_Data& data2 =
                intersection.data_intersection(index2);
              const int i_facette2 = data2.numero_facette_;
              // On ne teste que les couples (i_facette2 > i_facette)
              if (i_facette2 > i_facette)
                {
                  if (!facettes_testees.testsetbit(i_facette2))
                    {
                      // La facette n'a pas encore ete testee
                      liste_facettes_testees.append_array(i_facette2);
                      int resu;
                      if (Objet_U::dimension == 3)
                        resu = test_intersection_facettes_3D(i_facette, i_facette2, maillage);
                      else
                        resu = test_intersection_facettes_2D(i_facette, i_facette2, maillage);
                      if (resu > 0)
                        {
                          // Intersection si le resultat est positif.
                          collision_trouvee = resu > 0;
                          liste_elements_collision.append_array(i_elem);
                        }
                    }
                }
              index2 = data2.index_facette_suivante_;
            }
          // ----------------------------------------------------
          index = data.index_element_suivant_;
        }
      // Efface les marqueurs des facettes testees
      {
        int i;
        const int n = liste_facettes_testees.size_array();
        for (i = 0; i < n; i++)
          {
            int i_facette2 = liste_facettes_testees[i];
            facettes_testees.clearbit(i_facette2);
          }
        liste_facettes_testees.resize_array(0);
      }
    }
  array_trier_retirer_doublons(liste_elements_collision);
  collision_trouvee = (int)mp_max(collision_trouvee);
  return collision_trouvee;
}

// Description:
//  Computes eulerian connex components of the phase indicator function "indicatrice"
//  according to the get_phase_continue() property. This method must be used to
//  compute the input data for supprimer_interfaces.
//
//  For each eulerian mesh cell
//  num_compo[i] = -1 for all mesh cells containing phase "get_phase_continue()"
//  num_compo[i] = N  with 0 <= N < Nmax, for all other mesh cells.
//   N is a global connex component number. Two adjacent cells (by face) for
//   which indicatrice[i] != get_phase_continue() will have the same N.
// Parametre: num_compo
// Signification:
//  output: resized, md_vector setup, filled with -1 <= num_compo[i] < Nmax, and virtual space updated.
// Valeur de retour:
//  Returns Nmax, global number of connex components found.
int Topologie_Maillage_FT::calculer_composantes_connexes_pour_suppression(const Zone_VF& zone_vf,
                                                                          const DoubleTab& indicatrice,
                                                                          IntVect& num_compo) const
{
  zone_vf.zone().creer_tableau_elements(num_compo);

  double phase_cont = get_phase_continue();
  const int nb_elem = zone_vf.zone().nb_elem();
  // On marque les elements dont on veut avoir les composantes connexes:
  //  -1 pour la phase a conserver (ne pas ranger ces elements dans une compo connexe)
  //  1 pour la phase a supprimer
  int i;
  for (i = 0; i < nb_elem; i++)
    {
      const double indic = indicatrice[i];
      // Pas de epsilon pour ce test: on attend une egalite exacte avec 0. ou 1.
      // sinon c'est une maille diphasique qui contient une interface.
      num_compo[i] = (indic == phase_cont) ? -1 : 1;
    }
  num_compo.echange_espace_virtuel();
  const int nb_local_components = search_connex_components_local(zone_vf.elem_faces(), zone_vf.face_voisins(), num_compo);
  const int nb_compo_glob = compute_global_connex_components(num_compo, nb_local_components);

  return nb_compo_glob;
}

// Description: Removes all interfaces contained in eulerian elements marked
//  by the "flags_" array, and updates the "indicatrice" field by putting
//  "get_phase_condinue()" in those elements.
//  Virtual lagrangian elements are removed.
//  "flags" must be of size "nb_elem" and must be built with one or more
//  connex components (see calcul_composantes_connexes_eurleriennes).
// Return value: Integral of the "indicatrice" change (e.g. volume of phase changed,
//  positive or negative)
double Topologie_Maillage_FT::suppression_interfaces(const IntVect& num_compo,
                                                     const ArrOfInt& flags_compo_a_supprimer,
                                                     Maillage_FT_Disc& maillage,
                                                     DoubleTab& indicatrice)
{
  const int n = num_compo.size();
  double phase_cont = get_phase_continue();

  assert(n == indicatrice.dimension(0));

  const int nb_facettes = maillage.nb_facettes();
  ArrOfInt liste_facettes_a_supprimer;
  array_smart_allocate(liste_facettes_a_supprimer, nb_facettes);

  const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
  const ArrOfInt& index_elem = intersections.index_elem();
  ArrOfBit flags(nb_facettes);
  flags = 0;

  double variation_indicatrice = 0.;

  for (int i = 0; i < n; i++)
    {
      const int compo_connexe = num_compo[i];
      if (compo_connexe >= 0 && flags_compo_a_supprimer[compo_connexe])
        {
          // On vide la maille i:
          variation_indicatrice += phase_cont - indicatrice[i];
          indicatrice[i] = phase_cont;
          // Parcours des facettes traversant l'element i:
          int index = index_elem[i];
          while (index >= 0)
            {
              const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
              const int i_facette = data.numero_facette_;
              if (!flags.testsetbit(i_facette))
                liste_facettes_a_supprimer.append_array(i_facette);
              index = data.index_facette_suivante_;
            }
        }
    }
  maillage.supprimer_facettes(liste_facettes_a_supprimer);
  variation_indicatrice = mp_sum(variation_indicatrice);
  declare_espace_virtuel_invalide(indicatrice);
  return variation_indicatrice;
}

// Description:
//  Remaillage de l'interface:
//   - amelioration petites et grandes facettes,
//   - barycentrage,
//   - gestion des coalescences-fragmentations.
void Topologie_Maillage_FT::remailler_interface(const double temps,
                                                Maillage_FT_Disc& maillage,
                                                Champ_base& indicatrice,
                                                Remaillage_FT& algo_remaillage_local)
{
  // Le remaillage qui suit le traitement des coalescences peut a nouveau
  //  produire un maillage impropre (avec des facettes qui se coupent).
  // Dans ce cas, on s'arrete apres le traitement des coalescences.

  // Note B.M. : le remaillag systematique peut retarder la necessite de remailler
  //  localement (apparition de grandes ou petites aretes). Donc je le fais d'abord,
  //  et ensuite je teste s'il faut faire un remaillage local.

  // L'intervalle de temps entre deux lissages est-il ecoule ?
  if (algo_remaillage_local.a_lisser(temps))
    {
      algo_remaillage_local.barycentrer_lisser_systematique(temps,maillage);
    }

  // L'intervalle de temps entre deux remaillages locaux est-il ecoule:
  if (algo_remaillage_local.a_remailler(temps, maillage))
    {
      // Declanchement d'un remaillage local.
      // Pour que ca marche bien, les parametres de barycentrage et lissage "apres_remaillage"
      // doivent etre suffisants pour ramener le maillage dans un etat correct sans
      // appliquer de lissage systematique apres
      algo_remaillage_local.remaillage_local_interface(temps, maillage);
    }

  // Test de collision des interfaces (peut declancher un remaillage "global")
  if (active_)
    {
      ArrOfInt liste_elements_collision;
      for (int num_tentative = 0; num_tentative < 2; num_tentative++)
        {
          maillage.parcourir_maillage();
          const int collision_trouvee = test_collision_facettes(maillage, liste_elements_collision);
          if (!collision_trouvee)
            {
              // Pas de coalescence: on quitte la boucle tout de suite
              break;
            }

          if (Process::je_suis_maitre())
            Journal() << "Collision t= " << temps << finl;
          Remailleur_Collision_FT_base& remesh = remailleur_Collision_.valeur();
          if (juric_local_)
            {
              // Indicatrice de la partie du maillage a conserver (ie a ne pas remailler)
              // (COPIE !)
              Maillage_FT_Disc maillage2;
              maillage2.recopie(maillage, Maillage_FT_Disc::MINIMAL);
              // On fait une copie car on va modifier indicatrice et on veut avoir acces a la valeur precedente:
              DoubleTab indicatrice_partie_non_remaillee = indicatrice.valeurs();
              // Supprimer les triangles des composantes connexes ou se trouvent les collisions
              const Zone_VF& zone_vf = ref_cast(Zone_VF, indicatrice.zone_dis_base());
              {
                IntVect num_compo;
                const int nb_compo = calculer_composantes_connexes_pour_suppression(zone_vf, indicatrice_partie_non_remaillee, num_compo);
                // Selection des composantes a supprimer, initialisation a zero
                ArrOfInt flags_compo_a_supprimer(nb_compo);
                const int n = liste_elements_collision.size_array();
                for (int i = 0; i < n; i++)
                  {
                    const int elem = liste_elements_collision[i];
                    const int compo = num_compo[elem];
                    if (compo >= 0)
                      flags_compo_a_supprimer[compo] = 1;
                  }
                // Tous les procs d'accord sur les composantes a supprimer:
                mp_max_for_each_item(flags_compo_a_supprimer);
                suppression_interfaces(num_compo,
                                       flags_compo_a_supprimer,
                                       maillage,
                                       indicatrice_partie_non_remaillee);
              }
              // Constuire une fonction indicatrice contenant uniquement les interfaces a remailler
              DoubleTab& indic = indicatrice.valeurs();
              const int nb_elem = indic.dimension(0);
              int i;
              const double phase_cont = get_phase_continue();
              for (i = 0; i < nb_elem; i++)
                {
                  const double ind = indic[i];
                  const double indic2 = indicatrice_partie_non_remaillee[i];
                  if (ind == indic2)
                    indic[i] = phase_cont;
                }
              indic.echange_espace_virtuel();
              // Recontruire la portion de maillage dans une structure separee
              maillage2.parcourir_maillage();
              if (num_tentative == 0)
                remesh.traite_RuptureCoalescenceInterfaces_Conservatif(maillage2, indicatrice);
              else
                remesh.traite_RuptureCoalescenceInterfaces(maillage2, indicatrice);
              // Remaillage local du maillage2:
              algo_remaillage_local.remaillage_local_interface(temps, maillage2);
              // Concatener les deux maillages
              maillage.ajouter_maillage(maillage2);
              // Concatener les indicatrices
              for (i = 0; i < nb_elem; i++)
                {
                  const double indic2 = indicatrice_partie_non_remaillee[i];
                  if (indic2 != phase_cont)
                    indic[i] = indic2;
                }
            }
          else
            {
              if (num_tentative == 0)
                {
                  remesh.traite_RuptureCoalescenceInterfaces_Conservatif(maillage, indicatrice);
                  // Remaillage local:
                  algo_remaillage_local.remaillage_local_interface(temps, maillage);
                }
              else
                {
                  remesh.traite_RuptureCoalescenceInterfaces(maillage, indicatrice);
                  // Le premier remaillage local a reproduit une coalescence, ne pas remailler localement
                }
            }
        }
    }
}
