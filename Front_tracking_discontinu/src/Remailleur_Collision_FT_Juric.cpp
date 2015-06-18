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
// File:        Remailleur_Collision_FT_Juric.cpp
// Directory:   $TRUST_ROOT/src/Front_tracking_discontinu
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Remailleur_Collision_FT_Juric.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Probleme_base.h>
#include <Domaine.h>
#include <Param.h>
#include <Scatter.h>

Implemente_instanciable_sans_constructeur(Remailleur_Collision_FT_Juric,"Remailleur_Collision_FT_Juric",Remailleur_Collision_FT_base);

Remailleur_Collision_FT_Juric::Remailleur_Collision_FT_Juric() :
  source_isovaleur_(INDICATRICE) // Valeur par defaut: la plus robuste.
{
}

Sortie& Remailleur_Collision_FT_Juric::printOn(Sortie& os) const
{
  return os;
}

// Description: Lecture du parametre Source_Isovaleur.
//  Format attendu :
//    { [ source_isovaleur indicatrice | fonction_distance ] }
Entree& Remailleur_Collision_FT_Juric::readOn(Entree& is)
{
  source_isovaleur_ = INDICATRICE;
  Param param(que_suis_je());
  int src = INDICATRICE;
  param.ajouter("source_isovaleur", & src);
  param.dictionnaire("indicatrice", INDICATRICE);
  param.dictionnaire("fonction_distance", FONCTION_DISTANCE);
  param.lire_avec_accolades_depuis(is);
  if (src == INDICATRICE)
    {
      source_isovaleur_ = INDICATRICE;
      Cerr << " Using indicator function to reconstruct interfaces" << finl;
    }
  else
    {
      source_isovaleur_ = FONCTION_DISTANCE;
      Cerr << " Using distance function to reconstruct interfaces (DANGEROUS: not robust)" << finl;
    }
  return is;
}

static void calculer_iso_avec_distance(const Domaine&     domaine,
                                       const Champ_base& indicatrice,
                                       const DoubleTab&   distance,
                                       DoubleTab&         fonction_iso)
{
  fonction_iso = distance;

  // Remplissage des valeurs non calculees de la distance :
  //  on met une valeur positive si l'indicatrice vaut 1 et
  // negative si l'indicatrice vaut 0
  int i;
  const int nb_som = distance.dimension(0);
  DoubleTab indicatrice_sommets(nb_som,1);
  indicatrice.valeur_aux_sommets(domaine, indicatrice_sommets);
  const double valeur_invalide = -1.E30;
  for (i = 0; i < nb_som; i++)
    {
      double d = fonction_iso(i);
      if (d < valeur_invalide)
        {
          double indic = indicatrice_sommets(i,0);
          d = (indic < 0.5) ? -1.e6 : 1.e6;
          fonction_iso(i) = d;
        }
    }
  fonction_iso.echange_espace_virtuel();
}

static void calculer_iso_avec_indicatrice(const Domaine&     domaine,
                                          const Champ_base& indicatrice,
                                          DoubleTab&         fonction_iso)
{
  const int nb_som = fonction_iso.dimension(0);
  // Creation d'un tableau temporaire a deux dimensions
  // (attendu par valeur_aux_sommets())
  DoubleTab f(nb_som,1);
  indicatrice.valeur_aux_sommets(domaine, f);
  int i;
  for (i = 0; i < nb_som; i++)
    {
      fonction_iso(i) = f(i, 0) - 0.5;
    }
  fonction_iso.echange_espace_virtuel();
}

// Description:
//  Remaillage complet des interfaces par la methode de Juric:
//  On calcule une fonction distance a l'interface evaluee aux sommets
//  du maillage eulerien et on construit une isovaleur de cette fonction distance.
int Remailleur_Collision_FT_Juric::traite_RuptureCoalescenceInterfaces(
  Maillage_FT_Disc& maillage,
  Champ_base& indicatrice) const
{
  int res = 1;
  Process::Journal()<<"DEB Remailleur_Collision_FT_Juric::traite_RuptureCoalescenceInterfaces"<<finl;

  const Transport_Interfaces_FT_Disc& eq_transport_FT = maillage.equation_transport();
  const Domaine& domaine = eq_transport_FT.zone_dis().zone().domaine();

  DoubleTab fonction_iso;

  switch(source_isovaleur_)
    {
    case INDICATRICE:
      {
        domaine.creer_tableau_sommets(fonction_iso);
        calculer_iso_avec_indicatrice(domaine, indicatrice, fonction_iso);
        break;
      }
    case FONCTION_DISTANCE:
      {
        const DoubleTab& distance = eq_transport_FT.get_update_distance_interface_sommets();
        calculer_iso_avec_distance(domaine, indicatrice, distance, fonction_iso);
        break;
      }
    default:
      Cerr << "Erreur dans Remailleur_Collision_FT_Juric::traite_RuptureCoalescenceInterfaces" << finl;
      Cerr << " option source_isovaleur_ non implementee" << finl;
      assert(0);
      exit();
    }

  const Marching_Cubes& marching_cubes = eq_transport_FT.marching_cubes();

  marching_cubes.construire_iso(fonction_iso,         // Le champ dont on prend l'isovaleur
                                0.,                   // La valeur
                                maillage,
                                indicatrice.valeurs(),  // L'indicatrice est remise a jour dans les phases
                                Maillage_FT_Disc::AJOUTE_TOUTES_PHASES);

  Process::Journal()<<"FIN Remailleur_Collision_FT_Juric::traite_RuptureCoalescenceInterfaces res="<<res<<finl;
  return res;
}

int Remailleur_Collision_FT_Juric::traite_RuptureCoalescenceInterfaces(Maillage_FT_Disc& maillage,
                                                                       int nb_fa7Intersectees,
                                                                       const IntTab& couples_fa7Intersectees,
                                                                       const DoubleTab& segmentsInter_fa7Intersectees,
                                                                       Champ_base& indicatrice)
{
  int res = traite_RuptureCoalescenceInterfaces(maillage, indicatrice);
  return res;
}
