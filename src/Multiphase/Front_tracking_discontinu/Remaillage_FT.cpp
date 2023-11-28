/****************************************************************************
* Copyright (c) 2022, CEA
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
// File:        Remaillage_FT.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Multiphase/Front_tracking_discontinu/src
// Version:     /main/patch_168/1
//
//////////////////////////////////////////////////////////////////////////////

#include <Remaillage_FT.h>
#include <TRUST_Deriv.h>
#include <Motcle.h>
#include <Domaine_VF.h>
#include <Domaine.h>
#include <Triangle.h>
#include <Rectangle.h>
#include <Rectangle_2D_axi.h>
#include <Hexaedre.h>
#include <Tetraedre.h>
#include <Comm_Group.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <EcritureLectureSpecial.h>
#include <Param.h>
#include <Array_tools.h>
#include <DebogFT.h>
#include <stat_counters.h>

Implemente_instanciable_sans_constructeur(Remaillage_FT,"Remaillage_FT",Objet_U);


Remaillage_FT::Remaillage_FT() :
  // Values initialized:
  temps_(0),
  temps_dernier_remaillage_(-1.e40),
  temps_dernier_lissage_(-1.e40),
  // Values initialized but which can be set in the data file:
  dt_remaillage_(-1),
  dt_lissage_(-1),
  nb_iter_remaillage_(1),
  nb_iter_barycentrage_(1),
  nb_iter_bary_volume_seul_(3),
  seuil_dvolume_residuel_(0.),
  relax_barycentrage_(1.),
  critere_arete_(0.35),
  impr_(-1),
  valeur_longueur_fixe_(-1.),
  facteur_longueur_ideale_(-1.),
  equilateral_(0), // Par defaut, on regarde
  //                  l'orientation de la facette pour ajuster sa taille a l'element eulerien.
  //                  alors que equilateral=1 utilise la diagonal de l'element comme longueur de reference.
  //                  Avec equilateral_ = 0, les facettes sont etirees comme le maillage et facteur_longueur vaut
  //                  1. si l'arete traverse tout l'element dans la direction donnee par l'arrete.
  lissage_courbure_coeff_(-0.05), // valeur typique pour stabilite
  lissage_courbure_iterations_systematique_(0),
  lissage_courbure_iterations_si_remaillage_(0),
  lissage_courbure_iterations_old_(-1),
  lissage_critere_(0) // Default value to 0, when lissage is applied, it is for the whole mesh
{
}

/*! @brief Lecture des parametres de remaillage
 *
 */
Entree& Remaillage_FT::readOn(Entree& is)
{
  Param p(que_suis_je());
  set_param(p);
  double critere_remaillage_compat=-1;
  p.ajouter("critere_remaillage",& critere_remaillage_compat);
  p.ajouter("equilateral", &equilateral_);
  facteur_longueur_ideale_ = -1.;
  valeur_longueur_fixe_ = -1.;
  if ( critere_remaillage_compat!=-1)
    {
      Cerr<<"Warning : critere_remaillage is obsolete "<<finl;
    }
  p.lire_avec_accolades_depuis(is);

  if (lissage_courbure_iterations_old_>=0)
    {
      Cerr << "Remaillage_FT: old keyword lissage_courbure_iterations:\n"
           << " setting lissage_courbure_iterations_systematique and lissage_courbure_iterations_si_remaillage to same value" << finl;
      lissage_courbure_iterations_si_remaillage_ = lissage_courbure_iterations_old_;
      lissage_courbure_iterations_systematique_ = lissage_courbure_iterations_old_;
    }

  if (facteur_longueur_ideale_ > 0. && valeur_longueur_fixe_ == -1.)
    {
      // ok on utilisera facteur_longueur_ideale_
    }
  else if (facteur_longueur_ideale_ == -1. && valeur_longueur_fixe_ > 0.)
    {
      // ok on utilisera longueur fixe
    }
  else if (facteur_longueur_ideale_ == -1. && valeur_longueur_fixe_ == -1.)
    {
      // Valeur par defaut:
      if (dimension == 3)
        facteur_longueur_ideale_ = 1.5;
      else
        facteur_longueur_ideale_ = 2.;
    }
  else
    {
      Cerr << "Error in Remaillage_FT::readOn: It is an error to specify both facteur_longueur_ideale"
           << " and critere_longueur_fixe parameters." << finl;
      exit();
    }
  return is;
}

/*! @brief Cette fonction ne doit pas etre appelee
 *
 */
Sortie& Remaillage_FT::printOn(Sortie& os) const
{
  Cerr << "Erreur : ::printOn n'est pas code." << finl;
  assert(0);
  return os;
}

int Remaillage_FT::reprendre(Entree& is)
{
  Nom motlu;
  is >> motlu;
  if (motlu != que_suis_je())
    {
      Cerr << "Erreur dans Remaillage_FT::reprendre\n";
      Cerr << " On attendait le motcle " << que_suis_je();
      Cerr << "\n On a trouve " << motlu << finl;
      assert(0);
      exit();
    }
  is >> temps_dernier_remaillage_;
  is >> temps_dernier_lissage_;
  return 1;
}

int Remaillage_FT::sauvegarder(Sortie& os) const
{
  int special, afaire;
  const int format_xyz = EcritureLectureSpecial::is_ecriture_special(special, afaire);
  if (format_xyz)
    {
      if (Process::je_suis_maitre())
        {
          os << que_suis_je() << finl;
          os << temps_dernier_remaillage_ << finl;
          os << temps_dernier_lissage_ << finl;
        }
    }
  else
    {
      os << que_suis_je() << finl;
      os << temps_dernier_remaillage_ << finl;
      os << temps_dernier_lissage_ << finl;
    }
  return 0;
}

/*! @brief Methode appelee par readOn.
 *
 * Declaration des parametres a lire dans le .data
 *
 */
void Remaillage_FT::set_param(Param& p)
{
  p.ajouter("pas", &dt_remaillage_);
  p.ajouter("pas_lissage", &dt_lissage_);
  p.ajouter("nb_iter_remaillage", &nb_iter_remaillage_);
  p.ajouter("nb_iter_barycentrage", &nb_iter_barycentrage_);
  p.ajouter("nb_iter_correction_volume",  &nb_iter_bary_volume_seul_);
  p.ajouter("seuil_dvolume_residuel", &seuil_dvolume_residuel_);
  p.ajouter("relax_barycentrage", &relax_barycentrage_);
  p.ajouter("critere_arete", &critere_arete_);
  p.ajouter("impr", &impr_);
  p.ajouter("critere_longueur_fixe", &valeur_longueur_fixe_);
  p.ajouter("facteur_longueur_ideale", &facteur_longueur_ideale_);
  p.ajouter("lissage_courbure_coeff", &lissage_courbure_coeff_);
  // lissage_courbure_iterations est un synonyme de lissage_courbure_iterations_systematique
  p.ajouter("lissage_courbure_iterations", &lissage_courbure_iterations_old_);
  p.ajouter("lissage_courbure_iterations_systematique", &lissage_courbure_iterations_systematique_);
  p.ajouter("lissage_courbure_iterations_si_remaillage", &lissage_courbure_iterations_si_remaillage_);
  p.ajouter("lissage_critere", &lissage_critere_);
}

/*! @brief Cette fonction stocke le domaine_dis dans refdomaine_dis_
 *
 * @param (domaine_dis) domaine discretisee de calcul
 * @return le flot d'entree
 */
void Remaillage_FT::associer_domaine(const Domaine_dis& domaine_dis)
{
  Cerr<<"Remaillage_FT::associer_domaine_dis"<<finl;
  refdomaine_VF_ = ref_cast(Domaine_VF,domaine_dis.valeur());
}

int Remaillage_FT::a_remailler(double temps, const Maillage_FT_Disc& maillage) const
{
  int res = 0;
  if (dt_remaillage_ > 0.)
    {
      if (temps > (temps_dernier_remaillage_ + dt_remaillage_) * (1.-1e-15))
        {
          res = 1;
        }
    }
  return res;
}

int Remaillage_FT::a_lisser(double temps) const
{
  int res = 0;
  if ((dt_lissage_ > 0.) && (temps > (temps_dernier_lissage_ + dt_lissage_) * (1.-1e-15)))
    {
      res = 1;
    }
  return res;
}

/*! @brief Verifie les criteres de remaillage locaux (longueur des aretes) et effectue les remaillages locaux si necessaires.
 *
 */
void Remaillage_FT::remaillage_local_interface(double temps, Maillage_FT_Disc& maillage)
{
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "Remaillage_local_interface", "FrontTracking");
  statistiques().begin_count(stat_counter);

  temps_dernier_remaillage_ = temps_dernier_lissage_ = temps_ = temps;

  maillage.nettoyer_elements_virtuels();
  maillage.check_mesh();
  //boucle sur les remaillages
  int iter;
  ArrOfDoubleFT varVolume;
  for (iter = 0; iter < nb_iter_remaillage_; iter++)
    {
      const int nb_sommets = maillage.nb_sommets();
      varVolume.resize_array(nb_sommets);
      varVolume = 0.;
      variation_volume_ = 0.;
      // n = nombre de sommets supprimes
      const int n = supprimer_petites_aretes(maillage,varVolume);
      if (Comm_Group::check_enabled())
        maillage.check_mesh();
      if (n > 0)
        {
          supprimer_doublons_facettes(maillage); // Marque les facettes mais ne les supprime pas effectivement.
          if (Comm_Group::check_enabled())
            maillage.check_mesh();
          // On a supprime les petites aretes en deplacant un noeud sur
          //  un autre, la variation de volume engendree a ete mise dans varVolume:
          // On barycentre et on lisse, notamment pour recuperer cette variation.
          barycentrer_lisser_apres_remaillage(maillage, varVolume);
          if (Comm_Group::check_enabled())
            maillage.check_mesh();
        }
      const int m = diviser_grandes_aretes(maillage);
      if (Comm_Group::check_enabled())
        maillage.check_mesh();
      if (m > 0)
        {
          varVolume.resize_array(maillage.nb_sommets());
          varVolume = 0.;
          barycentrer_lisser_apres_remaillage(maillage, varVolume);
          if (Comm_Group::check_enabled())
            maillage.check_mesh();
        }
      if (Process::je_suis_maitre())
        Journal() << "remaillage_local_interface t= " << temps << " suppressions: " << n << " divisions: " << m << finl;
    }
  nettoyer_maillage(maillage);
  statistiques().end_count(stat_counter);
}


/*! @brief Cette fonction calcule les connectivites sommet ->facettes voisines Les facettes voisines des sommets sont stockees dans le tableau fa7VoisinesSom_data
 *
 *     sous forme de liste chainee.
 *     Pour le sommet som, le premier index de la liste est index = fa7VoisinesSom_index[som]
 *       la premiere facette est alors fa7VoisinesSom_data(index,1)
 *       et l'indice suivant index = fa7VoisinesSom_data(index,0)
 *     La chaine est terminee lorsque index==-1.
 *
 * @param (maillage) maillage a barycentrer
 * @param (fa7VoisinesSom_index) premier index pour le sommet som
 * @param (fa7VoisinesSom_data) premier liste des index et facettes voisines pour le sommet som
 * @return (int) le nombre de connectivites trouvees
 */
int Remaillage_FT::calculer_connectivites_sommetFacettes(const Maillage_FT_Disc& maillage,
                                                         ArrOfInt& fa7VoisinesSom_index,
                                                         IntTab& fa7VoisinesSom_data) const
{
  int compteur = 0;

  const int nb_som = maillage.nb_sommets();
  const int nb_facettes = maillage.nb_facettes();
  const IntTab& facettes = maillage.facettes();
  const int nb_som_par_facette = facettes.dimension(1);

  int fa7, isom, som, trouve, index, index0 = -1;

  //initialisation
  for (som=0 ; som<nb_som ; som++)
    {
      fa7VoisinesSom_index[som] = -1;
    }

  //on balaye les facettes
  for (fa7=0 ; fa7<nb_facettes ; fa7++)
    {
      //puis on balayes les sommets des facettes
      for (isom=0 ; isom<nb_som_par_facette ; isom++)
        {
          som = facettes(fa7,isom);
          //on va ensuite balayer la liste des facettes voisines du sommet som
          //pour voir si la facette fa7 est deja stockee
          index = fa7VoisinesSom_index[som];
          if (index==-1)
            {
              //le sommet n'a pas encore de fa7 voisine stockee :
              //l'ajoute a la premiere place dispo = compteur
              fa7VoisinesSom_index[som] = compteur;
              fa7VoisinesSom_data(compteur,0) = -1;
              fa7VoisinesSom_data(compteur,1) = fa7;
              compteur++;
            }
          else
            {
              //sinon : balaye  la liste des facettes voisines du sommet som
              trouve = -1;
              while (index>=0 && trouve==-1)
                {
                  index0 = index;
                  if (fa7VoisinesSom_data(index,1) == fa7)
                    {
                      trouve = index;
                    }
                  else
                    {
                      index = fa7VoisinesSom_data(index,0);
                    }
                }
              if (trouve==-1)
                {
                  //fa7 non encore stockee
                  assert(index==-1 && index0!=-1);
                  //l'ajoute a la premiere place dispo = compteur
                  fa7VoisinesSom_data(index0,0) = compteur;
                  fa7VoisinesSom_data(compteur,0) = -1;
                  fa7VoisinesSom_data(compteur,1) = fa7;
                  compteur++;
                }
            }
        }
    }

  return compteur;
}

/*! @brief Calcul de la differentielle du volume de phase 0 par rapport au deplacement de chaque sommet de l'interface.
 *
 * C'est une forme
 *   lineaire qu'on exprime sous la forme d'un vecteur v tel que
 *   differentielle(deplacement) = deplacement scalaire v.
 *   En un certain sens, v est la normale a l'interface evaluee aux sommets.
 *   Si on deplace le sommet dans un plan orthogonal au vecteur, le volume
 *   des phases est conserve a l'ordre 1 par rapport a l'amplitude du
 *   deplacement.
 *
 * @param (maillage)
 * @param (differentielle_volume)
 * @return (int (1))
 */
int Remaillage_FT::calculer_differentielle_volume(
  const Maillage_FT_Disc& maillage,
  DoubleTab& differentielle_volume) const
{
  const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();
  const int nb_sommets = maillage.nb_sommets();
  const int nb_facettes = maillage.nb_facettes();
  const IntTab& facettes = maillage.facettes();
  const int nb_som_par_facette = facettes.dimension(1);
  int isom, fa7, k;

  const int dim = Objet_U::dimension;
  const double angle_bidim_axi = Maillage_FT_Disc::angle_bidim_axi();

  differentielle_volume.resize(nb_sommets, dim);
  differentielle_volume = 0.;

  const double facteur = 1. / (double) dim;

  // on balaye les facettes reelles : la differentielle du volume
  // par rapport au deplacement d'un sommet est la somme des
  // differentielles de volume engendre par chaque facette voisine
  // du sommet. Pour une facette la differentielle de volume
  // est :
  if (!bidim_axi)
    {
      const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
      //  normale_unitaire * surface / nb_sommets_par_facette
      double normale[3] = { 0., 0., 0. };
      for (fa7 = 0; fa7 < nb_facettes; fa7++)
        {
          if (maillage.facette_virtuelle(fa7))
            continue;
          //on balaye ses sommets
          double surf = surface_facettes[fa7];
          double x = surf * facteur;
          for (k = 0; k < dim; k++)
            normale[k] = normale_facettes(fa7, k);
          // Ajout de la contribution de la facette a chacun des trois sommets.
          for (isom=0 ; isom<nb_som_par_facette ; isom++)
            {
              int som = facettes(fa7,isom);
              for (k = 0; k < dim; k++)
                differentielle_volume(som,k) += x * normale[k];
            }
        }
    }
  else
    {
      const double un_tiers = 1. / 3.;
      const double un_sixieme = 1. / 6.;
      const DoubleTab& sommets = maillage.sommets();
      // Cas bidim_axi
      for (int facette = 0; facette < nb_facettes; facette++)
        {
          if (maillage.facette_virtuelle(facette))
            continue;
          const int s1 = facettes(facette, 0);
          const int s2 = facettes(facette, 1);
          const double r1 = sommets(s1, 0);
          const double y1 = sommets(s1, 1);
          const double r2 = sommets(s2, 0);
          const double y2 = sommets(s2, 1);
          const double L2 = (r2-r1)*(r2-r1) + (y2-y1)*(y2-y1);
          const double L  = sqrt(L2);
          const double nx = normale_facettes(facette, 0);
          const double ny = normale_facettes(facette, 1);
          const double dv_dn1 = L * (r1 * un_tiers + r2 * un_sixieme) * angle_bidim_axi;
          const double dv_dn2 = L * (r2 * un_tiers + r1 * un_sixieme) * angle_bidim_axi;
          differentielle_volume(s1, 0) += nx * dv_dn1;
          differentielle_volume(s1, 1) += ny * dv_dn1;
          differentielle_volume(s2, 0) += nx * dv_dn2;
          differentielle_volume(s2, 1) += ny * dv_dn2;
        }
    }
  // Collection des donnees des sommets virtuels
  const Desc_Structure_FT& desc = maillage.desc_sommets();
  desc.collecter_espace_virtuel(differentielle_volume, MD_Vector_tools::EV_SOMME);
  desc.echange_espace_virtuel(differentielle_volume);

  return 1;
}

/*! @brief Cette fonction calcule la difference de volume au niveau d'une facette par rapport a une position initiale des sommets
 *
 *   version 2D
 *
 * @param (fa7) indice de la facette de calcul
 * @param (maillage) maillage a barycentrer
 * @param (position_initiale) position initiale des sommets
 * @return (double) la variation de volume
 */
double Remaillage_FT::calculer_variation_volume_facette_2D(int fa7, const Maillage_FT_Disc& maillage,
                                                           const DoubleTab& position_initiale) const
{
  const ArrOfInt& pe_owner = maillage.sommet_PE_owner();
  const IntTab& facettes = maillage.facettes();
  const DoubleTab& sommets = maillage.sommets();

#if DEBUG_CONSERV_VOLUME
#if VERBEUX
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const DoubleTab& normales_facettes = maillage.get_update_normale_facettes();
#endif
#endif
  double coord_som0[2], coord_som1[2], coord_som_opp[2];
  const int som0 = facettes(fa7,0);
  const int som1 = facettes(fa7,1);

  const double angle_bidim_axi =  bidim_axi ? Maillage_FT_Disc::angle_bidim_axi() : 0.;
  const double un_tiers = 1. / 3.;

  if (som0==som1)
    {
      //facette de surface nulle : sort
      return 0.;
    }

  //on decompose la variation de volume en 2 triangles, lies aux deplacements successifs des sommets
  int ordre = 0;
  if ( (pe_owner[som0]==pe_owner[som1] && som0>som1) || (pe_owner[som0]>pe_owner[som1]) )
    {
      //l'ordre est inverse si :
      // - les 2 sommets sont sur le meme proc, mais l'indice du sommet1 est plus petit que celui du sommet0
      // - ou si le numero du proc du sommet1 est plus petit que celui du proc du sommet0
      ordre = 1;
    }
  //premier triangle
  coord_som0[0] = position_initiale(som0,0);
  coord_som0[1] = position_initiale(som0,1);
  coord_som1[0] = position_initiale(som1,0);
  coord_som1[1] = position_initiale(som1,1);
  if (ordre==0)
    {
      coord_som_opp[0] = sommets(som0,0);
      coord_som_opp[1] = sommets(som0,1);
    }
  else
    {
      coord_som_opp[0] = sommets(som1,0);
      coord_som_opp[1] = sommets(som1,1);
    }
  double v1 = FTd_calculer_aire_triangle(coord_som0,coord_som1,coord_som_opp);
  if (bidim_axi)
    // x_G : x coordinate of the centre of gravity is located at a third of the position of all vertex.
    // Guldin's Theorem :
    v1 *= (coord_som0[0] + coord_som1[0] + coord_som_opp[0]) * un_tiers * angle_bidim_axi;

  //second triangle
  if (ordre==0)
    {
      coord_som0[0] = sommets(som0,0);
      coord_som0[1] = sommets(som0,1);
      coord_som1[0] = position_initiale(som1,0);
      coord_som1[1] = position_initiale(som1,1);
      coord_som_opp[0] = sommets(som1,0);
      coord_som_opp[1] = sommets(som1,1);
    }
  else
    {
      coord_som0[0] = position_initiale(som0,0);
      coord_som0[1] = position_initiale(som0,1);
      coord_som1[0] = sommets(som1,0);
      coord_som1[1] = sommets(som1,1);
      coord_som_opp[0] = sommets(som0,0);
      coord_som_opp[1] = sommets(som0,1);
    }
  double v2 = FTd_calculer_aire_triangle(coord_som0,coord_som1,coord_som_opp);
  if (bidim_axi)
    v2 *= (coord_som0[0] + coord_som1[0] + coord_som_opp[0]) * un_tiers * angle_bidim_axi;
#if DEBUG_CONSERV_VOLUME
#if VERBEUX
// On calcule pour comparaison la variation de volume dv par difference simple entre les
// volumes initiaux et finaux (avant vs apres mouvement).
// Je n'ai pas verifie si la somme est la meme in fine, sur le maillage.
// En tout cas, sur un seul segment, les deux methodes different car lorsque la projection
// sur dy du segment change, le volume genere par cette facette change, au detriment de sa voisine.
// Mais in fine, ca n'a pas d'effet au global.
// La methode standard, originelle, ne considere pas le mouvement d'un sommet selon la direction de la
// facette comme generatrice de variation de volume. C'est mieux ainsi, de voir le transfert d'une facette
// a l'autre comme globalement conservatif, donc pas comme une variation de volume.
// Donc v1+v2 c'est mieux :D
  if (bidim_axi)
    {
      const double un_sixieme = 1. / 6.;
      double dy = position_initiale(som0,1)-position_initiale(som1,1);
      double x0 = position_initiale(som0,0);
      double x1 = position_initiale(som1,0);
      // On n'a pas la vieille normale et la vieille surface pour calculer dy par s*nx
      // ON S'APPUIE SUR ORDRE en esperant que ca marche.
      const double vp1 = angle_bidim_axi * un_sixieme *dy*(x1*x1+x0*x0+x1*x0);

      // Apres deplacement :
      dy = sommets(som0,1)-sommets(som1,1);
      x0 = sommets(som0,0);
      x1 = sommets(som1,0);

      //const double vp2 = angle_bidim_axi * un_sixieme *s*normale_scalaire_direction_x*(x1*x1+x0*x0+x1*x0);
      const double vp2 = angle_bidim_axi * un_sixieme *dy*(x1*x1+x0*x0+x1*x0);

      // Si tout va bien, on devrait avoir le meme signe :
      const double s = surface_facettes[fa7];
      const double normale_scalaire_direction_x = normales_facettes(fa7, 0);
      Cerr << "dy=" << dy << " s.nx " << s*normale_scalaire_direction_x << finl;
      assert(normale_scalaire_direction_x*dy>=0.);

      const double dv = vp2 - vp1;

      Cerr << "Remaillage_FT::calculer_variation_volume_facette_2D" << finl;
      Cerr << "Positions ini(x;y) -> fin(x;y) "
           << "som0(" << position_initiale(som0,0)<< ";"<< position_initiale(som0,1)<< ") "
           "som1(" << position_initiale(som1,0)<< ";"<< position_initiale(som1,1)<< ") "
           << " -> "
           "som0(" << sommets(som0,0)<< ";"<< sommets(som0,1)<< ") "
           "som1(" << sommets(som1,0)<< ";"<< sommets(som1,1)<< ") "
           << finl;
      Cerr << "v1= " << v1 << finl;
      Cerr << "v2= " << v2 << finl;
      Cerr << "v1+v2= " << v1+v2 << finl;
      Cerr << "vp1= " << vp1 << finl;
      Cerr << "vp2= " << vp2 << finl;
      Cerr << "dv= " << dv << finl;
      return v1+v2; // quand meme !
      //return dv;
    }
#endif
#endif
  return v1 + v2;
}

/*! @brief Cette fonction calcule la difference de volume au niveau d'une facette par rapport a une position initiale des sommets
 *
 *     Methode corrigee le 11 juillet 2013 par B.M: dans la version precedente, l'implementation
 *      du tri des trois sommets en fonction de l'indice global des sommets etait faux
 *      d'ou de petites et rares erreurs de conservation du volume en parallele.
 *   version 3D
 *
 * @param (fa7) indice de la facette de calcul
 * @param (maillage) maillage a barycentrer
 * @param (position_initiale) position initiale des sommets
 * @return (double) la variation de volume
 */
static inline double calculer_volume_facette_3D_avec_ordre(const DoubleTab& position_initiale,
                                                           const DoubleTab& position_finale,
                                                           int facette[3],
                                                           int ordre_sommets[3])
{
  double v = 0.;
  // On va couper le prisme en 3 tetraedres.
  // 4 sommets du tetraedre courant:
  FTd_vecteur3 sommets_base[3], sommet4;
  // Initialisation des sommets de la base avec la facette en position initiale:
  {
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        sommets_base[i][j] = position_initiale(facette[i], j);
  }
  for (int tetra = 0; tetra < 3; tetra++)
    {
      // Indice du 4e sommet dans la facette courante (sur le maillage deplace)
      const int indice_sommet4 = ordre_sommets[tetra];
      // Position du 4e sommet:
      {
        for (int i = 0; i < 3; i++)
          sommet4[i] = position_finale(facette[indice_sommet4], i);
      }
      v += FTd_calculer_volume_tetraedre(sommets_base[0], sommets_base[1], sommets_base[2], sommet4);
      {
        for (int i = 0; i < 3; i++)
          sommets_base[indice_sommet4][i] = sommet4[i];
      }
    }
  return v;
}

double Remaillage_FT::calculer_variation_volume_facette_3D(int fa7, const Maillage_FT_Disc& maillage,
                                                           const DoubleTab& position_initiale) const
{
  double varVolume = 0.;
  int facette[3];
  {
    const IntTab& facettes = maillage.facettes();
    for (int i = 0; i < 3; i++)
      facette[i] = facettes(fa7, i);
  }
  const DoubleTab& position_finale = maillage.sommets();
#ifdef ALGO_NON_PARALLELE
  int ordre_sommets[3] = { 0, 1, 2 };
  // Calcul de l'ordre dans lequel on va construire les 3 tetraedres (ordre croissant des indices globaux
  // des sommets (indice global = numero du pe_owner, numero sommet sur le pe_owner)
  // Ce tri est important pour que le volume construit en extrudant le triangle de l'ancienne position a la
  // nouvelle, approxime par des tetraedres, soit jointif avec les volumes construits pour les triangles
  // voisins.
  // En triant, on s'assure que pour chaque arete du maillage en triangle, le premier sommet choisi comme sommet4
  // est le meme pour les deux triangles adjacents a l'arete (c'est le plus petit en indice global).
  {
    long long indice_global[3];
    for (int i = 0; i < 3; i++)
      {
        int isom = facette[i];
        indice_global[i] = ((long long) maillage.sommet_PE_owner_[isom]) << 32;
        indice_global[i] += maillage.sommet_num_owner_[isom];
      }
    // Tri dans l'ordre croissant des indices_globaux (tri a bulles immediat)
#define check_swap(i,j) \
if (indice_global[ordre_sommets[i]] > indice_global[ordre_sommets[j]]) \
  { int x = ordre_sommets[i]; ordre_sommets[i] = ordre_sommets[j]; ordre_sommets[j] = x; }
    check_swap(0,1);
    check_swap(1,2);
    check_swap(0,1);
#undef check_swap
  }
  varVolume = calculer_volume_facette_3D_avec_ordre(position_initiale, position_finale, facette, ordre_sommets);
#else
  {
    // Calcul le volume avec toutes les permutations possibles de l'ordre des sommets, devient independant de l'ordre
    // des sommets du maillage, donc du decoupage
    int ordre_sommets[6][3] =
    {
      { 0, 1, 2 },
      { 0, 2, 1 },
      { 1, 2, 0 },
      { 1, 0, 2 },
      { 2, 0, 1 },
      { 2, 1, 0 }
    };
    for (int i = 0; i < 6; i++)
      varVolume += calculer_volume_facette_3D_avec_ordre(position_initiale, position_finale, facette, ordre_sommets[i]);
    varVolume /= 6.;
  }
#endif
  return varVolume;
}

/*! @brief Cette fonction calcule le volume de phase 0 engendre par chaque sommet lors du deplacement de l'interface entre la position initiale et la position actuelle
 *
 *   du maillage (le volume engendre est positif si l'interface se deplace dans la
 *   direction de la normale aux facettes).
 *   On definit d'abord le volume engendre par le deplacement de chaque facette :
 *    C'est le volume d'un polyedre construit comme reunion de trois tetraedres,
 *    engendres par deplacement successifs des trois sommets
 *    dans un ordre conventionnel (deux sommets d'une arete sont toujours deplaces
 *    dans le meme ordre pour que les deux triangles voisins engendrent des
 *    tetraedres conformes, sans laisser de trous et sans se chevaucher)
 *   Ensuite, le volume engendre par chaque facette est divise en trois parts egales
 *   et reparti sur les trois sommets de la facette.
 *
 * @param (maillage) maillage
 * @param (position_initiale) position initiale des sommets (doit avoir la meme taille que maillage.sommets(), et l'espace virtuel doit etre a jour)
 * @param (varVolume) la variation de volume pour chaque sommet (on lui donne la bonne taille et on met a jour l'espace virtuel)
 * @return (double) Variation de volume totale sur l'ensemble du domaine.
 */
double Remaillage_FT::calculer_variation_volume(const Maillage_FT_Disc& maillage,
                                                const DoubleTab& position_initiale,
                                                ArrOfDouble& varVolume) const
{
  double dvolume_total = 0.;

  const int nb_facettes = maillage.nb_facettes();
  const int nb_sommets = maillage.nb_sommets();
  const IntTab& facettes = maillage.facettes();
  const int nb_som_par_facette = facettes.dimension(1);
  const double inv_som_facette = 1. / (double) nb_som_par_facette;

  varVolume.resize_array(nb_sommets);
  varVolume = 0.;

  int fa7;
  const int dimension3 = (Objet_U::dimension==3);

  for (fa7 = 0; fa7 < nb_facettes; fa7++)
    {
      // Traitement des facettes reelles uniquement:
      if (maillage.facette_virtuelle(fa7))
        continue;

      double dv;
      if (dimension3)
        dv = calculer_variation_volume_facette_3D(fa7,maillage,position_initiale);
      else
        dv = calculer_variation_volume_facette_2D(fa7,maillage,position_initiale);
      dvolume_total += dv;
      // Redistributes the volume variation to the ndim vertices (2 or 3).
      // Some volume will end on virtual vertices; hence, we'll have to "collect" contribution in the end.
      dv *= inv_som_facette;
      int isom;
      for (isom = 0 ; isom < nb_som_par_facette ; isom++)
        {
          const int sommet = facettes(fa7, isom);
          varVolume[sommet] += dv;
        }
    }
  // Collecting varVolume contributions on virtual space.
  // and finally add it to the real contribution.
  const Desc_Structure_FT& desc = maillage.desc_sommets();
  desc.collecter_espace_virtuel(varVolume, MD_Vector_tools::EV_SOMME);
  desc.echange_espace_virtuel(varVolume);

  dvolume_total = mp_sum(dvolume_total);
  return dvolume_total;
}

/*! @brief Cette fonction calcule une correction sur un deplacement liee a une variation de volume imposee
 *  Utile pour IJK
 *
 * @param (deplacement) delacement a corriger
 * @param (varVolume) la variation de volume
 * @param (deplacement_varVolume) la normale au plan de conservation de volume
 * @param (norme2_deplacement_varVolume) carre de la normale au plan
 * @return (int) 1 si le barycentrage s'est deroule correctement, 0 sinon
 */
int Remaillage_FT::calculer_correction_deplacement(DoubleTab& deplacement,
                                                   const ArrOfDouble& varVolume,
                                                   const DoubleTab& deplacement_varVolume,
                                                   const ArrOfDouble& norme2_deplacement_varVolume) const
{
  int res = 1;
  const int nb_som = deplacement.dimension(0);
  int som,k;
  double scal, dvn;
#ifndef  NDEBUG
  if (impr_>1000)
    {
      Process::Journal()<<"Remaillage_FT::calculer_correction_deplacement  nb_som="<<nb_som<<"  variation_volume_= "<<variation_volume_;
      Process::Journal()<<"  surface_interface_= "<<surface_interface_<<"  ->depl= "<<variation_volume_/surface_interface_<<finl;
    }
#endif

  double coeff_depl,coeff_varV;
#ifndef NDEBUG
  double depl,depl_max = 0.;
#endif

  coeff_depl = 1.;
  coeff_varV = 1.;  //correction due a la variation de volume deconnectee

  for (som=0 ; som<nb_som ; som++)
    {
      if (norme2_deplacement_varVolume[som]!=0.)
        {
          //si norme2_deplacement_varVolume[som]==0., on ne modifie pas le deplacement,
          //ca peut arriver lors d'etape de remaillage (c'est sommet supprime en fait, et qui
          //n'est plus lie a aucune facette)
          scal = 0.;
          for (k=0 ; k<dimension ; k++)
            {
              scal += deplacement(som,k) * deplacement_varVolume(som,k);
            }

#ifndef NDEBUG
          if (impr_>9000)
            {
              if (varVolume[som]!=0.)
                Process::Journal()<<"correcDepl sommet "<<som<<" varVolume= "<<varVolume[som]<<" scal= "<<scal<<finl;
            }
#endif
#if 1
          // cas avec correction par sommet
          double norme = sqrt(norme2_deplacement_varVolume[som]);
          for (k=0 ; k<dimension ; k++)
            {
              dvn = deplacement_varVolume(som,k) / norme2_deplacement_varVolume[som];
              deplacement(som,k) = relax_barycentrage_ * (
                                     coeff_depl * ( deplacement(som,k) - scal * dvn ) -
                                     coeff_varV * varVolume[som] * deplacement_varVolume(som,k)/norme );
            }
#else
          //cas avec correction globale
          double norme = sqrt(norme2_deplacement_varVolume[som]);
          for (k=0 ; k<dimension ; k++)
            {
              dvn = deplacement_varVolume(som,k) / norme2_deplacement_varVolume[som];
              deplacement(som,k) = relax_barycentrage_ * (
                                     coeff_depl * ( deplacement(som,k) - scal * dvn )  -
                                     coeff_varV * variation_volume_/surface_interface_ * deplacement_varVolume(som,k)/norme );
            }
#endif
#ifndef NDEBUG
          if (impr_>5000)
            {
              depl = 0.;
              for (k=0 ; k<dimension ; k++)
                {
                  depl += deplacement(som,k) * deplacement(som,k);
                }
              depl = sqrt(depl);
              if (depl>depl_max)
                {
                  depl_max = depl;
                }
            }
#endif
        }
    }

#ifndef NDEBUG
  if (impr_>5000)
    {
      Process::Journal()<<" variation_volume_depl_max="<<depl_max<<finl;
    }
#endif

  return res;
}

/*! @brief Cette fonction calcule pour chaque sommet le barycentre de l'ensemble des facettes voisines du sommet.
 *
 * Les "dimension" premieres colonnes
 *      contiennent les coordonnees du barycentre, la colonne "dimension"
 *      contient la somme des surfaces des facettes voisines du sommet.
 *
 * @param (maillage) maillage a barycentrer
 * @param (barycentres) par sommet, barycentres de ses facettes voisines
 * @return (int) toujours 1...
 */
int Remaillage_FT::calculer_barycentre_facettes_voisines(const Maillage_FT_Disc& maillage,
                                                         DoubleTab& barycentres) const
{
  int res = 1;
  const int nb_facettes = maillage.nb_facettes();
  const int nb_sommets = maillage.nb_sommets();
  const IntTab& facettes = maillage.facettes();
  const DoubleTab& sommets = maillage.sommets();
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();

  // colonnes 1..dimension - 1 : coordonnees du barycentre
  // colonne dimension : surface des facettes voisines
  barycentres.resize(nb_sommets, dimension+1);
  barycentres = 0.;

  const int dim = Objet_U::dimension;
  const double inv_dim = 1. / (double) dim;

  int fa7, som, isom, k;
  for (fa7=0 ; fa7<nb_facettes ; fa7++)
    {
      if (maillage.facette_virtuelle(fa7))
        continue;

      double bary[3] = {0., 0., 0.};
      for (isom = 0; isom < dim; isom++)
        {
          const int sommet = facettes(fa7, isom);
          for (k = 0; k < dim; k++)
            bary[k] += sommets(sommet,k);
        }
      const double surface = surface_facettes[fa7];
      const double surface_inv_dim = surface * inv_dim;
      for (isom = 0; isom < dim; isom++)
        {
          const int sommet = facettes(fa7, isom);
          for (k = 0; k < dim; k++)
            barycentres(sommet, k) += bary[k] * surface_inv_dim;
          barycentres(sommet, dim) += surface;
        }
    }

  //sommation des surf*barycentre et surfaces des fa7 voisines
  const Desc_Structure_FT& desc = maillage.desc_sommets();
  desc.collecter_espace_virtuel(barycentres, MD_Vector_tools::EV_SOMME);

  //normalisation par la surface totale des facettes voisines
  for (som=0 ; som<nb_sommets ; som++)
    {
      if (maillage.sommet_virtuel(som))
        continue;

      const double surface = barycentres(som, dim);
      if (surface <= 0.)
        {
          // Surface nulle (pas impossible), le barycentre est le sommet
          for (k=0 ; k<dim ; k++)
            {
              barycentres(som,k) = sommets(som,k);
            }
        }
      else
        {
          const double inverse_surface = 1. / surface;
          for (k=0 ; k<dimension ; k++)
            barycentres(som,k) *= inverse_surface;
        }
    }
  desc.echange_espace_virtuel(barycentres);

  return res;
}

#if DEBUG_CONSERV_VOLUME
// Beware, this algorithm is only efficient
// IF the interfaces are closed in the direction considered
//    (for the normal as well as for the coordinate)!
// OR IF the opening face is at the value zero of the coordinate?
double Remaillage_FT::calculer_volume_mesh(const Maillage_FT_Disc& mesh) const
{
  const int n = mesh.nb_facettes();
  double volume = 0.;
  const ArrOfDouble& surfaces_facettes = mesh.get_update_surface_facettes();
  const DoubleTab& normales_facettes = mesh.get_update_normale_facettes();
  const IntTab& facettes = mesh.facettes();
  const DoubleTab& sommets = mesh.sommets();
  const int DIR_projection = 0; // On projette sur (0->x ou 1->y ou 2->z if DIM=3)
  for (int i = 0; i < n; i++)
    {
      if (mesh.facette_virtuelle(i))
        continue;
      const double s = surfaces_facettes[i];
      const double normale_scalaire_direction = normales_facettes(i, DIR_projection);
      // Coordonnee du centre de gravite de la facette
      const int i0 = facettes(i,0);
      const int i1 = facettes(i,1);
      double coord_centre_gravite_DIR = sommets(i0,DIR_projection) + sommets(i1,DIR_projection) ;
      if (Objet_U::dimension==3)
        {
          const int i2 = facettes(i,2);
          coord_centre_gravite_DIR += sommets(i2,DIR_projection);
        }
      coord_centre_gravite_DIR /= (float) Objet_U::dimension;
      //const double coord_centre_gravite_i = (sommets(i0,0) + sommets(i1,0) + sommets(i2,0)) / 3.;
      //const double coord_centre_gravite_k = (sommets(i0,2) + sommets(i1,2) + sommets(i2,2)) / 3.;
      double volume_prisme=0;
      if (Objet_U::bidim_axi)
        {
          const double angle_bidim_axi =  Objet_U::bidim_axi ? Maillage_FT_Disc::angle_bidim_axi() : 0.;
          const double x0 = sommets(i0,0);
          const double x1 = sommets(i1,0);
          // The volume can actually be negative when normale_scalaire_direction is negative.
          volume_prisme = angle_bidim_axi / 6. * (s * normale_scalaire_direction) * (x0*x0+x1*x1+x0*x1);
        }
      else
        {
          volume_prisme = coord_centre_gravite_DIR * s * normale_scalaire_direction; // coord_centre_gravite_DIR, where dir should be the
        }
      volume += volume_prisme;
    }
  volume = Process::mp_sum(volume);
  return volume;
}

double Remaillage_FT::calculer_somme_dvolume(const Maillage_FT_Disc& mesh, const ArrOfDouble& dvolume) const
{

  const int n = dvolume.size_array();
  assert(n == mesh.nb_sommets());
  double somme = 0.;
  for (int i = 0; i < n; i++)
    {
      if (!mesh.sommet_virtuel(i))
        {
          somme += dvolume[i];
        }
    }
  somme = Process::mp_sum(somme);
  return somme;
}
#endif

/*! @brief Deplace les sommets du maillage pour les redistribuer de facon homogene.
 *
 * Le deplacement est la somme d'une composante tangentielle et d'une composante
 *   normale a l'interface. La composante tangentielle ramene le sommet vers le
 *   barycentre des facettes voisines du sommet. La composante normale est calculee
 *   de sorte a produire la variation de volume imposee.
 *   Traitement des lignes de contact: le deplacement est decompose en une composante
 *   normale a la ligne de contact et une composante tangente.
 *
 * @param (maillage) le maillage a redistribuer (il retourne dans l'etat minimal)
 * @param (relaxation_direction_tangente) si ce parametre vaut 1, le sommet est deplace sur le barycentre. le parametre fixe la fraction du deplacement a realiser.
 * @param (relaxation_direction_normale)
 */
double Remaillage_FT::redistribuer_sommets(Maillage_FT_Disc&   maillage,
                                           const double        relaxation_direction_tangente,
                                           const double        relaxation_direction_normale,
                                           ArrOfDouble& var_volume_impose,
                                           ArrOfDouble& var_volume_obtenu) const
{
  assert(relaxation_direction_tangente >= 0. && relaxation_direction_tangente <= 1.);
  assert(relaxation_direction_normale >= 0. && relaxation_direction_normale <= 1.);
  assert(var_volume_impose.size_array() == maillage.sommets().dimension(0));

  const int nb_som = maillage.nb_sommets();
  const int dim = Objet_U::dimension;

  DoubleTabFT differentielle_volume;
  DoubleTabFT barycentres;
  DebogFT::verifier_maillage_ft("Remaillage_FT::redistribuer_sommets maillage initial", maillage);
  DebogFT::verifier_tableau_sommets("Remaillage_FT::redistribuer_sommets var_volume_impose", maillage, var_volume_impose);
  calculer_differentielle_volume       (maillage, differentielle_volume);
  DebogFT::verifier_tableau_sommets("Remaillage_FT::redistribuer_sommets differentielle volume", maillage, differentielle_volume);
  calculer_barycentre_facettes_voisines(maillage, barycentres);
  DebogFT::verifier_tableau_sommets("Remaillage_FT::redistribuer_sommets barycentres", maillage, barycentres);


  // Copie du tableau des sommets pour calculer la variation de volume a posteriori:
  DoubleTabFT sommets(maillage.sommets());
  DoubleTabFT deplacement(nb_som, dim);

  // Estimation du deplacement a realiser
  //  Dans la direction tangente, on va vers le barycentre des facettes voisines.
  //  Dans la direction normale, on va pour obtenir la variation de volume impose.
  // Le vecteur normal (direction de deplacement pour faire varier le volume)
  double normale[3] = { 0., 0., 0. };
  // Differentielle du volume par rapport au deplacement
  double diff[3] = {0., 0., 0.};
  // Le deplacement tangent (direction de deplacement qui conserve le volume)
  double dtangent[3] = { 0., 0., 0. };
  double depl[3] = { 0., 0., 0. };
  int sommet, k;
  int clipped_move = 0; // Compteur de deplacements tronques
  const int dim3 = (Objet_U::dimension==3);
  const Parcours_interface& parcours = maillage.refparcours_interface_.valeur();

  for (sommet = 0; sommet < nb_som; sommet++)
    {
      // Calcul des sommets reels uniquement :
      if (maillage.sommet_virtuel(sommet))
        continue;
      // Calcul de la direction du deplacement normal:
      //  Pour une ligne de contact c'est la projection de la differentielle
      //  du volume sur le bord (composante normale a la ligne de contact)
      //  Sinon c'est la differentielle du volume
      normale[0] = diff[0] = differentielle_volume(sommet, 0);
      normale[1] = diff[1] = differentielle_volume(sommet, 1);
      if (dim3)
        normale[2] = diff[2] = differentielle_volume(sommet, 2);
      const int face_bord = maillage.sommet_face_bord_[sommet];
      if (face_bord >= 0)
        parcours.projeter_vecteur_sur_face(face_bord, normale[0], normale[1], normale[2]);

      // Normaliser le vecteur "normale".
      double norme2_n = normale[0]*normale[0] + normale[1] * normale[1] + normale[2]*normale[2];
      if (norme2_n == 0.)
        {
          // Le maillage est mal foutu, on va vers le barycentre sans se
          // preoccuper de conserver le volume.
        }
      else
        {
          double x = 1. / sqrt(norme2_n);
          normale[0] *= x;
          normale[1] *= x;
          normale[2] *= x;
        }

      // Produit scalaire normale * differentielle:
      // (is actually equal to the norm of the differential IF we're not on a contact line)
      const double normale_scal_diff = normale[0] * diff[0] + normale[1] * diff[1] + normale[2] * diff[2];
      const double i_normale_scal_diff = (normale_scal_diff == 0.) ? 0. : 1./normale_scal_diff;
      // Produit scalaire de deplacement tangent avec la normale:
      double dtangent_scalaire_n = 0.;
      for (k = 0; k < dim; k++)
        {
          double x = barycentres(sommet,k) - sommets(sommet, k);
          dtangent[k] = x;
          dtangent_scalaire_n += x * normale[k];
        }
      const double a = dtangent_scalaire_n;
      // Amplitude du deplacement normal a realiser pour obtenir la variation de volume impose
      const double b = var_volume_impose[sommet] * i_normale_scal_diff;
      double norme2_deplacement = 0.;
      for (k = 0; k < dim; k++)
        {
          // Deplacement dans la direction tangente (projection du vecteur
          // "sommets -> barycentres" sur le plan orthogonal a la normale :
          double dt = (dtangent[k] - normale[k] * a) * relaxation_direction_tangente;
          // Deplacement dans la direction normale (tel que la variation
          // du volume soit egale a la variation imposee)
          double dn = normale[k] * b * relaxation_direction_normale;
          double x = dt + dn;
          depl[k] = x;
          norme2_deplacement += x * x;
        }
      {
        // Si l'amplitude du deplacement est superieure a la taille caracteristique
        // des facettes, alors on tronque le deplacement.
        double facteur = 1.;
        double surface_totale_sommet = barycentres(sommet, dim);
        if (norme2_deplacement > barycentres(sommet, dim))
          {
            clipped_move++;
            facteur = sqrt(surface_totale_sommet / norme2_deplacement);
          }
        for (k = 0; k < dim; k++)
          deplacement(sommet,k) = depl[k] * facteur;
      }
    }
  if (clipped_move > 0)
    {
      Journal() << "Remaillage_FT::redistribuer_sommets clipped_move "
                << clipped_move << " / " << nb_som << finl;
    }
  DebogFT::verifier_tableau_sommets("Remaillage_FT::redistribuer_sommets deplacement", maillage, deplacement);
  // Deplacement des sommets du maillage.
  const Desc_Structure_FT& desc = maillage.desc_sommets();
  desc.echange_espace_virtuel(deplacement);
  maillage.preparer_tableau_avant_transport(var_volume_impose, desc);
  maillage.preparer_tableau_avant_transport(sommets, desc);

#if DEBUG_CONSERV_VOLUME
  double volume_init = calculer_volume_mesh(maillage);
  double dvol_init   = calculer_somme_dvolume(maillage, var_volume_impose);
#endif
  maillage.transporter(deplacement);
#if DEBUG_CONSERV_VOLUME
  double volume_apres = calculer_volume_mesh(maillage);
#endif
  const int new_nb_som = maillage.nb_sommets();
  maillage.update_tableau_apres_transport(var_volume_impose, new_nb_som, desc);
  maillage.update_tableau_apres_transport(sommets, new_nb_som, desc);

  DebogFT::verifier_maillage_ft("Remaillage_FT::redistribuer_sommets maillage final", maillage);
  // Calcul de la variation de volume obtenue :
  double dvolume_tot = calculer_variation_volume(maillage, sommets, var_volume_obtenu);
#if DEBUG_CONSERV_VOLUME
  double  volume_apres2 = calculer_volume_mesh(maillage);
  double dvol_obtenu = calculer_somme_dvolume(maillage, var_volume_obtenu);
  Cerr << "Remaillage_FT::redistribuer_sommets BILAN "
       << " volume avt= " << volume_init
       << " apres= " << volume_apres
       << " apres2= " << volume_apres2
       << " dvolume demande= " << dvol_init
       << " obtenu= " << dvol_obtenu << finl;
#endif
  return dvolume_tot;
}

/*! @brief deplacement des sommets se sorte a produire la variation de volume prescrite a chaque sommet.
 *
 * Precondition: pas de facettes virtuelles
 */
void Remaillage_FT::corriger_volume(Maillage_FT_Disc& maillage, ArrOfDouble& var_volume)
{
  corriger_volume_(maillage, var_volume, nb_iter_bary_volume_seul_);
}
void Remaillage_FT::corriger_volume_(Maillage_FT_Disc& maillage, ArrOfDouble& var_volume, const int nb_iter_corrections_vol)
{
  regulariser_maillage(maillage, var_volume,
                       0., /* pas de deplacement tangent des sommets */
                       0.,
                       0, /* pas de barycentrage */
                       0, /* pas de lissage */
                       nb_iter_corrections_vol,
                       seuil_dvolume_residuel_);
  supprimer_facettes_bord(maillage);
  nettoyer_maillage(maillage);
}

/*! @brief applique barycentrage, lissage et correction de volume.
 *
 * On applique le nombre d'iterations de lissage systematique.
 *
 * Precondition: pas de facettes virtuelles
 */
void Remaillage_FT::barycentrer_lisser_systematique(double temps, Maillage_FT_Disc& maillage)
{
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "Barycentrer_lisser_sys");
  statistiques().begin_count(stat_counter);
  if (Process::je_suis_maitre())
    Journal() << "barycentrer_lisser_systematique" << finl;
  temps_dernier_lissage_ = temps;
  ArrOfDoubleFT var_volume(maillage.nb_sommets());
  var_volume = 0.;
  regulariser_maillage(maillage,
                       var_volume,
                       relax_barycentrage_,
                       lissage_courbure_coeff_,
                       nb_iter_barycentrage_,
                       lissage_courbure_iterations_systematique_,
                       nb_iter_bary_volume_seul_,
                       seuil_dvolume_residuel_);
  supprimer_facettes_bord(maillage);
  nettoyer_maillage(maillage);
  statistiques().end_count(stat_counter);
}

/*! @brief idem mais avec le nombre d'iterations de lissage si remaillage
 *
 * Precondition: pas d'elements virtuels (on doit avoir appele nettoyer_elements_virtuels())
 *
 */
void Remaillage_FT::barycentrer_lisser_apres_remaillage(Maillage_FT_Disc& maillage, ArrOfDouble& var_volume)
{
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "Barycentrer_lisser_apres_rem");
  statistiques().begin_count(stat_counter);
  if (Process::je_suis_maitre())
    Journal() << "barycentrer_lisser_apres_remaillage" << finl;
  regulariser_maillage(maillage, var_volume,
                       relax_barycentrage_,
                       lissage_courbure_coeff_,
                       nb_iter_barycentrage_,
                       lissage_courbure_iterations_si_remaillage_,
                       nb_iter_bary_volume_seul_,
                       seuil_dvolume_residuel_);
  supprimer_facettes_bord(maillage);
  nettoyer_maillage(maillage);
  statistiques().end_count(stat_counter);
}

/*! @brief Algorithme general de lissage du maillage.
 *
 * Permet de barycentrer, regulariser la courbure, et d'appliquer une correction de volume en deplacant
 *   les noeuds dans la direction normale.
 *
 * @param (maillage) maillage a barycentrer
 * @param (var_volume) pour chaque sommet du maillage initial, variation de volume a obtenir lors du deplacement. En retour, on y met le residu (difference entre var_volume initial et var_volume obtenu lors du deplacement).
 * @return (double) erreur sur la variation totale de volume obtenue (en m3)
 */
double Remaillage_FT::regulariser_maillage(Maillage_FT_Disc& maillage,
                                           ArrOfDouble& var_volume,
                                           const double facteur_barycentrage_tangent,
                                           const double coeff_lissage,
                                           const int nb_iter_barycentrage,
                                           const int nb_iter_lissage,
                                           const int max_nb_iter_correction_volume,
                                           const double seuil_dvolume) const
{
  ArrOfDoubleFT var_volume_obtenu;
  // on boucle les barycentrages tant que les sommets sont deplaces
  double dvolume = 0.;
  int iteration = 0;
  int iteration_correction_volume = 0;
  for(iteration = 0; 1; iteration++)
    {
#if DEBUG_CONSERV_VOLUME
      double volume_initial=0, volume_final=0, dvol_init=0, dvol_redistributed=0;
      double dvol_apres_regulariser=0, dvol_obtenu=0, dvol_futur=0;
      double volume_interm1=0, volume_interm2=0, volume_interm3=0;
      {
        double vol = calculer_volume_mesh(maillage);
        volume_initial = vol;
        double dvol = calculer_somme_dvolume(maillage, var_volume);
        dvol_init = dvol;
        Cerr.precision(16);
        Cerr << "Remaillage_FT::regulariser_maillage iteration " << iteration
             << " volume= " << vol << " dvolume= " << dvol << finl;

      }
#endif
      DebogFT::verifier_tableau_sommets("Remaillage_FT::regulariser var_vol initial", maillage, var_volume);
      if (iteration > 0)
        // Ne pas lisser la variation de volume demandee a la premiere iteration.
        // Ensuite, on repartit la correction de volume sur les noeuds voisins
        // pour eviter les singularites de deplacement des sommets:
        lisser_dvolume(maillage, var_volume, 1);
      DebogFT::verifier_tableau_sommets("Remaillage_FT::regulariser var_vol apres lisser", maillage, var_volume);

#if DEBUG_CONSERV_VOLUME
      {
        volume_interm1 = calculer_volume_mesh(maillage);
        double dvol = calculer_somme_dvolume(maillage, var_volume);
        dvol_redistributed = dvol;
        Cerr << "Remaillage_FT::regulariser_maillage apres lisser dvolume "
             <<  " dvolume= " << dvol << finl;

      }
#endif
      // A la premiere iteration, on ne regularise pas la courbure
      // (maillage risque d'etre vraiment pourri avec petites aretes partout)
      // regulariser_courbure() ne change pas le maillage, cette methode remplit
      // le tableau dv_courbure. Le deplacement est realise par redistribuer_sommets()
      if (iteration > 0 && iteration < nb_iter_lissage + 1)
        {
          ArrOfDoubleFT dv_courbure(maillage.nb_sommets());
          regulariser_courbure(maillage,
                               lissage_courbure_coeff_,
                               dv_courbure);
          const int n = var_volume.size_array();
          for (int i = 0; i < n; i++)
            var_volume[i] += dv_courbure[i];

          DebogFT::verifier_tableau_sommets("Remaillage_FT::regulariser var_vol apres ajout dv", maillage, var_volume);
#if DEBUG_CONSERV_VOLUME
          {
            volume_interm2 = calculer_volume_mesh(maillage);
            double dvol = calculer_somme_dvolume(maillage, dv_courbure);
            double dvolnew = calculer_somme_dvolume(maillage, var_volume);
            dvol_apres_regulariser = dvolnew;
            Cerr << "Remaillage_FT::regulariser_maillage apres calcul dvcourbure "
                 <<  " dv courbure= " << dvol << finl;

          }
#endif

        }

      // Faut-il barycentrer ?
      double facteur = facteur_barycentrage_tangent;
      if (iteration >= nb_iter_barycentrage)
        // non, on a fini les iterations de barycentrage
        facteur = 0.;

      dvolume =
        redistribuer_sommets(maillage,
                             facteur, // relaxation_direction_tangente
                             1.,    // relaxation_direction_normale
                             var_volume /* imposed */,
                             var_volume_obtenu /* resulting */);
#if DEBUG_CONSERV_VOLUME
      {
        volume_interm3 = calculer_volume_mesh(maillage);
        double dvol = calculer_somme_dvolume(maillage, var_volume);
        dvol_obtenu = calculer_somme_dvolume(maillage, var_volume_obtenu);
        Cerr << "Remaillage_FT::regulariser_maillage apres redistribuer_sommets " << iteration
             << " dvolume demande= " << dvol
             << " dvolume obtenu= " << dvol_obtenu << finl;
      }
#endif
      DebogFT::verifier_tableau_sommets("Remaillage_FT::regulariser var_vol_obtenu", maillage, var_volume_obtenu);

      // Calcul de la correction de volume a realiser a l'iteration suivante:
      var_volume -= var_volume_obtenu;
#if DEBUG_CONSERV_VOLUME
      {
        double vol = calculer_volume_mesh(maillage);
        volume_final = vol;
        double vvol = calculer_somme_dvolume(maillage, var_volume_obtenu);
        double dvol = calculer_somme_dvolume(maillage, var_volume);
        dvol_futur = dvol;
        Cerr << "Remaillage_FT::regulariser_maillage apres redistribuer_sommets " << iteration
             << " volume= " << vol
             << " var_vol obtenu= " << vvol
             << " dvolume= " << dvol
             << " variation calculee= " << volume_final-volume_initial << finl;

        Cerr << "Remaillage_FT::regulariser_maillage "
             << "iteration " << iteration << finl
             << "volume_initial " << volume_initial << finl
             << "volume_final " << volume_final << finl
             << "dvol_init " << dvol_init << finl
             << "dvol_redistributed " << dvol_redistributed  << finl//par lisser_dvol
             << "dvol_apres_regulariser " << dvol_apres_regulariser << finl // apres le lissage de la courbure
             << "dvol_obtenu " << dvol_obtenu  << finl                      // uniquement celui engendre par redistribuer_sommets
             << "dvol_futur " << dvol_futur << finl      // celui pour l'iteration future.
             << finl;

        Cerr << "Remaillage_FT::regulariser_maillage VOLUMES "
             << "iteration " << iteration << finl
             << "volume_initial " << volume_initial << finl
             << "volume_interm1 " << volume_interm1 << finl
             << "volume_interm2 " << volume_interm2 << finl
             << "volume_interm3 " << volume_interm3 << finl
             << "volume_final " << volume_final << finl
             << finl;

      }
#endif

      if (iteration >= nb_iter_barycentrage-1 && iteration >= nb_iter_lissage)
        {
          // On a termine le barycentrage et le lissage de courbure.
          iteration_correction_volume++;
          // On s'arrete apres avoir converge en volume ou apres avoir atteint le nombre
          // max d'iterations de correction.
          if (iteration_correction_volume > max_nb_iter_correction_volume
              || std::fabs(dvolume) < seuil_dvolume)
            break;
        }
    };

  return dvolume;
}

/*! @brief Regularise le champ scalaire "var_volume" defini aux sommets du "maillage".
 *
 * On regularise en faisant calculant une valeur moyenne par facette,
 *   puis en recalculant une moyenne aux sommets en fonction de la moyenne
 *   aux faces. Les moyennes sont ponderees par la surface des facettes...
 *   Le lissage est repete un nombre donne de fois.
 *   Cette methode est utilisee pour lisser les corrections de volume et eviter
 *   l'apparition de pointes sur les interfaces.
 *
 * @param (maillage) le support du champ a regulariser
 * @param (var_volume) le champ a regulariser (champ aux sommets du maillage)
 * @param (nb_iterations) nombre d'iterations de l'operation de regularisation.
 */
void Remaillage_FT::lisser_dvolume(const Maillage_FT_Disc& maillage,
                                   ArrOfDouble& var_volume,
                                   const int nb_iterations) const
{
  DoubleTabFT barycentres;
  calculer_barycentre_facettes_voisines(maillage, barycentres);
  const ArrOfDouble& surfaces = maillage.get_update_surface_facettes();
  const int nb_facettes = maillage.nb_facettes();
  const IntTab& facettes = maillage.facettes();
  ArrOfDoubleFT dv_facettes(nb_facettes);
  const int dim = Objet_U::dimension;
  const double inv_dim = 1. / (double) dim;
  int facette, j, iter;
  for (iter = 0; iter < nb_iterations; iter++)
    {
      // Repartition de dvolume sur les faces au prorata de la
      // surface de la face
      dv_facettes = 0.;
      for (facette = 0; facette < nb_facettes; facette++)
        {
          if (maillage.facette_virtuelle(facette))
            continue;
          const double surface = surfaces[facette];
          for (j = 0; j < dim; j++)
            {
              int sommet = facettes(facette, j);
              double surface_tot = barycentres(sommet, dim);
              if (surface_tot > 0.)
                dv_facettes[facette] += surface / surface_tot * var_volume[sommet];
            }
        }
      // Ne pas faire echange espace virtuel (inutile car on n'utilise pas la partie virtuelle)
      // Repartition de dv_facettes aux sommets
      var_volume = 0.;
      for (facette = 0; facette < nb_facettes; facette++)
        {
          if (maillage.facette_virtuelle(facette))
            continue;
          double dv = dv_facettes[facette] * inv_dim;
          for (j = 0; j < dim; j++)
            {
              int sommet = facettes(facette, j);
              var_volume[sommet] += dv;
            }
        }
      const Desc_Structure_FT& desc_s = maillage.desc_sommets();
      desc_s.collecter_espace_virtuel(var_volume, MD_Vector_tools::EV_SOMME);
      desc_s.echange_espace_virtuel(var_volume);
    }
}

/*! @brief Cette fonction permet de gerer le decollement de l'interface de la paroi Si un sommet de bord n'a pas de sommet voisin de bord
 *
 *      et que son deplacement tend a le faire rentrer dans le domaine,
 *      il est marque comme sommet interne
 *
 * @param (maillage) maillage a traiter
 * @param (deplacement) vecteur deplacement des sommets
 * @return (int) 1 si le remaillage s'est deroule correctement, 0 sinon
 */
int Remaillage_FT::traite_decollement(Maillage_FT_Disc& maillage, const DoubleTab& deplacement) const
{
  int res = 1;
  //fonction non utilisee en dimension 2
  if (dimension<3)
    {
      return 1;
    }
  const int nb_sommets = maillage.nb_sommets();
  const int nb_facettes = maillage.nb_facettes();
  const IntTab& facettes = maillage.facettes();
  const int nb_som_par_facette = facettes.dimension(1);
  ArrOfIntFT nbSomVois_bord(nb_sommets);
  nbSomVois_bord = 0;
  int fa7,som, nb_som_bord, isom;

  for (fa7=0 ; fa7<nb_facettes ; fa7++)
    {
      nb_som_bord = 0;
      for (isom=0 ; isom<nb_som_par_facette ; isom++)
        {
          som = facettes(fa7,isom);
          if (maillage.sommet_ligne_contact(som))
            {
              nb_som_bord++;
            }
        }
      if (nb_som_bord>1)
        {
          nb_som_bord--;
          for (isom=0 ; isom<nb_som_par_facette ; isom++)
            {
              som = facettes(fa7,isom);
              if (maillage.sommet_ligne_contact(som))
                {
                  nbSomVois_bord[som] += nb_som_bord;
                }
            }
        }
    }
  maillage.desc_sommets().collecter_espace_virtuel(nbSomVois_bord, MD_Vector_tools::EV_SOMME);

  const DoubleTab& sommets = maillage.sommets();
  const DoubleTab& xp = refdomaine_VF_->xp();
  double x,y,z=0., scal;
  for (som=0 ; som<nb_sommets ; som++)
    {
      const int elem = maillage.sommet_elem_[som];
      if (elem >= 0 /* sommet reel */
          && maillage.sommet_ligne_contact(som)
          && nbSomVois_bord[som]==0)
        {
          //sommet de bord, sans autre sommet voisin de bord
          //teste si le deplacement fait entrer le sommet dans le domaine
          //calcul de la normale entrante du domaine : vecteur sommet-cdg(elem)
          x = xp(elem,0) - sommets(som,0);
          y = xp(elem,1) - sommets(som,1);
          scal = x * deplacement(som,0) + y * deplacement(som,1);
          if (dimension==3)
            {
              z = xp(elem,2) - sommets(som,2);
              scal += z * deplacement(som,2);
            }
          if (scal>0.)
            {
#ifndef NDEBUG
              if (impr_>900)
                {
                  Process::Journal()<<"Remaillage_FT::traite_decollement : decollement detecte"<<finl;
                  maillage.printSom(som,Process::Journal());
                  Process::Journal()<<"  decoll= "<<deplacement(som,0)<<" "<<deplacement(som,1)<<" "<<deplacement(som,2)<<finl;
                }
#endif
              maillage.sommet_face_bord_[som] = -1;
            }
        }
    }
  maillage.desc_sommets().echange_espace_virtuel(maillage.sommet_face_bord_);

  maillage.maillage_modifie(Maillage_FT_Disc::MINIMAL);
  return res;
}
/*! @brief Cette fonction permet de gerer l'adherence de l'interface a la paroi Si une facette possede 3 sommets sur un paroi, elle est supprimee
 *
 * @param (maillage) maillage a traiter
 * @return (int) 1 si le remaillage s'est deroule correctement, 0 sinon
 */
int Remaillage_FT::traite_adherence(Maillage_FT_Disc& maillage) const
{
  int res = 1;
  if (supprimer_facettes_bord(maillage)==1)
    {
#ifndef NDEBUG
      Process::Journal()<<"Remaillage_FT::traite_adherence : adherence detectee"<<finl;
#endif
      nettoyer_maillage(maillage);
    }
  maillage.maillage_modifie(Maillage_FT_Disc::MINIMAL);
  return res;
}

inline double produit_scalaire(const FTd_vecteur3& a, const FTd_vecteur3& b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline void produit_vectoriel(const FTd_vecteur3& a, const FTd_vecteur3& b, FTd_vecteur3& resu)
{
  resu[0] = a[1]*b[2] - a[2]*b[1];
  resu[1] = a[2]*b[0] - a[0]*b[2];
  resu[2] = a[0]*b[1] - a[1]*b[0];
}

/*! @brief Cette fonction marque a supprimer les facettes ayant leurs 3 sommets de bord Marquer a supprimer = condenser les 3 sommets en un seul (le sommet 0)
 *
 *     MODIF GB 21/07/2015 : Pour le cas sloshing par exemple, si le domaine a des coins,
 *           l'interface peut contenir un element dont les 3 sommets sont sur le bord
 *           sans pour autant que la facette soit integralement sur le bord. Elle peut
 *           avoir son 3eme cote dans le domaine. Dans ce cas, on la conserve.
 *
 * @param (maillage) maillage a remailler
 * @return (int) 1 si le remaillage s'est deroule correctement, 0 sinon
 */
int Remaillage_FT::supprimer_facettes_bord(Maillage_FT_Disc& maillage) const
{

  int res = 0;
  const int nb_facettes = maillage.nb_facettes();
  IntTab& facettes = maillage.facettes_;
  const int nb_som_par_facette = facettes.dimension(1);
  // Raccourci vers les coordonnees des sommets du maillage eulerien
#define CONSERVER_FACETTES_COINS
#ifdef CONSERVER_FACETTES_COINS
  const int dim = Objet_U::dimension;
  const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();
  const DoubleTab& sommets = maillage.sommets();
  const Parcours_interface& parcours = maillage.refparcours_interface_.valeur();
#endif

  //DoubleTab xsom(3,3); // coords des sommets de la face: xs(isom_eul, direction)
  int fa7, isom, som =0, nb_bord;
  for (fa7=0 ; fa7<nb_facettes ; fa7++)
    {
      nb_bord = 0;
      for (isom=0 ; isom<nb_som_par_facette ; isom++)
        {
          som = facettes(fa7,isom);
          if (maillage.sommet_ligne_contact(som))
            {
              nb_bord++;
            }
        }
      if (nb_bord==dimension)
#ifndef CONSERVER_FACETTES_COINS
        {
          //facette de bord : supprimer
          res = 1;
#ifndef NDEBUG
          if (impr_>9000)
            {
              Process::Journal()<<"  fa7_bord -> Supp ";
              maillage.printFa7(fa7,0,Process::Journal());
            }
#endif

          for (isom=1 ; isom<nb_som_par_facette ; isom++)
            {
              facettes(fa7,isom) = facettes(fa7,0);
            }
        }
#else
        {
          //facette de bord : a supprimer ou pas? Pas si simple...
          // Pour le savoir, on regarde la normale a chacune des faces tour a tour.
          for (isom=0 ; isom<nb_som_par_facette ; isom++)
            {
              som = facettes(fa7,isom);
              const int face = maillage.sommet_face_bord_[som];
              FTd_vecteur3 v4= {0.,0.,0.};
              double s2 = 0.;
              if (dim==3) s2=sommets(som,2);
              parcours.calculer_normale_face_bord(face,
                                                  sommets(som,0), sommets(som,1),  s2,
                                                  v4[0], v4[1], v4[2]);
              double ps=0.;
              for (int direction = 0; direction < dim; direction++)
                {
                  ps +=  v4[direction]*normale_facettes(fa7,direction);
                }
              const double tol = 1.e-10;
              if (1. - std::fabs(ps)> tol)
                {
                  // facette et face_bord ne sont pas coplanaires, on conserve la facette!
                }
              else
                {
                  // facette et face_bord sont dans le meme plan, il faut supprimer cette facette
                  // car c'est bien une facette de bord.
                  res = 1;
#ifndef NDEBUG
                  if (impr_>9000)
                    {
                      Process::Journal()<<"  fa7_bord -> Supp ";
                      maillage.printFa7(fa7,0,Process::Journal());
                    }
#endif
                  for (isom=1 ; isom<nb_som_par_facette ; isom++)
                    {
                      facettes(fa7,isom) = facettes(fa7,0);
                    }
                  // Attention, il semble important de sortir sinon, on reinterroge une facette nulle...
                  // gain de temps et semble risque sinon..
                  break;
                }
            }
        }
#endif
    }

  maillage.check_mesh();

  if (mp_sum(res) > 0)
    res = 1;

  return res;
}

/*! @brief fonction outil utilisee par Remaillage_FT::supprimer_petites_aretes.
 *
 * (SPA = Supprimer Petites Aretes)
 *   Cette fonction parcourt le tableau "tab_aretesMarquees" qui est un resultat
 *   de "marquer_aretes". tab_aretesMarquees contient une liste d'aretes
 *   trop petites qu'on veut supprimer.
 *   On supprime ces aretes en remplacant l'un des sommets de l'arete par l'autre
 *   extremite. En general, deux triangles disparaissent (sauf arete
 *   de bord). Plusieurs aretes adjacentes au meme point peuvent etre supprimees.
 *   Pour chaque arete, on choisit le sommet a supprimer et le sommet a conserver.
 *   Puis on verifie qu'un sommet supprime n'est pas utilise par une autre arete
 *   pour etre conserve (non trivial en parallele).
 *
 * @param (sommets_remplacement) On resize le tableau a (maillage.nb_sommets(), 2), puis on y stocke pour chaque sommet: - sommets_remplacement(i,0) = -1, sommets_remplacement(i,1) = -1 (sommet non supprime) - sommets_remplacement(i,0) = PE, sommets_remplacement(i,1) = j (sommet a remplacer par le sommet reel d'indice local j sur PE) En sortie, sommets_remplacement a son espace_virtuel a jour.
 */
static void SPA_choisir_sommets_remplacement(const Maillage_FT_Disc& maillage,
                                             const IntTab& tab_aretesMarquees,
                                             IntTab& sommets_remplacement)
{
  const int       nb_sommets       = maillage.nb_sommets();
  const IntTab&    facettes         = maillage.facettes();
  const ArrOfInt& sommet_PE_owner  = maillage.sommet_PE_owner();
  const ArrOfInt& sommet_num_owner = maillage.sommet_num_owner();

  sommets_remplacement.resize(nb_sommets, 2);
  // On initialise le tableau a "no_PE": on mettra ensuite
  // les processeurs d'accord entre eux avec un
  //   "collecter(Descripteur_FT::MIN_COLONNE1)"
  // Si tout le monde contient "no_PE", alors le noeud ne devra pas etre remplace.
  const int no_PE = Process::nproc();
  sommets_remplacement = no_PE;

  // Pour chaque sommet, y a-t-il un processeur qui veut le conserver ?
  // (le sommet sert a remplacer d'autres sommets)
  // 0 => non  sinon => oui
  ArrOfIntFT sommets_conserves(nb_sommets);
  sommets_conserves = 0;

  // Parcours de la liste des aretes a supprimer. Pour chacune, on choisit
  // un sommet a conserver et un sommet a supprimer (le plus petit en numerotation
  // globale pour que le choix soit identique pour les deux triangles adjacents).
  const int n = tab_aretesMarquees.dimension(0);
  const int nb_som_par_facette = facettes.dimension(1);
  int i;
  for (i = 0; i < n; i++)
    {
      const int fa7    = tab_aretesMarquees(i,0);
      const int isom   = tab_aretesMarquees(i,1);
      assert(isom >= 0 && isom < nb_som_par_facette);
      const int isom_s = (isom + 1 < nb_som_par_facette) ? isom + 1 : 0;
      const int som    = facettes(fa7, isom);
      const int som_s  = facettes(fa7, isom_s);
      const int bord   = maillage.sommet_ligne_contact(som);
      const int bord_s = maillage.sommet_ligne_contact(som_s);

      int somRempl = -1; // Le sommet de remplacement (num_owner)
      int somSupp  = -1; // Le sommet supprime (indice local)
      if (bord == bord_s)
        {
          //les deux sommets sont sur le bord ou les 2 sont internes
          //on compare alors leur numero global
          const int pe         = sommet_PE_owner[som];
          const int pe_s       = sommet_PE_owner[som_s];
          const int numOwner   = sommet_num_owner[som];
          const int numOwner_s = sommet_num_owner[som_s];
          if (FTd_compare_sommets_global(pe, numOwner, pe_s, numOwner_s) > 0)
            {
              //som_s>som : on garde som, on supprime som_s
              somRempl = som;
              somSupp = som_s;
            }
          else
            {
              //som_s<=som : on garde som_s, on supprime som
              somRempl = som_s;
              somSupp = som;
            }
        }
      else if (bord==1)
        {
          //som est de bord : on garde som, on supprimer som_s
          somRempl = som;
          somSupp = som_s;
        }
      else
        {
          //som_s est de bord : on garde som_s, on supprime som
          somRempl = som_s;
          somSupp = som;
        }
      // Remplacement du sommet si
      // * le sommet a supprimer ne l'a pas encore ete
      // * le sommet de remplacement n'a pas encore ete supprime
      // * le sommet a supprimer n'est pas utilise pour en remplacer un autre
      if (   sommets_remplacement(somRempl, 0) == no_PE
             && sommets_remplacement(somSupp, 0) == no_PE
             && sommets_conserves[somSupp] == 0)
        {
          const int pe  = sommet_PE_owner[somRempl];
          const int le_som = sommet_num_owner[somRempl];
          sommets_remplacement(somSupp, 0) = pe;
          sommets_remplacement(somSupp, 1) = le_som;
          sommets_conserves[somRempl] = 1;
        }
    }

  // Si un sommet est conserve sur un processeur, il ne doit etre supprime
  // sur aucun processeur. On fait la somme pour tous les processeurs
  // de "sommets_conserves" : si un sommet est conserve sur l'un, il devra
  // etre conserve sur tous.
  const Desc_Structure_FT& desc_sommets = maillage.desc_sommets();
  desc_sommets.collecter_espace_virtuel(sommets_conserves,
                                        MD_Vector_tools::EV_SOMME);
  desc_sommets.echange_espace_virtuel(sommets_conserves);

  // On met aussi tous les processeurs d'accord sur le sommet de remplacement.
  // Choix arbitraire de l'un des sommets de remplacement parmi ceux
  // proposes par les differents processeurs.
  desc_sommets.collecter_espace_virtuel(sommets_remplacement,
                                        MD_Vector_tools::EV_MINCOL1);
  desc_sommets.echange_espace_virtuel(sommets_remplacement);

  // Mise a jour de sommets_remplacement : -1 pour les sommets a conserver
  for (i = 0; i < nb_sommets; i++)
    {
      const int a_conserver = sommets_conserves[i];
      const int non_remplace = (sommets_remplacement(i,0) == no_PE);
      if (a_conserver || non_remplace)
        {
          sommets_remplacement(i,0) = -1;
          sommets_remplacement(i,1) = -1;
        }
    }
}

/*! @brief A l'aide de "marquer_aretes", on determine les aretes trop petites du maillage.
 *
 * On les supprime en remplacant un sommet
 *   d'une extremite de l'arete par le sommet de l'autre extremite.
 *   Pour cela, on cree des sommets virtuels, des facettes nulles (deux ou
 *   trois sommets confondus) et on change le volume des phases.
 *   La variation de volume engendree est ajoutee a varVolume.
 *
 * @param (maillage) Le maillage a optimiser. Le maillage retourne a l'etat minimal, Le volume change, Certaines facettes sont "nulles" (plusieurs sommets confondus), On cree des sommets virtuels.
 * @param (varVolume) Un tableau de taille nb_sommets() contenant une valeur initiale de variation de volume. On augmente la taille du tableau (creation de sommets virtuels) et on ajoute la variation de volume du maillage due a la suppression des aretes. L'espace virtuel est a jour. Valeur de retour: nombre total (sur tous les procs) de sommets du maillage supprimes
 */
int Remaillage_FT::supprimer_petites_aretes(Maillage_FT_Disc& maillage,
                                            ArrOfDouble& varVolume) const
{
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "Supprimer_petites_aretes");
  statistiques().begin_count(stat_counter);

  int nb_sommets_supprimes_tot = 0;
  int nb_sommets_supprimes = 0;
  do
    {

      // Pour chaque sommet a remplacer:
      // colonne0: indice local du sommet a remplacer
      // colonne1: indice local du sommet de remplacement
      IntTabFT remplacement_ilocal(0,2);

      {
        // ******************************************************
        // Recherche des aretes trop petites
        IntTabFT tab_aretesMarquees;
        {
          ArrOfIntFT tab_somD;             // Inutilise
          DoubleTabFT tab_deplacement_somD; // Inutilise
          marquer_aretes(maillage,
                         tab_aretesMarquees,
                         tab_somD,
                         tab_deplacement_somD,
                         -1 /* marquage des aretes trop petites */);
        }

        // ******************************************************
        // On supprime les aretes en remplacant une extremite par l'autre
        // (genere des triangles dont deux sommets sont confondus)
        // Determination des sommets a remplacer et des sommets de remplacement
        // en numerotation globale.
        IntTabFT sommets_remplacement;
        SPA_choisir_sommets_remplacement(maillage,
                                         tab_aretesMarquees,
                                         sommets_remplacement);

        // ******************************************************
        // Certains sommets de remplacement n'existent pas encore sur le processeur local.
        // Creation de la liste des sommets de remplacement virtuels dont on aura besoin
        ArrOfIntFT request_sommets_pe;
        ArrOfIntFT request_sommets_num;
        int i;
        const int nb_sommets = maillage.nb_sommets();
        const int moi = Process::me();
        assert(nb_sommets == sommets_remplacement.dimension(0));
        for (i = 0; i < nb_sommets; i++)
          {
            const int pe = sommets_remplacement(i,0);
            const int som= sommets_remplacement(i,1);
            if (pe >= 0)
              {
                int j;
                if (pe != moi)   // Le sommet de remplacement est virtuel
                  {
                    request_sommets_pe.append_array(pe);
                    request_sommets_num.append_array(som);
                    j = -1;
                  }
                else
                  {
                    j = som;
                  }
                remplacement_ilocal.append_line(i, j);
              }
          }
        // ******************************************************
        // Creation des sommets virtuels qui remplaceront les sommets supprimes
        // (certains existent peut-etre deja, ils ne seront pas recrees)

        maillage.creer_sommets_virtuels_numowner(request_sommets_pe,
                                                 request_sommets_num);
        // On a cree de nouveaux sommets virtuels. Mise a jour de varVolume:
        // Inutile de mettre a jour l'espace virtuel, l'echange sera fait
        // apres le calcul de dvolume
        {
          const int old_size = varVolume.size_array();
          const int new_size = maillage.nb_sommets();
          varVolume.resize_array(new_size);
          int ii;
          for (ii = old_size; ii < new_size; ii++)
            varVolume[ii] = 0.;
        }
        // ******************************************************
        // Calcul de l'indice local des sommets de remplacement
        ArrOfIntFT request_sommets_ilocal;

        maillage.convertir_numero_distant_local(maillage.desc_sommets(),
                                                maillage.sommet_num_owner(),
                                                request_sommets_num,
                                                request_sommets_pe,
                                                request_sommets_ilocal);
        int j = 0;
        const int n_rempl = remplacement_ilocal.dimension(0);
        for (i = 0; i < n_rempl; i++)
          {
            int som_new = remplacement_ilocal(i,1);
            if (som_new < 0)   // Le sommet de remplacement est virtuel ?
              {
                som_new = request_sommets_ilocal[j];
                remplacement_ilocal(i,1) = som_new;
                j++;
              }
          }
      }
      // ******************************************************
      // Determination de la variation de volume liee a la suppression des aretes
      // On calcule la variation de volume correspondant au deplacement des
      // noeuds remplaces vers les noeuds de remplacement.
      {
        const DoubleTab& sommets = maillage.sommets();
        // Nombre de sommets a remplacer
        const int n_rempl = remplacement_ilocal.dimension(0);
        DoubleTabFT position_finale(sommets); // Copie du tableau.
        int i;
        const int dim = Objet_U::dimension;
        for (i = 0; i < n_rempl; i++)
          {
            const int som_old = remplacement_ilocal(i,0);
            const int som_new = remplacement_ilocal(i,1);
            int j;
            for (j = 0; j < dim; j++)
              position_finale(som_old, j) = position_finale(som_new, j);
          }
        {
          ArrOfDoubleFT dvolume;
          calculer_variation_volume(maillage,
                                    position_finale,
                                    dvolume);

          // Avant d'affecter la variation de volume des noeuds remplaces aux noeuds
          // de remplacement: certains sommets de remplacement sont virtuels,
          // il faut donc "collecter" les contributions a la fin. Donc il faut
          // annuler la valeur de varVolume pour les sommets virtuels avant
          // de commencer.
          const int n = maillage.nb_sommets();
          for (i = 0; i < n; i++)
            {
              if (maillage.sommet_virtuel(i))
                varVolume[i] = 0.;
              else
                varVolume[i] += dvolume[i];
            }
        }
        for (i = 0; i < n_rempl; i++)
          {
            const int som_old = remplacement_ilocal(i,0);
            const int som_new = remplacement_ilocal(i,1);
            const double dv = varVolume[som_old];
            varVolume[som_new] += dv;
            varVolume[som_old] = 0;
          }
        const Desc_Structure_FT& desc = maillage.desc_sommets();
        desc.collecter_espace_virtuel(varVolume, MD_Vector_tools::EV_SOMME);
        desc.echange_espace_virtuel(varVolume);
      }

      // ******************************************************
      // Remplacement des sommets dans le tableau des facettes,
      // On rend invalides les facettes qui disparaissent (deux sommets confondus)
      {
        const int nb_som = maillage.nb_sommets();
        // Creation d'une table de remplacement de sommets:
        ArrOfIntFT table_old_new(nb_som);
        table_old_new = -1;
        int i;
        const int nb_rempl = remplacement_ilocal.dimension(0);
        for (i = 0; i < nb_rempl; i++)
          {
            const int som_old = remplacement_ilocal(i,0);
            const int som_new = remplacement_ilocal(i,1);
            table_old_new[som_old] = som_new;
          }
        IntTab& facettes = maillage.facettes_;
        const int nb_facettes = facettes.dimension(0);
        const int nb_som_par_facette = facettes.dimension(1);
        for (i = 0; i < nb_facettes; i++)
          {
            int j;
            int flag = 0;
            for (j = 0; j < nb_som_par_facette; j++)
              {
                const int i_sommet = facettes(i,j);
                const int new_sommet = table_old_new[i_sommet];
                if (new_sommet >= 0)
                  {
                    facettes(i,j) = new_sommet;
                    flag = 1;
                  }
              }
            // En dimension2, une facette est invalide si les deux sommets
            // sont identiques. En dimension 3, si deux sommets sont confondus,
            // il faut invalider la facette (sommet0 == sommet1)
            if (flag && nb_som_par_facette == 3)
              {
                // Facette modifiee, on teste si elle disparait
                const int s0 = facettes(i,0);
                const int s1 = facettes(i,1);
                const int s2 = facettes(i,2);
                if (s0 == s2 || s1 == s2)
                  {
                    facettes(i,1) = s0;
                    facettes(i,2) = s0;
                  }
              }
          }

      }
      nb_sommets_supprimes = remplacement_ilocal.dimension(0);
      nb_sommets_supprimes = Process::mp_sum(nb_sommets_supprimes);
      nb_sommets_supprimes_tot += nb_sommets_supprimes;
    }
  while (nb_sommets_supprimes > 0);

  maillage.corriger_proprietaires_facettes();
  // corriger_proprietaires_facettes() cree eventuellement des sommets virtuels.
  // => Mise a jour de varVolume
  varVolume.resize_array(maillage.nb_sommets());
  maillage.desc_sommets().echange_espace_virtuel(varVolume);

  maillage.maillage_modifie(Maillage_FT_Disc::MINIMAL);
  statistiques().end_count(stat_counter);
  return nb_sommets_supprimes_tot;
}

/*! @brief Cette fonction marque a supprimer les facettes en double.
 *
 * Lorsque 2 facettes se trouvent avoir 3 de leur sommets en commun,
 *      il faut supprimer les 2 facettes?
 *     Marquer a supprimer = condenser les 3 sommets en un seul (le sommet 0)
 *
 * Precondition:
 *    Le maillage doit verifier la propriete "proprietaire facette = premier sommet"
 *
 * @param (maillage) maillage a remailler
 * @return (int) Le nombre de facettes supprimees
 */
int Remaillage_FT::supprimer_doublons_facettes(Maillage_FT_Disc& maillage) const
{
  if (Comm_Group::check_enabled()) maillage.check_mesh();

  const int nb_sommets = maillage.nb_sommets();
  const int nb_facettes = maillage.nb_facettes();
  IntTab& facettes = maillage.facettes_;
  const int nb_som_par_facette = facettes.dimension(1);

  // Marqueurs des facettes supprimees:
  ArrOfIntFT facettes_a_supprimer(nb_facettes);
  facettes_a_supprimer = 0;

  ArrOfIntFT fa7VoisinesSom_index;
  IntTabFT fa7VoisinesSom_data;

  fa7VoisinesSom_index.resize_array(nb_sommets);
  fa7VoisinesSom_data.resize(3*nb_facettes,2); //estimation nb connectivites = 3 x nb de facettes

  //on commence par recalculer les connectivites sommet -> facettes voisines
  calculer_connectivites_sommetFacettes(maillage, fa7VoisinesSom_index,fa7VoisinesSom_data);

  //puis on balaye les sommets
  int som, index, index2, fa7,fa72, compteur, isom,isom2, trouve;
  for (som=0 ; som<nb_sommets ; som++)
    {
      index = fa7VoisinesSom_index[som];
      while (index>=0)
        {
          //recupere la facette voisine
          fa7 = fa7VoisinesSom_data(index,1);
          if (facettes(fa7,0)!=facettes(fa7,1))
            {
              //facette valide
              index2 = fa7VoisinesSom_data(index,0);
              while (index2>=0)
                {
                  //recupere la facette voisine
                  fa72 = fa7VoisinesSom_data(index2,1);
                  if (facettes(fa72,0)!=facettes(fa72,1))
                    {
                      //facette valide
                      //compte le nombre de sommet communs entre les facettes
                      compteur = 0;
                      trouve = 0;
                      for (isom=0 ; isom<nb_som_par_facette && trouve==0 ; isom++)
                        {
                          for (isom2=0 ; isom2<nb_som_par_facette && trouve==0 ; isom2++)
                            {
                              if (facettes(fa7,isom) == facettes(fa72,isom2))
                                {
                                  trouve = 1;
                                }
                            }
                          if (trouve==1)
                            {
                              //on a bien trouve le sommet facettes(fa7,isom) dans fa72
                              compteur++;
                              trouve = 0; //pour pouvoir continuer le test pour le sommet suivant de fa7
                            }
                          else
                            {
                              //on n'a pas trouve le sommet facettes(fa7,isom) dans fa72
                              //on peut donc arreter le test ici
                              trouve = -1;
                            }
                        }
                      if (compteur==nb_som_par_facette)
                        {
#ifndef NDEBUG
                          if (impr_>5000)
                            {
                              Process::Journal()<<"Facettes superposees : doublon a supprimer"<<finl;
                              Process::Journal()<<"  ";
                              maillage.printFa7(fa7,0,Process::Journal());
                              Process::Journal()<<"  ";
                              maillage.printFa7(fa72,0,Process::Journal());
                            }
#endif
                          index2 = -1;
                          //marque les facettes a supprimer, ie ses 3 sommets sont identiques
                          facettes_a_supprimer[fa7] = 1;
                          facettes_a_supprimer[fa72] = 1;
                        }
                      else
                        {
                          //index suivant
                          index2 = fa7VoisinesSom_data(index2,0);
                        }
                    }
                  else
                    {
                      //index suivant
                      index2 = fa7VoisinesSom_data(index2,0);
                    }
                }
            }
          //index suivant
          index = fa7VoisinesSom_data(index,0);
        }
    }
  // Mise a jour des espaces virtuels pour supprimer les facettes sur tous les procs
  {
    const Desc_Structure_FT& desc_facettes = maillage.desc_facettes();
    desc_facettes.echange_espace_virtuel(facettes_a_supprimer);
  }
  // On rend les facettes a supprimer invalides (trois sommets identiques)
  {
    int i, j;
    for (i = 0; i < nb_facettes; i++)
      {
        if (facettes_a_supprimer[i])
          {
            const int sommet0 = facettes(i, 0);
            for (j = 1; j < nb_som_par_facette; j++)
              facettes(i, j) = sommet0;
          }
      }
  }
  maillage.maillage_modifie(Maillage_FT_Disc::MINIMAL);

  return 1;
}


/*! @brief Cette fonction calcule le carre de la longueur ideale d'une arete Dans un premier, cette longueur ideale est la racine cubique (en 3D) de la moyenne
 *
 *      des volumes des elements euleriens contenant les sommets
 *
 * @param (som0) indice du 1er sommet de l'arete
 * @param (som1) indice du 2e sommet de l'arete
 * @return (double) carre de la longueur idealede l'arete
 */
double Remaillage_FT::calculer_longueurIdeale2_arete(const Maillage_FT_Disc& maillage,
                                                     int som0,
                                                     double x, double y, double z) const
{
  double lgrId2 = 0.;
  if (facteur_longueur_ideale_ > 0.)
    {
      // La longueur ideale est un multiple de la taille locale des
      // elements euleriens.
      const int elem0 = maillage.sommet_elem_[som0];

      const Elem_geom_base& elem_geom = refdomaine_VF_->domaine().type_elem().valeur();

      if (sub_type(Triangle, elem_geom))
        {
          double vol = refdomaine_VF_->volumes(elem0);
          lgrId2 = 2. * vol;

        }
      else if (sub_type(Tetraedre, elem_geom))
        {

          double vol = refdomaine_VF_->volumes(elem0);
          lgrId2 = 2. * pow(vol,2./dimension);

        }
      else if (sub_type(Rectangle, elem_geom)
               || sub_type(Hexaedre, elem_geom))
        {

          const DoubleTab& sommets = maillage.sommets();
          const DoubleTab& xv = refdomaine_VF_->xv();
          const IntTab& elem_faces = refdomaine_VF_->elem_faces();
          int k;
          FTd_vecteur3 v = {0., 0., 0.};
          FTd_vecteur3 xyz = {x, y, z};
          FTd_vecteur3 delta_xv = {0., 0., 0.};
          const int dim = Objet_U::dimension;
          double norme2 = 0.;
          for (k = 0; k < dim; k++)
            {
              v[k] = xyz[k] - sommets(som0, k);
              norme2 += v[k] * v[k];
              const int face0 = elem_faces(elem0,k);
              const int face1 = elem_faces(elem0,k+dim);
              delta_xv[k] = xv(face1,k) - xv(face0,k);
            }
          if (equilateral_)
            {
              // On calcul la diagonale de l'element eulerien. Puis longueur au carre.
              norme2 = 0.;
              for (k = 0; k < dim; k++)
                norme2 += delta_xv[k] * delta_xv[k];
            }
          else
            {
              // On calcul la projection de la diagonale de l'element eulerien
              // sur l'arete de la facette.  Puis longueur au carre.
              if (norme2 == 0)
                {
                  v[0] = 1.;
                  v[1] = 1.;
                  v[2] = dim==3 ? 1. : 0.;
                  norme2 = dim;
                }
              double f = 1. / sqrt(norme2);
              norme2 = 0.;
              for (k = 0; k < dim; k++)
                {
                  v[k] *= f * delta_xv[k];
                  norme2 += v[k] * v[k];
                }
            }
          lgrId2 = norme2;

        }
      else
        {
          Cerr << "Remaillage_FT::calculer_longueurIdeale2_arete non implemente pour les elements "
               << elem_geom.que_suis_je() << finl;
          assert(0==1);
          exit();
        }

      lgrId2 *= facteur_longueur_ideale_ * facteur_longueur_ideale_;
    }
  else if (valeur_longueur_fixe_ > 0.)
    {
      // La longueur ideale est imposee en valeur absolue dans le jeu de donnees
      // (valeur en metres)
      lgrId2 = valeur_longueur_fixe_ * valeur_longueur_fixe_;
    }
  else
    {
      Cerr << "Erreur Remaillage_FT::calculer_longueurIdeale2_arete" << finl;
      exit();
    }
  return lgrId2;
}

/*! @brief Cette fonction calcule la variation de volume liee a la suppression de sommets
 *
 * @param (maillage) maillage a barycentrer
 * @param (tab_somSupp) tableau des sommets a supprimer tab[som] < 0  : sommet a conserver tab[som] >=0  : sommet a remplacer par tab[som]
 * @param (varVolume) la variation de volume
 * @return (double) volume total supprime
 */
double Remaillage_FT::calculer_volume_sommets_supprimes(const Maillage_FT_Disc& maillage, const ArrOfInt& tab_somSupp,
                                                        ArrOfDouble& varVolume) const
{
  double res = 0.;

  const int nb_facettes = maillage.nb_facettes();
  const DoubleTab& sommets = maillage.sommets();
  const IntTab& facettes = maillage.facettes();
  const int nb_som_par_facette = facettes.dimension(1);

  int isom,som, som_rempl,fa7, k;
  double volume;
  FTd_vecteur3 som0,som1,som2,som3;

  //const int nb_sommets = maillage.nb_sommets();
  //for (som=0 ; som<nb_sommets ; som++) {
  //  varVolume[som] = 0.;
  //}

  //on balaye les facettes
  for (fa7=0 ; fa7<nb_facettes ; fa7++)
    {
      //on balaye ses sommets
      for (isom=0 ; isom<nb_som_par_facette ; isom++)
        {
          som = facettes(fa7,isom);
          som_rempl = tab_somSupp[som];
          if (som_rempl>=0)
            {
              //sommet som a remplacer par som_rempl
              //la variation de volume est alors egale au volume genere par les facettes voisines
              //de som (et le 4e sommet des tetraedres = som_rempl)
              //on recupere les coordonnees des sommets de la facette et du sommet de remplacement
              for (k=0 ; k<dimension ; k++)
                {
                  som0[k] = sommets(facettes(fa7,0),k);
                  som1[k] = sommets(facettes(fa7,1),k);
                  som2[k] = sommets(facettes(fa7,2),k);
                  som3[k] = sommets(som_rempl,k);
                }
              //on calcule le volume du tetraedre som0/som1/som2/som_rempl
              volume = FTd_calculer_volume_tetraedre(som0,som1,som2,som3);
              res += volume;
              varVolume[som] = volume;
            }
        }
    }

#ifndef NDEBUG
  if (impr_>9000)
    {
      const int nb_sommets = maillage.nb_sommets();
      for (som=0 ; som<nb_sommets ; som++)
        {
          if (tab_somSupp[som]!=-1)
            Process::Journal()<<"Som_apresSupp "<<tab_somSupp[som]<<" "<<som<<"  varVolume= "<<varVolume[som]<<finl;
        }
    }
#endif
#ifndef NDEBUG
  if (impr_>5000)
    {
      Process::Journal()<<"calculer_volume_sommets_supprimes : varVolume_global= "<<res<<finl;
    }
#endif
  return res;
}

/*! @brief Cette fonction effectue des permutations d'aretes afin d'ameliorer la qualite du maillage
 *
 * @param (maillage) maillage a remailler
 * @return (int) 1 si le remaillage s'est deroule correctement, 0 sinon
 */
int Remaillage_FT::permuter_aretes(Maillage_FT_Disc& maillage) const
{
  int res = 1;
  if (dimension!=3)
    {
      //ceci n'est valable qu'en 3D
      return res;
    }


  maillage.maillage_modifie(Maillage_FT_Disc::MINIMAL);
  return res;
}

//=============================================================================================
// FONCTIONS DE DIVISION DES GRANDES ARETES

static const IntTabFT* static_tab_sort;
// Fonction de comparaison lexicographique de deux aretes a scinder:
// Note : la fonction utilise non pas les donnees du tableau passe en parametre,
//   mais les donnees du tableau static_tab_sort pour la comparaison.
// Tri par ordre croissant de 2e colonne, puis 3e, puis 4e colonne.
True_int fct_compare_tab_aretes(const void *pt1, const void *pt2)
{
  const int index1 = *(const int *) pt1;
  const int index2 = *(const int *) pt2;

  True_int x = FTd_compare_sommets_global((*static_tab_sort)(index1,0), (*static_tab_sort)(index1,1), (*static_tab_sort)(index2,0), (*static_tab_sort)(index2,1));
  if (x==0)
    {
      x = FTd_compare_sommets_global((*static_tab_sort)(index1,2), (*static_tab_sort)(index1,3), (*static_tab_sort)(index2,2), (*static_tab_sort)(index2,3));
    }

  return x;
}

/*! @brief Cette fonction divise les grandes aretes en 2
 *
 * @param (maillage) maillage a remailler
 * @return (int) 1 si le remaillage s'est deroule correctement, 0 sinon
 */
int Remaillage_FT::diviser_grandes_aretes(Maillage_FT_Disc& maillage) const
{
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "Diviser_grandes_aretes");
  statistiques().begin_count(stat_counter);

  static int compteur = 0;
  static int test_val = -1;

  Process::Journal()<<"Remaillage_FT::diviser_grandes_aretes "<<temps_<<"  nb_som="<<maillage.nb_sommets()<<finl;
  Process::Journal()<<" Compteur = " << compteur << finl;
  compteur++;
  if (compteur == test_val)
    {
      Process::Journal() << " STOP." << finl;
    }
  //  int res = 1;

  maillage.nettoyer_elements_virtuels();

  //tableaux de stockage
  IntTabFT tab_aretesMarquees;
  ArrOfIntFT tab_somD;
  DoubleTabFT tab_deplacement_somD;

  //on commence par marquer les grandes aretes
  marquer_aretes(maillage,
                 tab_aretesMarquees,
                 tab_somD,
                 tab_deplacement_somD,
                 1 /* marquage des aretes trop grandes */);
  // resultat =  tab_aretesMarquees : [ fa7 iarete pe som ]

  const int nb_aretes_divis = tab_aretesMarquees.dimension(0);


  int nb_facettes = maillage.nb_facettes();
  const int nb_facettes0 = nb_facettes;
  IntTab& facettes = maillage.facettes_;
  const int nb_som_par_facette = facettes.dimension(1);

  const ArrOfInt& sommet_num_owner = maillage.sommet_num_owner_;
  int fa7,iarete, isom,som,isom_s,som_s,isom_ss,som_ss, pe_somD,numOwner_somD,somD,somD_s,somD_ss;

  const int dimension3 = (dimension==3);
  const int nb_aretes_par_facette = (dimension3)?3:1;
  //tableau stockant le sommet servant a decouper l'arete (ou -1 si arete non divisee)
  IntTabFT tab_fa7Divis(nb_facettes,nb_aretes_par_facette);
  for (fa7=0 ; fa7<nb_facettes ; fa7++)
    {
      for (isom=0 ; isom<nb_aretes_par_facette ; isom++)
        {
          tab_fa7Divis(fa7,isom) = -1;
        }
    }

  //on va balayer les aretes a memoriser l'ensemble des aretes a scinder
  for (iarete=0 ; iarete<nb_aretes_divis ; iarete++)
    {
      fa7 = tab_aretesMarquees(iarete,0);
      isom = tab_aretesMarquees(iarete,1);
      isom_s = (isom+1)%nb_som_par_facette;
      //sommet qui va rester dans fa7
      som = facettes(fa7,isom);
      //sommet qui va aller dans une nouvelle facette
      som_s = facettes(fa7,isom_s);
      //sommet a inserer dans l'arete
      pe_somD = tab_aretesMarquees(iarete,2);
      numOwner_somD = tab_aretesMarquees(iarete,3);

      if (pe_somD!=me())
        {
          //je ne connais pas l'indice du sommet a inserer dans me()
          maillage.convertir_numero_distant_local(maillage.desc_sommets(),sommet_num_owner,numOwner_somD,pe_somD,somD);
          assert(somD >= 0);
        }
      else
        {
          //je suis le proprietaire de somD : je connais donc le bon indice du sommet a inserer
          somD = numOwner_somD;
        }

      tab_fa7Divis(fa7,isom) = somD;
    }

  //on va ensuite balayer les facettes et les scinder
  //la configuration depend du nb d'aretes a scinder par facette
  int nb_areteScinder, isom0=-1,isom1=-1;
  for (fa7=0 ; fa7<nb_facettes0 ; fa7++)
    {
      //on compte le nombre d'aretes a scinder
      nb_areteScinder = 0;
      for (isom=0 ; isom<nb_aretes_par_facette ; isom++)
        {
          if (tab_fa7Divis(fa7,isom)>=0)
            {
              if (nb_areteScinder==0)
                {
                  isom0 = isom;
                }
              else
                {
                  isom1 = isom;
                }
              nb_areteScinder++;
            }
        }
      if (nb_areteScinder==1)
        {
          //s'il n'y a qu'une arete a scinder
          somD = tab_fa7Divis(fa7,isom0);
          //on modifie la facettes d'origine
          isom_s = (isom0+1)%nb_som_par_facette;
          som_s = facettes(fa7,isom_s);
          facettes(fa7,isom_s) = somD;
          //on cree une nouvelle facette
          if (nb_facettes>=facettes.dimension(0))
            {
              //Redimensionnement
              facettes.resize(nb_facettes+10,facettes.dimension(1));
            }
          isom_ss = (isom_s+1)%nb_som_par_facette;
          som_ss = facettes(fa7,isom_ss);
          if (dimension==2)
            {
              facettes(nb_facettes,isom0) = somD;
              facettes(nb_facettes,isom_s) = som_s;
            }
          else
            {
              facettes(nb_facettes,0) = som_ss;
              facettes(nb_facettes,1) = somD; //sommet insere
              facettes(nb_facettes,2) = som_s;
            }
          nb_facettes++;
        }
      else if (nb_areteScinder==2)
        {
          //si 2 aretes sont a scinder
          //on cree deux nouvelles facettes
          if (nb_facettes+2>=facettes.dimension(0))
            {
              //Redimensionnement
              facettes.resize(nb_facettes+12,facettes.dimension(1));
            }
          //on positionne isom et isom_s tq elles soient les aretes a scinder
          if (isom1==(isom0+1)%nb_som_par_facette)
            {
              isom = isom0;
              isom_s = isom1;
            }
          else
            {
              isom = isom1;
              isom_s = isom0;
            }
          isom_ss = ((isom_s+1)%nb_som_par_facette);
          som = facettes(fa7,isom);
          som_s = facettes(fa7,isom_s);
          som_ss = facettes(fa7,isom_ss);
          somD = tab_fa7Divis(fa7,isom);
          somD_s = tab_fa7Divis(fa7,isom_s);
          somD_ss = tab_fa7Divis(fa7,isom_ss);
          //on modifie la facette existante
          facettes(fa7,0) = som;
          facettes(fa7,1) = somD;
          facettes(fa7,2) = som_ss;
          //premiere nouvelle facette
          facettes(nb_facettes,0) = somD;
          facettes(nb_facettes,1) = som_s;
          facettes(nb_facettes,2) = somD_s;
          nb_facettes++;
          //seconde nouvelle facette
          facettes(nb_facettes,0) = somD;
          facettes(nb_facettes,1) = somD_s;
          facettes(nb_facettes,2) = som_ss;
          nb_facettes++;
        }
      else if (nb_areteScinder==3)
        {
          //si toutes les aretes sont a scinder
          //on cree trois nouvelles facettes
          if (nb_facettes+3>=facettes.dimension(0))
            {
              //Redimensionnement
              facettes.resize(nb_facettes+13,facettes.dimension(1));
            }
          isom = 0;
          isom_s = 1;
          isom_ss = 2;
          som = facettes(fa7,isom);
          som_s = facettes(fa7,isom_s);
          som_ss = facettes(fa7,isom_ss);
          somD = tab_fa7Divis(fa7,isom);
          somD_s = tab_fa7Divis(fa7,isom_s);
          somD_ss = tab_fa7Divis(fa7,isom_ss);
          //on modifie la facette existante
          facettes(fa7,0) = som;
          facettes(fa7,1) = somD;
          facettes(fa7,2) = somD_ss;
          //premiere nouvelle facette
          facettes(nb_facettes,0) = somD;
          facettes(nb_facettes,1) = som_s;
          facettes(nb_facettes,2) = somD_s;
          nb_facettes++;
          //seconde nouvelle facette
          facettes(nb_facettes,0) = somD_s;
          facettes(nb_facettes,1) = som_ss;
          facettes(nb_facettes,2) = somD_ss;
          nb_facettes++;
          //troisieme nouvelle facette
          facettes(nb_facettes,0) = somD;
          facettes(nb_facettes,1) = somD_s;
          facettes(nb_facettes,2) = somD_ss;
          nb_facettes++;
        }
    }
  //Redimensionnement
  facettes.resize(nb_facettes,facettes.dimension(1));
  maillage.desc_facettes_.calcul_schema_comm(nb_facettes);
  maillage.corriger_proprietaires_facettes();

  ArrOfIntFT liste_sommets_sortis;
  ArrOfIntFT numero_face_sortie;
  maillage.deplacer_sommets(tab_somD,tab_deplacement_somD,liste_sommets_sortis,numero_face_sortie);
  maillage.corriger_proprietaires_facettes();

  Process::Journal()<<"FIN Remaillage_FT::diviser_grandes_aretes "<<temps_<<"  nb_som="<<maillage.nb_sommets()
                    <<"  nb_aretes_divisees="<< nb_aretes_divis<<finl;
  int nb_aretes_divis_tot = Process::mp_sum(nb_aretes_divis);

  Process::Journal()<<"FIN Remaillage_FT::diviser_grandes_aretes " <<"  nb_aretes_divisees_tot="<< nb_aretes_divis_tot<<finl;

  maillage.maillage_modifie(Maillage_FT_Disc::MINIMAL);
  statistiques().end_count(stat_counter);
  //  return res;
  return nb_aretes_divis_tot;
}

//cette fonction insere une ligne de "requete arete" dans le tableau
//renvoie l'indice ou cela a ete mis dans le tableau
int Remaillage_FT::inserer_tab_aretes(int& nb_tab_aretes,IntTab& tab_aretes, DoubleTab& tab_criteres,
                                      int pe0, int numOwner0, int pe1, int numOwner1, int sommet_face_bord1,
                                      int peRequete, int fa7_peR, int iarete_fa7_peR) const
{
  int res = nb_tab_aretes;
  //Redimensionnement eventuel
  if (nb_tab_aretes>=tab_aretes.dimension(0))
    {
      tab_aretes.resize(nb_tab_aretes+10,tab_aretes.dimension(1));
      tab_criteres.resize(nb_tab_aretes+10,tab_criteres.dimension(1));
    }
  //insertion des donnees a la fin
  tab_aretes(nb_tab_aretes,0) = pe0;
  tab_aretes(nb_tab_aretes,1) = numOwner0;
  tab_aretes(nb_tab_aretes,2) = pe1;
  tab_aretes(nb_tab_aretes,3) = numOwner1;
  tab_aretes(nb_tab_aretes,4) = sommet_face_bord1;
  tab_aretes(nb_tab_aretes,5) = peRequete;
  tab_aretes(nb_tab_aretes,6) = fa7_peR;
  tab_aretes(nb_tab_aretes,7) = iarete_fa7_peR;

  nb_tab_aretes++;

  return res;
}

//Partant de l'indice index (dans tab_index), cette fonction renvoie l'index suivant correspondant a l'arete (pe0/som0,pe1/som1) passee
//Premiere utilisation : passer index = -1
//Renvoie tmp=-1 s'il n'y a pas d'arete correspondante
int Remaillage_FT::chercher_arete_tab(int tmp, const ArrOfInt& tab_index, const IntTab& tab_aretes,
                                      int pe0, int numOwner0, int pe1, int numOwner1) const
{
  int pe0_=-1,numOwner0_=-1, pe1_=-1,numOwner1_=-1;
  int iarete;
  const int nb_tab_aretes = tab_index.size_array();
  tmp++;
  if (tmp>=nb_tab_aretes)
    {
      return -1;
    }
  int trouve = -1;
  while (tmp<nb_tab_aretes && trouve==-1)
    {
      iarete = tab_index[tmp];
      pe0_ = tab_aretes(iarete,0);
      numOwner0_ = tab_aretes(iarete,1);
      pe1_ = tab_aretes(iarete,2);
      numOwner1_ = tab_aretes(iarete,3);
      if ( pe0_<pe0 || (pe0_==pe0 && numOwner0_<numOwner0)
           || (pe0_==pe0 && numOwner0_==numOwner0 && pe1_<pe1)
           || (pe0_==pe0 && numOwner0_==numOwner0 && pe1_==pe1 && numOwner1_<numOwner1) )
        {
          tmp++;
        }
      else
        {
          trouve = 0;
        }
    }
  if (! (pe0_==pe0 && numOwner0_==numOwner0 && pe1_==pe1 && numOwner1_==numOwner1) )
    {
      tmp = -1;
    }

  return tmp;
}

//Methode qui permet de collecter la liste des aretes marquees.
// ATTENTION : creation de sommets supplementaires au milieu des aretes trop longues
//
//Si drap > 0 : les facettes marquees sont celles qui sont trop grandes,
//   ie elles ne verifient pas le critere suivant:
//
//             (L_arete)^2
//     -------------------------  < 1.
//      (L_ideale^2 x (1 + C)^2
//
//Si drap < 0 : les facettes marquees sont celles qui sont trop petites,
//   ie elles ne verifient pas le critere suivant:
//
//             (L_arete)^2
//     -------------------------  > 1.
//      (L_ideale^2 x (1 - C)^2
//
//La methode remplit le tableau tab_aretesMarquees, dont les colonnes sont
//        [ fa7 iarete pe som ]
// ou
//  -fa7 : indice de la facette contenant l'arete marquee
//  -iarete : indice localc de l'arete marquee
//  -peD/somD : numerotation globale du sommet a utiliser lors d'une division d'arete
//
//Rq : les deplacement des sommets peD/somD est a realiser en-dehors de cette fonction (liste des sommets dans tab_somD)
//en utilisant le tableau tab_deplacement_somD, et apres avoir realise la subdivision des aretes
int Remaillage_FT::marquer_aretes(Maillage_FT_Disc& maillage, IntTab& tab_aretesMarquees,
                                  ArrOfInt& tab_somD, DoubleTab& tab_deplacement_somD, int drap) const
{
  int res = 1;
  int nb_aretesMarquees = 0;
  tab_aretesMarquees.resize(10,4); //tableau des aretes marquees [fa7 iarete pe som]
  //ou pe som est le sommet a utiliser lors d'une division d'arete

  int nb_somD = 0;
  tab_somD.resize_array(10);
  tab_deplacement_somD.resize(10,dimension);
  static IntTabFT tab_arete_somD(10,4);
  static IntTabFT tab_aretes(10,10); //tableau [0=pe0 1=ind0 2=pe1 3=ind1 4=face_bord1 5=peReq 6=fa7 7=iarete 8=critere 9=som]
  //ou (pe0) som est le sommet a utiliser lors d'une division d'arete
  int nb_tab_aretes = 0;
  static DoubleTabFT tab_criteres(10,6); //tableau [pox posy posz L^2 LI0^2 LI1^2]

  const int dimension3 = (dimension==3);

  const int nb_sommets = maillage.nb_sommets();
  const int nb_facettes = maillage.nb_facettes();
  const IntTab& facettes = maillage.facettes();
  DoubleTab& sommets = maillage.sommets_;
  ArrOfInt& sommet_elem = maillage.sommet_elem_;
  ArrOfInt& sommet_face_bord = maillage.sommet_face_bord_;
  ArrOfInt& sommet_PE_owner = maillage.sommet_PE_owner_;
  ArrOfInt& sommet_num_owner = maillage.sommet_num_owner_;
  ArrOfInt& drapeaux_sommets = maillage.drapeaux_sommets_;
  const int nb_som_par_facette = facettes.dimension(1);
  const int nb_aretes_par_facette = (dimension3)?3:1;
  const int _MOI_ = me();

  //on commence par recuperer le schema de communication
  //on prend le schema inverse du descripteur de sommets car
  const Schema_Comm_FT& comm = maillage.desc_sommets_.schema_comm_inverse();
  //initialise la communication
  comm.begin_comm();

  int fa7, isom,som,pe,numOwner, isom_s,som_s,pe_s,numOwner_s, indice, tmp;
  double x,y,z=0.;
  //on balaye les facettes reelles
  for (fa7=0 ; fa7<nb_facettes ; fa7++)
    {
      //on balaye les aretes de la facette
      if (!maillage.facette_virtuelle(fa7) && facettes(fa7,0)!=facettes(fa7,1))
        {
          //facette non virtuelle et non supprimee
          //on balaye les aretes de la facette
          for (isom=0 ; isom<nb_aretes_par_facette ; isom++)
            {
              //on recupere les sommets de l'arete
              isom_s = (isom+1)%nb_som_par_facette;
              som = facettes(fa7,isom);
              som_s = facettes(fa7,isom_s);
              pe = sommet_PE_owner[som];
              numOwner = sommet_num_owner[som];
              pe_s = sommet_PE_owner[som_s];
              numOwner_s = sommet_num_owner[som_s];
              tmp = FTd_compare_sommets_global(pe,numOwner,pe_s,numOwner_s);
              assert(tmp!=0);
              if (tmp>0)
                {
                  //on a en fait som>som_s : on permute les sommets
                  tmp = som;
                  som = som_s;
                  som_s = tmp;
                  pe = sommet_PE_owner[som];
                  numOwner = sommet_num_owner[som];
                  pe_s = sommet_PE_owner[som_s];
                  numOwner_s = sommet_num_owner[som_s];
                }
              if (pe!=_MOI_)
                {
                  //on envoie les informations au processeur possedant som
                  //pour qu'il puisse calculer sa partie du critere
                  Sortie& buffer = comm.send_buffer(pe);
                  //-les numerotations globales des sommets
                  buffer << pe << numOwner << pe_s << numOwner_s << sommet_face_bord[som_s];
                  //-les indices de la  fa7 et de 'larete
                  buffer << fa7 << isom;
                  //-la position du sommet som_s
                  buffer << sommets(som_s,0) << sommets(som_s,1);
                  if (dimension3)
                    buffer << sommets(som_s,2);
                }
              else
                {
                  //je suis pe le proprietaire de som
                  //on ajoute la ligne de donnees aux tableaux
                  //insere le tableau en tant que requete, avec moi comme requerant
                  tmp = inserer_tab_aretes(nb_tab_aretes,tab_aretes,tab_criteres,pe,numOwner,pe_s,numOwner_s, sommet_face_bord[som_s], _MOI_,fa7,isom);
                  tab_criteres(tmp,0) = sommets(som_s,0);
                  tab_criteres(tmp,1) = sommets(som_s,1);
                  tab_criteres(tmp,2) = (dimension3) ? sommets(som_s,2) : 0;
                  tab_criteres(tmp,3) = tab_criteres(tmp,4) = tab_criteres(tmp,5) = -1.;
                }
              if (pe_s!=pe)
                {
                  //test pour eviter de refaire le bloc ci-dessus
                  if (pe_s!=_MOI_)
                    {
                      //on envoie aussi les informations au processeur possedant som_s
                      //pour qu'il puisse calculer sa partie du critere
                      Sortie& buffer = comm.send_buffer(pe_s);
                      //-les numerotations globales des sommets
                      buffer << pe << numOwner << pe_s << numOwner_s << sommet_face_bord[som_s];
                      //-les indices de la  fa7 et de l'arete
                      buffer << fa7 << isom;
                      //-la position du sommet som
                      buffer << sommets(som,0) << sommets(som,1);
                      if (dimension3)
                        buffer << sommets(som,2);
                    }
                  else
                    {
                      //je suis pe_s le proprietaire de som_s (et la ligne n'a pas deja ete ajoutee ci-dessu)
                      //on ajoute la ligne de donnees aux tableaux
                      //insere le tableau en tant que requete, avec moi comme requerant
                      tmp =inserer_tab_aretes(nb_tab_aretes,tab_aretes,tab_criteres,pe,numOwner,pe_s,numOwner_s, sommet_face_bord[som_s], -1,-1,-1);
                      tab_criteres(tmp,0) = sommets(som,0);
                      tab_criteres(tmp,1) = sommets(som,1);
                      tab_criteres(tmp,2) = (dimension3) ? sommets(som,2) : 0;
                      tab_criteres(tmp,3) = tab_criteres(tmp,4) = tab_criteres(tmp,5) = -1.;
                    }
                }
            } //fin du balayage des aretes de la facette
        } //fin facette non virtuelle
    } //fin du balayage des facettes

  int peV, pe0,numOwner0, pe1,numOwner1, iarete, face_bord1;
  //on lance l'echange des messages
  comm.echange_taille_et_messages();

  //on recupere les messages
  //il s'agit de message que l'on m'a envoye si je suis
  //  -le proprietaire du sommet 0 de l'arete ou
  //  -le proprietaire du sommet 1 de l'arete
  //le message est constitue de pe0 / numOwner0 / pe1 / numOwner1 / fa7 / isom / posx / posy / posz
  //ou
  // -fa7 est l'indice dans le pe envoyant de la facette contenant l'arete
  // -isom est l'indeice de l'arete dans cette facette
  // -(posX, posy, posz) est la position
  //   -du sommet 1 si je suis pe0
  //   -du sommet 0 si je suis pe1
  {
    const ArrOfInt& pe_voisins = comm.get_recv_pe_list();
    const int nb_pe_voisins = pe_voisins.size_array();
    z = 0;
    for (int indice_pe = 0; indice_pe < nb_pe_voisins; indice_pe++)
      {
        peV = pe_voisins[indice_pe];
        Entree& buffer = comm.recv_buffer(peV);
        buffer >> pe0;
        while (!buffer.eof())
          {
            buffer >> numOwner0 >> pe1 >> numOwner1 >> face_bord1;
            buffer >> fa7 >> isom;
            buffer >> x >> y;
            if (dimension3)
              buffer >> z;

            assert(pe0==_MOI_ || pe1==_MOI_);
            //on insere alors la ligne de requete
            tmp = inserer_tab_aretes(nb_tab_aretes,tab_aretes,tab_criteres,pe0,numOwner0,pe1,numOwner1,face_bord1, peV,fa7,isom);
            tab_criteres(tmp,0) = x;
            tab_criteres(tmp,1) = y;
            tab_criteres(tmp,2) = z;
            tab_criteres(tmp,3) = tab_criteres(tmp,4) = tab_criteres(tmp,5) = -1.;
            buffer >> pe0;
          }
      }
  }
  //finit cette communication
  comm.end_comm();

  //Redimensionnements
  tab_aretes.resize(nb_tab_aretes,tab_aretes.dimension(1));
  tab_criteres.resize(nb_tab_aretes,tab_criteres.dimension(1));

  //a la fin de cela, tout processeur proprietaire d'un sommet connait
  // -toutes les aretes partant de ce sommet
  // -la position de tous les seconds sommets de l'arete

  //on va maintenant trier les tableaux tab_aretes et tab_criteres en fonction des numeros globaux des sommets de l'arete
  //pour cela on va utiliser un tableau d'indices annexe (pour ne pas modifier les tableaux eux-memes)
  ArrOfIntFT tab_index(nb_tab_aretes);
  for (indice=0 ; indice<nb_tab_aretes ; indice++)
    {
      tab_index[indice] = indice;
    }
  //on initialise le tableau servant au tri de la liste d'index
  static_tab_sort = &tab_aretes;
  //on recupere un tableau d'indices triant tab_fa7AScinder selon ses colonnes [1 2 3]
  qsort(tab_index.addr(),nb_tab_aretes,sizeof(int), fct_compare_tab_aretes);

  double lgr2=-1.,lgr_ideale02=-1.,lgr_ideale12,d;
  int pe0p,numOwner0p,pe1p,numOwner1p;
  pe0p = numOwner0p = pe1p = numOwner1p = -1;

  //on balaye alors les lignes des tableaux tab_*, pour calculer
  // -la longueur de l'arete si je suis le proprietaire du sommet 0
  // -la partie 0 de longueur ideale de l'arete si je suis le proprietaire du sommet 0
  // -la partie 1 de longueur ideale de l'arete si je suis le proprietaire du sommet 1
  //Une partie de longueur ideale est la longueur ideale calculee via l'element eulerien contenant le sommet
  for (indice=0 ; indice<nb_tab_aretes ; indice++)
    {
      iarete = tab_index[indice];
      pe0 = tab_aretes(iarete,0);
      numOwner0 = tab_aretes(iarete,1);
      pe1 = tab_aretes(iarete,2);
      numOwner1 = tab_aretes(iarete,3);
      if (FTd_compare_sommets_global(pe1,numOwner1,pe1p,numOwner1p)!=0 || FTd_compare_sommets_global(pe0,numOwner0,pe0p,numOwner0p)!=0)
        {
          //on est sur une arete differente de la precedente
          //on va donc recalculer les longueurs
          lgr2 = lgr_ideale02 = lgr_ideale12 = -1.;
          if (pe0==_MOI_)
            {
              //je suis le proprietaire du sommet 0
              //je dois donc calculer la longueur reelle de l'arete
              d = tab_criteres(iarete,0) - sommets(numOwner0,0);
              lgr2 = d*d;
              d = tab_criteres(iarete,1) - sommets(numOwner0,1);
              lgr2 += d*d;
              if (dimension3)
                {
                  d = tab_criteres(iarete,2) - sommets(numOwner0,2);
                  lgr2 += d*d;
                }
              //je dois aussi calculer ma longueur ideale
              lgr_ideale02 = calculer_longueurIdeale2_arete(maillage,numOwner0,tab_criteres(iarete,0),tab_criteres(iarete,1),tab_criteres(iarete,2));
            }
          if (pe1==_MOI_)
            {
              //je suis le proprietaire du sommet 1
              //je dois donc calculer ma longueur ideale
              lgr_ideale12 = calculer_longueurIdeale2_arete(maillage,numOwner1,tab_criteres(iarete,0),tab_criteres(iarete,1),tab_criteres(iarete,2));
            }
        }
      //et on stocke
      tab_criteres(iarete,3) = lgr2;
      tab_criteres(iarete,4) = lgr_ideale02;
      tab_criteres(iarete,5) = lgr_ideale12;
      pe0p = pe0;
      numOwner0p = numOwner0;
      pe1p = pe1;
      numOwner1p = numOwner1;
    }

  //on construit un nouveau schema de communication :
  // -un processeur proprietaire du sommet 0 saura qu'il devra recevoir des infos venant du proprietaire du sommet 1
  // -un processeur proprietaire du sommet 1 saura qu'il devra envoyer des infos au proprietaire du sommet 0
  iarete = 0;
  static Schema_Comm_FT comm2;
  {
    static ArrOfIntFT recv_pe_list;    // liste des pe qui recoivent des infos
    static ArrOfIntFT send_pe_list;    // liste des pe qui envoient des infos
    ArrOfIntFT drap_per(Process::nproc());
    ArrOfIntFT drap_pes(Process::nproc());
    drap_per = 0;
    drap_pes = 0;
    send_pe_list.resize_array(0);
    recv_pe_list.resize_array(0);
    //pour cela, on va alors balayer les aretes
    for (iarete=0 ; iarete<nb_tab_aretes ; iarete++)
      {
        pe0 = tab_aretes(iarete,0);
        pe1 = tab_aretes(iarete,2);
        if (pe0!=pe1)
          {
            if (pe0==_MOI_)
              {
                //je suis pe0 : je dois donc m'attendre a recevoir des infos de pe1
                if (drap_per[pe1]==0)
                  recv_pe_list.append_array(pe1);
                else
                  drap_per[pe1] = 1;
              }
            else
              {
                //je suis pe1 : je dois donc m'attendre a devoir envoyer des infos a pe0
                if (drap_pes[pe0]==0)
                  send_pe_list.append_array(pe0);
                else
                  drap_pes[pe0] = 1;
              }
          }
      }
    array_trier_retirer_doublons(send_pe_list);
    array_trier_retirer_doublons(recv_pe_list);
    comm2.set_send_recv_pe_list(send_pe_list, recv_pe_list);
  }
  comm2.begin_comm();

  //Les proprietaires de sommet 1 doivent alors envoyer leur partie de longueur ideale
  // aux proprietaires des somemts 0 (pour que celui-ci puisse calculer la longueur ideale finale)
  pe0p = numOwner0p = pe1p = numOwner1p = -1;
  for (indice=0 ; indice<nb_tab_aretes ; indice++)
    {
      iarete = tab_index[indice];
      pe0 = tab_aretes(iarete,0);
      pe1 = tab_aretes(iarete,2);
      if (pe1==_MOI_ && pe0!=_MOI_)
        {
          //je suis le proprietaire de numOwner1, mais pas de numOwner0
          //je dois donc envoyer a pe0 mon calcul de la longueur ideale
          numOwner0 = tab_aretes(iarete,1);
          numOwner1 = tab_aretes(iarete,3);
          if (FTd_compare_sommets_global(pe1,numOwner1,pe1p,numOwner1p)!=0 || FTd_compare_sommets_global(pe0,numOwner0,pe0p,numOwner0p)!=0)
            {
              //on est sur une arete differente de la precedente : on doit donc envoyer le message
              lgr_ideale12 = tab_criteres(iarete,5);
              Sortie& buffer = comm2.send_buffer(pe0);
              buffer << pe0 << numOwner0 << pe1 << numOwner1 << lgr_ideale12;
              pe0p = pe0;
              numOwner0p = numOwner0;
              pe1p = pe1;
              numOwner1p = numOwner1;
            }
        }
    }
  //on lance l'echange des messages
  comm2.echange_taille_et_messages();

  //on recupere les messages :
  //ce sont des messages venant des proprietaires de sommets 1
  //et m'envoyant leur partie de calcul de la longueur ideale
  {
    const ArrOfInt& pe_voisins = comm2.get_recv_pe_list();
    const int nb_pe_voisins = pe_voisins.size_array();
    int peVoi, pe0bis,numOwner0bus, pe1bis,numOwner1bis, trouve;
    for (int indice_pe = 0; indice_pe < nb_pe_voisins; indice_pe++)
      {
        peVoi = pe_voisins[indice_pe];
        Entree& buffer = comm2.recv_buffer(peVoi);
        buffer >> pe0bis;
        while (!buffer.eof())
          {
            buffer >> numOwner0bus >> pe1bis >> numOwner1bis >> lgr_ideale12;
            assert(pe1bis==peVoi && pe0bis!=peVoi);
            tmp = -1;
            trouve = -1;
            do
              {
                tmp = chercher_arete_tab(tmp,tab_index,tab_aretes, pe0bis,numOwner0bus,pe1bis,numOwner1bis);
                if (tmp>=0)
                  {
                    trouve = 1;
                    iarete = tab_index[tmp];
                    tab_criteres(iarete,5) = lgr_ideale12;
                  }
              }
            while (tmp>=0);
            if (trouve==-1)
              {
                Process::Journal()<<"PB : n'a pu insere l'arete ";
                Process::Journal()<<"  arete=("<<tab_aretes(iarete,0)<<"-"<<tab_aretes(iarete,1)<<" / "<<tab_aretes(iarete,2)<<"-"<<tab_aretes(iarete,3);
                Process::Journal()<<"  crit=("<<tab_criteres(iarete,3)<<" / "<<tab_criteres(iarete,4)<<" / "<<tab_criteres(iarete,5)<<finl;
                tmp = chercher_arete_tab(tmp,tab_index,tab_aretes, pe0bis,numOwner0bus,pe1bis,numOwner1bis);
                assert(0==1);
              }
            buffer >> pe0bis;
          }
      }
  }
  //finit la communication
  comm2.end_comm();

  //a la fin de cela, je connais toutes les infos necessaires pour calculer le critere sur les aretes
  //dont je suis proprietaire du sommet 0
  //on va donc calculer le critere pour chacune des aretes dont je suis le proprietaire du sommet 0
  int NV_nb_sommets = maillage.nb_sommets();
  iarete = 0;
  {
    ArrOfIntFT liste_pe, liste_sommets;
    // Rapport entre longueur minimale ou maximale de l'arete et la longueur ideale:
    const double tolerance_longueur_arete = 1 + critere_arete_ * ((drap>0) ? 1. : -1);
    // Carre de ce rapport (on va comparer le rapport des carres des longueurs):
    const double coeff2 = tolerance_longueur_arete * tolerance_longueur_arete;
    double critere, lgrIdeale2;
    int test=-1;
    //pour cela, on va alors balayer les aretes, dans l'ordre trie
    pe0p = numOwner0p = pe1p = numOwner1p = -1;
    for (indice=0 ; indice<nb_tab_aretes ; indice++)
      {
        iarete = tab_index[indice];
        pe0 = tab_aretes(iarete,0);
        if (pe0==_MOI_)
          {
            //je suis le proprietaire de som0
            numOwner0 = tab_aretes(iarete,1);
            pe1 = tab_aretes(iarete,2);
            numOwner1 = tab_aretes(iarete,3);
            if (FTd_compare_sommets_global(pe1,numOwner1,pe1p,numOwner1p) != 0
                || FTd_compare_sommets_global(pe0,numOwner0,pe0p,numOwner0p) != 0)
              {
                //on est sur une arete differente de la precedente
                //on va donc recalculer le critere
                //on calcule la moyenne des carres des longueurs ideales calculees en chacun des sommets
                lgrIdeale2 = 0.5 * (tab_criteres(iarete,4) + tab_criteres(iarete,5));
                assert(lgrIdeale2>0);
                //calcul du critere  : lrg2/ (lgrIdeale2 * coeff2)  ou
                // -coeff2 = (1+C)^2 si on veut tester les aretes a diviser (drap>0)
                // -coeff2 = (1-C)^2 si on veut tester les aretes a supprimer (drap<0)
                lgr2 = tab_criteres(iarete,3);
                assert(lgr2>=0.);
                critere = lgr2 / (lgrIdeale2 * coeff2);
                //le test a verifier est
                // (critere<1) si on veut tester les aretes a diviser (drap>0)
                // (critere>1) si on veut tester les aretes a supprimer (drap<0)
                test = (drap>0)? (critere<1) : (critere>1);
              }
            tab_aretes(iarete,8) = test;
            if (drap>0 && !test)
              {
                // le test n'est pas verifie pour une arete trop grande :
                // on doit creer un sommet
                // verifie d'abord si le sommet n'a pas deja ete cree
                tmp = -1;
                for (isom=0 ; isom<nb_somD && tmp==-1 ; isom++)
                  {
                    if (FTd_compare_sommets_global(pe0,numOwner0,tab_arete_somD(isom,0),tab_arete_somD(isom,1))==0
                        && FTd_compare_sommets_global(pe1,numOwner1,tab_arete_somD(isom,2),tab_arete_somD(isom,3))==0)
                      {
                        tmp = isom;
                      }
                  }
                if (tmp==-1)
                  {
                    //sommet non trouve : le cree
                    if (NV_nb_sommets>=sommets.dimension(0))
                      {
                        tmp = NV_nb_sommets+10;
                        sommets.resize(tmp, sommets.dimension(1));
                        sommet_elem.resize_array(tmp);
                        sommet_face_bord.resize_array(tmp);
                        sommet_PE_owner.resize_array(tmp);
                        sommet_num_owner.resize_array(tmp);
                        drapeaux_sommets.resize_array(tmp);
                      }
                    sommets(NV_nb_sommets,0) = sommets(numOwner0,0);
                    sommets(NV_nb_sommets,1) = sommets(numOwner0,1);
                    if (dimension3)
                      sommets(NV_nb_sommets,2) = sommets(numOwner0,2);
                    sommet_elem[NV_nb_sommets] = sommet_elem[numOwner0];
                    face_bord1 = tab_aretes(iarete,4);
                    maillage.convertir_numero_distant_local(maillage.desc_sommets_,sommet_num_owner,numOwner1,pe1,som_s);
                    if (face_bord1!=-1)
                      {
                        sommet_face_bord[NV_nb_sommets] = sommet_face_bord[numOwner0];
                      }
                    else
                      {
                        sommet_face_bord[NV_nb_sommets] = -1;
                      }
                    sommet_PE_owner[NV_nb_sommets]  = sommet_PE_owner[numOwner0];
                    sommet_num_owner[NV_nb_sommets] = NV_nb_sommets;
                    drapeaux_sommets[NV_nb_sommets] = drapeaux_sommets[numOwner0];
                    //on stocke les infos du nouveau sommet
                    if (nb_somD>=tab_somD.size_array())
                      {
                        tmp = nb_somD+10;
                        tab_somD.resize_array(tmp);
                        tab_deplacement_somD.resize(tmp,dimension);
                        tab_arete_somD.resize(tmp,tab_arete_somD.dimension(1));
                      }
                    //indice
                    tab_somD[nb_somD] = NV_nb_sommets;
                    //deplacement
                    tab_deplacement_somD(nb_somD,0) = 0.5 * ( tab_criteres(iarete,0) - sommets(NV_nb_sommets,0) );
                    tab_deplacement_somD(nb_somD,1) = 0.5 * ( tab_criteres(iarete,1) - sommets(NV_nb_sommets,1) );
                    if (dimension3)
                      tab_deplacement_somD(nb_somD,2) = 0.5 * ( tab_criteres(iarete,2) - sommets(NV_nb_sommets,2) );
                    //arete
                    tab_arete_somD(nb_somD,0) = tab_aretes(iarete,0);
                    tab_arete_somD(nb_somD,1) = tab_aretes(iarete,1);
                    tab_arete_somD(nb_somD,2) = tab_aretes(iarete,2);
                    tab_arete_somD(nb_somD,3) = tab_aretes(iarete,3);
                    nb_somD++;
                    tmp = NV_nb_sommets;
                    NV_nb_sommets++;
                  }
                else
                  {
                    tmp = tab_somD[tmp];
                  }
                //on stocke aussi l'indice du sommet dans les infos des aretes
                tab_aretes(iarete,9) = tmp;

                int peReq = tab_aretes(iarete,5);
                if (peReq!=_MOI_)
                  {
                    liste_pe.append_array(peReq);
                    liste_sommets.append_array(tmp);
                  }
              }
          }
      }
    test = (NV_nb_sommets!=nb_sommets);
    test = Process::mp_sum(test);
    if (test>0)
      {
        //le nb de sommets a ete modifie :
        //redimensionnement
        tmp = NV_nb_sommets;
        sommets.resize(tmp, sommets.dimension(1));
        sommet_elem.resize_array(tmp);
        sommet_face_bord.resize_array(tmp);
        sommet_PE_owner.resize_array(tmp);
        sommet_num_owner.resize_array(tmp);
        drapeaux_sommets.resize_array(tmp);
        maillage.desc_sommets_.calcul_schema_comm(tmp);
        //creation des sommets virtuels dans les processeurs le necessitant
        const Schema_Comm_FT& comm_ = maillage.desc_sommets().schema_comm();
        maillage.creer_sommets_virtuels(liste_sommets,liste_pe,comm_);
      }
  }

  // Les proprietaires des sommets envoient les resultats du test aux
  // PEs qui ont demande le test (ce sont les proprietaires des aretes,
  // chez qui les sommets sont virtuels). C'est donc le schema de comm direct:
  //  proprietaire du sommet envoie des donnees aux PEs qui ont des sommets virtuels.
  const Schema_Comm_FT& comm_real_to_virt = maillage.desc_sommets_.schema_comm();
  comm_real_to_virt.begin_comm();
  //on va maintenant balayer les aretes pour stocker celles ne verifiant pas le test
  int peReq, test;
  pe0p = numOwner0p = pe1p = numOwner1p = -1;
  for (indice=0 ; indice<nb_tab_aretes ; indice++)
    {
      iarete = tab_index[indice];
      pe0 = tab_aretes(iarete,0);
      if (pe0==_MOI_)
        {
          //je suis le proprietaire de som0
          numOwner0 = tab_aretes(iarete,1);
          pe1 = tab_aretes(iarete,2);
          numOwner1 = tab_aretes(iarete,3);
          if (FTd_compare_sommets_global(pe1,numOwner1,pe1p,numOwner1p)!=0 || FTd_compare_sommets_global(pe0,numOwner0,pe0p,numOwner0p)!=0)
            {
              //on est sur une arete differente de la precedente
              test = tab_aretes(iarete,8);
              if (!test)
                {
                  //le critere n'est pas verifie
                  peReq = tab_aretes(iarete,5);
                  fa7 = tab_aretes(iarete,6);
                  isom = tab_aretes(iarete,7);
                  if (drap>0)
                    {
                      //on divise des aretes : un sommet a ajouter
                      som = tab_aretes(iarete,9);
                    }
                  else
                    {
                      //on supprime des aretes : pas de sommet a ajouter
                      som = -1;
                    }
                  if( peReq==_MOI_)
                    {
                      //je suis aussi celui qui avait besoin de l'info
                      //je stocke les infos dans mon tableau
                      if (nb_aretesMarquees>=tab_aretesMarquees.dimension(0))
                        {
                          tab_aretesMarquees.resize(nb_aretesMarquees+10,tab_aretesMarquees.dimension(1));
                        }
                      tab_aretesMarquees(nb_aretesMarquees,0) = fa7;
                      tab_aretesMarquees(nb_aretesMarquees,1) = isom;
                      tab_aretesMarquees(nb_aretesMarquees,2) = _MOI_;
                      tab_aretesMarquees(nb_aretesMarquees,3) = som;
                      nb_aretesMarquees++;
                    }
                  else
                    {
                      //je dois envoyer l'info a celui qui me l'a demandee
                      Sortie& buffer = comm_real_to_virt.send_buffer(peReq);
                      buffer << fa7 << isom << som;
                    }
                }
            }
        }
    }

  //on lance l'echange des messages
  comm_real_to_virt.echange_taille_et_messages();

  //on recupere les messages venant des proprietaires des sommets 0
  //lorsque ceux-ci ont trouve une arete ne verifiant pas le critere
  {
    const ArrOfInt& pe_voisins = comm_real_to_virt.get_recv_pe_list();
    const int nb_pe_voisins = pe_voisins.size_array();
    int peVbis;
    for (int indice_pe = 0; indice_pe < nb_pe_voisins; indice_pe++)
      {
        peVbis = pe_voisins[indice_pe];
        Entree& buffer = comm_real_to_virt.recv_buffer(peVbis);
        buffer >> fa7;
        while (!buffer.eof())
          {
            buffer >> isom >> som;
            assert(_MOI_!=peVbis);
            //je stocke les infos dans mon tableau
            if (nb_aretesMarquees>=tab_aretesMarquees.dimension(0))
              {
                tab_aretesMarquees.resize(nb_aretesMarquees+10,tab_aretesMarquees.dimension(1));
              }
            tab_aretesMarquees(nb_aretesMarquees,0) = fa7;
            tab_aretesMarquees(nb_aretesMarquees,1) = isom;
            tab_aretesMarquees(nb_aretesMarquees,2) = peVbis;
            tab_aretesMarquees(nb_aretesMarquees,3) = som;
            nb_aretesMarquees++;
            buffer >> fa7;
          }
      }
  }
  //termine cette communication
  comm_real_to_virt.end_comm();

  //Redimensionnement
  tab_aretesMarquees.resize(nb_aretesMarquees,tab_aretesMarquees.dimension(1));
  tab_somD.resize_array(nb_somD);
  tab_deplacement_somD.resize(nb_somD,tab_deplacement_somD.dimension(1));

  maillage.maillage_modifie(Maillage_FT_Disc::MINIMAL);

  return res;
}


/*! @brief Cette methode calcule, pour un triangle donne, sa qualite : celle-ci est comprise dans ]0,1], et vaut 1 pour un triangle equilateral.
 *
 * @param (som0) id du premier sommet
 * @param (som1) id du second sommet
 * @param (som2) id du troisieme sommet
 * @param (CoordSom) coordonnees des sommets
 * @return (double) le coefficient de qualite du triangle.
 */
double Remaillage_FT::qualiteTriangle(const FTd_vecteur3& som0, const FTd_vecteur3& som1, const FTd_vecteur3& som2, double& aire) const
{
  int k;

  //calcul de perimetre2 = somme des carres de aretes
  double perimetre2 = 0.,d;
  for (k=0 ; k<dimension ; k++)
    {
      d = som1[k] - som0[k];
      perimetre2 += d*d;
    }
  for (k=0 ; k<dimension ; k++)
    {
      d = som2[k] - som1[k];
      perimetre2 += d*d;
    }
  for (k=0 ; k<dimension ; k++)
    {
      d = som0[k] - som2[k];
      perimetre2 += d*d;
    }

  assert(perimetre2>0.);

  //calcul de l'aire du triangle
  aire = 0.;
  FTd_vecteur3 vect;
  vect[0] = (som1[1]-som0[1]) * (som2[2]-som0[2]) - (som1[2]-som0[2]) * (som2[1]-som0[1]);
  vect[1] = (som1[2]-som0[2]) * (som2[0]-som0[0]) - (som1[0]-som0[0]) * (som2[2]-som0[2]);
  vect[2] = (som1[0]-som0[0]) * (som2[1]-som0[1]) - (som1[1]-som0[1]) * (som2[0]-som0[0]);
  for (k=0 ; k<dimension ; k++)
    {
      aire += vect[k]*vect[k];
    }
  aire = sqrt(aire)/2.;

#ifndef NDEBUG
  if (impr_>9000 && aire<=0.)
    {
      Cerr<<"PB qualiteTriangle :"<<finl;
      Cerr<<"  som0 "<<" coord= "<<som0[0]<<" "<<som0[1]<<" "<<som0[2]<<finl;
      Cerr<<"  som1 "<<" coord= "<<som1[0]<<" "<<som1[1]<<" "<<som1[2]<<finl;
      Cerr<<"  som2 "<<" coord= "<<som2[0]<<" "<<som2[1]<<" "<<som2[2]<<finl;
      Cerr<<" perimetre= "<<perimetre2<<finl;
      Cerr<<"  vect= "<<vect[0]<<" "<<vect[1]<<" "<<vect[2]<<"  aire= "<<aire<<finl;
    }
#endif
  // sqrt(3) : erreur sur xlC
  const double sq3x4 = sqrt(3.) * 4.;

  double qual = sq3x4 * aire / perimetre2;

  return qual;
}

/*! @brief Cette fonction "supprime" les facettes de surface nulle : Elle les reduit a 1 sommet (le sommet 0, pour ne pas changer de processeur)
 *
 * @param (maillage) maillage a remailler
 * @return (int) 1 si le remaillage s'est deroule correctement, 0 sinon
 */
int Remaillage_FT::supprimer_facettes_nulles(Maillage_FT_Disc& maillage) const
{
  int res = 1;

  Process::Journal()<<"Remaillage_FT::supprimer_facettes_nulles nb_fa7="<<maillage.nb_facettes()<<finl;
  maillage.check_mesh();

  IntTab& facettes = maillage.facettes_;
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const int nb_facettes = maillage.nb_facettes();
  const int nb_som_par_facette = facettes.dimension(1);

  int fa7,isom, nb_fa70 = 0;
  //on balaye les facettes pour marquer celles a supprimer
  for (fa7=0 ; fa7<nb_facettes ; fa7++)
    {
      if (surface_facettes[fa7]<=0.)
        {
          for (isom=1 ; isom<nb_som_par_facette ; isom++)
            {
              facettes(fa7,isom) = facettes(fa7,0);
            }
          nb_fa70++;
        }
    }

  // On a modifie le maillage
  maillage.maillage_modifie(Maillage_FT_Disc::MINIMAL);

  Process::Journal()<<"Remaillage_FT::supprimer_facettes_nulles nb_fa70="<<nb_fa70<<finl;
  maillage.check_mesh();

  //lance le nettoyage (statut minimal remis a cet endroit)
  nettoyer_maillage(maillage);

  Process::Journal()<<"FIN Remaillage_FT::supprimer_facettes_nulles nb_fa7="<<maillage.nb_facettes()<<finl;
  return res;
}
/*! @brief Cette fonction nettoie le maillage : Elle supprime les facettes reduites a 1 sommet
 *
 *     Elle supprime les sommets non utilises
 *
 * @param (maillage) maillage a remailler
 * @return (int) 1 si le remaillage s'est deroule correctement, 0 sinon
 */
int Remaillage_FT::nettoyer_maillage(Maillage_FT_Disc& maillage) const
{
  int res = 1;

  Process::Journal()<<"Remaillage_FT::nettoyer_maillage nb_fa7="<<maillage.nb_facettes()<<finl;

  //lance le nettoyage (statut minimal remis a cet endroit)
  maillage.nettoyer_maillage();

  Process::Journal()<<"FIN Remaillage_FT::nettoyer_maillage nb_fa7="<<maillage.nb_facettes()<<finl;
  return res;
}

/*! @brief Regularise le maillage en deplacant les sommets pour reduire les gradients de courbure.
 *
 *   L'equation d'evolution du maillage est une diffusion de masse le long
 *   de l'interface proportionnelle au gradient de courbure.
 *
 */
void Remaillage_FT::regulariser_courbure(Maillage_FT_Disc& maillage,
                                         const double coeff,
                                         ArrOfDouble& dvolume) const
{
  const IntTab& facettes = maillage.facettes();
  const DoubleTab& sommets = maillage.sommets();
  const int nb_facettes = facettes.dimension(0);
  const int dim = facettes.dimension(1);
  const ArrOfDouble& surface_facettes =
    maillage.get_update_surface_facettes();

  const double angle_bidim_axi =  bidim_axi ? Maillage_FT_Disc::angle_bidim_axi() : 0.;
  // Calcul de la courbure de l'interface
  const ArrOfDouble& courbure = maillage.get_update_courbure_sommets();

  // Calcul de la longueur de la plus petite arete
  double l_min = 1e10;
  int count = 0;
  int total_count = 0;
  const double one_third = 1. / 3.;

  dvolume = 0.;
  int i_facette;
  for (i_facette = 0; i_facette < nb_facettes; i_facette++)
    {
      // Ne pas calculer de flux pour les facettes virtuelles:
      if (maillage.facette_virtuelle(i_facette))
        continue;
      if (dim == 2)
        {
          const int s1 = facettes(i_facette, 0);
          const int s2 = facettes(i_facette, 1);
          const double x1 = sommets(s1, 0);
          const double x2 = sommets(s2, 0);
          const double dx = x2 - x1;
          const double dy = sommets(s2, 1) - sommets(s1, 1);
          const double l2 = dx * dx + dy * dy;
          if (l2 > 0.)
            {
              const double inv_l = 1. / sqrt(l2);
              const double inv_l2 = 1. / (l2);
              if (l2 < l_min)
                l_min = l2;

              const double c1 = courbure[s1];
              const double c2 = courbure[s2];
              const double gradient_c = (c2 - c1) / sqrt(l2);
              // We assume that the highest curvature achievable is 1/l (corresponds to an hexagon)
              // and that the maximal gradient is a 50%
              const double criterion1 = lissage_critere_*0.5*inv_l2;
              // Or we assume the variation around mean curvature is at maximum 50% of c_moy
              const double criterion2 = lissage_critere_*0.25*std::fabs(c1+c2)*inv_l;
              if ((std::fabs(gradient_c)>= criterion1)||(std::fabs(gradient_c)>= criterion2))
                {
                  double h = 1.;
                  if (bidim_axi)
                    h = (x1+x2) * 0.5 * angle_bidim_axi;
                  const double flux = gradient_c * h;
                  dvolume[s1] += flux;
                  dvolume[s2] -= flux;
                  count++;
                }
              total_count++;
            }
        }
      else
        {
          int som[3]; // Indices des trois sommets
          double c[3];   // Courbure des trois sommets
          double coord[3][3]; // Coordonnees
          int i, j;
          for (i = 0; i < 3; i++)
            {
              const int s = facettes(i_facette, i);
              som[i] = s;
              c[i] = courbure[s];
              for (j = 0; j < 3; j++)
                coord[i][j] = sommets(s, j);
            }
          const double surface = surface_facettes[i_facette];
          int i_suiv = 2;
          for (i = 0; i < 3; i++)
            {
              const double dx = coord[i_suiv][0] - coord[i][0];
              const double dy = coord[i_suiv][1] - coord[i][1];
              const double dz = coord[i_suiv][2] - coord[i][2];
              // Longueur de l'arete au carre
              const double l2 = dx * dx + dy * dy + dz * dz;
              if (l2 < l_min)
                l_min = l2;
              const double inv_l = (l2 == 0.) ? 1. : 1. / sqrt(l2);
              const double inv_l2 = (l2 == 0.) ? 1. : 1. / (l2);
              const double c1 = c[i];
              const double c2 = c[i_suiv];
              const double gradient_c = (c2 - c1) * inv_l;
              // We assume that the highest curvature achievable is 1/l (corresponds to an hexagon)
              // and that the maximal gradient is a 50%
              const double criterion1 = lissage_critere_*0.5*inv_l2;
              // Or we assume the variation around mean curvature is at maximum 50% of c_moy
              const double criterion2 = lissage_critere_*0.25*std::fabs(c1+c2)*inv_l;
              if ((std::fabs(gradient_c)>= criterion1)||(std::fabs(gradient_c)>= criterion2))
                {
                  // Hauteur d'un tiers du triangle (facette, la base du triangle etant [i,i_suiv])
                  const double h = (surface * inv_l) * one_third;
                  // Integrale du flux de masse sur le morceau de volume de controle autour des
                  // sommets (au coefficient de diffusion pres) :
                  const double flux = gradient_c * h;

                  const int s1 = som[i]; // Numero du premier sommet
                  const int s2 = som[i_suiv]; // Numero du deuxieme sommet
                  dvolume[s1] += flux;
                  dvolume[s2] -= flux;
                  count++;
                }
              total_count++;
              i_suiv = i;
            }
        }
    }
  l_min = Process::mp_min(l_min);
  if ((total_count)&&(lissage_critere_>DMINFLOAT))
    Journal() << "Proportion of smoothed aretes (similar to sommets) " << 100.*count/total_count << " %" << finl;
  // Calcul du coefficient de diffusion * pas de temps pour stabilite:
  // proportionnel a longueur^4 :
  const double coeff_dt = l_min * l_min * coeff;
  dvolume *= coeff_dt;

  maillage.desc_sommets().collecter_espace_virtuel(dvolume, MD_Vector_tools::EV_SOMME);
  maillage.desc_sommets().echange_espace_virtuel(dvolume);
}
