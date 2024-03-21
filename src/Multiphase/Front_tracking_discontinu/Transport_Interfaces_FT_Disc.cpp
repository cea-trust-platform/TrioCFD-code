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
// File:        Transport_Interfaces_FT_Disc.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/patch_168/1
//
//////////////////////////////////////////////////////////////////////////////

#include <Transport_Interfaces_FT_Disc.h>
#include <Probleme_FT_Disc_gen.h>
#include <Operateur.h>
#include <Motcle.h>
#include <EcritureLectureSpecial.h>
#include <Discret_Thyd.h>
#include <Schema_Temps_base.h>
#include <time.h>
#include <Parser.h>
#include <Domaine_VF.h>
#include <Domaine_VDF.h>
#include <Sous_Domaine.h>
#include <Champ_Face_VDF.h>
#include <Champ_P1NC.h>
#include <Champ_Fonc_P1NC.h>
#include <Fluide_Diphasique.h>
#include <Fluide_Incompressible.h>
#include <Statistiques.h>
#include <Paroi_FT_disc.h>
#include <Champ_Fonc_Face_VDF.h>
#include <Scatter.h>
#include <LecFicDiffuse.h>
#include <Sauvegarde_Reprise_Maillage_FT.h>
#include <SFichier.h>
#include <Debog.h>
#include <Format_Post_Lata.h>
#include <Connex_components.h>
#include <Connex_components_FT.h>
#include <Loi_horaire.h>
#include <Transport_Marqueur_FT.h>
#include <communications.h>
#include <Convection_Diffusion_Temperature_FT_Disc.h>
#include <Champ.h>
#include <Champ_Fonc_P0_base.h>
#include <EFichierBin.h>
#include <Param.h>
#include <sys/stat.h>
#include <TRUSTList.h>
#include <Domaine.h>
#include <TRUSTTrav.h>
#include <stat_counters.h>
#include <Dirichlet_paroi_fixe.h>
#include <Dirichlet_paroi_defilante.h>
#include <Echange_contact_VDF_FT_Disc.h>

Implemente_instanciable_sans_constructeur_ni_destructeur(Transport_Interfaces_FT_Disc,"Transport_Interfaces_FT_Disc",Transport_Interfaces_base);

/*! @brief Classe outil ou on stocke tout le bazar qui sert au fonctionnement de l'equation de transport.
 *
 * Permet de sauvegarder/reprendre les donnees d'interface
 *   (probleme des reprises avec pas de temps multiples),
 *   et de reduire les dependances en n'incluant presque rien dans Transport_Interfaces_FT_Disc.h
 *
 */
static void eval_vitesse(double x, double y, double z, double t,
                         Parser& px, Parser& py, Parser& pz,
                         double& vx, double& vy, double& vz)
{
  int i0=0, i1=1, i2=2, i3=3;
  px.setVar(i0,x);
  px.setVar(i1,y);
  px.setVar(i2,z);
  px.setVar(i3,t);

  vx = px.eval();
  py.setVar(i0,x);
  py.setVar(i1,y);
  py.setVar(i2,z);
  py.setVar(i3,t);

  vy = py.eval();
  pz.setVar(i0,x);
  pz.setVar(i1,y);
  pz.setVar(i2,z);
  pz.setVar(i3,t);
  vz = pz.eval();
}
static void integrer_vitesse_imposee(
  Parser& parser_vx, Parser& parser_vy, Parser& parser_vz,
  double temps, double dt, double& x, double& y, double& z)
{
  // Runge Kutta ordre 3:
  double vx, vy, vz;
  eval_vitesse(x,
               y,
               z,
               temps,
               parser_vx,parser_vy,parser_vz,vx,vy,vz);
  double vx1,vy1,vz1;
  eval_vitesse(x + vx * dt * 0.5,
               y + vy * dt * 0.5,
               z + vz * dt * 0.5,
               temps + dt * 0.5,
               parser_vx,parser_vy,parser_vz,vx1,vy1,vz1);
  double vx2,vy2,vz2;
  eval_vitesse(x + vx1 * dt * 0.5,
               y + vy1 * dt * 0.5,
               z + vz1 * dt * 0.5,
               temps + dt * 0.5,
               parser_vx,parser_vy,parser_vz,vx2,vy2,vz2);
  double vx3,vy3,vz3;
  eval_vitesse(x + vx2 * dt,
               y + vy2 * dt,
               z + vz2 * dt,
               temps + dt,
               parser_vx,parser_vy,parser_vz,vx3,vy3,vz3);

  x += (vx + 2.*vx1 + 2.*vx2 + vx3) / 6. * dt;
  y += (vy + 2.*vy1 + 2.*vy2 + vy3) / 6. * dt;
  z += (vz + 2.*vz1 + 2.*vz2 + vz3) / 6. * dt;
}

#if 0
static void normer_vecteurs(DoubleTab& tab)
{
  const int n = tab.dimension(0);
  int i;
  switch (tab.dimension(1))
    {
    case 2:
      for (i = 0; i < n; i++)
        {
          const double nx = tab(i,0);
          const double ny = tab(i,1);
          double n = nx*nx+ny*ny;
          if (n > 0.)
            {
              n = 1. / sqrt(n);
              tab(i,0) = nx * n;
              tab(i,1) = ny * n;
            }
        }
      break;
    case 3:
      for (i = 0; i < n; i++)
        {
          const double nx = tab(i,0);
          const double ny = tab(i,1);
          const double nz = tab(i,2);
          double n = nx*nx+ny*ny+nz*nz;
          if (n > 0.)
            {
              n = 1. / sqrt(n);
              tab(i,0) = nx * n;
              tab(i,1) = ny * n;
              tab(i,2) = nz * n;
            }
        }
      break;
    default:
      assert(0);
    }
}
#endif

/*! @brief calcul du vecteur normal a l'interface, aux sommets du maillage d'interface.
 *
 * Le tableau "normale" est efface et resize.
 *   La normale est la moyenne des normales des facettes voisines, ponderees par
 *   la surface de la facette.
 *   La norme du vecteur normal n'est pas unitaire !
 *   L'espace virtuel n'est pas a jour !
 *
 */
static void calculer_normale_sommets_interface(const Maillage_FT_Disc& maillage,
                                               DoubleTab& normale)
{
  const int nsom = maillage.nb_sommets();
  const int nfaces = maillage.nb_facettes();
  const int dim  = maillage.sommets().line_size();
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();
  const IntTab& facettes         = maillage.facettes();
  const int nsommets_faces = facettes.line_size();

  normale.resize(nsom, dim);
  normale = 0.;
  double n[3]= {0,0,0};

  for (int i = 0; i < nfaces; i++)
    {

      // Somme pour les faces reelles:
      if (maillage.facette_virtuelle(i))
        continue;

      const double surface = surface_facettes[i];
      for (int k = 0; k < dim; k++)
        n[k] = normale_facettes(i,k) * surface;

      for (int j = 0; j < nsommets_faces; j++)
        {
          const int sommet = facettes(i,j);
          for (int k = 0; k < dim; k++)
            normale(sommet, k) += n[k];
        }
    }

  // Sommer les contributions pour les sommets partages par plusieurs processeurs:
  maillage.desc_sommets().collecter_espace_virtuel(normale, MD_Vector_tools::EV_SOMME);
}

Implemente_instanciable_sans_constructeur_ni_destructeur(Transport_Interfaces_FT_Disc_interne,"Transport_Interfaces_FT_Disc_interne",Objet_U);

Entree& Transport_Interfaces_FT_Disc_interne::readOn(Entree& is)
{
  assert(0);
  return is;
}
Sortie& Transport_Interfaces_FT_Disc_interne::printOn(Sortie& os) const
{
  assert(0);
  return os;
}

int Transport_Interfaces_FT_Disc_interne::sauvegarder(Sortie& os) const
{
  // Il faut sauvegarder l'indicatrice_cache car elle ne peut pas toujours etre
  // correctement reconstruite a partir de l'interface (on tolere qu'il y ait
  // des inconsistances) :
  int bytes=0;
  bytes += indicatrice_cache.sauvegarder(os);
  int special, afaire;
  const int format_xyz = EcritureLectureSpecial::is_ecriture_special(special, afaire);
  if (format_xyz)
    {
      if (Process::je_suis_maitre())
        {
          os << indicatrice_cache_tag << finl;
        }
    }
  else
    {
      os << indicatrice_cache_tag << finl;
    }
  bytes += maillage_interface.sauvegarder(os);
  bytes += remaillage_interface_.sauvegarder(os); // Sauvegarde du temps du dernier remaillage...
  return bytes;
}

int Transport_Interfaces_FT_Disc_interne::reprendre(Entree& is)
{
  Nom ident, type;
  is >> ident >> type;
  if (! indicatrice_cache.non_nul())
    {
      // Le champ n'est pas discretise, on lit ceci pour sauter le bloc
      indicatrice_cache.typer(type);
    }
  indicatrice_cache.reprendre(is);
  is >> indicatrice_cache_tag;
  maillage_interface.reprendre(is);
  remaillage_interface_.reprendre(is);
  return 1;
}

/*! @brief constructeur par defaut
 *
 */
Transport_Interfaces_FT_Disc::Transport_Interfaces_FT_Disc()
{
  variables_internes_ = new Transport_Interfaces_FT_Disc_interne;
  temps_debut_=-1.e38;
  /*
    Noms& nom=champs_compris_.liste_noms_compris();
    nom.dimensionner(5);
    nom[0]="INDICATRICE_";
    nom[1]="INDICATRICE_CACHE_";
    nom[2]="VITESSE_FILTREE_";
    nom[3]="DISTANCE_INTERFACE_ELEM_";
    nom[4]="NORMALE_INTERFACE";
  */
  interpolation_repere_local_ = 0;
  force_.resize(dimension);
  moment_.resize((dimension==2?1:dimension));
}

/*! @brief le destructeur qui va avec
 *
 */
Transport_Interfaces_FT_Disc::~Transport_Interfaces_FT_Disc()
{
  delete variables_internes_;
}

Sortie& Transport_Interfaces_FT_Disc::printOn(Sortie& os) const
{
  return Transport_Interfaces_base::printOn(os);
}

Entree& Transport_Interfaces_FT_Disc::readOn(Entree& is)
{

  /*  Noms& nom=champs_compris_.liste_noms_compris();
      nom[0] += le_nom();
      nom[1] += le_nom();
      nom[2] += le_nom();
      nom[3] += le_nom();
  */

  Transport_Interfaces_base::readOn(is);
  if (suppression_interfaces_sous_domaine_!="??")
    {
      Cerr << "Name of subdomaine for interfaces deletion: " << suppression_interfaces_sous_domaine_ << finl;
      // Juste un test pour verifier que le nom existe:
      domaine_dis().domaine().ss_domaine(suppression_interfaces_sous_domaine_);
    }
  return is;
}

void Transport_Interfaces_FT_Disc::set_param(Param& param)
{
  Transport_Interfaces_base::set_param(param);
  param.ajouter_non_std("maillage",(this));
  param.ajouter("remaillage",&remaillage_interface());
  param.ajouter("collisions",&topologie_interface());
  param.ajouter_non_std("methode_transport",(this),Param::REQUIRED);
  param.ajouter_non_std("n_iterations_distance",(this));
  param.ajouter_non_std("iterations_correction_volume",(this)); // Former Keyword, Obsolete
  param.ajouter("VOFlike_correction_volume",&variables_internes_->VOFlike_correction_volume); // a flag: 0 or 1.
  param.ajouter("nb_lissage_correction_volume",&variables_internes_->nb_lissage_correction_volume);
  param.ajouter("nb_iterations_correction_volume",&variables_internes_->nb_iterations_correction_volume);
  param.ajouter_non_std("volume_impose_phase_1",(this));
  param.ajouter_non_std("methode_interpolation_v",(this));
  param.ajouter("parcours_interface",&variables_internes_->parcours_interface_);
  param.ajouter_non_std("injecteur_interfaces",(this));
  param.ajouter("suppression_sous_domaine",&suppression_interfaces_sous_domaine_);
  param.ajouter("sous_domaine_volume_impose",&variables_internes_->nom_domaine_volume_impose_);
  param.ajouter_flag("interpolation_repere_local", &interpolation_repere_local_);
  param.ajouter_non_std("interpolation_champ_face",(this));
  param.ajouter("n_iterations_interpolation_ibc",&variables_internes_->n_iterations_interpolation_ibc);
  param.ajouter_non_std("type_vitesse_imposee",(this));
  param.ajouter("nombre_facettes_retenues_par_cellule",&variables_internes_->nomb_fa7_accepted);
  param.ajouter("seuil_convergence_uzawa",&variables_internes_->seuil_uzawa) ;
  param.ajouter("nb_iteration_max_uzawa",&variables_internes_->nb_iter_uzawa) ;
  param.ajouter_non_std("distance_IBC_faces",(this));
  param.ajouter_non_std("distance_projete_faces",(this));
  param.ajouter("temps_debut",&temps_debut_);
  param.ajouter("vitesse_imposee_regularisee", &variables_internes_->vimp_regul) ;
  //param.ajouter("indic_faces_modifiee", &variables_internes_->indic_faces_modif) ;
  //param.ajouter_non_std("indic_faces_modifiee", (this)) ;
  param.ajouter_non_std("type_indic_faces", (this)) ;
}

int Transport_Interfaces_FT_Disc::lire_motcle_non_standard(const Motcle& un_mot, Entree& is)
{
  if (un_mot=="maillage")
    {
      //lecture des parametres du remaillage
      maillage_interface().lire_param_maillage(is);
      return 1;
    }
  else if (un_mot=="methode_transport")
    {
      Motcles motcles2(3);
      motcles2[0] = "vitesse_imposee";
      motcles2[1] = "loi_horaire";
      motcles2[2] = "vitesse_interpolee";
      Motcle motlu;
      is >> motlu;
      Cerr << un_mot << " " << motlu << finl;
      int rang = motcles2.search(motlu);
      switch(rang)
        {
        case 0:
          {
            if (Process::je_suis_maitre())
              Cerr << " Transport velocity imposed : reading " << Objet_U::dimension
                   << " analytical expressions vx(x,y,[z,]t) vy(...) [vz(...)]" << finl;
            // vitesse de deplacement de la phase 1 imposee
            variables_internes_->methode_transport =
              Transport_Interfaces_FT_Disc_interne::VITESSE_IMPOSEE;
            // lecture de Objet_U::dimension expressions analytiques
            //  contenant les composantes du champ de vitesse de transport
            //  (voir calculer_vitesse_imposee)
            // Ce sont des fonctions de x,y,[z,]t
            for (int i = 0; i < Objet_U::dimension; i++)
              {
                Nom& expression = variables_internes_->expression_vitesse_imposee[i];
                is >> expression;
                if (Process::je_suis_maitre())
                  Cerr << " Component " << i << " of the velocity : " << expression << finl;
              }
            break;
          }
        case 1:
          {
            Cerr << " Imposed motion by the time scheme : ";
            variables_internes_->methode_transport = Transport_Interfaces_FT_Disc_interne::LOI_HORAIRE;
            Nom nom_loi;
            is >> nom_loi;
            Cerr << nom_loi << finl;
            variables_internes_->loi_horaire_ = ref_cast(Loi_horaire, Interprete::objet(nom_loi));
            break;
          }
        case 2:
          {
            if (Process::je_suis_maitre())
              Cerr << " Transport velocity interpolated from the velocity \n"
                   << " of the following Navier_Stokes equation :" << finl;
            // vitesse de deplacement interpolee a partir du champ de vitesse fluide
            variables_internes_->methode_transport =
              Transport_Interfaces_FT_Disc_interne::VITESSE_INTERPOLEE;
            // lecture du nom de l'equation de navier_stokes :
            Motcle nom_eq;
            is >> nom_eq;
            Cerr << nom_eq << finl;
            variables_internes_->refequation_vitesse_transport =
              ref_cast(Probleme_FT_Disc_gen,
                       probleme_base_.valeur()).equation_hydraulique(nom_eq);
            break;
          }
        default:
          Cerr << "Transport_Interfaces_FT_Disc::lire\n"
               << "The options for methode_transport are :\n"
               << motcles2;
          Process::exit();
        }
      return 1;
    }

  else if (un_mot=="n_iterations_distance")
    {
      int n;
      is >> n;
      variables_internes_->n_iterations_distance = n;
      if (Process::je_suis_maitre())
        Cerr << "Iterations number of the smoothing tool \"distance_interface\": " << n << finl;
      return 1;
    }
  else if (un_mot=="iterations_correction_volume")
    {
      int n;
      is >> n;
      variables_internes_->iterations_correction_volume = n;
      if (n>0)
        variables_internes_->VOFlike_correction_volume = 1;
      variables_internes_->nb_lissage_correction_volume = n; // Historical behavior was to set it to the same value as the nb of iterations
      int m = remaillage_interface().get_nb_iter_bary_volume_seul();
      variables_internes_->nb_iterations_correction_volume = n;
      if (Process::je_suis_maitre())
        {
          Cerr << "Obsolete Keyword is read Iterations_correction_volume : " << n << finl;
          Cerr << "For future compatibility, it is recommended to switch to the new syntax : " << finl;
          Cerr << "   VOFlike_correction_volume 1 # a flag (0 or 1) for activation #" << finl;
          Cerr << "   nb_lissage_correction_volume " << n << " # Select 0 or N the number of smothing to apply to avoid potential spikes due to volume correction. #" << finl;
          Cerr << "   nb_iterations_correction_volume " << m << " # to iterate on the volume correction until seuil is reached. #" << finl;
          Cerr << "The value of is read nb_iterations_correction_volume from bloc remaillage at keyword: nb_iter_correction_volume" << finl;
        }
      return 1;
    }
  else if (un_mot=="volume_impose_phase_1")
    {
      double x;
      is >> x;
      variables_internes_->volume_impose_phase_1 = x;
      if (Process::je_suis_maitre())
        Cerr << "volume_impose_phase_1 : " << x << finl;
      return 1;
    }
  else if (un_mot=="methode_interpolation_v")
    {
      Motcles motcles2(2);
      motcles2[0] = "valeur_a_elem";
      motcles2[1] = "vdf_lineaire";
      Motcle motlu;
      is >> motlu;
      if (Process::je_suis_maitre())
        Cerr << un_mot << " " << motlu << finl;
      int rang = motcles2.search(motlu);
      switch(rang)
        {
        case 0:
          variables_internes_->methode_interpolation_v =
            Transport_Interfaces_FT_Disc_interne::VALEUR_A_ELEM;
          break;
        case 1:
          variables_internes_->methode_interpolation_v =
            Transport_Interfaces_FT_Disc_interne::VDF_LINEAIRE;
          break;
        default:
          Cerr << "Transport_Interfaces_FT_Disc::lire\n"
               << "The options for " << un_mot << " are :\n"
               << motcles2;
          exit();
        }
      return 1;
    }
  else if (un_mot=="parcours_interface")
    {
      is >> variables_internes_->parcours_interface_;
      return 1;
    }
  else if (un_mot=="injecteur_interfaces")
    {
      Nom nom_fichier;
      is >> nom_fichier;
      if (Process::je_suis_maitre())
        {
          Cerr << "Reading of the interfaces to inject in the file " << nom_fichier << finl;
          EFichier fic(nom_fichier);
          ArrOfDouble& temps = variables_internes_->injection_interfaces_temps_;
          ArrOfInt& phase = variables_internes_->injection_interfaces_phase_;
          Noms& expr = variables_internes_->injection_interfaces_expressions_;
          Nom une_expr;
          while(1)
            {
              double t;
              fic >> t;
              if (fic.eof())
                break;
              temps.append_array(t);
              int ph;
              fic >> ph;
              fic >> une_expr;
              phase.append_array(ph);
              expr.add(une_expr);
              Cerr << "Time=" << t << " phase=" << ph << " expression=" << une_expr << finl;
            }
        }
      envoyer_broadcast(variables_internes_->injection_interfaces_temps_, 0);
      envoyer_broadcast(variables_internes_->injection_interfaces_phase_, 0);
      envoyer_broadcast(variables_internes_->injection_interfaces_expressions_, 0);
      barrier();
      return 1;
    }
  else if (un_mot=="interpolation_champ_face")
    {
      Motcles motcles2(2);
      motcles2[0] = "base";
      motcles2[1] = "lineaire";
      Motcle motbis;
      is >> motbis;
      int rang = motcles2.search(motbis);
      switch(rang)
        {
        case 0:
          variables_internes_->interpolation_champ_face =
            Transport_Interfaces_FT_Disc_interne::BASE;
          break;
        case 1:
          {
            variables_internes_->interpolation_champ_face =
              Transport_Interfaces_FT_Disc_interne::LINEAIRE;

            Motcles mots;
            mots.add("vitesse_fluide_explicite");        // vitesse utilisee pour calculer l'interpolation du schema lineaire
            mots.add("extrapolation_diphasique_solide"); // le gradient est etendu dans les mailles diphasiques solides
            mots.add("extrapolation_solide");            // le gradient est etendu dans les mailles solides
            mots.add("distance_projete_face");           // utilisation de la distance reelle entre le barycentre de la face
            // et le point projete (au lieu d'utiliser dSigma)

            Motcle mot2;
            is >> mot2;
            Motcle accouverte = "{" , accfermee = "}" ;
            if (mot2 == accouverte)
              {
                is >> mot2;

                while (mot2 != accfermee)
                  {
                    int rang2 = mots.search(mot2);
                    switch(rang2)
                      {
                      case 0:
                        {
                          variables_internes_->vf_explicite = 1 ;
                          Cerr << "Lecture du type de vitesse fluide pour le calcul de v_imp  : vf_explicite " << finl;
                          break;
                        }
                      case 1:
                        {
                          variables_internes_->is_extra_diphasique_solide = 1 ;
                          variables_internes_->is_extra_solide = 0 ;
                          Cerr << "Lecture du type d'extrapolation considere : diphasique_solide " << finl;
                          break;
                        }
                      case 2:
                        {
                          variables_internes_->is_extra_solide = 1 ;
                          variables_internes_->is_extra_diphasique_solide = 1 ;
                          Cerr << "Lecture du type d'extrapolation considere : solide " << finl;
                          break;
                        }
                      case 3:
                        {
                          variables_internes_->is_distance_projete_face = 1 ;
                          Cerr << "Lecture du type de distance consideree dans l'interpolation : distance reelle entre projete et face " << finl;
                          break;
                        }
                      }
                    is >> mot2;
                  }
              }
            else
              {
                Cerr << "Erreur a la lecture des parametres de l'interpolation lineaire : " << finl;
                Cerr << "On attendait : " << accouverte << " , on a eu " << mot2 << finl;
                exit();
              }
            break;
          }
        default:
          Cerr << "Transport_Interfaces_FT_Disc::lire\n"
               << " les options de interpolation_champ_face sont :\n"
               << motcles2;
          exit();
        }
    }
  else if (un_mot=="type_vitesse_imposee")
    {
      Motcles motcles2(2);
      motcles2[0] = "uniforme";
      motcles2[1] = "analytique";
      Motcle motbis;
      is >> motbis;
      int rang = motcles2.search(motbis);
      switch(rang)
        {
        case 0:
          variables_internes_-> type_vitesse_imposee = Transport_Interfaces_FT_Disc_interne::UNIFORME;
          break;
        case 1:
          variables_internes_-> type_vitesse_imposee = Transport_Interfaces_FT_Disc_interne::ANALYTIQUE;
          break;
        default:
          Cerr << "Transport_Interfaces_FT_Disc::lire\n"
               << " les options de type_vitesse_imposee sont :\n"
               << motcles2;
          exit();
        }
    }
  else if (un_mot=="distance_IBC_faces")
    {
      Motcles motcles2(2);
      motcles2[0] = "initiale";
      motcles2[1] = "modifiee";
      Motcle motbis;
      is >> motbis;

      int rang = motcles2.search(motbis);
      switch(rang)
        {
        case 0:
          variables_internes_-> type_distance_calculee = Transport_Interfaces_FT_Disc_interne::DIST_INITIALE;
          break;
        case 1:
          {
            variables_internes_-> type_distance_calculee = Transport_Interfaces_FT_Disc_interne::DIST_MODIFIEE;
            break;
          }
        default:
          Cerr << "Transport_Interfaces_FT_Disc::lire\n"
               << " les options de distance_IBC_faces sont :\n"
               << motcles2;
          exit();
        }
    }
  else if (un_mot=="distance_projete_faces")
    {
      Motcles motcles2(3);
      motcles2[0] = "initiale";
      motcles2[1] = "modifiee";
      motcles2[2] = "simplifiee";
      Motcle motbis;
      is >> motbis;

      int rang = motcles2.search(motbis);
      switch(rang)
        {
        case 0:
          variables_internes_-> type_projete_calcule = Transport_Interfaces_FT_Disc_interne::PROJETE_INITIAL;
          break;
        case 1:
          variables_internes_-> type_projete_calcule = Transport_Interfaces_FT_Disc_interne::PROJETE_MODIFIE;
          break;
        case 2:
          variables_internes_-> type_projete_calcule = Transport_Interfaces_FT_Disc_interne::PROJETE_SIMPLIFIE;
          break;
        default:
          Cerr << "Transport_Interfaces_FT_Disc::lire\n"
               << " les options de distance_projete_faces sont :\n"
               << motcles2;
          Process::exit();
        }
    }
  else if (un_mot=="type_indic_faces")
    {
      Motcles motcles2(3);
      motcles2[0] = "standard";
      motcles2[1] = "modifiee";
      motcles2[2] = "ai_based";
      Motcle motlu;
      is >> motlu;
      Cerr << un_mot << " " << motlu << finl;
      int rang = motcles2.search(motlu);
      switch(rang)
        {
        case Transport_Interfaces_FT_Disc_interne::STANDARD:
          {
            variables_internes_->type_indic_faces_ = Transport_Interfaces_FT_Disc_interne::STANDARD;
            if (Process::je_suis_maitre())
              Cerr << " Standard interpolation of indicatrice to faces." << finl;
            return 1;
          }
        case Transport_Interfaces_FT_Disc_interne::MODIFIEE:
          {
            if (Process::je_suis_maitre())
              Cerr << " Transport velocity imposed : reading " << Objet_U::dimension
                   << " analytical expressions vx(x,y,[z,]t) vy(...) [vz(...)]" << finl;

            variables_internes_->type_indic_faces_ = Transport_Interfaces_FT_Disc_interne::MODIFIEE;
            // modified_indic_faces_* : Valeurs de la position et de l epaisseur de l indicatrice modifiee,
            // exprimees en multiples de taille de maille. Cette nouvelle indicatrice est calculee
            // a partir de la distance a l interface. La position represente une iso-ligne
            // de distance a l interface, et l epaisseur, le domaine de variation lineaire vis a vis de la distance,
            // qui fait passer l indicatrice de 0 a 1 (le domaine est centre sur l iso-ligne representee par la position).
            // Par defaut :
            variables_internes_->modified_indic_faces_position = 0. ;
            variables_internes_->modified_indic_faces_thickness= 1. ;
            Motcles mots;
            mots.add("position");
            mots.add("thickness");

            Motcle mot;
            is >> mot;
            Motcle accouverte = "{" , accfermee = "}" ;
            if (mot == accouverte)
              {
                is >> mot;

                while (mot != accfermee)
                  {
                    int rang2 = mots.search(mot);
                    switch(rang2)
                      {
                      case 0:
                        {
                          is >> variables_internes_->modified_indic_faces_position ;
                          Cerr << "Lecture de la position de l interface modifiee " << finl;
                          break;
                        }
                      case 1:
                        {
                          is >> variables_internes_->modified_indic_faces_thickness ;
                          if ( variables_internes_->modified_indic_faces_thickness < 0.)
                            {
                              Cerr << "L epaisseur de l interface doit etre positive ou nulle!!" << finl;
                              Process::exit();
                            }
                          Cerr << "Lecture de l epaisseur de l interface modifiee " << finl;
                          break;
                        }
                      }
                    is >> mot;
                  }
              }
            else
              {
                Cerr << "Erreur a la lecture des parametres de l'indicatrice modifiee " << finl;
                Cerr << "On attendait : " << accouverte << finl;
                exit();
              }
            if (Process::je_suis_maitre())
              {
                Cerr << "L indicatrice face sera calculee a partir de la distance. Position : d=" << variables_internes_->modified_indic_faces_position << "h ; Epaisseur : "<< variables_internes_->modified_indic_faces_thickness << "h" <<finl;
              }
            return 1;
          }
        case Transport_Interfaces_FT_Disc_interne::AI_BASED:
          {
            variables_internes_->type_indic_faces_ = Transport_Interfaces_FT_Disc_interne::AI_BASED;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces is based on the interfacial area"
                   << " and on the normal to the interface." << finl;
            return 1;
          }
        default:
          Cerr << "Transport_Interfaces_FT_Disc::lire\n"
               << " les options de type_indic_faces sont :\n"
               << motcles2;
          Process::exit();
        }
    }
  else
    return Transport_Interfaces_base::lire_motcle_non_standard(un_mot,is);
  return 1;
}

/*! @brief Methode appelee par Equation_base::readOn On verifie que toutes les cl sont de type Paroi_FT_disc.
 *
 *   Fait exit() si erreur.
 *
 */
int Transport_Interfaces_FT_Disc::verif_Cl() const
{
  const Conds_lim& les_cl = le_dom_Cl_dis.valeur().les_conditions_limites();
  const int n = les_cl.size();
  int i;
  for (i = 0; i < n; i++)
    {
      if (! sub_type(Paroi_FT_disc, les_cl[i].valeur()))
        break;
    }
  if (i < n)
    {
      Cerr << "Error in Transport_Interfaces_FT_Disc::verif_Cl():\n"
           << " Boundary conditions for Transport_Interfaces_FT_Disc must be\n"
           << " of type Paroi_FT_disc (example : \"Paroi_FT_disc symetrie\")" << finl;
      exit();
    }
  return 1;
}

static void fct_tri_sommet_fa7(const int* in, int* out)
{
  const int dim = Objet_U::dimension;
  out[0]=in[0];
  out[1]=in[1];
  for (int i=2; i<dim; i++)
    out[i]=in[i];

  if(out[1]<out[0])
    {
      const int temp=out[1];
      out[1]=out[0];
      out[0]=temp;
    }
  if(dim==3)
    {
      out[2]=in[2];
      if(out[2]<out[1])
        {
          const int temp=out[2];
          out[2]=out[1];
          out[1]=temp;
        }
      if(out[1]<out[0])
        {
          const int temp=out[1];
          out[1]=out[0];
          out[0]=temp;
        }
    }
  else
    out[2]=-123;
}

static True_int fct_tri_facettes(const void *pt1, const void *pt2)
{
  const int *a = (const int *) pt1;
  const int *b = (const int *) pt2;

  int i, x = 0;
  const int dim = Objet_U::dimension;

  int A[3],B[3];
  fct_tri_sommet_fa7(a,A);
  fct_tri_sommet_fa7(b,B);

  for (i = 0; i < dim; i++)
    {
      x = A[i] - B[i];
      if (x != 0)
        break;
    }
  return x;
}

void Transport_Interfaces_FT_Disc::lire_maillage_ft_cao(Entree& is)
{
  Domaine dom;
  REF(Domaine) ref_dom;
  Nom filename;
  ArrOfInt phase_specifiee;
  DoubleTab points;
  int default_phase = -1;
  Nom lata_file;
  Nom domain_name;
  int reverse_normal = 0;

  Motcle motlu;
  is >> motlu;
  if (motlu != "{")
    {
      Cerr << "Error in Transport_Interfaces_FT_Disc::lire_cond_init, expected { after fichier_geom" << finl;
      exit();
    }

  Motcles motcles(6);
  motcles[0] = "fichier_geom";
  motcles[1] = "point_phase";
  motcles[2] = "default_phase";
  motcles[3] = "lata_dump";
  motcles[4] = "nom_domaine";
  motcles[5] = "reverse_normal";
  do
    {
      is >> motlu;
      if (motlu == "}")
        break;
      const int rang = motcles.search(motlu);
      switch(rang)
        {
        case 0:
          is >> filename;
          Cerr << " Reading .geom file in following file: " << filename << finl;
          break;
        case 1:
          {
            const int i = phase_specifiee.size_array();
            phase_specifiee.resize_array(i+1);
            points.resize(i+1, dimension);
            is >> phase_specifiee[i];
            Cerr << " Reading point in phase " << phase_specifiee[i];
            if (phase_specifiee[i] != 0 && phase_specifiee[i] != 1)
              {
                Cerr << " Error reading point_phase : expected 0 or 1" << finl;
                exit();
              }
            for (int j = 0; j < dimension; j++)
              {
                is >> points(i, j);
                Cerr << " " << points(i,j);
              }
            Cerr << finl;
            break;
          }
        case 2:
          is >> default_phase;
          if (default_phase != 0 && default_phase != 1)
            {
              Cerr << " Error reading default_phase : expected 0 or 1" << finl;
              exit();
            }
          Cerr << " Default phase : " << default_phase << finl;
          break;
        case 3:
          is >> lata_file;
          Cerr << " Dumping connex components and indicator function in lata file : "
               << lata_file << finl;
          break;
        case 4:
          is >> domain_name;
          Cerr << " Reading interface data in domain: " << domain_name << finl;
          break;
        case 5:
          reverse_normal = 1;
          Cerr << " Reverse_normal : swap nodes to reverse the normal vector of the surface mesh" << finl;
          break;
        default:
          Cerr << " Unknown keyword " << motlu
               << "\n Known keywords are " << motcles << finl;
          exit();
        }
    }
  while (1);

  if (filename != "??" && domain_name != "??")
    {
      Cerr << "Error: you specified both a .geom file name AND a domain name" << finl;
      exit();
    }
  else if (filename == "??" && domain_name == "??")
    {
      Cerr << "Error: you must specify a FICHIER_GEOM or a NOM_DOMAINE" << finl;
      exit();
    }
  if (filename != "??")
    {
      Cerr << "Reading geometry file and building interface mesh" << finl;
      if (is_a_binary_file(filename))
        {
          EFichierBin fic(filename);
          fic >> dom;
        }
      else
        {
          LecFicDiffuse fic(filename);
          fic >> dom;
        }
      ref_dom = dom;
    }
  else
    {
      Cerr << "Reading interface in existing domain: " << domain_name << finl;
      if (!sub_type(Domaine, Interprete::objet(domain_name)))
        {
          Cerr << "Error : object " << domain_name << " is not a domain" << finl;
          exit();
        }
      ref_dom = ref_cast(Domaine, Interprete::objet(domain_name));
    }
  if (reverse_normal)
    {
      Cerr << "Reversing mesh normal vectors" << finl;
      Domaine& domaine = ref_dom.valeur();
      IntTab& elems = domaine.les_elems();
      const int nb_elem = elems.dimension(0);
      const int nb_som_elem = elems.line_size();
      if (nb_som_elem != 2 && nb_som_elem != 3)
        {
          Cerr << "Error: mesh has wrong dimension (must be segments in 2D, triangles in 3D)" << finl;
          exit();
        }
      const int j0 = nb_som_elem - 2;
      for (int i = 0; i < nb_elem; i++)
        {
          const int s0 = elems(i, j0);
          const int s1 = elems(i, j0+1);
          elems(i, j0) = s1;
          elems(i, j0+1) = s0;
        }
    }

  //MODIF SM 01/12/08
  //Effacer les fa7 doublons dans la cas d'un reprise d'un .MED parallele
  // Verification qu'il n'existe pas deux fois la meme facette
  {
    Scatter::uninit_sequential_domain(ref_dom.valeur());
    Domaine& domaine = ref_dom.valeur();
    IntTab& fa7 = domaine.les_elems();

    // tri du tableau
    int * data = fa7.addr();
    const int nb_facettes = fa7.dimension(0);
    assert(Objet_U::dimension == fa7.line_size());
    qsort(data, nb_facettes, fa7.line_size()*sizeof(int),
          fct_tri_facettes);

    // recherche et suppression des doublons
    int i;
    int count = 0;
    const int nb_som_facettes = Objet_U::dimension;
    for (i = 1; i < nb_facettes; i++)
      {
        if(fct_tri_facettes(&fa7(i,0),&fa7(i-1,0))!=0)
          {
            count++;
            for(int j=0; j<nb_som_facettes; j++)
              fa7(count,j)=fa7(i,j);

          }
      }
    fa7.resize(count+1,dimension);
    Scatter::init_sequential_domain(ref_dom.valeur());
    Cerr<<"End correction SM 12/08 the removed faces number is: "<<  nb_facettes <<" - "<< count+1 << " = " << nb_facettes-count-1 << finl;
  }
  //Fin modif SM 01/12/08


  Cerr << "Building interface mesh" << finl;
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
  Maillage_FT_Disc& maillage = maillage_interface();
  Sauvegarde_Reprise_Maillage_FT::lire_xyz(maillage, &domaine_vf, 0, &(ref_dom.valeur()));
  Cerr << "Extracting connex components and assigning indicator function." << finl;

  maillage.parcourir_maillage();
  const int nb_elem = domaine_vf.domaine().nb_elem();
  IntVect num_compo;
  domaine_vf.domaine().creer_tableau_elements(num_compo);
  // On marque avec -1 les elements traverses par une interface:
  const ArrOfInt& index_elem = maillage.intersections_elem_facettes().index_elem();
  {
    for (int i = 0; i < nb_elem; i++)
      {
        int n = 0;
        if (index_elem[i] >= 0)
          n = -1; // L'element est traverse par une interface
        num_compo[i] = n;
      }
  }
  num_compo.echange_espace_virtuel();
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const IntTab& faces_elem = domaine_vf.face_voisins();
  const int nb_local_connex_components = search_connex_components_local(elem_faces, faces_elem, num_compo);

  const int nb_connex_components = compute_global_connex_components(num_compo, nb_local_connex_components);

  //const int nb_connex_components = nb_local_connex_components;
  Cerr << " found " << nb_connex_components << " connex components" << finl;
  // Calcul de la phase a associer a chaque composante:
  ArrOfInt phase_of_component(nb_connex_components);
  phase_of_component= default_phase;
  // Pour chaque point ou la phase a ete specifiee, trouver l'element dans
  // lequel il se trouve, puis sa composante connexe et associer cette phase
  // a la composante connexe:
  const int nb_pts = phase_specifiee.size_array();
  if (nb_pts!=0)
    {
      ArrOfInt elem_points;
      domaine_vf.domaine().chercher_elements(points, elem_points);
      for (int i = 0; i < nb_pts; i++)
        {
          int composante_connexe = -2;
          if (elem_points[i] >= 0)
            composante_connexe = num_compo[elem_points[i]];
          // En general un seul processeur trouve le point chez lui, il distribue sa
          // trouvaille aux autres :
          composante_connexe = (int) mp_max(composante_connexe);
          if (composante_connexe == -2)
            {
              Cerr << "Error : point " << i << " is not in the domain" << finl;
              barrier();
              exit();
            }
          else if (composante_connexe == -1)
            {
              Cerr << "Error : point " << i << " is in connex_component -1"
                   << " (this point is in a mesh cell containing the interface)."
                   << " Use the dump_lata keyword to find a point away from the interface" << finl;
              barrier();
              exit();
            }
          else
            {
              const int phase = phase_specifiee[i];
              Cerr << "Assigning phase " << phase << " to component " << composante_connexe << finl;
              phase_of_component[composante_connexe] = phase;
            }
        }
    }
  // Met a jour l'indicatrice de phase:
  DoubleTab& indic = variables_internes_->indicatrice_cache.valeur().valeurs();
  for (int i = 0; i < nb_elem; i++)
    {
      const int compo = num_compo[i];
      if (compo >= 0)
        {
          const int phase = phase_of_component[compo];
          indic(i) = phase;
        }
    }
  get_update_indicatrice();
  // Postraitement des composantes connexes et de l'indicatrice
  if (lata_file != "??")
    {
      Cerr << "Writing lata file" << finl;
      Format_Post_Lata lata;
      const Domaine& un_dom = domaine_vf.domaine();
      constexpr double TEMPS = 0.;
      constexpr int FIRST_POST = 1;
      lata.initialize(lata_file, Format_Post_Lata::BINAIRE, "SIMPLE");
      lata.ecrire_entete(TEMPS, 0 /*reprise*/, FIRST_POST);
      lata.ecrire_domaine(un_dom, FIRST_POST);
      lata.ecrire_temps(TEMPS);
      DoubleTab data(nb_elem);
      for (int i = 0; i < nb_elem; i++)
        data(i) = num_compo(i);
      Noms unites;
      unites.add("-");
      Noms noms_compo;
      noms_compo.add("");
      Nom nom_dom(un_dom.le_nom());

      Nom nom_champ("connex_component");
      lata.ecrire_champ(un_dom, unites, noms_compo, 1, TEMPS, nom_champ, nom_dom, "elem", "scalar",
                        data);
      nom_champ = "indicatrice";
      lata.ecrire_champ(un_dom, unites, noms_compo, 1, TEMPS, nom_champ, nom_dom, "elem", "scalar",
                        indic);
      nom_champ = "distance";
      lata.ecrire_champ(un_dom, unites, noms_compo, 1, TEMPS, nom_champ, nom_dom, "elem", "scalar",
                        get_update_distance_interface().valeurs());
    }

  if (phase_of_component.size_array() > 0 && min_array(phase_of_component) < 0)
    {
      Cerr << "Error: some connex components of the indicator function were not initialized\n"
           << " you must either specify more points or specify a default_phase.\n"
           << " With the lata_dump keyword, you can watch connex components and find where\n"
           << " points should be added." << finl;
      exit();
    }
}


/*! @brief Lecture des conditions initiales.
 *
 * On s'attend a trouver ceci : { fonction EXPRESSION }
 *  ou expression depend de x, y et z et sera interpretee par le
 *  parser de TRUST. L'expression est envoyee a Marching_Cubes
 *  pour construire l'interface.
 *  Voir aussi parser et Marching_Cubes.
 *
 */
Entree& Transport_Interfaces_FT_Disc::lire_cond_init(Entree& is)
{
  if (Process::je_suis_maitre())
    Cerr << "Reading initial condition" << finl;
  Motcles motcles(4);
  motcles[0] = "fonction";
  motcles[1] = "fichier_geom";
  motcles[2] = "fonction_ignorer_collision";
  motcles[3] = "reprise";
  Motcle motlu;
  is >> motlu;
  if (motlu != "{")
    {
      Cerr << "Erreur dans Transport_Interfaces_FT_Disc::lire_cond_init\n";
      Cerr << " On attendait {\n On a trouve " << motlu << finl;
      exit();
    }
  const Motcle virgule = ",";
  int init = 1;
  do
    {
      Maillage_FT_Disc maillage_tmp;
      maillage_tmp.associer_equation_transport(*this);
      //boucle sur la lecture des conditions intiales
      is >> motlu;
      int rang = motcles.search(motlu);
      switch (rang)
        {
        case 0:
        case 2:
          {
            Maillage_FT_Disc::AjoutPhase phase = Maillage_FT_Disc::AJOUTE_TOUTES_PHASES;
            if (init)
              {
                init = 0;
              }
            else
              {
                //a partir de la 2e interface ajoutee, il faut indiquer quelle phase on ajoute
                Motcle nom_phase;
                is >> nom_phase;
                if (nom_phase=="ajout_phase1")
                  {
                    phase = Maillage_FT_Disc::AJOUTE_PHASE1;
                  }
                else if (nom_phase=="ajout_phase0")
                  {
                    phase = Maillage_FT_Disc::AJOUTE_PHASE0;
                  }
                else
                  {
                    Cerr << "Transport_Interfaces_FT::lire_cond_init :\n";
                    Cerr << " The keyword " << nom_phase << " is not understood.\n";
                    Cerr << " Allowed keywords are : ajout_phase0 or ajout_phase1" << finl;
                    assert(0);
                    exit();
                  }
              }
            Nom expression;
            is >> expression;;
            if (probleme().reprise_effectuee())
              {
                Cerr << " Interface not build since a restarting is expected." << finl;
              }
            else
              {
                Cerr << " Interface construction : " << expression << finl;
                // Construction de l'interface comme l'isovaleur zero de la fonction
                // La valeur de la fonction "expression" aux sommets est temporairement stockee
                // dans distance_interface_sommets qui a la bonne structure d'espace virtuel
                // (items communs correctement initialises).
                int ignorer_collision = (rang==2);
                const int ok = marching_cubes().construire_iso(expression, 0., maillage_tmp,
                                                               variables_internes_->indicatrice_cache.valeur().valeurs(),
                                                               phase,
                                                               variables_internes_->distance_interface_sommets,
                                                               ignorer_collision);
                if (!ok)
                  {
                    if (ignorer_collision)
                      Cerr << "Warning: the interface is colliding with an existing interface (ignorer_collision=1)." << finl;
                    else
                      {
                        Cerr << "Error:  the interface is colliding with an existing interface." << finl;
                        barrier();
                        exit();
                      }
                  }
                maillage_interface().ajouter_maillage(maillage_tmp);
              }
            break;
          }
        case 1:
          {
            lire_maillage_ft_cao(is);
            break;
          }
        case 3:
          {
            Nom file;
            is >> file;
            if (!file.finit_par(".xyz"))
              {
                Cerr << "A .xyz file is waited to restart the calculation for the interface equation." << finl;
                Process::exit();
              }
            int mode_lec_sa = EcritureLectureSpecial::mode_lec;
            EcritureLectureSpecial::mode_lec=1;
            EFichierBin is2(file);
            Motcle tmp;
            int version;
            is2 >> tmp >> version; // FORMAT_SAUVEGARDE: VERSION
            reprendre(is2);
            EcritureLectureSpecial::mode_lec = mode_lec_sa;
            break;
          }
        default:
          Cerr << "Transport_Interfaces_FT::lire_cond_init :\n";
          Cerr << " The keyword " << motlu << " is not understood here.\n";
          Cerr << " Allowed keywords are :\n" << motcles << finl;
          exit();
        }
      is >> motlu;
      Cerr<<"lectureF "<<motlu<<finl;
    }
  while (motlu==virgule);

  if (motlu != "}")
    {
      Cerr << "Error for method Transport_Interfaces_FT_Disc::lire_cond_init\n";
      Cerr << " A } was expected\n It has been found " << motlu << finl;
      exit();
    }

  return is;
}

void Transport_Interfaces_FT_Disc::associer_pb_base(const Probleme_base& un_probleme)
{
  assert(! probleme_base_.non_nul());
  if ( (! sub_type(Probleme_FT_Disc_gen, un_probleme))
       && !(un_probleme.is_dilatable()))
    {
      Cerr << "Error for method Transport_Interfaces_FT_Disc::associer_pb_base\n";
      Cerr << " The Transport_Interfaces_FT_Disc equation must be associated to a problem\n";
      Cerr << " of type Probleme_FT_Disc_gen" << finl;
      Cerr << " or type Pb_Thermohydraulique_Especes_QC" << finl;
      Cerr << " or type Pb_Thermohydraulique_Especes_Turbulent_QC" << finl;
      Process::exit();
    }
  probleme_base_ = un_probleme;
  Equation_base::associer_pb_base(un_probleme);
}

/*! @brief Discretisation des champs: - indicatrice_ : champ scalaire discretise aux elements
 *
 *  - typage du maillage et de l'algorithme marching cubes
 *
 */
void Transport_Interfaces_FT_Disc::discretiser(void)
{
  // Le nom des differents champs est un identifiant (indicatrice, vitesse, ...)
  // suivi de "_nom_de_l_equation" :
  Nom suffix("_");
  suffix += le_nom();

  const Discretisation_base& dis = discretisation();
  const double temps = schema_temps().temps_courant();
  const Domaine_dis_base& mon_dom_dis = domaine_dis().valeur();
  const int nb_valeurs_temps = schema_temps().nb_valeurs_temporelles();

  Nom fieldname;
  fieldname = "INDICATRICE";
  fieldname += suffix;
  dis.discretiser_champ("champ_elem", mon_dom_dis,
                        fieldname, "-",
                        1 /* composantes */, nb_valeurs_temps,
                        temps,
                        indicatrice_);
  indicatrice_.associer_eqn(*this);
  //Nouvelle formulation
  champs_compris_.ajoute_champ(indicatrice_);
  //champs_compris_.liste_noms_compris()[0]+le_nom();

  fieldname = "INDICATRICE_CACHE";
  fieldname += suffix;
  dis.discretiser_champ("champ_elem", mon_dom_dis,
                        fieldname, "-",
                        1 /* composantes */, 1 /* valeur temporelle */,
                        temps,
                        variables_internes_->indicatrice_cache);
  variables_internes_->indicatrice_cache.associer_eqn(*this);
  champs_compris_.ajoute_champ(variables_internes_->indicatrice_cache);
  //champs_compris_.liste_noms_compris()[1]+le_nom();

  fieldname = "INDICATRICE_FACES";
  fieldname += suffix;
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        fieldname, "-",
                        1 /* composantes */, nb_valeurs_temps,
                        temps,
                        indicatrice_faces_);
  indicatrice_faces_.associer_eqn(*this);
  champs_compris_.ajoute_champ(indicatrice_faces_);

  fieldname = "VITESSE_FILTREE";
  fieldname += suffix;
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        fieldname, "m/s",
                        Objet_U::dimension /* composantes */, 1, /* valeur temporelle */
                        temps,
                        variables_internes_->vitesse_filtree);
  variables_internes_->vitesse_filtree.associer_eqn(*this);
  champs_compris_.ajoute_champ(variables_internes_->vitesse_filtree);
  //champs_compris_.liste_noms_compris()[2]+le_nom();

  fieldname = "FLUX_TMP";
  fieldname += suffix;
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        fieldname, "m/s",
                        1 /* composantes */,
                        temps,
                        variables_internes_->tmp_flux);
  champs_compris_.ajoute_champ(variables_internes_->tmp_flux);
  //champs_compris_.liste_noms_compris()[2]+le_nom();

  fieldname = "INDEX_ELEM";
  dis.discretiser_champ("champ_elem", mon_dom_dis,
                        fieldname, "m",
                        Objet_U::dimension /* composantes */,
                        temps,
                        variables_internes_->index_element);
  champs_compris_.ajoute_champ(variables_internes_->index_element);

  fieldname = "NELEM_PAR_DIRECTION";
  dis.discretiser_champ("champ_elem", mon_dom_dis,
                        fieldname, "m",
                        1 /* composantes */,
                        temps,
                        variables_internes_->nelem_par_direction);
  champs_compris_.ajoute_champ(variables_internes_->nelem_par_direction);

  fieldname = "DISTANCE_INTERFACE_ELEM";
  fieldname += suffix;
  dis.discretiser_champ("champ_elem", mon_dom_dis,
                        fieldname, "m",
                        1 /* composantes */,
                        temps,
                        variables_internes_->distance_interface);
  champs_compris_.ajoute_champ(variables_internes_->distance_interface);

  fieldname = "DISTANCE_INTERFACE_FACE";
  fieldname += suffix;
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        fieldname, "m",
                        1 /* composantes */,
                        temps,
                        variables_internes_->distance_interface_faces);
  champs_compris_.ajoute_champ(variables_internes_->distance_interface_faces);

  fieldname = "DISTANCE_INTERFACE_FACE_COR";
  fieldname += suffix;
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        fieldname, "m",
                        1 /* composantes */,
                        temps,
                        variables_internes_->distance_interface_faces_corrigee);
  champs_compris_.ajoute_champ(variables_internes_->distance_interface_faces_corrigee);

  fieldname = "DISTANCE_INTERFACE_FACE_DIF";
  fieldname += suffix;
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        fieldname, "m",
                        1 /* composantes */,
                        temps,
                        variables_internes_->distance_interface_faces_difference);
  champs_compris_.ajoute_champ(variables_internes_->distance_interface_faces_difference);

  fieldname = "VITESSE_IMP_INTERP";
  fieldname += suffix;
  dis.discretiser_champ("vitesse", mon_dom_dis,
                        fieldname, "m",
                        -1 /* composantes */,
                        temps,
                        vitesse_imp_interp_);
  champs_compris_.ajoute_champ(vitesse_imp_interp_);

  //champs_compris_.liste_noms_compris()[3]+le_nom();

  // Construction de la structure du tableau avec l'espace virtuel:
  {
    DoubleTab& d = variables_internes_->distance_interface_sommets;
    const Domaine& dom = domaine_dis().domaine();
    d.resize(0);
    dom.creer_tableau_sommets(d);
  }
  variables_internes_->distance_interface_sommets = 0.;

  fieldname = "NORMALE_INTERFACE";
  fieldname += suffix;
  dis.discretiser_champ("champ_elem", mon_dom_dis,
                        fieldname, "-",
                        Objet_U::dimension /* composantes */,
                        temps,
                        variables_internes_->normale_interface);
  champs_compris_.ajoute_champ(variables_internes_->normale_interface);
  //champs_compris_.liste_noms_compris()[4]+le_nom();

  //GB : Ajout du post-traitement de l'aire interfaciale par cellule :
  fieldname = "SURFACE_INTERFACE";
  fieldname += suffix;
  dis.discretiser_champ("champ_elem", mon_dom_dis,
                        fieldname, "m2",
                        1 /* composantes */,
                        temps,
                        variables_internes_->surface_interface);
  champs_compris_.ajoute_champ(variables_internes_->surface_interface);
  //GB : Fin de l'ajout.

  if (calculate_time_derivative())
    {
      //TF : Ajout de la gestion de la derivee en temps
      fieldname = "TIME_DERIVATIVE_INTERFACE";
      fieldname += suffix;
      dis.discretiser_champ(Motcle("champ_elem"), mon_dom_dis,
                            fieldname, Nom(""),
                            1 /* composantes */,
                            schema_temps().nb_valeurs_temporelles(),
                            temps,
                            derivee_en_temps());
      //TF : Fin de l'ajout
    }

  maillage_interface().changer_temps(0.);
  {
    // On modifie le domaine ici => on a besoin d'une reference non constante
    Domaine_dis_base& le_dom_dis2 = domaine_dis().valeur();
    le_dom_dis2.domaine().construire_elem_virt_pe_num();
  }
  {
    Domaine_VF& domaine_vf = ref_cast_non_const(Domaine_VF, mon_dom_dis);
    domaine_vf.construire_face_virt_pe_num();
    variables_internes_->connectivite_frontieres_.associer_domaine_vf(domaine_vf);
  }
  {
    Parcours_interface& parcours = variables_internes_->parcours_interface_;
    parcours.associer_domaine_dis(domaine_dis());
    parcours.associer_connectivite_frontieres(connectivite_frontieres());
  }
  {
    const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, mon_dom_dis);
    marching_cubes().associer_domaine_vf(domaine_vf);
  }
  maillage_interface().associer_equation_transport(*this);
  remaillage_interface().associer_domaine(domaine_dis());

  variables_internes_->algorithmes_transport_.typer("Algorithmes_Transport_FT_Disc");
  // On n'appelle pas Equation_base::discretiser car on ne veut pas
  // de solveur masse.
  discretisation().domaine_Cl_dis(domaine_dis(), le_dom_Cl_dis);
  le_dom_Cl_dis->associer_eqn(*this);
  le_dom_Cl_dis->associer_inconnue(inconnue());
}

/*! @brief Remaillage de l'interface : - amelioration petites et grandes facettes,
 *
 *    - barycentrage,
 *    - gestion des coalescences-fragmentations.
 *
 */
void Transport_Interfaces_FT_Disc::remailler_interface(void)
{
  Journal() << "Transport_Interfaces_FT_Disc::remailler_interface " << le_nom() << finl;
  const double temps = schema_temps().temps_courant();
  Maillage_FT_Disc& maillage = maillage_interface();
  // Ce champ sera mis a jour si on realise un remaillage global...
  Champ_base& indicatrice = variables_internes_->indicatrice_cache.valeur();
  Remaillage_FT& algo_remaillage_local = remaillage_interface();
  Topologie_Maillage_FT& algo_topologie = topologie_interface();
  algo_topologie.remailler_interface(temps, maillage, indicatrice, algo_remaillage_local);
}

int Transport_Interfaces_FT_Disc::preparer_calcul(void)
{
  Process::Journal()<<"Transport_Interfaces_FT_Disc::preparer_calcul"<<finl;

  const double temps = schema_temps().temps_courant();
  // La ligne suivante doit figurer avant le premier remaillage
  // car le remaillage utilise les angles de contact (lissage courbure)
  le_dom_Cl_dis.valeur().initialiser(temps);

  if (probleme().reprise_effectuee())
    {
      Cerr << "Calculation restarting : no initial meshing" << finl;
    }
  else
    {
      Cerr << "Initial condition remeshing." << finl;
      //maillage_interface().nettoyer_maillage(); // GB : Pour supprimer toutes ces facettes nulles...
      remailler_interface();
      //maillage_interface().nettoyer_maillage(); // GB : Pour supprimer toutes ces facettes nulles...
    }

  // Ajout pour la sauvegarde au premier pas de temps si reprise
  indicatrice_.changer_temps(temps);
  variables_internes_->indicatrice_cache.changer_temps(temps);
  indicatrice_faces_.changer_temps(temps);
  variables_internes_->vitesse_filtree.changer_temps(temps);
  variables_internes_->tmp_flux.changer_temps(temps);
  variables_internes_->index_element.changer_temps(temps);
  variables_internes_->nelem_par_direction.changer_temps(temps);
  variables_internes_->distance_interface.changer_temps(temps);
  variables_internes_->distance_interface_faces.changer_temps(temps);
  variables_internes_->distance_interface_faces_corrigee.changer_temps(temps);
  variables_internes_->distance_interface_faces_difference.changer_temps(temps);
  vitesse_imp_interp_.changer_temps(temps);
  variables_internes_->normale_interface.changer_temps(temps);
  variables_internes_->surface_interface.changer_temps(temps);

  //calcul de l'indicatrice
  indicatrice_.valeurs() = get_update_indicatrice().valeurs();
  get_update_distance_interface();
  get_update_normale_interface();

  // On verifie que la methode de transport a bien ete fournie dans le jeu
  // de donnees:
  if (variables_internes_->methode_transport == Transport_Interfaces_FT_Disc_interne::INDEFINI)
    {
      Cerr << "Error for the method Transport_Interfaces_FT_Disc::preparer_calcul\n"
           << "The transport method has not indicated for the equation\n"
           << le_nom() << finl;
      assert(0);
      Process::exit();
    }

  // Pour le premier postraitement :
  //variables_internes_->maillage_pour_post.recopie(maillage,
  //                                                  Maillage_FT_Disc::MINIMAL);

  //Ajout TF : gestion de la derivee en temps
  if (calculate_time_derivative()) derivee_en_temps().changer_temps(temps);
  //Fin TF
  return 1;
}

void Transport_Interfaces_FT_Disc::preparer_pas_de_temps(void)
{

}

double Transport_Interfaces_FT_Disc::calculer_pas_de_temps(void) const
{
  // TODO:
// We should think of implementing it as in eq. 16 for instance:
  // https://hal.archives-ouvertes.fr/hal-02304125/document
  // and use the Lagrangian min edge length instead of Delta_x.
  //
  // Ou encore :
  // Par ailleurs du fait de la contrainte sur le pas de temps capillaire (Popinet 2009) cette approche
  // connait egalement des limites (Pierson 2021). En effet, quand les films deviennent suffisamment minces,
  // et que les cellules composant le film sont suffisamment petites, le pas de temps necessaire au calcul
  // decroit drastiquement, comme la taille des cellules a la puissance trois-demi.
  // Pierson, J. (2021). Numerical study of drop bouncing on a fluid-fluid interface. ICTAM, Milan.
  // Popinet, S. (2020). A vertically-Lagrangian, non-hydrostatic, multilayer model for multiscale free-surface flows. Journal of Computational Physics, 109609
  //
  return DMAXFLOAT;
}

/*! @brief Recalcul du champ variables_internes_->indicatrice_cache a partir de la position des interfaces.
 *
 *   ATTENTION, ce n'est pas l'inconnue du probleme. L'inconnue est mise a jour
 *   a partir de ce champ dans mettre_a_jour.
 *
 */
const Champ_base& Transport_Interfaces_FT_Disc::get_update_indicatrice()
{
  const int tag = maillage_interface().get_mesh_tag();
  if (tag != variables_internes_->indicatrice_cache_tag)
    {
      DoubleVect& valeurs_indicatrice = variables_internes_->indicatrice_cache.valeur().valeurs();
      maillage_interface().parcourir_maillage();
      maillage_interface().calcul_indicatrice(valeurs_indicatrice,
                                              valeurs_indicatrice);
      variables_internes_->indicatrice_cache_tag = tag;
    }
  return variables_internes_->indicatrice_cache.valeur();
}

/*! @brief Interpolation lineaire d'un champ de vitesse VDF aux faces en un point de coordonnees coord_som.
 *
 * Le point coord_som est suppose se trouver dans
 *   l'element "element".
 *   Pour chaque composante, on cherche le cube contenant coord_som et dont
 *   les sommets sont des noeuds de vitesse pour cette composante et on
 *   interpole lineairement dans ce cube. Au bord du domaine, la vitesse
 *   tangentielle dans le demi-element colle a la paroi est la vitesse
 *   discrete de l'element (la condition aux limites n'est PAS utilisee).
 *
 */
void interpoler_vitesse_point_vdf(const Champ_base& champ_vitesse,
                                  const FTd_vecteur3& coord_som,
                                  const int element,
                                  FTd_vecteur3& vitesse)
{
  const DoubleTab& valeurs_v = champ_vitesse.valeurs();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, champ_vitesse.domaine_dis_base());
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const DoubleTab& xv = domaine_vf.xv();
  const DoubleTab& xp = domaine_vf.xp();
  const IntTab& face_voisins = domaine_vf.face_voisins();

  // La vitesse est la moyenne ponderee des 4 ou 8 vitesses
  // aux faces qui encadrent le sommet.
  // On designe par i,j,k l'une des 8 faces, avec 0<=i<=1, 0<=j<=1, 0<=k<=1.
  // La composante vitesse au sommet est
  //  vitesse = SOMME(vitesse(face_i,j,k)*coef[0]*coef[1]*coef[2]);
  // Exemple : interpolation de la composante horizontale de vitesse
  //  pour le sommet "x". On fait une interpolation bilineaire
  //  entre les quatre composantes horizontales de vitesses v (8 vitesses en 3d)
  //  On a besoin de l'element voisin situe au dessus car le sommet
  //  est dans la moitie superieure de l'element :
  //   --------
  //  |        |
  //  | v10    | v11
  //  |-->     |-->
  //  |        |
  //  |        |
  //   --------
  //  |     x  |
  //  | v00    | v01
  //  |-->     |-->
  //  |        |
  //  |        |
  //   --------

  // Boucle sur les 2 ou 3 composantes de la vitesse
  int compo;
  const int dim = Objet_U::dimension;
  for (compo = 0; compo < dim; compo++)
    {
      // ****
      // Calcul des 2 ou 3 coefficients de ponderation coef[x]
      // ****
      double coef[3];
      // Numeros des elements voisins dans chaque direction
      // (le numero de l'element voisins dans la direction "compo" n'est pas
      //  utilise et reste a -1)
      int elem_voisins[3] = { -1, -1, -1 };
      // Indice local de la face de "element" qui est commune avec elem_voisin[i]:
      int indice_local_face_voisine[3] = { -1, -1, -1 };
      // Numero du quatrieme element utilise pour interpoler la vitesse
      // (les elements utilises en 3D sont "element", "elem_voisins[(compo+1)%3]",
      //  "elem_voisins[(compo+2)%3]" et "elem_diagonal".
      int elem_diagonal = -1;

      // Boucle sur les 2 ou 3 dimensions du cube dans lequel on interpole la valeur
      int direction;
      for (direction = 0; direction < dim; direction++)
        {
          // Coordonnee du point ou il faut interpoler la vitesse
          const double x = coord_som[direction];
          double a; // Le coef d'interpolation dans la direction
          if (direction == compo)
            {
              // Dans la direction de la composante traitee on interpole
              // entre les deux faces opposees de l'element. Coordonnees
              // de ces faces :
              const int face_inf_compo = elem_faces(element, compo);
              const int face_sup_compo = elem_faces(element, compo + dim);
              const double xmin = xv(face_inf_compo, direction);
              const double xmax = xv(face_sup_compo, direction);
              a = (xmax - x) / (xmax - xmin);
            }
          else
            {
              // Dans les autres directions, on cherche les faces des
              // elements voisins si elles existent.
              // Calcul de l'element voisin a utiliser: le plus proche du point
              // a interpoler
              const double centre_elem = xp(element, direction);
              // Indice local dans l'element de la face la plus proche dans la
              // direction i
              const int i_face_voisine = direction + ((x < centre_elem) ? 0 : dim);
              // Indice de la face dans le domaine
              const int face_voisine = elem_faces(element, i_face_voisine);
              // Indice de l'element voisin par cette face
              const int elem_voisin =
                face_voisins(face_voisine, 0) + face_voisins(face_voisine, 1) - element;
              if (elem_voisin >= 0)
                {
                  // Il y a un voisin: coordonnee du centre de la face voisine
                  const double centre_voisin = xp(elem_voisin, direction);
                  a = (centre_voisin - x) / (centre_voisin - centre_elem);
                  elem_voisins[direction] = elem_voisin;
                  indice_local_face_voisine[direction] = i_face_voisine;
                }
              else
                {
                  // Pas de voisin (face de bord). La vitesse est supposee constante
                  // dans la direction "i".
                  a = 1.;
                }
            }
          if (a < 0.)
            a = 0.;
          if (a > 1.)
            a = 1.;
          coef[direction] = a;
        }

      if (dim == 3)
        {
          const int direction1 = (compo + 1) % 3;
          const int direction2 = (compo + 2) % 3;
          const int elem_voisin1 = elem_voisins[direction1];
          const int elem_voisin2 = elem_voisins[direction2];
          // Recherche de l'element diagonal:
          if (elem_voisin1 >= 0 && elem_voisin2 >= 0)
            {
              // On cherche l'element voisin de elem_voisin1 dans la direction direction2:
              int i_face_voisine, face_voisine, elem_diagonal1, elem_diagonal2;

              i_face_voisine = indice_local_face_voisine[direction2];
              face_voisine = elem_faces(elem_voisin1, i_face_voisine);
              // element diagonal obtenu par le voisin de elem_voisin1:
              elem_diagonal1 = face_voisins(face_voisine, 0)
                               + face_voisins(face_voisine, 1) - elem_voisin1;
              if (elem_diagonal1 >= 0)
                {
                  // Ok, on a un element diagonal, on verifie qu'on l'a aussi par
                  // l'autre cote (voisin de elem_voisin2 dans la direction direction1) :
                  i_face_voisine = indice_local_face_voisine[direction1];
                  face_voisine = elem_faces(elem_voisin2, i_face_voisine);
                  // element diagonal obtenu par le voisin de elem_voisin2:
                  elem_diagonal2 = face_voisins(face_voisine, 0)
                                   + face_voisins(face_voisine, 1) - elem_voisin2;
                  if (elem_diagonal1 == elem_diagonal2)
                    {
                      // Les deux elements existent et sont identiques:
                      elem_diagonal = elem_diagonal1;
                    }
                  else
                    {
                      // L'un des deux elements n'existe pas (il y a un bord)
                      // => pas d'element diagonal
                      elem_diagonal = -1;
                    }
                }
              else
                {
                  // Pas de voisin => pas d'element diagonal
                  elem_diagonal = -1;
                }
              if (elem_diagonal < 0)
                {
                  // Pas d'element diagonal => il y a un coin => pas d'interpolation.
                  elem_voisins[direction1] = -1;
                  elem_voisins[direction2] = -1;
                }
            }
          else
            {
              // Il n'y a qu'un voisin au maximum, donc pas d'element diagonal.
              // Rien a faire.
            }
        }
      if (dim == 2)
        {
          // Numeros des faces utilisees pour interpoler le champ:
          const int direction1 = 1 - compo;
          int element_voisin = elem_voisins[direction1];
          if (element_voisin < 0)
            element_voisin = element;
          const int f00 = elem_faces(element, compo);
          const int f10 = elem_faces(element, compo + dim);
          const int f01 = elem_faces(element_voisin, compo);
          const int f11 = elem_faces(element_voisin, compo + dim);
          // Coefficient d'interpolation dans la direction de la composante traitee
          double ci = coef[compo];
          // Coefficient dans l'autre direction
          double cj = coef[direction1];
          vitesse[compo] =
            ci      * cj      * valeurs_v(f00)
            + (1.-ci) * cj      * valeurs_v(f10)
            + ci      * (1.-cj) * valeurs_v(f01)
            + (1.-ci) * (1.-cj) * valeurs_v(f11);
          if (Objet_U::bidim_axi && (compo==0) && (ci>Objet_U::precision_geom) && (ci<1.-Objet_U::precision_geom)
              && (xv(f00,0) <DMINFLOAT)
              && ((fabs(valeurs_v(f00)-valeurs_v(f10))>DMINFLOAT) || (fabs(valeurs_v(f01)-valeurs_v(f11))>DMINFLOAT)))
            {
              Cerr << "In bidim_axi, when interpolating u_r within the first cell, we use the value on the symetry axis u_r(r=0)=0." << finl;
              Cerr << "We take a simple mean on that and the value at the other face. But for a divergence-free field, neglecting dv/dy, " << finl;
              Cerr << "it would be better to assume a velocity as u(x) = x_1/x * u_1 (if x!=0). GB 2020/03/05." << finl;
              const double x = coord_som[0];
              Cerr << "Here, the difference is "<< (1.-ci) << " vs. " << xv(f10, 0)/x << finl;
              Cerr << "u1= " << valeurs_v(f10) << " direction1 : " << direction1 << finl;
              Cerr << "interpoler_vitesse_point_vdf of Transport_Interface..cpp not exiting but interpolation is adapted" << finl;
              Cerr << "Former: " << vitesse[compo];
              vitesse[compo] = xv(f10, 0)/x * (
                                 cj      * valeurs_v(f10)
                                 + (1.-cj) * valeurs_v(f11));
              Cerr << " New: " << vitesse[compo] << finl;
              //Process::exit();
            }
        }
      else if (dim == 3)
        {
          // 4 elements utilises pour interpoler la valeur:
          const int direction1 = (compo + 1) % 3;
          const int direction2 = (compo + 2) % 3;
          int elem00 = element;
          int elem10 = elem_voisins[direction1];
          int elem01 = elem_voisins[direction2];
          int elem11 = elem_diagonal;
          if (elem10 < 0) elem10 = element;
          if (elem01 < 0) elem01 = element;
          if (elem11 < 0) elem11 = element;
          // Numeros des faces utilisees pour interpoler le champ:
          const int f000 = elem_faces(elem00, compo);
          const int f100 = elem_faces(elem00, compo + dim);
          const int f010 = elem_faces(elem10, compo);
          const int f110 = elem_faces(elem10, compo + dim);
          const int f001 = elem_faces(elem01, compo);
          const int f101 = elem_faces(elem01, compo + dim);
          const int f011 = elem_faces(elem11, compo);
          const int f111 = elem_faces(elem11, compo + dim);
          // Coefficients d'interpolation
          double ci = coef[compo];
          double cj = coef[direction1];
          double ck = coef[direction2];
          // Calcul de la valeur interpolee
          vitesse[compo] =
            ( ci      * cj      * valeurs_v(f000)
              + (1.-ci) * cj      * valeurs_v(f100)
              + ci      * (1.-cj) * valeurs_v(f010)
              + (1.-ci) * (1.-cj) * valeurs_v(f110)) * ck
            +(  ci      * cj      * valeurs_v(f001)
                + (1.-ci) * cj      * valeurs_v(f101)
                + ci      * (1.-cj) * valeurs_v(f011)
                + (1.-ci) * (1.-cj) * valeurs_v(f111)) * (1. - ck);

        }
      else
        {
          // Hey, what are you doing here ?
          Process::exit();
        }
    }
}


// Description:
//  Interpolation uni-lineaire d'un champ de vitesse VDF aux faces en un point
//  de coordonnees coord_som. Le point coord_som est suppose se trouver dans
//  l'element "element". Contrairement a la methode interpoler_vitesse_point_vdf
//  dont elle s'inspire, cette methode n'interpole chaque composante que dans sa propre
//  direction
//  Pour chaque composante, on cherche le cube contenant coord_som et dont
//  les sommets sont des noeuds de vitesse pour cette composante et on
//  interpole lineairement dans ce cube. Au bord du domaine, la vitesse
//  tangentielle dans le demi-element colle a la paroi est la vitesse
//  discrete de l'element (la condition aux limites n'est PAS utilisee).
void interpoler_simple_vitesse_point_vdf(const Champ_base& champ_vitesse,
                                         const FTd_vecteur3& coord_som,
                                         const int element,
                                         FTd_vecteur3& vitesse)
{
  const DoubleTab& valeurs_v = champ_vitesse.valeurs();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, champ_vitesse.domaine_dis_base());
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const DoubleTab& xv = domaine_vf.xv();

  // Chaque composante de vitesse est la moyenne ponderee des 2 vitesses
  // aux faces de l'element qui encadrent le sommet.
  // On designe par i,j,k l'une des 8 faces, avec 0<=i<=1, 0<=j<=1, 0<=k<=1.
  // La composante vitesse au sommet est
  //  vitesse = SOMME(vitesse(face_i,j,k)*coef[0]*coef[1]*coef[2]);
  // Exemple : interpolation de la composante horizontale de vitesse
  //  pour le sommet "x". On fait une interpolation bilineaire
  //  entre les DEUX composantes horizontales de vitesses v.
  //  Pas besoin des voisins!
  //   --------
  //  |     x  |
  //  |  v0    | v1
  //  |-->     |-->
  //  |        |
  //  |        |
  //   --------

  // Boucle sur les 2 ou 3 composantes de la vitesse
  int compo;
  const int dim = Objet_U::dimension;
  for (compo = 0; compo < dim; compo++)
    {
      // ****
      // Calcul d'un unique coefficient de ponderation coef
      // ****
      double coef;

      // Boucle sur les 2 ou 3 dimensions du cube dans lequel on interpole la valeur
      int direction = compo;
      // Coordonnee du point ou il faut interpoler la vitesse
      const double x = coord_som[direction];
      // Dans la direction de la composante traitee on interpole
      // entre les deux faces opposees de l'element. Coordonnees
      // de ces faces :
      const int face_inf_compo = elem_faces(element, compo);
      const int face_sup_compo = elem_faces(element, compo + dim);
      const double xmin = xv(face_inf_compo, direction);
      const double xmax = xv(face_sup_compo, direction);
      coef = (xmax - x) / (xmax - xmin); // Le coef d'interpolation dans la direction

      if (coef < 0.)
        coef = 0.;
      if (coef > 1.)
        coef = 1.;

      vitesse[compo] = coef * valeurs_v(face_inf_compo) + (1.-coef) * valeurs_v(face_sup_compo);
    }
}

double Transport_Interfaces_FT_Disc::calculer_integrale_indicatrice(const DoubleVect& indicatrice, double& integrale_ph0) const
{
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
  const DoubleVect& volumes = domaine_vf.volumes();
  int nd, nb_nd = indicatrice.size();
  assert(nb_nd==domaine_vf.nb_elem());

  double integrale = 0.;
  integrale_ph0 = 0.;
  for (nd=0 ; nd<nb_nd ; nd++)
    {
      integrale += indicatrice(nd) * volumes(nd);
      integrale_ph0 += (1 - indicatrice(nd))*volumes(nd);
    }
  integrale_ph0 = mp_sum(integrale_ph0);
  integrale = mp_sum(integrale);
  return integrale;
}

/*! @brief Calcul de la vitesse de deplacement des noeuds du maillage a partir d'un champ eulerien par interpolation.
 *
 *   Le deplacement fourni n'a aucune propriete particuliere de conservation
 *   du volume.
 *   Les lignes de contact sont deplacees avec une vitesse qui n'a pas de
 *   propriete particuliere non plus...
 *
 *  ATTENTION : on evalue simplement la vitesse a l'endroit ou sont les sommets.
 *  Param nv_calc : si =1 : recalcule le champ eulerien de la vitesse par filtrage L2
 *    sinon, reutilise celui stocke dans variables_internes_
 *
 */
void Transport_Interfaces_FT_Disc::calculer_vitesse_transport_interpolee(
  const Champ_base&        champ_vitesse,
  const Maillage_FT_Disc& maillage,
  DoubleTab&             vitesse_noeuds,
  int                   nv_calc,
  int standard) const
{

  switch(variables_internes_->methode_interpolation_v)
    {
    case Transport_Interfaces_FT_Disc_interne::VALEUR_A_ELEM:
      {
        const Champ_base * champ_vitesse_interp = 0; // Le champ a utiliser pour interpoler
        if (sub_type(Champ_P1NC, champ_vitesse) || sub_type(Champ_Fonc_P1NC,champ_vitesse))
          {
            // Champ P1NC en VEF: on projette L2 sur un champ P1,
            // puis on interpole lineairement
            Champ_base& champ_filtre = variables_internes_->vitesse_filtree.valeur();
            champ_vitesse_interp = & champ_filtre;
            if (nv_calc)
              {
                // Premier jet :
                // Calcul d'un champ aux sommets par filtrage L2 (inversion d'une matrice)
                // c'est assez long...
                champ_filtre.valeurs() = champ_vitesse.valeurs();
                ////if (type == "Champ_P1NC")
                if (sub_type(Champ_P1NC, champ_vitesse))
                  ref_cast(Champ_P1NC, champ_vitesse).filtrer_L2(champ_filtre.valeurs());
                else
                  ref_cast(Champ_Fonc_P1NC, champ_vitesse).filtrer_L2(champ_filtre.valeurs());
              }
          }
        else if (sub_type(Champ_Face_VDF, champ_vitesse) || sub_type(Champ_Fonc_Face_VDF,champ_vitesse))
          {
            // On appelle 'valeur_aux_elems' du champ aux faces.
            champ_vitesse_interp = &champ_vitesse;
          }
        else
          {
            Cerr << "Error for the method Transport_Interfaces_FT_Disc::calculer_vitesse_transport_interpolee\n";
            Cerr << "The interpolation from a field " << champ_vitesse.que_suis_je();
            Cerr << "is not developped." << finl;
            assert(0); // a coder...
            exit();
          }

        // Calcul de la vitesse interpolee a partir du champ P1
        // On ne calcule que les vitesses des sommets reels :
        DoubleTab& les_positions = variables_internes_->doubletab_pos;
        IntVect& les_elements = variables_internes_->intvect_elements;
        DoubleTab& les_vitesses = variables_internes_->doubletab_vitesses;
        // Remplissage des tableaux :
        const DoubleTab& pos = maillage.sommets();
        const ArrOfInt& elem = maillage.sommet_elem();
        const int nb_pos_tot = pos.dimension(0);
        const int dim = Objet_U::dimension;
        les_positions.resize(nb_pos_tot, dim);
        les_elements.resize(nb_pos_tot);
        int i, j;
        int nb_positions = 0;
        for (i = 0; i < nb_pos_tot; i++)
          {
            const int num_elem = elem[i];
            if (num_elem >= 0)
              {
                for (j = 0; j < dim; j++)
                  les_positions(nb_positions, j) = pos(i, j);
                les_elements(nb_positions) = num_elem;
                nb_positions++;
              }
          }
        les_positions.resize(nb_positions, dim);
        les_elements.resize(nb_positions);
        les_vitesses.resize(nb_positions, dim);

        champ_vitesse_interp->valeur_aux_elems(les_positions, les_elements, les_vitesses);

        // Copie des vitesses :
        vitesse_noeuds.resize(nb_pos_tot, dim);
        nb_positions = 0;
        for (i = 0; i < nb_pos_tot; i++)
          {
            if (elem[i] >= 0)
              {
                for (j = 0; j < dim; j++)
                  vitesse_noeuds(i, j) = les_vitesses(nb_positions, j);
                nb_positions++;
              }
          }
        maillage.desc_sommets().echange_espace_virtuel(vitesse_noeuds);

        break;
      }
    case Transport_Interfaces_FT_Disc_interne::VDF_LINEAIRE:
      {
        if (!sub_type(Champ_Face_VDF, champ_vitesse) && !sub_type(Champ_Fonc_Face_VDF,champ_vitesse) )
          {
            Cerr << "Error for the method Transport_Interfaces_FT_Disc::calculer_vitesse_transport_interpolee\n"
                 << "the interpolation VDF_LINEAIRE is valid only for a VDF discretization with a Champ_face field\n"
                 << " (type for the current field: " << champ_vitesse.que_suis_je() << finl;
            Process::exit();
          }
        const DoubleTab& pos = maillage.sommets();
        const ArrOfInt& elem = maillage.sommet_elem();
        const int nb_pos_tot = pos.dimension(0);
        const int dim        = Objet_U::dimension;
        int i;
        FTd_vecteur3 coord;
        FTd_vecteur3 vitesse;
        vitesse_noeuds.resize(nb_pos_tot, dim);
        for (i = 0; i < nb_pos_tot; i++)
          {
            const int element = elem[i];
            if (element >= 0)   // sommet reel ?
              {
                int j;
                for (j = 0; j < dim; j++)
                  coord[j] = pos(i,j);
                if (standard)
                  {
                    // Interpolation M-lineaire dans toutes les M-directions pour chaque compo
                    interpoler_vitesse_point_vdf(champ_vitesse, coord, element, vitesse);
                  }
                else
                  {
                    // Interpolation uni-lineaire dans chaque direction de chaque compo :
                    interpoler_simple_vitesse_point_vdf(champ_vitesse, coord, element, vitesse);
                  }
                for (j = 0; j < dim; j++)
                  vitesse_noeuds(i,j) = vitesse[j];
              }
          }
        maillage.desc_sommets().echange_espace_virtuel(vitesse_noeuds);
        break;
      }
    default:
      {
        Cerr << "Transport_Interfaces_FT_Disc::calculer_vitesse_transport_interpolee\n"
             << " interpolation method not developped" << finl;
        Process::exit();
      }
    }
}


void Transport_Interfaces_FT_Disc::calculer_scalaire_interpole(
  const Champ_base&        champ_scal,
  const Maillage_FT_Disc& maillage,
  DoubleTab&               val_scal_noeuds,
  int                   nv_calc) const
{

  switch(variables_internes_->methode_interpolation_v)
    {
    case Transport_Interfaces_FT_Disc_interne::VALEUR_A_ELEM:
      {

        Champ champ_scal_interp;
        champ_scal_interp=champ_scal;
        if (sub_type(Champ_P1NC, champ_scal) || sub_type(Champ_Fonc_P1NC,champ_scal))
          {

            if (nv_calc)
              {
                // Premier jet :
                // Calcul d'un champ aux sommets par filtrage L2 (inversion d'une matrice)
                // c'est assez long...

                DoubleTab scal_filtre_val(champ_scal.valeurs());

                if (sub_type(Champ_P1NC, champ_scal))
                  ref_cast(Champ_P1NC, champ_scal).filtrer_L2(scal_filtre_val);
                else if (sub_type(Champ_Fonc_P1NC, champ_scal))
                  ref_cast(Champ_Fonc_P1NC, champ_scal).filtrer_L2(scal_filtre_val);
                champ_scal_interp.valeurs() = scal_filtre_val;
              }

          }
        else if (sub_type(Champ_Fonc_P0_base, champ_scal))
          champ_scal_interp.valeurs() = champ_scal.valeurs();

        else
          {
            Cerr<<"champ_scal type ="<<champ_scal.que_suis_je()<<finl;
            Cerr << "Erreur dans Transport_Interfaces_FT_Disc::calculer_scalaire_interpole\n";
            Cerr << " L'interpolation a partir d'un champ " << champ_scal.que_suis_je();
            Cerr << " n'est pas codee." << finl;

            Cerr<<"champ_scal type ="<<champ_scal.que_suis_je()<<finl;
            Cerr << "Error for the method Transport_Interfaces_FT_Disc::calculer_scalaire_interpole\n";
            Cerr << " The interpolation from a field " << champ_scal.que_suis_je();
            Cerr << " is not developped." << finl;
            assert(0); // a coder...
            Process::exit();
          }

        // Calcul de la vitesse interpolee a partir du champ P1
        // On ne calcule que les vitesses des sommets reels :
        DoubleTab& les_positions = variables_internes_->doubletab_pos;
        IntVect& les_elements = variables_internes_->intvect_elements;
        ///DoubleTab & les_vitesses = variables_internes_->doubletab_vitesses;
        // Remplissage des tableaux :
        const DoubleTab& pos = maillage.sommets();
        const ArrOfInt& elem = maillage.sommet_elem();
        const int nb_pos_tot = pos.dimension(0);
        const int dim = Objet_U::dimension;
        les_positions.resize(nb_pos_tot, dim);
        les_elements.resize(nb_pos_tot);
        int i, j;
        int nb_positions = 0;
        for (i = 0; i < nb_pos_tot; i++)
          {
            const int num_elem = elem[i];
            if (num_elem >= 0)
              {
                for (j = 0; j < dim; j++)
                  les_positions(nb_positions, j) = pos(i, j);
                les_elements(nb_positions) = num_elem;
                nb_positions++;
              }
          }
        les_positions.resize(nb_positions, dim);
        les_elements.resize(nb_positions);

        champ_scal_interp->valeur_aux_elems(les_positions, les_elements,val_scal_noeuds);

        break;
      }
    case Transport_Interfaces_FT_Disc_interne::VDF_LINEAIRE:
      {

        exit();
        break;
      }
    default:
      {
        Cerr << "Transport_Interfaces_FT_Disc::calculer_scalaire_interpole\n"
             << " the interpolation method is not coded" << finl;
        exit();
      }
    }
}


/*! @brief Calcul de la derivee en temps de l'inconnue : zero.
 *
 */
DoubleTab& Transport_Interfaces_FT_Disc::derivee_en_temps_inco(DoubleTab& derivee)
{
  derivee = 0.;
  return derivee;
}

void Transport_Interfaces_FT_Disc::assembler( Matrice_Morse& mat_morse, const DoubleTab& present, DoubleTab& secmem)
{
  Cerr<<" On ne doit pas resoudre  Transport_Interfaces_FT_Disc en implicite "<<finl;
  abort();
  return;
}

static void init_parser_v_impose(const Noms& expression_vitesse, Parser& parser_x, Parser& parser_y, Parser& parser_z, double temps)
{
  const int dimension3 = (Objet_U::dimension==3);
  // Preparation des parsers...
  std::string sx(expression_vitesse[0]);
  parser_x.setString(sx);
  parser_x.setNbVar(4);
  parser_x.addVar("x");
  parser_x.addVar("y");
  parser_x.addVar("z");
  parser_x.addVar("t");
  parser_x.parseString();
  parser_x.setVar("z", 0.);
  parser_x.setVar("t", temps);

  std::string sy(expression_vitesse[1]);
  parser_y.setString(sy);
  parser_y.setNbVar(4);
  parser_y.addVar("x");
  parser_y.addVar("y");
  parser_y.addVar("z");
  parser_y.addVar("t");
  parser_y.parseString();
  parser_y.setVar("z", 0.);
  parser_y.setVar("t", temps);

  Nom unused_expr("0");
  std::string sz(dimension3
                 ? expression_vitesse[2]
                 : unused_expr /* inutilise */);

  parser_z.setString(sz);
  parser_z.setNbVar(4);
  parser_z.addVar("x");
  parser_z.addVar("y");
  parser_z.addVar("z");
  parser_z.addVar("t");
  parser_z.parseString();
  parser_z.setVar("z", 0.);
  parser_z.setVar("t", temps);
}


//La presence d un interface solide dans l ecoulement implique d appliquer un terme source
//dans l equation de quantite de mouvement pour imposer au fluide la vitesse de l interface.
//Ce terme source n est pas standard mais consiste en une modification de vpoint qui depend
//de la valeur de vpoint elle meme.
//Contrainte :  methode_transport==VITESSE_IMPOSEE
//
//Les etapes de la methode sont :
//-Calcul de la vitesse imposee a l interface a un temps donne: calcul_vitesse(...)
//-On determine le type du domaine discretisee (is_VDF vaut 1 pour VDF et 0 pour VEF)
// et le type  dequation traitee (is_QC vaut 0 pour  Navier_Stokes_FT_Disc et 1 sinon)
//-Estimation de l increment de quantite de mouvement a ajouter a vpoint
// puis ajout de cet increment : calcul_source_et_modifie_vpoint(...)
// Rq : la methode est parallele pour vpoint
//  Le terme vpoint*dt permet d'impliciter un peu le terme source, mais
//  semble poser quelques problemes avec la gravite. Un solide fixe dans un
//  liquide sous gravite avec l'option terme_gravite rho_g met en mouvement
//  le fluide dans le sens inverse a la gravite !!!
void Transport_Interfaces_FT_Disc::modifier_vpoint_pour_imposer_vit(const DoubleTab& inco_val,DoubleTab& vpoint0,
                                                                    DoubleTab& vpoint,const DoubleTab& rho_faces,
                                                                    DoubleTab& source_val,const double temps,
                                                                    const double dt,const int is_explicite,
                                                                    const double eta)
{
  if (!(temps>temps_debut_))
    return;
  if (variables_internes_->methode_transport == Transport_Interfaces_FT_Disc_interne::VITESSE_IMPOSEE
      || variables_internes_->methode_transport == Transport_Interfaces_FT_Disc_interne::LOI_HORAIRE)
    {
      // Etape 1 : calcul de la vitesse imposee sur l'interface a un temps donne
      DoubleTab vit_imposee;
      calcul_vitesse(vit_imposee,inco_val,vpoint0,temps,dt); // vpoint0 au lieu de vpoint
      vit_imposee.echange_espace_virtuel() ;

      // Etape 2.1 : determination du domaine de discretisation
      const Domaine_dis_base& mon_dom_dis = domaine_dis().valeur();
      const IntTab& face_voisins = mon_dom_dis.face_voisins();
      assert(inco_val.dimension(0) == face_voisins.dimension(0));

      // Etape 2.2 : determiniation du systeme d'equations a resoudre
      int is_QC=0;
      const Equation_base& eq = probleme_base_->equation(0);
      if (sub_type(Navier_Stokes_FT_Disc,eq))
        is_QC=0;
      else
        is_QC=1;

      // Etape 3 : calcul et ajout du terme de penalisation dans la qdm
      const DoubleTab& indicatrice = get_update_indicatrice().valeurs();

      calcul_indicatrice_faces(indicatrice,face_voisins);
      const DoubleTab& indicatrice_faces = get_indicatrice_faces().valeurs();

      // Si on souhaite regulariser l interface IBC/fluide :
      // combinaison convexe de la vitesse imposee et la vitesse predite.
      // Ceci est realise seulement pendant une etape explicite.
      // En penalise, pas de regularisation pendant l etape implicite
      // (pour regulariser en penalise, on ajoute une etape explicite)
      const int n = vpoint.dimension(0);
      const int m = vpoint.line_size();
      if(variables_internes_->vimp_regul && is_explicite)
        {

          DoubleTab vitesse(vpoint0);
          vitesse *= dt ;
          vitesse += inco_val ;

          for (int i=0 ; i < n; i++)
            {
              if (indicatrice_faces(i) > 0. )
                {
                  double f = indicatrice_faces(i);
                  for (int j = 0; j < m; j++)
                    {
                      if (!is_QC)
                        vit_imposee(i,j) = f*vit_imposee(i,j) + (1.-f)*vitesse(i,j);
                      else
                        vit_imposee(i,j) = f*vit_imposee(i,j) + (1.-f)*vitesse(i,j)/rho_faces(i);

                      vitesse_imp_interp_.valeur().valeurs()(i,j)= vit_imposee( i,j ) ;
                    }

                }
            }
        }

      calcul_source(inco_val,vpoint0,rho_faces,source_val,vit_imposee,indicatrice_faces,
                    is_QC,dt,is_explicite,eta); // vpoint0 au lieu de vpoint
      source_val.echange_espace_virtuel();

      const DoubleVect& volumes_entrelaces = ref_cast(Domaine_VF,mon_dom_dis).volumes_entrelaces();
      const Solveur_Masse& le_solveur_masse = eq.solv_masse();
      int i, j;

      DoubleTab termes_sources_face(vpoint);
      termes_sources_face=0.;
      modifie_source(termes_sources_face,source_val,rho_faces,n,m,is_QC,volumes_entrelaces,le_solveur_masse);
      termes_sources_face.echange_espace_virtuel() ; // CI

      for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
          vpoint(i, j) += termes_sources_face(i,j);
    }
  else if (sub_type(Transport_Marqueur_FT,*this))
    {

    }
  else
    {
      Cerr << "Error for the method Transport_Interfaces_FT_Disc::modifier_vpoint_pour_imposer_vit\n"
           << " The transport equation is not of type \"methode_transport vitesse_imposee\"."
           << " or \"methode_transport loi_horaire\"."
           << finl;
      exit();
    }
  vpoint.echange_espace_virtuel();
  Debog::verifier("Transport_Interfaces_FT_Disc::calculer_vitesse_imposee vpoint",vpoint);
}


void Transport_Interfaces_FT_Disc::calcul_indicatrice_faces(const DoubleTab& indicatrice,
                                                            const IntTab& face_voisins)
{
  DoubleTab& indicatrice_faces = indicatrice_faces_.valeurs();
  const int nfaces = face_voisins.dimension_tot(0);
  for (int i = 0; i < nfaces; i++)
    {
      const int elem0 = face_voisins(i, 0);
      const int elem1 = face_voisins(i, 1);
      indicatrice_faces(i)= 0.;

      if (elem0 >= 0)
        indicatrice_faces(i) = indicatrice(elem0);
      if (elem1 >= 0)
        indicatrice_faces(i) += indicatrice(elem1);
      if (elem0 >= 0 && elem1 >= 0)
        indicatrice_faces(i) *= 0.5;
      double bmax=1.-1e-9;
      if(indicatrice_faces(i) <= 1.0e-9) indicatrice_faces(i) = 0. ;
      else if(indicatrice_faces(i) >= bmax) indicatrice_faces(i) = 1. ;
    }

  switch(variables_internes_->type_indic_faces_)
    {
    case Transport_Interfaces_FT_Disc_interne::STANDARD:
      {
        // Does nothing more.
        break;
      }
    case Transport_Interfaces_FT_Disc_interne::MODIFIEE:
      {
        // si on souhaite calculer l'indicatrice a partir de la distance :
        const DoubleTab& dist_face = get_update_distance_interface_faces().valeurs();
        const Domaine_dis_base& domaine_dis_base = domaine_dis().valeur();
        const Domaine_VDF&   domaine_vdf       = ref_cast(Domaine_VDF, domaine_dis_base);
        const DoubleVect& face_surfaces = domaine_vdf.face_surfaces();
        const DoubleVect& volumes_entrelaces = domaine_vdf.volumes_entrelaces();
        double& position  = variables_internes_->modified_indic_faces_position;
        double& thickness = variables_internes_->modified_indic_faces_thickness;
        for (int i = 0; i < nfaces; i++)
          {
            double h=volumes_entrelaces(i)/face_surfaces(i);
            if (dist_face(i) > (position+thickness/2.)*h || indicatrice_faces(i)==1.)
              indicatrice_faces(i)=1.;
            else if (dist_face(i) >= (position-thickness/2.)*h && thickness !=0.)
              indicatrice_faces(i) = (dist_face(i) - position*h)/(thickness*h) + 0.5 ;
            else
              indicatrice_faces(i)=0.;
          }
        break;
      }
    case Transport_Interfaces_FT_Disc_interne::AI_BASED:
      {
        // This method should be very much like case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED in Navier_Stokes_FT_Disc::calculer_dI_dt
        // (but just the part to compute indicatrices faces...
        // WARNING : contrary to what is done in Navier_Stokes_FT_Disc::calculer_dI_dt, we compute chi_1 (ie same as indicatrice),
        // not the opposite chi_0 = 1-chi_1
        //
        const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
        const Equation_base& eqn_hydraulique = variables_internes_->refequation_vitesse_transport.valeur();
        if (sub_type(Navier_Stokes_FT_Disc, eqn_hydraulique))
          {
            // On recupere le saut de vitesse a l'interface (changement de phase)
            const Navier_Stokes_FT_Disc& ns = ref_cast(Navier_Stokes_FT_Disc, eqn_hydraulique);
            const DoubleTab& interfacial_area = ns.get_interfacial_area();
            const DoubleTab& normale_elements = get_update_normale_interface().valeurs();

            const int dim = ns.inconnue().valeurs().line_size();
            const int vef = (dim == 2);
            if (vef)
              {
                Cerr << "Code never applied or checked in VEF. You should read the algo first and assess it!" << finl;
                Process::exit();
              }
            // On fait la moyenne des 2 valeurs calculees sur les voisins
            // ATTENTION, ici on veut la valeur de chiv (cad chi_0) a la face.
            for (int face = 0; face < nfaces; face++)
              {
                double indic_face = 0.;
                int v;
                for (v = 0; v < 2; v++)
                  {
                    const int elem = face_voisins(face, v);
                    if (elem >=0)
                      {
                        // If a neighbour is pure, we use that value at the face and stop further calculation.
                        const double indic = indicatrice[elem]; // This is the value of chi_1 (ie =1 in phase 1!)
                        //if (indic == 0. || indic == 1.)
                        if (indic <=5e-3 || indic >= 1.-5e-3)
                          {
                            indic_face = indic; // Obviously, We want chi of phase_1
                            break;
                          }
                        else
                          {
                            const double surface=domaine_vf.face_surfaces(face);
                            const double ai= interfacial_area(elem); // nx pointe vers le liquide (sortant de phase 0)
                            if (fabs(ai)>DMINFLOAT)
                              {
                                double x = 0.;
                                if (vef)
                                  {
                                    for (int j = 0; j < dim; j++)
                                      {
                                        const double nf = domaine_vf.face_normales(face , j);
                                        const double nx = normale_elements(elem, j);
                                        // produit scalaire :
                                        x +=  nf*nx;
                                      }
                                    x *= ai/surface;
                                    // Que/comment Choisir?
                                    indic_face += x;
                                    Cerr << "Never tested. To be verified. It should depend on a scalar product with the vect (xp-xv)" << finl;
                                    Process::exit();
                                  }
                                else
                                  {
                                    // En VDF, l'acces a orientation permet d'eviter le calcul du produit scalaire.
                                    const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, domaine_dis().valeur());
                                    const IntVect& orientation = zvdf.orientation();
                                    const int dir = orientation[face];
                                    const double nx = normale_elements(elem, dir);
                                    // Assumes a cube, nx larger than diag means we can use the method rather safely
                                    if (nx>0.707)
                                      {
                                        x = ai/surface*nx;
                                        // On suppose que v0 est a gauche et v1 a droite!!!
                                        if (v==0)
                                          indic_face += x; // This way, we build chi_1 because normale points towards chi_1
                                        else
                                          indic_face += 1-x;
                                      }
                                    else
                                      {
                                        // L'interface croise probablement la face d'en face et la methode ne marche plus.
                                        // We revert back to the standard method :
                                        double tmp = 0.;
                                        if (elem >= 0)
                                          tmp = indicatrice(elem);

                                        double bmax=1.-1e-9;
                                        if(tmp <= 1.0e-9) tmp = 0. ;
                                        else if(tmp >= bmax) tmp = 1. ;

                                        indic_face += tmp;
                                      }
                                  }
                              }
                            else
                              {
                                Cerr <<" WTF, c'est impossible" << finl;
                                Process::exit();
                              }
                          }
                      }
                    else
                      {
                        // The only neighbour to the face :
                        const int elem_voisin = face_voisins(face, 1-v); // The other one is accessed by 1-v
                        const double indic = indicatrice[elem_voisin]; // This is the value of chi_1 (ie =1 in phase 1!)
                        indic_face = indic; // We want chi of phase_1
                        break; // c'est important pour le if d'apres.
                      }
                  }
                if (v==2)
                  // On n'a pas touche le break, on est donc passe 2 fois. donc :
                  indic_face*=0.5;

                // assert((indic_face >=0) && (indic_face<=1.));
                // ca arrive des petits derapages..
                if (indic_face <0)
                  indic_face=0.;
                if (indic_face >1.)
                  indic_face=1.;

                indicatrice_faces_(face) = indic_face;
              }
          }
        else
          {
            Cerr << "Interpolation option AI_BASED in Transport_Interfaces_FT_Disc is not available "
                 << "in Transport_Interfaces_FT_Disc::calcul_indicatrice_faces if we do not "
                 << "have a Navier_Stokes_FT_Disc equation..." << finl;
            Process::exit();
          }
        break;
      }
    default:
      Cerr << "Transport_Interfaces_FT_Disc::calcul_indicatrice_faces\n"
           << " unknown case?" << finl;
      Process::exit();
    }

  indicatrice_faces.echange_espace_virtuel();
  Debog::verifier("Transport_Interfaces_FT_Disc::calcul_indicatrice_faces indicatrice_faces",indicatrice_faces);
}

const int& Transport_Interfaces_FT_Disc::get_vimp_regul() const
{
  return variables_internes_->vimp_regul;
}
const Champ_base& Transport_Interfaces_FT_Disc::get_indicatrice_faces()
{
  return indicatrice_faces_;
}

const Champ_base& Transport_Interfaces_FT_Disc::get_compute_indicatrice_faces()
{
  const DoubleTab& indicatrice = get_update_indicatrice().valeurs();
  const Domaine_dis_base& mon_dom_dis = domaine_dis().valeur();
  const IntTab& face_voisins = mon_dom_dis.face_voisins();
  calcul_indicatrice_faces(indicatrice,face_voisins);
  return indicatrice_faces_;
}

//Methode outil qui estime la valeur d increment de qdm a ajouter a
//la derivee temporelle de la vitesse
// L'expression calculee
// Cas Front-Tracking  (is_QC vaut 0):
//  (v_impose - (champ_vitesse + vpoint*dt)) / dt * indicatrice * rho_face;
// Cas Quasi-compressible (is_QC vaut 1) :
//  (rho_np1 v_impose - (rho_vitesse_n + vpoint*dt)) / dt * indicatrice

void Transport_Interfaces_FT_Disc::calcul_source(const DoubleTab& inco_val,
                                                 const DoubleTab& vpoint,
                                                 const DoubleTab& rho_faces,
                                                 DoubleTab& source_val,const DoubleTab& vit_imposee,const DoubleTab& indicatrice_faces,
                                                 const int is_QC,
                                                 const double dt,
                                                 const int is_explicite,
                                                 const double eta)

{
  DoubleTab terme_explicite(vpoint);
  double c;
  const int nfaces = vpoint.dimension(0);
  c = 1. / eta;
  if (indicatrice_faces.dimension(0) == nfaces )
    {
      if (is_explicite)
        {
          terme_explicite *= dt;
        }
      else
        {
          terme_explicite *= 0.;
        }

      for (int i = 0; i < nfaces; i++)
        {

          double indic = (indicatrice_faces(i) > 0. ? 1.0 : 0.0) ;
          // Si on veut regulariser le passage de IBC/fluide et que l'on est penalise,
          // seules les faces telles que indic=1 sont penalisees, les autres sont DF.
          if (variables_internes_->vimp_regul && !is_explicite ) indic = (indicatrice_faces(i) ==1. ? 1.0 : 0.0);
          double tsource;
          double rho_face = rho_faces(i);
          const int dim = inco_val.line_size(); // > 1 => VEF, =1 => VDF

          for (int j = 0; j < dim; j++)
            {
              double increment_inco = inco_val(i,j) + terme_explicite(i,j);
              if (!is_QC)
                tsource = (vit_imposee(i,j) - increment_inco) / dt * indic * rho_face * c;
              else
                tsource = (vit_imposee(i,j) * rho_face - increment_inco) / dt * indic * c;

              source_val(i,j) = tsource;
            }
        }
    }
  else
    {

      Cerr << "Erreur dans Transport_Interfaces_FT_Disc::calcul_source\n"         << finl;
      Cerr << "Dimension de indicatrice_face differente de dimension de vpoint.\n" << finl;
      exit();
    }
}

void ouvrir_fichier(SFichier& os,const Nom& type, const int flag, const Transport_Interfaces_FT_Disc& equation)
{

  // flag nul on n'ouvre pas le fichier
  if (flag==0)
    return ;
  Nom fichier=Objet_U::nom_du_cas();
  if (type=="force")
    fichier+="_Force_totale_sur_";
  else if( type=="force_totale" )
    fichier+="_Friction_totale_sur_" ;
  else if( type=="Friction" )
    fichier+="_Friction_conv_diff_sur_" ;
  else if( type=="Pressure" )
    fichier+="_Friction_Pression_sur_" ;
  else
    fichier+="_Moment_total_sur_";
  fichier+=equation.le_nom();
  fichier+=".out";
  const Schema_Temps_base& sch=equation.probleme().schema_temps();
  const int precision=sch.precision_impr();
  // On cree le fichier a la premiere impression avec l'en tete ou si le fichier n'existe pas
  struct stat f;
  if (stat(fichier,&f) || (sch.nb_impr()==1 && !equation.probleme().reprise_effectuee()))
    {
      os.ouvrir(fichier);
      SFichier& fic=os;
      Nom espace="\t\t";
      fic << (Nom)"# Printing " << (type=="moment"?"of the drag moment exerted":"of the drag exerted");
      fic << " by the fluid on the interface " << equation.le_nom();
      fic << " " << (type=="moment"?"[N.m]":"[N]") << finl;
      int nb_compo=(type=="moment" && Objet_U::dimension==2?1:Objet_U::dimension);
      fic << "# Time";

      Nom ch=espace;
      if (type=="moment")
        {
          if (Objet_U::dimension==2) ch+="Mz";
          else
            {
              ch+="Mx";
              ch+=espace+"My";
              ch+=espace+"Mz";
            }
        }
      else
        {
          if (nb_compo>1) ch+="Fx";
          if (nb_compo>=2) ch+=espace+"Fy";
          if (nb_compo>=3) ch+=espace+"Fz";
        }
      fic << ch << finl;
    }
  // Sinon on l'ouvre
  else
    {
      os.ouvrir(fichier,ios::app);
    }
  os.precision(precision);
  os.setf(ios::scientific);
}

void Transport_Interfaces_FT_Disc::modifie_source(DoubleTab& termes_sources_face,const DoubleTab& source_val,const DoubleTab& rho_faces,
                                                  const int n,const int m, const int is_QC,
                                                  const DoubleVect& vol_entrelaces,const Solveur_Masse& un_solv_masse)
{

  for (int face=0; face<n; face++)
    for (int dim=0; dim<m; dim++)
      termes_sources_face(face,dim)=vol_entrelaces(face)*source_val(face,dim);

  termes_sources_face.echange_espace_virtuel() ; // CI
  un_solv_masse.appliquer(termes_sources_face);

  if (!is_QC)
    {
      for (int i = 0; i < n; i++)
        {
          const double rho_face = rho_faces(i);

          for (int j = 0; j < m; j++)
            termes_sources_face(i,j) = termes_sources_face(i,j) / rho_face;
        }
    }
}

void Transport_Interfaces_FT_Disc::impr_effort_fluide_interface( DoubleTab& source_val, DoubleTab& pressure_part, DoubleTab& friction_part  )
{
  const DoubleTab& indicatrice_faces = get_indicatrice_faces().valeurs();
  const int n = source_val.dimension(0);
  const int nbdim1 = source_val.line_size() == 1; // VDF
  const int m = source_val.line_size();

  const Domaine_dis_base& mon_dom_dis = domaine_dis().valeur();
  const Domaine_VDF * zvdf = 0;
  if (sub_type(Domaine_VDF, domaine_dis().valeur())) zvdf = &ref_cast(Domaine_VDF, domaine_dis().valeur());

  DoubleTab termes_sources_face(source_val);
  DoubleTab termes_pressure_face(pressure_part);
  DoubleTab termes_friction_face(friction_part);

  DoubleTrav values(3,dimension);
  values=0.;

  const DoubleVect& vol_entrelaces = ref_cast(Domaine_VF,mon_dom_dis).volumes_entrelaces();
  // Construction d'un tableau des items reels non communs
  ArrOfInt sequential_items_flags;
  MD_Vector_tools::get_sequential_items_flags(source_val.get_md_vector(), sequential_items_flags);


  for (int face=0; face<n; face++)
    {
      double indic = (indicatrice_faces(face) > 0. ? 1.0 : 0.0);
      double coef = vol_entrelaces(face)*indic;
      for (int dim=0; dim<m; dim++)
        {
          termes_sources_face(face,dim)=source_val(face,dim)*coef;
          termes_pressure_face(face,dim)=pressure_part(face,dim)*coef;
          termes_friction_face(face,dim)=friction_part(face,dim)*coef;
        }
      // Calcul de dforce contribution de la force du fluide sur la face i
      // si ce n'est pas une face commune a plusieurs processeurs
      if (sequential_items_flags[face])
        {
          if (nbdim1) // VDF
            {
              int j = zvdf->orientation(face);

              values(0,j) -= termes_sources_face(face,0);
              values(1,j) -= termes_pressure_face(face,0);
              values(2,j) -= termes_friction_face(face,0);
            }
          else // VEF
            {
              for (int j = 0; j < dimension; j++)
                {
                  values(0,j) -= termes_sources_face(face,j);
                  values(1,j) -= termes_pressure_face(face,j);
                  values(2,j) -= termes_friction_face(face,j);
                }
            }
        }
    }

  // Impression des efforts exerces par le fluide sur l'interface
  {
    // Ajout des differents processeurs en //
    mp_sum_for_each_item(values);

    // Impression dans les fichiers
    if (Process::je_suis_maitre())
      {
        SFichier Force;
        ouvrir_fichier(Force,"force_totale",1,*this);
        Nom espace=" \t";
        schema_temps().imprimer_temps_courant(Force);
        Force.precision(10) ;
        for(int k=0; k<dimension; k++)
//            Force << espace << dforce(k);
          Force << espace << values(0,k);
        Force << finl;
        const int impr_mom = 1 ;

        SFichier Pressure;
        ouvrir_fichier(Pressure,"Pressure",impr_mom,*this);
        schema_temps().imprimer_temps_courant(Pressure);
        Pressure.precision(10) ;
//          for(int k=0; k<pressure.size_array(); k++)
        for(int k=0; k<dimension ; k++)
//            Pressure << espace << pressure(k);
          Pressure << espace << values(1,k);
        Pressure << finl;

        SFichier Friction;
        ouvrir_fichier(Friction,"Friction",impr_mom,*this);
        schema_temps().imprimer_temps_courant(Friction);
        Friction.precision(10) ;
//          for(int k=0; k<friction.size_array(); k++)
        for(int k=0; k<dimension; k++)
//            Friction << espace << friction(k);
          Friction << espace << values(2,k);
        Friction << finl;
      }
  }
}


// Impression des forces et moment
int Transport_Interfaces_FT_Disc::impr(Sortie& os) const
{
  // Impression dans les fichiers
  if (Process::je_suis_maitre())
    {
      SFichier Force;
      ouvrir_fichier(Force,"force",1,*this);
      Nom espace=" \t";
      schema_temps().imprimer_temps_courant(Force);
      for(int k=0; k<dimension; k++)
        Force << espace << force_[k];
      Force << finl;
      const Domaine& domaine=domaine_dis().domaine();
      const int impr_mom = domaine.moments_a_imprimer();
      if (impr_mom)
        {
          SFichier Moment;
          ouvrir_fichier(Moment,"moment",impr_mom,*this);
          schema_temps().imprimer_temps_courant(Moment);
          for(int k=0; k<moment_.size_array(); k++)
            Moment << espace << moment_[k];
          Moment << finl;
        }
    }
  return 1;
}

//Cette methode actualise le critere de stationnnarite dI/dt (en derivee partielle)
//du fait que pour ce type d equation derivee_en_temps_inco fixe dI/dt a 0.
//
void Transport_Interfaces_FT_Disc::update_critere_statio()
{
  Schema_Temps_base& sch_tps = schema_temps();
  const DoubleTab& present = inconnue().valeurs();
  const DoubleTab& passe = inconnue().passe();
  const double dt = sch_tps.pas_de_temps();
  DoubleTab tab_critere(present);

  tab_critere = present;
  tab_critere -=passe;
  tab_critere /= dt;
  sch_tps.update_critere_statio(tab_critere,*this);
}

void Transport_Interfaces_FT_Disc::calcul_effort_fluide_interface(const DoubleTab& vpoint,
                                                                  const DoubleTab& rho_faces,DoubleTab& source_val,
                                                                  const int is_explicite,const double eta)
{
  const DoubleTab& indicatrice_faces = get_indicatrice_faces().valeurs();
  const int n = vpoint.dimension(0);
  const int m = vpoint.line_size();
  double c= 1./eta;
  const Domaine_dis_base& mon_dom_dis = domaine_dis().valeur();

  int is_QC=0;
  const Equation_base& eq = probleme_base_->equation(0);
  const Solveur_Masse& le_solveur_masse = eq.solv_masse();
  if (sub_type(Navier_Stokes_FT_Disc,eq))
    is_QC=0;
  else
    is_QC=1;

  if ( !is_explicite )
    {
      for (int i = 0; i<n; i++)
        {
          double indic = (indicatrice_faces(i) > 0. ? 1.0 : 0.0);
          for (int j = 0; j < m; j++)
            {
              if (!is_QC)
                source_val(i,j) -= vpoint(i,j) * indic * rho_faces(i) * c;
              else
                source_val(i,j) -= vpoint(i,j) * indic * c;
            }
        }
    }

  DoubleTab termes_sources_face(vpoint);
  const DoubleVect& vol_entrelaces = ref_cast(Domaine_VF,mon_dom_dis).volumes_entrelaces();

  for (int face=0; face<n; face++)
    for (int dim=0; dim<m; dim++)
      {
        double indic = (indicatrice_faces(face) > 0. ? 1.0 : 0.0);
        termes_sources_face(face,dim)=vol_entrelaces(face)*vol_entrelaces(face)*source_val(face,dim)*indic;
      }

  termes_sources_face.echange_espace_virtuel() ;
  le_solveur_masse.appliquer(termes_sources_face);

  // Impression des efforts exerces par le fluide sur l'interface
  {
    ArrOfDouble dforce(dimension);
    force_=0;
    moment_=0;
    const Domaine& domaine=domaine_dis().domaine();
    const int impr_mom = domaine.moments_a_imprimer();
    const ArrOfDouble& centre_gravite = domaine.cg_moments();
    const DoubleTab& centre_faces = ref_cast(Domaine_VF,domaine_dis().valeur()).xv();
    ArrOfDouble xgr(dimension);

    const Domaine_VDF * zvdf = 0;
    if (sub_type(Domaine_VDF, domaine_dis().valeur()))
      zvdf = &ref_cast(Domaine_VDF, domaine_dis().valeur());

    // Construction d'un tableau des items sequentiels
    ArrOfInt sequential_items_flags;
    MD_Vector_tools::get_sequential_items_flags(rho_faces.get_md_vector(), sequential_items_flags);

    // Calcul de la force et du moment en fonction de la discretisation
    for (int i = 0; i < n; i++)
      {
        // Calcul de dforce contribution de la force du fluide sur la face i
        // si c'est une face sequentielle
        if (sequential_items_flags[i])
          {
            dforce=0;
            if (zvdf)
              {
                int j = zvdf->orientation(i);
                dforce[j] = -termes_sources_face(i);
                force_[j] += dforce[j];
              }
            else
              for (int j = 0; j < dimension; j++)
                {
                  dforce[j] = -termes_sources_face(i,j);
                  force_[j] += dforce[j];
                }
            // Ajout de dforce au calcul eventuel du moment
            if (impr_mom)
              {
                for (int j = 0; j < dimension; j++)
                  xgr[j] = centre_faces(i,j) - centre_gravite[j];

                if (dimension==2)
                  moment_[0] += dforce[1]*xgr[0] - dforce[0]*xgr[1];
                else
                  {
                    moment_[0] += dforce[2]*xgr[1] - dforce[1]*xgr[2];
                    moment_[1] += dforce[0]*xgr[2] - dforce[2]*xgr[0];
                    moment_[2] += dforce[1]*xgr[0] - dforce[0]*xgr[1];
                  }
              }
          }
      }

    // Ajout des differents processeurs en //

    mp_sum_for_each_item(force_);
    mp_sum_for_each_item(moment_);
    // Impression dans les fichiers
    impr(Process::Journal());
  }
}

void Transport_Interfaces_FT_Disc::get_expression_vitesse_imposee(DoubleTab& vit_ibc)
{
  const double temps=le_schema_en_temps->temps_courant();
  const int dim = Objet_U::dimension;
  const int dimension3 = (dim==3);
  Parser parser_x, parser_y, parser_z;
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
  const IntTab& face_voisins = domaine_vf.face_voisins();
  const int nfaces = face_voisins.dimension(0);
  const DoubleTab& xv = domaine_vf.xv();
  const Domaine_VDF * zvdf = 0;
  if (sub_type(Domaine_VDF, domaine_vf))
    zvdf = &ref_cast(Domaine_VDF, domaine_vf);

  if (zvdf)
    vit_ibc.resize(nfaces,1);
  else
    vit_ibc.resize(nfaces,dim);

  init_parser_v_impose(variables_internes_->expression_vitesse_imposee,
                       parser_x, parser_y, parser_z,temps);
  for (int i = 0; i < nfaces; i++)
    {
      for (int j = 0; j < dim; j++)
        {
          const double coord = xv(i,j);
          parser_x.setVar(j, coord);
          parser_y.setVar(j, coord);
          if (dimension3)
            parser_z.setVar(j, coord);
        }
      if (zvdf)
        {

          switch(zvdf->orientation(i))
            {
            case 0:
              vit_ibc(i,0) = parser_x.eval();
              break;
            case 1:
              vit_ibc(i,0) = parser_y.eval();
              break;
            case 2:
              vit_ibc(i,0) = parser_z.eval();
              break;
            }
        }
      else
        {
          vit_ibc(i,0) = parser_x.eval();
          vit_ibc(i,1) = parser_y.eval();
          if (dimension3)
            vit_ibc(i,2) = parser_z.eval();
        }
    }
}

//Methode qui effectue le calcul de la vitesse imposee a l interface a l instant temps
//- Utilise l expression de la vitesse imposee specifiee par l utilisateur : expression_vitesse_imposee
//- La vitesse calculee est contenue dans vitesse_imp

void Transport_Interfaces_FT_Disc::calcul_vitesse(DoubleTab& vitesse_imp,
                                                  const DoubleTab& vitesse,
                                                  const DoubleTab& vpoint,
                                                  const double temps,
                                                  const double dt)
{
  const int dim = Objet_U::dimension;
  const int dimension3 = (dim==3);
  const Domaine_dis_base& mon_dom_dis = domaine_dis().valeur();
  const IntTab& face_voisins = mon_dom_dis.face_voisins();
  const int nfaces = face_voisins.dimension(0);
  const DoubleTab& xv = ref_cast(Domaine_VF, mon_dom_dis).xv();
  const Domaine_VDF * zvdf = 0;
  if (sub_type(Domaine_VDF, mon_dom_dis))
    zvdf = &ref_cast(Domaine_VDF, mon_dom_dis);

  if (zvdf)
    vitesse_imp.resize(nfaces,1);
  else
    vitesse_imp.resize(nfaces,dim);

  if (variables_internes_->methode_transport == Transport_Interfaces_FT_Disc_interne::LOI_HORAIRE)
    {
      ArrOfDouble coord(dim);
      ArrOfDouble v_imp(dim);
      for (int i = 0; i < nfaces; i++)
        {
          // Calcul de la vitesse au centre de la face par la loi horaire a temps+dt
          for (int j = 0; j < dim; j++)
            coord[j] = xv(i,j);
          v_imp = variables_internes_->loi_horaire_->vitesse(temps+dt,coord);

          if (zvdf)
            vitesse_imp(i,0) = v_imp[zvdf->orientation(i)];
          else
            for (int j = 0; j < dim; j++)
              vitesse_imp(i,j) = v_imp[j];
        }
    }
  else if (variables_internes_->methode_transport == Transport_Interfaces_FT_Disc_interne::VITESSE_IMPOSEE)
    {
      switch(variables_internes_->interpolation_champ_face)
        {
        case Transport_Interfaces_FT_Disc_interne::BASE:
          {
            Parser parser_x, parser_y, parser_z;
            init_parser_v_impose(variables_internes_->expression_vitesse_imposee,
                                 parser_x, parser_y, parser_z,temps);
            for (int i = 0; i < nfaces; i++)
              {

                for (int j = 0; j < dim; j++)
                  {
                    // Vitesse evaluee au centre de la face
                    const double coord = xv(i,j);
                    parser_x.setVar(j, coord);
                    parser_y.setVar(j, coord);
                    if (dimension3)
                      parser_z.setVar(j, coord);
                  }
                if (zvdf)
                  {

                    switch(zvdf->orientation(i))
                      {
                      case 0:
                        vitesse_imp(i,0) = parser_x.eval();
                        break;
                      case 1:
                        vitesse_imp(i,0) = parser_y.eval();
                        break;
                      case 2:
                        vitesse_imp(i,0) = parser_z.eval();
                        break;
                      }

                  }
                else
                  {

                    vitesse_imp(i,0) = parser_x.eval();
                    vitesse_imp(i,1) = parser_y.eval();
                    if (dimension3)
                      vitesse_imp(i,2) = parser_z.eval();
                  }
              }
            break;
          }
        case Transport_Interfaces_FT_Disc_interne::LINEAIRE:
          {
            double dt_loc = 0.;
            if ( !variables_internes_->vf_explicite )
              {
                //On calcule un dt_loc car diff et conv sont explicites pour cette estimation de vpoint
                const Equation_base& eq = probleme_base_->equation(0);
                const double nb_op = eq.nombre_d_operateurs();
                if(nb_op > 0) dt_loc += 1. / eq.operateur(0).calculer_pas_de_temps();
                if(nb_op > 1) dt_loc += 1. / eq.operateur(1).calculer_pas_de_temps();
                dt_loc = (1./dt_loc)*(le_schema_en_temps->facteur_securite_pas());
                dt_loc = Process::mp_min(dt_loc);
                // si le dt_loc obtenu est plus grand que le dt_min du jdd :
                dt_loc = std::min(dt_loc,dt);
                Cerr << "Transport_Interfaces_FT_Disc::calculer_vitesse_imposee avec dt_loc : "<<dt_loc<<finl;
              }
            if (zvdf)
              {
                const DoubleTab& dist_face = get_update_distance_interface_faces().valeurs();
                const int phase = 0;
                vitesse_imp = vpoint;
                vitesse_imp*= dt_loc;
                vitesse_imp += vitesse;
                vitesse_imp.echange_espace_virtuel() ; // CI
                DoubleTab gradient(vitesse);
                interpoler_vitesse_face(dist_face,phase,variables_internes_->n_iterations_interpolation_ibc,vitesse_imp,gradient,temps,dt);
                Debog::verifier("Transport_Interfaces_FT_Disc::calcul_vitesse vitesse_imp apres interp. ",vitesse_imp);
              }
            else
              {
                Cerr << "Transport_Interfaces_FT_Disc::calculer_vitesse_imposee\n"
                     << " interpolation lineaire non codee en VEF" << finl;
                exit();
              }
            break;
          }
        default:
          {
            Cerr << "Transport_Interfaces_FT_Disc::calculer_vitesse_imposee\n"
                 << " methode d'interpolation non codee" << finl;
            exit();
          }
        }
    }
  else
    {
      Cerr << "Error for the method Transport_Interfaces_FT_Disc::calcul_vitesse\n"
           << " The transport equation is not of type \"methode_transport vitesse_imposee\"."
           << " or \"methode_transport loi_horaire\"."
           << finl;
      exit();
    }
}


//Calcul d'une distance signee a l'interface aux centres des faces, interpolee a partir de la distance signee aux
//centres des elements (dist_elem) et de la normale a l'interface evaluee aux centres des elements euleriens (normale_elem).
//
//Pour une face, on evalue la distance entre le centre de la face et l'interface par :
//   dist_face = d1 + d2,
//   d1 = dist_elem,
//   d2 = normale_elem scalaire (position_centre_face - centre_element)
//Ensuite, la distance entre le centre d'une face et l'interface est la moyenne de toutes
//les distances calculee a l'aide des elements adjacents a cette face.
//La distance est invalide au-dela d'une certaine epaisseur autour de l'interface
//(voir iterations de lissage dans calculer_distance_interface).
void Transport_Interfaces_FT_Disc::calculer_distance_interface_faces(
  const DoubleTab& dist_elem,
  const DoubleTab& normale_elem,
  DoubleTab&        dist_face) const
{
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "Calculer_distance_interface");
  statistiques().begin_count(stat_counter);

  static const double distance_faces_invalides = -1.e30;

  const Domaine_VF&    domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const DoubleTab& xv = domaine_vf.xv();
  const DoubleTab& xp = domaine_vf.xp();
  const int nb_faces = dist_face.dimension_tot(0);
  ArrOfInt ncontrib(nb_faces);
  ncontrib = 0;
  dist_face = 0.;

  const int dim = Objet_U::dimension;
  double centre[3] = {0., 0., 0.};
  double normale[3] = {0., 0., 0.};
  // Calcul de SOMME(d1+d2) pour tous les elements voisins de chaque face :
  int elem, i;
  const int nb_elem_tot = dist_elem.dimension_tot(0);
  const int nb_faces_elem = elem_faces.line_size();

  for (elem = 0; elem < nb_elem_tot; elem++)
    {
      const double d1 = dist_elem(elem);
      // Si la distance est invalide, on passe:
      if (d1 < distance_faces_invalides)
        continue;
      // Centre de l'element et normale a l'interface pour cet element :
      for (i = 0; i < dim; i++)
        {
          centre[i] = xp(elem, i);
          normale[i] = normale_elem(elem, i);
        }
      // Boucle sur les faces de l'element
      for (i = 0; i < nb_faces_elem; i++)
        {
          const int face = elem_faces(elem, i);
          // dist_face ne contient que des faces reelles en general.
          // si la face n'est pas dans dist_face, on ne calcule pas.
          if (face < nb_faces)
            {
              double d2 = 0.;
              int j;
              for (j = 0; j < dim; j++)
                {
                  double position_centre_face = xv(face, j);
                  d2 += (position_centre_face - centre[j]) * normale[j];
                }
              dist_face(face) += d1 + d2;
              ++(ncontrib[face]);
            }
        }
    }
  // Division par le nombre d'elements voisins
  const double valeur_invalide = distance_faces_invalides * 1.1;
#ifdef __INTEL_COMPILER
#pragma novector // Desactive vectorisation sur Intel car crash sinon
#endif
  for (i = 0; i < nb_faces; i++)
    {
      const int n = ncontrib[i];
      if (n > 0)
        dist_face(i) /= n;
      else
        dist_face(i) = valeur_invalide;
    }
  dist_face.echange_espace_virtuel();
  statistiques().end_count(stat_counter);
}


const Champ_base& Transport_Interfaces_FT_Disc::get_update_distance_interface_faces() const
{
  // Si le tag du maillage et le tag du champ sont identiques, inutile de recalculer:
  const int tag = maillage_interface().get_mesh_tag();
  if (tag == variables_internes_->distance_faces_cache_tag)
    return variables_internes_->distance_interface_faces.valeur();

  variables_internes_->distance_faces_cache_tag = tag;

  const DoubleTab& dist_elem = get_update_distance_interface().valeurs();
  const DoubleTab& normale_elem = get_update_normale_interface().valeurs();
  DoubleTab&        dist_face = variables_internes_->distance_interface_faces.valeur().valeurs();

  calculer_distance_interface_faces(dist_elem, normale_elem, dist_face);
  return variables_internes_->distance_interface_faces.valeur();
}

// Verification que la vitesse imposee est uniforme
inline void check(const DoubleTab& v_imp, int& i_face, double v, int v_est_initialise)
{
  if (v_est_initialise && v_imp(i_face)!= v)
    {
      Cerr << "=====================================================================================" << finl;
      Cerr << "You defined a non-uniform function for the keyword: methode_transport vitesse_imposee" << finl;
      Cerr << "Please add:" << finl;
      Cerr << "type_vitesse_imposee analytique" << finl;
      Process::exit();
    }
  else
    {
      v = v_imp(i_face);
      v_est_initialise = 1;
    }
}
// Modifie les valeurs de la vitesse aux faces pres de l'interface (le stencil depend du nombre d'iterations) via une
// interpolation lineaire de la vitesse.
// A partir d'une indicatrice interpolee aux faces, on calcule le gradient de chaque composante de la vitesse pour les
// faces ou :
// - indicatrice=phase
// - distance face/interface valide (depend du nb iterations de lissage dans calculer_distance_interface).
// On en deduit le gradient aux faces telles que indicatrice!=phase, via un processus iteratif qui permet d'etendre
// l'interpolation aux points interieurs a l'IBC. L'interpolation s'effectue sur les quatre faces les plus proches de
// meme orientation (vitesses non colocalisees).
// Le resultat Champ est prealablement construit a partir du champ de vitesse Un+vpoint.dt dans la methode
// calcule_vitesse
//
// Attention :
// - phase=0 impose que l'indicatrice de la partie fluide soit egale a 0
// - l'interpolation est specifique au VDF
// - codage pour l'instant pour une IBC immobile ou a vitesse uniforme

void Transport_Interfaces_FT_Disc::interpoler_vitesse_face(
  const DoubleTab& distance_interface_faces,
  const int   phase,
  const int   stencil_width,
  DoubleTab& champ,
  DoubleTab& gradient,
  const double t, const double dt)
{
  const Domaine_dis_base& mon_dom_dis = domaine_dis().valeur();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, mon_dom_dis);
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, mon_dom_dis);
  const Domaine_VDF * zvdf = 0;
  if (sub_type(Domaine_VDF, mon_dom_dis))
    zvdf = &ref_cast(Domaine_VDF, mon_dom_dis);

  const IntTab& elem_faces = domaine_vf.elem_faces();
  const IntTab& faces_elem = domaine_vf.face_voisins();
  const int nfaces = faces_elem.dimension(0);
  const int nb_elem = domaine_vf.nb_elem() ;
  const int nb_elem_tot = domaine_vf.nb_elem_tot() ;
  const IntVect& orientation = domaine_vdf.orientation();
  const DoubleTab& xv = domaine_vf.xv();
  const DoubleTab& xp = domaine_vf.xp();
  assert(champ.dimension(0) == nfaces);
  assert(gradient.dimension(0) == nfaces);
  assert(distance_interface_faces.dimension(0) == nfaces);

  const int dim = Objet_U::dimension;
  const int dimension3 = (Objet_U::dimension==3);

  Parser parser_x, parser_y, parser_z;
  init_parser_v_impose(variables_internes_->expression_vitesse_imposee,parser_x, parser_y, parser_z, t);

  const double invalid_test = -1.e30;
  const double invalid_value = -2.e30;
  gradient = invalid_value;
  const double indic_phase = (phase == 0) ? 0. : 1.;
  int i_face;
  Maillage_FT_Disc& maillage = maillage_interface() ;
  const DoubleTab& indicatrice = get_update_indicatrice().valeurs();
  const DoubleTab& indicatrice_face = get_compute_indicatrice_faces().valeurs();
  // distance signee aux faces corrigee
  DoubleTab dist_face_cor(distance_interface_faces) ;
  // type de face : fluide (0), solide (1) et 0.5 si le barycentre de la face est proche de l'IBC
  DoubleTab typeface(champ) ;// CI
  typeface = -1.e30 ;

  if(variables_internes_-> type_distance_calculee == Transport_Interfaces_FT_Disc_interne::DIST_MODIFIEE )
    {
      //////////////////////////////////////////////////////////////////////////////////////////////
      //////  LANCE DE RAYON
      //////////////////////////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////////////////////////
      //------- Calcul du nombre de fois qu'un rayon local (par direction) traverse l'IBC dans une cellule
      IntTab trav(nb_elem,dim) ;
      if( nb_elem < nb_elem_tot )
        domaine_vf.domaine().creer_tableau_elements(trav, RESIZE_OPTIONS::NOCOPY_NOINIT) ;
      trav = 0 ;

      DoubleTab xe(dim) ;
      xe = 0. ;
      for(int i_elem=0 ; i_elem<nb_elem ; i_elem++)
        {
          for(int i=0 ; i<dim ; i++ )
            xe(i)=xp(i_elem,i) ;
          // Cas des elements diphasiques
          if( indicatrice(i_elem) != 0. && indicatrice(i_elem) != 1. )
            {
              for(int dir = 0 ; dir < dim ; dir++ )
                {
                  int a = elem_faces(i_elem, dir) ;
                  int b = elem_faces(i_elem, dir+dim) ;
                  const double pas = std::fabs(domaine_vdf.distance_face(a,b,dir)) ;
                  int traverse = 0 ;
                  calcul_nb_traverse(xe,pas,dim,dir,maillage,i_elem,traverse) ;
                  trav(i_elem,dir) = traverse ;
                }
            }
          else if( indicatrice(i_elem) == 0. )
            {
              // Cas des elements fluides
              for( int k=0; k<dim; k++)
                trav(i_elem,k) = -1 ;
            }
          else if( indicatrice(i_elem) == 1. )
            {
              // Cas des elements solides
              for( int k=0; k<dim; k++)
                trav(i_elem,k) = -2 ;
            }
        }
      trav.echange_espace_virtuel() ;
      //-----------------------------------------------------------------
      //-----------------------------------------------------------------
      // initialisation des compteurs de faces
      // faces monophasiques (de reference) : cpt_ref
      // faces diphasiques avec une cellule voisine de reference : cpt_vref
      // faces avec indicatrice_face = 0.5 et deux cellules voisines de reference : cpt_halfref
      int cpt_ref = 0 ;
      int cpt_vref = 0 ;
      int cpt_halfref = 0 ;
      DoubleTab typefacejoint ;
      IntList FACEJOINT  ;

      int nbjoints=domaine_vf.nb_joints();
      for(int njoint=0; njoint<nbjoints; njoint++)
        {
          const Joint& joint_temp = domaine_vf.joint(njoint);
          const IntTab& indices_faces_joint = joint_temp.joint_item(Joint::FACE).renum_items_communs();
          const int nb_faces = indices_faces_joint.dimension(0);
          for (int j = 0; j < nb_faces; j++)
            {
              int face_de_joint = indices_faces_joint(j, 1);
              FACEJOINT.add_if_not(face_de_joint) ;
            }
        }
      //-----------------------------------------------------------------
      // Parcourt des faces monophasiques sur le processeur courant
      for( i_face = 0 ; i_face < nfaces ; i_face++ )
        {
          if( indicatrice_face(i_face) == 0. || indicatrice_face(i_face) == 1. )
            {
              typeface(i_face) = indicatrice_face(i_face) ;
              cpt_ref++ ;
            }
        }
      typefacejoint = typeface ;
      typeface.echange_espace_virtuel() ;
      for( int i=0 ; i<FACEJOINT.size() ; i++ )
        {
          typeface(FACEJOINT[i]) = std::max(typefacejoint(FACEJOINT[i]),typeface(FACEJOINT[i])) ;
          if( typefacejoint(FACEJOINT[i]) < 0 && typeface(FACEJOINT[i]) > -1 )
            cpt_ref++;
        }

      // Traitement des faces diphasiques locales ayant une cellule voisine de reference
      for( i_face = 0 ; i_face < nfaces ; i_face++ )
        {
          // les faces avec une indicatrice de 0.5 sont traitees apres
          if( indicatrice_face(i_face) != 0. && indicatrice_face(i_face) != 1. && indicatrice_face(i_face) != 0.5 )
            {
              // Si l'un des deux voisins est monophasique alors i_face est du meme cote
              int reference_trouvee = 0 ;
              int i_voisin = 0 ;
              while( i_voisin < 2 && reference_trouvee == 0)
                {
                  const int elem_voisin = faces_elem(i_face,i_voisin) ;

                  if( elem_voisin > -1 &&
                      ( indicatrice(elem_voisin) == 0. || indicatrice(elem_voisin) == 1. ) )
                    {
                      reference_trouvee = 1 ;
                      typeface(i_face) = indicatrice(elem_voisin) ;
                      cpt_vref++;
                    }
                  i_voisin++ ;
                }
            }
        }
      typefacejoint = typeface ;
      typeface.echange_espace_virtuel() ;
      for( int i=0 ; i<FACEJOINT.size() ; i++ )
        {
          typeface(FACEJOINT[i]) = std::max(typefacejoint(FACEJOINT[i]),typeface(FACEJOINT[i])) ;
          if( typefacejoint(FACEJOINT[i]) < 0 && typeface(FACEJOINT[i]) > -1
              && indicatrice_face(FACEJOINT[i]) != 0. && indicatrice_face(FACEJOINT[i]) != 1.
              && indicatrice_face(FACEJOINT[i]) != 0.5 )
            cpt_vref++;
        }

      // Traitement des faces avec indicatrice_face = 0.5 et deux cellules voisines de reference
      for( i_face = 0 ; i_face < nfaces ; i_face++ )
        {
          if( indicatrice_face(i_face) == 0.5 )
            {
              int element = -1 ;
              if(  faces_elem(i_face,0) > -1 )
                element = faces_elem(i_face,0) ;
              else
                element = faces_elem(i_face,1) ;

              if( indicatrice(element) == 0. || indicatrice(element) == 1. )
                {
                  typeface(i_face) = 0.5 ;
                  cpt_halfref++;
                }
            }
        }
      typefacejoint = typeface ;
      typeface.echange_espace_virtuel() ;
      for( int i=0 ; i<FACEJOINT.size() ; i++ )
        {
          typeface(FACEJOINT[i]) = std::max(typefacejoint(FACEJOINT[i]),typeface(FACEJOINT[i])) ;
          if( typefacejoint(FACEJOINT[i]) < 0 && typeface(FACEJOINT[i]) > -1  && indicatrice_face(FACEJOINT[i]) == 0.5 )
            cpt_halfref++;
        }


      // A ce niveau, il reste a traiter les faces diphasiques ayant
      // deux cellules voisines traversees par l'IBC
      // Pour determiner leur type, on procede de proche en proche
      // a partir des faces deja determinees.
      // On itere le processus tant que l'ensemble des faces
      // n'a pas ete determine.
      // En parallele, il faut que chaque processeur procede
      // simultanement car a la fin de chaque processus
      // l'info est echangee. Ceci explique le "barrier()" (utile?)
      // On s'arrete des que chaque processeur a determine
      // toutes ces faces.
      //-----------------------------------------------------------
      // Initialisation des compteurs de faces
      // nombre de faces restant a determiner par processeur
      int cpt_rest = nfaces - cpt_ref - cpt_vref - cpt_halfref ;
      // nombre de face diphasique determinee au processus n
      int cpt_diph = 0 ;
      // critere d'arret pour le processeur courant
      double proc_stop = -1. ;
      // critere d'arret pour l'ensemble des processeurs
      double stop = -1. ;
      //////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////
      while( stop < 0 )
        {
          if( proc_stop < 0. )
            {

              for( i_face = 0 ; i_face<nfaces ; i_face++)
                {
                  if( typeface( i_face ) < 0 )
                    {
                      const int ori    = orientation[i_face] ;
                      const int elem0  = faces_elem(i_face,0) ;
                      const int elem1  = faces_elem(i_face,1) ;

                      int elem = 0 ;
                      if( elem0 > -1  && elem0 < nb_elem_tot )
                        elem = elem0 ;
                      else if( elem1 > -1  && elem1 < nb_elem_tot )
                        elem = elem1 ;
                      else
                        {
                          Process::Journal() << "erreur dans le choix de l'element "<<finl ;
                          exit() ;
                        }
                      int nb = trav(elem,ori) ;

                      // On regarde la face voisine que l'on va determiner. Elle peut etre virtuelle.
                      int j_face = elem_faces(elem, ori) + elem_faces(elem, ori+dim) - i_face ;

                      if( typeface(j_face) > -1 )
                        {
                          // On modifie typeface(i_face)
                          if( fmod( nb, 2. ) > 0  )
                            {
                              //Cas impair
                              // a) le rayon rencontre une face fluide alors la face est solide
                              // b) le rayon rencontre une face solide alors la face est fluide
                              typeface(i_face) = 1. - typeface(j_face) ;
                            }
                          else
                            {
                              //Cas pair
                              // a) le rayon rencontre une face fluide alors la face est fluide
                              // b) le rayon rencontre une face solide alors la face est solide
                              typeface(i_face) = typeface(j_face) ;
                            }
                          cpt_diph++ ;
                        }
                      else
                        {
                          elem = elem0 + elem1 - elem ;
                          if( elem > -1 && elem < nb_elem_tot )
                            {
                              nb = trav(elem,ori) ;
                              j_face = elem_faces(elem, ori) + elem_faces(elem, ori+dim) - i_face ;

                              if( typeface(j_face) > -1 )
                                {
                                  // On modifie typeface(i_face)
                                  if( fmod( nb, 2. ) > 0  )
                                    {
                                      //Cas impair
                                      // a) le rayon rencontre une face fluide alors la face est solide
                                      // b) le rayon rencontre une face solide alors la face est fluide
                                      typeface(i_face) = 1. - typeface(j_face) ;
                                    }
                                  else
                                    {
                                      //Cas pair
                                      // a) le rayon rencontre une face fluide alors la face est fluide
                                      // b) le rayon rencontre une face solide alors la face est solide
                                      typeface(i_face) = typeface(j_face) ;
                                    }
                                  cpt_diph++ ;
                                }
                            }
                        }
                    }
                }
            }
          typefacejoint = typeface ;
          typeface.echange_espace_virtuel() ;
          for( int i=0 ; i<FACEJOINT.size() ; i++ )
            {
              typeface(FACEJOINT[i]) = std::max(typefacejoint(FACEJOINT[i]),typeface(FACEJOINT[i])) ;
              if( typefacejoint(FACEJOINT[i]) < 0 && typeface(FACEJOINT[i]) > -1
                  && indicatrice_face(FACEJOINT[i]) != 0. && indicatrice_face(FACEJOINT[i]) != 1.
                  && indicatrice_face(FACEJOINT[i]) != 0.5 )
                cpt_diph++ ;
            }

          // A t on determiner toutes les faces reelles du processeurs
          if( cpt_diph == cpt_rest )
            proc_stop = 0. ;

          Process::barrier() ;
          stop = mp_sum( proc_stop ) ;
        }

      // Traitement secondaire des faces tres proche de l'IBC
      for( i_face = 0 ; i_face < nfaces ; i_face++ )
        {
          if( indicatrice_face(i_face) != 0. &&  indicatrice_face(i_face) != 1. )
            {
              int elem_voisin = 0 ;
              int j_face = 0 ;
              double dxa = 0 ;
              double dxb = 0 ;
              double Lref = 0. ;
              const int ori = orientation[i_face] ;
              if( faces_elem(i_face,0) > -1 )
                {
                  j_face = elem_faces(faces_elem(i_face,0), ori) + elem_faces(faces_elem(i_face,0), ori+dim) - i_face ;
                  dxa = std::fabs(domaine_vdf.distance_face(i_face,j_face,ori)) ;
                }
              if( faces_elem(i_face,1) > -1 )
                {
                  j_face = elem_faces(faces_elem(i_face,1), ori) + elem_faces(faces_elem(i_face,1), ori+dim) - i_face ;
                  dxb =  std::fabs(domaine_vdf.distance_face(i_face,j_face,ori)) ;
                }
              if( dxa < dxb )
                {
                  elem_voisin = faces_elem(i_face,1) ;
                  Lref = dxb*dxb ;
                }
              else
                {
                  elem_voisin = faces_elem(i_face,0) ;
                  Lref = dxa*dxa ;
                }

              for( int k=0 ; k<dim ; k++ )
                {
                  if( k != ori )
                    {
                      int a = elem_faces(elem_voisin, k) ;
                      int b = elem_faces(elem_voisin, k+dim) ;
                      double pas = std::fabs(domaine_vdf.distance_face(a,b,k)) ;
                      Lref += 0.25*pas*pas ;
                    }
                }
              const double ratio = std::fabs(distance_interface_faces(i_face))/sqrt(Lref) ;
              const double tol_distface = 1e-3 ;
              const int test1 = ( ratio < tol_distface ) ;
              if( test1 )
                typeface(i_face) = 0.5 ;
            }
        }
      typefacejoint = typeface ;
      typeface.echange_espace_virtuel() ;
      for( int i=0 ; i<FACEJOINT.size() ; i++ )
        {
          typeface(FACEJOINT[i]) = typefacejoint(FACEJOINT[i]) ;
        }


      //////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////
      // Fin du lance de rayon
      //////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////


      //////////////////////////////////////////////////////////////////////////////////////////////
      //---- Correction de la distance entre les faces et l'IBC (hyp: |distance_interface_faces| correcte)
      for( i_face=0 ; i_face<nfaces ; i_face++ )
        {
          if( distance_interface_faces(i_face) > invalid_test )
            dist_face_cor(i_face) = (2.*typeface(i_face) - 1.)*std::fabs(distance_interface_faces(i_face)) ;
          else
            dist_face_cor(i_face) = distance_interface_faces(i_face) ;
        }
    }
  else
    {
      for( i_face=0 ; i_face<nfaces ; i_face++ )
        {
          dist_face_cor(i_face) = distance_interface_faces(i_face) ;
        }
    }
  dist_face_cor.echange_espace_virtuel() ;

  //////////////////////////////////////////////////////////////////////////////////////////////
  ///----- variables pour post-traitement
  for( i_face=0 ; i_face<nfaces ; i_face++ )
    {
      variables_internes_->distance_interface_faces_corrigee.valeur().valeurs()(i_face) = dist_face_cor(i_face) ;
      variables_internes_->distance_interface_faces_difference.valeur().valeurs()(i_face) = dist_face_cor(i_face) - distance_interface_faces(i_face) ;
    }

  if( variables_internes_-> type_vitesse_imposee == Transport_Interfaces_FT_Disc_interne::ANALYTIQUE
      && variables_internes_-> type_distance_calculee == Transport_Interfaces_FT_Disc_interne::DIST_MODIFIEE )
    {
      int cpt_face_fluide_signefaux = 0 ;
      int cpt_face_solide_signefaux = 0 ;
      int cpt_face_diphasique_signefaux = 0 ;
      int cpt_face_diphasique_signemodifiee = 0 ;
      int cpt_face_diphasique_signemodifiee_dnulle = 0 ;
      int nfaces_diphasique = 0 ;
      int nfaces_fluide = 0 ;
      int nfaces_solide = 0 ;

      for( i_face=0 ; i_face<nfaces ; i_face++ )
        {
          if( indicatrice_face( i_face ) == 0. )
            nfaces_fluide++ ;
          else if( indicatrice_face( i_face ) == 1. )
            nfaces_solide++ ;
          else
            nfaces_diphasique++;

          if( indicatrice_face( i_face ) == 0. && dist_face_cor( i_face ) >= 0. )
            {
              cpt_face_fluide_signefaux++ ;
            }
          else if( indicatrice_face( i_face ) == 1. && dist_face_cor( i_face ) <= 0. && dist_face_cor( i_face ) > invalid_test )
            {
              cpt_face_solide_signefaux++ ;
            }
          else if( indicatrice_face( i_face ) != 1. && indicatrice_face( i_face ) != 0. )
            {
              int element = -1 ;
              if(  faces_elem(i_face,0) > -1 )
                element = faces_elem(i_face,0) ;
              else
                element = faces_elem(i_face,1) ;

              if( indicatrice( element ) == 0. && dist_face_cor( i_face ) >= 0. )
                {
                  cpt_face_diphasique_signefaux++ ;
                }
              else if( indicatrice( element ) == 1. && dist_face_cor( i_face ) <= 0. && dist_face_cor( i_face ) > invalid_test )
                {
                  cpt_face_diphasique_signefaux++ ;
                }
              else
                {
                  if( std::fabs(dist_face_cor(i_face) - distance_interface_faces(i_face)) != 0. )
                    {
                      cpt_face_diphasique_signemodifiee++ ;
                    }
                  if( dist_face_cor( i_face ) == 0. )
                    {
                      cpt_face_diphasique_signemodifiee_dnulle++ ;
                    }
                }
            }
        }
      Cerr <<" -----------------------------------------------------"<< finl ;
      Cerr <<" -----------------------------------------------------"<< finl ;
      Cerr <<" -----------------------------------------------------"<< finl ;
      Cerr <<" Nombre de faces : "<< nfaces <<finl ;
      Cerr <<" Nombre de faces fluides : "<<nfaces_fluide<<finl;
      Cerr <<" Nombre de faces fluides avec mauvais signe : "<<cpt_face_fluide_signefaux<<finl;
      Cerr <<" Nombre de faces solides : "<<nfaces_solide<<finl;
      Cerr <<" Nombre de faces solides avec mauvais signe : "<<cpt_face_solide_signefaux<<finl;
      Cerr <<" Nombre de faces diphasiques : "<<nfaces_diphasique<<finl;
      Cerr <<" Nombre de faces diphasiques avec mauvais signe : "<<cpt_face_diphasique_signefaux<<finl;
      Cerr <<" Nombre de faces diphasiques avec signe modifie : "<<cpt_face_diphasique_signemodifiee<<finl;
      Cerr <<" Nombre de faces diphasiques avec signe modifie et distance nulle : "<<cpt_face_diphasique_signemodifiee_dnulle<<finl;
      Cerr <<" -----------------------------------------------------"<< finl ;
      Cerr <<" -----------------------------------------------------"<< finl ;
      Cerr <<" -----------------------------------------------------"<< finl ;
    }

  //on stocke dans v_imp la composante de la vitesse de l'interface correspondant
  //au projete du centre de gravite d'une face sur l'interface
  DoubleTab v_imp(champ) ;
  v_imp = 0. ;
  Cerr << " INITIALISATION DE Vimp a Vsolide pour les faces solides "<<finl ;
  for ( i_face = 0; i_face < nfaces; i_face++)
    {
      const double d = dist_face_cor(i_face) ;
      if ( d>0. || (indicatrice_face(i_face) == 1. && d < -1.e20) )
        {
          for (int j = 0; j < dim; j++)
            {
              const double coord = xv(i_face,j);
              parser_x.setVar(j, coord) ;
              parser_y.setVar(j, coord) ;
              if (dim==3)
                parser_z.setVar(j, coord) ;
            }
          if (zvdf)
            {
              switch(zvdf->orientation(i_face))
                {
                case 0:
                  v_imp(i_face) = parser_x.eval() ;
                  break ;
                case 1:
                  v_imp(i_face) = parser_y.eval() ;
                  break ;
                case 2:
                  v_imp(i_face) = parser_z.eval() ;
                  break ;
                }
            }
        }
    }
  v_imp.echange_espace_virtuel() ; // CI

  if(variables_internes_-> type_vitesse_imposee == Transport_Interfaces_FT_Disc_interne::UNIFORME)
    {
      //si la vitesse est uniforme on peut calculer directement le gradient ici
      Cerr << "Transport_Interfaces_FT_Disc_interne::UNIFORME " << finl;
      double vx=0,vy=0,vz=0;
      int ix=0,iy=0,iz=0;
      for ( i_face = 0; i_face < nfaces; i_face++)
        {
          // Calcul de la composante normale du gradient du champ:
          // gradient normal = (champ - valeur_interface) / distance.

          double d = dist_face_cor(i_face) ;
          if (indicatrice_face(i_face) == indic_phase && d > invalid_test)
            {
              //Vitesse imposee evaluee au centre des faces
              for (int j = 0; j < dim; j++)
                {
                  const double coord = xv(i_face,j);
                  parser_x.setVar(j, coord);
                  parser_y.setVar(j, coord);
                  if (dimension3)
                    parser_z.setVar(j, coord);
                }

              if (zvdf)
                {
                  v_imp(i_face) = 0.;
                  // PL Ajout de la verification que la vitesse est bien uniforme
                  switch(zvdf->orientation(i_face))
                    {
                    case 0:
                      v_imp(i_face) = parser_x.eval() ;
                      check(v_imp,i_face,vx,ix);
                      break ;
                    case 1:
                      v_imp(i_face) = parser_y.eval() ;
                      check(v_imp,i_face,vy,iy);
                      break ;
                    case 2:
                      v_imp(i_face) = parser_z.eval() ;
                      check(v_imp,i_face,vz,iz);
                      break ;
                    }
                  gradient(i_face) = (champ(i_face) - v_imp(i_face)) / d ;
                }
            }
        }
    }
  else if (variables_internes_-> type_vitesse_imposee == Transport_Interfaces_FT_Disc_interne::ANALYTIQUE)
    {
      // Si la vitesse n'est pas uniforme, on projete les points suceptibles
      // d'etre interpoles sur l'interface de facon a calculer un gradient
      // correct
      Cerr << "Transport_Interfaces_FT_Disc_interne::ANALYTIQUE "<< finl ;

      if( variables_internes_->type_indic_faces_ == Transport_Interfaces_FT_Disc_interne::MODIFIEE && variables_internes_-> type_projete_calcule != Transport_Interfaces_FT_Disc_interne::PROJETE_SIMPLIFIE)
        {
          double indic_pos=variables_internes_->modified_indic_faces_position;
          double indic_thi=variables_internes_->modified_indic_faces_thickness;
          if (indic_pos-0.5*indic_thi < -1.)
            {
              Cerr << "-------------------------------------------------------------------------"  << finl;
              Cerr << "-------------------------------------------------------------------------"  << finl;
              Cerr << "The datafile ask to use type_vitesse_imposee=ANALYTIQUE and"                << finl;
              Cerr << "type_indic_faces modifiee with a too large transition to reach zero"        << finl;
              Cerr << "There is a big risk that fluid gradients cannot be computed"                << finl;
              Cerr << "because fluid points (Indic_faces=0) are too far from the interface"        << finl;
              Cerr << "To fix that, just do one of the following things:"                          << finl;
              Cerr << "If the imposed velocity field is uniform, use type_vitesse_imposee=UNIFORM" << finl;
              Cerr << "If not, use \"distance_projete_faces simplifiee\" "                         << finl;
              Cerr << "-------------------------------------------------------------------------"  << finl;
              Cerr << "-------------------------------------------------------------------------"  << finl;
              Process::exit();
            }

        }

      if(variables_internes_-> type_projete_calcule == Transport_Interfaces_FT_Disc_interne::PROJETE_SIMPLIFIE)
        {
          Cerr << "Transport_Interfaces_FT_Disc_interne::PROJETE_SIMPLIFIE, calcul a partir des distances et normales" << finl;

          // Ici, on ne calcule pas de projete avec uzawa, on utilise simplement les quantites deja calculees dans le code
          // on a la distance interface-faces d, donc en determinant une normale aux faces n_F, on peut en deduire
          // le projete x_P de la facon suivante :
          // x_P = x_F - d * N_F
          // on determine la normale aux faces comme la moyenne de la normale aux elements.

          const DoubleTab& normale_elem = get_update_normale_interface().valeurs();
          for (i_face = 0; i_face < nfaces; i_face++)
            {
              double d = dist_face_cor(i_face) ;
              const int diphasique_tmp = (indicatrice_face(i_face) != 0. &&
                                          (indicatrice_face(i_face) != 1. || (indicatrice_face(i_face) == 1. && d<=0. && d> invalid_test)));
              const int diphasique = ( variables_internes_->is_extra_diphasique_solide ? diphasique_tmp : (diphasique_tmp && d <= 0.) ) ;
              const int fluide = ( indicatrice_face(i_face) == 0. ) ;
              const int solide = ( variables_internes_->is_extra_solide && (d>0. || (indicatrice_face(i_face) == 1. && d<-1.e20)) ) ;
              const int monophasique = ( fluide || solide ) ;

              if (diphasique || (monophasique && d> invalid_test ))
                {
                  DoubleTab coord_projete(dim) ;
                  DoubleTab normale_face(dim) ;
                  DoubleTab coord_face(dim) ;
                  coord_projete = 0. ;
                  normale_face  = 0. ;
                  coord_face    = 0. ;
                  double cpt_normale = 0.;
                  // Coordonnee au centre de la face :
                  for (int i=0 ; i<dim ; i++)
                    coord_face(i) = xv(i_face, i) ;
                  // Normale a l IBC au centre de la face
                  // calculee avec la moyenne des normales aux elements voisins
                  for(int i=0 ; i<2; i++)
                    {
                      const int elem_voisin = faces_elem(i_face,i) ;
                      if (elem_voisin > -1)
                        {
                          for (int j = 0; j < dim; j++)
                            normale_face(j) += normale_elem(elem_voisin, j);
                          cpt_normale += 1. ;
                        }
                    }
                  assert(cpt_normale!=0.);
                  normale_face /=cpt_normale;
                  double norme_normale = 0.;
                  for (int i=0 ; i<dim ; i++)
                    norme_normale +=normale_face(i)*normale_face(i);
                  assert(norme_normale!=0.);
                  normale_face /= sqrt(norme_normale);
                  // Coordonnee du projete :
                  for (int i=0 ; i<dim ; i++)
                    coord_projete(i) = coord_face(i) - d*normale_face(i);
                  //v_imp au projete :
                  for (int j = 0; j < dim; j++)
                    {
                      parser_x.setVar(j, coord_projete(j)) ;
                      parser_y.setVar(j, coord_projete(j)) ;
                      if (dim==3)
                        parser_z.setVar(j, coord_projete(j)) ;
                    }
                  if (zvdf)
                    {
                      switch(zvdf->orientation(i_face))
                        {
                        case 0:
                          v_imp(i_face) = parser_x.eval() ;
                          break ;
                        case 1:
                          v_imp(i_face) = parser_y.eval() ;
                          break ;
                        case 2:
                          v_imp(i_face) = parser_z.eval() ;
                          break ;
                        }
                    }
                }
            }
        }
      else
        {

          //------------------------
          //Initialisation
          //------------------------
          // 'exec_planfa7existant = 1' permet d'activer la routine
          // 'plan_facette_existant' dans 'RenumFa7' et 'StockageFa7'.
          // Pour un element donne, cette routine permet d'eviter de
          // selectionner deux facettes de meme plan.
          // Basee sur des criteres geometriques, elle peut conduire dans
          // certains cas a des differences sequentiel/parallele.
          // Son activation est recommandee lorsque le nombre de facettes
          // par element est important et qu'une majorite d'entre elles
          // a le meme plan.
          const int exec_planfa7existant = 0 ;
          int nb_facettes = maillage.nb_facettes() ;
          // dimfa7 :
          // Surement gourmand en allocation memoire pour des maillages fins
          // +1 pour s'assurer qu'en cas de doublon num_max_facette
          // +nb_facettes_dim!=dimfa7
          // La somme des facettes par processeur n'est pas egale au nombre
          //  de facettes total en sequentiel a cause des doublons.
          int dimfa7 = int( Process::mp_sum( double(nb_facettes) ) )+1 ;
          int nb_facettes_dim = int( Process::mp_max( double(nb_facettes) ) ) ;
          // Nombre max de facettes acceptees et calculees par element
          // Le nombre de facette acceptees est choisi des que le maillage
          // Lagrangien est significativement plus fin que le maillage Eulerien.
          // Dans ce cas, on se retrouve avec un nombre consequent de facettes
          // par element (eg. maillage importe d'une CAO)
          int nb_fa7_accepted = variables_internes_->nomb_fa7_accepted;
          int max_fa7 = nb_fa7_accepted ;
          // OutElem : liste des eventuels elements qui contiennent un nombre
          // de facettes superieurs a nb_fa7_accepted.
          // OutFa7 : liste des facettes considerees comme etant les plus
          // representatives de l'IBC.
          IntList OutElem, OutFa7 ;

          // Remplissage des listes OutElem et OutFa7
          calcul_OutElemFa7(maillage,indicatrice,nb_elem,nb_fa7_accepted,OutElem,OutFa7) ;

          // Dans la mesure ou le nombre de facettes par element est inferieur a
          // nb_fa7_accepted alors on calcule le nombre max de facettes afin
          // d'optimiser la dimension des tableaux distribues (voir plus bas)
          if( OutElem.est_vide() )
            calcul_maxfa7(maillage,indicatrice,nb_elem,max_fa7,exec_planfa7existant) ;

          int nb_facettes_in_elem = int( Process::mp_max( double(max_fa7) ) ) ;

          // Tableaux distribues pour le stokage des facettes par elements
          // nb_elem : nb d'elements reels par processeur
          // Tab10i :
          // - entrees i= 0,1,2,3 : parametres du plan d'une facette (a,b,c,d)
          // Tab11i :
          // - entrees i= 0,1,2 : coordonnees du barycentre d'une facette
          // Tab12 : numero local de la facette ( < nb_facettes )
          // utile pour la renumerotation des facettes virtuelles
          DoubleTab Tab100( nb_elem, nb_facettes_in_elem ) ;
          DoubleTab Tab101( nb_elem, nb_facettes_in_elem ) ;
          DoubleTab Tab102( nb_elem, nb_facettes_in_elem ) ;
          DoubleTab Tab103( nb_elem, nb_facettes_in_elem ) ;
          DoubleTab Tab110( nb_elem, nb_facettes_in_elem ) ;
          DoubleTab Tab111( nb_elem, nb_facettes_in_elem ) ;
          DoubleTab Tab112( nb_elem, nb_facettes_in_elem ) ;
          IntTab Tab12( nb_elem, nb_facettes_in_elem ) ;
          // Tableau distribue stockant le nb de facettes par element
          IntTab CptFacette( nb_elem ) ;

          // Tableau distribue pour le calul du sommet de facette
          // le plus proche du barycentre d'une face
          DoubleTab Vertex(nfaces,dim) ;
          DoubleTab PPP(nfaces,dim) ;

          // Distribution des tableaux (uniquement en parallele)
          if( nb_elem < nb_elem_tot )
            {
              domaine_vf.domaine().creer_tableau_elements(Tab100, RESIZE_OPTIONS::NOCOPY_NOINIT) ;
              domaine_vf.domaine().creer_tableau_elements(Tab101, RESIZE_OPTIONS::NOCOPY_NOINIT) ;
              domaine_vf.domaine().creer_tableau_elements(Tab102, RESIZE_OPTIONS::NOCOPY_NOINIT) ;
              domaine_vf.domaine().creer_tableau_elements(Tab103, RESIZE_OPTIONS::NOCOPY_NOINIT) ;
              domaine_vf.domaine().creer_tableau_elements(Tab110, RESIZE_OPTIONS::NOCOPY_NOINIT) ;
              domaine_vf.domaine().creer_tableau_elements(Tab111, RESIZE_OPTIONS::NOCOPY_NOINIT) ;
              domaine_vf.domaine().creer_tableau_elements(Tab112, RESIZE_OPTIONS::NOCOPY_NOINIT) ;
              domaine_vf.domaine().creer_tableau_elements(Tab12, RESIZE_OPTIONS::NOCOPY_NOINIT) ;
              domaine_vf.domaine().creer_tableau_elements(CptFacette, RESIZE_OPTIONS::NOCOPY_NOINIT) ;
              domaine_vf.creer_tableau_faces(Vertex, RESIZE_OPTIONS::NOCOPY_NOINIT) ;
              domaine_vf.creer_tableau_faces(PPP, RESIZE_OPTIONS::NOCOPY_NOINIT) ;
            }
          CptFacette = 0 ;
          Tab100 = -1e+30 ;
          Tab101 = -1e+30 ;
          Tab102 = -1e+30 ;
          Tab103 = -1e+30 ;
          Tab110 = -1e+30 ;
          Tab111 = -1e+30 ;
          Tab112 = -1e+30 ;
          Tab12 = -1 ;
          Vertex = 1.e30 ;
          PPP = 1.e30 ;

          // Tableau des coordonnees barycentrique des facettes
          DoubleTab Barycentre(nb_facettes, 3) ;
          Barycentre = 0. ;

          Cerr << "STOKAGE DES FACETTES " << finl ;
          ArrOfBit fa7(maillage.nb_facettes()) ;
          fa7 = 0 ;
          StockageFa7(maillage,CptFacette,Tab100,Tab101,Tab102,Tab103,
                      Tab110,Tab111,Tab112,Tab12,Barycentre,indicatrice,
                      OutElem,fa7,exec_planfa7existant) ;

          // On traite ici des 'OutElem.size()' elements dont le nombre de
          // facettes est superieur a nb_fa7_accepted.
          if( !OutElem.est_vide() )
            {
              if( OutFa7.size() != OutElem.size()*nb_fa7_accepted )
                {
                  Cerr << "ERREUR DE DIMENSIONNEMENT " << finl ;
                  Cerr << " OutFa7.size() "<<OutFa7.size()<<finl ;
                  Cerr << " OutElem.size()*nb_fa7_accepted "<<OutElem.size()*nb_fa7_accepted<<finl ;
                  exit() ;
                }

              // Tableau des numeros de facettes retenues
              // dans les 'OutElem.size()' elements.
              IntTab TabOutFa7(OutElem.size(),nb_fa7_accepted) ;
              for( int i=0 ; i< OutElem.size() ; i++)
                {
                  // Chaque element contient nb_fa7_accepted facettes.
                  CptFacette( OutElem[i] ) = nb_fa7_accepted ;
                  for( int j=0 ; j< nb_fa7_accepted ; j++)
                    {
                      TabOutFa7(i,j) = OutFa7[j+i*nb_fa7_accepted] ;
                    }
                }
              StockageFa7(maillage,Tab100,Tab101,Tab102,Tab103,Tab110,Tab111,Tab112,
                          Tab12,Barycentre,OutElem,TabOutFa7,fa7) ;
            }
          // Echange entre les processeurs
          Tab100.echange_espace_virtuel() ;
          Tab101.echange_espace_virtuel() ;
          Tab102.echange_espace_virtuel() ;
          Tab103.echange_espace_virtuel() ;
          Tab110.echange_espace_virtuel() ;
          Tab111.echange_espace_virtuel() ;
          Tab112.echange_espace_virtuel() ;
          Tab12.echange_espace_virtuel() ;
          CptFacette.echange_espace_virtuel() ;
          // Numerotation locale des facettes virtuelles (uniquement en parallele)
          if( nb_elem < nb_elem_tot )
            {
              Cerr << " RENUMEROTATION DES FACETTES VIRTUELLES " << finl ;
              RenumFa7(Barycentre,Tab110,Tab111,Tab112,Tab12,CptFacette,
                       nb_facettes,nb_facettes_dim) ;
            }

          if(variables_internes_-> type_projete_calcule == Transport_Interfaces_FT_Disc_interne::PROJETE_MODIFIE)
            {
              // Demarche :
              // 1) chaque face diphasique appartient au moins a un element diphasique
              // 1) on identifie donc dans ces elements voisins diphasique
              // 1) le sommet de facette le plus proche et on le stocke dans
              // 1) Vertex.
              //
              // 2) On parcourt a nouveau les faces diphasiques mais cette fois en regardant les
              // 2) elements voisins des elements voisins. Le sommet le plus proche definitif
              // 2) pour ces faces est stocke dans PPP.
              //
              // 3) On parcourt enfin les faces fluides, on determine leurs faces voisines diphasiques
              // 3) et a partir de ces dernieres on determine le sommet le plus proche.
              // 3) qui est stocke dans PPP.
              //
              Cerr << " IDENTIFICATICATION DU VERTEX LE PLUS PROCHE POUR LES FACES CONCERNEES PAR LA PROJECTION " << finl ;
              Cerr << " CAS DES FACES APPARTENANT A AU MOINS UN ELEMENT TRAVERSE : ETAPE 1 " <<finl ;
              PPP_face_interface(maillage,indicatrice,indicatrice_face,Vertex) ;
              MD_Vector_tools::echange_espace_virtuel(Vertex, MD_Vector_tools::EV_MINCOL1) ;
              Cerr << " CAS DES FACES APPARTENANT A AU MOINS UN ELEMENT TRAVERSE : ETAPE 2 " <<finl ;
              PPP_face_interface_voisin(indicatrice,indicatrice_face,Vertex,PPP) ;
              PPP.echange_espace_virtuel() ;
              Cerr << " CAS DES FACES MONOPHASIQUE : ETAPE 3 " <<finl ;
              PPP_face_voisin(indicatrice,indicatrice_face,PPP) ;
              PPP.echange_espace_virtuel() ;
              Cerr << " FIN IDENTIFICATION " << finl ;
            }
          // Nombre de point projete modifie a partir d'un critere de distance
          int nb_proj_modif = 0 ;
          int nb_proj_modif_tot = 0 ;
          Cerr << "CALCUL DU PROJETE POUR LES FACES QUI APPARTIENNENT " << finl ;
          Cerr << "A UN ELEMENT TRAVERSE PAR L INTERFACE" << finl;
          projete_point_face_interface(nb_proj_modif,dimfa7,indicatrice_face,indicatrice,dist_face_cor,t,dt,Tab100,
                                       Tab101,Tab102,Tab103,Tab12,CptFacette,v_imp,PPP,parser_x,parser_y,parser_z) ;
          nb_proj_modif_tot = mp_sum(nb_proj_modif) ;
          if( Process::je_suis_maitre() )
            Cerr << " Nombre de projete modifie pour les faces diphasiques : "<<nb_proj_modif_tot<<" au temps "<<schema_temps().temps_courant()<<finl ;

          //----------------------------------------------------------------
          //----------------------------------------------------------------
          Cerr << "CALCUL DU PROJETE POUR LES FACES FLUIDE " << finl ;
          projete_point_face_fluide(nb_proj_modif,dimfa7,indicatrice_face,indicatrice,dist_face_cor,
                                    t,dt,Tab100,Tab101,Tab102,Tab103,Tab12,CptFacette,v_imp,PPP,parser_x,parser_y,parser_z) ;
          v_imp.echange_espace_virtuel() ;
          nb_proj_modif_tot = mp_sum(nb_proj_modif) ;
          if( Process::je_suis_maitre() )
            Cerr << " Nombre de projete modifie pour les faces fluides : "<<nb_proj_modif_tot<<" au temps "<<schema_temps().temps_courant()<<finl ;
          Debog::verifier("Transport_Interfaces_FT_Disc::interpoler_vitesse_face vimp",v_imp);
        }
    }

  for (i_face = 0; i_face < nfaces; i_face++)
    {
      double d = dist_face_cor(i_face) ;
      if (indicatrice_face(i_face) == indic_phase && d > invalid_test)
        {
          gradient(i_face) = (champ(i_face) - v_imp(i_face))/d  ;
        }
    }
  gradient.echange_espace_virtuel();
  Debog::verifier("Transport_Interfaces_FT_Disc::interpoler_vitesse_face gradient avant etalement",gradient);

  // On etale ce gradient par continuite sur une epaisseur de "stencil_value"
  // Les iterations de cet algorithme convergent vers une sorte de laplacien=0
  // avec condition aux limites de Dirichlet sur les elements de la phase
  // "phase".
  DoubleTab gradient_old;
  Cerr << "Nombre d'iterations de la methode interpoler_vitesse_face = " << stencil_width << finl;
  for (int iteration = 0; iteration < stencil_width; iteration++)
    {
      // Copie de la valeur du gradient: on ne veut pas utiliser les valeurs
      // calculees lors de l'iteration courante
      gradient_old = gradient;
      // La valeur sur une face est la moyenne des valeurs sur les faces
      // les plus proches (interpolation 'en croix' : a partir de 4 faces
      // voisines en 2D et 6 faces voisines en 3D dans le cas general)
      for (i_face = 0; i_face < nfaces; i_face++)
        {
          const double d = dist_face_cor(i_face) ;
          if (indicatrice_face(i_face) != indic_phase && d > invalid_test)
            {
              // Ne pas toucher au gradient de la phase "phase".
              // Iterer sur les autres valeurs.
              double somme = 0.;
              double coeff = 0.;
              const int ori=orientation(i_face);
              const int elem0 = faces_elem(i_face, 0);
              const int elem1 = faces_elem(i_face, 1);
              int elem=-1;
              if (elem0 >= 0)
                elem=elem0;
              else
                elem=elem1;

              for (int direction=0; direction < dimension; direction++)
                {
                  if (direction==ori)
                    {
                      //Dans la direction de la composante de i_face, on ajoute la contribution des 2 faces qui ont un element commun avec i_face
                      if (elem0 >= 0)
                        {
                          const int voisin0 = elem_faces(elem0, ori) + elem_faces(elem0, ori+dimension) - i_face;
                          if( (elem_faces(elem0, ori) != i_face) && (elem_faces(elem0, ori+dimension) != i_face) )
                            {
                              Cerr << "Attention, il y a un probleme sur la numerotation des faces. On calcule en effet" << finl ;
                              Cerr << "elem_faces(elem0, ori) = " <<elem_faces(elem0, ori)<< finl ;
                              Cerr << "elem_faces(elem0, ori+dimension) = " <<elem_faces(elem0, ori+dimension)<< finl ;
                              Cerr << "i_face = " <<i_face<< finl;
                              Cerr << "tandis qu'au moins un des deux calculs 'elem_faces(elem0, ori)' ou 'elem_faces(elem0, ori+dimension)'"<<finl ;
                              Cerr << "doit correspondre au nuemro de la face i_face"<<finl ;
                              Cerr << "Attention, verifier les conditions aux limites."<<finl ;
                              Cerr << "En effet, le codage suivant ne permet pas de prendre en compte des conditions" << finl ;
                              Cerr << "de periodicite avec une vitesse solide imposee analytique." << finl ;
                              exit() ;
                            }
                          const double grad0 = gradient_old(voisin0);
                          if (grad0 > invalid_test)
                            {
                              somme += grad0;
                              coeff++;
                            }
                        }
                      if (elem1 >= 0)
                        {
                          const int voisin1 = elem_faces(elem1, ori) + elem_faces(elem1, ori+dimension) - i_face;
                          if( (elem_faces(elem1, ori) != i_face) && (elem_faces(elem1, ori+dimension) != i_face) )
                            {
                              Cerr << "Attention, il y a un probleme sur la numerotation des faces. On calcule en effet" << finl ;
                              Cerr << "elem_faces(elem1, ori) = " <<elem_faces(elem1, ori)<< finl ;
                              Cerr << "elem_faces(elem1, ori+dimension) = " <<elem_faces(elem1, ori+dimension)<< finl ;
                              Cerr << "i_face = " <<i_face<< finl;
                              Cerr << "tandis qu'au moins un des deux calculs 'elem_faces(elem1, ori)' ou 'elem_faces(elem1, ori+dimension)'"<<finl ;
                              Cerr << "doit correspondre au nuemro de la face i_face"<<finl ;
                              Cerr << "Attention, verifier les conditions aux limites."<<finl ;
                              Cerr << "En effet, le codage suivant ne permet pas de prendre en compte des conditions" << finl ;
                              Cerr << "de periodicite avec une vitesse solide imposee analytique." << finl ;
                              exit() ;
                            }
                          const double grad1 = gradient_old(voisin1);
                          if (grad1 > invalid_test)
                            {
                              somme += grad1;
                              coeff++;
                            }
                        }
                    }
                  else
                    {
                      //Dans les autres directions, on ajoute les contributions des faces appartenant aux elements voisins, si elles existent:

                      //-on identifie les faces de elem (=un des deux elements auxquels appartient i_face) de direction differente de i_face
                      //-on identifie les elements auxquels appartiennent ces faces
                      //-pour chaque element voisin, on veut trouver la face la plus proche de i_face : on compare la distance entre les centres des faces
                      const int face0=elem_faces(elem,direction);
                      const int face1=elem_faces(elem,direction+dimension);
                      const int elem_voisin0=faces_elem(face0,0)+faces_elem(face0,1)-elem;
                      const int elem_voisin1=faces_elem(face1,0)+faces_elem(face1,1)-elem;

                      double a_carre=0.;
                      double b_carre=0.;
                      double grad_voisin0=0.;
                      double c_carre=0.;
                      double d_carre=0.;
                      double grad_voisin1=0.;

                      if (elem_voisin0 >= 0)
                        {
                          for (int compo = 0; compo<dimension; compo++)
                            {
                              double c_gravite=xv(i_face,compo);
                              double c_face_voisin1=xv(elem_faces(elem_voisin0,ori),compo);
                              double c_face_voisin2=xv(elem_faces(elem_voisin0,ori+dimension),compo);
                              a_carre+=(c_gravite-c_face_voisin1)*(c_gravite-c_face_voisin1);
                              b_carre+=(c_gravite-c_face_voisin2)*(c_gravite-c_face_voisin2);
                            }
                          if (a_carre<b_carre)
                            grad_voisin0=gradient_old(elem_faces(elem_voisin0,ori));
                          else
                            grad_voisin0=gradient_old(elem_faces(elem_voisin0,ori+dimension));
                          if (grad_voisin0 > invalid_test)
                            {
                              somme += grad_voisin0;
                              coeff++;
                            }
                        }
                      if (elem_voisin1 >= 0)
                        {
                          for (int compo = 0; compo<dimension; compo++)
                            {
                              double c_gravite=xv(i_face,compo);
                              double c_face_voisin3=xv(elem_faces(elem_voisin1,ori),compo);
                              double c_face_voisin4=xv(elem_faces(elem_voisin1,ori+dimension),compo);
                              c_carre+=(c_gravite-c_face_voisin3)*(c_gravite-c_face_voisin3);
                              d_carre+=(c_gravite-c_face_voisin4)*(c_gravite-c_face_voisin4);
                            }
                          if (c_carre<d_carre)
                            grad_voisin1=gradient_old(elem_faces(elem_voisin1,ori));
                          else
                            grad_voisin1=gradient_old(elem_faces(elem_voisin1,ori+dimension));
                          if (grad_voisin1 > invalid_test)
                            {
                              somme += grad_voisin1;
                              coeff++;
                            }
                        }
                    }
                }

              if (coeff > 0.)
                gradient(i_face) = somme / coeff;
            }
        }
      gradient.echange_espace_virtuel();
    }

  Debog::verifier("Transport_Interfaces_FT_Disc::interpoler_vitesse_face gradient",gradient);
  //    Debog::verifier("Transport_Interfaces_FT_Disc::interpoler_vitesse_face distance_interface_faces",distance_interface_faces);
  // On calcule la valeur extrapolee:
  for (i_face = 0; i_face < nfaces; i_face++)
    {
      if (indicatrice_face(i_face) != indic_phase)
        {
          double d = dist_face_cor(i_face) ;
          const double grad = gradient[i_face];
          const double indic = indicatrice_face(i_face) ;

          const int interpolation    = ( d < 0. ) ;
          const int extra_diphasique = ( variables_internes_->is_extra_diphasique_solide && indic != 0. && indic != 1. && d > 0. ) ;
          const int extra_solide     = ( variables_internes_->is_extra_solide && indic == 1. && (d>0. || d < -1.e20) ) ;
          const int extrapolation    = ( extra_diphasique || extra_solide ) ;
          if (variables_internes_-> type_vitesse_imposee == Transport_Interfaces_FT_Disc_interne::UNIFORME)
            {
              //Vitesse imposee evaluee au centre des faces
              for (int j = 0; j < dim; j++)
                {
                  const double coord = xv(i_face,j);
                  parser_x.setVar(j, coord);
                  parser_y.setVar(j, coord);
                  if (dimension3)
                    parser_z.setVar(j, coord);
                }

              if (zvdf)
                {
                  v_imp(i_face) = 0. ;
                  switch(zvdf->orientation(i_face))
                    {
                    case 0:
                      v_imp(i_face) = parser_x.eval();
                      break;
                    case 1:
                      v_imp(i_face) = parser_y.eval();
                      break;
                    case 2:
                      v_imp(i_face) = parser_z.eval();
                      break;
                    }
                  champ( i_face ) = v_imp(i_face) ;


                  if (grad > invalid_test && d > invalid_test  && ( interpolation || extrapolation ) )
                    {
                      champ(i_face) += d*grad ;
                    }
                }
            }
          else
            {
              // Vitesse non uniforme : il faut projeter utiliser le tableau v_imp
              // qui contient la vitesse du projete de i_face sur l'interface

              champ(i_face) = v_imp(i_face) ;
              if (grad > invalid_test && d > invalid_test  && ( interpolation || extrapolation ) )
                {
                  champ(i_face) += d * grad ;
                }
            }
          vitesse_imp_interp_.valeur().valeurs()(i_face)= champ( i_face ) ;
        }
      else
        {
          vitesse_imp_interp_.valeur().valeurs()(i_face)=-3e30;
        }
    }
  champ.echange_espace_virtuel() ;
  Debog::verifier("Transport_Interfaces_FT_Disc::interpoler_vitesse_face champ",champ );

}

// Calcul le nombre de fois que l'on traverse l'IBC
void Transport_Interfaces_FT_Disc::calcul_nb_traverse( const DoubleTab& xe, const double dx,
                                                       const int dim, const int ori,
                                                       Maillage_FT_Disc& maillage, int elem,
                                                       int& traverse )

{
  const Intersections_Elem_Facettes& intersection = maillage.intersections_elem_facettes();
  const ArrOfInt& index_elem = intersection.index_elem() ;
  const DoubleTab& normale_facettes = maillage.get_update_normale_facettes() ;
  const IntTab& facettes = maillage.facettes() ;
  //const int nb_som  = facettes.line_size() ;
  const DoubleTab& sommets = maillage.sommets() ;
  //IntTab Som( nb_som ) ;
  // Pour un element donne, on parcourt ces facettes
  int index = index_elem[elem] ;
  const double precision = Objet_U::precision_geom ;
  const double tol_theta = 1e-4 ;
  double I_sur_arete_fa7 = 0. ;
  double I_sur_arete_vertex = 0. ;
  int I_dans_fa7 = 0 ;

  while (index >= 0 )
    {
      const Intersections_Elem_Facettes_Data& data = intersection.data_intersection(index);
      const int i_facette = data.numero_facette_;

      // produit scalaire entre un vecteur directeur au rayon et la normale a la facette
      // le vecteur de reference est unitaire valant un sur la composante ori.
      // Si le produit scalaire est nul alors il n'existe pas de I dans la facette
      if( std::fabs( normale_facettes(i_facette,ori) ) > 0. )
        {
          DoubleTab A(dim), B(dim), C(dim), U(dim), I(dim);
          A = 0. ;
          B = 0. ;
          C = 0. ;
          U = 0.  ;
          I = 0. ;
          for(int i = 0 ; i<dim ; i++)
            {
              I(i) = xe(i) ;
            }
          I(ori) = 0. ;
          for(int i=0 ; i<dim ; i++)
            {
              U(i) = normale_facettes(i_facette,i) ;
              A(i) = sommets(facettes(i_facette, 0),i) ;
              B(i) = sommets(facettes(i_facette, 1),i) ;
              if( dim == 3 )
                C(i) = sommets(facettes(i_facette, 2),i) ;
            }

          for( int i =0 ; i<dim ; i++)
            {
              if( i!=ori )
                I(ori) += (A(i) - I(i))*U(i)/U(ori) ;
            }
          I(ori) += A(ori) ;


          double FIscalU = (I(ori)-xe(ori))*U(ori) ;
          double FI = std::fabs( I(ori)-xe(ori) ) ;
          double theta = acos(FIscalU/FI) ;
          double theta_rel = std::fabs( theta - (2.0*atan(1.)) ) / (2.0*atan(1.)) ;
          double ecart_elem = std::fabs(xe(ori)-I(ori))  ;

          if( theta_rel >= tol_theta && ecart_elem <= 0.5*dx  )
            {
              if( dim == 2 )
                {
                  double AIscalAB = (I(0)-A(0))*(B(0)-A(0)) + (I(1)-A(1))*(B(1)-A(1)) ;
                  double ABscalAB = (B(0)-A(0))*(B(0)-A(0)) + (B(1)-A(1))*(B(1)-A(1)) ;
                  double xx = AIscalAB / ABscalAB ;

                  const int test1 = ( std::fabs( xx ) <= precision ) ;
                  const int test2 = ( std::fabs( xx  - 1. ) <= precision ) ;
                  const int test3 = ( xx > precision ) ;
                  const int test4 = ( xx < 1.-precision ) ;
                  if( test1 || test2 )
                    {
                      // si I est dans un voisinage d'un vertex
                      I_sur_arete_vertex += 1./2. ;
                    }
                  else if( test3 && test4 )
                    {
                      // I est a l'interieur du triangle
                      I_dans_fa7++;
                    }
                }
              else if( dim == 3 )
                {
                  DoubleTab ABvectAC(dim) ;
                  ABvectAC(0) = (B(1)-A(1))*(C(2)-A(2))-(B(2)-A(2))*(C(1)-A(1)) ;
                  ABvectAC(1) = (B(2)-A(2))*(C(0)-A(0))-(B(0)-A(0))*(C(2)-A(2)) ;
                  ABvectAC(2) = (B(0)-A(0))*(C(1)-A(1))-(B(1)-A(1))*(C(0)-A(0)) ;
                  double ABvectACscalU = ABvectAC(0)*U(0) +  ABvectAC(1)*U(1) +  ABvectAC(2)*U(2) ;
                  double AireABC = ABvectACscalU / 2. ;

                  DoubleTab IBvectIC(dim) ;
                  IBvectIC(0) = (B(1)-I(1))*(C(2)-I(2))-(B(2)-I(2))*(C(1)-I(1)) ;
                  IBvectIC(1) = (B(2)-I(2))*(C(0)-I(0))-(B(0)-I(0))*(C(2)-I(2)) ;
                  IBvectIC(2) = (B(0)-I(0))*(C(1)-I(1))-(B(1)-I(1))*(C(0)-I(0)) ;
                  double IBvectICscalU = IBvectIC(0)*U(0) +  IBvectIC(1)*U(1) +  IBvectIC(2)*U(2) ;
                  double AireIBC = IBvectICscalU / 2. ;

                  DoubleTab ICvectIA(dim) ;
                  ICvectIA(0) = (C(1)-I(1))*(A(2)-I(2))-(C(2)-I(2))*(A(1)-I(1)) ;
                  ICvectIA(1) = (C(2)-I(2))*(A(0)-I(0))-(C(0)-I(0))*(A(2)-I(2)) ;
                  ICvectIA(2) = (C(0)-I(0))*(A(1)-I(1))-(C(1)-I(1))*(A(0)-I(0)) ;
                  double ICvectIAscalU = ICvectIA(0)*U(0) +  ICvectIA(1)*U(1) +  ICvectIA(2)*U(2) ;
                  double AireICA = ICvectIAscalU / 2. ;

                  DoubleTab IAvectIB(dim) ;
                  IAvectIB(0) = (A(1)-I(1))*(B(2)-I(2))-(A(2)-I(2))*(B(1)-I(1)) ;
                  IAvectIB(1) = (A(2)-I(2))*(B(0)-I(0))-(A(0)-I(0))*(B(2)-I(2)) ;
                  IAvectIB(2) = (A(0)-I(0))*(B(1)-I(1))-(A(1)-I(1))*(B(0)-I(0)) ;
                  double IAvectIBscalU = IAvectIB(0)*U(0) +  IAvectIB(1)*U(1) +  IAvectIB(2)*U(2) ;
                  double AireIAB = IAvectIBscalU / 2. ;

                  double alpha = AireIBC/AireABC ;
                  double beta  = AireICA/AireABC ;
                  double gamma = AireIAB/AireABC ;

                  int test1 = ( std::fabs(alpha-1.)<=precision && std::fabs(beta)<=precision && std::fabs(gamma)<=precision ) ;
                  int test2 = ( std::fabs(beta-1.)<=precision && std::fabs(alpha)<=precision && std::fabs(gamma)<=precision ) ;
                  int test3 = ( std::fabs(gamma-1.)<=precision && std::fabs(alpha)<=precision && std::fabs(beta)<=precision ) ;
                  int test4 = ( std::fabs(alpha)<=precision && precision<beta && beta<1.-precision && precision<gamma && gamma<1.-precision ) ;
                  int test5 = ( std::fabs(beta)<=precision && precision<alpha && alpha<1.-precision && precision<gamma && gamma<1.-precision ) ;
                  int test6 = ( std::fabs(gamma)<=precision && precision<beta && beta<1.-precision && precision<alpha && alpha<1.-precision ) ;
                  int test7 = ( precision<alpha && alpha<1.-precision) ;
                  int test8 = ( precision<beta && beta<1.-precision) ;
                  int test9 = ( precision<gamma && gamma<1.-precision) ;

                  if( test1 || test2 || test3 )
                    {
                      // si I est dans un voisinage d'un vertex
                      // Dans ce cas, une "infinite" (nombre max de facette dans elem) de facettes peuvent etre concernee.
                      I_sur_arete_vertex += 1e-5 ;
                    }
                  else if( test4 || test5 || test6 )
                    {
                      // si I est dans un voisinage d'une arete sans etre un vertex
                      I_sur_arete_fa7 += 1./2. ;
                    }
                  else if( test7 && test8 && test9 )
                    {
                      // I est a l'interieur
                      I_dans_fa7++;
                    }
                }
              else
                {
                  Cerr << " attention seules les dimension 2 et 3 sont traitees "<<finl ;
                  exit() ;
                }
            }
        }
      // Iteration des facettes de l'element
      index = data.index_facette_suivante_;
    }
  traverse = int(ceil(I_sur_arete_vertex)) + int(ceil(I_sur_arete_fa7)) + I_dans_fa7 ;
}

// Methode listant les elements contenant un nombre de facettes superieur a
// nb_fa7_accepeted et l'ensemble des facettes associees.
void Transport_Interfaces_FT_Disc::calcul_OutElemFa7( Maillage_FT_Disc& maillage, const DoubleTab& indicatrice,
                                                      const int nb_elem, int& nb_fa7_accepted,
                                                      IntList& OutElem, IntList& OutFa7 )
{
  const Intersections_Elem_Facettes& intersection = maillage.intersections_elem_facettes();
  const ArrOfInt& index_elem = intersection.index_elem() ;
  // Calcul de la surface des facettes contenues dans i_elem
  const ArrOfDouble& surfaces = maillage.get_update_surface_facettes();

  // Parcourt de l'ensemble des elements reels
  for (int i_elem = 0; i_elem < nb_elem ; i_elem++)
    {
      // Test si l'element est traverse par une IBC
      if( indicatrice(i_elem) != 0. && indicatrice(i_elem) != 1. )
        {
          int index = index_elem[i_elem] ;
          // Listes contenant la surface des facettes et le numero de la facette
          // Les 'nb_fa7_accepted' facettes de plus grande surface sont
          // considerees ici comme etant les plus representatives.
          DoubleList E ;
          IntList F ;

          // Parcours des facettes traversant l'element i_elem
          while (index >= 0 )
            {
              const Intersections_Elem_Facettes_Data& data = intersection.data_intersection(index);
              const int i_facette = data.numero_facette_;

              // Stokage de la surface de la facette i_facette dans la liste E
              E.add(surfaces[i_facette]) ;
              // Stokage du numero de la facette dans la liste F
              F.add(i_facette) ;
              index = data.index_facette_suivante_;
            }

          // Si le nombre de facette est superieur a nb_fa7_accepted
          // alors on choisit les facettes les plus representatives.
          // Ici, representative := de plus grande surface
          if( F.size() > nb_fa7_accepted )
            {
              // Stockage de l'element i_elem dans la liste OutElem
              OutElem.add_if_not( i_elem ) ;
              // Liste des facettes stockees dans i_elem
              IntList Fa7Elem ;

              // Listing des facettes de plus grande surface
              int i = 0 ;
              while( i<nb_fa7_accepted )
                {
                  // Recherche de la surface max
                  double max_data = -1e30 ;
                  int l=0 ;
                  int ll = 0 ;
                  while( l<E.size() )
                    {
                      if( E[l]>max_data && !Fa7Elem.contient(F[l]) )
                        {
                          max_data = E[l] ;
                          ll = l ;
                        }
                      l++ ;
                    }
                  // Ajout de la facette de plus grande surface dans
                  // la liste OutFa7
                  Fa7Elem.add( F[ll] ) ;
                  OutFa7.add( F[ll] ) ;
                  i++ ;
                }
            }
        }
    }
}

// Calcul du vertex le plus proche pour les faces diphasiques ( etape 1 )
void Transport_Interfaces_FT_Disc::PPP_face_interface( Maillage_FT_Disc& maillage, const DoubleTab& indicatrice,
                                                       const DoubleTab& indicatrice_face, DoubleTab& Vertex )

{
  const int dim = Objet_U::dimension;
  const Domaine_dis_base& mon_dom_dis = domaine_dis().valeur();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, mon_dom_dis);
  const IntTab& faces_elem = domaine_vf.face_voisins();
  const int nb_elem = domaine_vf.nb_elem() ;
  const int nfaces = faces_elem.dimension(0) ;
  const Intersections_Elem_Facettes& intersection = maillage.intersections_elem_facettes();
  const ArrOfInt& index_elem = intersection.index_elem() ;
  const IntTab& facettes = maillage.facettes() ;
  const DoubleTabFT& sommets = maillage.sommets() ;
  const DoubleTab& xv = domaine_vf.xv();
  DoubleTab coord(dim) ;
  coord = 0. ;

  for( int i_face=0 ; i_face<nfaces ; i_face++ )
    {
      // Distance face facette
      double dist_fa7 = 1e+30 ;
      int ppp_calcule = 0 ;

      // Parcourt des faces diphasiques
      if( indicatrice_face(i_face) == 0.5  )
        {
          int voisin = -1 ;
          if( faces_elem(i_face,0) > -1 )
            voisin = faces_elem(i_face,0) ;
          else
            voisin = faces_elem(i_face,1) ;

          // Dans ce cas, i_face coincide avec les facettes
          // Le point le plus proche est lui-meme
          if( indicatrice(voisin) == 0. || indicatrice(voisin) == 1. )
            {
              for( int j=0; j<dim; j++)
                {
                  Vertex(i_face,j) = xv(i_face,j) ;
                }
              ppp_calcule = 1 ;
            }
        }

      if (indicatrice_face(i_face) != 0. && indicatrice_face(i_face) != 1. && ppp_calcule == 0 )
        {
          // Parcourt des elements voisins traverses physiques
          for(int i=0 ; i<2; i++)
            {
              const int elem_voisin = faces_elem(i_face,i) ;

              if( elem_voisin<nb_elem && elem_voisin>-1 && indicatrice(elem_voisin) != 0. && indicatrice(elem_voisin) != 1. )
                {
                  int index = index_elem[elem_voisin] ;
                  // Parcourt des facettes de l'element
                  while (index >= 0 )
                    {
                      const Intersections_Elem_Facettes_Data& data = intersection.data_intersection(index);
                      const int i_facette = data.numero_facette_;
                      // Distance min par elements
                      double dist_min = 1e+30 ;
                      // Parcourt des sommets des facettes
                      for( int som=0; som<facettes.dimension(1); som++)
                        {
                          // Distance face sommet
                          double dist_som = 0. ;
                          for( int j=0; j<dim; j++)
                            {
                              // j-ieme coordonnees du sommets som
                              double sj = sommets(facettes(i_facette, som), j) ;
                              // j-ieme coordonnees du barycentre de la face
                              double xj = xv(i_face, j) ;
                              // Carre de la distance face sommet som
                              dist_som += (xj - sj)*(xj - sj) ;
                            }
                          // Distance face sommet som
                          dist_som = sqrt(dist_som) ;

                          // Distance face facette par rapport a distance face sommet de facette
                          if( dist_som < dist_min )
                            {
                              dist_min = dist_som ;
                              for( int j=0; j<dim; j++)
                                {
                                  coord(j) = sommets(facettes(i_facette, som), j) ;
                                }
                            }
                        }
                      // Distance face sommet minimale par rapport a distance face facette
                      if( dist_min < dist_fa7 )
                        {
                          dist_fa7 = dist_min ;
                          for( int j=0; j<dim; j++)
                            {
                              Vertex(i_face,j) = coord(j) ;
                            }
                        }
                      // Ieration des facettes de l'element
                      index = data.index_facette_suivante_;
                    }
                }
            }
        }
    }
}

// Calcul du vertex le plus proche des faces diphasiques (etape 2)
void Transport_Interfaces_FT_Disc::PPP_face_interface_voisin( const DoubleTab& indicatrice, const DoubleTab& indicatrice_face,
                                                              DoubleTab& Vertex, DoubleTab& PPP )

{
  const int dim = Objet_U::dimension;
  const Domaine_dis_base& mon_dom_dis = domaine_dis().valeur();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, mon_dom_dis);
  const IntTab& faces_elem = domaine_vf.face_voisins();
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const int nfaces = faces_elem.dimension(0) ;
  const DoubleTab& xv = domaine_vf.xv();

  for( int i_face=0 ; i_face<nfaces ; i_face++ )
    {
      double d = 0. ;
      double dmin = 0. ;
      int ppp_calcule = 0 ;

      // Parcourt des faces diphasiques
      if( indicatrice_face(i_face) == 0.5  )
        {
          int voisin = -1 ;
          if( faces_elem(i_face,0) > -1 )
            voisin = faces_elem(i_face,0) ;
          else
            voisin = faces_elem(i_face,1) ;

          // Dans ce cas, i_face coincide avec les facettes
          // Le point le plus proche est lui-meme
          // On ne modifie pas Vertex
          if( indicatrice(voisin) == 0. || indicatrice(voisin) == 1. )
            {
              for( int j=0; j<dim; j++)
                {
                  PPP(i_face,j) = xv(i_face,j) ;
                }
              ppp_calcule = 1 ;
            }
        }

      if (indicatrice_face(i_face) != 0. && indicatrice_face(i_face) != 1. && ppp_calcule == 0 )
        {
          for( int i=0 ; i<dim ; i++)
            {
              d += (Vertex(i_face,i)-xv(i_face,i))*(Vertex(i_face,i)-xv(i_face,i)) ;
            }
          d = sqrt(d) ;
          dmin = d ;
          for( int k=0 ; k<dim ; k++)
            {
              PPP(i_face,k) = Vertex(i_face,k) ;
            }
          // On parcourt les autres faces de mes deux voisins
          // On compare le PPP de i_face avec ceux des autres faces.
          // Parmi eux, on choisit le plus proche de i_face.
          for( int i=0; i<2; i++ )
            {
              const int elem = faces_elem(i_face,i) ;
              // Si cet element voisin existe
              if( elem >-1 )
                {
                  // alors on parcourt ses autres faces
                  for (int face_loc=0 ; face_loc<2*dim ; face_loc++)
                    {
                      const int autre_face = elem_faces(elem,face_loc) ;
                      if (indicatrice_face(autre_face) != 0. && (indicatrice_face(autre_face) != 1.
                                                                 || (indicatrice_face(autre_face) == 1. && d<0.))) //YG: si indicatrice modifiee
                        {
                          double dd = 0. ;
                          for( int j=0 ; j<dim ; j++)
                            {
                              dd += (Vertex(autre_face,j)-xv(i_face,j))*(Vertex(autre_face,j)-xv(i_face,j)) ;
                            }
                          dd = sqrt(dd) ;
                          if( dd <= dmin )
                            {
                              dmin = dd ;
                              for( int k=0 ; k<dim ; k++)
                                {
                                  PPP(i_face,k) = Vertex(autre_face,k) ;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// Calcul du vertex le plus proche des faces fluides
void Transport_Interfaces_FT_Disc::PPP_face_voisin( const DoubleTab& indicatrice, const DoubleTab& indicatrice_face, DoubleTab& PPP )

{
  const int dim = Objet_U::dimension;
  const Domaine_dis_base& mon_dom_dis = domaine_dis().valeur();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, mon_dom_dis);
  const IntTab& faces_elem = domaine_vf.face_voisins();
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const int nfaces = faces_elem.dimension(0) ;
  const DoubleTab& xv = domaine_vf.xv();

  for( int i_face=0 ; i_face<nfaces ; i_face++ )
    {
      // Distance min entre la face et sommet le plus proche
      double dmin = 1.e30 ;

      if (indicatrice_face(i_face) == 0. || indicatrice_face(i_face) == 1.)
        {
          // On parcourt les autres faces de mes deux voisins
          // On compare le PPP de i_face avec ceux des autres faces.
          // Parmi eux, on choisit le plus proche de i_face.
          for( int ii=0; ii<2; ii++ )
            {
              const int elem = faces_elem(i_face,ii) ;
              // Si cet element voisin existe
              if( elem >-1 )
                {
                  // alors on parcourt ses autres faces
                  for (int face_loc=0 ; face_loc<2*dim ; face_loc++)
                    {
                      const int face = elem_faces(elem,face_loc) ;

                      if (indicatrice_face(face) != 0. && indicatrice_face(face) != 1.)
                        {
                          // on identifie les elements auxquels appartient autre_face
                          for(int i=0 ; i<2; i++)
                            {
                              const int elem_voisin = faces_elem(face,i) ;

                              // On veut que l'element existe, qu'il soit lui-meme
                              // traverse par une interface
                              if( elem_voisin>-1 && indicatrice(elem_voisin) != 0. && indicatrice(elem_voisin) != 1. )
                                {
                                  // alors on parcourt ses autres faces
                                  for (int autre_face_loc=0 ; autre_face_loc<2*dim ; autre_face_loc++)
                                    {
                                      const int autre_face = elem_faces(elem_voisin,autre_face_loc) ;
                                      double dd = 0. ;
                                      for( int j=0 ; j<dim ; j++)
                                        {
                                          dd += (PPP(autre_face,j)-xv(i_face,j))*(PPP(autre_face,j)-xv(i_face,j)) ;
                                        }
                                      dd = sqrt(dd) ;
                                      if( dd <= dmin )
                                        {
                                          dmin = dd ;
                                          for( int k=0 ; k<dim ; k++)
                                            {
                                              PPP(i_face,k) = PPP(autre_face,k) ;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// Calcul du nombre max de facettes du maillage Lagrangien par element.
void Transport_Interfaces_FT_Disc::calcul_maxfa7( Maillage_FT_Disc& maillage,
                                                  const DoubleTab& indicatrice,
                                                  const int nb_elem, int& max_fa7,
                                                  const int exec_planfa7existant)
{
  const Intersections_Elem_Facettes& intersection = maillage.intersections_elem_facettes();
  const ArrOfInt& index_elem = intersection.index_elem() ;
  max_fa7 = -1 ;

  // Parcourt des elements reels
  for (int i_elem = 0; i_elem < nb_elem ; i_elem++)
    {
      // Test si l'element est traverse par une IBC
      if( indicatrice(i_elem) != 0. && indicatrice(i_elem) != 1. )
        {
          int index = index_elem[i_elem] ;
          int cpt = 0 ;
          // Liste des parametres du plan de la facette i_facette
          DoubleList A, B, C, D ;

          // Pour un element, on parcourt les facettes qui le traversent
          while (index >= 0 )
            {
              const Intersections_Elem_Facettes_Data& data = intersection.data_intersection(index);
              const int i_facette = data.numero_facette_;

              int test_liste = 0 ; // test_liste vaut 1 si le plan de i_facette existe
              // Calcul des parametres du plan de la facette
              // (sauf pour la premiere car les listes sont vides )
              if( !(A.est_vide()) && exec_planfa7existant == 1 )
                plan_facette_existant( maillage, A, B, C, D, i_facette, test_liste ) ;

              // Test si i_facette est la premiere facette parcourue
              // ou si le plan qui la traverse n'existe pas encore
              if( test_liste==0 )
                {
                  if( exec_planfa7existant == 1 )
                    {
                      double aa=0., bb=0., cc=0., dd=0. ;
                      calcul_eq_plan_facette( i_facette, aa, bb, cc, dd ) ;
                      // Stokage des parametres du plan de la facette.
                      // Listes utilisees dans plan_facette_existant
                      A.add(aa) ;
                      B.add(bb) ;
                      C.add(cc) ;
                      D.add(dd) ;
                    }
                  cpt++ ;
                }

              index = data.index_facette_suivante_;
            }
          // Calcul du nombre max de facettes retenues par element
          if( cpt > max_fa7 )
            max_fa7 = cpt ;
        }
    }
}

// Renumerotation des facettes virtuelles
void Transport_Interfaces_FT_Disc::RenumFa7( DoubleTab& Barycentre, DoubleTab& Tab110,
                                             DoubleTab& Tab111, DoubleTab& Tab112,
                                             IntTab& Tab12, IntTab& CptFacette,
                                             const int nb_facettes,
                                             const int nb_facettes_dim )
{
  //double eps = 1e-20 ;
  const double eps = Objet_U::precision_geom ;

  int nn = Tab110.dimension(0) ;
  int nn_tot = Tab110.dimension_tot(0) ;
  for (int i_elem = nn; i_elem < nn_tot ; i_elem++)
    {

      for( int cpt=0 ; cpt < CptFacette(i_elem) ; cpt++ )
        {
          int bary_trouve = 0 ;
          int ii = 0 ;
          while( ii < nb_facettes && bary_trouve==0 )
            {
              bool test_0 = (Barycentre(ii,0)-Tab110(i_elem,cpt))*(Barycentre(ii,0)-Tab110(i_elem,cpt))<eps ;
              bool test_1 = (Barycentre(ii,1)-Tab111(i_elem,cpt))*(Barycentre(ii,1)-Tab111(i_elem,cpt))<eps ;
              bool test_2 = (Barycentre(ii,2)-Tab112(i_elem,cpt))*(Barycentre(ii,2)-Tab112(i_elem,cpt))<eps ;

              if( test_0 && test_1 && test_2 )
                {
                  bary_trouve = 1 ;
                  Tab12( i_elem, cpt ) = ii ;
                }
              ii++ ;
            }
          if( bary_trouve==0 && ii==nb_facettes )
            {

              // On renumerote localement ssi le numero de la facette virtuelle
              // existe dans la numerotation des facettes reelles.
              if( Tab12( i_elem, cpt ) < nb_facettes )
                Tab12( i_elem, cpt )+=nb_facettes_dim ;
            }
        }
    }
}

void Transport_Interfaces_FT_Disc::StockageFa7( Maillage_FT_Disc& maillage,
                                                IntTab& CptFacette, DoubleTab& Tab100,
                                                DoubleTab& Tab101, DoubleTab& Tab102,
                                                DoubleTab& Tab103, DoubleTab& Tab110,
                                                DoubleTab& Tab111, DoubleTab& Tab112,
                                                IntTab& Tab12, DoubleTab& Barycentre,
                                                const DoubleTab& indicatrice,
                                                IntList& OutElem, ArrOfBit& fa7,
                                                const int exec_planfa7existant )
{
  const Intersections_Elem_Facettes& intersection = maillage.intersections_elem_facettes();
  const ArrOfInt& index_elem = intersection.index_elem() ;
  int nn = Tab100.dimension(0) ;
  int nn_tot = Tab100.dimension_tot(0) ;

  // Matrice de booleen evitant de calculer deux fois le meme barycentre
  fa7=0;
  for (int i_elem = 0; i_elem < nn ; i_elem++)
    {
      if( indicatrice(i_elem) != 0. && indicatrice(i_elem) != 1. && !OutElem.contient(i_elem) )
        {

          int index = index_elem[i_elem] ;
          CptFacette(i_elem) = 0 ;
          int cpt = 0 ;
          // Les listes sont vides pour la premieres facettes de chaque element
          DoubleList A, B, C, D ;
          while (index >= 0 )
            {
              const Intersections_Elem_Facettes_Data& data = intersection.data_intersection(index);
              const int i_facette = data.numero_facette_;


              int test_liste = 0 ; // si test_liste=1 le plan de i_facette existe
              if( !(A.est_vide()) && exec_planfa7existant == 1 )
                plan_facette_existant( maillage, A, B, C, D, i_facette, test_liste ) ;

              // Test si i_facette est la premiere facette parcourue
              // ou si le plan qui la traverse n'existe pas encore
              if( test_liste==0 )
                {
                  // Calcul des parametres du plan de la facette
                  calcul_eq_plan_facette(maillage,i_facette,Tab100(i_elem,cpt),Tab101(i_elem,cpt),
                                         Tab102(i_elem,cpt),Tab103(i_elem,cpt)) ;
                  if( exec_planfa7existant == 1 )
                    {
                      A.add(Tab100(i_elem,cpt)) ;
                      B.add(Tab101(i_elem,cpt)) ;
                      C.add(Tab102(i_elem,cpt)) ;
                      D.add(Tab103(i_elem,cpt)) ;
                    }
                  Tab12(i_elem,cpt) = i_facette ;

                  // Calcul de l'octree de la facette ( seulement en parallele )
                  if( nn < nn_tot )
                    {
                      // Calcul des coordonnees barycentriques de i_facette
                      if( !fa7.testsetbit(i_facette) )
                        {
                          BaryFa7(maillage,i_facette,Barycentre) ;
                        }
                      // Stokage des coordonnees barycentriques de i_facette
                      Tab110(i_elem,cpt) = Barycentre(i_facette,0) ;
                      Tab111(i_elem,cpt) = Barycentre(i_facette,1) ;
                      Tab112(i_elem,cpt) = Barycentre(i_facette,2) ;
                    }
                  cpt++ ;
                }
              index = data.index_facette_suivante_;
            }
          // Nombre de facettes dans l'element
          CptFacette(i_elem) = cpt ;
          if( cpt > Tab12.line_size() )
            {
              Cerr<<"Attention, seulement "<< cpt << " facettes ont ete stokees "<< finl ;
              Cerr<<"dans l'element "<< i_elem << " qui en contient plus."<< finl ;
              Cerr<<"Regarder la finesse du maillage Lagrangien par rapport au maillage Eulerien" <<finl ;
              exit() ;
            }
        }
    }
}

void Transport_Interfaces_FT_Disc::StockageFa7( Maillage_FT_Disc& maillage, DoubleTab& Tab100,
                                                DoubleTab& Tab101, DoubleTab& Tab102, DoubleTab& Tab103,
                                                DoubleTab& Tab110,DoubleTab& Tab111, DoubleTab& Tab112,
                                                IntTab& Tab12, DoubleTab& Barycentre, IntList& OutElem,
                                                IntTab& TabOutFa7, ArrOfBit& fa7 )
{
  int cpt_elem = 0 ;
  int nn = Tab100.dimension(0) ;
  int nn_tot = Tab100.dimension_tot(0) ;
  while( cpt_elem < OutElem.size() )
    {
      for( int k=0 ; k<TabOutFa7.line_size() ; k++ )
        {
          const int i_facette = TabOutFa7(cpt_elem,k) ;
          calcul_eq_plan_facette(maillage,i_facette,Tab100( OutElem[cpt_elem],k),
                                 Tab101(OutElem[cpt_elem],k),Tab102(OutElem[cpt_elem],k),
                                 Tab103(OutElem[cpt_elem],k)) ;
          Tab12(OutElem[cpt_elem],k) = i_facette ;

          // Calcul de l'octree de la facette ( seulement en parallele )
          if( nn < nn_tot )
            {
              // Calcul des coordonnees barycentriques de i_facette
              if( !fa7.testsetbit(i_facette) )
                {
                  BaryFa7(maillage,i_facette,Barycentre) ;
                }
              // Stokage des coordonnees barycentriques de i_facette
              Tab110(OutElem[cpt_elem],k) = Barycentre(i_facette,0) ;
              Tab111(OutElem[cpt_elem],k) = Barycentre(i_facette,1) ;
              Tab112(OutElem[cpt_elem],k) = Barycentre(i_facette,2) ;
            }
        }
      cpt_elem++ ;
    }
}

// Calcul des coordonees barycentriques des facettes
void Transport_Interfaces_FT_Disc::BaryFa7( Maillage_FT_Disc& maillage, const int i_facette, DoubleTab& Bary )
{
  const DoubleTabFT& sommets = maillage.sommets() ;
  const IntTab& facettes = maillage.facettes() ;
  const int dim = Objet_U::dimension ;
  IntTab Som( facettes.line_size() ) ;
  DoubleTab CoordSom( facettes.line_size(), dim ) ;

  // Calcul des coordonnees de chaque sommets
  for ( int i=0 ; i<facettes.line_size() ; i ++ )
    {
      Bary(i_facette,i) = 0. ;
      Som(i) =  facettes(i_facette, i) ;
      for ( int j = 0 ; j < dim ; j++ )
        {
          CoordSom(i,j) = sommets(Som(i),j) ;
        }
    }
  // Calcul des coordonees barycentriques
  for ( int k = 0 ; k< dim ; k++ )
    {
      for( int l = 0 ; l< facettes.line_size() ; l++ )
        {
          Bary(i_facette,k) += CoordSom(l,k) / double(facettes.line_size()) ;
        }
    }
  if( dim == 2 )
    Bary(i_facette,2) = 0. ;
}

// Methode determinant si le plan passant par la facette i_facette existe
void Transport_Interfaces_FT_Disc::plan_facette_existant( Maillage_FT_Disc& maillage,DoubleList A,
                                                          DoubleList B,DoubleList C,DoubleList D,
                                                          const int i_facette,int& test_liste )
{
  assert( A.size() == B.size() ) ;
  assert( A.size() == C.size() ) ;
  assert( A.size() == D.size() ) ;
  const DoubleTabFT& sommets = maillage.sommets() ;
  const IntTab& facettes = maillage.facettes() ;
  double x1, y1, z1 = 0. ;
  double x2, y2, z2 = 0. ;
  double x3 = 0., y3 = 0., z3 = 0. ;
  double plan3 = 0. ;
  const int s1 = facettes(i_facette, 0) ;
  x1 = sommets(s1,0) ;
  y1 = sommets(s1,1) ;
  const int s2 = facettes(i_facette, 1) ;
  x2 = sommets(s2,0) ;
  y2 = sommets(s2,1) ;
  double eps = 1e-20 ;

  // const double eps = Objet_U::precision_geom ;

  if (Objet_U::dimension == 3)
    {
      z1 = sommets(s1,2) ;
      z2 = sommets(s2,2) ;
      const int s3 = facettes(i_facette, 2) ;
      x3 = sommets(s3,0) ;
      y3 = sommets(s3,1) ;
      z3 = sommets(s3,2) ;
    }

  int i=0 ;
  while( i< A.size() && test_liste==0 )
    {
      double plan1 = D[i] + A[i] * x1 + B[i] * y1 + C[i] * z1 ;
      double plan2 = D[i] + A[i] * x2 + B[i] * y2 + C[i] * z2 ;
      if(Objet_U::dimension == 3)
        plan3 = D[i] + A[i] * x3 + B[i] * y3 + C[i] * z3 ;


      if( std::fabs(plan1)<eps && std::fabs(plan2)<eps && std::fabs(plan3)<eps )
        test_liste= 1 ;
      i++ ;
    }
}

// Calcul des parametres du plan de la facette. Utile pour l'algorithme d'Uzawa.
void Transport_Interfaces_FT_Disc::calcul_eq_plan_facette( Maillage_FT_Disc& maillage,const int i_facette,
                                                           double& a,double& b,double& c,double& d )
{
  const DoubleTabFT& normale_facettes = maillage.get_update_normale_facettes() ;
  const DoubleTabFT& sommets = maillage.sommets() ;
  const IntTab& facettes = maillage.facettes() ;

  c = 0. ;
  double x, y, z = 0. ;
  double inverse_norme = 0. ;

  const int s1 = facettes(i_facette, 0) ;

  a = normale_facettes(i_facette, 0) ;
  b = normale_facettes(i_facette, 1) ;
  x = sommets(s1,0) ;
  y = sommets(s1,1) ;

  if (Objet_U::dimension == 3)
    {
      c = normale_facettes(i_facette, 2) ;
      z = sommets(s1,2) ;
    }
  inverse_norme = 1. / sqrt( a*a + b*b + c*c ) ;

  a *= inverse_norme ;
  b *= inverse_norme ;
  c *= inverse_norme ;
  d = - (a * x + b * y + c * z) ;
}

// Calcul des parametres du plan de la facette. Utile pour l'algorithme d'Uzawa.
void Transport_Interfaces_FT_Disc::calcul_eq_plan_facette( const int i_facette,double& a,
                                                           double& b,double& c,double& d )
{
  Maillage_FT_Disc& maillage = maillage_interface() ;
  const DoubleTabFT& normale_facettes = maillage.get_update_normale_facettes() ;
  const DoubleTabFT& sommets = maillage.sommets() ;
  const IntTab& facettes = maillage.facettes() ;

  c = 0. ;
  double x, y, z = 0. ;
  double inverse_norme = 0. ;

  const int s1 = facettes(i_facette, 0) ;

  a = normale_facettes(i_facette, 0) ;
  b = normale_facettes(i_facette, 1) ;
  x = sommets(s1,0) ;
  y = sommets(s1,1) ;

  if (Objet_U::dimension == 3)
    {
      c = normale_facettes(i_facette, 2) ;
      z = sommets(s1,2) ;
    }
  inverse_norme = 1. / sqrt( a*a + b*b + c*c ) ;

  a *= inverse_norme ;
  b *= inverse_norme ;
  c *= inverse_norme ;
  d = - (a * x + b * y + c * z) ;
}

void Transport_Interfaces_FT_Disc::calcul_tolerance_projete_diphasique( const int i_face, const int ori, const int voisin0,
                                                                        const int voisin1, const DoubleTab& indicatrice, double& tol )
{
  const int dim = Objet_U::dimension;
  const Domaine_dis_base& mon_dom_dis = domaine_dis().valeur();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, mon_dom_dis);
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, mon_dom_dis) ;
  const IntTab& elem_faces = domaine_vf.elem_faces();
  int voisin = 0 ;
  DoubleTab L(dim) ;
  L = 0. ;
  tol = 0. ;
  if( voisin0 > -1 && indicatrice(voisin0) != 0. && indicatrice(voisin0) != 1. )
    {
      voisin = voisin0 ;
      const int face_voisine = elem_faces(voisin0,ori) + elem_faces(voisin0,ori+dim) - i_face ;
      L(ori) = std::fabs(domaine_vdf.distance_face(i_face,face_voisine,ori)) ;
    }
  if( voisin1 > -1  && indicatrice(voisin1) != 0. && indicatrice(voisin1) != 1. )
    {
      voisin = voisin1 ;
      const int face_voisine = elem_faces(voisin1,ori) + elem_faces(voisin1,ori+dim) - i_face ;
      double xx = std::fabs(domaine_vdf.distance_face(i_face,face_voisine,ori)) ;
      L(ori) = std::max( L(ori), xx ) ;
    }

  for( int k = 0 ; k<dim ; k++ )
    {
      if( k != ori )
        {
          const int face0 = elem_faces(voisin,k) ;
          const int face1 = elem_faces(voisin,k+dim) ;
          L(k) = std::fabs(domaine_vdf.distance_face(face0,face1,k))/2. ;
        }
    }
  for( int k = 0 ; k<dim ; k++ )
    {
      tol += L(k)*L(k) ;
    }
  tol = sqrt(tol) ;
}

void Transport_Interfaces_FT_Disc::calcul_tolerance_projete_monophasique( const int i_face, const int ori, const int voisin0,
                                                                          const int voisin1, const DoubleTab& indicatrice_face,
                                                                          const DoubleTab& indicatrice, double& tol )
{
  const int dim = Objet_U::dimension ;
  const Domaine_dis_base& mon_dom_dis = domaine_dis().valeur() ;
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, mon_dom_dis) ;
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, mon_dom_dis) ;
  const IntVect& orientation = domaine_vdf.orientation() ;
  const IntTab& elem_faces = domaine_vf.elem_faces() ;
  const IntTab& faces_elem = domaine_vf.face_voisins() ;

  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  // cas longitudinal
  DoubleTab L(dim) ;
  L = 0. ;
  double tol0 = 0. ;
  if( voisin0 > -1 )
    {
      const int face_voisine = elem_faces(voisin0,ori) + elem_faces(voisin0,ori+dim) - i_face ;
      if( indicatrice_face(face_voisine) != 0. && indicatrice_face(face_voisine) != 1. )
        {
          L(ori) = std::fabs(domaine_vdf.distance_face(i_face,face_voisine,ori)) ;

          const int voisin_voisin = faces_elem(face_voisine,0) + faces_elem(face_voisine,1) - voisin0 ;
          if( voisin_voisin > -1 )
            {
              const int face_voisine_voisine = elem_faces(voisin_voisin,ori) + elem_faces(voisin_voisin,ori+dim) - face_voisine ;
              L(ori) += std::fabs(domaine_vdf.distance_face(face_voisine,face_voisine_voisine,ori)) ;
            }
          for( int k = 0 ; k<dim ; k++ )
            {
              if( k != ori )
                {
                  const int face0 = elem_faces(voisin0,k) ;
                  const int face1 = elem_faces(voisin0,k+dim) ;
                  L(k) = std::fabs(domaine_vdf.distance_face(face0,face1,k))/2. ;
                }
            }
          for( int k = 0 ; k<dim ; k++ )
            {
              tol0 += L(k)*L(k) ;
            }
          tol0 = sqrt(tol0) ;
        }
    }
  double tol1 = 0. ;
  L = 0. ;
  if( voisin1 > -1 )
    {
      const int face_voisine = elem_faces(voisin1,ori) + elem_faces(voisin1,ori+dim) -i_face ;
      if( indicatrice_face(face_voisine) != 0. && indicatrice_face(face_voisine) != 1. )
        {
          L(ori) = std::fabs(domaine_vdf.distance_face(i_face,face_voisine,ori)) ;

          const int voisin_voisin = faces_elem(face_voisine,0) + faces_elem(face_voisine,1) - voisin1 ;
          if( voisin_voisin > -1 )
            {
              const int face_voisine_voisine = elem_faces(voisin_voisin,ori) + elem_faces(voisin_voisin,ori+dim) - face_voisine ;
              L(ori) += std::fabs(domaine_vdf.distance_face(face_voisine,face_voisine_voisine,ori)) ;
            }
          for( int k = 0 ; k<dim ; k++ )
            {
              if( k != ori )
                {
                  const int face0 = elem_faces(voisin1,k) ;
                  const int face1 = elem_faces(voisin1,k+dim) ;
                  L(k) = std::fabs(domaine_vdf.distance_face(face0,face1,k))/2. ;
                }
            }
          for( int k = 0 ; k<dim ; k++ )
            {
              tol1 += L(k)*L(k) ;
            }
          tol1 = sqrt(tol1) ;
        }
    }
  double tol_long = std::max(tol0,tol1) ;


  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  // cas transverse
  L = 0. ;
  int voisin = 1 ;
  if( voisin0 > -1 )
    {
      voisin = voisin0 ;
      const int face_voisine = elem_faces(voisin0,ori) + elem_faces(voisin0,ori+dim) - i_face ;
      L(ori) = std::fabs(domaine_vdf.distance_face(i_face,face_voisine,ori)) ;
    }
  if( voisin1 > -1 )
    {
      voisin = voisin1 ;
      const int face_voisine = elem_faces(voisin1,ori) + elem_faces(voisin1,ori+dim) - i_face ;
      double xx = std::fabs(domaine_vdf.distance_face(i_face,face_voisine,ori)) ;
      L(ori) = std::max( L(ori), xx ) ;
    }

  DoubleTab xx_max(dim) ;
  xx_max = 0. ;
  for (int i=0 ; i<2*dim ; i++)
    {
      const int autre_face = elem_faces(voisin,i) ;
      const int face_voisine = elem_faces(voisin,ori) + elem_faces(voisin,ori+dim) - i_face ;

      if( autre_face != i_face && autre_face != face_voisine
          && indicatrice_face(autre_face) != 0. && indicatrice_face(autre_face) != 1. )
        {
          const int voisin_voisin = faces_elem(autre_face,0) +  faces_elem(autre_face,1) - voisin ;
          const int ori_face = orientation[autre_face] ;
          const int autre_face_voisine = elem_faces(voisin_voisin,ori_face) + elem_faces(voisin_voisin,ori_face+dim) - autre_face ;
          double xx = std::fabs(domaine_vdf.distance_face(autre_face,autre_face_voisine,ori_face)) ;
          if( xx_max(ori_face) < xx )
            xx_max(ori_face) = xx ;
        }
    }

  DoubleTab TOL_TMP(dim-1) ;
  TOL_TMP = 0. ;
  int cpt = 0 ;
  for(int k=0 ; k<dim ; k++)
    {
      if( k != ori )
        {
          for(int l=0 ; l<dim ; l++)
            {
              if( l != k && l != ori )
                {
                  const int face0 = elem_faces(voisin,k) ;
                  const int face1 = elem_faces(voisin,k+dim) ;
                  double yy  = 0.5*std::fabs(domaine_vdf.distance_face(face0,face1,k)) ;
                  double a = 0.5*yy + xx_max(k) ;
                  const int face3 = elem_faces(voisin,l) ;
                  const int face4 = elem_faces(voisin,l+dim) ;
                  double b = 0.5*std::fabs(domaine_vdf.distance_face(face3,face4,l)) ;
                  TOL_TMP(cpt) = L(ori)*L(ori) + a*a + b*b ;
                  cpt++ ;
                }
            }
        }
    }
  double tol_trans = TOL_TMP(0) ;
  if( dim == 3 )
    tol_trans = std::max(tol_trans,TOL_TMP(1)) ;
  tol_trans = sqrt( tol_trans ) ;
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  tol = std::max(tol_long,tol_trans) ;
}

// Verification du projete
void Transport_Interfaces_FT_Disc::verifprojete( const int monophasique,const double Lref, double d, const DoubleTab& x,
                                                 const DoubleTab& V, DoubleTab& coord_projete, int& cpt )
{
  const int dim = coord_projete.dimension(0) ;
  const double dist_IBC = std::fabs(d) ;
  const double precision = Objet_U::precision_geom ;
  double dist_proj = 0. ;
  int projete_modifie = 0 ;
  cpt = 0 ;

  DoubleTab coord_projete_old(dim) ;
  coord_projete_old = coord_projete ;

  for (int j = 0; j < dim ; j++)
    {
      dist_proj += (x(j) - coord_projete(j))*(x(j) - coord_projete(j)) ;
    }
  dist_proj = sqrt( dist_proj ) ;
  double dist_proj_old = dist_proj ;

  const double facsec = 0.5 ;
  const int test0 = ( dist_proj < precision ) ;
  const int test1 = ( std::fabs( dist_IBC - dist_proj ) > facsec*dist_IBC ) ;
  const int test2 = ( dist_proj > (1.+facsec)*Lref ) ;

  if( monophasique )
    {
      if( test0 || (test1 && test2) )
        {
          // on modifie le projete en prenant le vertex lagrangien le plus proche
          for (int j = 0; j < dim ; j++)
            {
              coord_projete(j) = V(j) ;
            }
          projete_modifie = 1 ;
          cpt++ ;
        }
    }
  else
    {
      const int test = ( dist_IBC >= (1.-facsec)*Lref ) ;

      if( (test0 && test) || (test1 && test2) )
        {
          // on modifie le projete en prenant le vertex lagrangien le plus proche
          for (int j = 0; j < dim ; j++)
            {
              coord_projete(j) = V(j) ;
            }
          projete_modifie = 1 ;
          cpt++ ;
        }
    }

  if( projete_modifie == 1 )
    {
      dist_proj = 0. ;
      for (int j = 0; j < dim ; j++)
        {
          dist_proj += (x(j) - coord_projete(j))*(x(j) - coord_projete(j)) ;
        }
      dist_proj = sqrt( dist_proj ) ;

      Cerr << " ------------------------------------------------------- "<<finl;
      Cerr << " le projete "<<finl;
      for(int i=0; i<dim; i++)
        {
          Cerr <<" X("<<i<<") = "<<coord_projete_old(i)<<finl ;
        }
      Cerr << " situe a une distance "<<dist_proj_old<<" du barycentre de la face "<<finl ;
      for(int i=0; i<dim; i++)
        {
          Cerr <<" Y("<<i<<") = "<<x(i)<<finl ;
        }
      Cerr << " a ete modifie en "<<finl ;
      for(int i=0; i<dim; i++)
        {
          Cerr <<" X("<<i<<") = "<<coord_projete(i)<<finl ;
        }
      Cerr << "situe a une distance "<<dist_proj<<" du barycentre de face "<<finl ;
      Cerr << "distance a l'IBC : "<<dist_IBC<<finl ;
      Cerr << "Longueur de reference : "<<Lref<<finl ;
      Cerr << " ------------------------------------------------------- "<<finl;
    }
}

// Algorithme d'Uzawa calculant les coordonnees du projete du centre de gravite
// d'une face 'i_face' sur l'interface.
void Transport_Interfaces_FT_Disc::uzawa(const double d, const DoubleTab& C,
                                         const DoubleTab& x, const DoubleTab& f,
                                         DoubleTab& x_ibc) const
{
  double rho = 0. ;
  const int nb_lignes = C.dimension(0) ;
  const int nb_colonnes = C.line_size() ;
  assert(x_ibc.dimension(0) == nb_colonnes) ;
  assert(f.dimension(0) == nb_lignes) ;
  assert(x.dimension(0) == nb_colonnes) ;

  // Calcul du pas rho :
  // Definition et initialisation des normes induites ||.||_1 et||.||_inf
  double norme1 = 0. ;
  double norme_infty = 0. ;

  for (int j=0; j<nb_colonnes; j++)
    {
      double norme1_int = 0. ;
      for (int i=0; i<nb_lignes; i++)
        {
          norme1_int += std::fabs( C(i,j) ) ;
        }
      norme1 = std::max( norme1 , norme1_int ) ;
    }

  for (int i=0; i<nb_lignes; i++)
    {
      double norme_infty_int = 0.;
      for (int j=0; j<nb_colonnes; j++)
        {
          norme_infty_int += std::fabs( C(i,j) ) ;
        }
      norme_infty = std::max( norme_infty , norme_infty_int ) ;
    }
  rho = 1. / ( norme1 * norme_infty ) ; //car ||A||_2^2 <= ||A||_1 ||A||_infty

  // Definition et initialisation d'un multiplicateur de Lagrange
  DoubleTab lambda(nb_lignes) ;
  lambda = 0. ;

  double distance = 1. ;
  double distance_old = 0. ;
  double deplacement_relatif = 1.;
  int nb_iter = 0 ;

  double err_uzawa = variables_internes_->seuil_uzawa ;
  int nb_iter_max = variables_internes_->nb_iter_uzawa ;
  // On arrete l'algorithme des que le carre de la distance (projete,origine)
  // a l'etape k moins le carre de la distance (projete,origine) a l'etape k-1
  // est inferieur en valeur absolue au seuil 'err_uzawa'.
  // Pour des raisons de temps CPU, l'algorithme s'arrete egalement si
  // le nombre d'iterations est superieur a 'nb_iter_max'.
  // En pratique, les valeurs par defaut sont 'err_uzawa=1e-8' et 'nb_iter_max=30'.
  // Dans certains cas, par exemple pour des maillages tres grossiers, il a ete
  // observe que ce choix conduit a quelques differences sequentiel/parallele.
  // On recommande alors de diminuer le seuil et d'augmenter le nombre
  // d'iteration pour ameliorer la convergence de l'algorithme, au risque d'augmenter
  // le temps CPU.
  x_ibc = 0. ;
  while( deplacement_relatif > err_uzawa  && nb_iter <= nb_iter_max )
    {
      // Calcul de la solution a l'etape k
      x_ibc = 0. ;
      distance_old = distance ;
      for (int i=0; i<x_ibc.dimension(0); i++)
        {
          for (int j=0; j<lambda.dimension(0); j++)
            {
              x_ibc(i) -= 0.5 * C(j,i) * lambda(j) ;
            }
          x_ibc(i) += x[i] ;
        }

      // Calcul du multiplicateur de lagrange a l'etape k+1
      for (int i=0; i<lambda.dimension(0); i++)
        {
          //Parametre xx = C*x_ibc-f
          double xx = 0. ;

          for (int j=0; j<x_ibc.dimension(0); j++)
            {
              xx += C(i,j)*x_ibc(j) ;
            }
          xx -= f(i) ;
          if ( d > 0 )
            {
              // cote solide
              lambda(i) = std::max( lambda(i) + rho * xx , 0. ) ;
            }
          else
            {
              //cote fluide
              lambda(i) = std::min( lambda(i) + rho * xx , 0. ) ;
            }
        }
      // Calcul et mise a jour de la condition du while
      distance = ( x_ibc(0) - x[0] ) * ( x_ibc(0) - x[0] )
                 + ( x_ibc(1) - x[1] ) * ( x_ibc(1) - x[1] ) ;
      if (Objet_U::dimension==3)
        distance += ( x_ibc(2) - x[2] ) * ( x_ibc(2) - x[2] ) ;

      deplacement_relatif = std::fabs(distance - distance_old) ;
      nb_iter++ ;

    }

}

void Transport_Interfaces_FT_Disc::projete_point_face_fluide( int& nb_proj_modif,const int dimfa7,
                                                              const DoubleTab& indicatrice_face, const DoubleTab& indicatrice,
                                                              const DoubleTab& dist_face,const double t,const double dt,
                                                              DoubleTab& Tab100,DoubleTab& Tab101,DoubleTab& Tab102,DoubleTab& Tab103,
                                                              IntTab& Tab12,IntTab& CptFacette,DoubleTab& v_imp,DoubleTab& Vertex,Parser& parser_x, Parser& parser_y, Parser& parser_z  )
{
  const int dim = Objet_U::dimension;
  const Domaine_dis_base& mon_dom_dis = domaine_dis().valeur();
  const Domaine_VDF * zvdf = 0;
  if (sub_type(Domaine_VDF, mon_dom_dis))
    zvdf = &ref_cast(Domaine_VDF, mon_dom_dis) ;
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, mon_dom_dis);
  const IntVect& orientation = domaine_vdf.orientation();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, mon_dom_dis);
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const IntTab& faces_elem = domaine_vf.face_voisins();
  const int nfaces = faces_elem.dimension(0) ;
  const DoubleTab& xv = domaine_vf.xv();
  const double invalid_test = -1.e30;

  Cerr << "CALCUL DU PROJETE POUR LES FACES FLUIDE ET d > invalid_test" << finl ;

  for (int i_face = 0; i_face < nfaces; i_face++)
    {
      double d = dist_face(i_face) ;
      if (indicatrice_face(i_face) == 0. && d > invalid_test)
        {
          DoubleList mat_0, mat_1, mat_2, mat_3 ;
          ArrOfBit fa7(dimfa7);
          fa7=0;
          int nb_facettes_to_uzawa = 0 ;
          // On travaille sur les elements voisins de chaque face 'i_face'
          for( int ii=0; ii<2; ii++ )
            {
              const int elem = faces_elem(i_face,ii) ;
              // Si cet element voisin existe
              if( elem >-1 )
                {
                  // alors on parcourt ses autres faces
                  for (int face_loc=0 ; face_loc<2*dim ; face_loc++)
                    {
                      const int autre_face = elem_faces(elem,face_loc) ;

                      if (indicatrice_face(autre_face) != 0. && indicatrice_face(autre_face) != 1.)
                        {
                          // on identifie les elements auxquels appartient autre_face
                          // et qui sont traverses par l'interface
                          for(int ii2=0 ; ii2<2; ii2++)
                            {
                              const int elem_voisin = faces_elem(autre_face,ii2) ;

                              // On veut que l'element existe, qu'il soit lui-meme
                              // traverse par une interface
                              if( elem_voisin>-1 && indicatrice(elem_voisin) != 0.
                                  && indicatrice(elem_voisin) != 1. )
                                {

                                  for( int i=0 ; i< CptFacette(elem_voisin) ; i++)
                                    {
                                      if( !fa7.testsetbit( Tab12(elem_voisin,i) ) )
                                        {
                                          nb_facettes_to_uzawa++ ;
                                          mat_0.add( Tab100(elem_voisin,i) ) ;
                                          mat_1.add( Tab101(elem_voisin,i) ) ;
                                          mat_2.add( Tab102(elem_voisin,i) ) ;
                                          mat_3.add( Tab103(elem_voisin,i) ) ;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
          if( nb_facettes_to_uzawa > 0 )
            {
              DoubleTab matrice_plans(nb_facettes_to_uzawa, dim) ;
              DoubleTab x(dim) ;
              DoubleTab secmem(nb_facettes_to_uzawa) ;

              for( int i_facette=0 ; i_facette<nb_facettes_to_uzawa ; i_facette++ )
                {
                  matrice_plans(i_facette,0) = mat_0[i_facette] ;
                  matrice_plans(i_facette,1) = mat_1[i_facette] ;
                  if (dim==3)
                    matrice_plans(i_facette,2) = mat_2[i_facette] ;
                  secmem[i_facette] = -mat_3[i_facette] ;
                }

              for (int i_coord=0; i_coord<dim; i_coord++)
                {
                  x(i_coord) = xv(i_face,i_coord) ;
                }

              // Calcul des coordonnees du projete 'coord_projete'
              DoubleTab coord_projete(dim) ;
              coord_projete = 0. ;
              uzawa(d,matrice_plans,x,secmem,coord_projete) ;

              if (variables_internes_->type_projete_calcule == Transport_Interfaces_FT_Disc_interne::PROJETE_MODIFIE)
                {
                  DoubleTab V(dim) ;
                  for (int i_coord=0; i_coord<dim; i_coord++)
                    {
                      V(i_coord) = Vertex(i_face,i_coord) ;
                    }
                  int cpt = 0 ;
                  double Lref =0. ;
                  const int ori = orientation[i_face] ;
                  const int voisin0 = faces_elem(i_face,0) ;
                  const int voisin1 = faces_elem(i_face,1) ;
                  const int monophasique = 1 ;
                  calcul_tolerance_projete_monophasique(i_face,ori,voisin0,voisin1,indicatrice_face,indicatrice,Lref) ;
                  verifprojete(monophasique,Lref,d,x,V,coord_projete,cpt) ;
                  if( cpt > 0 )
                    nb_proj_modif++ ;
                }

              //Calcul de la vitesse du projete
              if (variables_internes_->methode_transport == Transport_Interfaces_FT_Disc_interne::VITESSE_IMPOSEE)
                {
                  for (int j = 0; j < dim; j++)
                    {
                      parser_x.setVar(j, coord_projete(j)) ;
                      parser_y.setVar(j, coord_projete(j)) ;
                      if (dim==3)
                        parser_z.setVar(j, coord_projete(j)) ;
                    }
                  if (zvdf)
                    {
                      switch(zvdf->orientation(i_face))
                        {
                        case 0:
                          v_imp(i_face) = parser_x.eval() ;
                          break ;
                        case 1:
                          v_imp(i_face) = parser_y.eval() ;
                          break ;
                        case 2:
                          v_imp(i_face) = parser_z.eval() ;
                          break ;
                        }
                    }
                }
              else if (variables_internes_->methode_transport == Transport_Interfaces_FT_Disc_interne::LOI_HORAIRE )
                {
                  ArrOfDouble v(dim) ;
                  v = variables_internes_->loi_horaire_->vitesse(t+dt,coord_projete);
                  if (zvdf)
                    v_imp(i_face)=v[zvdf->orientation(i_face)];
                }
            }
        }
    }
}

void Transport_Interfaces_FT_Disc::projete_point_face_interface( int& nb_proj_modif, const int dimfa7, const DoubleTab& indicatrice_face,
                                                                 const DoubleTab& indicatrice, const DoubleTab& dist_face,
                                                                 const double t, const double dt, DoubleTab& Tab100,
                                                                 DoubleTab& Tab101, DoubleTab& Tab102,DoubleTab& Tab103,
                                                                 IntTab& Tab12,IntTab& CptFacette, DoubleTab& v_imp,
                                                                 DoubleTab& Vertex,Parser& parser_x, Parser& parser_y,
                                                                 Parser& parser_z )
{
  const int dim = Objet_U::dimension;
  const Domaine_dis_base& mon_dom_dis = domaine_dis().valeur();
  const Domaine_VDF * zvdf = 0 ;
  if (sub_type(Domaine_VDF, mon_dom_dis))
    zvdf = &ref_cast(Domaine_VDF, mon_dom_dis);
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, mon_dom_dis);
  const IntVect& orientation = domaine_vdf.orientation();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, mon_dom_dis);
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const IntTab& faces_elem = domaine_vf.face_voisins();
  const int nfaces = faces_elem.dimension(0) ;
  const DoubleTab& xv = domaine_vf.xv() ;
  const double invalid_test = -1.e30;

  for( int i_face=0 ; i_face<nfaces ; i_face++ )
    {
      double d = dist_face(i_face) ;
      const int diphasique_tmp = (indicatrice_face(i_face) != 0. &&
                                  (indicatrice_face(i_face) != 1. || (indicatrice_face(i_face) == 1. && d<=0. && d> invalid_test)));
      const int diphasique = ( variables_internes_->is_extra_diphasique_solide ? diphasique_tmp : (diphasique_tmp && d <= 0.) ) ;
      const int fluide = ( indicatrice_face(i_face) == 0. ) ;
      const int solide = ( variables_internes_->is_extra_solide && (d>0. || (indicatrice_face(i_face) == 1. && d<-1.e20)) ) ;
      const int monophasique = ( fluide || solide ) ;

      if (diphasique || (monophasique && d> invalid_test ))
        {
          DoubleList mat_0, mat_1, mat_2, mat_3 ;
          ArrOfBit fa7(dimfa7);
          fa7=0;
          int nb_facettes_to_uzawa = 0 ;

          // Coordonnes du projet 'coord_projete'
          DoubleTab coord_projete(dim) ;
          coord_projete = 0. ;
          int projete_calcule = 0 ;

          // Si la face a une indicatrice de 0.5 celui implique
          // qu'elle est parfaitement traverse par au moins une facette.
          // Dans ce cas, le projet du barycentre de la face sur le
          // plan des facettes est lui-mme.
          if( indicatrice_face(i_face) == 0.5 )
            {
              int voisin = -1 ;
              if( faces_elem(i_face,0) > -1 )
                voisin = faces_elem(i_face,0) ;
              else
                voisin = faces_elem(i_face,1) ;

              if( indicatrice(voisin) == 0. || indicatrice(voisin) == 1. )
                {
                  for( int j=0 ; j<dim ; j++)
                    {
                      coord_projete(j)=xv(i_face,j) ;
                    }
                  projete_calcule = 1 ;
                }
            }

          if( projete_calcule == 0 )
            {
              if( d == 0. )
                {
                  for( int j=0 ; j<dim ; j++)
                    {
                      coord_projete(j)=xv(i_face,j) ;
                    }
                }
              else
                {
                  // on identifie les elements auxquels appartient i_face
                  // et qui sont traverses par l'interface
                  for(int ii=0 ; ii<2; ii++)
                    {
                      const int elem_voisin = faces_elem(i_face,ii) ;


                      // On veut que l'element existe, qu'il soit lui-meme
                      // traversepar une interface
                      // YG : ce test a ete scinde en deux. Dans le cas d une indicatrice
                      // modifiee, il peut arriver qu aucun des voisins ait une indicatrice
                      // element !=0. Pour autant, on a envie de regarder les voisins des
                      // voisins pour chercher des facettes. Idem pour le cas monophasique!
                      if( elem_voisin>-1)
                        {
                          if (indicatrice(elem_voisin) != 0. && indicatrice(elem_voisin) != 1. )
                            {
                              for( int i=0 ; i< CptFacette(elem_voisin) ; i++)
                                {
                                  if( !fa7.testsetbit( Tab12(elem_voisin,i) ) )
                                    {
                                      nb_facettes_to_uzawa++ ;
                                      mat_0.add( Tab100(elem_voisin,i) ) ;
                                      mat_1.add( Tab101(elem_voisin,i) ) ;
                                      mat_2.add( Tab102(elem_voisin,i) ) ;
                                      mat_3.add( Tab103(elem_voisin,i) ) ;
                                    }
                                }
                            }
                          // On identifie ensuite les elements voisins (elem_voisin_voisin)
                          // traverses par l'IBC.
                          for(int j=0 ; j<2*dim ; j++)
                            {
                              const int elem_voisin_voisin = faces_elem(elem_faces(elem_voisin,j),0)
                                                             + faces_elem(elem_faces(elem_voisin,j),1) - elem_voisin ;

                              if( elem_voisin_voisin>-1 && indicatrice(elem_voisin_voisin) != 0.
                                  && indicatrice(elem_voisin_voisin) != 1.
                                  && elem_voisin_voisin!=faces_elem(i_face,0)
                                  && elem_voisin_voisin!=faces_elem(i_face,1) )
                                {
                                  for( int i=0 ; i< CptFacette(elem_voisin_voisin) ; i++)
                                    {
                                      if( !fa7.testsetbit( Tab12(elem_voisin_voisin,i) ) )
                                        {
                                          nb_facettes_to_uzawa++ ;
                                          mat_0.add( Tab100(elem_voisin_voisin,i) ) ;
                                          mat_1.add( Tab101(elem_voisin_voisin,i) ) ;
                                          mat_2.add( Tab102(elem_voisin_voisin,i) ) ;
                                          mat_3.add( Tab103(elem_voisin_voisin,i) ) ;
                                        }
                                    }
                                }
                            }
                        }
                    }

                  if( nb_facettes_to_uzawa == 0 && diphasique )
                    {
                      Cerr << "erreur concernant la face : "<< i_face << " I_f = "<<indicatrice_face(i_face)<<finl ;
                      Cerr << "aucune facette detecte "<< finl ;
                      exit() ;
                    }
                  else if (nb_facettes_to_uzawa >0)
                    {

                      DoubleTab matrice_plans(nb_facettes_to_uzawa, dim) ;
                      DoubleTab x(dim) ;
                      DoubleTab secmem(nb_facettes_to_uzawa) ;

                      for( int i_facette=0 ; i_facette<nb_facettes_to_uzawa ; i_facette++ )
                        {
                          matrice_plans(i_facette,0) = mat_0[i_facette] ;
                          matrice_plans(i_facette,1) = mat_1[i_facette] ;
                          if (dim==3)
                            matrice_plans(i_facette,2) = mat_2[i_facette] ;
                          secmem[i_facette] = -mat_3[i_facette] ;
                        }
                      for (int i_coord=0; i_coord<dim; i_coord++)
                        {
                          x(i_coord) = xv(i_face,i_coord) ;
                        }
                      coord_projete = 0. ;
                      uzawa(d,matrice_plans,x,secmem,coord_projete) ;

                      if(variables_internes_-> type_projete_calcule == Transport_Interfaces_FT_Disc_interne::PROJETE_MODIFIE)
                        {
                          DoubleTab V(dim) ;
                          for (int i_coord=0; i_coord<dim; i_coord++)
                            {
                              V(i_coord) = Vertex(i_face,i_coord) ;
                            }
                          int cpt = 0 ;
                          double Lref = 0. ;
                          const int ori = orientation[i_face] ;
                          const int voisin0 = faces_elem(i_face,0) ;
                          const int voisin1 = faces_elem(i_face,1) ;
                          const int monophasique2 = 0 ;
                          calcul_tolerance_projete_diphasique(i_face,ori,voisin0,voisin1,indicatrice,Lref) ;
                          verifprojete(monophasique2,Lref,d,x,V,coord_projete,cpt) ;
                          if( cpt > 0 )
                            nb_proj_modif++ ;
                        }
                    }
                }
            }

          const int test = (diphasique || fluide || solide ) ;


          if(test)
            {
              // Calcul de la vitesse du projete
              if (variables_internes_->methode_transport == Transport_Interfaces_FT_Disc_interne::VITESSE_IMPOSEE)
                {
                  for (int j = 0; j < dim; j++)
                    {
                      parser_x.setVar(j, coord_projete(j)) ;
                      parser_y.setVar(j, coord_projete(j)) ;
                      if (dim==3)
                        parser_z.setVar(j, coord_projete(j)) ;
                    }
                  if (zvdf)
                    {
                      switch(zvdf->orientation(i_face))
                        {
                        case 0:
                          v_imp(i_face) = parser_x.eval() ;
                          break ;
                        case 1:
                          v_imp(i_face) = parser_y.eval() ;
                          break ;
                        case 2:
                          v_imp(i_face) = parser_z.eval() ;
                          break ;
                        }
                    }
                }
              else if (variables_internes_->methode_transport == Transport_Interfaces_FT_Disc_interne::LOI_HORAIRE )
                {
                  ArrOfDouble v(dim);
                  v = variables_internes_->loi_horaire_->vitesse(t+dt,coord_projete);
                  if (zvdf)
                    v_imp(i_face)=v[zvdf->orientation(i_face)];
                }
            }
        }
    }
}

/*! @brief Deplace le maillage a l'aide du champ de vitesse impose entre l'instant maillage.
 *
 * temps() et l'instant "temps".
 *   La nouvelle position des sommets est obtenue par integration des lignes de courant de ce champ,
 *   par un schema RK3 (en un pas de temps) qui prend en compte la dependance en temps des vitesses imposees.
 *   On "nettoie" et on change le temps du maillage. Pas de remaillage. maillage MINIMAL en sortie.
 *
 */
static void deplacer_maillage_ft_v_impose(Noms expression_vitesse,
                                          Maillage_FT_Disc& maillage, double temps)
{
  const double dt = temps - maillage.temps();
  Parser parser_x, parser_y, parser_z;
  init_parser_v_impose(expression_vitesse,
                       parser_x, parser_y, parser_z, temps);

  const DoubleTab& sommets = maillage.sommets();
  const int dim = sommets.dimension(1);
  const int dimension3 = (dim == 3);
  int i;
  const int nb_sommets = maillage.nb_sommets();
  DoubleTabFT deplacement(nb_sommets, dim);
  for (i = 0; i < nb_sommets; i++)
    {
      const double x0 = sommets(i, 0);
      const double y0 = sommets(i, 1);
      const double z0 = dimension3 ? sommets(i, 2) : 0.;
      double x = x0, y = y0, z = z0;
      integrer_vitesse_imposee(parser_x, parser_y, parser_z,
                               temps, dt, x, y, z);
      deplacement(i, 0) = x - x0;
      deplacement(i, 1) = y - y0;
      if (dimension3)
        deplacement(i, 2) = z - z0;
    }
  maillage.transporter(deplacement);
  maillage.nettoyer_elements_virtuels();
  maillage.changer_temps(temps);
}

// Calcul d'une vitesse de deplacement minimisant les deplacements relatifs
// des sommets d'interface dans le repere local (donc minimisant les remaillages locaux)
// On veut que, dans le repere mobile associe au centre de gravite de l'interface,
//  le deplacement tangentiel soit nul.
void Transport_Interfaces_FT_Disc::calculer_vitesse_repere_local(const Maillage_FT_Disc& maillage,  DoubleTab& deplacement,DoubleTab& Positions,DoubleTab& Vitesses) const
{
  // const Maillage_FT_Disc & maillage = maillage_interface();
  //vitesse des centres d'aire
  const int nb_facettes=maillage.nb_facettes();
  const IntTab& facettes=maillage.facettes ();
  const int dim3 = (deplacement.line_size() == 3);
  ArrOfInt compo_connexe_facettes(nb_facettes); // Init a zero
  int n = search_connex_components_local_FT(maillage, compo_connexe_facettes);
  int nb_compo_tot=compute_global_connex_components_FT(maillage, compo_connexe_facettes, n);

  DoubleTabFT normale;
  calculer_normale_sommets_interface(maillage, normale);
  Positions.resize(nb_compo_tot,dimension);
  Vitesses.resize(nb_compo_tot,dimension);

  calculer_vmoy_composantes_connexes(maillage, compo_connexe_facettes, nb_compo_tot,
                                     deplacement, Vitesses, Positions);

  // Calcul de la composante connexe des sommets
  // Attention un sommet peut n'etre rattache a aucune facette sur le meme processeur !
  // (il faudrait calculer les compo connexes sur les sommets, et ensuite passer aux faces
  //  ce serait plus simple, voire calculer les deux en meme temps)
  IntVect compo_sommets;
  maillage.creer_tableau_sommets(compo_sommets, RESIZE_OPTIONS::NOCOPY_NOINIT);
  compo_sommets = -1;
  {
    const int dim = deplacement.dimension(1);
    for (int iface = 0; iface < nb_facettes; iface++)
      {
        const int compo = compo_connexe_facettes[iface];
        for (int j = 0; j < dim; j++)
          compo_sommets[facettes(iface, j)] = compo;
      }
    // On prend le max sur tous les processeurs qui partagent le sommet pour les sommets isoles
    // (max calcule uniquement pour les sommets reels, sommets virtuels faux)
    MD_Vector_tools::echange_espace_virtuel(compo_sommets, MD_Vector_tools::EV_MAX);
    // Inutile de synchroniser, on utilise uniquement les sommets reels
  }

  const int nb_sommets_tot = compo_sommets.size_totale();
  for (int som = 0; som < nb_sommets_tot; som++)
    {
      if (maillage.sommet_virtuel(som))
        {
          // Valeur pipo pour dire qu'on initialise pas
          deplacement(som, 0) = 100.;
          deplacement(som, 1) = 100.;
          continue;
        }

      const int compo = compo_sommets[som];

      // (v-vmoy) doit etre normal a l'interface
      // Donc on fait v_corrige = v_initial - composante_tangentielle_de(v_initial-vmoy)
      // Demonstration que (v_corrige - vmoy) est normal a l'interface :
      // On note ct() = composante_tangentielle_de()
      // ct(v_corrige - vmoy)
      //  = ct(v_corrige)                            - ct(vmoy)    car ct est lineaire
      //  = ct(v_initial - ct(v_initial) + ct(vmoy)) - ct(vmoy)    car ct est lineaire
      //  = ct(v_initial) - ct(v_initial) + ct(vmoy) - ct(vmoy)    cat ct(ct(x)) = ct(x) et linearite
      //  = 0

      // v_corrige = cn(v_initial - vmoy) + vmoy
      // v_corrige = cn(v_initial) - cn(vmoy) + vmoy
      // v_corrige = cn(v_initial) + ct(vmoy)
      // v_corrige = v_initial - ct(v_inital) + ct(vmoy)
      // v_corrige = v_initial - ct(v_initial - vmoy)

      double nx = 0., ny = 0., nz = 0.;
      double vi_x = 0., vi_y = 0., vi_z = 0.;

      nx = normale(som, 0);
      ny = normale(som, 1);

      vi_x = deplacement(som, 0);
      vi_y = deplacement(som, 1);
      if (dim3)
        {
          vi_z = deplacement(som, 2);
          nz=normale(som, 2);
        }

      double prodscal =
        (vi_x - Vitesses(compo, 0)) * nx
        + (vi_y - Vitesses(compo, 1)) * ny;
      if (dim3)
        prodscal += (vi_z - Vitesses(compo, 2)) * nz;
      double norme_carre = nx * nx + ny * ny + nz * nz;
      if (norme_carre != 0.)
        {
          prodscal /= norme_carre;
          deplacement(som, 0) = nx * prodscal + Vitesses(compo, 0);
          deplacement(som, 1) = ny * prodscal + Vitesses(compo, 1);
          if (dim3)
            deplacement(som, 2) = nz * prodscal + Vitesses(compo, 2); // BugFix reported from baltik TCL on 2020/10/26
        }
    }
}

void Transport_Interfaces_FT_Disc::deplacer_maillage_ft_v_fluide(const double temps)
{
  DoubleTab& deplacement = variables_internes_->deplacement_sommets;
  const Equation_base& eqn_hydraulique = variables_internes_->refequation_vitesse_transport.valeur();
  const Champ_base& champ_vitesse = eqn_hydraulique.inconnue().valeur();
  Maillage_FT_Disc& maillage = maillage_interface();
  // Calcul de la vitesse de deplacement des sommets par interpolation
  // (deplacement contient en fait la vitesse en m/s)
  int flag = 1;
  if (sub_type(Navier_Stokes_FT_Disc, eqn_hydraulique))
    {
      const Navier_Stokes_FT_Disc& ns = ref_cast(Navier_Stokes_FT_Disc, eqn_hydraulique);
      if (ns.get_delta_vitesse_interface())
        flag=0;
    }
  calculer_vitesse_transport_interpolee(champ_vitesse,
                                        maillage,
                                        deplacement, 1 /* recalculer le champ de vitesse L2 */,
                                        flag /* Interpolation Multi-lineaire en VDF */);

#if DEBUG_CONSERV_VOLUME
  int n = maillage.nb_sommets();
  double dmin=0.;
  double dmax=0.;
  double dmean = 0.;
  DoubleTab norm_deplacement(n);
  if (n>0)
    {
      for (int i=0; i<n; i++)
        {
          double x = 0.;
          for (int j=0; j<deplacement.dimension(1); j++)
            x += deplacement(i,j)*deplacement(i,j);
          norm_deplacement[i]  = sqrt(x);
          dmean +=sqrt(x);
        }
      dmean /=n;
      dmin = min_array(norm_deplacement);
      dmax = max_array(norm_deplacement);
    }
#endif
  // On recupere et ajoute a deplacement le saut de vitesse a l'interface si changement de phase
  ajouter_contribution_saut_vitesse(deplacement);

#if DEBUG_CONSERV_VOLUME
  if (n>0)
    {
      double dmean_wpch = 0.;
      for (int i=0; i<n; i++)
        {
          double x = 0.;
          for (int j=0; j<deplacement.dimension(1); j++)
            x += deplacement(i,j)*deplacement(i,j);
          norm_deplacement[i]  = sqrt(x);
          dmean_wpch +=sqrt(x);
        }
      dmean_wpch /=n;
      double dmin_wpch = min_array(norm_deplacement);
      double dmax_wpch = max_array(norm_deplacement);
      Cerr << "Interfacial_velocity before/after pch [min/mean/max]. Time: "
           << schema_temps().temps_courant()
           << " Without-pch: " << dmin << " " << dmean<< " "  << dmax
           << " With-pch: "  << dmin_wpch << " " << dmean_wpch << " "  << dmax_wpch
           << " [Warning, values invalid in //] " << finl;
    }
#endif

  if (interpolation_repere_local_)
    {
      DoubleTab Positions,Vitesses;

      calculer_vitesse_repere_local(maillage, deplacement,Positions,Vitesses);
      assert(Positions.dimension(0)==Vitesses.dimension(0));
      if(Process::je_suis_maitre())
        {
          ofstream fout;
          fout.open("composantes_connexes.txt",ios::app);
          fout << "TEMPS: " << temps << std::endl;
          const int nb_bulles=Positions.dimension(0);

          for(int bulle=0; bulle<nb_bulles; bulle++)
            {
              fout << "composante " << bulle << " position " ;
              for (int i=0; i <Positions.dimension(1); i++) fout << " " << Positions(bulle,i);
              fout << " vitesse " ;
              for (int i=0; i <Vitesses.dimension(1); i++) fout << " " << Vitesses(bulle,i);
              fout << std::endl;
            }
          fout.close();
        }
    }

  const double delta_t = temps - maillage.temps();
  // Calcul du deplacement :
  deplacement *= delta_t;

// Debug GB 2020.03.20 Conservation de volume
#if DEBUG_CONSERV_VOLUME
  double  volume_avt = remaillage_interface().calculer_volume_mesh(maillage);
#endif
  if (variables_internes_->VOFlike_correction_volume > 0)
    {
      // Transport avec correction du volume des phases :
      // Sauvegarde de la position actuelle des sommets :
      DoubleTabFT position_precedente = maillage.sommets();
      maillage.preparer_tableau_avant_transport(position_precedente,
                                                maillage.desc_sommets());
      // Calcul de la variation de volume exacte a obtenir a partir du
      // champ de vitesse eulerien :
      ArrOfDoubleFT var_volume;
      {
        // Desole pour le ref_cast_non_const, il y a probablement plus propre mais je ne sais pas faire.
        // En l'etat, 'ns' doit etre non-const car il fait appel a get_update...normale et get_compute_indicatrice_faces
        // qui modifient l'objet Transport_... dans lequel nous sommes!
        Navier_Stokes_FT_Disc& ns = ref_cast_non_const(Navier_Stokes_FT_Disc, eqn_hydraulique);

        DoubleVect dI_dt;
        domaine_dis().valeur().domaine().creer_tableau_elements(dI_dt);
        ns.calculer_dI_dt(dI_dt);
        dI_dt.echange_espace_virtuel();
#if DEBUG_CONSERV_VOLUME
        const int nb_elem = domaine_dis().valeur().nb_elem();
        double sum_before_rm = 0.;
        double sum_before_rm_dvol = 0.;
        for (int i = 0; i < nb_elem; i++)
          sum_before_rm +=-dI_dt[i]; // It's already homogeneous to *volumes[i];

        sum_before_rm_dvol +=sum_before_rm*delta_t;
#endif
        ramasse_miettes(maillage, variables_internes_->tmp_flux->valeurs(), dI_dt);
        transfert_conservatif_eulerien_vers_lagrangien_sommets(maillage, dI_dt, var_volume);

#if DEBUG_CONSERV_VOLUME
        double sum = 0.;
        for (int i = 0; i < nb_elem; i++)
          sum +=-dI_dt[i]; // It's already homogeneous to *volumes[i];

        const double dvoldt_totale = remaillage_interface().calculer_somme_dvolume(maillage, var_volume);
        const double sum_dvol =sum*delta_t;
        Cerr << " time " << temps << " sum_dI_dt " << sum  << " sum_dvol " << sum_dvol
             << " sum_before_rm_dI_dt " << sum_before_rm  << " sum_before_rm_dvol " << sum_before_rm_dvol
             << " V_lagrangien= " << dvoldt_totale <<finl;
#endif
      }
      // var_volume est une derivee par rapport au temps.
      // calcul de l'integrale pendant le pas de temps ...
      // et changement de signe car on veut la variation de volume de la phase 0
      // (et non celle de la phase 1)
      var_volume *= -delta_t;
// Debug GB 2019.02.08 Conservation de volume
#if DEBUG_CONSERV_VOLUME
      double dvol_theo_depl = remaillage_interface().calculer_somme_dvolume(maillage, var_volume);
#endif
      maillage.preparer_tableau_avant_transport(var_volume,
                                                maillage.desc_sommets());
      // Transport avec le deplacement interpole :
      remaillage_interface().traite_decollement(maillage_interface(), deplacement);
      maillage.transporter(deplacement);
      maillage.update_tableau_apres_transport(position_precedente,
                                              maillage.nb_sommets(),
                                              maillage.desc_sommets());
      maillage.update_tableau_apres_transport(var_volume,
                                              maillage.nb_sommets(),
                                              maillage.desc_sommets());

      // Calcul de la variation de volume obtenue par ce deplacement :
      ArrOfDoubleFT var_volume_deplacement;
      remaillage_interface().calculer_variation_volume(maillage,
                                                       position_precedente,
                                                       var_volume_deplacement);
#if DEBUG_CONSERV_VOLUME
      double  volume_apres = remaillage_interface().calculer_volume_mesh(maillage);
      double dvol_reel_depl = remaillage_interface().calculer_somme_dvolume(maillage, var_volume_deplacement);
      Cerr << "Transport_Interfaces_FT_Disc::calculer_vitesse_repere_local " << finl
           << " volume avt= " << volume_avt << finl
           << " apres= " << volume_apres << finl
           << " dvol_theo_depl= " << dvol_theo_depl << finl
           << " dvol_reel_depl= " << dvol_reel_depl << finl
           << finl;
#endif
      // Calcul de la variation de volume de la phase 0 a imposer lors de la correction
      // de volume :
      var_volume -= var_volume_deplacement;

      // Si volume de phase_1 imposee : on calcule une deuxieme correction
      if (variables_internes_->volume_impose_phase_1 > 0.)
        {
          DoubleVect values(3);
          values=0.;
//        volume_phase_1     ->   values(0)
//        volume_sous_domaine   ->   values(1)
//        volume_phase_0     ->   values(2)
          if (variables_internes_->nom_domaine_volume_impose_ == "??")
            values(0) = calculer_integrale_indicatrice(indicatrice_.valeurs(), values(2));
          else
            {
              const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
              const DoubleVect& volumes = domaine_vf.volumes();
              const Sous_Domaine& sous_domaine = domaine_dis().valeur().domaine().ss_domaine(variables_internes_->nom_domaine_volume_impose_);
              const int nb_elem_sous_domaine = sous_domaine.nb_elem_tot();
              const DoubleTab& indic = indicatrice_.valeurs();
              const int nb_elem = domaine_vf.nb_elem();
              for (int i = 0; i < nb_elem_sous_domaine; i++)
                {
                  const int elem = sous_domaine[i];
                  if (elem < nb_elem)
                    {
                      values(0) += indic(elem) * volumes(elem);
                      values(1) += volumes(elem);
                      values(2) += (1.-indic(elem)) * volumes(elem);
                    }
                }
              mp_sum_for_each_item(values);
              /*              Cerr << "Volume_sous_domaine " << variables_internes_->nom_domaine_volume_impose_ <<
                                 mp_sum(volume_sous_domaine) << finl;
                            volume_phase_1 = mp_sum(volume_phase_1); */
              Cerr << "Volume_sous_domaine " << variables_internes_->nom_domaine_volume_impose_ <<
                   values(1) << finl;

            }
          const double erreur_volume = variables_internes_->volume_impose_phase_1 - values(0);
          Journal() << "Transport_Interfaces_FT_Disc::mettre_a_jour "
                    << "correction_volume_impose_phase_1= " << erreur_volume << finl;
          // Calcul de la correction de volume a imposer a chaque sommet
          // Surface totale d'interface:
          double surface_tot = 0.;
          {
            int i;
            const int nb_facettes = maillage.nb_facettes();
            const ArrOfDouble& surfaces = maillage.get_update_surface_facettes();

            for (i = 0; i < nb_facettes; i++)
              surface_tot += surfaces[i];

            surface_tot = mp_sum(surface_tot) * 3.;
            if (surface_tot > 0.)
              {
                const double inv_surface_tot = 1. / surface_tot;
                const IntTab& facettes = maillage.facettes();
                const int nb_som_facette = facettes.dimension(1);
                for (i = 0; i < nb_facettes; i++)
                  {
                    int j;
                    const double dvolume = erreur_volume * surfaces[i] * inv_surface_tot;
                    for (j = 0; j < nb_som_facette; j++)
                      {
                        int isom = facettes(i,j);
                        // var_volume est la variation de volume de phase 0, donc signe "moins"
                        var_volume[isom] -= dvolume;
                      }
                  }
              }
          }
        }
      // Correction de volume :
      // GB. 2019.05.13. Pourquoi utiliser iterations_correction_volume pour le nombre de lissage de dvol?
      // remaillage_interface().lisser_dvolume(maillage, var_volume,
      //                                      variables_internes_->iterations_correction_volume);
      // 2019.05.13. Nouveau code :
      remaillage_interface().lisser_dvolume(maillage, var_volume,
                                            variables_internes_->nb_lissage_correction_volume);

#if DEBUG_CONSERV_VOLUME
      double  volume_before = remaillage_interface().calculer_volume_mesh(maillage);
      double dvol_total_before = remaillage_interface().calculer_somme_dvolume(maillage, var_volume);
      Cerr << " dvol_error_before= " << dvol_total_before << " Volume_before= " << volume_before << " time " << temps << finl;
#endif
      // GB. 2019.05.13. Pourquoi utiliser nb_iter_bary_volume_seul_ plutot que iterations_correction_volume?
      //remaillage_interface().corriger_volume(maillage, var_volume);
      // 2019.05.13. Nouveau code :
      remaillage_interface().corriger_volume_(maillage, var_volume,
                                              variables_internes_->nb_iterations_correction_volume);
#if DEBUG_CONSERV_VOLUME
      double  volume_after = remaillage_interface().calculer_volume_mesh(maillage);
      double dvol_total_after = remaillage_interface().calculer_somme_dvolume(maillage, var_volume);
      Cerr << " dvol_error_after= " << dvol_total_after << " Volume_after= " << volume_after << " time " << temps << finl;
#endif
    }
  else
    {
#if DEBUG_CONSERV_VOLUME
      // Sauvegarde de la position actuelle des sommets :
      DoubleTabFT position_precedente = maillage.sommets();
      maillage.preparer_tableau_avant_transport(position_precedente,
                                                maillage.desc_sommets());
#endif
      // Transport par interpolation de vitesse seule :
      remaillage_interface().traite_decollement(maillage_interface(), deplacement);
      maillage.transporter(deplacement);
#if DEBUG_CONSERV_VOLUME
      // Calcul de la variation de volume obtenue par ce deplacement :
      ArrOfDoubleFT var_volume_deplacement;
      remaillage_interface().calculer_variation_volume(maillage,
                                                       position_precedente,
                                                       var_volume_deplacement);

      maillage.update_tableau_apres_transport(position_precedente,
                                              maillage.nb_sommets(),
                                              maillage.desc_sommets());
      double  volume_apres = remaillage_interface().calculer_volume_mesh(maillage);
      double dvol_reel_depl = remaillage_interface().calculer_somme_dvolume(maillage, var_volume_deplacement);
      Cerr << "Transport_Interfaces_FT_Disc::calculer_vitesse_repere_local " << finl
           << " volume avt= " << volume_avt << finl
           << " apres= " << volume_apres << finl;
// pas calcule          Cerr << " dvol_theo_depl= " << dvol_theo_depl << finl;
      Cerr     << " dvol_reel_depl= " << dvol_reel_depl << finl
               << " dt= " << delta_t << finl;
#endif
    }
  remaillage_interface().traite_adherence(maillage_interface());
  maillage.changer_temps(temps);
}

void Transport_Interfaces_FT_Disc::ajouter_contribution_saut_vitesse(DoubleTab& deplacement) const
{
  const Equation_base& eqn_hydraulique = variables_internes_->refequation_vitesse_transport.valeur();
  const Champ_base *u0_ptr = 0;

  if (sub_type(Navier_Stokes_FT_Disc, eqn_hydraulique))
    {
      // On recupere le saut de vitesse a l'interface (changement de phase)
      const Navier_Stokes_FT_Disc& ns = ref_cast(Navier_Stokes_FT_Disc, eqn_hydraulique);
      u0_ptr = ns.get_delta_vitesse_interface();
      if (u0_ptr)
        {
          const Champ_base& u0 = *u0_ptr;
          DoubleTabFT d2(deplacement);
          // If ns.get_new_mass_source() == 0, we use a standard multi-linear interpolation.
          // If ns.get_new_mass_source() == 1, we use the new method (1D-interpolation of each velocity component in its direction)
          calculer_vitesse_transport_interpolee(u0,
                                                maillage_interface(),
                                                d2, 1 /* recalculer le champ de vitesse L2 */,
                                                1-ns.get_new_mass_source());

          const int n = d2.dimension(0);
          const int dim = d2.line_size();
          for (int i = 0; i < n; i++)
            {
              for (int j = 0; j < dim; j++)
                {
                  const double depl2 = d2(i,j);
                  deplacement(i,j) -= depl2;
                }
            }
        }
    }

}

int Transport_Interfaces_FT_Disc::calculer_composantes_connexes_pour_suppression(IntVect& num_compo)
{
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
  const DoubleTab& indicatrice = get_update_indicatrice().valeurs();
  const int nb_compo = topologie_interface().calculer_composantes_connexes_pour_suppression(domaine_vf, indicatrice, num_compo);
  return nb_compo;
}

double Transport_Interfaces_FT_Disc::suppression_interfaces(const IntVect& num_compo, const ArrOfInt& flags_compo_a_supprimer)
{
  Maillage_FT_Disc& maillage = maillage_interface();
  const double volume = topologie_interface().suppression_interfaces(num_compo,
                                                                     flags_compo_a_supprimer, maillage,
                                                                     variables_internes_->indicatrice_cache.valeur().valeurs());
  return volume;
}

void Transport_Interfaces_FT_Disc::test_suppression_interfaces_sous_domaine()
{
  if (suppression_interfaces_sous_domaine_ == "??")
    return;

  const DoubleTab& indicatrice = get_update_indicatrice().valeurs();
  const Sous_Domaine& sous_domaine = domaine_dis().domaine().ss_domaine(suppression_interfaces_sous_domaine_);
  // Construction de la liste des elements de la sous-domaine contenant la phase a supprimer
  ArrOfInt liste_elems_sous_domaine;
  int i;
  const double phase_continue = topologie_interface().get_phase_continue();

  const int nb_elems_sous_domaine = sous_domaine.nb_elem_tot();
  for (i = 0; i < nb_elems_sous_domaine; i++)
    {
      const int i_elem = sous_domaine[i];
      const double indic = indicatrice[i_elem];
      if (indic != phase_continue)
        liste_elems_sous_domaine.append_array(i_elem);
    }
  const int sz_liste_elems_sous_domaine = liste_elems_sous_domaine.size_array();
  if (mp_sum(sz_liste_elems_sous_domaine) > 0)
    {
      IntVect num_compo;
      const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
      const int nb_compo = topologie_interface().calculer_composantes_connexes_pour_suppression(domaine_vf, indicatrice, num_compo);
      // Tableau initialise a zero:
      ArrOfInt flags_compo_a_supprimer(nb_compo);
      for (i = 0; i < sz_liste_elems_sous_domaine; i++)
        {
          const int elem = liste_elems_sous_domaine[i];
          const int compo = num_compo[elem];
          flags_compo_a_supprimer[compo] = 1;
        }
      mp_sum_for_each_item(flags_compo_a_supprimer);
      Maillage_FT_Disc& maillage = maillage_interface();
      topologie_interface().suppression_interfaces(num_compo, flags_compo_a_supprimer, maillage,
                                                   variables_internes_->indicatrice_cache.valeur().valeurs());

      // Parcours de toutes les equations du probleme,
      // Pour les equations "temperature FT" on appelle la methode "suppression_interfaces"
      Probleme_base& pb = probleme();
      const int n = pb.nombre_d_equations();
      for (int ii = 0; ii < n; ii++)
        {
          Equation_base& eq = pb.equation(ii);
          if (sub_type(Convection_Diffusion_Temperature_FT_Disc, eq))
            {
              Convection_Diffusion_Temperature_FT_Disc& eq_temp =
                ref_cast(Convection_Diffusion_Temperature_FT_Disc, eq);
              const int i_phase_continue = (phase_continue < 0.5) ? 0 : 1;
              eq_temp.suppression_interfaces(num_compo, flags_compo_a_supprimer,
                                             i_phase_continue /* phase qui remplace l'ancienne */);
            }
        }
    }
}

//Integration des trajectoires de particules
// La methode est une version simplifiee de deplacer_maillage_ft_v_fluide(...)
// -Interpolation du champ de vitesse de l equation de qdm
// -Calcul du deplacement a appliquer aux particules
// -Transport des particules en fonction du deplacement calcule
// -Mise a jour du temps de l objet ensemble de points meme si t<t_debut_integr
// Rq : t_debut_integr est fixe par defaut a t_init mais peut etre specifie par l utilisateur
void Transport_Interfaces_FT_Disc::integrer_ensemble_lagrange(const double temps)
{
  Cerr<<"The method Transport_Interfaces_FT_Disc::integrer_ensemble_lagrange"<<finl;
  Cerr<<"must be overloaded."<<finl;
  exit();
}

void Transport_Interfaces_FT_Disc::mettre_a_jour(double temps)
{
  Process::Journal() << "Transport_Interfaces_FT_Disc::mettre_a_jour " << le_nom() << " temps= "<< temps << finl;
  Maillage_FT_Disc& maillage = maillage_interface();

  switch (variables_internes_->methode_transport)
    {
    case Transport_Interfaces_FT_Disc_interne::VITESSE_IMPOSEE:
      {
        deplacer_maillage_ft_v_impose(variables_internes_->expression_vitesse_imposee, maillage, temps);
        break;
      }
    case Transport_Interfaces_FT_Disc_interne::VITESSE_INTERPOLEE:
      {
        deplacer_maillage_ft_v_fluide(temps);
        break;
      }
    case Transport_Interfaces_FT_Disc_interne::LOI_HORAIRE:
      {
        // Le deplacement de l'interface est donne par l'expression de la loi horaire
        const double dt = temps - maillage.temps();
        // Le temps auquel est evalue la vitesse:
        //  pour obtenir un deplacement precis a l'ordre 2, on evalue la vitesse au temps suivant:
        const double t = maillage.temps() + dt * 0.5;
        if (Process::je_suis_maitre())
          {
            Cerr << "Transport_Interfaces_FT_Disc::mettre_a_jour imposed motion by a time law at t:" << t << finl;
          }
        const int nb_sommets = maillage.nb_sommets();
        DoubleTabFT deplacement(nb_sommets, Objet_U::dimension);

        const DoubleTab& sommets = maillage.sommets();

        const int dim = Objet_U::dimension;
        ArrOfDouble coord(dim);
        ArrOfDouble coord_barycentre(dim); // Coordonnes du barycentre des points de l'interface
        int nb_sommets_reels=0;
        for (int i = 0; i < nb_sommets; i++)
          {
            // Determination si le sommet de l'interface est reel ou virtuel
            int est_un_sommet_reel=0;
            if (!maillage.sommet_virtuel(i))
              {
                nb_sommets_reels++;
                est_un_sommet_reel=1;
              }
            // Coordonnees coord du ieme sommet de l'interface a l'instant t
            for (int j = 0; j < dim; j++)
              {
                double pos = sommets(i, j);
                coord[j] = pos;
                if (est_un_sommet_reel)
                  coord_barycentre[j] += pos;
              }
            // Calcul de ses nouvelles coordonnees a t+dt
            coord = variables_internes_->loi_horaire_->position(t+dt,t,coord);

            // Calcul du tableau deplacement
            for (int j = 0; j < dim; j++)
              deplacement(i, j) = coord[j] - sommets(i,j);
          }
        // Calcul du barycentre des noeuds de l'interface
        // en tenant compte du parallelisme
        nb_sommets_reels=mp_sum(nb_sommets_reels);
        /*        for (int j = 0; j < dim; j++)
                  {
                    coord_barycentre(j)=mp_sum(coord_barycentre(j))/nb_sommets_reels;
                  } */
        mp_sum_for_each_item(coord_barycentre);
        for (int j = 0; j < dim; j++)
          {
            coord_barycentre[j]=coord_barycentre[j]/nb_sommets_reels;
          }
        // Impression eventuelle de la loi horaire et du centre de gravite "discret" de l'interface
        if (schema_temps().limpr())
          variables_internes_->loi_horaire_->imprimer(schema_temps(),coord_barycentre);

        maillage.transporter(deplacement);
        maillage.nettoyer_elements_virtuels();
        maillage.completer_maillage();
        maillage.changer_temps(temps);
        break;
      }
    default:
      Cerr << "Error for the method Transport_Interfaces_FT_Disc::mettre_a_jour :\n"
           << " the requested transport mehod is not developped" << finl;
      exit();
    }


  // injection des interfaces
  if (variables_internes_->injection_interfaces_temps_.size_array() > 0)
    {
      const ArrOfDouble& tps = variables_internes_->injection_interfaces_temps_;
      const int n = tps.size_array();
      const Noms& expr = variables_internes_->injection_interfaces_expressions_;
      assert(expr.size() == n);
      // Recherche du prochain temps d'injection:
      const double last_time = variables_internes_->injection_interfaces_last_time_;
      int i;
      for (i = 0; i < n; i++)
        if (tps[i] > last_time)
          break;
      for (; i < n && tps[i] <= temps; i++)
        {
          variables_internes_->injection_interfaces_last_time_=tps[i];
          // On essaye d'injecter l'interface
          Maillage_FT_Disc maillage_tmp;
          maillage_tmp.associer_equation_transport(*this);
          Maillage_FT_Disc::AjoutPhase phase = variables_internes_->injection_interfaces_phase_[i]
                                               ? Maillage_FT_Disc::AJOUTE_PHASE1 : Maillage_FT_Disc::AJOUTE_PHASE0;
          DoubleTab sauvegarde(variables_internes_->indicatrice_cache.valeur().valeurs());
          const int ok = marching_cubes().construire_iso(expr[i],
                                                         0., maillage_tmp,
                                                         variables_internes_->indicatrice_cache.valeur().valeurs(),
                                                         phase,
                                                         variables_internes_->distance_interface_sommets);



          Cerr << "Injection_interface time " << temps << " " << expr[i];
          if (ok)
            {
              maillage_interface().ajouter_maillage(maillage_tmp);
              get_update_indicatrice();
              double unused_vol_phase_0 = 0.;
              const double volume_phase_1_old = calculer_integrale_indicatrice(sauvegarde, unused_vol_phase_0);
              unused_vol_phase_0= 0.;
              const double volume_phase_1 = calculer_integrale_indicatrice(variables_internes_->indicatrice_cache.valeur().valeurs(), unused_vol_phase_0);
              double volume = volume_phase_1-volume_phase_1_old;
              // pow(-1,1-phase) ne compile pas avec xlC sur AIX car n'a que pow(double,int)
              volume*=pow(-1.,1-phase);
              Cerr << " volume " << volume << finl;
            }
          else
            {
              Cerr << " failure: collision" << finl;
              variables_internes_->indicatrice_cache.valeur().valeurs() = sauvegarde;
            }
        }
      if (i == n)
        {
          Cerr << "The injection list is fulfilled !"
               << " (if it is wished, add an interface at a large time in the file)"
               << finl;
          barrier();
          exit();
        }
    }

  const Probleme_base&  pb = get_probleme_base();

  if (sub_type(Probleme_FT_Disc_gen,pb))
    {

      Probleme_FT_Disc_gen& pb_ft = ref_cast_non_const(Probleme_FT_Disc_gen, pb);
      // injection des interfaces with temperature of activation
      if (pb_ft.tcl().reinjection_tcl() && pb_ft.tcl().ready_inject_tcl())
        {
          const double thetac = pb_ft.tcl().thetaC_tcl();
          const double Rc = pb_ft.tcl().Rc_inject();

          // const Nom expr = "x^2+(y-8e-05*cos(50.0*pi/180.0))^2-8e-05^2";
          const Nom expr = Nom("x^2+(y-") + Nom(Rc)
                           + Nom("*cos(") + Nom(thetac)
                           + Nom("*pi/180.0))^2-") + Nom(Rc)
                           + Nom("^2");

          // On essaye d'injecter l'interface
          Maillage_FT_Disc maillage_tmp;
          maillage_tmp.associer_equation_transport (*this);

          // By default, inject vapeur, phase 0
          Maillage_FT_Disc::AjoutPhase phase =
            0 ?
            Maillage_FT_Disc::AJOUTE_PHASE1 : Maillage_FT_Disc::AJOUTE_PHASE0;

          DoubleTab sauvegarde (
            variables_internes_->indicatrice_cache.valeur ().valeurs ());

          const int ok = marching_cubes ().construire_iso (
                           expr, 0., maillage_tmp,
                           variables_internes_->indicatrice_cache.valeur ().valeurs (), phase,
                           variables_internes_->distance_interface_sommets);

          Cerr << "Injection_interface time " << temps << " " << expr;
          if (ok)
            {
              maillage_interface ().ajouter_maillage (maillage_tmp);
              get_update_indicatrice ();
              double unused_vol_phase_0 = 0.;
              const double volume_phase_1_old = calculer_integrale_indicatrice (
                                                  sauvegarde, unused_vol_phase_0);
              unused_vol_phase_0 = 0.;
              const double volume_phase_1 = calculer_integrale_indicatrice (
                                              variables_internes_->indicatrice_cache.valeur ().valeurs (),
                                              unused_vol_phase_0);
              double volume = volume_phase_1 - volume_phase_1_old;
              // pow(-1,1-phase) ne compile pas avec xlC sur AIX car n'a que pow(double,int)
              volume *= pow (-1., 1 - phase);
              Cerr << " volume " << volume << finl;
            }
          else
            {
              Cerr << " failure: collision" << finl;
              variables_internes_->indicatrice_cache.valeur ().valeurs () =
                sauvegarde;
            }
        }
    }


  // Traitement des domaines de suppression
  test_suppression_interfaces_sous_domaine();

  // Remaillage de l'interface:
  remailler_interface();

  // Attention: get_update_indicatrice renvoie une ref a indicatrice_cache.
  //  C'est ici qu'on copie le contenu de indicatrice_cache dans indicatrice :
  indicatrice_.valeurs() = get_update_indicatrice().valeurs();

  variables_internes_->indicatrice_cache.changer_temps(temps);
  indicatrice_.changer_temps(temps);

  update_critere_statio();

  get_update_distance_interface();
  get_update_normale_interface();

  {
    double volume_phase_0 = 0.;
    const double volume_phase_1 = calculer_integrale_indicatrice(indicatrice_.valeurs(), volume_phase_0);
    if (Process::je_suis_maitre())
      {
        Cerr << "Volume_phase_0 " << Nom(volume_phase_0, "%20.14g") << " time " << temps << finl;
        Cerr << "Volume_phase_1 " << Nom(volume_phase_1, "%20.14g") << " time " << temps << finl;
      }
  }

  // Affichage de la surface totale d'interfaces dans le fichier .err
  // Affichage du centre de gravite des phases 0 et 1
  {
    // Calcul de la somme des surfaces des facettes reelles:
    //const Maillage_FT_Disc& maillage = maillage_interface();
    const ArrOfDouble& surfaces = maillage.get_update_surface_facettes();
    const int nb_facettes = maillage.nb_facettes();
    double s = 0.;
    int i;
    for (i = 0; i < nb_facettes; i++)
      {
        const int virt = maillage.facette_virtuelle(i);
        if (! virt)
          s += surfaces[i];
      }
    const double s_tot = mp_sum(s);
    if (Process::je_suis_maitre())
      {
        Cerr << "Surface_Totale_Interface " << le_nom()
             << " t= " << temps << " surface= " << s_tot << finl;
      }
    // Calcul du centre de gravite des phases 0 et 1 a partir de l'indicatrice
    // indicatrice de phase:
    const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
    const DoubleTab& indic = indicatrice_.valeurs();
    // centre de gravite des elements euleriens:
    const DoubleTab& xp = domaine_vf.xp();
    // volumes des elements euleriens:
    const DoubleVect& volumes = domaine_vf.volumes();

    const int nb_elem = domaine_vf.nb_elem();
    const int dim = xp.line_size();
    DoubleTrav values(3,dim);
    values=0.;
    // somme des volumes des phases:
//  v0_tot -> values(2,0)
//  v1_tot -> values(2,1)
    // centre de gravite de la phase 0:
//  g0[j] -> values(0,j)
    // phase 1:
//  g1[j] -> values(1,j)
    assert(indic.dimension(0) == nb_elem); // indicatrice aux elements
    // Boucle sur les elements reels
    for (i = 0; i < nb_elem; i++)
      {
        const double ind = indic(i);
        const double v = volumes(i);
        const double v1 = ind * v;
        const double v0 = (1. - ind) * v;
        values(2,0) += v0;
        values(2,1) += v1;
        int j;
        for (j = 0; j < dim; j++)
          {
            const double x = xp(i,j);
            values(0,j) += v0 * x;
            values(1,j) += v1 * x;
          }
      }
    // Somme sur tous les processeurs et calcul du centre de gravite
    double mp_g0[3] = {0., 0., 0.};
    double mp_g1[3] = {0., 0., 0.};

    mp_sum_for_each_item(values);
    {
      int j;
      /*      mp_vtot0 = mp_sum(v0_tot);
            mp_vtot1 = mp_sum(v1_tot);
            if (mp_vtot0 > 0.)
              for (j = 0; j < dim; j++)
                mp_g0[j] = mp_sum(g0[j]) / mp_vtot0;
            if (mp_vtot1 > 0.)
              for (j = 0; j < dim; j++)
                mp_g1[j] = mp_sum(g1[j]) / mp_vtot1; */
      if (values(2,0) > 0.)
        {
          {
            for ( j = 0; j<dim; j++)
              mp_g0[j] = values(0,j) / values(2,0);
          }
          if (values(2,1) > 0.)
            for ( j = 0; j<dim; j++)
              mp_g1[j] = values(1,j) / values(2,1);
        }
    }
    // Affichage
    if (je_suis_maitre())
      {
        Cerr << "Centre_gravite_phases " << le_nom() << " t= " << temps;
        Cerr << " phase0: " << mp_g0[0] << " " << mp_g0[1] << " " << mp_g0[2] ;
        Cerr << " phase1: " << mp_g1[0] << " " << mp_g1[1] << " " << mp_g1[2] << finl;
      }
  }
  //GB : Calcul de la surface d'interfaces par element
  {
    const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
    const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
    const ArrOfInt& index_elem = intersections.index_elem();
    DoubleTab& surface = variables_internes_->surface_interface.valeur().valeurs();
    const int nb_elements = surface.dimension(0);
    for (int element = 0; element < nb_elements; element++)
      {
        int index = index_elem[element];
        double surface_totale = 0.;
        // Boucle sur les faces qui traversent l'element:
        while (index >= 0)
          {
            const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
            surface_totale += data.fraction_surface_intersection_ * surface_facettes[data.numero_facette_];
            index = data.index_facette_suivante_;
          }
        surface[element] = surface_totale;
      }
    surface.echange_espace_virtuel();
    variables_internes_->surface_interface.mettre_a_jour(temps);
  }
  // Fin de GB

  //TF : Gestion de l avancee en temps de la derivee
  if (calculate_time_derivative()) derivee_en_temps().changer_temps(temps);
  //Fin de TF
}

// Deplace les sommets de l'interface du deplacement prescrit (vitesse * coeff)
// ou coeff est homogene a un temps en secondes.
// On met a jour le tableau vitesse en fonction des sommets qui sont
// passes d'un processeur a l'autre.
void Transport_Interfaces_FT_Disc::transporter_sans_changement_topologie(DoubleTab& vitesse, const double coeff,const double temps)
{
  Maillage_FT_Disc& maillage = maillage_interface();
  DoubleTabFT deplacement(vitesse);
  deplacement *= coeff;

  maillage.preparer_tableau_avant_transport(vitesse,
                                            maillage.desc_sommets());

  maillage.transporter(deplacement);

  // attention, creation de sommets virtuels... on fait update_...
  // seulement a la fin.
  maillage.parcourir_maillage();

  maillage.update_tableau_apres_transport(vitesse,
                                          maillage.nb_sommets(),
                                          maillage.desc_sommets());
  //ajout pour postraiter ss pas de tps RK3_FT
  maillage.changer_temps(temps);
  indicatrice_.valeurs()=get_update_indicatrice().valeurs();
  variables_internes_->indicatrice_cache.changer_temps(temps);
  indicatrice_.changer_temps(temps);
  get_update_distance_interface();
  get_update_normale_interface();
}


// ************************************************************************
// ************************************************************************
// Un tas de routines vides ou presque qui ne font rien ou pas grand chose
// ************************************************************************
// ************************************************************************

void Transport_Interfaces_FT_Disc::associer_milieu_base(const Milieu_base& un_milieu)
{
  // Si on veut creer un probleme contenant uniquement une equation Transport,
  // il faut un milieu sinon ca plante dans pb_base::discretiser !
  ref_milieu_ = un_milieu;
}

void Transport_Interfaces_FT_Disc::associer_equation_ns(const Navier_Stokes_FT_Disc& ns)
{
  equation_ns_=ns;
}

Milieu_base& Transport_Interfaces_FT_Disc::milieu()
{
  if (!ref_milieu_.non_nul())
    {
      Cerr << "You forgot to associate the diphasic fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return ref_milieu_.valeur();
}

const Milieu_base& Transport_Interfaces_FT_Disc::milieu() const
{
  if (!ref_milieu_.non_nul())
    {
      Cerr << "You forgot to associate the diphasic fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return ref_milieu_.valeur();
}

int Transport_Interfaces_FT_Disc::nombre_d_operateurs(void) const
{
  return 0;
}

const Operateur& Transport_Interfaces_FT_Disc::operateur(int i) const
{
  assert(0);
  exit();
  throw;
}

Operateur& Transport_Interfaces_FT_Disc::operateur(int i)
{
  assert(0);
  exit();
  throw;
}

const Champ_Inc& Transport_Interfaces_FT_Disc::inconnue(void) const
{
  return indicatrice_;
}

Champ_Inc& Transport_Interfaces_FT_Disc::inconnue(void)
{
  return indicatrice_;
}

Marching_Cubes& Transport_Interfaces_FT_Disc::marching_cubes()
{
  return variables_internes_->marching_cubes_;
}

const Marching_Cubes& Transport_Interfaces_FT_Disc::marching_cubes() const
{
  return variables_internes_->marching_cubes_;
}

Maillage_FT_Disc& Transport_Interfaces_FT_Disc::maillage_interface()
{
  return variables_internes_->maillage_interface;
}

const Maillage_FT_Disc& Transport_Interfaces_FT_Disc::maillage_interface() const
{
  return variables_internes_->maillage_interface;
}
Remaillage_FT& Transport_Interfaces_FT_Disc::remaillage_interface()
{
  return variables_internes_->remaillage_interface_;
}

const Remaillage_FT& Transport_Interfaces_FT_Disc::remaillage_interface() const
{
  return variables_internes_->remaillage_interface_;
}
Topologie_Maillage_FT& Transport_Interfaces_FT_Disc::topologie_interface()
{
  return variables_internes_->topologie_interface_;
}

const Topologie_Maillage_FT& Transport_Interfaces_FT_Disc::topologie_interface() const
{
  return variables_internes_->topologie_interface_;
}

DoubleTab& Transport_Interfaces_FT_Disc::tableaux_positions()
{
  return variables_internes_->doubletab_pos;
}

IntVect& Transport_Interfaces_FT_Disc::vecteur_elements()
{
  return variables_internes_->intvect_elements;
}

DoubleTab& Transport_Interfaces_FT_Disc::deplacement_som()
{
  return variables_internes_->deplacement_sommets;
}

Proprietes_part_vol& Transport_Interfaces_FT_Disc::proprietes_particules()
{
  return variables_internes_->proprietes_particules_;
}

const Proprietes_part_vol& Transport_Interfaces_FT_Disc::proprietes_particules() const
{
  return variables_internes_->proprietes_particules_;
}

Maillage_FT_Disc& Transport_Interfaces_FT_Disc::maillage_inject()
{
  return variables_internes_->maillage_inject_;
}

const Maillage_FT_Disc& Transport_Interfaces_FT_Disc::maillage_inject() const
{
  return variables_internes_->maillage_inject_;
}

Proprietes_part_vol& Transport_Interfaces_FT_Disc::proprietes_inject()
{
  return variables_internes_->proprietes_inject_;
}

const Proprietes_part_vol& Transport_Interfaces_FT_Disc::proprietes_inject() const
{
  return variables_internes_->proprietes_inject_;
}

void Transport_Interfaces_FT_Disc::nettoyer_proprietes_particules(const ArrOfInt& som_utilises)
{
  proprietes_particules().nettoyer(som_utilises);
}


/*! @brief
 *
 */
int Transport_Interfaces_FT_Disc::sauvegarder(Sortie& os) const
{
  int bytes = Equation_base::sauvegarder(os);
  {
    int special, afaire;
    const int format_xyz = EcritureLectureSpecial::is_ecriture_special(special, afaire);
    double temps=inconnue().temps();
    Nom mon_ident("variables_internes_transport");
    mon_ident += Nom(temps,"%e");
    if (format_xyz)
      {
        if (Process::je_suis_maitre())
          {
            os << mon_ident << finl;
            os << variables_internes_->que_suis_je() << finl;
          }
      }
    else
      {
        os << mon_ident << finl;
        os << variables_internes_->que_suis_je() << finl;
      }
    bytes += variables_internes_->sauvegarder(os);
  }
  os.flush();

  return bytes;
}
int Transport_Interfaces_FT_Disc::reprendre(Entree& is)
{
  Equation_base::reprendre(is);
  {
    Nom id, type_name;
    is >> id >> type_name;
    if ( (! id.debute_par("variables_internes_transport"))
         || type_name != variables_internes_->que_suis_je())
      {
        Cerr << "Error for the method Transport_Interfaces_FT_Disc::reprendre" << finl;
        Cerr << variables_internes_->que_suis_je() <<" was expected."<< finl;
        exit();
      }
    variables_internes_->maillage_interface.associer_equation_transport(*this);
    variables_internes_->reprendre(is);
    variables_internes_->injection_interfaces_last_time_ = schema_temps().temps_courant();
  }
  return 1;
}

const Probleme_base& Transport_Interfaces_FT_Disc::get_probleme_base() const
{
  return probleme_base_.valeur();
}

const Parcours_interface& Transport_Interfaces_FT_Disc::parcours_interface() const
{
  return variables_internes_->parcours_interface_;
}

const Connectivite_frontieres& Transport_Interfaces_FT_Disc::connectivite_frontieres() const
{
  return variables_internes_->connectivite_frontieres_;
}

const Algorithmes_Transport_FT_Disc& Transport_Interfaces_FT_Disc::algorithmes_transport() const
{
  return variables_internes_->algorithmes_transport_.valeur();
}

/*! @brief Cherche le champ discret aux interfaces dont le nom est "champ", et verifie qu'il peut etre postraite a la localisation demandee (loc).
 *
 *   Si oui on renvoie 1 et, si ftab est non nul, on remplit le champ ftab
 *   avec le champ demande.
 *   Si non, on renvoie 0.
 *   (la fonction est appelee avec ftab=0 lors de la lecture du postraitement,
 *    car on n'a pas besoin de la valeur du champ, on veut seulement verifier
 *    qu'il existe).
 *
 */
int Transport_Interfaces_FT_Disc::get_champ_post_FT(const Motcle& champ, Postraitement_base::Localisation loc, DoubleTab *ftab) const
{
  int res = 1;

  const Motcle som = "sommets";            //postraitement possible uniquement aux sommets
  const Motcle elem = "elements";          //postraitement possible uniquement aux elements
  const Motcle bi = "elements et sommets"; //postraitement possible aux sommets et aux elements
  const int nb_champs = 5;
  Motcles les_champs(nb_champs);
  {
    les_champs[0] = Postraitement_base::demande_description;
    les_champs[1] = "courbure";
    les_champs[2] = "vitesse";
    les_champs[3] = "vitesse_repere_local";
    les_champs[4] = "normale_unitaire";
  }
  Motcles localisations(nb_champs);
  {
    localisations[0] = bi;
    localisations[1] = som;
    localisations[2] = som;
    localisations[3] = som;
    localisations[4] = elem;
  }

  int rang=les_champs.search(champ), i;

  if (rang==0)
    {
      Cerr<<"The real fields to be post-processed are :"<<finl;
      for (i=1 ; i<nb_champs ; i++)
        {
          Cerr << " Fields("<<i<<") : "<< les_champs[i] << " # Localisations : " << localisations[i] << finl;
        }
      res = 0;
    }
  else if (rang==-1)     //test champ existe ?
    {
      //champ inexistant
      res = 0;
    }
  else if (! (localisations[rang]==bi
              || (localisations[rang]==som && loc==Postraitement_base::SOMMETS)
              || (localisations[rang]==elem && loc==Postraitement_base::ELEMENTS)) )   //test localisation autorisee ?
    {
      //localisation non autorisee
      res = 0;
    }
  else
    {
      switch(rang)
        {
        case 1:
          {
            if (!ftab) break; // Pointeur nul : ne pas calculer la valeur du champ.
            const Maillage_FT_Disc& mesh = maillage_interface_pour_post();
            const ArrOfDouble& valeurs = mesh.get_update_courbure_sommets();
            const int n = valeurs.size_array();
            ftab->resize(n,1);
            for (int ii = 0; ii < n; ii++)
              (*ftab)(ii,0) = valeurs[ii];
            break;
          }
        case 2:
          {
            if (!ftab) break; // Pointeur nul : ne pas calculer la valeur du champ.
            if (variables_internes_->refequation_vitesse_transport.non_nul())
              {
                DoubleTabFT vit;
                // Calcul de la vitesse de deplacement des sommets par interpolation
                // (deplacement contient en fait la vitesse en m/s)
                const Equation_base& eqn_hydraulique = variables_internes_->refequation_vitesse_transport.valeur();
                const Champ_base& champ_vitesse = eqn_hydraulique.inconnue().valeur();
                int flag = 1;
                if (sub_type(Navier_Stokes_FT_Disc, eqn_hydraulique))
                  {
                    const Navier_Stokes_FT_Disc& ns = ref_cast(Navier_Stokes_FT_Disc, eqn_hydraulique);
                    if (ns.get_delta_vitesse_interface())
                      flag=0;
                  }
                calculer_vitesse_transport_interpolee(champ_vitesse, maillage_interface_pour_post(), vit,
                                                      0 /* ne pas recalculer le champ de vitesse L2 */,  // GB 2020/03/20 -> Je suis surpris par cette option 0 en desaccord avec le moment du calcul
                                                      // mais je ne prends pas la responsabilite de changer car je ne sais pas trop ce que ca represente (c'est du VEF uniquement?)
                                                      flag /* Interpolation Multi-lineaire en VDF */);
                // Pour ajouter le saut de vitesse a l'interface :
                ajouter_contribution_saut_vitesse(vit); // ici, l'interpolation depend de ns.get_new_mass_source()
                const int nb_noeuds = vit.dimension(0);
                const int nb_compo = vit.line_size();
                ftab->resize(nb_noeuds, nb_compo);
                int som2,k;
                for (som2=0 ; som2<nb_noeuds ; som2++)
                  for (k=0 ; k<nb_compo ; k++)
                    (*ftab)(som2,k) = vit(som2,k);
              }
            else
              {
                // Coder le postraitement d'une vitesse imposee
                Cerr << "Error : velocity nodes post-processing : to be developped." << finl;
                assert(0);
                Process::exit();
              }
            break;
          }
        case 3:
          {
            if (!ftab) break; // Pointeur nul : ne pas calculer la valeur du champ.
            if (variables_internes_->refequation_vitesse_transport.non_nul())
              {
                DoubleTabFT vit;
                DoubleTab Positions,Vitesses;
                // Calcul de la vitesse de deplacement des sommets par interpolation
                // (deplacement contient en fait la vitesse en m/s)
                const Equation_base& eqn_hydraulique = variables_internes_->refequation_vitesse_transport.valeur();
                const Champ_base& champ_vitesse = eqn_hydraulique.inconnue().valeur();
                calculer_vitesse_transport_interpolee(champ_vitesse, maillage_interface_pour_post(), vit,
                                                      0 /* ne pas recalculer le champ de vitesse L2 */,  // GB 2020/03/20 -> Je suis surpris par cette option 0 en desaccord avec le moment du calcul
                                                      // mais je ne prends pas la responsabilite de changer car je ne sais pas trop ce que ca represente (c'est du VEF uniquement?)
                                                      1 /* Interpolation Multi-lineaire en VDF */);
                // Pour ajouter le saut de vitesse a l'interface :
                ajouter_contribution_saut_vitesse(vit);// ici, l'interpolation depend de ns.get_new_mass_source()
                calculer_vitesse_repere_local( maillage_interface_pour_post(), vit,Positions,Vitesses);

                const int nb_noeuds = vit.dimension(0);
                const int nb_compo = vit.line_size();
                ftab->resize(nb_noeuds, nb_compo);
                int som2,k;
                for (som2=0 ; som2<nb_noeuds ; som2++)
                  for (k=0 ; k<nb_compo ; k++)
                    (*ftab)(som2,k) = vit(som2,k);
              }
            else
              {
                // Coder le postraitement d'une vitesse imposee
                Cerr << "Error : vitesse_repere_local nodes post-processing : to be developped." << finl;
                assert(0);
                exit();
              }
            break;
          }
        case 4:
          {
            if (!ftab) break; // Pointeur nul : ne pas calculer la valeur du champ.
            const Maillage_FT_Disc& mesh = maillage_interface_pour_post();
            const DoubleTab& valeurs = mesh.get_update_normale_facettes();
            const int nb_fa7 = valeurs.dimension(0);
            const int nb_compo = valeurs.line_size();
            ftab->resize(nb_fa7, nb_compo);
            int fa7,k;
            for (fa7=0 ; fa7<nb_fa7 ; fa7++)
              for (k=0 ; k<nb_compo ; k++)
                (*ftab)(fa7,k) = valeurs(fa7,k);
            break;
          }
        default:
          Cerr << "Error for the method Transport_Interfaces_FT_Disc::get_champ_post_FT" << finl;
          assert(0);
          exit();
        }
      res = 1;
    }

  return res;
}

/*! @brief Voir l'autre get_champ_post_FT.
 *
 * Cette fonction est specifique aux champs d'entiers.
 *
 */
int Transport_Interfaces_FT_Disc::get_champ_post_FT(const Motcle& champ, Postraitement_base::Localisation loc, IntTab *itab) const
{
  int res = 1;

  const Motcle som = "sommets";            //postraitement possible uniquement aux sommets
  const Motcle elem = "elements";          //postraitement possible uniquement aux elements
  const Motcle bi = "elements et sommets"; //postraitement possible aux sommets et aux elements
  const int nb_champs = 5;
  Motcles les_champs(nb_champs);
  {
    les_champs[0] = Postraitement_base::demande_description;
    les_champs[1] = "pe";        // PE owner
    les_champs[2] = "numero";    // numero local du sommet/element
    les_champs[3] = "pe_local";  // PE local
    les_champs[4] = "compo_connexe";
  }
  Motcles localisations(nb_champs);
  {
    localisations[0] = bi;
    localisations[1] = bi;
    localisations[2] = bi;
    localisations[3] = bi;
    localisations[4] = elem;
  }

  int rang=les_champs.search(champ), i;

  if (rang==0)
    {
      Cerr<<"The integer fields to be post-processed are:"<<finl;
      for (i=1 ; i<nb_champs ; i++)
        {
          Cerr << " Fields("<<i<<") : "<< les_champs[i] << " # Localisations : " << localisations[i] << finl;
        }
      res = 0;
    }
  else if (rang==-1)
    {
      //champ inexistant
      res = 0;
    }
  else if (! (localisations[rang]==bi
              || (localisations[rang]==som && loc==Postraitement_base::SOMMETS)
              || (localisations[rang]==elem && loc==Postraitement_base::ELEMENTS)) )   //test localisation autorisee ?
    {
      //localisation non autorisee
      res = 0;
    }
  else
    {
      if (itab)   // Pointeur non nul : calculer le champ
        {

          const Maillage_FT_Disc& maillage = maillage_interface_pour_post();
          const int n =
            (loc==Postraitement_base::SOMMETS)
            ? maillage.nb_sommets()
            : maillage.nb_facettes();
          int i2;
          itab->resize(n);

          switch (rang )
            {
            case 1:
              {
                if (loc==Postraitement_base::SOMMETS)
                  {
                    //TMP : tant que IntTabFt et IntTab n'ont pas fusionne :
                    //*itab = maillage_interface_->sommet_PE_owner();
                    const ArrOfInt& pe_som = maillage.sommet_PE_owner();
                    for (i2=0 ; i2<n ; i2++)
                      {
                        (*itab)(i2) = pe_som[i2];
                      }
                  }
                else
                  {
                    ArrOfIntFT pe_fac;
                    maillage.facette_PE_owner(pe_fac);
                    for (i2=0 ; i2<n ; i2++)
                      {
                        (*itab)(i2) = pe_fac[i2];
                      }
                  }
                break;
              }
            case 2:
              {
                for (i2=0 ; i2<n ; i2++)
                  {
                    (*itab)(i2) = i2;
                  }
                break;
              }
            case 3:
              {
                (*itab) = Process::me();
                break;
              }
            case 4:
              {
                maillage.intersections_elem_facettes();
                ArrOfIntFT compo(maillage.nb_facettes());
                compo = 0;
                int n2 = search_connex_components_local_FT(maillage, compo);
                compute_global_connex_components_FT(maillage, compo, n2);
                const int nbf = maillage.nb_facettes();
                for (int ii = 0; ii < nbf; ii++)
                  (*itab)[ii] = compo[ii];

                break;
              }
            default:
              Cerr << "Transport_Interfaces_FT_Disc::get_champ_post_FT : unexpected case" << finl;
              assert(0);
              exit();
            }
        }
    }
  return res;
}

/*! @brief Renvoie le maillage stocke specialement pour le postraitement (si on veut postraiter un etat intermediaire.
 *
 * ..)
 *
 */
const Maillage_FT_Disc& Transport_Interfaces_FT_Disc::maillage_interface_pour_post() const
{
  //return variables_internes_->maillage_pour_post;
  return maillage_interface();
}

const int& Transport_Interfaces_FT_Disc::get_n_iterations_distance() const
{
  return variables_internes_->n_iterations_distance;
}

/*! @brief Calcule et renvoie la distance a l'interface, evaluee sur une epaisseur egale a n_iterations_distance aux elements et discretisee aux elements
 *
 */
const Champ_base& Transport_Interfaces_FT_Disc::get_update_distance_interface() const
{
  // Si le tag du maillage et le tag du champ sont identiques, inutile de recalculer:
  const int tag = maillage_interface().get_mesh_tag();
  if (tag != variables_internes_->distance_normale_cache_tag)
    {
      DoubleTab& distance = variables_internes_->distance_interface.valeur().valeurs();
      DoubleTab& normale  = variables_internes_->normale_interface.valeur().valeurs();
      calculer_distance_interface(maillage_interface(),
                                  distance,
                                  normale,
                                  variables_internes_->n_iterations_distance);
      variables_internes_->distance_normale_cache_tag = tag;
    }
  return variables_internes_->distance_interface.valeur();
}

/*! @brief Calcule et renvoie la normale a l'interface, evaluee sur une epaisseur egale a n_iterations_distance aux elements et discretisee aux elements.
 *
 */
const Champ_base& Transport_Interfaces_FT_Disc::get_update_normale_interface() const
{
  const int tag = maillage_interface().get_mesh_tag();
  if (tag != variables_internes_->distance_normale_cache_tag)
    {
      DoubleTab& distance = variables_internes_->distance_interface.valeur().valeurs();
      DoubleTab& normale  = variables_internes_->normale_interface.valeur().valeurs();
      calculer_distance_interface(maillage_interface(),
                                  distance,
                                  normale,
                                  variables_internes_->n_iterations_distance);
      variables_internes_->distance_normale_cache_tag = tag;
    }
  return variables_internes_->normale_interface.valeur();
}

/*! @brief Renvoi de la distance signee entre l'interface et les sommets du maillage eulerien.
 *
 * Si cette distance n'a pas encore ete calculee, appel a calculer_distance_interface_sommets.
 *  C'est un DoubleTab parce qu'il n'existe pas (encore) de champ aux sommets en VDF ...
 *
 */
const DoubleTab&   Transport_Interfaces_FT_Disc::get_update_distance_interface_sommets() const
{
  // Si le tag du maillage et le tag du champ sont identiques, inutile de recalculer:
  const int tag = maillage_interface().get_mesh_tag();
  if (tag == variables_internes_->distance_sommets_cache_tag)
    {
      return variables_internes_->distance_interface_sommets;
    }
  variables_internes_->distance_sommets_cache_tag = tag;

  const DoubleTab& dist_elem = get_update_distance_interface().valeurs();
  const DoubleTab& normale_elem = get_update_normale_interface().valeurs();
  DoubleTab&        dist_som = variables_internes_->distance_interface_sommets;

  calculer_distance_interface_sommets(dist_elem, normale_elem, dist_som);
  return dist_som;
}

/*! @brief Calcule dist_som, la distance entre l'interface et les sommets du maillage eulerien a partir de dist_elem et normale_elem,
 *
 *   distance et normale a l'interface aux centres des elements euleriens.
 *   Pour un element, on evalue la distance entre chaque sommet de l'element et l'interface
 *   comme :
 *    d = d1 + d2,
 *    d2 = normale scalaire (position_sommet - centre_element)
 *    d1 est la distance entre l'interface et le centre de l'element,
 *    normale est la normale a l'interface evaluee au centre de l'element
 *   Ensuite, la distance entre un sommet et l'interface est la moyenne de toutes
 *   les distances calculee a l'aide des elements adjacents a ce sommet.
 *  La distance est invalide au-dela d'une certaine epaisseur autour de l'interface
 *  (voir iterations de lissage dans calculer_distance_interface).
 *  Dans ce cas on met une distance de +1e30 si l'indicatrice est >0.5,
 *  sinon on met -1e30 (ce choix permet d'utiliser la fonction
 *  distance dans les marching-cubes sans avoir a calculer une vraie distance partout).
 *  Parametre : dist_elem
 *  Signification : tableau contenant pour chaque element reel et virtuel la distance
 *                  entre l'interface et le centre de l'element (calculee par
 *                  calculer_distance_interface). L'espace virtuel doit etre a jour.
 *  Parametre : normale_elem
 *  Signification : idem pour la normale a l'interface
 *  Parametre : dist_som
 *  Signification : tableau ou on stocke le resultat du calcul. Le tableau doit
 *                  avoir la bonne taille et un descripteur adequat (voir "discretiser",
 *                  a priori un tableau avec une epaisseur de joint de zero et uniquement
 *                  des items communs).
 *
 */
void Transport_Interfaces_FT_Disc::calculer_distance_interface_sommets(
  const DoubleTab& dist_elem,
  const DoubleTab& normale_elem,
  DoubleTab&        dist_som) const
{
  static const double distance_sommets_invalides = -1.e30;

  const Domaine_VF&    domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
  const IntTab&     elem_som = domaine_vf.domaine().les_elems();
  const DoubleTab& xp = domaine_vf.xp();
  const DoubleTab& coord_som = domaine_vf.domaine().les_sommets();

  const int nb_sommets = dist_som.dimension_tot(0);
  ArrOfInt ncontrib(nb_sommets);
  ncontrib = 0;
  dist_som = 0.;

  const int dim = Objet_U::dimension;
  double centre[3] = {0., 0., 0.};
  double normale[3] = {0., 0., 0.};
  // Calcul de SOMME(d1+d2) pour tous les elements voisins de chaque sommet :
  int elem, i;
  const int nb_elem_tot = dist_elem.dimension_tot(0);
  const int nb_som_elem = elem_som.line_size();

  for (elem = 0; elem < nb_elem_tot; elem++)
    {
      const double d1 = dist_elem(elem);
      // Si la distance est invalide, on passe:
      if (d1 < distance_sommets_invalides)
        continue;
      // Centre de l'element et normale a l'interface pour cet element :
      for (i = 0; i < dim; i++)
        {
          centre[i] = xp(elem, i);
          normale[i] = normale_elem(elem, i);
        }
      // Boucle sur les sommets de l'element
      for (i = 0; i < nb_som_elem; i++)
        {
          const int som = elem_som(elem, i);
          // dist_som ne contient que des sommets reels en general.
          // si le sommet n'est pas dans dist_som, on ne calcule pas.
          if (som < nb_sommets)
            {
              double d2 = 0.;
              int j;
              for (j = 0; j < dim; j++)
                {
                  double position_sommet = coord_som(som, j);
                  d2 += (position_sommet - centre[j]) * normale[j];
                }
              dist_som(som) += d1 + d2;
              ++(ncontrib[som]);
            }
        }
    }
  // Division par le nombre d'elements voisins
  const double valeur_invalide = distance_sommets_invalides * 1.1;
#ifdef __INTEL_COMPILER
#pragma novector // Desactive vectorisation sur Intel car crash sinon
#endif
  for (i = 0; i < nb_sommets; i++)
    {
      const int n = ncontrib[i];
      if (n > 0)
        dist_som(i) /= n;
      else
        dist_som(i) = valeur_invalide;
    }
  dist_som.echange_espace_virtuel();
  Debog::verifier("Transport_Interfaces_FT_Disc::calculer_distance_interface_sommets",dist_som);
}

/*! @brief Calcul d'un champ scalaire aux elements contenant une distance signee entre le centre de l'element et l'interface.
 *
 * La distance est positive dans
 *   la phase 1 et negative dans la phase 0.
 *   On calcule aussi un champ vectoriel aux elements contenant une normale
 *   a l'interface. Ce champ est evalue en resolvant moralement
 *    laplacien(normale) = gradient(indicatrice)
 *   ou gradient(indicatrice) est le gradient de l'indicatrice continue
 *   c'est a dire un dirac localise a la surface de l'interface.
 *   Pour l'instant, cette normale est calculee de facon approchee avec quelques
 *   iterations d'un lisseur. Le support est donc limite au voisinage de l'interface.
 *   Pour les autres elements, la distance vaut -1.e30
 *  Precondition : le maillage doit etre parcouru
 *
 */
void Transport_Interfaces_FT_Disc::calculer_distance_interface(
  const Maillage_FT_Disc& maillage,
  DoubleTab& distance_elements,
  DoubleTab& normale_elements,
  const int n_iter) const
{
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "Calculer_distance_interface");
  statistiques().begin_count(stat_counter);

  static const double distance_sommets_invalides = -1.e30;

  // Coordonnees des sommets du maillage eulerien:
  const Domaine_dis_base& mon_dom_dis = domaine_dis().valeur();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, mon_dom_dis);
  const DoubleTab& centre_element = domaine_vf.xp();

  // Tableau contenant une approximation de la normale aux sommets du maillage
  // eulerien :
  distance_elements = distance_sommets_invalides * 1.1;
  normale_elements  = 0.;

  const int nb_elem = mon_dom_dis.domaine().nb_elem();
  const int dim = Objet_U::dimension;

  // Calcul de la distance pour l'epaisseur 0 (sommets des elements traverses par
  // l'interface). Pour chaque element, on calcule le plan passant par
  // le centre de gravite des portions de facettes et dont la normale est
  // la moyenne des normales aux portions de facettes, ponderees par la surface
  // des portions de facettes. La distance interface/element est la distance
  // entre ce plan et le centre de l'element.
  {
    const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
    const ArrOfInt& index_elem = intersections.index_elem();
    const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();
    const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
    const IntTab& facettes = maillage.facettes();
    const DoubleTab& sommets = maillage.sommets();
    // Boucle sur les elements
    for (int elem = 0; elem < nb_elem; elem++)
      {
        int index = index_elem[elem];
        // Moyenne ponderee des normales aux facettes qui traversent l'element
        double normale[3] = {0., 0., 0.};
        // Centre de gravite de l'intersection facettes/element
        double centre[3] = {0., 0., 0.};
        // Somme des poids
        double surface_tot = 0.;
        // Boucle sur les facettes qui traversent cet element
        while (index >= 0)
          {
            const Intersections_Elem_Facettes_Data& data =
              intersections.data_intersection(index);

            const int num_facette = data.numero_facette_;
#ifdef AVEC_BUG_SURFACES
            const double surface = data.surface_intersection_;
#else
            const double surface = data.fraction_surface_intersection_ * surface_facettes[num_facette];
#endif
            surface_tot += surface;
            for (int i = 0; i < dim; i++)
              {
                normale[i] += surface * normale_facettes(num_facette, i);
                // Calcul du centre de gravite de l'intersection facette/element
                double g_i = 0.; // Composante i de la coordonnee du centre de gravite
                for (int j = 0; j < dim; j++)
                  {
                    const int som   = facettes(num_facette, j);
                    const double coord = sommets(som, i);
                    const double coeff = data.barycentre_[j];
                    g_i += coord * coeff;
                  }
                centre[i] += surface * g_i;
              }
            index = data.index_facette_suivante_;
          }
        if (surface_tot > 0.)
          {
            // La normale stockee n'est pas normee : norme a peu pres egale a la surface
            // d'interface qui a servi a la calculer.
            // centre = somme(centre[facette] * surface) / somme(surface)
            const double inverse_surface_tot = 1. / surface_tot;
            double norme = 0.;
            int j;
            for (j = 0; j < dim; j++)
              {
                norme += normale[j] * normale[j];
                normale_elements(elem, j) = normale[j];
                centre[j] *= inverse_surface_tot;
              }
            if (norme > 0)
              {
                double i_norme = 1./sqrt(norme);
                double distance = 0.;
                for (j = 0; j < dim; j++)
                  {
                    double n_j = normale[j] * i_norme; // normale normee
                    distance += (centre_element(elem, j) - centre[j]) * n_j;
                  }
                distance_elements(elem) = distance;
              }
          }
      }
    normale_elements.echange_espace_virtuel();
    distance_elements.echange_espace_virtuel();
    Debog::verifier("Transport_Interfaces_FT_Disc::calculer_distance_interface",normale_elements);
    Debog::verifier("Transport_Interfaces_FT_Disc::calculer_distance_interface",distance_elements);
  }

  DoubleTab terme_src(normale_elements);
  DoubleTab tmp(normale_elements);

  const IntTab& face_voisins = domaine_vf.face_voisins();
  const IntTab& elem_faces   = domaine_vf.elem_faces();
  const int nb_elem_voisins = elem_faces.line_size();

  // Calcul d'une normale aux elements :
  int iteration;
  for (iteration = 0; iteration < n_iter; iteration++)
    {
      // Iteration du lisseur : moralement on fait
      //  normale = normale + (terme_source - laplacien(normale)) * facteur
      // qui converge vers
      //  laplacien(normale) = terme_source
      // mais le laplacien est pipo : en pratique on fait juste :
      //  normale = moyenne(normale sur les elem voisins) + terme_source

      const double un_sur_ncontrib = 1. / (1. + nb_elem_voisins);
      int elem, i, k;
      for (elem = 0; elem < nb_elem; elem++)
        {
          // La moyenne du vecteur normal sur les voisins:
          double n[3] = {0., 0., 0.};
          for (i = 0; i < dim; i++)
            n[i] = normale_elements(elem, i);
          for (k = 0; k < nb_elem_voisins; k++)
            {
              // On cherche l'element voisin par la face k
              const int face = elem_faces(elem, k);
              const int e_voisin = face_voisins(face, 0) + face_voisins(face, 1) - elem;
              if (e_voisin >= 0) // Si on n'est pas au bord...
                for (i = 0; i < dim; i++)
                  n[i] += normale_elements(e_voisin, i);
            }
          for (i = 0; i < dim; i++)
            tmp(elem, i) = terme_src(elem, i) + n[i] * un_sur_ncontrib;
        }
      normale_elements = tmp;
      normale_elements.echange_espace_virtuel();
    }
  // On normalise la normale, et on cree une liste des elements pour lesquels
  // la normale est connue :
  ArrOfIntFT liste_elements;
  {
    int elem;
    for (elem = 0; elem < nb_elem; elem++)
      {
        double nx = normale_elements(elem, 0);
        double ny = normale_elements(elem, 1);
        double nz = (dim==3) ? normale_elements(elem, 2) : 0.;
        double norme2 = nx*nx + ny*ny + nz*nz;
        if (norme2 > 0.)
          {
            double i_norme = 1. / sqrt(norme2);
            normale_elements(elem, 0) = nx * i_norme;
            normale_elements(elem, 1) = ny * i_norme;
            if (dim==3)
              normale_elements(elem, 2) = nz * i_norme;
            liste_elements.append_array(elem);
          }
      }
    normale_elements.echange_espace_virtuel();//Manque cet echange (MB)
  }
  // Calcul d'une distance a l'interface :
  terme_src = distance_elements;
  tmp = distance_elements;
  for (iteration = 0; iteration < n_iter; iteration++)
    {
      int i_elem, elem;
      const int liste_elem_size = liste_elements.size_array();
      for (i_elem = 0; i_elem < liste_elem_size; i_elem++)
        {
          elem = liste_elements[i_elem];
          if (terme_src(elem) > distance_sommets_invalides)
            {
              // Pour les elements traverses par l'interface, la distance n'est pas recalculee.
              tmp(elem) = distance_elements(elem);
            }
          else
            {
              // Pour les autres, on calcule une distance pour chaque element voisin
              double ncontrib = 0.;
              double somme_distances = 0.;
              int k;
              for (k = 0; k < nb_elem_voisins; k++)
                {
                  // On cherche l'element voisin par la face k
                  const int face = elem_faces(elem, k);
                  const int e_voisin = face_voisins(face, 0) + face_voisins(face, 1) - elem;
                  if (e_voisin >= 0) // Si on n'est pas au bord...
                    {
                      const double distance_voisin = distance_elements(e_voisin);
                      if (distance_voisin > distance_sommets_invalides)
                        {
                          // Calcul d'une normale moyenne entre l'element et le voisin
                          double nx = normale_elements(elem, 0) + normale_elements(e_voisin, 0);
                          double ny = normale_elements(elem, 1) + normale_elements(e_voisin, 1);
                          double nz = (dim==3)
                                      ? normale_elements(elem, 2) + normale_elements(e_voisin, 2) : 0.;
                          double norm2 = nx*nx + ny*ny + nz*nz;
                          if (norm2 > 0.)
                            {
                              double i_norm = 1./sqrt(norm2);
                              nx *= i_norm;
                              ny *= i_norm;
                              nz *= i_norm;
                            }
                          // Calcul du vecteur (element - element_voisin)
                          double dx = centre_element(elem, 0) - centre_element(e_voisin, 0);
                          double dy = centre_element(elem, 1) - centre_element(e_voisin, 1);
                          double dz = (dim==3)
                                      ? centre_element(elem, 2) - centre_element(e_voisin, 2) : 0.;
                          double d = nx * dx + ny * dy + nz * dz + distance_voisin;
                          somme_distances += d;
                          ncontrib++;
                        }
                    }
                }
              // Moyenne des distances obtenues avec les elements voisins.
              if (ncontrib > 0.)
                {
                  double d = somme_distances / ncontrib;
                  tmp(elem) = d;
                }
            }
        }
      distance_elements = tmp;
      distance_elements.echange_espace_virtuel();
    }
  statistiques().end_count(stat_counter);
}

/*! @brief Calcul de la derivee par rapport au temps du volume de phase 1 aux sommets du maillage lagrangien a partir du champ de vitesse
 *
 *   eulerien. On utilise le fait que le champ eulerien est a divergence
 *   nulle
 *   Cette grandeur permet de corriger le deplacement des sommets pour
 *   conserver le volume des phases (voir these B.M. paragraphe 3.2.10).
 *   On a I = rho_0 + (rho_1-rho_0) * I. Donc:
 *    drho/dt = (rho_1-rho_0) * dI/dt  (d'une part)
 *            = div(rho*u) = div( (rho_0 + (rho_1-rho_0)*I) * u)  (d'autre part)
 *   Donc, si div(u) = 0
 *    dI/dt = div(I * u)
 *   Si non (changement de phase):
 *    dI/dt = div((rho_0/(rho_1-rho_0) + I) * u)
 *  Parametre : vitesse
 *  Signification: le champ de vitesse eulerien aux faces du maillage.
 *                 si vitesse.dimension(1)==1, on suppose que c'est la
 *                 composante normale a la face de la vitesse, sinon
 *                 on suppose que c'est le vecteur vitesse 2d ou 3d.
 *  Parametre : indicatrice
 *  Signification: indicatrice de phase aux elements euleriens.
 *                 doit avoir son espace virtuel a jour et correspondre
 *                 au maillage suivant...
 *  Parametre : maillage
 *  Signification : le maillage de l'interface,
 *                  doit etre parcouru.
 *  Parametre : rho_0_sur_delta_rho_div_u
 *  Signification: resultat de l'operateur div(u) aux elements (integrale de div_u
 *                 sur les elements. Si ce tableau est de taille non nulle, on ajoute
 *                 le terme en div_u
 *  Parametre : var_volume
 *  Signification : Tableau ou on stocke la variation de volume de phase 1
 *                  pour chaque sommet du maillage lagrangien.
 *
 */
#if 0
void Transport_Interfaces_FT_Disc::calculer_derivee_volume_phase1(
  const DoubleTab& vitesse,
  const DoubleTab& indicatrice,
  const Maillage_FT_Disc& maillage,
  const DoubleTab& rho_0_sur_delta_rho_div_u,
  ArrOfDouble& var_volume) const
{
  const int nb_sommets = maillage.nb_sommets();
  var_volume.resize_array(nb_sommets);
  // Initialisation a zero car on va ajouter des contributions dans le desordre :
  var_volume = 0.;

  // Ancien codage (ne supporte pas le changement de phase)

  // Recuperation des tableaux et autres raccourcis pour le calcul :
  const Domaine_dis_base& domaine_dis_base = domaine_dis().valeur();
  const Domaine_VF&        domaine_vf       = ref_cast(Domaine_VF, domaine_dis_base);
  const Domaine&           domaine          = domaine_dis_base.domaine();
  const int nb_faces_element = domaine.nb_faces_elem();
  const int dim = Objet_U::dimension;
  // En vef : 2 ou 3 composantes de vitesse,
  // en vdf : composante normale a la face uniquement.
  const int   vitesse_n_composantes = (vitesse.line_size() > 1) ? 1 : 0; // 1 si VEF et 0 si VDF
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const IntTab& face_voisins = domaine_vf.face_voisins();
  // Pour le vdf, on a besoin de la surface des faces...
  const DoubleVect * face_surfaces_ptr = 0;
  if (sub_type(Domaine_VDF, domaine_vf))
    {
      const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_vf);
      face_surfaces_ptr = & (domaine_vdf.face_surfaces());
    }
  const IntTab& facettes = maillage.facettes();
  const Intersections_Elem_Facettes& intersections =
    maillage.intersections_elem_facettes();
  const ArrOfInt& index_elem = intersections.index_elem();
  // Verification: soit on est en vdf avec une seule composante,
  //  soit on n'est pas en vdf, avec "dimension" composantes :
  assert(vitesse_n_composantes!=0 || face_surfaces_ptr!=0);
  assert(vitesse_n_composantes==0 || vitesse.line_size()==dim);
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();

  // Boucle sur les elements euleriens
  const int nb_elements = domaine.nb_elem();
  int element;
  for (element = 0; element < nb_elements; element++)
    {

      // Si cet element n'est pas traverse par l'interface, ne pas calculer
      // la variation de volume.
      const int index_premiere_intersection = index_elem[element];
      if (index_premiere_intersection < 0) // Pas de facette dans cet element
        continue;

      // Calcul du bilan de volume de la phase 1 sur cet element :
      // integrale sur les faces de vitesse * normale * indicatrice_face.
      // indicatrice_face est l'indicatrice a la face: 1 si l'element voisin
      // a une indicatrice = 1, 0 si l'element voisin a une indicatrice 0,
      // indicatrice moyenne sinon (voir these B.M. paragraphe 3.2.10)

      // Calcul de la derivee par rapport au temps du volume de phase1 dans l'element :
      double derivee_volume_phase1 = 0.;
      int face_locale;
      for (face_locale = 0; face_locale < nb_faces_element; face_locale++)
        {
          const int face = elem_faces(element, face_locale);
          const int elem_voisin_0 = face_voisins(face, 0);
          const int elem_voisin_1 = face_voisins(face, 1);
          const int elem_voisin = elem_voisin_0 + elem_voisin_1 - element;
          double indicatrice_face;
          if (elem_voisin < 0)
            {
              // Face de bord
              // Si c'est un bord ferme (vitesse normale nulle), alors
              // la valeur n'a pas d'importance. Si c'est un bord ouvert,
              // alors je ne connais pas de formulation qui permette de
              // conserver le volume des phases.
              // Choix d'une valeur pas debile :
              indicatrice_face = indicatrice(element);
            }
          else
            {
              const double indicatrice_voisin = indicatrice(elem_voisin);
              if (indicatrice_voisin == 0.)
                indicatrice_face = 0.;
              else if (indicatrice_voisin == 1.)
                indicatrice_face = 1.;
              else
                {
                  double indic = indicatrice(element);
                  indicatrice_face = (indic + indicatrice_voisin) * 0.5;
                }
            }
          // Composante normale de la vitesse a la face :
          // (la normale pointe de l'elem_voisin_0 vers l'elem_voisin_1,
          //  sa norme est egale a la surface de la face).
          double vitesse_normale = 0;
          if (vitesse_n_composantes)
            {
              // c'est le vef en general :
              // Produit scalaire normale*vitesse
              int i;
              for (i = 0; i < dim; i++)
                {
                  double n = domaine_vf.face_normales(face, i);
                  double v = vitesse(face, i);
                  vitesse_normale += v * n;
                }
            }
          else
            {
              // c'est le vdf : facile :
              double surface = (*face_surfaces_ptr)(face);
              vitesse_normale = vitesse(face) * surface;
            }
          // Calcul du debit volumique de phase 1 entrant dans l'element par cette face
          // Si la normale est dirigee vers l'exterieur de l'element, on change de signe :
          double debit_entrant_phase1 =
            (elem_voisin_0 == element) ? -vitesse_normale : vitesse_normale;
          debit_entrant_phase1 *= indicatrice_face;
          // Contribution au bilan de volume de phase 1 dans l'element :
          derivee_volume_phase1 += debit_entrant_phase1;
        }

      // Repartition conservative de cette derivee de volume sur les noeuds
      // du maillage lagrangien :

      // Premier passage : calcul de la surface totale d'intersection entre le
      // maillage lagrangien et l'element eulerien "element".
      double surface_totale = 0.;
      // Boucle sur les faces qui traversent l'element:
      int index = index_premiere_intersection;
      while (index >= 0)
        {
          const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
          const int facette = data.numero_facette_;
          const double surface_facette = surface_facettes[facette];
#ifdef AVEC_BUG_SURFACES
          surface_totale += data.surface_intersection_;
#else
          surface_totale += data.fraction_surface_intersection_ * surface_facette;
#endif
          index = data.index_facette_suivante_;
        }
      // Deuxieme passage : repartition de la derivee de volume sur les sommets,
      // proportionnelle a la surface et a la coordonnee barycentrique du sommet :
      index = index_premiere_intersection;
      if (surface_totale > 0.)
        {
          const double inv_surface_tot = 1. / surface_totale;
          while (index >= 0)
            {
              const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
              const int facette = data.numero_facette_;
              const double surface_facette = surface_facettes[facette];
              // Fraction de la surface d'intersection de cette facette a la surface totale
              // d'intersection avec l'element :
#ifdef AVEC_BUG_SURFACES
              const double fraction_surface = data.surface_intersection_ * inv_surface_tot;
#else
              const double fraction_surface = data.fraction_surface_intersection_ * surface_facette * inv_surface_tot;
#endif
              int i;
              for (i = 0; i < dim; i++)
                {
                  const int sommet = facettes(facette, i);
                  const double coeff  = data.barycentre_[i];
                  // La somme des "coeff" pour une intersection vaut 1
                  // et la sommet des fraction_surface pour toutes les intersections
                  // vaut 1, donc on distribue "derivee_volume_phase1" de facon conservative
                  // sur les sommets du maillage lagrangien :
                  const double valeur = coeff * fraction_surface * derivee_volume_phase1;
                  // Certaines contributions sont ajoutees a des sommets virtuels,
                  // on collecte le tout a la fin :
                  var_volume[sommet] += valeur;
                }
              index = data.index_facette_suivante_;
            }
        }
    }

  // Collection des donnees des sommets virtuels :
  maillage.desc_sommets().collecter_espace_virtuel(var_volume, Descripteur_FT::SOMME);
  // Mise a jour des espaces virtuels :
  maillage.desc_sommets().echange_espace_virtuel(var_volume);
}

#endif


void Transport_Interfaces_FT_Disc::calculer_vmoy_composantes_connexes(const Maillage_FT_Disc& maillage,
                                                                      const ArrOfInt& compo_connexes_facettes,
                                                                      const int nb_compo_tot,
                                                                      const DoubleTab& vitesse_sommets,
                                                                      DoubleTab& vitesses,
                                                                      DoubleTab& positions) const
{
  assert(nb_compo_tot == vitesses.dimension(0));
  assert(nb_compo_tot == positions.dimension(0));

  const int dim = vitesses.line_size();
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();
  const IntTab& facettes = maillage.facettes();
  const DoubleTab& sommets = maillage.sommets();
  assert(facettes.line_size() == dim);

  // Surface totale de chaque composante connexe, initialise a zero
  ArrOfDouble surfaces_compo(nb_compo_tot);
  positions = 0.;

  // Calcul du centre de gravite de la composante connexe
  //  (centre de gravite de la surface, pas du volume)
  const int nb_facettes_tot = facettes.dimension_tot(0);
  {
    for (int i = 0; i < nb_facettes_tot; i++)
      {
        if (maillage.facette_virtuelle(i))
          continue;
        const int compo = compo_connexes_facettes[i];
        const double surface = surface_facettes[i];
        surfaces_compo[compo] += surface;
        // Centre de gravite de la facette, pondere par la surface
        for (int j = 0; j < dim; j++)
          {
            // Indice du sommet
            const int s = facettes(i, j);
            for (int k = 0; k < dim; k++)
              // On divisera par dim a la fin:
              positions(compo, k) += surface * sommets(s, k);
          }
      }
    mp_sum_for_each_item(surfaces_compo);
    mp_sum_for_each_item(positions);

    positions *= (1. / dim);

    DoubleVect s; // tab_divide prend DoubleVect, pas ArrOfDouble...
    s.ref_array(surfaces_compo);
    tab_divide_any_shape(positions, s);
  }
  ArrOfInt som(3);
  // Calcul de la vitesse de deplacement moyenne
  vitesses = 0.;
  //calcul de la vitesse moyenne de deplacement de l'interface
  for (int fa7 = 0; fa7 < nb_facettes_tot; fa7++)
    {

      if (maillage.facette_virtuelle(fa7))
        continue;
      const double si=surface_facettes[fa7];
      som=0;

      for (int d=0; d<dim; d++)
        som[d]=facettes(fa7,d);

      //calcul de ds_dt
      double ds_dt=0.;

      if(dim==2)
        {
          const double n0=normale_facettes(fa7,0);
          const double n1=normale_facettes(fa7,1);
          ds_dt=-n1*vitesse_sommets(som[0],0)+n0*vitesse_sommets(som[0],1)+n1*vitesse_sommets(som[1],0)-n0*vitesse_sommets(som[1],1);
        }
      else
        {
          const double n0=normale_facettes(fa7,0)*0.5;
          const double n1=normale_facettes(fa7,1)*0.5;
          const double n2=normale_facettes(fa7,2)*0.5;
          double s2s1[3],d_surface[3];

          for (int i = 0; i < 3; i++)
            {
              // La differentielle de surface pour un deplacement du sommet i
              // est le produit vectoriel de la normale par le vecteur
              // s2s1 = (sommet[(i+1)%3] - sommet[(i+2)%3]) * 0.5
              // (vecteur de norme "base du triangle * 0.5" et de direction
              //  la hauteur du triangle)
              const int s0 = som[i];
              const int s1 = som[ (i+1)%3 ];
              const int s2 = som[ (i+2)%3 ];

              s2s1[0] = sommets(s1,0) - sommets(s2,0);
              s2s1[1] = sommets(s1,1) - sommets(s2,1);
              s2s1[2] = sommets(s1,2) - sommets(s2,2);

              d_surface[0] = s2s1[1] * n2 - s2s1[2] * n1;
              d_surface[1] = s2s1[2] * n0 - s2s1[0] * n2;
              d_surface[2] = s2s1[0] * n1 - s2s1[1] * n0;

              ds_dt+=d_surface[0]*vitesse_sommets(s0,0) +        d_surface[1]*vitesse_sommets(s0,1) +d_surface[2]*vitesse_sommets(s0,2);
            }
        }

      const int compo = compo_connexes_facettes[fa7];
      for (int d=0; d<dim; d++)
        {
          double V=0.;
          double p=0.;
          for (int d1=0; d1<dim; d1++)
            {
              V+=vitesse_sommets(som[d1],d);
              p+=sommets(som[d1],d);
            }
          V /= dim;
          p /= dim;
          vitesses(compo, d) += si * V + (p-positions(compo, d)) * ds_dt;
        }
    }

  mp_sum_for_each_item(vitesses);
  {
    DoubleVect s; // tab_divide prend DoubleVect, pas ArrOfDouble...
    s.ref_array(surfaces_compo);
    tab_divide_any_shape(vitesses, s);
  }
}

void Transport_Interfaces_FT_Disc::ramasse_miettes(const Maillage_FT_Disc& maillage,
                                                   DoubleVect& flux,
                                                   DoubleVect& valeurs)
{
  // Calcul d'un flux a travers chaque face, proportionnel a
  //  Surface_face * (normale_interface scalaire normale_face) * grandeur_amont_a_transporter
  const int dim = Objet_U::dimension;
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());
  const Domaine& domaine = domaine_vf.domaine();
  const IntTab&   face_voisins = domaine_vf.face_voisins();
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const int nb_faces_tot = domaine_vf.nb_faces_tot();
  //const int nb_elem = domaine.nb_elem();
  const int nb_elem_tot = domaine.nb_elem_tot();
  const int nb_faces_elem = domaine_vf.domaine().nb_faces_elem();
  const DoubleVect& indic = get_update_indicatrice().valeurs();
  const DoubleTab& normale_interface =  get_update_normale_interface().valeurs();
  for (int i_face = 0; i_face < nb_faces_tot; i_face++)
    {
      double f;
      const int elem0 = face_voisins(i_face, 0);
      const int elem1 = face_voisins(i_face, 1);
      if (elem0 < 0 || elem1 < 0)
        f = 0.;
      else
        {
          double prod_scal = 0.;
          // Le produit scalaire est fait avec la normale a l'interface evaluee
          // a la face (moyenne des normales n0 et n1 des deux elements voisins.
          // faces_normales() donne une normale de norme egale a la surface de la face
          for (int i = 0; i < dim; i++)
            {
              double n0 = normale_interface(elem0, i);
              double n1 = normale_interface(elem1, i);
              prod_scal += (n0 + n1) * 0.5 * domaine_vf.face_normales(i_face, i);
            }
          if (indic(elem0) == 1. || indic(elem1) == 1.)
            prod_scal = - prod_scal;
          // Si la face appartient a deux elements diphasiques alors
          // on annule le flux pour eviter les echanges entre ces cellules
          if (indic(elem0) != 0. && indic(elem0) != 1.)
            {
              if (indic(elem1) != 0. && indic(elem1) != 1.)
                {
                  prod_scal=0.;
                }
            }
          f = prod_scal;
        }
      flux(i_face) = f;
    }
  //
  flux.echange_espace_virtuel();
  //
  // Renormalisation des flux de sorte que la somme des flux sortants de chaque
  //  element vide completement cet element.
  // On travaille sur une copie du tableau flux pour eviter les interactions
  // entre les deux boucles qui suivent
  DoubleVect tmp_flux = flux;
  ArrOfInt flag(nb_faces_elem);
  for (int i_elem = 0; i_elem < nb_elem_tot; i_elem++)
    {
      // Somme des flux sortants de l'element:
      double somme = 0.;
      for (int j = 0; j < nb_faces_elem; j++)
        {
          int i_face = elem_faces(i_elem, j);
          const int elem0 = face_voisins(i_face, 0);
          const int elem1 = face_voisins(i_face, 1);
          if (elem0 < 0 || elem1 < 0)
            continue;
          // Le flux est-il sortant ?
          const double signe = (elem0 == i_elem) ? 1. : -1.;
          const double f = flux(i_face);
          if (f * signe > 0.)
            {
              somme += f * signe; // oui
              flag[j] = 1;//
            }
        }
      // renormalisation des flux sortants et multiplication par la valeur a transporter:
      if (somme > 0.)
        {
          const double facteur = valeurs(i_elem) / somme;
          for (int j = 0; j < nb_faces_elem; j++)
            {
              int i_face = elem_faces(i_elem, j);
              if (flag[j])
                {
                  tmp_flux[i_face] = tmp_flux[i_face] * facteur;
                  flag[j] = 0;
                }
            }
        }
    }
  // Recopie de tmp_flux dans flux
  for (int i_elem = 0; i_elem < nb_elem_tot; i_elem++)
    {
      for (int j = 0; j < nb_faces_elem; j++)
        {
          int i_face = elem_faces(i_elem, j);
          flux[i_face]=tmp_flux[i_face];
        }
    }
  // Application du flux calcule
  // nb_elem_tot a la place de nb_elem
  for (int i_elem = 0; i_elem < nb_elem_tot; i_elem++)
    {
      double x = valeurs(i_elem);
      for (int j = 0; j < nb_faces_elem; j++)
        {
          int i_face = elem_faces(i_elem, j);
          double f = flux(i_face);
          double signe = (i_elem == face_voisins(i_face, 1)) ? 1. : -1.;
          x += f * signe;
        }
      valeurs(i_elem) = x;

    }
  valeurs.echange_espace_virtuel();
}

// Calcule des valeurs aux sommets du maillage ft a partir des valeurs euler qui sont supposees etre aux elements.
// Le transfert est conservatif: on suppose que les valeurs_euler sont des integrales sur chaque element.
// La somme de ces valeurs sera egale a la somme des valeurs associees aux sommets lagrangiens.
void Transport_Interfaces_FT_Disc::transfert_conservatif_eulerien_vers_lagrangien_sommets(const Maillage_FT_Disc& maillage,
                                                                                          const DoubleVect& valeurs_euler,
                                                                                          ArrOfDouble& valeurs_lagrange)
{
  const int nb_sommets = maillage.nb_sommets();
  valeurs_lagrange.resize_array(nb_sommets);
  // Initialisation a zero car on va ajouter des contributions dans le desordre :
  valeurs_lagrange = 0.;

  const int dim = Objet_U::dimension;
  const IntTab& facettes = maillage.facettes();
  const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
  const ArrOfInt& index_elem = intersections.index_elem();
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const int nb_elements = valeurs_euler.size();
  assert(nb_elements == index_elem.size_array());

  // Boucle sur les elements euleriens
  int element;
  for (element = 0; element < nb_elements; element++)
    {

      // Si cet element n'est pas traverse par l'interface, ne pas calculer
      // la variation de volume.
      const int index_premiere_intersection = index_elem[element];
      if (index_premiere_intersection < 0) // Pas de facette dans cet element
        continue;

      // Repartition conservative de la valeur_euler sur les noeuds
      // du maillage lagrangien :

      // Premier passage : calcul de la surface totale d'intersection entre le
      // maillage lagrangien et l'element eulerien "element".
      double surface_totale = 0.;
      // Boucle sur les faces qui traversent l'element:
      int index = index_premiere_intersection;
      while (index >= 0)
        {
          const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
          const int facette = data.numero_facette_;
          const double surface_facette = surface_facettes[facette];
          surface_totale += data.fraction_surface_intersection_ * surface_facette;
          index = data.index_facette_suivante_;
        }

      const double valeur_euler = valeurs_euler[element];

      // Deuxieme passage : repartition de la valeur_euler sur les sommets,
      // proportionnelle a la surface et a la coordonnee barycentrique du sommet :
      index = index_premiere_intersection;
      if (surface_totale > 0.)
        {
          const double inv_surface_tot = 1. / surface_totale;
          while (index >= 0)
            {
              const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
              const int facette = data.numero_facette_;
              const double surface_facette = surface_facettes[facette];
              // Fraction de la surface d'intersection de cette facette a la surface totale
              // d'intersection avec l'element :
              const double fraction_surface = data.fraction_surface_intersection_ * surface_facette * inv_surface_tot;
              int i;
              for (i = 0; i < dim; i++)
                {
                  const int sommet = facettes(facette, i);
                  const double coeff  = data.barycentre_[i];
                  // La somme des "coeff" pour une intersection vaut 1
                  // et la sommet des fraction_surface pour toutes les intersections
                  // vaut 1, donc on distribue "valeur_euler" de facon conservative
                  // sur les sommets du maillage lagrangien :
                  const double valeur = coeff * fraction_surface * valeur_euler;
                  // Certaines contributions sont ajoutees a des sommets virtuels,
                  // on collecte le tout a la fin :
                  valeurs_lagrange[sommet] += valeur;
                }
              index = data.index_facette_suivante_;
            }
        }
    }

  maillage.desc_sommets().collecter_espace_virtuel(valeurs_lagrange, MD_Vector_tools::EV_SOMME);
  // Mise a jour des espaces virtuels :
  maillage.desc_sommets().echange_espace_virtuel(valeurs_lagrange);
}
