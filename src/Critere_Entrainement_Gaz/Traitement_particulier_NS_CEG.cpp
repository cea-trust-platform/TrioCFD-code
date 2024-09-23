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
// File:        Traitement_particulier_NS_CEG.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Critere_Entrainement_Gaz/src
// Version:     /main/23
//
//////////////////////////////////////////////////////////////////////////////

#include <Traitement_particulier_NS_CEG.h>
#include <Navier_Stokes_Turbulent.h>
#include <Domaine_VF.h>
#include <Probleme_base.h>
#include <Schema_Temps_base.h>
#include <sys/stat.h>
#include <communications.h>
#include <Param.h>
#include <Milieu_base.h>
#include <SFichier.h>
#include <Domaine.h>
#include <TRUSTList.h>
#include <Statistiques.h>
#include <stat_counters.h>
#include <Fluide_Incompressible.h>
#include <Modele_turbulence_hyd_RANS_K_Eps_base.h>
#include <Pb_Hydraulique_Turbulent.h>
#include <Champ_P1NC.h>
#include <Postraitement.h>

Implemente_base_sans_constructeur_ni_destructeur(Traitement_particulier_NS_CEG,"Traitement_particulier_NS_CEG",Traitement_particulier_NS_base);

/*! @brief
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Traitement_particulier_NS_CEG::printOn(Sortie& is) const
{
  return is;
}


/*! @brief
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Traitement_particulier_NS_CEG::readOn(Entree& is)
{
  return is;
}

inline void error(const Nom& message)
{
  Cerr << finl;
  Cerr << "=========================================================" << finl;
  Cerr << "Error in the use of Traitement_particulier CEG !" << finl;
  if (message!="")
    Cerr << "-> " << message << finl;
  else
    Cerr << "-> Case not expected, contact pierre.ledac@c-s.fr" << finl;
  Process::exit();
}

Entree& Traitement_particulier_NS_CEG::lire(Entree& is)
{
  haspi_=0;
  C_=125.;
  calculer_critere_areva_=0;
  calculer_critere_cea_jaea_=0;
  nb_mailles_mini_=0;
  critere_cea_jaea_normalise_=0;
  debug_=1;
  dt_post_=-1;
  t_deb_=DMAXFLOAT;
  t_fin_=DMAXFLOAT;
  min_critere_Q_sur_max_critere_Q_=0;
  dernier_temps_=-1e9;
  // Verification gravite
  if ( !mon_equation->milieu().a_gravite() ||
       mon_equation->milieu().gravite().valeurs()(0,0)!=0 ||
       mon_equation->milieu().gravite().valeurs()(0,1)!=0 ||
       mon_equation->milieu().gravite().valeurs()(0,2)>=0 ) error("Error ! Gravity should be defined and oriented parallel to Z, downwards.");
  // Verification VEF 3D
  const DoubleTab& vitesse = mon_equation->inconnue().valeurs();
  if (vitesse.nb_dim()!=2 || vitesse.dimension(1)!=3) error("Error ! Only works in VEF 3D.");

  // XD traitement_particulier_ceg traitement_particulier_base ceg 1  Keyword for a CEG ( Gas Entrainment Criteria)  calculation. An objective is deepening gas entrainment on the free surface. Numerical analysis can be performed to predict the hydraulic and geometric conditions that can   handle gas entrainment from the free surface.

  Param param("lire_ceg");
  param.ajouter("frontiere",&la_surface_libre_nom_,Param::REQUIRED); // XD_ADD_P chaine  To specify the boundaries conditions representing the free surfaces
  param.ajouter("t_deb",&t_deb_,Param::REQUIRED); // XD_ADD_P double value of the CEG's initial calculation time
  param.ajouter("t_fin",&t_fin_); // XD_ADD_P double not_set time during which the CEG's calculation was stopped
  param.ajouter("dt_post",&dt_post_); // XD_ADD_P double periode refers to the printing period, this value is expressed in seconds
  param.ajouter("haspi",&haspi_,Param::REQUIRED); // XD_ADD_P double The suction height  required to calculate AREVA's criterion
  param.ajouter("debug",&debug_); // XD_ADD_P int not_set
  Param& param_areva=param.ajouter_param("AREVA"); // XD_ADD_P ceg_areva AREVA's criterion
  // 2XD ceg_areva objet_lecture nul -1 not_set
  param_areva.ajouter("C",&C_);  // 2XD_ADD_P double not_set
  Param& param_cea_jaea =param.ajouter_param("CEA_JAEA"); // XD_ADD_P ceg_cea_jaea CEA_JAEA's criterion
  // 2XD ceg_cea_jaea objet_lecture nul -1 not_set
  param_cea_jaea.ajouter("normalise",&critere_cea_jaea_normalise_);   // 2XD_ADD_P int renormalize (1) or not (0) values alpha and gamma
  param_cea_jaea.ajouter("nb_mailles_mini",&nb_mailles_mini_);  // 2XD_ADD_P int Sets the minimum number of cells for the detection of a vortex.
  param_cea_jaea.ajouter("min_critere_Q_sur_max_critere_Q",&min_critere_Q_sur_max_critere_Q_);  // 2XD_ADD_P double Is an optional keyword used to correct the minimum values of Q's  criterion taken into account in the detection of a vortex
  param.lire_avec_accolades_depuis(is);


  calculer_critere_areva_=(param.get_list_mots_lus().rang("AREVA")!=-1);
  calculer_critere_cea_jaea_=(param.get_list_mots_lus().rang("CEA_JAEA")!=-1);
  Motcle mot;
  is >> mot;
  if (mot!="}") error("On attendait une }");
  Cerr << "haspi=" << haspi_ << finl;
  Cerr << "calculer_critere_areva="<< calculer_critere_areva_ << finl;
  Cerr << "calculer_critere_cea_jaea="<<calculer_critere_cea_jaea_<<finl;
  Cerr << "nb_mailles_mini=" << nb_mailles_mini_ << finl;
  Cerr << "t_deb_=" << t_deb_ << finl;
  Cerr << "t_fin_=" << t_deb_ << finl;
  Cerr << "min_critere_Q_sur_max_critere_Q_=" << min_critere_Q_sur_max_critere_Q_ << finl;
  if (min_critere_Q_sur_max_critere_Q_>0) error("Error ! min_critere_Q_sur_max_critere_Q_ should less or equal to 0.");
  return is;
}


void Traitement_particulier_NS_CEG::preparer_calcul_particulier()
{
  // Recherche de la frontiere surface libre:
  int trouve=0;
  const Conds_lim& les_cls=mon_equation->domaine_Cl_dis().les_conditions_limites();
  for (int num_cl=0; num_cl<les_cls.size(); num_cl++)
    {
      // Surface libre trouvee
      if (les_cls[num_cl]->frontiere_dis().le_nom()==la_surface_libre_nom_)
        {
          const Domaine_VF& domaine_VF = ref_cast(Domaine_VF, mon_equation->inconnue().domaine_dis_base());
          la_surface_libre_ = ref_cast(Front_VF,les_cls[num_cl]->frontiere_dis());
          trouve=1;
          int nb_faces=la_surface_libre_->nb_faces();
          for (int ind_face=0; ind_face<nb_faces; ind_face++)
            {
              int face = la_surface_libre_->num_face(ind_face);
              if (domaine_VF.face_normales(face,0)!=0 || domaine_VF.face_normales(face,1)) error("The free surface should be normal to Z on all faces.");
            }
        }
    }
  if (!trouve)
    {
      Cerr << "The free surface named " << la_surface_libre_nom_ << " has not been found." << finl;
      exit();
    }
  // KEps si critere AREVA
  if (calculer_critere_areva_)
    if (!sub_type(Pb_Hydraulique_Turbulent,mon_equation->probleme()) || !sub_type(Modele_turbulence_hyd_RANS_K_Eps_base,ref_cast(Navier_Stokes_Turbulent,ref_cast(Pb_Hydraulique_Turbulent,mon_equation->probleme()).equation(0)).modele_turbulence()))
      error("AREVA criterion can only be calculated with a RANS K-eps simulation.");

  // Vorticite dans le jeu de donnees si AREVA
  if (calculer_critere_areva_)
    mon_equation->creer_champ("vorticite");
  // && !ref_cast(Navier_Stokes_std,mon_equation.valeur()).vorticite().non_nul()) error("The vorticity should be postreated in the datafile for the AREVA criterion.");

  // CritereQ dans le jeu de donnees si CEA_JAEA
  if (calculer_critere_cea_jaea_ )
    mon_equation->creer_champ("critere_Q");
  //&& !ref_cast(Navier_Stokes_std,mon_equation.valeur()).critereQ().non_nul()) error("The Q criterion should be postreated in the datafile for the CEA/JAEA criterion.");
}

void Traitement_particulier_NS_CEG::post_traitement_particulier()
{
  if (mon_equation->probleme().schema_temps().temps_courant()>=t_deb_ && mon_equation->probleme().schema_temps().temps_courant()<t_fin_)
    {
      Cerr << "Beginning calculation of the criteria..." << finl;
      statistiques().begin_count(m1_counter_);
      // Calcul des 2 criteres
      if (calculer_critere_areva_) critere_areva();
      statistiques().end_count(m1_counter_);
      if (debug_) Cout << "CPU AREVA criterion " << statistiques().last_time(m1_counter_) << " s" << finl << finl;
      if (calculer_critere_cea_jaea_) critere_cea_jaea();
      dernier_temps_ = mon_equation->probleme().schema_temps().temps_courant();
      Cerr << "End of the calculation of the criteria." << finl;
    }
}

void Traitement_particulier_NS_CEG::critere_areva()
{
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF, mon_equation->inconnue().domaine_dis_base());
  const DoubleTab& vitesse = mon_equation->inconnue().valeurs();
  const DoubleTab& vorticite = mon_equation->get_champ("vorticite").valeurs();
  const Navier_Stokes_Turbulent& eqn = ref_cast(Navier_Stokes_Turbulent,ref_cast(Pb_Hydraulique_Turbulent,mon_equation->probleme()).equation(0));
  const DoubleTab& KEps = ref_cast(Modele_turbulence_hyd_RANS_K_Eps_base,eqn.modele_turbulence()).equation_k_eps(0).inconnue().valeurs();

  double gz = mon_equation->milieu().gravite().valeurs()(0,2);
  double K_max_local=0;
  ArrOfDouble centre_vortex(3);
  int nb_face_par_elem = domaine_VF.elem_faces().dimension(1);
  int nb_faces=la_surface_libre_->nb_faces();

  // On boucle sur les faces reelles
  for (int ind_face=0; ind_face<nb_faces; ind_face++)
    {
      int face = la_surface_libre_->num_face(ind_face);
      // On recupere chaque maille adjacente
      int elem = domaine_VF.face_voisins(face,0);
      if (elem<0) elem = domaine_VF.face_voisins(face,1);

      // Interpolation des champs aux elements:
      double vz=0;
      double k=0;
      double epsilon=0;
      for (int i=0; i<nb_face_par_elem; i++)
        {
          face = domaine_VF.elem_faces(elem,i);
          vz+=vitesse(face,2);
          k+=KEps(face,0);
          epsilon+=KEps(face,1);
        }
      vz*=0.25;
      k*=0.25;
      epsilon*=0.25;

      // Calcul de K
      double wz = std::fabs(vorticite(elem,2));
      double omega;
      if (epsilon>0)
        omega = k * wz / epsilon;			// Rotation de l'ecoulement
      else
        omega = 0;
      double lambda = std::max(-vz,0.) / sqrt(-gz*haspi_); 	// Aspiration de l'ecoulement
      double K = C_ * omega * lambda;			// C_ pour une normalisation

      // On cherche le vortex (std::max(K) et son centre)
      if (K>K_max_local)
        {
          K_max_local=K;
          for (int i=0; i<3; i++) centre_vortex[i] = domaine_VF.xp(elem,i);
        }
    }
  // Parallelisme : Impression par le processeur le plus grand qui possede le vortex
  double K_max_global=mp_max(K_max_local);
  int pe=(K_max_global==K_max_local?Process::me():-1);
  int pe_max=(int)mp_max(pe);
  if (pe==pe_max) imprimer(K_max_global, "AREVA", centre_vortex, 0);
  Cerr << "-> AREVA criterion   : Biggest vortex at (x,y,z)=(" << centre_vortex[0] << "," << centre_vortex[1] << "," << centre_vortex[2] << ")" << finl;
}

int Traitement_particulier_NS_CEG::lpost(double temps_courant, double dt_post) const
{
  double epsilon = 1.e-8;
  if (dt_post<=temps_courant - dernier_temps_)
    return 1;
  else
    {
      // Voir Schema_Temps_base::limpr pour information sur epsilon et modf
      double i, j;
      modf(temps_courant/dt_post + epsilon, &i);
      modf(dernier_temps_/dt_post + epsilon, &j);
      return ( i>j );
    }
}

void Traitement_particulier_NS_CEG::critere_cea_jaea()
{
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF, mon_equation->inconnue().domaine_dis_base());
  int nb_faces=la_surface_libre_->nb_faces();
  int nb_elem=domaine_VF.nb_elem();

  IntVect elements_surface_libre(nb_faces);	// Tableau donnant les elements au contact des faces de la surface libre
  elements_surface_libre=-1;
  IntVect vortex_potentiel(nb_elem); // Tableau donnant pour chaque element si le vortex a ete evalue
  vortex_potentiel=0;
  // On calcule le critereQ

  // GF on laisse faire NS...
  // DoubleTab critereQ(nb_elem);
  //ref_cast(Champ_P1NC,mon_equation->inconnue()).calcul_critere_Q(critereQ);
  const DoubleTab& critereQ = mon_equation->get_champ("critere_Q").valeurs();


  // On boucle sur les faces reelles
  // ou Q>0 (domaine potentielle de vortex)
  for (int ind_face=0; ind_face<nb_faces; ind_face++)
    {
      int face = la_surface_libre_->num_face(ind_face);
      // On recupere chaque maille adjacente
      int elem = domaine_VF.face_voisins(face,0);
      if (elem<0) elem = domaine_VF.face_voisins(face,1);
      // On marque:
      if (critereQ(elem)>0)
        {
          elements_surface_libre(ind_face)=elem;
          vortex_potentiel(elem)=1;
        }
    }
  // Dimensionnement des tableaux recueillant les infos sur les vortexes
  int nb_vortex=0;
  double R0=0;
  double Pi = 4.*atan(1);
  DoubleTab Critere(nb_vortex,3);
  DoubleTab Centre(nb_vortex,3);
  ArrOfDouble Rayon(nb_vortex);
  ArrOfDouble Max_Critere_Q(nb_vortex);
  ArrOfInt Taille(nb_vortex);
  ArrOfDouble centre_vortex(3);  // Tableau de travail

  int nb_dtheta=36; // Angle dtheta=10deg
  double dtheta=2*Pi/nb_dtheta;
  ArrOfInt elements(nb_dtheta);
  DoubleTab points(nb_dtheta,3);
  DoubleTab u(nb_dtheta,3);
  int taille_maxi = 0;
  // Boucle tant que la liste d'elements n'est pas vide(tous les vortexes seront alors evalues)
  while (elements_surface_libre.mp_max_vect()>=0)
    {
      statistiques().begin_count(m1_counter_);
      double critereQ_max=0;
      int ind_face_centre_vortex=-1;
      // On boucle sur les elements non evalues pour trouver le plus grand critere Q
      for (int ind_face=0; ind_face<nb_faces; ind_face++)
        {
          int elem=elements_surface_libre(ind_face);
          if (elem>=0 && critereQ(elem)>critereQ_max)
            {
              critereQ_max=critereQ(elem);
              ind_face_centre_vortex=ind_face;
            }
        }
      // Plus grand Q sur l'ensemble des processeurs:
      double critereQ_mp_max=mp_max(critereQ_max);
      // On traite le cas ou meme Q (on prend alors le processeur de rang le plus eleve)
      int pe=(critereQ_mp_max==critereQ_max?Process::me():-1);
      int pe_mp_max=(int)mp_max(pe);
      int taille_vortex=0;
      // Le centre du vortex est le centre de elem_centre_vortex
      // on le communique a tous les processeurs
      if (pe==pe_mp_max)
        {
          int face_centre_fortex = la_surface_libre_->num_face(ind_face_centre_vortex);
          int elem_centre_vortex = elements_surface_libre(ind_face_centre_vortex);
          vortex_potentiel(elem_centre_vortex)=0; // Ne sera plus evalue plus tard
          taille_vortex=1;

          // On calcule le rayon initial (rayon du cercle inscrit) dans le triangle:
          // https://fr.wikipedia.org/wiki/Cercles_inscrit_et_exinscrits_d%27un_triangle
          // Rayon d'un cercle inscrit dans un triangle: R=2*Surface(triangle)/perimetre(triangle)
          // Calcul du perimetre
          double surface=domaine_VF.face_surfaces(face_centre_fortex);
          const DoubleTab& coord=mon_equation->probleme().domaine().coord_sommets();
          int nb_sommet_face=domaine_VF.face_sommets().dimension(1);
          double perimetre=0;
          for (int i=0; i<nb_sommet_face; i++)
            {
              int som1 = domaine_VF.face_sommets(face_centre_fortex,i);
              int som2 = (i!=nb_sommet_face-1 ? domaine_VF.face_sommets(face_centre_fortex,i+1) : domaine_VF.face_sommets(face_centre_fortex,0));
              double dx=coord(som2,0)-coord(som1,0);
              double dy=coord(som2,1)-coord(som1,1);
              perimetre+=sqrt(dx*dx+dy*dy);
            }
          // On multiplie par 3/4 par homothetie
          R0 = 0.6*0.75*(2*surface/perimetre);

          // Centre du vortex
          for (int i=0; i<3; i++) centre_vortex[i]=domaine_VF.xp(elem_centre_vortex,i);
          // Envoi des donnees tous les autres processes:
          for (int p=0; p<nproc(); p++)
            if (p!=pe_mp_max)
              {
                envoyer(centre_vortex,pe_mp_max,p,p);
                envoyer(R0,pe_mp_max,p,p);
              }
        }
      else
        {
          // Receptions des donnees
          recevoir(centre_vortex,pe_mp_max,me(),me());
          recevoir(R0,pe_mp_max,me(),me());
        }
      /**************************************************************/
      // Calcul des criteres autour du vortex de centre centre_vortex
      /**************************************************************/
      // On va cherche le cercle le plus large pour lequel Q>0 et
      // ensuite on calcule les circulations
      // Processus de dichotomie sur R
      statistiques().end_count(m1_counter_);
      double R=R0;
      double dR=R;
      int niter=0;
      int limite_vortex_atteinte=0;
      int points_trouves=0;
      statistiques().begin_count(m2_counter_);
      while (dR>0.01*R0)
        {
          double inside_vortex=1;
          points_trouves=0;
          for (int i_theta=0; i_theta<nb_dtheta; i_theta++)
            {
              double theta=i_theta*dtheta;
              points(i_theta,0)=centre_vortex[0]+R*cos(theta);
              points(i_theta,1)=centre_vortex[1]+R*sin(theta);
              points(i_theta,2)=centre_vortex[2];
            }
          // On cherche les elements contenant le point x,y,z:
          mon_equation->domaine_dis().domaine().chercher_elements(points,elements);
          for (int i_theta=0; i_theta<nb_dtheta; i_theta++)
            {
              int elem=elements[i_theta];
              // Test pour verifier que le cercle R0 est bien inscrit au tetraedre:
              if (niter==0 && pe==pe_mp_max)
                {
                  int elem_centre_vortex=elements_surface_libre(ind_face_centre_vortex);
                  if (elem!=elem_centre_vortex)
                    error("Error, R0 is too large. Algorithm problem, contact pierre.ledac@c-s.fr");
                }
              if (elem>=0 && elem<nb_elem)
                {
                  // Point trouve dans le domaine reel
                  points_trouves++;
                  if (critereQ(elem)<=min_critere_Q_sur_max_critere_Q_*critereQ_mp_max)
                    {
                      // On est en dehors du tourbillon
                      inside_vortex=0;
                      break;
                    }
                }
              //Cerr << "elem=" << elem << " Q=" << critereQ(elem) << finl;
            }
          inside_vortex=mp_min(inside_vortex);
          if (inside_vortex)
            {
              if (limite_vortex_atteinte) dR*=0.5;
              //if (debug_) Cerr << "[" << niter << "] R=" << R << " dR=" << dR << finl;
              R+=dR; // On continue d'avancer
              // On supprime les elements du vortex actuel de la liste
              for (int i_theta=0; i_theta<nb_dtheta; i_theta++)
                {
                  int elem=elements[i_theta];
                  if (elem>=0 && elem<nb_elem)
                    {
                      if (vortex_potentiel(elem)!=0) taille_vortex++;
                      vortex_potentiel(elem)=0; // Ne seront plus evalues plus tard
                    }
                }
            }
          else
            {
              limite_vortex_atteinte=1;
              dR*=0.5;
              //if (debug_) Cerr << "[" << niter << "] R=" << R << " dR=" << -dR << finl;
              R-=dR; // On recule
            }
          niter++;
        }
      R-=dR;
      statistiques().end_count(m2_counter_);
      statistiques().begin_count(m3_counter_);
      // On ne prend que des vortex superieres a N mailles dont toutes les
      // mailles sont dans le domaine fluide (points_trouves==nb_dtheta)
      points_trouves=mp_sum(points_trouves);
      int taille_vortex_mailles=mp_sum(taille_vortex);
      taille_maxi = (taille_vortex_mailles > taille_maxi ? taille_vortex_mailles : taille_maxi);
      if (points_trouves==nb_dtheta && taille_vortex_mailles>nb_mailles_mini_)
        {
          //Cerr << printf("-> Critere CEA_JAEA: Vortex de %d mailles en (x,y,z)=(%.4f,%.4f,%.4f) Rayon=%.2f std::max(critere_Q)=%.2f\n",taille_vortex_mailles,centre_vortex(0),centre_vortex(1),centre_vortex(2),R,critereQ_mp_max);
          if (debug_)
            {
              Cerr << "-> CEA_JAEA criterion : Vortex " << nb_vortex << " de " << taille_vortex_mailles << " elements ";
              Cerr << "in (x,y,z)=(" << centre_vortex[0] << "," << centre_vortex[1] << "," << centre_vortex[2] << ") ";
              Cerr << "Radius=" << R << " std::max(critere_Q)=" << critereQ_mp_max;
              //int e = mon_equation->domaine_dis().domaine().chercher_elements(centre_vortex(0),centre_vortex(1),centre_vortex(2));
              //if (e>=0) Cerr << "Sur process " << Process::me() << " critere_Q=" << critereQ(e) << finl;
            }
          // On evalue la vitesse aux points:
          mon_equation->inconnue().valeur_aux(points,u);
          double alpha=0;
          double gamma=0;
          // On integre sur le cercle pour calcul alpha et gamma:
          for (int i_theta=0; i_theta<nb_dtheta; i_theta++)
            {
              double theta=i_theta*dtheta;
              int elem=elements[i_theta];
              if (elem>=0 && elem<nb_elem)
                {
                  alpha+=(u(i_theta,0)*cos(theta)+u(i_theta,1)*sin(theta))*R*dtheta;		// Integration vitesse normale (dans le plan)
                  gamma+=(u(i_theta,1)*cos(theta)-u(i_theta,0)*sin(theta))*R*dtheta; 	// Integration vitesse tangentielle
                }
            }

          alpha=mp_sum(alpha)/(Pi*R*R);
          gamma=mp_sum(gamma);
          // On adimensionnalise eventuellement selon jeu de donnees (conseille):
          if (critere_cea_jaea_normalise_)
            {
              double nu = ref_cast(Fluide_Incompressible,mon_equation->milieu()).viscosite_cinematique()->valeurs()(0,0);
              double gz = mon_equation->milieu().gravite().valeurs()(0,2);
              alpha*=nu/(gz*haspi_);  	// alpha*=alpha*nu/(gh)
              gamma/=nu;	    		// gamma*=gamma/nu
            }
          // On redimensionne:
          Critere.resize(nb_vortex+1,3);
          Centre.resize(nb_vortex+1,3);
          Rayon.resize(nb_vortex+1);
          Max_Critere_Q.resize(nb_vortex+1);
          Taille.resize(nb_vortex+1);
          // On remplit
          Critere(nb_vortex,0)=alpha;				// alpha
          Critere(nb_vortex,1)=Critere(nb_vortex,0)*gamma*gamma;	// alpha*gamma^2
          Critere(nb_vortex,2)=alpha*Critere(nb_vortex,1);		// (alpha*gamma)^2
          Rayon[nb_vortex]=R;
          Max_Critere_Q[nb_vortex]=critereQ_mp_max;
          Taille[nb_vortex]=taille_vortex_mailles;
          for (int i=0; i<3; i++) Centre(nb_vortex,i)=centre_vortex[i];
          if (debug_) Cerr << " alpha=" << Critere(nb_vortex,0) << " alpha*gamma^2=" << Critere(nb_vortex,1) << " (alpha*gamma)^2=" << Critere(nb_vortex,2) << finl;
          R+=dR;
          nb_vortex++;
        }
      // Mise a jour de elements_surface_libre par rapport a vortex_potentiel:
      for (int ind_face=0; ind_face<nb_faces; ind_face++)
        {
          int elem=elements_surface_libre(ind_face);
          if (elem>=0 && vortex_potentiel(elem)==0) elements_surface_libre(ind_face)=-1;
        }
    }
  statistiques().end_count(m3_counter_);
  if (debug_) Cout << "CPU AREVA criterion m1= " << statistiques().last_time(m1_counter_) << " s m2= " << statistiques().last_time(m2_counter_) << " s m3= " << statistiques().last_time(m3_counter_) << finl;
  const Schema_Temps_base& sch=mon_equation->probleme().schema_temps();
  if (nb_vortex==0)
    {
      Cerr << "-> CEA_JAEA criterion : No vortex detected of size larger than " << nb_mailles_mini_ << " elements. A potential vortex was " << taille_maxi << " elements large." << finl;
    }
  else
    {
      Cerr << "-> CEA_JAEA criterion : we found " << nb_vortex << " vortices near the free surface." << finl;
      // Recherche du vortex donnant la valeur maximale pour chaque critere:
      ArrOfDouble Critere_max(3);
      ArrOfInt Vortex_max(3);
      Vortex_max=-1;
      for (int critere=0; critere<3; critere++)
        for (int vortex=0; vortex<nb_vortex; vortex++)
          if (Critere(vortex,critere)>Critere_max[critere])
            {
              Critere_max[critere]=Critere(vortex,critere);
              Vortex_max[critere]=vortex;
            }

      // On prend comme vortex principal celui indique par alpha*gamma^2:
      int vortex_principal = Vortex_max[1];
      if (vortex_principal>=0)
        {
          if (Process::je_suis_maitre())
            {
              const Postraitement& post = ref_cast(Postraitement,mon_equation->probleme().postraitements()[0].valeur());
              // Si dt_post pas lu dans le jeu de donnees, on prend celui du postraitement
              if (dt_post_<0) dt_post_ = post.dt_post_ch();
              // Impression des criteres
              Noms Critere_nom(3);
              Critere_nom[0]="CEA_JAEA_alpha";
              Critere_nom[1]="CEA_JAEA_alphaXgamma2";
              Critere_nom[2]="CEA_JAEA_alpha2Xgamma2";

              for (int critere=0; critere<3; critere++)
                {
                  //int vortex = Vortex_max(critere);
                  for (int i=0; i<3; i++) centre_vortex[i]=Centre(vortex_principal,i);
                  //imprimer(Critere_max(critere), Critere_nom(critere), centre_vortex, Rayon(vortex));
                  imprimer(Critere(vortex_principal,critere), Critere_nom[critere], centre_vortex, Rayon[vortex_principal]);

                  // Impression periodique des caracteristiques du vortex pour chaque critere
                  if (lpost(sch.temps_courant(),dt_post_))
                    {
                      Nom filename("Vortex_");
                      filename+=Critere_nom[critere];
                      filename+="_";
                      filename+=Objet_U::nom_du_cas();
                      filename+=".";
                      filename+=(Nom)(sch.temps_courant());
                      filename+=".csv";
                      SFichier file(filename);
                      file << "# " << filename << finl;
                      file << "# Temps";
                      for (int i_theta=0; i_theta<nb_dtheta; i_theta++)
                        {
                          double theta=i_theta*dtheta;
                          double x=centre_vortex[0]+Rayon[vortex_principal]*cos(theta);
                          double y=centre_vortex[1]+Rayon[vortex_principal]*sin(theta);
                          double z=centre_vortex[2]+2*R0;
                          file << " x= " << x << " y= " << y << " z= " << z;
                        }
                      file << finl;
                      file << "# Champ CRITERE" << finl;
                      file << "# Type POINTS" << finl;
                      file << sch.temps_courant();
                      //for (int i_theta=0;i_theta<nb_dtheta;i_theta++) file << " " << Critere_max(critere);
                      for (int i_theta=0; i_theta<nb_dtheta; i_theta++) file << " " << Critere(vortex_principal,critere);
                      file << finl;
                    }
                }
              // Impression periodique des centres de tous les vortex trouves
              if (lpost(sch.temps_courant(),dt_post_))
                {
                  Nom filename("Centres_vortex.");
                  filename+=(Nom)(sch.temps_courant());
                  filename+=".csv";
                  SFichier file(filename);
                  file << "# " << filename << finl;
                  file << "# Temps";
                  for (int vortex=0; vortex<nb_vortex; vortex++) file << " x= " << Centre(vortex,0) << " y= " << Centre(vortex,1) << " z= " << Centre(vortex,2)+2*R0;
                  file << finl;
                  file << "# Champ RAYON" << finl;
                  file << "# Type POINTS" << finl;
                  file << sch.temps_courant();
                  for (int vortex=0; vortex<nb_vortex; vortex++) file << " " << Rayon[vortex];
                  file << finl;
                }
            }
          Cerr << "-> CEA_JAEA criterion : The main vortex (largest alpha*gamma^2) was vortex number " << vortex_principal << finl;
        }
    }
}

void Traitement_particulier_NS_CEG::imprimer(const double valeur_critere, const Nom& critere, const ArrOfDouble& centre_vortex, const double rayon_vortex)
{
  Nom filename=Objet_U::nom_du_cas();
  filename+="_";
  filename+=mon_equation->probleme().le_nom()+"_";
  filename+=critere;
  filename+=".csv";
  SFichier fic;
  struct stat f;
  const Schema_Temps_base& sch=mon_equation->probleme().schema_temps();
  if (stat(filename,&f) || (sch.nb_impr()==0 && !mon_equation->probleme().reprise_effectuee()))
    {
      fic.ouvrir(filename);
      // Ecriture en tete
      fic << "# Time \tVortex center (x,y,z) \t\t\tVortex radius \t" << critere << " criterion" << finl;
    }
  else
    fic.ouvrir(filename,ios::app);
  fic.precision(sch.precision_impr());
  fic.setf(ios::scientific);
  // Ecriture
  fic << sch.temps_courant() << " \t" << centre_vortex[0] << " \t" << centre_vortex[1] << " \t" << centre_vortex[2] << " \t" << rayon_vortex << " \t" << valeur_critere << finl;
  fic.close();
}
