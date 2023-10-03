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
// File      : Trait_part_NS_plan.cpp
// Directory : $NEW_ALGO_QC_ROOT/src/modif_QC/Trait_part_NS_plan
//
/////////////////////////////////////////////////////////////////////////////

#include <Trait_part_NS_plan.h>
// Lectures ecriture dans les fichiers
#include <LecFicDistribue.h>
#include <EcrFicCollecte.h>
// coordonees des points
#include <Zone_VF.h>
#include <Zone_VDF.h>
// types specifiques de tableaux.
#include <DoubleTrav.h>
// communications mpi
#include <communications.h>
#include <Schema_Comm.h>
// Equations ou donnees necessaires au post traitement
#include <Probleme_base.h>
#include <Loi_Etat_GP.h>
#include <Convection_Diffusion_Chaleur_Turbulent_QC.h>
#include <Schema_Temps_base.h>  // Utile pour les pas de temps (doublons avec Schema_Temps.h) 
#include <Schema_Temps.h> // Utile pour les pas de temps
#include <Navier_Stokes_Turbulent.h>
#include <Fluide_Quasi_Compressible.h>

#include <Param.h>
#include <IJK_Field.h>
#include <communications.h>

Implemente_base(Traitement_particulier_NS_plan,"Traitement_particulier_NS_plan",Traitement_particulier_NS_base);


void test_convergence(const DoubleTab& Moyennes, const DoubleTab& sauv_moyennes,
                      const double tps, const double temps_de_moyenne,
                      const char *filename_suivi_convergence)
{
  ///parfait ici je vais tous faire
  /// pour faire les test on va utiliser un tableau qui ne va exister que sur le maitre.
  int j = 0 ;
  double var = 0 , MaxE =0 ;
  double usdt = 0 ;
  if(temps_de_moyenne)
    usdt =  1. / (temps_de_moyenne) ;

  for (j = 0 ; j < Moyennes.dimension(0); j++ )
    {
      var =  (Moyennes(j,22) -  Moyennes(j,5) * Moyennes(j,5)) * usdt ;

      if (fabs(var) >= 1e-10)
        {
          var -= sauv_moyennes(j,0,0) -  sauv_moyennes(j,1,0) * sauv_moyennes(j,1,0);
          var /= (Moyennes(j,22) -  Moyennes(j,5) * Moyennes(j,5))*usdt ;
        }
      var = fabs(var) ;
      if ( var > MaxE )
        MaxE = var ;
    }

  SFichier fic(filename_suivi_convergence,ios::app);
  fic.setf(ios::scientific);
  fic << tps << "  " << MaxE << finl;
  fic.close();
}

static void stocker_moy(const DoubleTab& Moyennes, DoubleTab& sauv_moyennes, double temps_de_moyenne)
{

  double usdt =0;
  if(temps_de_moyenne)
    usdt =  1. / (temps_de_moyenne) ;

  const int n = Moyennes.dimension(0);
  for (int j = 0 ; j < n; j++ )
    {
      sauv_moyennes(j,0,0) =  sauv_moyennes(j,0,1) ;
      sauv_moyennes(j,1,0) =  sauv_moyennes(j,1,1) ;
      sauv_moyennes(j,0,1) =  sauv_moyennes(j,0,2) ;
      sauv_moyennes(j,1,1) =  sauv_moyennes(j,1,2) ;
      sauv_moyennes(j,0,2) =  Moyennes(j,22)* usdt ;
      sauv_moyennes(j,1,2) =  Moyennes(j,5 )* usdt ;
    }
}

static void ecriture_fichiers_moyennes(const IJK_Splitting& split,
                                       const DoubleTab& Moyennes, const Nom& fichier,
                                       double temps_courant, double temps_de_moyenne)
{
  double tps = temps_courant;

  int i,j;
  SFichier fic (fichier) ;
  fic.setf(ios::scientific);
  const int ni = split.get_grid_geometry().get_nb_elem_tot(DIRECTION_I);
  const int nj = split.get_grid_geometry().get_nb_elem_tot(DIRECTION_J);
  const int nk = split.get_grid_geometry().get_nb_elem_tot(DIRECTION_K);


  fic << "# Temps : "  <<  tps << " Nombre de points :" << ni * nj  << finl;
  fic << "# (1) y " << finl;
  fic << "# (2) = <U1>           (13) = U1rms " << finl;
  fic << "# (3) = <U2>           (14) = U2rms " << finl;
  fic << "# (4) = <U3>           (15) = U3rms " << finl;
  fic << "# (5) = <V1>           (16) = V1rms "  << finl;
  fic << "# (6) = <V2>           (17) = V2rms "  << finl;
  fic << "# (7) = <V3>           (18) = V3rms "  << finl;
  fic << "# (8) = <rho>          (19) = <u1'v2'> " << finl;
  fic << "# (9) = <T>            (20) = Trms " << finl;
  fic << "# (10) = <C>           (21) = Crms " << finl;
  fic << "# (11) = <mu>        (22) = <u1'.T'> " << finl;
  fic << "# (12) = <lambdaeff>   " << finl;
  if (temps_de_moyenne == 0)
    temps_de_moyenne = 1 ; //si on post traite en espace on evite de diviser par zero.

  const ArrOfDouble& coord_node = split.get_grid_geometry().get_node_coordinates(DIRECTION_K);

  for (j = 0 ; j < nk ; j++ )
    {
      fic << (coord_node[j] + coord_node[j+1]) * 0.5;
      for (i = 0 ; i <6  ; i++ )
        fic << "  " <<  Moyennes(j,i)/temps_de_moyenne ;
      fic << "  " <<  Moyennes(j,35)/temps_de_moyenne ;
      fic << "  " <<  Moyennes(j,7)/temps_de_moyenne ;
      fic << "  " <<  Moyennes(j,8)/temps_de_moyenne ;
      fic << "  " <<  Moyennes(j,9)/temps_de_moyenne ;
      fic << "  " <<  Moyennes(j,34)/temps_de_moyenne ;
      for (i = 17 ; i <23  ; i++ )
        fic << "  " <<  Moyennes(j,i)/temps_de_moyenne ;
      fic << "  " <<  Moyennes(j,16)/temps_de_moyenne ;
      fic << "  " <<  Moyennes(j,24)/temps_de_moyenne ;
      fic << "  " <<  Moyennes(j,25)/temps_de_moyenne ;
      fic << "  " <<  Moyennes(j,36)/temps_de_moyenne << finl ;
    }

  fic.flush();
  fic.close();
}

// Fonction qui fait la reprise a partir d'un fichier de type Stat_post_traitement_temps
static void reprendre_stat_plan(const IJK_Splitting& splitting,
                                const Nom& fichier,
                                DoubleTab& Moyennes,
                                DoubleTab& sauv_moyennes,
                                double& temps_de_moyenne)
{
  if (!Process::je_suis_maitre())
    {
      Cerr << "reprendre_stat_plan ne doit etre appele que sur le maitre" << finl;
      Process::exit();
    }
  Cerr << "Traitement_particulier_NS_canal::reprendre_stat" << flush << finl;
  int i,j,k ;
  double temps_courant, y;
  EFichier fic(fichier); // constructor exits if open fails
  fic.setf(ios::scientific);
  fic >> temps_courant; // valeur du temps ou la sauvegarde a ete faite
  fic >> temps_de_moyenne; // valeur du temps de moyennes

  const int nj_tot = splitting.get_grid_geometry().get_nb_elem_tot(DIRECTION_K);
  const ArrOfDouble& coord_node = splitting.get_grid_geometry().get_node_coordinates(DIRECTION_K);
  for (j = 0 ; j < nj_tot ; j++ )
    {
      const double coord_y = (coord_node[j] + coord_node[j+1]) * 0.5;
      fic >> y;
      if (fabs(y - coord_y) > 1e-6)
        {
          Cerr << "Erreur dans reprendre_stat_plan: y["<<i<<"] ne colle pas: "
               << y << " != " << coord_y << finl;
          Process::exit();
        }
      for (i = 0 ; i < Moyennes.dimension(1); i++ )
        fic >> Moyennes(j,i); // on recupere le tableau.
    }
  for (j = 0 ; j < nj_tot ; j++ )
    for (i = 0 ; i < sauv_moyennes.dimension(1); i++ )
      for (k = 0 ; k < sauv_moyennes.dimension(2); k++ )
        fic >> sauv_moyennes(j,i,k);

}


static void sauver_stat_plan(const Nom& fichier,
                             const ArrOfDouble& coord_y,  // coordonnees Y des mailles
                             const DoubleTab& Moyennes,
                             const DoubleTab& sauv_moyennes,
                             const double tps,
                             const double temps_de_moyenne)
{
  if (!Process::je_suis_maitre())
    {
      Cerr << "reprendre_stat_plan ne doit etre appele que sur le maitre" << finl;
      Process::exit();
    }
  int i,j,k;
  SFichier fic2(fichier); // constructor exits if open fails
  fic2.setf(ios::scientific);

  fic2 << tps << finl;
  fic2 << temps_de_moyenne << finl;
  const int nj_tot = coord_y.size_array();
  for (j = 0 ; j < nj_tot ; j++ )
    {
      fic2 << coord_y(j) ;
      for (i = 0 ; i < Moyennes.dimension(1)	; i++ )
        fic2 << "  " <<  Moyennes(j,i) ;

      fic2 << finl;
    }
  for (j = 0 ; j < nj_tot ; j++ )
    {
      for (i = 0 ; i < sauv_moyennes.dimension(1)	; i++ )
        for (k = 0 ; k < sauv_moyennes.dimension(2)	; k++ )
          fic2 << "  " <<   sauv_moyennes(j,i,k) ;
      fic2 << finl;
    }
  fic2.close();

}

// on va post traiter les grandeurs parrietales calculer a partir des grandeurs spatiales (la calcule a partir de la moyennes globale interviendra plus tard.
static void ecrire_reynolds_tau(const DoubleTab& moyennes, const double temps_courant,
                                const ArrOfDouble& Y_tot, const char *nom_fichier)
{
  double tps = temps_courant; //

  int Ny = Y_tot.size_array();
  double h;      // demi-hauteur
  if (Ny%2) //si Ny est impaire ;
    h = Y_tot(Ny/2); // alors h est la position de l'element central
  else /// si NY est paire alors la demie hauteur est la coordonnees de face entre les deux mailles centrales
    h = (Y_tot(Ny/2) + Y_tot(Ny/2-1))/2. ;

  // definition des grandeurs interveannt dans la suite des calculs
  //////////////////////////////////////////////////////
  double y1,		y2				;
  double mu_b, 	mu_h			;    // viscosite dynamique
  double rho_b, 	rho_h			;	 // masse volumique
  double U1,		U2				;	 // norme de la vitesse tangente a la paroi
  double tauwb,	tauwh, 	tauwm	;	 // cisaillement a la paroi
  double utaub, 	utauh, 	utaum	;	 // vitesse de frottement
  double retaub, 	retauh, retaum	;	 // Reynolds de frottement

  mu_b  	= moyennes(0   ,9) ;
  mu_h 	= moyennes(moyennes.dimension(0)-1,9) ;

  rho_b  	= moyennes(0   ,35);
  rho_h 	= moyennes(moyennes.dimension(0)-1,35);

  U1  	= moyennes(0    ,0) ;
  U2	 	= moyennes(1    ,0) ;

  y1		=Y_tot(0);
  y2		=Y_tot(1);

  // calcul du cisaillement a la paroi suivant la vitesse tangentielle avec non glissement a la paroi.
  //  U0 | U1 | U2 |
  //  0  | y1 | y2 |
  // approximation d'ordre 2 : Tau_w = (U1*y2/y1-U2*y1/y2)/(y2-y1))
  // calcul et ecritures des differentes grandeurs parietales
  //////////////////////////////////////////////////////
  tauwb= mu_b * (U1*(y2/y1)-U2*(y1/y2)) /(y2-y1);
  //tauwb= mu_b * (fabs(U1)/y1);
  U1  	= moyennes(moyennes.dimension(0)-1,0);
  U2	 	= moyennes(moyennes.dimension(0)-2,0);

  tauwh= mu_h * (U1*(y2/y1)-U2*(y1/y2)) /(y2-y1);
  //tauwh= mu_h * (fabs(U1)/y1);
  tauwh= fabs(tauwh);
  utaub =sqrt(tauwb/rho_b);
  utauh =sqrt(tauwh/rho_h);

  retaub=rho_b*utaub*h/mu_b;
  retauh=rho_h*utauh*h/mu_h;

  tauwm =(tauwh+tauwb)/2.;
  retaum=(retauh+retaub)/2.;
  utaum =(utauh+utaub)/2.;

  SFichier fic1("Parietal.dat",ios::app);

  fic1 << tps   << " " << retaum << " " << retaub << " " << retauh << " ";
  fic1 << utaum << " " << utaub  << " " << utauh  << " " ;
  fic1 << tauwm << " " << tauwb  << " " << tauwh  << finl;

  fic1.close();

}



Sortie& Traitement_particulier_NS_plan::printOn(Sortie& is) const
{
  return is;
}

Entree& Traitement_particulier_NS_plan::readOn(Entree& is)
{
  return is;
}

Entree& Traitement_particulier_NS_plan::lire(Entree& is)
{
  start_time_ = 1e6; // default: never
  dt_write_spatial_avg_ = 1e6;
  fich_repr_ = "??";
  oui_calcul_ = false;

  Param param(que_suis_je());
  param.ajouter("dt_impr", &dt_write_spatial_avg_);
  param.ajouter("debut", &start_time_);
  param.ajouter("plan", &indices_planes_);
  param.ajouter("reprise", &fich_repr_);
  param.ajouter_flag("converger", &oui_calcul_);
  param.lire_avec_accolades(is);
  Motcle test ;
  is >> test ;
  if (oui_calcul_ && fich_repr_=="??")
    {
      Cerr << "Erreur Traitement_particulier_NS_plan: il faut une reprise pour l'option converger"<< endl;
      exit();
    }
  oui_repr_ = (fich_repr_ != "??");
  return is;
}

void Traitement_particulier_NS_plan::preparer_calcul_particulier()
{
  const Zone_dis_base& zdisbase = mon_equation->inconnue().zone_dis_base();
  const Zone_VDF& zone_VDF = ref_cast(Zone_VDF, zdisbase);
  IJK_Grid_Geometry ijk_grid;
  // swap y and z in conversion:
  ijk_grid.initialize_from_unstructured(zone_VDF.zone(), DIRECTION_I, DIRECTION_K, DIRECTION_J,
                                        true, true, false /* not periodic in k */);

  // Create a splitting of the ijk mesh such that whole z planes are on each processor
  int nslices_k = nproc();
  if (nslices_k > ijk_grid.get_nb_elem_tot(DIRECTION_K))
    nslices_k = ijk_grid.get_nb_elem_tot(DIRECTION_K);
  ijk_splitting_.initialize(ijk_grid,
                            1 /* 1 slice in direction i */,
                            1 /* 1 slice in direction j */,
                            nslices_k);

  vdf_to_ijk_.initialize(zone_VDF, ijk_splitting_, IJK_Splitting::ELEM,
                         DIRECTION_I, DIRECTION_K, DIRECTION_J);


  // Preparation des tableaux "Moyennes_temporelles_"
  if (Process::je_suis_maitre())
    {
      const int nk = ijk_grid.get_nb_elem_tot(DIRECTION_K);
      Moyennes_temporelles_.resize(nk, nb_moyennes_temporelles());
      Moyennes_temporelles_ = 0.;
      // on va sauvez 3 moyennes en temps de <W^2> et <W>. il n'y a que le maitre qui en a besoin pour le moment.
      sauv_moyennes_.resize(nk,2,3);
      sauv_moyennes_ = 0.;
    }
  temps_de_moyenne_ = 0;
  if (oui_repr_)   // si on fait une reprise des stats
    {
      if (Process::je_suis_maitre())
        {
          reprendre_stat_plan(ijk_splitting_, fich_repr_, Moyennes_temporelles_, sauv_moyennes_, temps_de_moyenne_ ); // recupere les statistiques
          if (oui_calcul_)
            {
              Moyennes_finale_ = Moyennes_temporelles_;
              Moyennes_finale_ /= temps_de_moyenne_ ; // on calcule la norme,
            }
        }
      envoyer_broadcast(temps_de_moyenne_, 0);
      envoyer_broadcast(Moyennes_finale_, 0);
    }

  ref_champ_temperature_ = mon_equation.valeur().probleme().equation(1).get_champ("temperature_qc");
  const Fluide_Quasi_Compressible& fluide_QC = ref_cast(Fluide_Quasi_Compressible, mon_equation->milieu());
  Pth_passe_ = fluide_QC.pression_th();

  // cette portion de code sera a retirer apres traduction du code en tout ijk:

  zone_VDF.zone().creer_tableau_elements(racinederho);
  zone_VDF.creer_tableau_faces(Rrhouf);

  {
    const IntTab& elem_faces = zone_VDF.elem_faces();

    int face; //recepteur des faces
    int nb_elems = zone_VDF.zone().nb_elem_tot(); // pour pouvoir faire les derivation le tableau doit comporter les adresses des fictifs.
    int num_elem;

    // inutile apres tout ijk.
    Tab_recap.resize(nb_elems,6); // on dimentionne le tableau.
    //Le tableau donne le numero de l'element voisin de la face 0,3,1,4,2,5

    for (num_elem=0; num_elem<nb_elems; num_elem++) //Remplissage du tableau des voisin et des vecteur X et Z
      {
        // on remplit le tableau des voisins
        // elem_voisin(numerode de l'elem, numero de la face, 0 dans la direction oposee du repere ou 1 dans la direction au repere)
        // si le voisin n'existe pas la valeur retourner est -1
        //faces X
        face                   =  elem_faces(num_elem,0);
        Tab_recap(num_elem,0)  =  zone_VDF.elem_voisin(num_elem,face,0);

        face                   =  elem_faces(num_elem,0+dimension);
        Tab_recap(num_elem,1)  =  zone_VDF.elem_voisin(num_elem,face,1);

        //faces Y
        face                   =  elem_faces(num_elem,1); //face inferieure
        Tab_recap(num_elem,2)  =  zone_VDF.elem_voisin(num_elem,face,0);

        face                   =  elem_faces(num_elem,1+dimension); //face superieure
        Tab_recap(num_elem,3)  =  zone_VDF.elem_voisin(num_elem,face,1);

        //face Z
        face                   =  elem_faces(num_elem,2);
        Tab_recap(num_elem,4)  =  zone_VDF.elem_voisin(num_elem,face,0);

        face                   =  elem_faces(num_elem,2+dimension);
        Tab_recap(num_elem,5)  =  zone_VDF.elem_voisin(num_elem,face,1);

      } //fin  for (num_elem=0;num_elem<nb_elems;num_elem++)

  }
  // cette portion de code sera a retirer apres traduction du code en tout ijk:
  {
    IJK_Field_double indice_k;
    indice_k.allocate(ijk_splitting_, IJK_Splitting::ELEM, 0 /* ghost size */);
    const int offset = ijk_splitting_.get_offset_local(DIRECTION_K);
    for (int k = 0; k < indice_k.nk(); k++)
      for (int j = 0; j < indice_k.nj(); j++)
        for (int i = 0; i < indice_k.ni(); i++)
          indice_k(i,j,k) = k + offset;
    DoubleVect ref_y_double;
    mon_equation->inconnue().zone_dis_base().zone().creer_tableau_elements(ref_y_double);
    vdf_to_ijk_.convert_from_ijk(indice_k, ref_y_double);
    mon_equation->inconnue().zone_dis_base().zone().creer_tableau_elements(Ref_Y);
    const int nb_elem = Ref_Y.size();
    for (int i = 0; i < nb_elem; i++)
      {
        Ref_Y[i] = (int) round(ref_y_double[i]);
      }
  }
}

// Cette methode calcule si on doit post-traiter le temps courant ou non, renvoie 1 si oui, 0 sinon.
int Traitement_particulier_NS_plan::calculer_critere_post() const
{
  int post = 0;
  if (je_suis_maitre())
    {
      double tps = mon_equation->inconnue().temps();
      double stationnaire_atteint = mon_equation->schema_temps().stationnaire_atteint();
      double tps_passe		  = mon_equation->schema_temps().temps_courant();
      int nb_pas_dt_max	   	  = mon_equation->schema_temps().nb_pas_dt_max();
      int nb_pas_dt			    = mon_equation->schema_temps().nb_pas_dt();
      double temps_max		    = mon_equation->schema_temps().temps_max();
      int impr_inst;

      if ( dt_write_spatial_avg_<=(tps-tps_passe) )
        {
          impr_inst=1;
        }
      else
        {
          double i, j, epsilon = 1.e-8; // Voir Schema_Temps_base::limpr pour information sur epsilon et modf
          modf(tps/dt_write_spatial_avg_ + epsilon, &i);
          modf(tps_passe/dt_write_spatial_avg_ + epsilon, &j);
          impr_inst=(i>j);
        }
      post = ( (nb_pas_dt+1<=1) || impr_inst || (temps_max <= tps) || (nb_pas_dt_max <= nb_pas_dt+1) || stationnaire_atteint  );
    }
  // Pour s'assurer que le resultat du test est identique sur tous les processeurs:
  envoyer_broadcast(post, 0);
  return post;
}

void Traitement_particulier_NS_plan::calculer_moyennes(DoubleTab& Moyennes, DoubleTab& rrhouf, DoubleTab& Racinederho) const
{
  const DoubleTab&            						temperature	= ref_champ_temperature_.valeur().valeurs();
  const Zone_dis_base&        						zdisbase	= mon_equation->inconnue().zone_dis_base();
  const Zone_VDF&            	 					zone_VDF	= ref_cast(Zone_VDF, zdisbase);
  const IntTab&              	  					elem_faces 	= zone_VDF.elem_faces();
  const Fluide_Incompressible&  					le_fluide 	= ref_cast(Fluide_Incompressible,mon_equation->milieu());
  const DoubleTab& 									vitesse 	= mon_equation->inconnue().valeurs();
  const DoubleTab& 									visco_dyn 	= le_fluide.viscosite_dynamique();
  const DoubleTab& 									tab_rho_elem= le_fluide.masse_volumique();
  const Fluide_Quasi_Compressible& 					fluide_QC	= ref_cast(Fluide_Quasi_Compressible,le_fluide);
  const DoubleTab& 									lambda 		= le_fluide.conductivite();
  const DoubleTab& 							   tab_rho_face = fluide_QC.rho_face_n();// rho discretisee aux faces.

  int REFY =0 ;
  int dimension  = Objet_U::dimension;
  int nb_elems   = zone_VDF.zone().nb_elem();
  int num_elem;
  int face_x_0, face_x_1, face_y_0, face_y_1, face_z_0, face_z_1;
  int Nyt = ijk_splitting_.get_grid_geometry().get_nb_elem_tot(DIRECTION_K);;
  double total_elem = 0. ; // variable qui recupere le nombre d'element dans un plan divise par Nyt.
  // on declare les variables qui vont stocker les derivees.

  Moyennes.resize(Nyt, nb_moyennes_temporelles(), Array_base::NOCOPY_NOINIT);
  Moyennes = 0.;

  double U1 , U2 , U3 ; // vitesses au centre des mailles
  double V1 , V2 , V3 ; // vitesses au centre des mailles * rho^(1/2)
  double b;// b = rho^(-1/2) C = b^(-1) * T
  double C;


  for (num_elem=0; num_elem<nb_elems; num_elem++)
    {
      REFY = Ref_Y(num_elem) ;
      face_x_0 = elem_faces(num_elem,0);
      face_x_1 = elem_faces(num_elem,dimension);
      face_y_0 = elem_faces(num_elem,1);
      face_y_1 = elem_faces(num_elem,1+dimension);
      face_z_0 = elem_faces(num_elem,2);
      face_z_1 = elem_faces(num_elem,2+dimension);

      rrhouf[face_x_0] = sqrt(tab_rho_face[face_x_0]) * vitesse[face_x_0] ;
      rrhouf[face_x_1] = sqrt(tab_rho_face[face_x_1]) * vitesse[face_x_1] ;
      rrhouf[face_y_0] = sqrt(tab_rho_face[face_y_0]) * vitesse[face_y_0] ;
      rrhouf[face_y_1] = sqrt(tab_rho_face[face_y_1]) * vitesse[face_y_1] ;
      rrhouf[face_z_0] = sqrt(tab_rho_face[face_z_0]) * vitesse[face_z_0] ;
      rrhouf[face_z_1] = sqrt(tab_rho_face[face_z_1]) * vitesse[face_z_1] ;

      // FA 1/09/11 : pour eviter de moyenner localement la vitesse u et w
      //
      // on est obliger de ramener la vitesse au centre des elements ou pas.
      U1 = .5*(vitesse[face_x_0]+vitesse[face_x_1]);
      U2 = .5*(vitesse[face_y_0]+vitesse[face_y_1]);
      U3 = .5*(vitesse[face_z_0]+vitesse[face_z_1]);

      Racinederho[num_elem] = sqrt(tab_rho_elem[num_elem]);

      // on calcule au centre de l'element meme si cela de doit pas etre utile.
      V1 = 0.5 * ( rrhouf[face_x_0] + rrhouf[face_x_1] ) ;
      V2 = 0.5 * ( rrhouf[face_y_0] + rrhouf[face_y_1] ) ;
      V3 = 0.5 * ( rrhouf[face_z_0] + rrhouf[face_z_1] ) ;

      C = racinederho[num_elem]*temperature(num_elem);
      b = 1/ racinederho[num_elem];

      // debut du codage de recuperation des moyennes pour le calcul des fluctuations.

      Moyennes(REFY,0) 	+= U1;
      Moyennes(REFY,1) 	+= U2;
      Moyennes(REFY,2) 	+= U3;
      Moyennes(REFY,3) 	+= V1;
      Moyennes(REFY,4) 	+= V2;
      Moyennes(REFY,5) 	+= V3;
      Moyennes(REFY,6)	+= b;
      Moyennes(REFY,7) 	+= temperature[num_elem];
      Moyennes(REFY,8) 	+= C;
      Moyennes(REFY,9)	+= visco_dyn[num_elem];
      Moyennes(REFY,10) 	+= vitesse[face_x_0];
      Moyennes(REFY,11) 	+= vitesse[face_x_1];
      Moyennes(REFY,12) 	+= vitesse[face_y_0];
      Moyennes(REFY,13) 	+= vitesse[face_y_1];
      Moyennes(REFY,14) 	+= vitesse[face_z_0];
      Moyennes(REFY,15) 	+= vitesse[face_z_1];

      Moyennes(REFY,37) 	+= rrhouf[face_x_0] ;
      Moyennes(REFY,38) 	+= rrhouf[face_x_1] ;
      Moyennes(REFY,39) 	+= rrhouf[face_y_0] ;
      Moyennes(REFY,40) 	+= rrhouf[face_y_1] ;
      Moyennes(REFY,41) 	+= rrhouf[face_z_0] ;
      Moyennes(REFY,42) 	+= rrhouf[face_z_1] ;

      Moyennes(REFY,34)	+= lambda[num_elem] ;
      Moyennes(REFY,35)	+= tab_rho_elem[num_elem] ;
      Moyennes(REFY,36)	+= U1 * temperature[num_elem] ;

      // Ici on calcule les RMS.
      Moyennes(REFY,16)	+= U1*U2; // on remplit 16 pour calculer 16 comme <U.V>-<U>.<V>
      Moyennes(REFY,17) 	+= U1*U1;
      Moyennes(REFY,18) 	+= U2*U2;
      Moyennes(REFY,19) 	+= U3*U3;
      Moyennes(REFY,20) 	+= V1*V1;
      Moyennes(REFY,21) 	+= V2*V2;
      Moyennes(REFY,22) 	+= V3*V3;
      Moyennes(REFY,23)	+= b*b;
      Moyennes(REFY,24) 	+= temperature[num_elem]*temperature[num_elem];//	- Moyennes_temporelle(REFY,7);
      Moyennes(REFY,25) 	+= C*C;
      Moyennes(REFY,26)	+= visco_dyn[num_elem]*visco_dyn[num_elem];//- Moyennes_temporelle(REFY,9);
      Moyennes(REFY,27) 	+= vitesse[face_x_0]*vitesse[face_x_0];
      Moyennes(REFY,28) 	+= vitesse[face_x_1]*vitesse[face_x_1];
      Moyennes(REFY,29) 	+= vitesse[face_y_0]*vitesse[face_y_0];
      Moyennes(REFY,30) 	+= vitesse[face_y_1]*vitesse[face_y_1];//	- Moyennes_temporelle(REFY,13);
      Moyennes(REFY,31) 	+= vitesse[face_z_0]*vitesse[face_z_0];
      Moyennes(REFY,32) 	+= vitesse[face_z_1]*vitesse[face_z_1];//	- Moyennes_temporelle(REFY,15);
    }

  //on calcule les moyennes globales de toutes les variables.
  // a calculer en preparation avec un argument e passer sera plus efficace.
  total_elem = (int)(mp_sum(nb_elems)/double(Nyt)); // nombre total d'element

  // pour chaque case du tableau, somme sur tous les processeurs.
  mp_sum_for_each_item(ref_cast(ArrOfDouble, Moyennes));

  Moyennes /= total_elem ;
}

void Traitement_particulier_NS_plan::post_traitement_particulier()
{
  const double current_time = mon_equation->inconnue().temps();
  int test_start_time = (current_time > start_time_);
  envoyer_broadcast(test_start_time, 0);
  if (! test_start_time)
    return;

  const int post = calculer_critere_post();
  int nb_pas_dt = mon_equation->schema_temps().nb_pas_dt();
  if (je_suis_maitre() && nb_pas_dt%30 == 0)
    {
      test_convergence(Moyennes_temporelles_, sauv_moyennes_,
                       current_time, temps_de_moyenne_,
                       "suivi_convergence.dat");
      stocker_moy(Moyennes_temporelles_, sauv_moyennes_, temps_de_moyenne_);
    }

  DoubleTab moyennes_spatiales;
  calculer_moyennes(moyennes_spatiales, Rrhouf, racinederho);

  double dt = mon_equation->schema_temps().pas_de_temps();

  if (post && je_suis_maitre())
    {
      Moyennes_temporelles_.ajoute(dt, moyennes_spatiales); // Moyenne_temporelle += dt * moyennes_spatiales
      temps_de_moyenne_ += dt;

      Nom filename;
      filename = Nom("Spatiale_") + Nom(current_time);
      ecriture_fichiers_moyennes(ijk_splitting_, moyennes_spatiales, filename, current_time, 0.);

      const ArrOfDouble& znodes = ijk_splitting_.get_grid_geometry().get_node_coordinates(DIRECTION_K);
      ArrOfDouble coord_k(znodes.size_array()-1);
      for (int i = 0; i < znodes.size_array()-1; i++)
        coord_k[i] = (znodes[i] + znodes[i+1]) * 0.5;
      ecrire_reynolds_tau(moyennes_spatiales, current_time, coord_k, "Parietal.dat");

      filename = Nom("Moyennes_") + Nom(current_time);
      ecriture_fichiers_moyennes(ijk_splitting_, Moyennes_temporelles_, filename, current_time, temps_de_moyenne_);

      filename = Nom("Stat_post_traitement_") + Nom(current_time);
      sauver_stat_plan(filename, coord_k, Moyennes_temporelles_, sauv_moyennes_, current_time, temps_de_moyenne_);
    }
  const Fluide_Quasi_Compressible& fluide_QC = ref_cast(Fluide_Quasi_Compressible, mon_equation->milieu());
  double Pth = fluide_QC.pression_th();
  double dpth  = Pth - Pth_passe_; // coefficent  R /(pth * cp ) * dpth/dt
  Pth_passe_ = dpth;
  dpth /=  dt; // ceci est calculer maintenant pour avoir l'evolution a chaque pas de temps.
  if (post && oui_calcul_)
    {
      DoubleTab val_post(Ref_Y.size(), nb_valeurs_spatiales_x()); // trop grand ???

      // on donne bel et bien moyenne spatiale qui a deja ete divisee par Delta t
      calculer_valeur_spatiale_vitesse_rho_mu(val_post, Moyennes_finale_, dpth); // ca j'ai verifier ca marche
      ecriture_val_post(val_post, current_time);
    }
}

void Traitement_particulier_NS_plan::ecriture_val_post(const DoubleTab& val_post, double current_time) const
{
  // Ecriture des valeurs sur les differents plans:
  ArrOfInt indices_plans_locaux;
  // Determine quels plans sont sur ce processeur:
  indices_plans_locaux.set_smart_resize(1);
  // accolade sans struture ???
  {
    for (int i = 0; i < indices_planes_.size_array(); i++)
      {
        const int local_k_index = indices_planes_[i] - ijk_splitting_.get_offset_local(DIRECTION_K);
        if (local_k_index >= 0 && local_k_index < ijk_splitting_.get_nb_elem_local(DIRECTION_K))
          {
            indices_plans_locaux.append_array(local_k_index);
          }
      }
  }
  // Ouverture des fichiers une seule fois pour toutes les variables (sur lustre, l'ouverture de fichier
  // est longue)
  // Les vecteurs de fichiers n'existent pas... oblige de faire une allocation dynamique
  SFichier *vecteur_fichiers = new SFichier[indices_plans_locaux.size_array()];
  SFichier *vecteur_fichiers_bin = new SFichier[indices_plans_locaux.size_array()];

  IJK_Field_float ijk_field; // on ouvre en float car on a pas besoin de plus de precision
  ijk_field.allocate(ijk_splitting_, IJK_Splitting::ELEM, 0 /* ghost size */);

  DoubleVect vdf_tmp;
  const Zone_dis_base& zdisbase = mon_equation->inconnue().zone_dis_base();
  const Zone_VDF& zone_VDF = ref_cast(Zone_VDF, zdisbase);
  zone_VDF.zone().creer_tableau_elements(vdf_tmp);

  const int ni = ijk_field.ni();
  const int nj = ijk_field.nj();
  const double delta_x = ijk_field.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double delta_y = ijk_field.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);


  long long binary_file_offset_in_bytes = 0;
  int data_size = ni * nj * sizeof(float);

  for (int ival = 0; ival < val_post.dimension(1); ival++)
    {
      // Copier les donnees pour une valeur post-traitee dans un champ et transformer en ijk:
      const int nb_elem_vdf = val_post.dimension(0);
      for (int i = 0; i < nb_elem_vdf; i++)
        vdf_tmp[i] = val_post(i, ival);

      // Methode parallele a appeler simultanement sur tous les processeurs:
      vdf_to_ijk_.convert_to_ijk(vdf_tmp, ijk_field);

      // Ecrire les donnees pour les plans a post-traiter situes sur ce processeur
      // ok mais moi j'ai seulement ce qui se trouve sur le 0.
      // quand et comment les autres ecrivent ?
      // je pige pas !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //
      for (int iplan = 0; iplan < indices_plans_locaux.size_array(); iplan++)
        {
          int local_k_layer = indices_plans_locaux[iplan];
          int global_k_layer = local_k_layer + ijk_splitting_.get_offset_local(DIRECTION_K);
          const ArrOfDouble& znodes = ijk_splitting_.get_grid_geometry().get_node_coordinates(DIRECTION_K);  //  recupere liste des noeuds
          double coord_z = (znodes[global_k_layer] + znodes[global_k_layer+1]) * 0.5; // calcul coords en z du plan a post traite
          Nom filename = Objet_U::nom_du_cas() + Nom("_") + Nom(current_time) + Nom("_y_") + Nom(coord_z); // nomme les fichiers
          Nom filename_bin = filename + Nom(".bin");

          SFichier& fic = vecteur_fichiers[iplan];
          SFichier& fic_bin = vecteur_fichiers_bin[iplan];
          if (ival == 0)
            {
              //journal() << "j'ouvre le fichier pour le proc " << me();
              fic.ouvrir(filename);
              fic_bin.set_bin(1);
              fic_bin.ouvrir(filename_bin);
              fic << "# Temps y Nval Nx Nz dx dz"    << finl;
              fic << current_time << " " <<  coord_z << " " << val_post.dimension(1) << " " ;
              fic << ni << "  " << nj << " " << delta_x << " " << delta_y <<  finl;
            }
          // On ecrit: la position en k en plan dans le maillage, le nom du fichier binaire, la position en octets dans le fichier, la taille en octets des valeurs
          fic << global_k_layer << " " << ival << " "
              << filename_bin << " " << (int)binary_file_offset_in_bytes << " " << data_size << finl;
          // les valeurs sont contigues dans la direction i, il peut y avoir des trous entre une rangee j et j+1
          // donc on n'ecrit pas tout d'un bloc:
          for (int i = 0; i < ni; i++)
            for (int j = 0; j < nj; j++)
              fic_bin.put(&ijk_field(i,j,local_k_layer), 1);
        }
      binary_file_offset_in_bytes += data_size;
    }

  delete [] vecteur_fichiers; // le destructeur de SFichier ferme les fichiers automatiquement
  delete [] vecteur_fichiers_bin;
}

#if 0
void Traitement_particulier_NS_plan::preparer_calcul_particulier()
{
  const Fluide_Incompressible& 		le_fluide 	= ref_cast(Fluide_Incompressible,mon_equation->milieu());
  const Fluide_Quasi_Compressible& 	fluide_QC   = ref_cast(Fluide_Quasi_Compressible,le_fluide);
  const Zone_dis_base& 			zdisbase	= mon_equation->inconnue().zone_dis_base();
  const Zone_VF& 			zone_VF		= ref_cast(Zone_VF, zdisbase);
  const Zone_VDF&            	 	zone_VDF	= ref_cast(Zone_VDF, zdisbase);
  const DoubleTab& 			xp 		= zone_VF.xp();
  const DoubleTab& 			vitesse 	= mon_equation->inconnue().valeurs();
  int 					nb_elems 	= zone_VF.zone().nb_elem();

  int num_elem	= 0; // compteur des elements
  int i,j,k,l;	  // compteurs des boucles
  int Nmoyenne = 43 ; // nombre de moyennes calculees
  Pth_passe = fluide_QC.pression_th();
  remplir_XYZ(X,Y,Z,Nx,Ny,Nz,Tab_recap); // renvoie vers Traitement_particulier_NS_plan_VDF ou Traitement_particulier_NS_plan_VEF

  remplir_reordonne_Z_tot(Z,Z_tot);
  remplir_reordonne_Y_tot(Y,Y_tot);
  remplir_reordonne_X_tot(X,X_tot);
  if ( oui_calcul )
    {
      if (!oui_repr)
        {
          Cerr << "L'option converger de Trait-part-plan e besoin d'une moyenne a partir de laquelle demarrer" << finl;
          Cerr << "Utiliser 'reprise nom_du_fichier_de_reprise' dans le jeu de donnee" << finl;
          Process::exit();
        }
      else
        envoyer_broadcast(oui_calcul,0);
    }
  int nb_elem_tot = zone_VDF.zone().nb_elem_tot(); // nombre total d'elements (reel + fict)

  racinederho.resize(nb_elem_tot) ;
  Rrhouf = vitesse ;

  remplir_Ref_Y(Ref_Y) ;
  Moyennes_temporelle.resize(Y_tot.size(),Nmoyenne); // on dimentionne le tableau des moyennes par plan
  Moyennes_temporelle = 0 ;
  if (je_suis_maitre())
    sauv_moyennes.resize(Y_tot.size(),2,3); // on va sauvez 3 moyennes en temps de <W^2> et <W>. il n'y a que le maitre qui en a besoin pour le moment.

  if (oui_repr) // si on fait une reprise des stats
    {
      reprendre_stat_plan(fich_repr,Moyennes_temporelle); // recupere les statistiques
      envoyer_broadcast(Moyennes_temporelle,0);
      envoyer_broadcast(sauv_moyennes,0);
      envoyer_broadcast(temps_de_moyenne,0);
      if (oui_calcul)
        {
          Moyennes_finale_.resize(Y_tot.size(),Nmoyenne);
          Moyennes_finale_ = Moyennes_temporelle;
          Moyennes_finale_ /= temps_de_moyenne ; // on calcule la norme,
        }
    }
  else // sinon tout e zero et on enregistre le temps.
    {
      Moyennes_temporelles_ = 0 ;
      sauv_moyennes_ = 0 ;
      temps_de_moyenne_ = 0 ;
    }

  Tab_post.resize(Np,X.size(),Z.size()); // fixe la taille du tableau
  Envoyer_val.resize(Np); // tableau qui va chercher a savoir si il faut envoyer les donnees

  Yp.resize(Np);
  int maxX = mp_max( X.size() ) ;
  int maxZ = mp_max( Z.size() ) ;
  DoubleTrav ptoTX(maxX);
  DoubleTrav ptoTZ(maxZ);
  ptoTX = -1 ;
  ptoTZ = -1 ;
  for (i=0; i<Np; i++)
    {
      if (yt[i]>Y_tot.size())
        {
          Cerr << "Erreur dans le choix des numeros de plan " <<  yt[i] << " est trop grand par rapport au domaine " << finl;
          Process::exit();
        }
      Yp[i]=Y_tot[yt[i]];
      Envoyer_val[i] = 0 ;// on en profite pour l'initialise
    }

  Yp.ordonne_array();

  // on initialise a -1 le tableau ceci sera utile apres.
  Tab_post = -1 ;
  for ( num_elem=0; num_elem<nb_elems; num_elem++) // boucle sur les elements
    for (i = 0 ; i <Np; i++) // boucles sur les plans a post traite
      if(xp(num_elem,1)==Yp[i]) // si c'est un plan a post traite
        {
          Envoyer_val[i] = 1 ; // marqueur globale qui permet d'envoyer l'info uniquement si elle est utile.

          for(j=0; j<X.size(); j++) // boucle sur les X
            if(xp(num_elem,0)==X(j))
              {
                break; // on cherche l'element correspondant dans X_tot
              }

          for(k=0; k<Z.size(); k++) //Boucle sur les Z
            if(xp(num_elem,2)==Z(k))
              {
                break; // Idem pour Z_tot
              }

          Tab_post(i,j,k) = num_elem ;

          for(l=0; l<X_tot.size(); l++) // boucle sur les X
            if(xp(num_elem,0)==X_tot(l))
              {
                break; // on cherche l'element correspondant dans X_tot
              }
          ptoTX(j)= l;
          for(l=0; l<Z_tot.size(); l++) //Boucle sur les Z
            if(xp(num_elem,2)==Z_tot(l))
              {
                break; // Idem pour Z_tot
              }
          ptoTZ(k) = l;

          break; // casse la boucle sur les Np .
        }
  envoyer(ptoTX ,Process::me(),0,0); // on dis si on va post traite.
  envoyer(ptoTZ ,Process::me(),0,0); // on dis si on va post traite.

  if(je_suis_maitre())
    {
      PtoTX.resize(maxX,Process::nproc());
      PtoTZ.resize(maxZ,Process::nproc());
      for(int p=0; p<Process::nproc(); p++)
        {
          recevoir(ptoTX,p,0,0);
          recevoir(ptoTZ,p,0,0);
          for (j = 0 ; j <ptoTX.size(); j++)
            PtoTX(j,p) = ptoTX(j) ;
          for (j = 0 ; j <ptoTZ.size(); j++)
            PtoTZ(j,p) = ptoTZ(j) ;
        }
    }

  if (!je_suis_maitre())
    {
      int dummy = 0;
      recevoir(dummy, 0, 0); // attendre que le processeur 0 reclame le message pour ne pas engorger le reseau
      envoyer(Envoyer_val, 0, 0); // envoyer tout le tableau
    }
  else
    {
      for_each_plane_recv_flag_.resize(Np, Process::nproc());
      for (int i = 0; i < Process::nproc(); i++)
        {
          IntTab tmp; // il faut recevoir le meme type que Envoyer_val
          if (i == 0)
            {
              tmp = Envoyer_val;
            }
          else
            {
              int dummy = 1;
              envoyer(dummy, i, 0); // signaler au processeur i que j'attends son message
              recevoir(tmp, i, 0);
            }
          for (int plan = 0; plan < Np; plan++)
            for_each_plane_recv_flag_(plan, i) = tmp[plan];
        }
    }

}


void Traitement_particulier_NS_plan::post_traitement_particulier()
{
  const Fluide_Incompressible&     le_fluide = ref_cast(Fluide_Incompressible,mon_equation->milieu());
  const Fluide_Quasi_Compressible& fluide_QC = ref_cast(Fluide_Quasi_Compressible,le_fluide);
  double 							dt		= mon_equation->schema_temps().pas_de_temps();
  double 			   				Pth 	= fluide_QC.pression_th();
  int i,j,k,p, val; // compteur de boucles.
  Nom fichier1;

  double dpth  = Pth - Pth_passe; // coefficent  R /(pth * cp ) * dpth/dt
  dpth /=  dt; // ceci est calculer maintenant pour avoir l'evolution a chaque pas de temps.

  double tps = mon_equation->inconnue().temps();
  Nom tampon ;
  Nom temps = Nom(tps);

  DoubleTrav moyennes_spatiales(Moyennes_temporelle);
  DoubleTrav Val_plan(X.size(),Z.size());
  DoubleTrav Val_plan_tot ;// on redimensionnera le tableau plus tard.
  if (tps>temps_deb)
    {
      double stationnaire_atteint = mon_equation->schema_temps().stationnaire_atteint();
      double tps_passe		  = mon_equation->schema_temps().temps_courant();
      int nb_pas_dt_max	   	  = mon_equation->schema_temps().nb_pas_dt_max();
      int nb_pas_dt			    = mon_equation->schema_temps().nb_pas_dt();
      double temps_max		    = mon_equation->schema_temps().temps_max();
      int impr_inst;

      if ( dt_impr_moy<=(tps-tps_passe) )
        {
          impr_inst=1;
        }
      else
        {
          double i, j, epsilon = 1.e-8; // Voir Schema_Temps_base::limpr pour information sur epsilon et modf
          modf(tps/dt_impr_moy + epsilon, &i);
          modf(tps_passe/dt_impr_moy + epsilon, &j);
          impr_inst=(i>j);
        }
      if ( !(nb_pas_dt%30) )
        {
          if (je_suis_maitre()) // il faut regler une conditions sur les pas de temps.
            {
              test_convergence(Moyennes_temporelle) ; // on test la convergence.
              stocker_moy(Moyennes_temporelle) ;
            }
        }
      int post = ( (nb_pas_dt+1<=1) || impr_inst || (temps_max <= tps) || (nb_pas_dt_max <= nb_pas_dt+1) || stationnaire_atteint  );
      // on calcule dans un nouveau tableau la moyenne spatiale
      calculer_moyennes(moyennes_spatiales);

      if (post)
        {
          if(je_suis_maitre())
            {
              fichier1 = "Spatiale_" ;
              fichier1 += Nom(tps) ;
              ecriture_fichiers_moyennes(moyennes_spatiales,fichier1,0); // on anticipe la division pour qu'elle soit par un
              calcul_reynolds_tau(moyennes_spatiales); //calcule les grandeurs parietalles.
            }
        }
      moyennes_spatiales *= dt ; // on la pondere par dt
      Moyennes_temporelle += moyennes_spatiales ; // on ajoute a Moyennes temporelles
      temps_de_moyenne += dt;

      if ( post )
        {
          const Zone_dis_base& zdisbase=mon_equation->inconnue().zone_dis_base();
          const Zone_VDF& zone_VDF=ref_cast(Zone_VDF, zdisbase);

          int nb_elem_p = zone_VDF.zone().nb_elem();
          int Nxt  =X_tot.size() ;
          int Nzt  =Z_tot.size() ;

          if (oui_calcul)
            {
              DoubleTrav val_post(nb_elem_p,Nval);
              calculer_valeur_spatiale_vitesse_rho_mu(val_post,Moyennes_finale,dpth); // on donne bel et bien moyenne spatiale qui a deja ete divise par Delta t

              if (je_suis_maitre())
                {
                  Val_plan_tot.resize(Nxt,Nzt); //le proc 0 cree le gros tableau
                }

              for (i = 0 ; i < Np ; i++)/// on fait le traitement pour chaque plan.
                {
                  Schema_Comm schema; // cree un schema de communication
                  ArrOfInt send_pe_list; // liste des processeurs a qui "moi", je vais envoyer un message.
                  ArrOfInt recv_pe_list; // liste des processeurs de qui je vais recevoir

                  // on remplit les tableaux de communication.
                  if (je_suis_maitre())
                    {
                      recv_pe_list.set_smart_resize(1);
                      for (p = 1; p < Process::nproc(); p++)
                        if (for_each_plane_recv_flag_(i,p))
                          recv_pe_list.append_array(p);
                    }
                  else
                    {
                      if (Envoyer_val[i])
                        {
                          send_pe_list.resize_array(1);
                          send_pe_list[0] = 0; // j'envoie un seul message, au processeur 0
                        }
                    }

                  schema.set_send_recv_pe_list(send_pe_list, recv_pe_list); // on donne les infos au schema.

                  for (val=0; val<Nval; val++)
                    {
                      Val_plan = 0 ;
                      Val_plan_tot = 0 ; // remise a zero des tableaux au changement de variable.

                      if (Envoyer_val[i]) // si il existe un element a post traite sur ce proc.
                        for(j=0; j<X.size(); j++)
                          for(k=0; k<Z.size(); k++)
                            if (Tab_post(i,j,k)+ 1 ) // faux si Tab_post(i,j,k) = -1 (soit le point n'est pas a remplir)
                              Val_plan(j,k) = val_post( Tab_post(i,j,k), val ) ;

                      schema.begin_comm(); // on demarre les communications.
                      if (!je_suis_maitre() && Envoyer_val[i] ) // on charge les donnees a envoye Val_plan
                        schema.send_buffer(0) << Val_plan;

                      schema.echange_taille_et_messages(); // on envoi effectivement les donnees.

                      if (je_suis_maitre())
                        {
                          for (p = 0; p < Process::nproc(); p++)
                            {
                              if (for_each_plane_recv_flag_(i,p))  // si des donnees existent
                                {
                                  if (p)
                                    {
                                      schema.recv_buffer(p) >> Val_plan;  // pour p == 0 les donnees sont deja dans val_plan.
                                    }

                                  for(j=0; j<Val_plan.dimension(0); j++) // on copie le tableau a sa place
                                    for(k=0; k<Val_plan.dimension(1); k++)
                                      if ( ( PtoTX(j,p) + 1 ) && ( PtoTZ(k,p)+ 1 ) )
                                        Val_plan_tot( PtoTX(j,p) , PtoTZ(k,p) ) = Val_plan(j,k);
                                }
                            }
                        }
                      schema.end_comm();

                      if(je_suis_maitre())// Ecriture des termes spectraux.
                        {
                          fichier1  = nom_du_cas() + "_" ;
                          fichier1 += temps;
                          fichier1 += "_y_";
                          fichier1 += Nom(Yp[i]);
                          ecriture_fichiers_spatial(Val_plan_tot,fichier1,i,val);
                        }

                    } // Boucle sur les grandeurs
                } // Boucle sur les plan.
            }
          if(je_suis_maitre())
            {
              fichier1  = "Moyennes_";
              fichier1 += temps ;
              ecriture_fichiers_moyennes(Moyennes_temporelle,fichier1,temps_de_moyenne);
              sauver_stat_plan(Moyennes_temporelle);
            }

        }// if ( (nb_pas_dt+1<=1) || impr_inst || (temps_max <= tps) || (nb_pas_dt_max <= nb_pas_dt+1) || stationnaire_atteint)
    } //  if (tps>temps_deb)

  Pth_passe  = Pth;
}

#endif
