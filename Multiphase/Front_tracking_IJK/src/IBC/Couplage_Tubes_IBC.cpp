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
// File      : Couplage_Tubes_IBC.cpp
// Directory : $IJK_ROOT/src/IBC
//
/////////////////////////////////////////////////////////////////////////////
#include <Couplage_Tubes_IBC.h>
#include <Param.h>
#include <Interprete_bloc.h>
#include <Linear_algebra_tools.h>
#include <communications.h>
#include <SFichier.h>
#include <EFichier.h>
#include <IJK_Navier_Stokes_tools.h>

Implemente_instanciable(Couplage_Tubes_IBC, "Couplage_Tubes_IBC", Objet_U);
// Fonction qui permet d'afficher toutes les composantes d'un Vecteur3
static Sortie& operator<<(Sortie& os, Vecteur3 v)
{
  os << "[ " << v[0] << " " << v[1] << " " << v[2] << " ]";
  return os;
}
// Pour ecrire l'objet dans un fichier ou a l'ecran (pour l'instant ne fait rien)
Sortie& Couplage_Tubes_IBC::printOn(Sortie& os) const
{
  return os;
}
#define select(a,x,y,z) ((a==0)?(x):((a==1)?(y):(z)))

// Lecture de l'objet dans le jeu de donnees
Entree& Couplage_Tubes_IBC::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("vitesse_pour_adim", &vitesse_pour_adim_, Param::REQUIRED);
  param.ajouter("rho_fluide_pour_adim", &rho_fluide_pour_adim_, Param::REQUIRED);
  param.ajouter("lissage", &lissage_, Param::REQUIRED);
  param.dictionnaire("etalement", ETALEMENT); // Attention le dictionnaire s'applique sur le dernier param ajoute !!!!
  param.dictionnaire("pas_etalement", PAS_ETALEMENT);
  param.ajouter("epaisseur_lissage", &epaisseur_lissage_, Param::REQUIRED);
  Nom nomobjet_faisceau_tubes;
  param.ajouter("faisceau_tubes", &nomobjet_faisceau_tubes);
  param.ajouter("champ_miroir", &champ_miroir_, Param::REQUIRED);
  param.dictionnaire("miroir", MIROIR);
  param.dictionnaire("symetrie_plane", SYMETRIE_PLANE);
  param.dictionnaire("pas_miroir", PAS_MIROIR);
  param.ajouter("methode_IBC", &methode_IBC_, Param::REQUIRED);
  param.dictionnaire("ibc0", IBC0);
  param.dictionnaire("ibc_diffuse", IBC_DIFFUSE);
  param.dictionnaire("ibc_localisee", IBC_LOCALISEE);
  param.dictionnaire("ibc_localisee_qdm", IBC_LOCALISEE_QDM);
  param.ajouter("solide", &solide_, Param::REQUIRED);
  param.dictionnaire("cylindre", CYLINDRE);
  param.dictionnaire("cube", CUBE);
  param.ajouter("n_P", &n_P_, Param::REQUIRED);
  param.ajouter("L_cube", &L_cube_, Param::REQUIRED);
  param.lire_avec_accolades(is);


  if (nomobjet_faisceau_tubes != "??")
    {
      faisceau_ = ref_cast(Faisceau_Tubes, Interprete_bloc::objet_global(nomobjet_faisceau_tubes));
      Cout << "On a lu les donnees dans le fichier " << nomobjet_faisceau_tubes << ":" << finl;
      Cout << faisceau_; // affiche les valeurs contenues dans le tableau donnes_tubes_ pr verifier qu'on a bien copier les valeurs donnees en entre
    }
  return is;
}

void Couplage_Tubes_IBC::initialize(const IJK_Splitting& splitting)
{
  ref_splitting_ = splitting; // copie de splitting pour l'utiliser dans la class Couplage_Tubes_IBC

  const int n_tubes_total = faisceau_.size();
  integrale_force_.resize(n_tubes_total , 3);
  integrale_force_post_proj_.resize(n_tubes_total , 3);
  integrale_force_pression_.resize(n_tubes_total , 2);
  pression_teta_.resize(n_tubes_total , n_P_);
  d_integrale_force_.resize(n_tubes_total , 3);
  d_integrale_force_post_proj_.resize(n_tubes_total , 3);
  integrale_force_N_moins_1_.resize(n_tubes_total , 3);
  integrale_force_N_moins_1_post_proj_.resize(n_tubes_total , 3);
  volume_cylindres_.resize(n_tubes_total , 3);
  masse_fluide_cylindres_.resize(n_tubes_total, 3);

  for (int i = 0; i < n_tubes_total; i++)
    faisceau_[i].valeur().initialize(splitting); // la fonction initialize() appartient a la class Tube_base, on l'appelle ici pour pouvoir accader a splitting dans cette class
}

void Couplage_Tubes_IBC::sauvegarder_probleme(SFichier& f)
{
  if (Process::je_suis_maitre())
    {
      f << "{\n"
        << "  faisceau_tubes " << faisceau_ << "\n"
        << "  integrale_force_ " << integrale_force_ << "\n"
        << "  integrale_force_N_moins_1_ " << integrale_force_N_moins_1_ << "\n"
        << "}\n";
    }
}

void Couplage_Tubes_IBC::sauvegarder_pression(SFichier& f)
{
  const int n_tubes_total = faisceau_.size();
  const double teta = 360./ n_P_;
  const double Adim_P = 1. / ( 0.5 * rho_fluide_pour_adim_ * vitesse_pour_adim_ * vitesse_pour_adim_);
  DoubleTab P_teta_Adim = pression_teta_;
  P_teta_Adim *= Adim_P;

  if (Process::je_suis_maitre())
    {
      f.precision(17);
      f.setf(std::ios_base::scientific);
      //f << "P_Pa teta_Rad" << endl;
      for (int itube = 0; itube < n_tubes_total; itube++)
        for (int l = 0; l < n_P_; l++)
          f  << P_teta_Adim(itube, l) << " " << 180 - teta * l  <<"\n" ;
    }
}

void Couplage_Tubes_IBC::reprendre_probleme(Entree& fichier)
{
  Param param(que_suis_je());
  param.ajouter("faisceau_tubes", &faisceau_);
  param.ajouter("integrale_force_", &integrale_force_);
  param.ajouter("integrale_force_N_moins_1_", &integrale_force_N_moins_1_);
  param.lire_avec_accolades(fichier);
}

void Couplage_Tubes_IBC::force_ibc(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                   const IJK_Field_double& rho_fluide_field,
                                   const double timestep, const IJK_Field_double& pressure, double current_time)
{
  // ATTENTION update() est appele dans le run() avant force_ibc_velocity()
  // donc quand on appelle force_ibc_velocity(), le timestep a change, c'est celui d'apres utilise dans update() au pas de temps d'avant
  // Ex : on demarre au pas de temps n-1, on rentre ds le run(), on appelle d'abord update() donc on utilise timestep=t(n)-t(n-1)
  // donc quand on entre dans force_ibc() on va calculer F_IBC(N) mais le pas de temps est toujours timestep=t(n)-t(n-1)
  // ce qui permet de calculer d_integrale_force_ = ( F_IBC(N)-F_IBC(N-1) ) / timestep
  // ensuite on passe au pas de temps suivant donc au pas de temps n et donc timestep=t(n+1)-t(n)
  // on entre d'abord dans update() on peut donc calculer F_IBC_extrapole car on stocke d_integrale_force_ precedemment!

  integrale_force_N_moins_1_ = integrale_force_;

  switch(champ_miroir_)
    {
    case PAS_MIROIR:

      switch(solide_)
        {
        case CYLINDRE:

          switch(methode_IBC_)
            {
            case IBC_DIFFUSE:

              force_ibc_velocity_frac_vol(vx, vy, vz,
                                          rho_fluide_field,
                                          masse_fluide_cylindres_,
                                          volume_cylindres_,
                                          integrale_force_,
                                          faisceau_);

              break;
            case IBC0:

              force_ibc_velocity(vx, vy, vz,
                                 rho_fluide_field,
                                 masse_fluide_cylindres_,
                                 volume_cylindres_,
                                 integrale_force_,
                                 faisceau_);


              break;
            default:
              Cerr << "Erreur dans Couplage_Tubes_IBC::force_ibc: methode_IBC " << methode_IBC_ << " non code" << finl;
              exit();
            }

          break;
        case CUBE:

          switch(methode_IBC_)
            {

            case IBC0:
              ibc0_velocity_cube(vx, vy, vz,
                                 rho_fluide_field,
                                 masse_fluide_cylindres_,
                                 volume_cylindres_,
                                 integrale_force_,
                                 faisceau_);


              break;

            case IBC_DIFFUSE:

              ibc_diffuse_velocity_cube(vx, vy, vz,
                                        rho_fluide_field,
                                        masse_fluide_cylindres_,
                                        volume_cylindres_,
                                        integrale_force_,
                                        faisceau_);


              break;

            case IBC_LOCALISEE:
              ibc_localisee_velocity_cube(vx, vy, vz,
                                          rho_fluide_field,
                                          masse_fluide_cylindres_,
                                          volume_cylindres_,
                                          integrale_force_,
                                          faisceau_, current_time);

              break;

            case IBC_LOCALISEE_QDM:
              ibc_localisee_velocity_cube_qdm(vx, vy, vz,
                                              rho_fluide_field,
                                              masse_fluide_cylindres_,
                                              volume_cylindres_,
                                              integrale_force_,
                                              faisceau_, current_time);

              break;


            default:
              Cerr << "Erreur dans Couplage_Tubes_IBC::force_ibc: methode_IBC " << methode_IBC_ << " non code" << finl;
              exit();
            }
          break;
        default:
          Cerr << "Erreur dans Couplage_Tubes_IBC::force_ibc: solide " << solide_ << " non code" << finl;
          exit();
        }

      break;
    case MIROIR:
      force_ibc_velocity_miroir(vx, vy, vz,
                                rho_fluide_field,
                                masse_fluide_cylindres_,
                                volume_cylindres_,
                                integrale_force_,
                                faisceau_);
      break;
    case SYMETRIE_PLANE:
      force_ibc_velocity_symetrie_plane(vx, vy, vz,
                                        rho_fluide_field,
                                        masse_fluide_cylindres_,
                                        volume_cylindres_,
                                        integrale_force_,
                                        faisceau_);
      break;
    default:
      Cerr << "Erreur dans Couplage_Tubes_IBC::force_ibc: champ_miroir " << champ_miroir_ << " non code" << finl;
      exit();
    }


  // Calcul des forces de pression
  //calcul_F_pression(pressure, vx, integrale_force_pression_, pression_teta_, faisceau_);
  calcul_F_pression2(pressure, vx, integrale_force_pression_, pression_teta_, faisceau_);

  double inv_timestep = (1./timestep);
  integrale_force_ *= inv_timestep;
  // on calcule : ( F_IBC(N)-F_IBC(N-1) ) / timestep, ou timestep = t(n) - t(n-1)
  // d_integrale_force_ = (integrale_force_ - integrale_force_N_moins_1_) * inv_timestep; mais si on ecrit comme ca en c++ ca plante
  d_integrale_force_ = integrale_force_;
  d_integrale_force_ -= integrale_force_N_moins_1_;
  d_integrale_force_ *= inv_timestep;

}


/////////////////////////////////////////////////
/// Calcul force ibc post projection
//////////////////////////////////////////////////


void Couplage_Tubes_IBC::calcul_force_post_projection(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                                      const IJK_Field_double& rho_fluide_field,
                                                      const double timestep, double current_time)
{
  // ATTENTION update() est appele dans le run() avant force_ibc_velocity()
  // donc quand on appelle force_ibc_velocity(), le timestep a change, c'est celui d'apres utilise dans update() au pas de temps d'avant
  // Ex : on demarre au pas de temps n-1, on rentre ds le run(), on appelle d'abord update() donc on utilise timestep=t(n)-t(n-1)
  // donc quand on entre dans force_ibc() on va calculer F_IBC(N) mais le pas de temps est toujours timestep=t(n)-t(n-1)
  // ce qui permet de calculer d_integrale_force_ = ( F_IBC(N)-F_IBC(N-1) ) / timestep
  // ensuite on passe au pas de temps suivant donc au pas de temps n et donc timestep=t(n+1)-t(n)
  // on entre d'abord dans update() on peut donc calculer F_IBC_extrapole car on stocke d_integrale_force_ precedemment!

  integrale_force_N_moins_1_post_proj_ = integrale_force_post_proj_;

  switch(methode_IBC_)
    {

    case IBC0:
      ibc0_force_cube(vx, vy, vz,
                      rho_fluide_field,
                      masse_fluide_cylindres_,
                      volume_cylindres_,
                      integrale_force_post_proj_,
                      faisceau_);


      break;

    case IBC_DIFFUSE:

      ibc_diffuse_force_cube(vx, vy, vz,
                             rho_fluide_field,
                             masse_fluide_cylindres_,
                             volume_cylindres_,
                             integrale_force_post_proj_,
                             faisceau_);


      break;
# if 0
    case IBC_LOCALISEE:
      ibc_localisee_force_cube(vx, vy, vz,
                               rho_fluide_field,
                               masse_fluide_cylindres_,
                               volume_cylindres_,
                               integrale_force_post_proj_,
                               faisceau_, current_time);

      break;
    case IBC_LOCALISEE_QDM:
      ibc_localisee_velocity_cube_qdm(vx, vy, vz,
                                      rho_fluide_field,
                                      masse_fluide_cylindres_,
                                      volume_cylindres_,
                                      integrale_force_,
                                      faisceau_, current_time);

      break;


#endif
    default:
      Cerr << "Erreur dans Couplage_Tubes_IBC::force_ibc_post_proj: methode_IBC " << methode_IBC_ << " non code" << finl;
      exit();
    }

  double inv_timestep = (1./timestep);
  integrale_force_post_proj_ *= inv_timestep;
  // on calcule : ( F_IBC(N)-F_IBC(N-1) ) / timestep, ou timestep = t(n) - t(n-1)
  // d_integrale_force_ = (integrale_force_ - integrale_force_N_moins_1_) * inv_timestep; mais si on ecrit comme ca en c++ ca plante
  d_integrale_force_post_proj_ = integrale_force_post_proj_;
  d_integrale_force_post_proj_ -= integrale_force_N_moins_1_post_proj_;
  d_integrale_force_post_proj_ *= inv_timestep;

}



/////////////////////////////
//// calcul du champ miroir pour permettre de mieux calculer le RHS ensuite
/////////////////////////////



void Couplage_Tubes_IBC::champ_miroir(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                      const IJK_Field_double& rho_fluide_field,
                                      const double timestep, const IJK_Field_double& pressure)
{

  DoubleTab tmp_integrale_force = integrale_force_;
  force_ibc_velocity_miroir(vx, vy, vz,
                            rho_fluide_field,
                            masse_fluide_cylindres_,
                            volume_cylindres_,
                            tmp_integrale_force,
                            faisceau_);
  // le champ miroir est stocke dans velocity_tmp_

}




/////////////////////////
//// forcage direct ordre 0
////////////////////////////////

void Couplage_Tubes_IBC::force_ibc_velocity(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                            const IJK_Field_double& rho_fluide_field,
                                            DoubleTab& masse_fluide_cylindres,
                                            DoubleTab& volume_cylindres,
                                            DoubleTab& integrale_force,
                                            const Faisceau_Tubes& faisceau) const
{
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  ArrOfDouble nodes_pos[3]; //defini un tableau dynamique de 3 lignes pour la position des nodes
  ArrOfDouble elem_pos[3];  //defini un tableau dynamique de 3 lignes pour la position des elements
  {
    for (int dir = 0; dir < 3; dir++)
      {
        nodes_pos[dir] = geom.get_node_coordinates(dir); // node=noeud, les noeuds sont a l'intersection des mailles; on ne connait que la position des noeuds
        const int n = nodes_pos[dir].size_array()-1; // n= nombre de noeuds-1 = nombre d'elements dans la direction dir; les alements sont au centre des mailles
        elem_pos[dir].resize_array(n); // donne la taille du tableau contenant les elements
        for (int i = 0; i < n; i++)
          elem_pos[dir][i] = (nodes_pos[dir][i] + nodes_pos[dir][i+1]) * 0.5; // la position de l'element selon dir est moyenne de la positions des noeuds de part et d'autre
      }
  }
  const int offset_i = splitting.get_offset_local(DIRECTION_I); // donne l'offset selon la direction i
  const int offset_k = splitting.get_offset_local(DIRECTION_K);  // donne l'offset selon la direction k

  const int ntubes = faisceau.size();
  const double epsilon = epaisseur_lissage_ * geom.get_constant_delta(DIRECTION_I); // definition de la demi largeur de la zone d'etalement du terme de forcage comme un nombre de fois la largeur de maille selon la direction x

  // Boucle sur les tubes du faisceau
  for (int itube = 0; itube < ntubes; itube++)
    {
      const Tube_base& tube = faisceau[itube].valeur();
      const double tube_r = tube.get_rayon();
      const double omega = tube.get_omega();
      const double r_tube_min = tube_r - epsilon;
      const double r_tube_max = tube_r + epsilon; // rayon inf et max de la zone d'etalement du terme de forcage
      const double r_tube_min_carre = r_tube_min * r_tube_min;
      const double r_tube_max_carre = r_tube_max * r_tube_max;

      const Vecteur3 position_cylindre = tube.get_current_pos();
      const Vecteur3 vitesse_cylindre = tube.get_current_velocity();

      // Boucle sur les trois composantes de qdm
      for (int direction = 0; direction < 3; direction++)   // boucle sur les 3 directions de la vitesse
        {

          const int ivoisin = (direction == 0) ? -1 : 0;
          const int jvoisin = (direction == 1) ? -1 : 0;
          const int kvoisin = (direction == 2) ? -1 : 0;

          const double volume_maille = geom.get_constant_delta(0) * geom.get_constant_delta(1) * geom.get_constant_delta(2); //donne le volume de la maille
          IJK_Field_double& velocity = select(direction, vx, vy, vz);  // velocity est egale a la vitesse selon la direction 0,1 ou 2
          int imin = 0;
          int jmin = 0;
          int kmin = 0;
          int imax = velocity.ni(); // donne le nombre de valeur de la vitesse selon i
          int jmax = velocity.nj();
          int kmax = velocity.nk();
          {
            // restreint le domaine ijk balaye a la zone couverte par le cylindre:
            const ArrOfDouble& x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I] : elem_pos[DIRECTION_I];
            const double xmin = position_cylindre[0] - (tube_r + epsilon); // on n'oublie pas que l'on doit chercher dans la zone courverte par le cylindre + etalement
            const double xmax = position_cylindre[0] + (tube_r + epsilon);
            while (imin < velocity.ni() && x_coord[imin + offset_i] < xmin)
              imin++;
            while (imax > 0 && x_coord[imax + offset_i - 1] > xmax)
              imax--;
          }
          {
            const ArrOfDouble& z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K] : elem_pos[DIRECTION_K];
            const double zmin = position_cylindre[2] - (tube_r + epsilon);
            const double zmax = position_cylindre[2] + (tube_r + epsilon);
            while (kmin < velocity.nk() && z_coord[kmin + offset_k] < zmin)
              kmin++;
            while (kmax > 0 && z_coord[kmax + offset_k - 1] > zmax)
              kmax--;
          }

          double delta_qdm_cylindre = 0.;
          double masse_fluide = 0.; // Masse de fluide dans le cylindre
          double volume = 0.;

          for (int k = kmin; k < kmax; k++)
            {
              // si direction==DIRECTION_I alors x_ccord = coordonnee du noeud dans la direction I donc x_noeud et z_coord = z_element
              // si direction==DIRECTION_K alors x_coord = x_element et z_coord = z_noeud
              // x_coord et z_coord sont les abcisses des vitesse qui sont situee au centre des faces
              const double z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K][k+offset_k] : elem_pos[DIRECTION_K][k+offset_k];
              for (int j = jmin; j < jmax; j++)
                {
                  //const double y_coord = (direction==DIRECTION_J) ? nodes_pos[DIRECTION_J][j+offset_j] : elem_pos[DIRECTION_J][j+offset_j];
                  for (int i = imin; i < imax; i++)
                    {
                      const double x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I][i+offset_i] : elem_pos[DIRECTION_I][i+offset_i];
                      const double delta_x = x_coord - position_cylindre[0];
                      const double delta_z = z_coord - position_cylindre[2];
                      const double d2 = delta_x * delta_x + delta_z * delta_z;
                      const double d = sqrt(d2); // distance du point M au centre du tube
                      const double rayon_tube_carre = tube_r * tube_r;


                      double critere_max_carre, critere_min_carre;
                      switch(lissage_)
                        {
                        case PAS_ETALEMENT:
                          critere_max_carre = rayon_tube_carre;
                          critere_min_carre = rayon_tube_carre;
                          break;
                        case ETALEMENT:
                          critere_max_carre = r_tube_max_carre;
                          critere_min_carre = r_tube_min_carre;
                          break;
                        default:
                          critere_max_carre = rayon_tube_carre; // pour initialiser, pour le compilateur
                          critere_min_carre = rayon_tube_carre;
                          Cerr << "Erreur dans Couplage_Tubes_IBC::update: lissage " << lissage_ << " non code" << finl;
                          exit();
                        }

                      if (d2 < critere_max_carre)
                        {
                          Vecteur3 rotation;
                          rotation[0] = omega * (- delta_z);
                          rotation[1] = 0.;
                          rotation[2] = omega * delta_x;
                          double extrapolated_v = vitesse_cylindre[direction] + rotation[direction];
                          double f_etalement;
                          if (d2 < critere_min_carre)
                            {
                              f_etalement = 1.;
                            }
                          else   // si on a pas d2 < r_tube_min_carre alors d2 > r_tube_min_carre
                            {
                              //f_etalement = 0.5 - 0.5 * (d - tube_r) / epsilon - 0.5 / M_PI * sin( M_PI * (d - tube_r)/ epsilon); //  pour d= tube_r la derive de f_etalement est de 1
                              f_etalement = 0.5 * ( 1 - tanh(4* (d - tube_r)/ epsilon) ) ; // pour d= tube_r la derive de f_etalement est de 10
                            }

                          double delta_v = f_etalement * (extrapolated_v - velocity(i,j,k)); // qdm ds la zone d'etalement

                          // Calcul de rho sur la face: moyenne des rho sur les elements voisins
                          double rho_fluide = (rho_fluide_field(i,j,k) + rho_fluide_field(i+ivoisin, j+jvoisin, k+kvoisin)) * 0.5; // rho sur la face est la moy des deux cotes
                          delta_qdm_cylindre += delta_v * volume_maille * rho_fluide;
                          velocity(i, j, k) = f_etalement * extrapolated_v + (1-f_etalement) * velocity(i, j, k); // forcage etale dans la zone d'etalement de largeur 2 * epsilon
                          masse_fluide += volume_maille * rho_fluide;
                          volume += volume_maille;
                        }
                    }
                }
            }

          integrale_force(itube, direction) = delta_qdm_cylindre;
          masse_fluide_cylindres(itube, direction) = masse_fluide;
          volume_cylindres(itube, direction) = volume;
        }
    }
  mp_sum_for_each_item(integrale_force);
  mp_sum_for_each_item(masse_fluide_cylindres);
  mp_sum_for_each_item(volume_cylindres);
}





/////////////////////////
//// forcage direct cylindre frac vol simplifie
////////////////////////////////

void Couplage_Tubes_IBC::force_ibc_velocity_frac_vol(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                                     const IJK_Field_double& rho_fluide_field,
                                                     DoubleTab& masse_fluide_cylindres,
                                                     DoubleTab& volume_cylindres,
                                                     DoubleTab& integrale_force,
                                                     const Faisceau_Tubes& faisceau) const
{
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  ArrOfDouble nodes_pos[3]; //defini un tableau dynamique de 3 lignes pour la position des nodes
  ArrOfDouble elem_pos[3];  //defini un tableau dynamique de 3 lignes pour la position des elements
  {
    for (int dir = 0; dir < 3; dir++)
      {
        nodes_pos[dir] = geom.get_node_coordinates(dir); // node=noeud, les noeuds sont a l'intersection des mailles; on ne connait que la position des noeuds
        const int n = nodes_pos[dir].size_array()-1; // n= nombre de noeuds-1 = nombre d'elements dans la direction dir; les alements sont au centre des mailles
        elem_pos[dir].resize_array(n); // donne la taille du tableau contenant les elements
        for (int i = 0; i < n; i++)
          elem_pos[dir][i] = (nodes_pos[dir][i] + nodes_pos[dir][i+1]) * 0.5; // la position de l'element selon dir est moyenne de la positions des noeuds de part et d'autre
      }
  }
  const int offset_i = splitting.get_offset_local(DIRECTION_I); // donne l'offset selon la direction i
  const int offset_k = splitting.get_offset_local(DIRECTION_K);  // donne l'offset selon la direction k

  const int ntubes = faisceau.size();


  // Boucle sur les tubes du faisceau
  for (int itube = 0; itube < ntubes; itube++)
    {
      const Tube_base& tube = faisceau[itube].valeur();
      const double tube_r = tube.get_rayon();
      const double omega = tube.get_omega();


      const Vecteur3 position_cylindre = tube.get_current_pos();
      const Vecteur3 vitesse_cylindre = tube.get_current_velocity();

      // Boucle sur les trois composantes de qdm
      for (int direction = 0; direction < 3; direction++)   // boucle sur les 3 directions de la vitesse
        {

          const int ivoisin = (direction == 0) ? -1 : 0;
          const int jvoisin = (direction == 1) ? -1 : 0;
          const int kvoisin = (direction == 2) ? -1 : 0;

          const double volume_maille = geom.get_constant_delta(0) * geom.get_constant_delta(1) * geom.get_constant_delta(2); //donne le volume de la maille
          IJK_Field_double& velocity = select(direction, vx, vy, vz);  // velocity est egale a la vitesse selon la direction 0,1 ou 2
          int imin = 0;
          int jmin = 0;
          int kmin = 0;
          int imax = velocity.ni(); // donne le nombre de valeur de la vitesse selon i
          int jmax = velocity.nj();
          int kmax = velocity.nk();
          //      const double d_x = geom.get_constant_delta(0);
          const double d_z = geom.get_constant_delta(2);
          {
            // restreint le domaine ijk balaye a la zone couverte par le cylindre:
            const ArrOfDouble& x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I] : elem_pos[DIRECTION_I];
            const double xmin = position_cylindre[0] - (tube_r + d_z); // on n'oublie pas que l'on doit chercher dans la zone courverte par le cylindre + etalement
            const double xmax = position_cylindre[0] + (tube_r + d_z);
            //Cout << "xmin : " << xmin << " xmax : " <<  xmax<<finl;
            while (imin < velocity.ni() && x_coord[imin + offset_i] < xmin)
              imin++;
            while (imax > 0 && x_coord[imax + offset_i - 1] > xmax)
              imax--;
          }
          {
            const ArrOfDouble& z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K] : elem_pos[DIRECTION_K];
            const double zmin = position_cylindre[2] - (tube_r + d_z);
            const double zmax = position_cylindre[2] + (tube_r + d_z);
            //Cout << "zmin : " << zmin << " zmax : " <<  zmax<<finl;
            while (kmin < velocity.nk() && z_coord[kmin + offset_k] < zmin)
              kmin++;
            while (kmax > 0 && z_coord[kmax + offset_k - 1] > zmax)
              kmax--;
          }

          double delta_qdm_cylindre = 0.;
          double masse_fluide = 0.; // Masse de fluide dans le cylindre
          double volume = 0.;
          //Cout << "imin : " << imin << " imax : " << imax << " kmin : " << kmin << " kmax : " << kmax << " x_O : " << position_cylindre[0] << " z_O : " << position_cylindre[2] << finl;

          for (int k = kmin; k <= kmax; k++)
            {
              // si direction==DIRECTION_I alors x_ccord = coordonnee du noeud dans la direction I donc x_noeud et z_coord = z_element
              // si direction==DIRECTION_K alors x_coord = x_element et z_coord = z_noeud
              // x_coord et z_coord sont les abcisses des vitesse qui sont situee au centre des faces
              const double z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K][k+offset_k] : elem_pos[DIRECTION_K][k+offset_k];
              for (int j = jmin; j < jmax; j++)
                {
                  //const double y_coord = (direction==DIRECTION_J) ? nodes_pos[DIRECTION_J][j+offset_j] : elem_pos[DIRECTION_J][j+offset_j];
                  for (int i = imin; i <= imax; i++)
                    {
                      const double x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I][i+offset_i] : elem_pos[DIRECTION_I][i+offset_i];
                      const double delta_x = x_coord - position_cylindre[0];
                      const double delta_z = z_coord - position_cylindre[2];
                      const double d2 = delta_x * delta_x + delta_z * delta_z;
                      const double d = sqrt(d2); // distance du point M au centre du tube


                      // attention on suppose que d_x = d_z !!
                      const double Borne_min = tube_r - d_z /2;
                      const double Borne_max = tube_r + d_z /2;

                      //Cout << "i : " << i << " j : " << j << " k : " << k << " dir : " << direction << " x_coord : " << x_coord << " z_coord : " ;
                      //Cout << z_coord <<  " d_OM : " << d << " velocity(i,j,k) : " << velocity(i,j,k) << " Borne_max : " << Borne_max << " Borne_min : " << Borne_min <<  " d_z/2 : " << d_z/2 << finl;


                      if (d <= Borne_max)
                        {
                          Vecteur3 rotation;
                          rotation[0] = omega * (- delta_z);
                          rotation[1] = 0.;
                          rotation[2] = omega * delta_x;
                          double extrapolated_v = vitesse_cylindre[direction] + rotation[direction];
                          double f;
                          const double d_occupe = Borne_max - d;
                          if (d <= Borne_min)
                            f = 1.;
                          else
                            f = d_occupe / d_z ;
                          double delta_v_non_etale = extrapolated_v - velocity(i,j,k);
                          double delta_v = f * delta_v_non_etale; // qdm ds la zone d'etalement
                          //Cout << " delta_v : " << delta_v << " extrapolated_v : " << extrapolated_v << " velocity(i,j,k) : " << velocity(i,j,k) << " f : " << f << " d_occupe : " << d_occupe << " delta_v_non_etale : " << delta_v_non_etale << finl;

                          // Calcul de rho sur la face: moyenne des rho sur les elements voisins
                          double rho_fluide = (rho_fluide_field(i,j,k) + rho_fluide_field(i+ivoisin, j+jvoisin, k+kvoisin)) * 0.5; // rho sur la face est la moy des deux cotes
                          delta_qdm_cylindre += delta_v * volume_maille * rho_fluide;
                          velocity(i, j, k) = f * extrapolated_v + (1-f) * velocity(i, j, k); // forcage etale dans la zone d'etalement de largeur 2 * epsilon
                          //Cout << " delta_qdm_cylindre : " << delta_qdm_cylindre << " volume_maille : " << volume_maille << " rho_fluide : " << rho_fluide << " velocity(i, j, k) : " << velocity(i, j, k) << finl;
                          masse_fluide += volume_maille * rho_fluide;
                          volume += volume_maille;
                        }
                    }
                }
            }

          integrale_force(itube, direction) = delta_qdm_cylindre;
          masse_fluide_cylindres(itube, direction) = masse_fluide;
          volume_cylindres(itube, direction) = volume;
          //Cout << " integrale_force(itube, direction) " << integrale_force(itube, direction) << " volume_cylindres(itube, direction) : " << volume_cylindres(itube, direction) << finl;
        }
    }
  mp_sum_for_each_item(integrale_force);
  mp_sum_for_each_item(masse_fluide_cylindres);
  mp_sum_for_each_item(volume_cylindres);
}




/////////////////////////////
//// forcer les mailles qui vont passer solide avt la conv + diff afin de ne pas avoir de saut en vitesse trop violent quand on fera le forcage derriere
/////////////////////////////



void Couplage_Tubes_IBC::forcage_anticipe(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                          const IJK_Field_double& rho_fluide_field,
                                          const double timestep, const IJK_Field_double& pressure)
{

  DoubleTab tmp_integrale_force = integrale_force_;
  DoubleTab tmp_masse_fluide_cylindres = masse_fluide_cylindres_;
  DoubleTab tmp_volume_cylindres = volume_cylindres_;
  force_ibc_velocity_anticipe_cube(vx, vy, vz,
                                   rho_fluide_field,
                                   tmp_masse_fluide_cylindres,
                                   tmp_volume_cylindres,
                                   tmp_integrale_force,
                                   faisceau_,
                                   timestep);
  // le champ miroir est stocke dans velocity_tmp_

}

//////////////
/// forcage anticpe dans le cas du cube
////////////////////
void Couplage_Tubes_IBC::force_ibc_velocity_anticipe_cube(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                                          const IJK_Field_double& rho_fluide_field,
                                                          DoubleTab& masse_fluide_cylindres,
                                                          DoubleTab& volume_cylindres,
                                                          DoubleTab& integrale_force,
                                                          const Faisceau_Tubes& faisceau,
                                                          const double timestep) const
{
  Cout << "Fonction_anticipe" << finl;
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  ArrOfDouble nodes_pos[3]; //defini un tableau dynamique de 3 lignes pour la position des nodes
  ArrOfDouble elem_pos[3];  //defini un tableau dynamique de 3 lignes pour la position des elements
  {
    for (int dir = 0; dir < 3; dir++)
      {
        nodes_pos[dir] = geom.get_node_coordinates(dir); // node=noeud, les noeuds sont a l'intersection des mailles; on ne connait que la position des noeuds
        const int n = nodes_pos[dir].size_array()-1; // n= nombre de noeuds-1 = nombre d'elements dans la direction dir; les alements sont au centre des mailles
        elem_pos[dir].resize_array(n); // donne la taille du tableau contenant les elements
        for (int i = 0; i < n; i++)
          elem_pos[dir][i] = (nodes_pos[dir][i] + nodes_pos[dir][i+1]) * 0.5; // la position de l'element selon dir est moyenne de la positions des noeuds de part et d'autre
      }
  }
  const int offset_i = splitting.get_offset_local(DIRECTION_I); // donne l'offset selon la direction i
  const int offset_k = splitting.get_offset_local(DIRECTION_K);  // donne l'offset selon la direction k

  const int ntubes = faisceau.size();

  // Boucle sur les tubes du faisceau
  for (int itube = 0; itube < ntubes; itube++)
    {
      const Tube_base& tube = faisceau[itube].valeur();
      const double omega = tube.get_omega();


      const Vecteur3 position_cylindre = tube.get_current_pos();
      const Vecteur3 vitesse_cylindre = tube.get_current_velocity();

      // Boucle sur les trois composantes de qdm
      for (int direction = 0; direction < 3; direction++)   // boucle sur les 3 directions de la vitesse
        {

          const int ivoisin = (direction == 0) ? -1 : 0;
          const int jvoisin = (direction == 1) ? -1 : 0;
          const int kvoisin = (direction == 2) ? -1 : 0;

          const double volume_maille = geom.get_constant_delta(0) * geom.get_constant_delta(1) * geom.get_constant_delta(2); //donne le volume de la maille
          IJK_Field_double& velocity = select(direction, vx, vy, vz);  // velocity est egale a la vitesse selon la direction 0,1 ou 2
          int imin = 0;
          int jmin = 0;
          int kmin = 0;
          int imax = velocity.ni(); // donne le nombre de valeur de la vitesse selon i
          int jmax = velocity.nj();
          int kmax = velocity.nk();
          const double delta_x = geom.get_constant_delta(0);
          const double delta_z = geom.get_constant_delta(2);

          {
            // restreint le domaine ijk balaye a la zone couverte par le cylindre:
            const ArrOfDouble& x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I] : elem_pos[DIRECTION_I];
            const double xmin = position_cylindre[0] - (L_cube_ + delta_x); // on n'oublie pas que l'on doit chercher dans la zone courverte par le cylindre + etalement
            const double xmax = position_cylindre[0] + (L_cube_ + delta_x);
            Cout << "xmin : " << xmin << " xmax : " <<  xmax << finl;
            while (imin < velocity.ni() && x_coord[imin + offset_i] < xmin)
              imin++;
            while (imax > 0 && x_coord[imax + offset_i - 1] > xmax)
              imax--;
          }
          {
            const ArrOfDouble& z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K] : elem_pos[DIRECTION_K];
            const double zmin = position_cylindre[2] - (L_cube_ + delta_x);
            const double zmax = position_cylindre[2] + (L_cube_ + delta_x);
            Cout  << "zmin : " << zmin << " zmax : " << zmax <<finl;
            while (kmin < velocity.nk() && z_coord[kmin + offset_k] < zmin)
              kmin++;
            while (kmax > 0 && z_coord[kmax + offset_k - 1] > zmax)
              {
                //Cout  << "kmax : " << kmax << " z_coord[kmax + offset_k - 1] : " << z_coord[kmax + offset_k - 1] <<finl;
                kmax--;
                //Cout  << "kmax : " << kmax << " z_coord[kmax + offset_k - 1] : " << z_coord[kmax + offset_k - 1] <<finl;
              }
          }
          double delta_qdm_cylindre = 0.;
          double masse_fluide = 0.; // Masse de fluide dans le cylindre
          double volume = 0.;
          Cout << "imin : " << imin << " imax : " << imax << " kmin : " << kmin << " kmax : " << kmax << " x_O : " << position_cylindre[0] << " z_O : " << position_cylindre[2] << finl;
          for (int k = kmin; k <= kmax; k++)
            {
              // si direction==DIRECTION_I alors x_ccord = coordonnee du noeud dans la direction I donc x_noeud et z_coord = z_element
              // si direction==DIRECTION_K alors x_coord = x_element et z_coord = z_noeud
              // x_coord et z_coord sont les abcisses des vitesse qui sont situee au centre des faces
              const double z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K][k+offset_k] : elem_pos[DIRECTION_K][k+offset_k];
              for (int j = jmin; j < jmax; j++)
                {
                  //const double y_coord = (direction==DIRECTION_J) ? nodes_pos[DIRECTION_J][j+offset_j] : elem_pos[DIRECTION_J][j+offset_j];
                  for (int i = imin; i <= imax; i++)
                    {
                      const double x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I][i+offset_i] : elem_pos[DIRECTION_I][i+offset_i];
                      const double d_OM_x = x_coord - position_cylindre[0];
                      const double d_OM_z = z_coord - position_cylindre[2];
                      const double ABS_dOM_x = fabs(d_OM_x);
                      const double ABS_dOM_z = fabs(d_OM_z);
                      Cout << "i : " << i << " j : " << j << " k : " << k << " dir : " << direction << " x_coord : " << x_coord << " z_coord : ";
                      Cout << z_coord << " ABS_d_OM_x : " << ABS_dOM_x << " ABS_d_OM_z : " << ABS_dOM_z << " velocity(i,j,k) : " << velocity(i,j,k) <<  finl;
                      const double Borne_max = L_cube_;
                      Vecteur3 rotation;
                      rotation[0] = omega * (- delta_z);
                      rotation[1] = 0.;
                      rotation[2] = omega * delta_x;
                      const double vs_z = vitesse_cylindre[2] + rotation[2];
                      if ( vs_z >0)
                        {

                          if ((d_OM_z > Borne_max) && (ABS_dOM_x <= Borne_max+1e-10))
                            {
                              Cout << "vs_z > 0" << " i : " << i << " j : " << j << " k : " << k << " dir : " << direction << " x_coord : " << x_coord << " z_coord : " ;
                              Cout << z_coord <<" ABS_d_OM_x : " << ABS_dOM_x << " d_OM_z : " << d_OM_z <<  " velocity(i,j,k) : " << velocity(i,j,k)<<  finl;

                              const double d_OM_L = d_OM_z - L_cube_;
                              const double dz   = vs_z * timestep;
                              const double ABS_dz = fabs(dz);
                              Cout << "d_OM_L : " << d_OM_L << " dz : " << dz << " vs_z : " << vs_z << finl;
                              if  ( d_OM_L <= ABS_dz  )
                                {
                                  double extrapolated_v = vitesse_cylindre[direction] + rotation[direction];
                                  //double extrapolated_v = velocity(i, j, k-1);
                                  double delta_v = (extrapolated_v - velocity(i,j,k));
                                  Cout << "On force la vitesse en prevision " << " delta_v : " << delta_v << " extrapolated_v : " << extrapolated_v << " velocity(i,j,k) : " << velocity(i,j,k) << finl;

                                  // Calcul de rho sur la face: moyenne des rho sur les elements voisins
                                  double rho_fluide = (rho_fluide_field(i,j,k) + rho_fluide_field(i+ivoisin, j+jvoisin, k+kvoisin)) * 0.5; // rho sur la face est la moy des deux cotes
                                  delta_qdm_cylindre += delta_v * volume_maille * rho_fluide;
                                  Cout << " delta_qdm_cylindre : " << delta_qdm_cylindre << " volume_maille : " << volume_maille << finl;
                                  velocity(i, j, k) = extrapolated_v; // forcage etale dans la zone d'etalement de largeur 2 * epsilon
                                  masse_fluide += volume_maille * rho_fluide;
                                  volume += volume_maille;

                                }
                            }
                        }

                      if ( vs_z <0)
                        {

                          if ((-d_OM_z > Borne_max) && (ABS_dOM_x <= Borne_max+1e-10))
                            {
                              Cout << "vs_z < 0" << " i : " << i << " j : " << j << " k : " << k << " dir : " << direction << " x_coord : " << x_coord << " z_coord : ";
                              Cout<< z_coord <<" ABS_d_OM_x : " << ABS_dOM_x << " d_OM_z : " << d_OM_z << " velocity(i,j,k) : " << velocity(i,j,k) <<  finl;

                              const double d_OM_L = d_OM_z - L_cube_;
                              const double dz   = vs_z * timestep;
                              const double ABS_dz = fabs(dz);
                              Cout << "d_OM_L : " << d_OM_L << " dz : " << dz << " vs_z : " << vs_z << finl;
                              if  ( d_OM_L <= ABS_dz  )
                                {
                                  double extrapolated_v = vitesse_cylindre[direction] + rotation[direction];
                                  //double extrapolated_v = velocity(i, j, k-1);
                                  double delta_v = (extrapolated_v - velocity(i,j,k));
                                  Cout << "On force la vitesse en prevision " << " delta_v : " << delta_v << " extrapolated_v : " << extrapolated_v << " velocity(i,j,k) : " << velocity(i,j,k) << finl;

                                  // Calcul de rho sur la face: moyenne des rho sur les elements voisins
                                  double rho_fluide = (rho_fluide_field(i,j,k) + rho_fluide_field(i+ivoisin, j+jvoisin, k+kvoisin)) * 0.5; // rho sur la face est la moy des deux cotes
                                  delta_qdm_cylindre += delta_v * volume_maille * rho_fluide;
                                  Cout << " delta_qdm_cylindre : " << delta_qdm_cylindre << " volume_maille : " << volume_maille << finl;
                                  velocity(i, j, k) = extrapolated_v; // forcage etale dans la zone d'etalement de largeur 2 * epsilon
                                  masse_fluide += volume_maille * rho_fluide;
                                  volume += volume_maille;

                                }
                            }

                        }


                    }
                }
            }

          integrale_force(itube, direction) = delta_qdm_cylindre;
          masse_fluide_cylindres(itube, direction) = masse_fluide;
          volume_cylindres(itube, direction) = volume;
          Cout << " integrale_force(itube, direction) " << integrale_force(itube, direction) << " volume_cylindres(itube, direction) : " << volume_cylindres(itube, direction) << finl;
        }
    }
  mp_sum_for_each_item(integrale_force);
  mp_sum_for_each_item(masse_fluide_cylindres);
  mp_sum_for_each_item(volume_cylindres);
}



///////////////////////////////////////////////////////
//////////////////////////////////////////////////////
//         CUBE
//         CUBE
//         CUBE
/////////////////////////////////////////////////////
////////////////////////////////////////////////////



/////////////////////////
//// IBC0 cube
////////////////////////////////

void Couplage_Tubes_IBC::ibc0_velocity_cube(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                            const IJK_Field_double& rho_fluide_field,
                                            DoubleTab& masse_fluide_cylindres,
                                            DoubleTab& volume_cylindres,
                                            DoubleTab& integrale_force,
                                            const Faisceau_Tubes& faisceau) const
{
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  ArrOfDouble nodes_pos[3]; //defini un tableau dynamique de 3 lignes pour la position des nodes
  ArrOfDouble elem_pos[3];  //defini un tableau dynamique de 3 lignes pour la position des elements
  {
    for (int dir = 0; dir < 3; dir++)
      {
        nodes_pos[dir] = geom.get_node_coordinates(dir); // node=noeud, les noeuds sont a l'intersection des mailles; on ne connait que la position des noeuds
        const int n = nodes_pos[dir].size_array()-1; // n= nombre de noeuds-1 = nombre d'elements dans la direction dir; les alements sont au centre des mailles
        elem_pos[dir].resize_array(n); // donne la taille du tableau contenant les elements
        for (int i = 0; i < n; i++)
          elem_pos[dir][i] = (nodes_pos[dir][i] + nodes_pos[dir][i+1]) * 0.5; // la position de l'element selon dir est moyenne de la positions des noeuds de part et d'autre
      }
  }
  const int offset_i = splitting.get_offset_local(DIRECTION_I); // donne l'offset selon la direction i
  const int offset_k = splitting.get_offset_local(DIRECTION_K);  // donne l'offset selon la direction k

  const int ntubes = faisceau.size();

  // Boucle sur les tubes du faisceau
  for (int itube = 0; itube < ntubes; itube++)
    {
      const Tube_base& tube = faisceau[itube].valeur();
      const double omega = tube.get_omega();


      const Vecteur3 position_cylindre = tube.get_current_pos();
      const Vecteur3 vitesse_cylindre = tube.get_current_velocity();


      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      // Affichage de la vitesse predite
      Cout << "v_predite : " << vx(64, 0, 73) << " " << vx(64, 0, 74) << " " <<vx(64, 0, 75) << " " <<vx(64, 0, 76) << " z_interface : " <<  position_cylindre[2] + L_cube_ << finl;
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      // Boucle sur les trois composantes de qdm
      for (int direction = 0; direction < 3; direction++)   // boucle sur les 3 directions de la vitesse
        {

          const int ivoisin = (direction == 0) ? -1 : 0;
          const int jvoisin = (direction == 1) ? -1 : 0;
          const int kvoisin = (direction == 2) ? -1 : 0;

          const double volume_maille = geom.get_constant_delta(0) * geom.get_constant_delta(1) * geom.get_constant_delta(2); //donne le volume de la maille
          IJK_Field_double& velocity = select(direction, vx, vy, vz);  // velocity est egale a la vitesse selon la direction 0,1 ou 2
          int imin = 0;
          int jmin = 0;
          int kmin = 0;
          int imax = velocity.ni(); // donne le nombre de valeur de la vitesse selon i
          int jmax = velocity.nj();
          int kmax = velocity.nk();
          const double delta_x = geom.get_constant_delta(0);
          const double delta_z = geom.get_constant_delta(2);
          {
            // restreint le domaine ijk balaye a la zone couverte par le cylindre:
            const ArrOfDouble& x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I] : elem_pos[DIRECTION_I];
            const double xmin = position_cylindre[0] - (L_cube_ + delta_x); // on n'oublie pas que l'on doit chercher dans la zone courverte par le cube
            const double xmax = position_cylindre[0] + (L_cube_ + delta_x);

            while (imin < velocity.ni() && x_coord[imin + offset_i] < xmin)
              imin++;
            while (imax > 0 && x_coord[imax + offset_i - 1] > xmax)
              imax--;
          }
          {
            const ArrOfDouble& z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K] : elem_pos[DIRECTION_K];
            const double zmin = position_cylindre[2] - (L_cube_ + delta_z);
            const double zmax = position_cylindre[2] + (L_cube_ + delta_z);

            while (kmin < velocity.nk() && z_coord[kmin + offset_k] < zmin)
              kmin++;
            while (kmax > 0 && z_coord[kmax + offset_k - 1] > zmax)
              {
                kmax--;
              }
          }
          double delta_qdm_cylindre = 0.;
          double masse_fluide = 0.; // Masse de fluide dans le cylindre
          double volume = 0.;

# if 0
          Cout << "                               " << finl;
          Cout << "dir : "<< direction << " x_O : " << position_cylindre[0] << " z_O : " << position_cylindre[2] << " V_maille : " << volume_maille << finl;
          Cout << "                               " << finl;
# endif
          for (int k = kmin; k <= kmax; k++)
            {
              // si direction==DIRECTION_I alors x_ccord = coordonnee du noeud dans la direction I donc x_noeud et z_coord = z_element
              // si direction==DIRECTION_K alors x_coord = x_element et z_coord = z_noeud
              // x_coord et z_coord sont les abcisses des vitesse qui sont situee au centre des faces
              const double z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K][k+offset_k] : elem_pos[DIRECTION_K][k+offset_k];
              for (int j = jmin; j < jmax; j++)
                {
                  //const double y_coord = (direction==DIRECTION_J) ? nodes_pos[DIRECTION_J][j+offset_j] : elem_pos[DIRECTION_J][j+offset_j];
                  for (int i = imin; i <= imax; i++)
                    {
                      const double x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I][i+offset_i] : elem_pos[DIRECTION_I][i+offset_i];
                      const double dOMx = x_coord - position_cylindre[0]; // distance du centre au point M selon x
                      const double dOMz = z_coord - position_cylindre[2]; // distance du centre au point M selon z
                      const double ABS_dOMx = fabs(dOMx);
                      const double ABS_dOMz = fabs(dOMz);
                      const double Borne_z = L_cube_;
                      const double Borne_x = L_cube_;


                      if (  (ABS_dOMx <= Borne_x)  && (ABS_dOMz <= Borne_z))
                        {
                          Vecteur3 rotation;
                          rotation[0] = omega * (- dOMz);
                          rotation[1] = 0.;
                          rotation[2] = omega * dOMx;
                          double extrapolated_v = vitesse_cylindre[direction] + rotation[direction];
                          // on force la vitesse a la vitesse du cylindre car : d_OM <= R
                          double delta_v = (extrapolated_v - velocity(i,j,k));

# if 0
                          Cout << "Forcage_interieur "<< "i : " << i << " j : " << j << " k : " << k <<  " x_M : " << x_coord << " z_M : " << z_coord << " v(i,j,k) : " << velocity(i,j,k) << finl;
# endif
                          // Calcul de rho sur la face: moyenne des rho sur les elements voisins
                          double rho_fluide = (rho_fluide_field(i,j,k) + rho_fluide_field(i+ivoisin, j+jvoisin, k+kvoisin)) * 0.5; // rho sur la face est la moy des deux cotes
                          delta_qdm_cylindre += delta_v * volume_maille * rho_fluide;

# if 0
                          Cout <<  "delta_v : "<< delta_v << " v_S : " << extrapolated_v << " qdm : " << delta_v * volume_maille * rho_fluide << " qdm_cube : " << delta_qdm_cylindre << finl;
# endif

                          //# if 0
                          // Sortie pour le cas ref de validation : cube de cote 22 mailles a Re10
                          // evolution de la vitesse et qdm dans deux mailles : une sortante et une entrante
                          if ( (direction==0) && ((j==0) && (i==64)) && ((k==54)|(k==55)|(k==76)|(k==77)) )
                            {
                              Cout << " " << finl;
                              Cout << "qdm_ijk" << i << j << k << " : " << delta_v * volume_maille * rho_fluide << " x_M : " << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2] + L_cube_ << finl;
                              Cout << "v_ijk" << i << j << k << " : " << velocity(i,j,k) << " x_M : " << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2] + L_cube_ << finl;
                              Cout << " " << finl;
                            }
                          //# endif

                          velocity(i, j, k) = extrapolated_v; // v_F = v_S car on est dans le cube
                          masse_fluide += volume_maille * rho_fluide;
                          volume += volume_maille;
                        }
                    } // fin de la boucle en i
                }
            }

          integrale_force(itube, direction) = delta_qdm_cylindre;
          masse_fluide_cylindres(itube, direction) = masse_fluide;
          volume_cylindres(itube, direction) = volume;
# if 0
          Cout << "integrale_force(itube, direction) " << integrale_force(itube, direction) << " V_S : " << volume_cylindres(itube, direction) << finl;
# endif
        }

      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      // Affichage de la vitesse imposee
      Cout << "v_imposee : " << vx(64, 0, 73) << " " << vx(64, 0, 74) << " " <<vx(64, 0, 75) << " " << vx(64, 0, 76) << " z_interface : " <<  position_cylindre[2] + L_cube_ << finl;
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    }
  mp_sum_for_each_item(integrale_force);
  mp_sum_for_each_item(masse_fluide_cylindres);
  mp_sum_for_each_item(volume_cylindres);
}

//////////////////////////////
//// IBC ibc0 FORCE cube
///////////////////////////////////



void Couplage_Tubes_IBC::ibc0_force_cube(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                         const IJK_Field_double& rho_fluide_field,
                                         DoubleTab& masse_fluide_cylindres,
                                         DoubleTab& volume_cylindres,
                                         DoubleTab& integrale_force,
                                         const Faisceau_Tubes& faisceau) const
{
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  ArrOfDouble nodes_pos[3]; //defini un tableau dynamique de 3 lignes pour la position des nodes
  ArrOfDouble elem_pos[3];  //defini un tableau dynamique de 3 lignes pour la position des elements
  {
    for (int dir = 0; dir < 3; dir++)
      {
        nodes_pos[dir] = geom.get_node_coordinates(dir); // node=noeud, les noeuds sont a l'intersection des mailles; on ne connait que la position des noeuds
        const int n = nodes_pos[dir].size_array()-1; // n= nombre de noeuds-1 = nombre d'elements dans la direction dir; les alements sont au centre des mailles
        elem_pos[dir].resize_array(n); // donne la taille du tableau contenant les elements
        for (int i = 0; i < n; i++)
          elem_pos[dir][i] = (nodes_pos[dir][i] + nodes_pos[dir][i+1]) * 0.5; // la position de l'element selon dir est moyenne de la positions des noeuds de part et d'autre
      }
  }
  const int offset_i = splitting.get_offset_local(DIRECTION_I); // donne l'offset selon la direction i
  const int offset_k = splitting.get_offset_local(DIRECTION_K);  // donne l'offset selon la direction k

  const int ntubes = faisceau.size();

  // Boucle sur les tubes du faisceau
  for (int itube = 0; itube < ntubes; itube++)
    {
      const Tube_base& tube = faisceau[itube].valeur();
      const double omega = tube.get_omega();


      const Vecteur3 position_cylindre = tube.get_current_pos();
      const Vecteur3 vitesse_cylindre = tube.get_current_velocity();

      // Boucle sur les trois composantes de qdm
      for (int direction = 0; direction < 3; direction++)   // boucle sur les 3 directions de la vitesse
        {

          const int ivoisin = (direction == 0) ? -1 : 0;
          const int jvoisin = (direction == 1) ? -1 : 0;
          const int kvoisin = (direction == 2) ? -1 : 0;

          const double volume_maille = geom.get_constant_delta(0) * geom.get_constant_delta(1) * geom.get_constant_delta(2); //donne le volume de la maille
          IJK_Field_double& velocity = select(direction, vx, vy, vz);  // velocity est egale a la vitesse selon la direction 0,1 ou 2
          int imin = 0;
          int jmin = 0;
          int kmin = 0;
          int imax = velocity.ni(); // donne le nombre de valeur de la vitesse selon i
          int jmax = velocity.nj();
          int kmax = velocity.nk();
          const double delta_x = geom.get_constant_delta(0);
          const double delta_z = geom.get_constant_delta(2);
          {
            // restreint le domaine ijk balaye a la zone couverte par le cylindre:
            const ArrOfDouble& x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I] : elem_pos[DIRECTION_I];
            const double xmin = position_cylindre[0] - (L_cube_ + delta_x); // on n'oublie pas que l'on doit chercher dans la zone courverte par le cube
            const double xmax = position_cylindre[0] + (L_cube_ + delta_x);

            while (imin < velocity.ni() && x_coord[imin + offset_i] < xmin)
              imin++;
            while (imax > 0 && x_coord[imax + offset_i - 1] > xmax)
              imax--;
          }
          {
            const ArrOfDouble& z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K] : elem_pos[DIRECTION_K];
            const double zmin = position_cylindre[2] - (L_cube_ + delta_z);
            const double zmax = position_cylindre[2] + (L_cube_ + delta_z);

            while (kmin < velocity.nk() && z_coord[kmin + offset_k] < zmin)
              kmin++;
            while (kmax > 0 && z_coord[kmax + offset_k - 1] > zmax)
              {
                kmax--;
              }
          }
          double delta_qdm_cylindre = 0.;
          double masse_fluide = 0.; // Masse de fluide dans le cylindre
          double volume = 0.;

          for (int k = kmin; k <= kmax; k++)
            {
              // si direction==DIRECTION_I alors x_ccord = coordonnee du noeud dans la direction I donc x_noeud et z_coord = z_element
              // si direction==DIRECTION_K alors x_coord = x_element et z_coord = z_noeud
              // x_coord et z_coord sont les abcisses des vitesse qui sont situee au centre des faces
              const double z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K][k+offset_k] : elem_pos[DIRECTION_K][k+offset_k];
              for (int j = jmin; j < jmax; j++)
                {
                  //const double y_coord = (direction==DIRECTION_J) ? nodes_pos[DIRECTION_J][j+offset_j] : elem_pos[DIRECTION_J][j+offset_j];
                  for (int i = imin; i <= imax; i++)
                    {
                      const double x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I][i+offset_i] : elem_pos[DIRECTION_I][i+offset_i];
                      const double dOMx = x_coord - position_cylindre[0]; // distance du centre au point M selon x
                      const double dOMz = z_coord - position_cylindre[2]; // distance du centre au point M selon z
                      const double ABS_dOMx = fabs(dOMx);
                      const double ABS_dOMz = fabs(dOMz);
                      const double Borne_z = L_cube_;
                      const double Borne_x = L_cube_;


                      if (  (ABS_dOMx <= Borne_x)  && (ABS_dOMz <= Borne_z))
                        {
                          Vecteur3 rotation;
                          rotation[0] = omega * (- dOMz);
                          rotation[1] = 0.;
                          rotation[2] = omega * dOMx;
                          double extrapolated_v = vitesse_cylindre[direction] + rotation[direction];
                          // on force la vitesse a la vitesse du cylindre car : d_OM <= R
                          double delta_v = (extrapolated_v - velocity(i,j,k));

                          // Calcul de rho sur la face: moyenne des rho sur les elements voisins
                          double rho_fluide = (rho_fluide_field(i,j,k) + rho_fluide_field(i+ivoisin, j+jvoisin, k+kvoisin)) * 0.5; // rho sur la face est la moy des deux cotes
                          delta_qdm_cylindre += delta_v * volume_maille * rho_fluide;

                          //# if 0
                          // Sortie pour le cas ref de validation : cube de cote 22 mailles a Re10
                          // evolution de la vitesse et qdm dans deux mailles : une sortante et une entrante
                          if ( (direction==0) && ((j==0) && (i==64)) && ((k==55)|(k==76)) )
                            {
                              Cout << " " << finl;
                              Cout << "qdm_ijk_post_proj" << i << j << k << " : " << delta_v * volume_maille * rho_fluide << " x_M : " << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2] + L_cube_ << finl;
                              Cout << " " << finl;
                            }
                          //# endif

                          masse_fluide += volume_maille * rho_fluide;
                          volume += volume_maille;
                        }
                    } // fin de la boucle en i
                }
            }

          integrale_force(itube, direction) = delta_qdm_cylindre;
          masse_fluide_cylindres(itube, direction) = masse_fluide;
          volume_cylindres(itube, direction) = volume;
        }
    }
  mp_sum_for_each_item(integrale_force);
  mp_sum_for_each_item(masse_fluide_cylindres);
  mp_sum_for_each_item(volume_cylindres);
}

/////////////////////////
//// IBC interface diffuse, une maille, Cube
////////////////////////////////

void Couplage_Tubes_IBC::ibc_diffuse_velocity_cube(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                                   const IJK_Field_double& rho_fluide_field,
                                                   DoubleTab& masse_fluide_cylindres,
                                                   DoubleTab& volume_cylindres,
                                                   DoubleTab& integrale_force,
                                                   const Faisceau_Tubes& faisceau) const
{
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  ArrOfDouble nodes_pos[3]; //defini un tableau dynamique de 3 lignes pour la position des nodes
  ArrOfDouble elem_pos[3];  //defini un tableau dynamique de 3 lignes pour la position des elements
  {
    for (int dir = 0; dir < 3; dir++)
      {
        nodes_pos[dir] = geom.get_node_coordinates(dir); // node=noeud, les noeuds sont a l'intersection des mailles; on ne connait que la position des noeuds
        const int n = nodes_pos[dir].size_array()-1; // n= nombre de noeuds-1 = nombre d'elements dans la direction dir; les alements sont au centre des mailles
        elem_pos[dir].resize_array(n); // donne la taille du tableau contenant les elements
        for (int i = 0; i < n; i++)
          elem_pos[dir][i] = (nodes_pos[dir][i] + nodes_pos[dir][i+1]) * 0.5; // la position de l'element selon dir est moyenne de la positions des noeuds de part et d'autre
      }
  }
  const int offset_i = splitting.get_offset_local(DIRECTION_I); // donne l'offset selon la direction i
  const int offset_k = splitting.get_offset_local(DIRECTION_K);  // donne l'offset selon la direction k

  const int ntubes = faisceau.size();

  // Boucle sur les tubes du faisceau
  for (int itube = 0; itube < ntubes; itube++)
    {
      const Tube_base& tube = faisceau[itube].valeur();
      const double omega = tube.get_omega();


      const Vecteur3 position_cylindre = tube.get_current_pos();
      const Vecteur3 vitesse_cylindre = tube.get_current_velocity();

      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      // Affichage de la vitesse predite
      Cout << "v_predite : " << vx(64, 0, 76) << " z_interface : " <<  position_cylindre[2] + L_cube_ << finl;
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      // Boucle sur les trois composantes de qdm
      for (int direction = 0; direction < 3; direction++)   // boucle sur les 3 directions de la vitesse
        {

          const int ivoisin = (direction == 0) ? -1 : 0;
          const int jvoisin = (direction == 1) ? -1 : 0;
          const int kvoisin = (direction == 2) ? -1 : 0;

          const double volume_maille = geom.get_constant_delta(0) * geom.get_constant_delta(1) * geom.get_constant_delta(2); //donne le volume de la maille
          IJK_Field_double& velocity = select(direction, vx, vy, vz);  // velocity est egale a la vitesse selon la direction 0,1 ou 2
          int imin = 0;
          int jmin = 0;
          int kmin = 0;
          int imax = velocity.ni(); // donne le nombre de valeur de la vitesse selon i
          int jmax = velocity.nj();
          int kmax = velocity.nk();
          const double delta_x = geom.get_constant_delta(0);
          const double delta_z = geom.get_constant_delta(2);
          {
            // restreint le domaine ijk balaye a la zone couverte par le cylindre:
            const ArrOfDouble& x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I] : elem_pos[DIRECTION_I];
            const double xmin = position_cylindre[0] - (L_cube_ + delta_x); // on n'oublie pas que l'on doit chercher dans la zone courverte par le cube + une maille
            const double xmax = position_cylindre[0] + (L_cube_ + delta_x);

            while (imin < velocity.ni() && x_coord[imin + offset_i] < xmin)
              imin++;
            while (imax > 0 && x_coord[imax + offset_i - 1] > xmax)
              imax--;
          }
          {
            const ArrOfDouble& z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K] : elem_pos[DIRECTION_K];
            const double zmin = position_cylindre[2] - (L_cube_ + delta_z);
            const double zmax = position_cylindre[2] + (L_cube_ + delta_z);

            while (kmin < velocity.nk() && z_coord[kmin + offset_k] < zmin)
              kmin++;
            while (kmax > 0 && z_coord[kmax + offset_k - 1] > zmax)
              {
                kmax--;
              }
          }
          double delta_qdm_cylindre = 0.;
          double masse_fluide = 0.; // Masse de fluide dans le cylindre
          double volume = 0.;
# if 0
          Cout << "                               " << finl;
          Cout << "dir : "<< direction << " x_O : " << position_cylindre[0] << " z_O : " << position_cylindre[2] << " V_maille : " << volume_maille << finl;
          Cout << "                               " << finl;
# endif
          for (int k = kmin; k <= kmax; k++)
            {
              // si direction==DIRECTION_I alors x_ccord = coordonnee du noeud dans la direction I donc x_noeud et z_coord = z_element
              // si direction==DIRECTION_K alors x_coord = x_element et z_coord = z_noeud
              // x_coord et z_coord sont les abcisses des vitesse qui sont situee au centre des faces
              const double z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K][k+offset_k] : elem_pos[DIRECTION_K][k+offset_k];
              for (int j = jmin; j < jmax; j++)
                {
                  //const double y_coord = (direction==DIRECTION_J) ? nodes_pos[DIRECTION_J][j+offset_j] : elem_pos[DIRECTION_J][j+offset_j];
                  for (int i = imin; i <= imax; i++)
                    {
                      const double x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I][i+offset_i] : elem_pos[DIRECTION_I][i+offset_i];
                      const double dOMx = x_coord - position_cylindre[0]; // distance du centre au point M selon x
                      const double dOMz = z_coord - position_cylindre[2]; // distance du centre au point M selon z
                      const double ABS_dOMx = fabs(dOMx);
                      const double ABS_dOMz = fabs(dOMz);
                      const double Borne_min_z = L_cube_ - delta_z/2;
                      const double Borne_max_z = L_cube_ + delta_z/2;
                      const double Borne_min_x = L_cube_ - delta_x/2;
                      const double Borne_max_x = L_cube_ + delta_x/2;

                      if (  (ABS_dOMx <= Borne_max_x)  && (ABS_dOMz <= Borne_max_z))
                        {
                          Vecteur3 rotation;
                          rotation[0] = omega * (- dOMz);
                          rotation[1] = 0.;
                          rotation[2] = omega * dOMx;
                          double extrapolated_v = vitesse_cylindre[direction] + rotation[direction];
                          double f_k;
                          double f_i;
                          double f_ik;
                          const double d_k = Borne_max_z - ABS_dOMz; // distance occupee selon z dans la maille par le cube
                          const double d_i = Borne_max_x - ABS_dOMx; // distance occupae selon x dans la maille par le cube
                          if (( ABS_dOMx <= Borne_min_x) && ( ABS_dOMz <= Borne_min_z))  // alors le forcage est total : d_OM <= R - delta_z/2
                            {
                              f_i = 1.;
                              f_k = 1.;
                            }
                          else if ( ABS_dOMx <= Borne_min_x)   // vol de controle occupe en k
                            {
                              f_i = 1.;
                              f_k = d_k / delta_z; // fraction volumique occupe
                            }
                          else if( ABS_dOMz <= Borne_min_z)    // vol de controle occupe en i
                            {
                              f_i = d_i / delta_x; // fraction volumique occupe
                              f_k = 1.;
                            }
                          else   // vol de controle occupe en i et k
                            {
                              f_i = d_i / delta_x;
                              f_k = d_k / delta_z;
                            }
                          f_ik = f_i * f_k; // frac volumique occupee

# if 0
                          Cout << "Forcage "<< "i : " << i << " j : " << j << " k : " << k <<  " x_M : " << x_coord << " z_M : " << z_coord << " v* : " ;
                          Cout << velocity(i,j,k) << " f_k : " << f_k <<" f_i : " << f_i << " f_ik : " << f_ik << finl;
#endif
                          double delta_v = f_ik * (extrapolated_v - velocity(i,j,k)); // qdm ds la zone d'etalement
                          // Calcul de rho sur la face: moyenne des rho sur les elements voisins
                          double rho_fluide = (rho_fluide_field(i,j,k) + rho_fluide_field(i+ivoisin, j+jvoisin, k+kvoisin)) * 0.5; // rho sur la face est la moy des deux cotes
                          delta_qdm_cylindre += delta_v * volume_maille * rho_fluide;
                          //# if 0
                          // Sortie pour le cas ref de validation avec translation selon z : cube de cote 22 mailles a Re10
                          // evolution de la qdm et dans deux mailles : une sortante et une entrante
                          if ( (direction==0) && ((j==0) && (i==64)) && ((k==54)|(k==55)|(k==76)|(k==77))  )
                            {
                              Cout << " " << finl;
                              Cout << "qdm_ijk" << i << j << k << " : " << delta_v * volume_maille * rho_fluide << " x_M : " << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2] + L_cube_<< " frac_vol : " << f_ik * volume_maille << finl;
                              Cout << "v_ijk" << i << j << k << " : " << velocity(i,j,k) << " x_M : "                        << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2] + L_cube_<< finl;
                              Cout << " " << finl;
                            }
                          //# endif
# if 0
                          Cout <<  "delta_v : "<< delta_v << " v_S : " << extrapolated_v << " qdm : " << delta_v * volume_maille * rho_fluide << " qdm_cube : " << delta_qdm_cylindre << " v_force : " << velocity(i, j, k) << finl;
# endif
                          velocity(i, j, k) = f_ik * extrapolated_v + (1-f_ik) * velocity(i, j, k); // forcage proportionnel au vol de controle occupe
                          masse_fluide += volume_maille * rho_fluide;
                          volume += volume_maille;
                        }
                    } // fin de la boucle en i
                }
            }

          integrale_force(itube, direction) = delta_qdm_cylindre;
          masse_fluide_cylindres(itube, direction) = masse_fluide;
          volume_cylindres(itube, direction) = volume;
# if 0
          Cout << "integrale_force(itube, direction) " << integrale_force(itube, direction) << " V_S : " << volume_cylindres(itube, direction) << finl;
# endif
        }

      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      // Affichage de la vitesse imposee
      Cout << "v_imposee : " << vx(64, 0, 76) << " z_interface : " <<  position_cylindre[2] + L_cube_ << finl;
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    } // fin boucle sur itube
  mp_sum_for_each_item(integrale_force);
  mp_sum_for_each_item(masse_fluide_cylindres);
  mp_sum_for_each_item(volume_cylindres);
}

///////////////
///// force ibc_diffuse post proj
///////////////////////


void Couplage_Tubes_IBC::ibc_diffuse_force_cube(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                                const IJK_Field_double& rho_fluide_field,
                                                DoubleTab& masse_fluide_cylindres,
                                                DoubleTab& volume_cylindres,
                                                DoubleTab& integrale_force,
                                                const Faisceau_Tubes& faisceau) const
{
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  ArrOfDouble nodes_pos[3]; //defini un tableau dynamique de 3 lignes pour la position des nodes
  ArrOfDouble elem_pos[3];  //defini un tableau dynamique de 3 lignes pour la position des elements
  {
    for (int dir = 0; dir < 3; dir++)
      {
        nodes_pos[dir] = geom.get_node_coordinates(dir); // node=noeud, les noeuds sont a l'intersection des mailles; on ne connait que la position des noeuds
        const int n = nodes_pos[dir].size_array()-1; // n= nombre de noeuds-1 = nombre d'elements dans la direction dir; les alements sont au centre des mailles
        elem_pos[dir].resize_array(n); // donne la taille du tableau contenant les elements
        for (int i = 0; i < n; i++)
          elem_pos[dir][i] = (nodes_pos[dir][i] + nodes_pos[dir][i+1]) * 0.5; // la position de l'element selon dir est moyenne de la positions des noeuds de part et d'autre
      }
  }
  const int offset_i = splitting.get_offset_local(DIRECTION_I); // donne l'offset selon la direction i
  const int offset_k = splitting.get_offset_local(DIRECTION_K);  // donne l'offset selon la direction k

  const int ntubes = faisceau.size();

  // Boucle sur les tubes du faisceau
  for (int itube = 0; itube < ntubes; itube++)
    {
      const Tube_base& tube = faisceau[itube].valeur();
      const double omega = tube.get_omega();


      const Vecteur3 position_cylindre = tube.get_current_pos();
      const Vecteur3 vitesse_cylindre = tube.get_current_velocity();

      // Boucle sur les trois composantes de qdm
      for (int direction = 0; direction < 3; direction++)   // boucle sur les 3 directions de la vitesse
        {

          const int ivoisin = (direction == 0) ? -1 : 0;
          const int jvoisin = (direction == 1) ? -1 : 0;
          const int kvoisin = (direction == 2) ? -1 : 0;

          const double volume_maille = geom.get_constant_delta(0) * geom.get_constant_delta(1) * geom.get_constant_delta(2); //donne le volume de la maille
          IJK_Field_double& velocity = select(direction, vx, vy, vz);  // velocity est egale a la vitesse selon la direction 0,1 ou 2
          int imin = 0;
          int jmin = 0;
          int kmin = 0;
          int imax = velocity.ni(); // donne le nombre de valeur de la vitesse selon i
          int jmax = velocity.nj();
          int kmax = velocity.nk();
          const double delta_x = geom.get_constant_delta(0);
          const double delta_z = geom.get_constant_delta(2);
          {
            // restreint le domaine ijk balaye a la zone couverte par le cylindre:
            const ArrOfDouble& x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I] : elem_pos[DIRECTION_I];
            const double xmin = position_cylindre[0] - (L_cube_ + delta_x); // on n'oublie pas que l'on doit chercher dans la zone courverte par le cube + une maille
            const double xmax = position_cylindre[0] + (L_cube_ + delta_x);

            while (imin < velocity.ni() && x_coord[imin + offset_i] < xmin)
              imin++;
            while (imax > 0 && x_coord[imax + offset_i - 1] > xmax)
              imax--;
          }
          {
            const ArrOfDouble& z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K] : elem_pos[DIRECTION_K];
            const double zmin = position_cylindre[2] - (L_cube_ + delta_z);
            const double zmax = position_cylindre[2] + (L_cube_ + delta_z);

            while (kmin < velocity.nk() && z_coord[kmin + offset_k] < zmin)
              kmin++;
            while (kmax > 0 && z_coord[kmax + offset_k - 1] > zmax)
              {
                kmax--;
              }
          }
          double delta_qdm_cylindre = 0.;
          double masse_fluide = 0.; // Masse de fluide dans le cylindre
          double volume = 0.;

          for (int k = kmin; k <= kmax; k++)
            {
              // si direction==DIRECTION_I alors x_ccord = coordonnee du noeud dans la direction I donc x_noeud et z_coord = z_element
              // si direction==DIRECTION_K alors x_coord = x_element et z_coord = z_noeud
              // x_coord et z_coord sont les abcisses des vitesse qui sont situee au centre des faces
              const double z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K][k+offset_k] : elem_pos[DIRECTION_K][k+offset_k];
              for (int j = jmin; j < jmax; j++)
                {
                  //const double y_coord = (direction==DIRECTION_J) ? nodes_pos[DIRECTION_J][j+offset_j] : elem_pos[DIRECTION_J][j+offset_j];
                  for (int i = imin; i <= imax; i++)
                    {
                      const double x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I][i+offset_i] : elem_pos[DIRECTION_I][i+offset_i];
                      const double dOMx = x_coord - position_cylindre[0]; // distance du centre au point M selon x
                      const double dOMz = z_coord - position_cylindre[2]; // distance du centre au point M selon z
                      const double ABS_dOMx = fabs(dOMx);
                      const double ABS_dOMz = fabs(dOMz);
                      const double Borne_min_z = L_cube_ - delta_z/2;
                      const double Borne_max_z = L_cube_ + delta_z/2;
                      const double Borne_min_x = L_cube_ - delta_x/2;
                      const double Borne_max_x = L_cube_ + delta_x/2;

                      if (  (ABS_dOMx <= Borne_max_x)  && (ABS_dOMz <= Borne_max_z))
                        {
                          Vecteur3 rotation;
                          rotation[0] = omega * (- dOMz);
                          rotation[1] = 0.;
                          rotation[2] = omega * dOMx;
                          double extrapolated_v = vitesse_cylindre[direction] + rotation[direction];
                          double f_k;
                          double f_i;
                          double f_ik;
                          const double d_k = Borne_max_z - ABS_dOMz; // distance occupee selon z dans la maille par le cube
                          const double d_i = Borne_max_x - ABS_dOMx; // distance occupae selon x dans la maille par le cube
                          if (( ABS_dOMx <= Borne_min_x) && ( ABS_dOMz <= Borne_min_z))  // alors le forcage est total : d_OM <= R - delta_z/2
                            {
                              f_i = 1.;
                              f_k = 1.;
                            }
                          else if ( ABS_dOMx <= Borne_min_x)   // vol de controle occupe en k
                            {
                              f_i = 1.;
                              f_k = d_k / delta_z; // fraction volumique occupe
                            }
                          else if( ABS_dOMz <= Borne_min_z)    // vol de controle occupe en i
                            {
                              f_i = d_i / delta_x; // fraction volumique occupe
                              f_k = 1.;
                            }
                          else   // vol de controle occupe en i et k
                            {
                              f_i = d_i / delta_x;
                              f_k = d_k / delta_z;
                            }
                          f_ik = f_i * f_k; // frac volumique occupee
                          double delta_v = f_ik * (extrapolated_v - velocity(i,j,k)); // qdm ds la zone d'etalement
                          // Calcul de rho sur la face: moyenne des rho sur les elements voisins
                          double rho_fluide = (rho_fluide_field(i,j,k) + rho_fluide_field(i+ivoisin, j+jvoisin, k+kvoisin)) * 0.5; // rho sur la face est la moy des deux cotes
                          delta_qdm_cylindre += delta_v * volume_maille * rho_fluide;
                          //# if 0
                          // Sortie pour le cas ref de validation avec translation selon z : cube de cote 22 mailles a Re10
                          // evolution de la qdm et dans deux mailles : une sortante et une entrante
                          if ( (direction==0) && ((j==0) && (i==64)) && ((k==55)|(k==76)) )
                            {
                              Cout << " " << finl;
                              Cout << "qdm_ijk_post_proj" << i << j << k << " : " << delta_v * volume_maille * rho_fluide << " x_M : " << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2] + L_cube_<< " frac_vol : " << f_ik * volume_maille << finl;
                              Cout << " " << finl;
                            }
                          //# endif
                          masse_fluide += volume_maille * rho_fluide;
                          volume += volume_maille;
                        }
                    } // fin de la boucle en i
                }
            }

          integrale_force(itube, direction) = delta_qdm_cylindre;
          masse_fluide_cylindres(itube, direction) = masse_fluide;
          volume_cylindres(itube, direction) = volume;
        }
    } // fin boucle sur itube
  mp_sum_for_each_item(integrale_force);
  mp_sum_for_each_item(masse_fluide_cylindres);
  mp_sum_for_each_item(volume_cylindres);
}



/////////////////////////
//// IBC interface localisee,forcage dans la couche limite sur une maille, Cube
////////////////////////////////

void Couplage_Tubes_IBC::ibc_localisee_velocity_cube(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                                     const IJK_Field_double& rho_fluide_field,
                                                     DoubleTab& masse_fluide_cylindres,
                                                     DoubleTab& volume_cylindres,
                                                     DoubleTab& integrale_force,
                                                     const Faisceau_Tubes& faisceau, double current_time) const
{
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  ArrOfDouble nodes_pos[3]; //defini un tableau dynamique de 3 lignes pour la position des nodes
  ArrOfDouble elem_pos[3];  //defini un tableau dynamique de 3 lignes pour la position des elements
  {
    for (int dir = 0; dir < 3; dir++)
      {
        nodes_pos[dir] = geom.get_node_coordinates(dir); // node=noeud, les noeuds sont a l'intersection des mailles; on ne connait que la position des noeuds
        const int n = nodes_pos[dir].size_array()-1; // n= nombre de noeuds-1 = nombre d'elements dans la direction dir; les alements sont au centre des mailles
        elem_pos[dir].resize_array(n); // donne la taille du tableau contenant les elements
        for (int i = 0; i < n; i++)
          elem_pos[dir][i] = (nodes_pos[dir][i] + nodes_pos[dir][i+1]) * 0.5; // la position de l'element selon dir est moyenne de la positions des noeuds de part et d'autre
      }
  }
  const int offset_i = splitting.get_offset_local(DIRECTION_I); // donne l'offset selon la direction i
  const int offset_k = splitting.get_offset_local(DIRECTION_K);  // donne l'offset selon la direction k

  const int ntubes = faisceau.size();

  // Boucle sur les tubes du faisceau
  for (int itube = 0; itube < ntubes; itube++)
    {
      const Tube_base& tube = faisceau[itube].valeur();
      const double omega = tube.get_omega();


      const Vecteur3 position_cylindre = tube.get_current_pos();
      const Vecteur3 vitesse_cylindre = tube.get_current_velocity();

      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      // Affichage de la vitesse predite
      Cout << "v_predite : " << vx(64, 0, 76) << " z_interface : " <<  position_cylindre[2] + L_cube_ << finl;
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      DoubleTab v_predite;

      // Boucle sur les trois composantes de qdm
      for (int direction = 0; direction < 3; direction++)   // boucle sur les 3 directions de la vitesse
        {

          const int ivoisin = (direction == 0) ? -1 : 0;
          const int jvoisin = (direction == 1) ? -1 : 0;
          const int kvoisin = (direction == 2) ? -1 : 0;

          const double volume_maille = geom.get_constant_delta(0) * geom.get_constant_delta(1) * geom.get_constant_delta(2); //donne le volume de la maille
          IJK_Field_double& velocity = select(direction, vx, vy, vz);  // velocity est egale a la vitesse selon la direction 0,1 ou 2
          int imin = 0;
          int jmin = 0;
          int kmin = 0;
          int imax = velocity.ni(); // donne le nombre de valeur de la vitesse selon i
          int jmax = velocity.nj();
          int kmax = velocity.nk();
          const double delta_x = geom.get_constant_delta(0);
          const double delta_z = geom.get_constant_delta(2);
          //const double d_x = geom.get_constant_delta(0);
          //const double d_z = geom.get_constant_delta(2);


          {
            // restreint le domaine ijk balaye a la zone couverte par le cylindre:
            const ArrOfDouble& x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I] : elem_pos[DIRECTION_I];
            const double xmin = position_cylindre[0] - (L_cube_ + 2 * delta_x); // on n'oublie pas que l'on doit chercher dans la zone courverte par le cube + deux mailles car pr le calcul de v ds la couche lim on a besoin de la vitesse au dessus
            const double xmax = position_cylindre[0] + (L_cube_ + 2 * delta_x);

            while (imin < velocity.ni() && x_coord[imin + offset_i] < xmin)
              imin++;
            while (imax > 0 && x_coord[imax + offset_i - 1] > xmax)
              imax--;
          }
          {
            const ArrOfDouble& z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K] : elem_pos[DIRECTION_K];
            const double zmin = position_cylindre[2] - (L_cube_ + 2 * delta_z);
            const double zmax = position_cylindre[2] + (L_cube_ + 2 * delta_z);

            while (kmin < velocity.nk() && z_coord[kmin + offset_k] < zmin)
              kmin++;
            while (kmax > 0 && z_coord[kmax + offset_k - 1] > zmax)
              {
                kmax--;
              }
          }
          double delta_qdm_cylindre = 0.;
          double masse_fluide = 0.; // Masse de fluide dans le cylindre
          double volume = 0.;
# if 0
          Cout << "                               " << finl;
          Cout << "dir : "<< direction << " x_O : " << position_cylindre[0] << " z_O : " << position_cylindre[2] << " V_maille : " << volume_maille << finl;
          Cout << "                               " << finl;
# endif

          // Calcul des bornes du tableau contenant les copies de v* dans le domaine ibc
          const int N_i = imax - imin + 1;  // imin   <= i <= imax, donc il y a imax - imin +1 valeurs de la vitesse selon i qu'on veut stocker
          const int N_j = jmax;            // jmin = 0 <=j < jmax, donc il y a jmax -1 - 0 + 1 = jmax valeurs de la vitesse selon j que l'on veut stocker
          const int N_k = kmax - kmin + 1;
          v_predite.resize(3, N_i, N_j, N_k);

          for (int k = kmin; k <= kmax; k++)
            {
              // si direction==DIRECTION_I alors x_ccord = coordonnee du noeud dans la direction I donc x_noeud et z_coord = z_element
              // si direction==DIRECTION_K alors x_coord = x_element et z_coord = z_noeud
              // x_coord et z_coord sont les abcisses des vitesse qui sont situee au centre des faces
              const double z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K][k+offset_k] : elem_pos[DIRECTION_K][k+offset_k];
              for (int j = jmin; j < jmax; j++)
                {
                  //const double y_coord = (direction==DIRECTION_J) ? nodes_pos[DIRECTION_J][j+offset_j] : elem_pos[DIRECTION_J][j+offset_j];
                  for (int i = imin; i <= imax; i++)
                    {
                      const double x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I][i+offset_i] : elem_pos[DIRECTION_I][i+offset_i];
                      const double dOMx = x_coord - position_cylindre[0];
                      const double dOMz = z_coord - position_cylindre[2];
                      const double ABS_dOMx = fabs(dOMx);
                      const double ABS_dOMz = fabs(dOMz);
                      const double Borne_min_z               = L_cube_ - delta_z / 2; //                  d_OM <= L - delta_z/2 : v~* = v_s            , qdm = rho * delta_v * V_maille
                      const double Borne_intermediaire1_z    = L_cube_              ; // L - delta_z/2 <= d_OM <= L             : v~* = v_s            , qdm = rho * delta_v * V_maille * f_k
                      const double Borne_intermediaire2_z    = L_cube_ + delta_z / 2; // L             <= d_OM <= L + delta_z/2 : v~* = v_couche_limite, qdm = rho * delta_v * V_maille * f_k
                      const double Borne_max_z               = L_cube_ + delta_z    ; // L + delta_z/2 <= d_OM <= L + delta_z   : v~* = v_couche_limite, qdm = 0
                      const double Borne_max_x = L_cube_;

                      v_predite(direction, i-imin, j, k-kmin) = velocity(i, j, k); // copie de la vitesse predite

                      if (  (ABS_dOMx <= Borne_max_x)  && (ABS_dOMz <= Borne_max_z))
                        {
                          Vecteur3 rotation;
                          rotation[0] = omega * (- dOMz);
                          rotation[1] = 0.;
                          rotation[2] = omega * dOMx;
                          double extrapolated_v = vitesse_cylindre[direction] + rotation[direction];
                          double delta_v = 0.;
                          double rho_fluide =0.;
                          double d_k; // distance occupe dans le vol de controle
                          double f_k; // frac vol occuope ds le vol de controle
                          double v_k_exterieure; // necessaire au calcul de v_couche_limite
                          double v_couche_limite;
                          double d_Cube_M; // distance entre le cube et le point M, necessaire pr le calcul de v_couche_limite

                          if ( ABS_dOMz <= Borne_min_z)  // d_OM <= L - delta_z/2 : v~* = v_s, qdm = rho * delta_v * V_maille
                            {
                              delta_v = extrapolated_v - velocity(i,j,k);
                              rho_fluide = (rho_fluide_field(i,j,k) + rho_fluide_field(i+ivoisin, j+jvoisin, k+kvoisin)) * 0.5; // Calcul de rho sur la face: moyenne des rho sur les elements voisins (moy des deux cotes)
                              delta_qdm_cylindre += delta_v * volume_maille * rho_fluide;

# if 0
                              Cout << "Forcage_complet "<< "i : " << i << " j : " << j << " k : " << k <<  " x_M : " << x_coord << " z_M : " << z_coord << " v*(i,j,k) : " << velocity(i,j,k) << finl;
#endif
                              //# if 0
                              // Sortie pour le cas ref de validation avec translation selon z : cube de cote 22 mailles a Re10
                              // evolution de la qdm et de la vitesse dans deux mailles : une sortante et une entrante
                              if ( (direction==0) && ((j==0) && (i==64)) && ((k==55)|(k==76)) )
                                {
                                  Cout << " " << finl;
                                  Cout << "qdm_ijk" << i << j << k << " : " << delta_v * volume_maille * rho_fluide << " x_M : " << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2]+L_cube_ << " forcage_complet " << finl;
                                  Cout << "v_ijk" << i << j << k << " : " << velocity(i,j,k) << " x_M : "                        << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2]+L_cube_ << " forcage_complet " << finl;
                                  Cout << " " << finl;
                                }
                              //# endif

                              velocity(i, j, k) = extrapolated_v; // v~* = v_S car d_OM <=L
                              masse_fluide += volume_maille * rho_fluide;
                              volume += volume_maille;

# if 0
                              Cout <<  " delta_v : "<< delta_v << " v_S : " << extrapolated_v << " qdm : " << delta_v * volume_maille * rho_fluide << " qdm_cube : " << delta_qdm_cylindre << " v~*(i,j,k) : " << velocity(i, j, k) << finl;
# endif
                            }
                          else if ( ABS_dOMz <= Borne_intermediaire1_z)  // L - delta_z/2 <= d_OM <= L : v~* = v_s, qdm = rho * delta_v * V_maille * f_k
                            {
                              d_k = L_cube_ + delta_z / 2 - ABS_dOMz; // distance occupee selon z dans le vol de controle par le cube
                              f_k = d_k / delta_z; // fraction volumique occupe selon z, on a :     0.5 <= f_k <=1
                              delta_v = f_k * (extrapolated_v - velocity(i,j,k)); // qdm prop au vol occupe dans le vol de controle
                              rho_fluide = (rho_fluide_field(i,j,k) + rho_fluide_field(i+ivoisin, j+jvoisin, k+kvoisin)) * 0.5; // Calcul de rho sur la face: moyenne des rho sur les elements voisins (moy des deux cotes)
                              delta_qdm_cylindre += delta_v * volume_maille * rho_fluide;

# if 0
                              Cout << "Forcage_complet_qdm_partielle "<< "i : " << i << " j : " << j << " k : " << k;
                              Cout <<  " x_M : " << x_coord << " z_M : " << z_coord << " v*(i,j,k) : " << velocity(i,j,k) << " f_k : " << f_k << finl;
#endif
                              //# if 0
                              // Sortie pour le cas ref de validation avec translation selon z : cube de cote 22 mailles a Re10
                              // evolution de la qdm dans deux mailles : une sortante et une entrante
                              if ( (direction==0) && ((j==0) && (i==64)) && ((k==55)|(k==76)) )
                                {
                                  Cout << " " << finl;
                                  Cout << "Forcage_complet_qdm_partielle "<< " v*(i,j,k) : " << velocity(i,j,k) << " f_k : " << f_k << " v_S : " << extrapolated_v << finl;
                                  Cout << "qdm_ijk" << i << j << k << " : " << delta_v * volume_maille * rho_fluide << " x_M : " << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2]+L_cube_ << finl;
                                  Cout << "v_ijk" << i << j << k << " : " << velocity(i,j,k) << " x_M : "                        << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2]+L_cube_ << finl;
                                  Cout << " " << finl;
                                }
                              //# endif

                              velocity(i, j, k) = extrapolated_v; // v~* = v_S car d_OM <=L
                              masse_fluide += volume_maille * rho_fluide;
                              volume += volume_maille;

# if 0
                              Cout <<  " delta_v : "<< delta_v << " v_S : " << extrapolated_v << " qdm : " << delta_v * volume_maille * rho_fluide << " qdm_cube : " << delta_qdm_cylindre << " v~*(i,j,k) : " << velocity(i, j, k) << finl;
# endif
                            }
                          else if ( ABS_dOMz <= Borne_intermediaire2_z)    // L<= d_OM <= L + delta_z/2 : v~* = v_couche_limite, qdm = rho * delta_v * V_maille * f_k
                            {
                              d_k = L_cube_ + delta_z / 2 - ABS_dOMz; // distance occupee selon z dans le vol de controle par le cube
                              f_k = d_k / delta_z; // fraction volumique occupe selon z, on a :     0. <= f_k <= 0.5
                              double v_k_interieure; // necessaire au calcul de v_predite_interpollee pr le calcul de la qdm
                              if ( dOMz > 0)
                                {
                                  v_k_exterieure  = velocity(i, j, k+1);
                                  v_k_interieure  = v_predite(direction, i-imin, j, k- kmin - 1); // la vitesse predite en k-1 a ete force a v_s precedemment on doit donc prendre la copie stockee dans v_predite, pour i=imin on ne peut pas etre dans ce else if

                                }
                              else
                                {
                                  v_k_exterieure  = velocity(i, j, k-1);
                                  v_k_interieure  = v_predite(direction, i-imin, j, k- kmin + 1); // pour k=kmax, on ne peut pas se retrouver dans ce else if
                                }
                              double v_predite_interpollee_sur_la_frontiere = ((d_k + delta_z / 2) / delta_z) * ( velocity(i,j,k) - v_k_interieure ) + v_k_interieure; // On calcule un v* par interpol lineaire au niveau de la frontiere du cube
                              delta_v = f_k * (extrapolated_v - v_predite_interpollee_sur_la_frontiere); // qdm prop au vol occupe dans le vol de controle
                              rho_fluide = (rho_fluide_field(i,j,k) + rho_fluide_field(i+ivoisin, j+jvoisin, k+kvoisin)) * 0.5; // Calcul de rho sur la face: moyenne des rho sur les elements voisins (moy des deux cotes)
                              delta_qdm_cylindre += delta_v * volume_maille * rho_fluide;

# if 0
                              Cout << "Forcage_couche_limite_qdm_partielle "<< "i : " << i << " j : " << j << " k : " << k <<  " x_M : " << x_coord << " z_M : " << z_coord << " v*_int :" << v_k_interieure << " v*_S : ";
                              Cout<< v_predite_interpollee_sur_la_frontiere << "  v*(i,j,k) : " << velocity(i,j,k) << " v*_ext :" << v_k_exterieure << " f_k : " << f_k << " d_Cube_M : " << d_Cube_M << finl;
#endif
                              d_Cube_M = ABS_dOMz - L_cube_; // distance entre M et le cube,    0 <= d_Cube_M <= delta_z / 2
                              v_couche_limite = (d_Cube_M / (d_Cube_M + delta_z )) * ( v_k_exterieure - extrapolated_v ) + extrapolated_v;

                              //# if 0
                              // Sortie pour le cas ref de validation avec translation selon z : cube de cote 22 mailles a Re10
                              // evolution de la qdm dans deux mailles : une sortante et une entrante
                              if ( (direction==0) && ((j==0) && (i==64)) && ((k==55)|(k==76)) )
                                {
                                  Cout << " " << finl;
                                  Cout << "Forcage_couche_limite_qdm_partielle "<< " v*_int :" << v_k_interieure << " v*_S : " << v_predite_interpollee_sur_la_frontiere << "  v*(i,j,k) : " << velocity(i,j,k) << " v*_ext :" << v_k_exterieure;
                                  Cout << " f_k : " << f_k << " v_couche_limite : " << v_couche_limite << " v_S : " << extrapolated_v << " d_Cube_M : " << d_Cube_M << finl;
                                  Cout << "qdm_ijk" << i << j << k << " : " << delta_v * volume_maille * rho_fluide << " x_M : " << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2]+L_cube_ << finl;
                                  Cout << "v_ijk" << i << j << k << " : " << velocity(i, j, k)                      << " x_M : " << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2]+L_cube_ << finl;
                                  Cout << " " << finl;
                                }
                              //# endif

                              velocity(i, j, k) = v_couche_limite; // v~* = v_couche_limite car  L<= d_OM <= L + delta_z/2
                              masse_fluide += volume_maille * rho_fluide;
                              volume += volume_maille;

# if 0
                              Cout <<  " delta_v : "<< delta_v << " v_S : " << extrapolated_v << " v_couche_limite : " << v_couche_limite << " qdm : " << delta_v * volume_maille * rho_fluide << " qdm_cube : " << delta_qdm_cylindre << " v~*(i,j,k) : " << velocity(i, j, k) << finl;
# endif
                            }

                          else    // L + delta_z/2 <= d_OM <= L + delta_z   : v~* = v_couche_limite, qdm = 0
                            {
                              d_Cube_M = ABS_dOMz - L_cube_;
                              if ( dOMz > 0)
                                {
                                  v_k_exterieure  = velocity(i, j, k+1);
                                }
                              else
                                {
                                  v_k_exterieure  = velocity(i, j, k-1);
                                }

# if 0
                              Cout << "Forcage_couche_limite_qdm_0 "<< "i : " << i << " j : " << j << " k : " << k <<  " x_M : " << x_coord << " z_M : " << z_coord;
                              Cout << "  v*(i,j,k) : " << velocity(i,j,k) << " v*_ext :" << v_k_exterieure << "  d_Cube_M: " << d_Cube_M << finl;
#endif
                              v_couche_limite = (d_Cube_M / (d_Cube_M + delta_z )) * ( v_k_exterieure - extrapolated_v ) + extrapolated_v;

                              //# if 0
                              // Sortie pour le cas ref de validation avec translation selon z : cube de cote 22 mailles a Re10
                              // evolution de la qdm dans deux mailles : une sortante et une entrante
                              if ( (direction==0) && ((j==0) && (i==64)) && ((k==55)|(k==76)) )
                                {
                                  Cout << " " << finl;
                                  Cout << "Forcage_couche_limite_qdm0 "<< "  v*(i,j,k) : " << velocity(i,j,k) << " v*_ext :" << v_k_exterieure << " v_couche_limite : " << v_couche_limite << " v_S : " << extrapolated_v << " d_Cube_M : " << d_Cube_M << finl;
                                  Cout << "qdm_ijk" << i << j << k << " : " << delta_v * volume_maille * rho_fluide << " x_M : " << x_coord << " z_M : " << z_coord << " z_O : " << position_cylindre[2] << finl;
                                  Cout << "v_ijk" << i << j << k << " : " << velocity(i, j, k)                      << " x_M : " << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2]+L_cube_ << finl;
                                  Cout << " " << finl;
                                }
                              //# endif

                              velocity(i, j, k) = v_couche_limite; // v~* = v_couche_limite car  L + delta_z/2<= d_OM <= L
                              masse_fluide += volume_maille * rho_fluide;
                              volume += volume_maille;
# if 0
                              Cout << " v_S : " << extrapolated_v << " v_couche_limite : " << v_couche_limite << " qdm_cube : " << delta_qdm_cylindre << " v~*(i,j,k) : " << velocity(i, j, k) << finl;
# endif
                            }
                        } // fin de la boucle sur les 4 tests
                    } // fin de la boucle en i
                }
            }

          integrale_force(itube, direction) = delta_qdm_cylindre;
          masse_fluide_cylindres(itube, direction) = masse_fluide;
          volume_cylindres(itube, direction) = volume;
# if 0
          Cout << "integrale_force(itube, direction) " << integrale_force(itube, direction) << " V_S : " << volume_cylindres(itube, direction) << finl;
# endif
        }
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      // Affichage de la vitesse imposee
      Cout << "v_imposee : " << vx(64, 0, 76) << " z_interface : " <<  position_cylindre[2] + L_cube_ << finl;
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    }
  mp_sum_for_each_item(integrale_force);
  mp_sum_for_each_item(masse_fluide_cylindres);
  mp_sum_for_each_item(volume_cylindres);
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Variante calcul qdm
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/////////////////////////
//// IBC interface localisee,forcage dans la couche limite sur une maille, Cube
////////////////////////////////

void Couplage_Tubes_IBC::ibc_localisee_velocity_cube_qdm(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                                         const IJK_Field_double& rho_fluide_field,
                                                         DoubleTab& masse_fluide_cylindres,
                                                         DoubleTab& volume_cylindres,
                                                         DoubleTab& integrale_force,
                                                         const Faisceau_Tubes& faisceau, double current_time) const
{
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  ArrOfDouble nodes_pos[3]; //defini un tableau dynamique de 3 lignes pour la position des nodes
  ArrOfDouble elem_pos[3];  //defini un tableau dynamique de 3 lignes pour la position des elements
  {
    for (int dir = 0; dir < 3; dir++)
      {
        nodes_pos[dir] = geom.get_node_coordinates(dir); // node=noeud, les noeuds sont a l'intersection des mailles; on ne connait que la position des noeuds
        const int n = nodes_pos[dir].size_array()-1; // n= nombre de noeuds-1 = nombre d'elements dans la direction dir; les alements sont au centre des mailles
        elem_pos[dir].resize_array(n); // donne la taille du tableau contenant les elements
        for (int i = 0; i < n; i++)
          elem_pos[dir][i] = (nodes_pos[dir][i] + nodes_pos[dir][i+1]) * 0.5; // la position de l'element selon dir est moyenne de la positions des noeuds de part et d'autre
      }
  }
  const int offset_i = splitting.get_offset_local(DIRECTION_I); // donne l'offset selon la direction i
  const int offset_k = splitting.get_offset_local(DIRECTION_K);  // donne l'offset selon la direction k

  const int ntubes = faisceau.size();

  // Boucle sur les tubes du faisceau
  for (int itube = 0; itube < ntubes; itube++)
    {
      const Tube_base& tube = faisceau[itube].valeur();
      const double omega = tube.get_omega();


      const Vecteur3 position_cylindre = tube.get_current_pos();
      const Vecteur3 vitesse_cylindre = tube.get_current_velocity();

      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      // Affichage de la vitesse predite
      Cout << "v_predite : " << vx(64, 0, 76) << " z_interface : " <<  position_cylindre[2] + L_cube_ << finl;
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      DoubleTab v_predite;

      // Boucle sur les trois composantes de qdm
      for (int direction = 0; direction < 3; direction++)   // boucle sur les 3 directions de la vitesse
        {

          const int ivoisin = (direction == 0) ? -1 : 0;
          const int jvoisin = (direction == 1) ? -1 : 0;
          const int kvoisin = (direction == 2) ? -1 : 0;

          const double volume_maille = geom.get_constant_delta(0) * geom.get_constant_delta(1) * geom.get_constant_delta(2); //donne le volume de la maille
          IJK_Field_double& velocity = select(direction, vx, vy, vz);  // velocity est egale a la vitesse selon la direction 0,1 ou 2
          int imin = 0;
          int jmin = 0;
          int kmin = 0;
          int imax = velocity.ni(); // donne le nombre de valeur de la vitesse selon i
          int jmax = velocity.nj();
          int kmax = velocity.nk();
          const double delta_x = geom.get_constant_delta(0);
          const double delta_z = geom.get_constant_delta(2);
          //const double d_x = geom.get_constant_delta(0);
          //const double d_z = geom.get_constant_delta(2);


          {
            // restreint le domaine ijk balaye a la zone couverte par le cylindre:
            const ArrOfDouble& x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I] : elem_pos[DIRECTION_I];
            const double xmin = position_cylindre[0] - (L_cube_ + 2 * delta_x); // on n'oublie pas que l'on doit chercher dans la zone courverte par le cube + deux mailles car pr le calcul de v ds la couche lim on a besoin de la vitesse au dessus
            const double xmax = position_cylindre[0] + (L_cube_ + 2 * delta_x);

            while (imin < velocity.ni() && x_coord[imin + offset_i] < xmin)
              imin++;
            while (imax > 0 && x_coord[imax + offset_i - 1] > xmax)
              imax--;
          }
          {
            const ArrOfDouble& z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K] : elem_pos[DIRECTION_K];
            const double zmin = position_cylindre[2] - (L_cube_ + 2 * delta_z);
            const double zmax = position_cylindre[2] + (L_cube_ + 2 * delta_z);

            while (kmin < velocity.nk() && z_coord[kmin + offset_k] < zmin)
              kmin++;
            while (kmax > 0 && z_coord[kmax + offset_k - 1] > zmax)
              {
                kmax--;
              }
          }
          double delta_qdm_cylindre = 0.;
          double masse_fluide = 0.; // Masse de fluide dans le cylindre
          double volume = 0.;
# if 0
          Cout << "                               " << finl;
          Cout << "dir : "<< direction << " x_O : " << position_cylindre[0] << " z_O : " << position_cylindre[2] << " V_maille : " << volume_maille << finl;
          Cout << "                               " << finl;
# endif

          // Calcul des bornes du tableau contenant les copies de v* dans le domaine ibc
          const int N_i = imax - imin + 1;  // imin   <= i <= imax, donc il y a imax - imin +1 valeurs de la vitesse selon i qu'on veut stocker
          const int N_j = jmax;            // jmin = 0 <=j < jmax, donc il y a jmax -1 - 0 + 1 = jmax valeurs de la vitesse selon j que l'on veut stocker
          const int N_k = kmax - kmin + 1;
          v_predite.resize(3, N_i, N_j, N_k);


          double qdm_76=0;
          double qdm_77=0;
          double qdm_78=0;
          double qdm_tot=0;
          for (int k = kmin; k <= kmax; k++)
            {
              // si direction==DIRECTION_I alors x_ccord = coordonnee du noeud dans la direction I donc x_noeud et z_coord = z_element
              // si direction==DIRECTION_K alors x_coord = x_element et z_coord = z_noeud
              // x_coord et z_coord sont les abcisses des vitesse qui sont situee au centre des faces
              const double z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K][k+offset_k] : elem_pos[DIRECTION_K][k+offset_k];
              for (int j = jmin; j < jmax; j++)
                {
                  //const double y_coord = (direction==DIRECTION_J) ? nodes_pos[DIRECTION_J][j+offset_j] : elem_pos[DIRECTION_J][j+offset_j];
                  for (int i = imin; i <= imax; i++)
                    {
                      const double x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I][i+offset_i] : elem_pos[DIRECTION_I][i+offset_i];
                      const double dOMx = x_coord - position_cylindre[0];
                      const double dOMz = z_coord - position_cylindre[2];
                      const double ABS_dOMx = fabs(dOMx);
                      const double ABS_dOMz = fabs(dOMz);
                      const double Borne_min_z               = L_cube_ - delta_z / 2; //                  d_OM <= L - delta_z/2 : v~* = v_s            , qdm = rho * delta_v * V_maille
                      const double Borne_intermediaire1_z    = L_cube_              ; // L - delta_z/2 <= d_OM <= L             : v~* = v_s            , qdm = rho * delta_v * V_maille * f_k
                      const double Borne_intermediaire2_z    = L_cube_ + delta_z / 2; // L             <= d_OM <= L + delta_z/2 : v~* = v_couche_limite, qdm = rho * delta_v * V_maille * f_k
                      const double Borne_max_z               = L_cube_ + delta_z    ; // L + delta_z/2 <= d_OM <= L + delta_z   : v~* = v_couche_limite, qdm = 0
                      const double Borne_max_x = L_cube_;

                      v_predite(direction, i-imin, j, k-kmin) = velocity(i, j, k); // copie de la vitesse predite

                      if (  (ABS_dOMx <= Borne_max_x)  && (ABS_dOMz <= Borne_max_z))
                        {
                          Vecteur3 rotation;
                          rotation[0] = omega * (- dOMz);
                          rotation[1] = 0.;
                          rotation[2] = omega * dOMx;
                          double extrapolated_v = vitesse_cylindre[direction] + rotation[direction];
                          double delta_v = 0.;
                          double rho_fluide =0.;
                          double d_k; // distance occupe dans le vol de controle
                          double f_k; // frac vol occuope ds le vol de controle
                          double v_k_exterieure; // necessaire au calcul de v_couche_limite
                          double v_couche_limite;
                          double d_Cube_M; // distance entre le cube et le point M, necessaire pr le calcul de v_couche_limite

                          if ( ABS_dOMz <= Borne_min_z)  // d_OM <= L - delta_z/2 : v~* = v_s, qdm = rho * delta_v * V_maille
                            {
                              delta_v = extrapolated_v - velocity(i,j,k);
                              rho_fluide = (rho_fluide_field(i,j,k) + rho_fluide_field(i+ivoisin, j+jvoisin, k+kvoisin)) * 0.5; // Calcul de rho sur la face: moyenne des rho sur les elements voisins (moy des deux cotes)
                              delta_qdm_cylindre += delta_v * volume_maille * rho_fluide;

# if 0
                              Cout << "Forcage_complet "<< "i : " << i << " j : " << j << " k : " << k <<  " x_M : " << x_coord << " z_M : " << z_coord << " v*(i,j,k) : " << velocity(i,j,k) << finl;
#endif
                              //# if 0
                              // Sortie pour le cas ref de validation avec translation selon z : cube de cote 22 mailles a Re10
                              // evolution de la qdm et de la vitesse dans deux mailles : une sortante et une entrante
                              if ( (direction==0) && ((j==0) && (i==64)) && ((k==55)|(k==76)|(k==77)|(k==78)) )
                                {
                                  Cout << " " << finl;
                                  if (k==76)
                                    {
                                      qdm_76 = delta_v* volume_maille * rho_fluide;
                                    }
                                  else if (k==77)
                                    {
                                      qdm_77 = delta_v* volume_maille * rho_fluide;
                                    }
                                  else if (k==78)
                                    {
                                      qdm_78 = delta_v* volume_maille * rho_fluide;
                                    }
                                  Cout << "qdm_ijk" << i << j << k << " : " << delta_v * volume_maille * rho_fluide << " x_M : " << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2]+L_cube_ << " forcage_complet " << finl;
                                  Cout << "v_ijk" << i << j << k << " : " << velocity(i,j,k) << " x_M : "                        << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2]+L_cube_ << " forcage_complet " << finl;
                                  Cout << " " << finl;
                                }
                              //# endif

                              velocity(i, j, k) = extrapolated_v; // v~* = v_S car d_OM <=L
                              masse_fluide += volume_maille * rho_fluide;
                              volume += volume_maille;

# if 0
                              Cout <<  " delta_v : "<< delta_v << " v_S : " << extrapolated_v << " qdm : " << delta_v * volume_maille * rho_fluide << " qdm_cube : " << delta_qdm_cylindre << " v~*(i,j,k) : " << velocity(i, j, k) << finl;
# endif
                            }
                          else if ( ABS_dOMz <= Borne_intermediaire1_z)  // L - delta_z/2 <= d_OM <= L : v~* = v_s, qdm = rho * delta_v * V_maille * f_k
                            {
                              d_k = L_cube_ + delta_z / 2 - ABS_dOMz; // distance occupee selon z dans le vol de controle par le cube
                              f_k = d_k / delta_z; // fraction volumique occupe selon z, on a :     0.5 <= f_k <=1
                              delta_v = (extrapolated_v - velocity(i,j,k)); // qdm totale
                              rho_fluide = (rho_fluide_field(i,j,k) + rho_fluide_field(i+ivoisin, j+jvoisin, k+kvoisin)) * 0.5; // Calcul de rho sur la face: moyenne des rho sur les elements voisins (moy des deux cotes)
                              delta_qdm_cylindre += delta_v * volume_maille * rho_fluide;

# if 0
                              Cout << "Forcage_complet_qdm_partielle "<< "i : " << i << " j : " << j << " k : " << k;
                              Cout <<  " x_M : " << x_coord << " z_M : " << z_coord << " v*(i,j,k) : " << velocity(i,j,k) << " f_k : " << f_k << finl;
#endif
                              //# if 0
                              // Sortie pour le cas ref de validation avec translation selon z : cube de cote 22 mailles a Re10
                              // evolution de la qdm dans deux mailles : une sortante et une entrante
                              if ( (direction==0) && ((j==0) && (i==64)) && ((k==55)|(k==76)|(k==77)|(k==78)) )
                                {
                                  Cout << " " << finl;
                                  if (k==76)
                                    {
                                      qdm_76 = delta_v* volume_maille * rho_fluide;
                                    }
                                  else if (k==77)
                                    {
                                      qdm_77 = delta_v* volume_maille * rho_fluide;
                                    }
                                  else if (k==78)
                                    {
                                      qdm_78 = delta_v* volume_maille * rho_fluide;
                                    }
                                  Cout << "Forcage_complet_qdm_partielle "<< " v*(i,j,k) : " << velocity(i,j,k) << " f_k : " << f_k << " v_S : " << extrapolated_v << finl;
                                  Cout << "qdm_ijk" << i << j << k << " : " << delta_v * volume_maille * rho_fluide << " x_M : " << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2]+L_cube_ << finl;
                                  Cout << "v_ijk" << i << j << k << " : " << velocity(i,j,k) << " x_M : "                        << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2]+L_cube_ << finl;
                                  Cout << " " << finl;
                                }
                              //# endif

                              velocity(i, j, k) = extrapolated_v; // v~* = v_S car d_OM <=L
                              masse_fluide += volume_maille * rho_fluide;
                              volume += volume_maille;

# if 0
                              Cout <<  " delta_v : "<< delta_v << " v_S : " << extrapolated_v << " qdm : " << delta_v * volume_maille * rho_fluide << " qdm_cube : " << delta_qdm_cylindre << " v~*(i,j,k) : " << velocity(i, j, k) << finl;
# endif
                            }
                          else if ( ABS_dOMz <= Borne_intermediaire2_z)    // L<= d_OM <= L + delta_z/2 : v~* = v_couche_limite, qdm = rho * delta_v * V_maille * f_k
                            {
                              d_k = L_cube_ + delta_z / 2 - ABS_dOMz; // distance occupee selon z dans le vol de controle par le cube
                              f_k = d_k / delta_z; // fraction volumique occupe selon z, on a :     0. <= f_k <= 0.5

                              if ( dOMz > 0)
                                {
                                  v_k_exterieure  = velocity(i, j, k+1);
                                }
                              else
                                {
                                  v_k_exterieure  = velocity(i, j, k-1);
                                }
# if 0
                              Cout << "Forcage_couche_limite_qdm_partielle "<< "i : " << i << " j : " << j << " k : " << k <<  " x_M : " << x_coord << " z_M : " << z_coord << " v*_int :" << v_k_interieure << " v*_S : ";
                              Cout << v_predite_interpollee_sur_la_frontiere << "  v*(i,j,k) : " << velocity(i,j,k) << " v*_ext :" << v_k_exterieure << " f_k : " << f_k << " d_Cube_M : " << d_Cube_M << finl;
#endif
                              d_Cube_M = ABS_dOMz - L_cube_; // distance entre M et le cube,    0 <= d_Cube_M <= delta_z / 2
                              v_couche_limite = (d_Cube_M / (d_Cube_M + delta_z )) * ( v_k_exterieure - extrapolated_v ) + extrapolated_v;
                              delta_v = (v_couche_limite -velocity(i, j, k) ); // qdm totale
                              rho_fluide = (rho_fluide_field(i,j,k) + rho_fluide_field(i+ivoisin, j+jvoisin, k+kvoisin)) * 0.5; // Calcul de rho sur la face: moyenne des rho sur les elements voisins (moy des deux cotes)
                              delta_qdm_cylindre += delta_v * volume_maille * rho_fluide;

                              //# if 0
                              // Sortie pour le cas ref de validation avec translation selon z : cube de cote 22 mailles a Re10
                              // evolution de la qdm dans deux mailles : une sortante et une entrante
                              if ( (direction==0) && ((j==0) && (i==64)) && ((k==55)|(k==76)|(k==77)|(k==78)) )
                                {
                                  Cout << " " << finl;
                                  if (k==76)
                                    {
                                      qdm_76 = delta_v* volume_maille * rho_fluide;
                                    }
                                  else if (k==77)
                                    {
                                      qdm_77 = delta_v* volume_maille * rho_fluide;
                                    }
                                  else if (k==78)
                                    {
                                      qdm_78 = delta_v* volume_maille * rho_fluide;
                                    }
                                  Cout << "Forcage_couche_limite_qdm_partielle "<< "  v*(i,j,k) : " << velocity(i,j,k) << " v*_ext :" << v_k_exterieure << " f_k : " << f_k << " v_couche_limite : " << v_couche_limite << " v_S : " << extrapolated_v << " d_Cube_M : " << d_Cube_M << finl;
                                  Cout << "qdm_ijk" << i << j << k << " : " << delta_v * volume_maille * rho_fluide << " x_M : " << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2]+L_cube_ << finl;
                                  Cout << "v_ijk" << i << j << k << " : " << velocity(i, j, k)                      << " x_M : " << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2]+L_cube_ << finl;
                                  Cout << " " << finl;
                                }
                              //# endif

                              velocity(i, j, k) = v_couche_limite; // v~* = v_couche_limite car  L<= d_OM <= L + delta_z/2
                              masse_fluide += volume_maille * rho_fluide;
                              volume += volume_maille;

# if 0
                              Cout <<  " delta_v : "<< delta_v << " v_S : " << extrapolated_v << " v_couche_limite : " << v_couche_limite << " qdm : " << delta_v * volume_maille * rho_fluide << " qdm_cube : " << delta_qdm_cylindre << " v~*(i,j,k) : " << velocity(i, j, k) << finl;
# endif
                            }

                          else    // L + delta_z/2 <= d_OM <= L + delta_z   : v~* = v_couche_limite, qdm = 0
                            {
                              d_Cube_M = ABS_dOMz - L_cube_;
                              if ( dOMz > 0)
                                {
                                  v_k_exterieure  = velocity(i, j, k+1);
                                }
                              else
                                {
                                  v_k_exterieure  = velocity(i, j, k-1);
                                }

# if 0
                              Cout << "Forcage_couche_limite_qdm_0 "<< "i : " << i << " j : " << j << " k : " << k;
                              Cout <<  " x_M : " << x_coord << " z_M : " << z_coord << "  v*(i,j,k) : " << velocity(i,j,k) << " v*_ext :" << v_k_exterieure << "  d_Cube_M: " << d_Cube_M << finl;
#endif
                              v_couche_limite = (d_Cube_M / (d_Cube_M + delta_z )) * ( v_k_exterieure - extrapolated_v ) + extrapolated_v;
                              delta_v = (v_couche_limite - velocity(i, j, k)); // qdm totale
                              rho_fluide = (rho_fluide_field(i,j,k) + rho_fluide_field(i+ivoisin, j+jvoisin, k+kvoisin)) * 0.5; // Calcul de rho sur la face: moyenne des rho sur les elements voisins (moy des deux cotes)
                              delta_qdm_cylindre += delta_v * volume_maille * rho_fluide;

                              //# if 0
                              // Sortie pour le cas ref de validation avec translation selon z : cube de cote 22 mailles a Re10
                              // evolution de la qdm dans deux mailles : une sortante et une entrante
                              if ( (direction==0) && ((j==0) && (i==64)) && ((k==55)|(k==76)|(k==77)|(k==78)) )
                                {
                                  Cout << " " << finl;
                                  if (k==76)
                                    {
                                      qdm_76 = delta_v* volume_maille * rho_fluide;
                                    }
                                  else if (k==77)
                                    {
                                      qdm_77 = delta_v* volume_maille * rho_fluide;
                                    }
                                  else if (k==78)
                                    {
                                      qdm_78 = delta_v* volume_maille * rho_fluide;
                                    }
                                  Cout << "Forcage_couche_limite_qdm0 "<< "  v*(i,j,k) : " << velocity(i,j,k) << " v*_ext :" << v_k_exterieure << " v_couche_limite : " << v_couche_limite << " v_S : " << extrapolated_v << " d_Cube_M : " << d_Cube_M << finl;
                                  Cout << "qdm_ijk" << i << j << k << " : " << delta_v * volume_maille * rho_fluide << " x_M : " << x_coord << " z_M : " << z_coord << " z_O : " << position_cylindre[2]+L_cube_ << finl;
                                  Cout << "v_ijk" << i << j << k << " : " << velocity(i, j, k)                      << " x_M : " << x_coord << " z_M : " << z_coord << " z_interface : " << position_cylindre[2]+L_cube_ << finl;
                                  Cout << " " << finl;
                                }
                              //# endif

                              velocity(i, j, k) = v_couche_limite; // v~* = v_couche_limite car  L + delta_z/2<= d_OM <= L
                              masse_fluide += volume_maille * rho_fluide;
                              volume += volume_maille;
# if 0
                              Cout << " v_S : " << extrapolated_v << " v_couche_limite : " << v_couche_limite << " qdm_cube : " << delta_qdm_cylindre << " v~*(i,j,k) : " << velocity(i, j, k) << finl;
# endif
                            }
                        } // fin de la boucle sur les 4 tests
                    } // fin de la boucle en i
                }
            }
          if (direction==0)
            {
              qdm_tot = qdm_76 + qdm_77 + qdm_78;
              Cout << "qdm_tot : " << qdm_tot << " z_interface : " << position_cylindre[2]+L_cube_ << finl;
            }
          integrale_force(itube, direction) = delta_qdm_cylindre;
          masse_fluide_cylindres(itube, direction) = masse_fluide;
          volume_cylindres(itube, direction) = volume;
# if 0
          Cout << "integrale_force(itube, direction) " << integrale_force(itube, direction) << " V_S : " << volume_cylindres(itube, direction) << finl;
# endif
        }
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      // Affichage de la vitesse imposee
      Cout << "v_imposee : " << vx(64, 0, 76) << " z_interface : " <<  position_cylindre[2] + L_cube_ << finl;
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    }
  mp_sum_for_each_item(integrale_force);
  mp_sum_for_each_item(masse_fluide_cylindres);
  mp_sum_for_each_item(volume_cylindres);
}


void interpoler_vitesses_xz(const DoubleTab& coords, DoubleTab& vitesses_interpolees,
                            const IJK_Field_double& vx,
                            const IJK_Field_double& vz,
                            ArrOfDouble& result)
{
  int i;
  const int n = coords.dimension(0);
  vitesses_interpolees.resize(n,3);
  ijk_interpolate(vx, coords, result);
  for (i = 0; i < n; i++)
    vitesses_interpolees(i,0)=result[i];
  for (i = 0; i < n; i++)
    vitesses_interpolees(i,1)=0.;
  ijk_interpolate(vz, coords, result);
  for (i = 0; i < n; i++)
    vitesses_interpolees(i,2)=result[i];
}


//////////////////////////////////////////////////////////////////////////
/// fonction_ibc pour faire le champ miroir
//////////////////////////////////////////////////////////////////////////
void Couplage_Tubes_IBC::force_ibc_velocity_symetrie_plane(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                                           const IJK_Field_double& rho_fluide_field,
                                                           DoubleTab& masse_fluide_cylindres,
                                                           DoubleTab& volume_cylindres,
                                                           DoubleTab& integrale_force,
                                                           const Faisceau_Tubes& faisceau) const
{
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  ArrOfDouble nodes_pos[3]; //defini un tableau dynamique de 3 lignes pour la position des nodes
  ArrOfDouble elem_pos[3];  //defini un tableau dynamique de 3 lignes pour la position des elements
  {
    for (int dir = 0; dir < 3; dir++)
      {
        nodes_pos[dir] = geom.get_node_coordinates(dir); // node=noeud, les noeuds sont a l'intersection des mailles; on ne connait que la position des noeuds
        const int n = nodes_pos[dir].size_array()-1; // n= nombre de noeuds-1 = nombre d'elements dans la direction dir; les alements sont au centre des mailles
        elem_pos[dir].resize_array(n); // donne la taille du tableau contenant les elements
        for (int i = 0; i < n; i++)
          elem_pos[dir][i] = (nodes_pos[dir][i] + nodes_pos[dir][i+1]) * 0.5; // la position de l'element selon dir est moyenne de la positions des noeuds de part et d'autre
      }
  }
  const int offset_i = splitting.get_offset_local(DIRECTION_I); // donne l'offset selon la direction i
  const int offset_j = splitting.get_offset_local(DIRECTION_J);  // donne l'offset selon la direction k
  const int offset_k = splitting.get_offset_local(DIRECTION_K);  // donne l'offset selon la direction k

  const int ntubes = faisceau.size();
  const double h_maillage = geom.get_constant_delta(DIRECTION_I);

  const double epsilon = epaisseur_lissage_ * h_maillage; // definition de la demi largeur de la zone d'etalement du terme de forcage comme un nombre de fois la largeur de maille selon la direction x

  DoubleTab coords(1,3); // tableau temporaire ou on stocke les coordonnees des points ou on veut interpoler la vitesse
  DoubleTab vitesses_interpolees(1,3);
  ArrOfDouble result(1); // tableau temporaire ou la fonction interpolation stocke le resultat
  // Boucle sur les tubes du faisceau
  for (int itube = 0; itube < ntubes; itube++)
    {
      const Tube_base& tube = faisceau[itube].valeur();
      const double tube_r = tube.get_rayon();
      const double r_tube_min = tube_r - epsilon;
      //const double r_tube_min_carre = r_tube_min * r_tube_min;

      const Vecteur3 position_cylindre = tube.get_current_pos();
      const Vecteur3 vitesse_cylindre = tube.get_current_velocity();

      // Boucle sur les trois composantes de qdm
      for (int direction = 0; direction < 3; direction++)   // boucle sur les 3 directions de la vitesse
        {

          const int ivoisin = (direction == 0) ? -1 : 0;
          const int jvoisin = (direction == 1) ? -1 : 0;
          const int kvoisin = (direction == 2) ? -1 : 0;
          const double volume_maille = geom.get_constant_delta(0) * geom.get_constant_delta(1) * geom.get_constant_delta(2); //donne le volume de la maille
          IJK_Field_double& velocity = select(direction, vx, vy, vz);  // velocity est egale a la vitesse selon la direction 0,1 ou 2
          int imin = 0;
          int jmin = 0;
          int kmin = 0;
          int imax = velocity.ni(); // donne le nombre de valeur de la vitesse selon i
          int jmax = velocity.nj();
          int kmax = velocity.nk();
          {
            // restreint le domaine ijk balaye a la zone couverte par le cylindre:
            const ArrOfDouble& x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I] : elem_pos[DIRECTION_I];
            const double xmin = position_cylindre[0] - (tube_r + epsilon); // on n'oublie pas que l'on doit chercher dans la zone couverte par le cylindre + etalement
            const double xmax = position_cylindre[0] + (tube_r + epsilon);
            while (imin < velocity.ni() && x_coord[imin + offset_i] < xmin)
              imin++;
            while (imax > 0 && x_coord[imax + offset_i - 1] > xmax)
              imax--;
          }
          {
            const ArrOfDouble& z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K] : elem_pos[DIRECTION_K];
            const double zmin = position_cylindre[2] - (tube_r + epsilon);
            const double zmax = position_cylindre[2] + (tube_r + epsilon);
            while (kmin < velocity.nk() && z_coord[kmin + offset_k] < zmin)
              kmin++;
            while (kmax > 0 && z_coord[kmax + offset_k - 1] > zmax)
              kmax--;
          }

          double delta_qdm_cylindre = 0.;
          double masse_fluide = 0.; // Masse de fluide dans le cylindre
          double volume = 0.;

          for (int k = kmin; k < kmax; k++)
            {
              // si direction==DIRECTION_I alors x_coord = coordonnee du noeud dans la direction I donc x_noeud et z_coord = z_element
              // si direction==DIRECTION_K alors x_coord = x_element et z_coord = z_noeud
              // x_coord et z_coord sont les abcisses des vitesse qui sont situee au centre des faces
              const double z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K][k+offset_k] : elem_pos[DIRECTION_K][k+offset_k];
              for (int j = jmin; j < jmax; j++)
                {
                  const double y_coord = (direction==DIRECTION_J) ? nodes_pos[DIRECTION_J][j+offset_j] : elem_pos[DIRECTION_J][j+offset_j];
                  for (int i = imin; i < imax; i++)
                    {
                      const double x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I][i+offset_i] : elem_pos[DIRECTION_I][i+offset_i];
                      const double delta_x = x_coord - position_cylindre[0]; // x_M - x_0
                      const double delta_z = z_coord - position_cylindre[2]; // z_M - z_0
                      const double d2 = delta_x * delta_x + delta_z * delta_z;
                      const double d = sqrt(d2); // distance du point M au centre du tube

                      // foraage si on est dans le cylindre
                      if (d < tube_r && d > r_tube_min)
                        {
                          // deuxieme essai: vitesse tangentielle egale a la vitesse tangentielle a la surface du cylindre
                          // vitesse normale nulle a l'interface, extrapolation a l'interieur du cylindre pour que
                          // div(u) = 0 (en fonction du gradient de vitesse tangentielle)

                          if (direction == DIRECTION_J)
                            {
                              // pas de forcage dans la direction j
                            }
                          else
                            {
                              //const double beta = tube_r / d;
                              //const double delta_x_P = beta * delta_x; // x_P - x_0
                              //const double delta_z_P = beta * delta_z; // z_P _ z_0
                              Vecteur3 normale;
                              normale[0] = delta_x / d;
                              normale[1] = 0.;
                              normale[2] = delta_z / d;
                              Vecteur3 tangente;
                              tangente[0] = - normale[2];
                              tangente[1] = 0.;
                              tangente[2] = normale[0];
#if 0
                              coords.resize(2,3);
                              coords(0,0) = position_cylindre[0] + delta_x_P - tangente_x * h_maillage;
                              coords(0,1) = y_coord;
                              coords(0,2) = position_cylindre[2] + delta_z_P - tangente_z * h_maillage;
                              coords(1,0) = position_cylindre[0] + delta_x_P + tangente_x * h_maillage;
                              coords(1,1) = y_coord;
                              coords(1,2) = position_cylindre[2] + delta_z_P + tangente_z * h_maillage;
                              Vecteur3 vitesse_interpolee;
                              Vecteur3 grad_v_theta;
                              ijk_interpolate(vx, coords, result); // la fonction remplit le tableau result
                              result[0] -=  vitesse_cylindre[0];
                              result[1] -=  vitesse_cylindre[0];
                              vitesse_interpolee[0]=(result[0]+result[1]) * 0.5;
                              grad_v_theta[0] = (result[1]-result[0])  / (2.*h_maillage);
                              /* ijk_interpolate(vy, coords, result); */
                              vitesse_interpolee[1]=0.;
                              grad_v_theta[1] = 0.;
                              ijk_interpolate(vz, coords, result);
                              result[0] -=  vitesse_cylindre[2];
                              result[1] -=  vitesse_cylindre[2];
                              vitesse_interpolee[2]=(result[0]+result[1]) * 0.5;
                              grad_v_theta[2] = (result[1]-result[0])  / (2.*h_maillage);
                              Vecteur3 vitesse_forcee;
                              // Gradient selon theta de la composante theta de vitesse:
                              const double grad_v_theta_theta =
                                grad_v_theta[0]*tangente_x + grad_v_theta[2]*tangente_z;
                              Vecteur3 normale;
                              normale[0] = delta_x / d;
                              normale[1] = 0.;
                              normale[2] = delta_z / d;
                              double v_scalaire_n =
                                vitesse_interpolee[0] *  normale[0]
                                /* + vitesse_interpolee[1] *  normale[1] */
                                + vitesse_interpolee[2] *  normale[2];
                              // Symetrique
                              for (int l=0; l<3; l++)
                                vitesse_forcee[l] = vitesse_cylindre[l] // referentiel cylindre->ref. fixe
                                                    + vitesse_interpolee[l] - normale[l] * v_scalaire_n // composante tangentielle
                                                    + 0.*normale[l] * grad_v_theta_theta * (tube_r/d-1); // composante normale
                              velocity(i, j, k) = vitesse_forcee[direction];
#else
                              coords.resize(3,3);
                              const double rr = tube_r * tube_r / d;
                              coords(0,0) = position_cylindre[0] + normale[0] * rr;
                              coords(0,1) = y_coord;
                              coords(0,2) = position_cylindre[2] + normale[2] * rr;

                              coords(1,0) = x_coord;
                              coords(1,1) = y_coord;
                              coords(1,2) = z_coord;

                              coords(2,0) = position_cylindre[0] + normale[0] * tube_r;
                              coords(2,1) = y_coord;
                              coords(2,2) = position_cylindre[2] + normale[2] * tube_r;

                              int l, m;
                              interpoler_vitesses_xz(coords, vitesses_interpolees, vx, vz, result);
                              for (l = 0; l < 3; l++)
                                for (m = 0; m < 3; m++)
                                  vitesses_interpolees(l,m) -= vitesse_cylindre[m];

                              // Composante normale = oppose de la composante normale au symetrique
                              // Composante tangentielle: on freine, un peu
                              double compo_normale = 0;
                              for (l = 0; l < 3; l++)
                                compo_normale += (- vitesses_interpolees(0,l)*rr/d - vitesses_interpolees(1,l)) * normale[l];
                              double compo_tan = 0;
                              for (l = 0; l < 3; l++)
                                compo_tan += (vitesses_interpolees(2,l) * (1 - 4*(tube_r - d) / tube_r) - vitesses_interpolees(1,l)) * tangente[l];
                              // Regularisation du terme de forcage:
                              double facteur = 1;
                              if (tube_r - d <  h_maillage)
                                facteur = (tube_r - d) / h_maillage;
                              if (d - r_tube_min < h_maillage)
                                facteur = (d - r_tube_min) / h_maillage;

                              velocity(i, j, k) += (compo_normale * normale[direction] + compo_tan * tangente[direction]) * facteur;

#endif
                            }
                        }
                      else
                        {
                          // pas de champ miroir forcage direct classique
                        }

                      // Calcul de rho sur la face: moyenne des rho sur les elements voisins
                      double rho_fluide = (rho_fluide_field(i,j,k) + rho_fluide_field(i+ivoisin, j+jvoisin, k+kvoisin)) * 0.5; // rho sur la face est la moy des deux cotes
                      //delta_qdm_cylindre += delta_v * volume_maille * rho_fluide;
                      masse_fluide += volume_maille * rho_fluide;
                      volume += volume_maille;
                    }
                }
            }

          integrale_force(itube, direction) = delta_qdm_cylindre;
          masse_fluide_cylindres(itube, direction) = masse_fluide;
          volume_cylindres(itube, direction) = volume;
        }
    }
  mp_sum_for_each_item(integrale_force);
  mp_sum_for_each_item(masse_fluide_cylindres);
  mp_sum_for_each_item(volume_cylindres);
}

void Couplage_Tubes_IBC::force_ibc_velocity_miroir(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                                   const IJK_Field_double& rho_fluide_field,
                                                   DoubleTab& masse_fluide_cylindres,
                                                   DoubleTab& volume_cylindres,
                                                   DoubleTab& integrale_force,
                                                   const Faisceau_Tubes& faisceau) const
{
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  ArrOfDouble nodes_pos[3]; //defini un tableau dynamique de 3 lignes pour la position des nodes
  ArrOfDouble elem_pos[3];  //defini un tableau dynamique de 3 lignes pour la position des elements
  {
    for (int dir = 0; dir < 3; dir++)
      {
        nodes_pos[dir] = geom.get_node_coordinates(dir); // node=noeud, les noeuds sont a l'intersection des mailles; on ne connait que la position des noeuds
        const int n = nodes_pos[dir].size_array()-1; // n= nombre de noeuds-1 = nombre d'elements dans la direction dir; les alements sont au centre des mailles
        elem_pos[dir].resize_array(n); // donne la taille du tableau contenant les elements
        for (int i = 0; i < n; i++)
          elem_pos[dir][i] = (nodes_pos[dir][i] + nodes_pos[dir][i+1]) * 0.5; // la position de l'element selon dir est moyenne de la positions des noeuds de part et d'autre
      }
  }
  const int offset_i = splitting.get_offset_local(DIRECTION_I); // donne l'offset selon la direction i
  const int offset_k = splitting.get_offset_local(DIRECTION_K);  // donne l'offset selon la direction k

  const int ntubes = faisceau.size();
  const double epsilon = epaisseur_lissage_ * geom.get_constant_delta(DIRECTION_I); // definition de la demi largeur de la zone d'etalement du terme de forcage comme un nombre de fois la largeur de maille selon la direction x

  // Boucle sur les tubes du faisceau
  for (int itube = 0; itube < ntubes; itube++)
    {
      const Tube_base& tube = faisceau[itube].valeur();
      const double tube_r = tube.get_rayon();
      const double omega = tube.get_omega();
      const double r_tube_min = tube_r - epsilon;
      const double r_tube_min_carre = r_tube_min * r_tube_min;

      const Vecteur3 position_cylindre = tube.get_current_pos();
      const Vecteur3 vitesse_cylindre = tube.get_current_velocity();

      // Boucle sur les trois composantes de qdm
      for (int direction = 0; direction < 3; direction++)   // boucle sur les 3 directions de la vitesse
        {

          const int ivoisin = (direction == 0) ? -1 : 0;
          const int jvoisin = (direction == 1) ? -1 : 0;
          const int kvoisin = (direction == 2) ? -1 : 0;

          const double volume_maille = geom.get_constant_delta(0) * geom.get_constant_delta(1) * geom.get_constant_delta(2); //donne le volume de la maille
          const double inv_delta_x = 1. / geom.get_constant_delta(0); // inverse de la taille de la maille selon x
          const double inv_delta_z = 1. / geom.get_constant_delta(2);
          IJK_Field_double& velocity = select(direction, vx, vy, vz);  // velocity est egale a la vitesse selon la direction 0,1 ou 2
          int imin = 0;
          int jmin = 0;
          int kmin = 0;
          int imax = velocity.ni(); // donne le nombre de valeur de la vitesse selon i
          int jmax = velocity.nj();
          int kmax = velocity.nk();
          {
            // restreint le domaine ijk balaye a la zone couverte par le cylindre:
            const ArrOfDouble& x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I] : elem_pos[DIRECTION_I];
            const double xmin = position_cylindre[0] - (tube_r + epsilon); // on n'oublie pas que l'on doit chercher dans la zone couverte par le cylindre + etalement
            const double xmax = position_cylindre[0] + (tube_r + epsilon);
            while (imin < velocity.ni() && x_coord[imin + offset_i] < xmin)
              imin++;
            while (imax > 0 && x_coord[imax + offset_i - 1] > xmax)
              imax--;
          }
          {
            const ArrOfDouble& z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K] : elem_pos[DIRECTION_K];
            const double zmin = position_cylindre[2] - (tube_r + epsilon);
            const double zmax = position_cylindre[2] + (tube_r + epsilon);
            while (kmin < velocity.nk() && z_coord[kmin + offset_k] < zmin)
              kmin++;
            while (kmax > 0 && z_coord[kmax + offset_k - 1] > zmax)
              kmax--;
          }

          double delta_qdm_cylindre = 0.;
          double masse_fluide = 0.; // Masse de fluide dans le cylindre
          double volume = 0.;

          for (int k = kmin; k < kmax; k++)
            {
              // si direction==DIRECTION_I alors x_coord = coordonnee du noeud dans la direction I donc x_noeud et z_coord = z_element
              // si direction==DIRECTION_K alors x_coord = x_element et z_coord = z_noeud
              // x_coord et z_coord sont les abcisses des vitesse qui sont situee au centre des faces
              const double z_coord = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K][k+offset_k] : elem_pos[DIRECTION_K][k+offset_k];
              for (int j = jmin; j < jmax; j++)
                {
                  //const double y_coord = (direction==DIRECTION_J) ? nodes_pos[DIRECTION_J][j+offset_j] : elem_pos[DIRECTION_J][j+offset_j];
                  for (int i = imin; i < imax; i++)
                    {
                      const double x_coord = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I][i+offset_i] : elem_pos[DIRECTION_I][i+offset_i];
                      const double delta_x = x_coord - position_cylindre[0]; // x_M - x_0
                      const double delta_z = z_coord - position_cylindre[2]; // z_M - z_0
                      const double d2 = delta_x * delta_x + delta_z * delta_z;
                      const double d = sqrt(d2); // distance du point M au centre du tube
                      const double rayon_tube_carre = tube_r * tube_r;

                      // foraage si on est dans le cylindre
                      if (d2 < rayon_tube_carre)
                        {
                          Vecteur3 rotation;
                          rotation[0] = omega * (- delta_z);
                          rotation[1] = 0.;
                          rotation[2] = omega * delta_x;
                          double extrapolated_v = vitesse_cylindre[direction] + rotation[direction];
                          double g_etalement;
                          double delta_v;
                          if (d2 > r_tube_min_carre )   // champ mirroir si d*d > (R-epsilon)*(R-epsilon)
                            {
                              g_etalement = 0.5 + 0.5 * (d - (tube_r - epsilon/2) ) / (epsilon * 0.5) - 0.5 / M_PI * sin( M_PI * (d - (tube_r - epsilon/2) )/ (epsilon*0.5) ); // fonction d'etalement
                              // on calcule les coordonnees de P qui est le sym de M par rapport a l'interface : OP=beta*OM en notation vectorielle
                              const double d_interface = tube_r - d; // distance du pt M a l'interface fluide solide
                              const double beta = (tube_r + d_interface) / (tube_r - d_interface);
                              const double delta_x_P = beta * delta_x; // x_P - x_0
                              const double delta_z_P = beta * delta_z; // z_P _ z_0

                              // determination des indices i_min, i_max, k_min et k_max pour connaitre le cube dans lequel est P afin d'interpoler la vitesse en P
                              // attention !!! : il s'agit d'indice locaux, ils permettent de se reperer que dans le proc considere
                              const int i_min = int(floor(i + (delta_x_P - delta_x) * inv_delta_x ));
                              const int i_max = i_min +1;
                              const int k_min = int(floor(k + (delta_z_P - delta_z) * inv_delta_z ));
                              const int k_max = k_min +1;
                              // a partir des indices locaux on determine les coordonnes globales des vitesses en ces points :
                              // x_11, z_11 correspond a i_min, k_min ; x_12, z_12 a i_max, k_min ; x_22, z_22 a i_max, k_max et x_21, z_21 a i_min, k_max.

                              //const double z_11 = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K][k_min+offset_k] : elem_pos[DIRECTION_K][k_min+offset_k];
                              //const double x_11 = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I][i_min+offset_i] : elem_pos[DIRECTION_I][i_min+offset_i];

                              const double z_12 = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K][k_max+offset_k] : elem_pos[DIRECTION_K][k_max+offset_k];
                              //const double x_12 = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I][i_min+offset_i] : elem_pos[DIRECTION_I][i_min+offset_i];

                              //const double z_22 = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K][k_max+offset_k] : elem_pos[DIRECTION_K][k_max+offset_k];
                              //const double x_22 = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I][i_max+offset_i] : elem_pos[DIRECTION_I][i_max+offset_i];

                              //const double z_21 = (direction==DIRECTION_K) ? nodes_pos[DIRECTION_K][k_min+offset_k] : elem_pos[DIRECTION_K][k_min+offset_k];
                              const double x_21 = (direction==DIRECTION_I) ? nodes_pos[DIRECTION_I][i_max+offset_i] : elem_pos[DIRECTION_I][i_max+offset_i];
                              // on peut alors calculer alpha et beta pour realiser l'interpolation bilineaire :

                              const double alpha = (z_12 - (  delta_z_P + position_cylindre[2] ) ) * inv_delta_z;
                              const double gamma  = (x_21 - (  delta_x_P + position_cylindre[0] ) ) * inv_delta_x;
                              // on fait l'interpolation bilineaire : mais comment on obtient la vitesse dans le repere global dc en i_min + offset_i par ex ????
                              // j'ai l'impression que velocity(i,j,k) c'est dans le repere local pas ds le global
                              // pour l'instant je le fais en monoproc il faudra dc mettre i_min + offset_i, j , k_min + offset_k
                              // Valeur des vitesses aux quatre points du carre ds lequel est P

                              // Test si on est bien a l'interieur de la zone locale sur ce processeur:
                              {
                                const int gh = velocity.ghost();
                                const int ni = velocity.ni();
                                const int nk = velocity.nk();
                                if (i_min < -gh || k_min < -gh || i_max >= ni + gh || k_max >= nk + gh)
                                  {
                                    Cerr << "Erreur, l'epaisseur ghost est insuffisante pour le calcul du champ mirroir ibc" << finl;
                                    Process::exit();
                                  }
                              }
                              const double v_11 = velocity(i_min, j, k_min);
                              const double v_21 = velocity(i_max, j, k_min);
                              const double v_22 = velocity(i_max, j, k_max);
                              const double v_12 = velocity(i_min, j, k_max);
                              // calcule de v_P par interpolation bilineaire
                              const double v_P = (1-alpha) * (1-gamma) * v_11 + alpha * gamma * v_22 + alpha * (1 - gamma) * v_12 + gamma * (1 - alpha) * v_21;

                              // On fait une symetrie axial de la vitesse de P en M par rapport a l'interface fluide solide : (v(M) - v(S)) = - (v(P) - v(S))
                              delta_v = (g_etalement * (-v_P + 2 * extrapolated_v) + (1-g_etalement) * extrapolated_v) - velocity(i, j, k); // qdm
                              velocity(i, j, k) = g_etalement * (-v_P + 2 * extrapolated_v) + (1-g_etalement) * extrapolated_v; // ne pas oublier de mettre la fonction g_etalement
# if 0
                              Cout << " direction : " << direction << " i : " << i << " j : "<< j << " k : " << k << " x_coord : " << x_coord << " z_coord : " << z_coord << endl;
                              Cout << " i_min : " << i_min << " i_max : "<< i_max << " k_min : " << k_min << " k_max : " << k_max << " x_P " << delta_x_P + position_cylindre[0] << " z_P " << delta_z_P + position_cylindre[2] << endl;
                              Cout << " alpha : " << alpha << " gamma : " << gamma << " x_11 " << x_11 << " z_11 " << z_11 << " x_21 " << x_21 ;
                              Cout << " z_21 " << z_21 << " x_22 " << x_22 << " z_22 " << z_22 <<  " x_12 " << x_12 << " z_12 " << z_12 <<endl;
                              Cout << " g_etalement : " << g_etalement << " v_P : " << v_P << " v(i,j,k) : " <<  velocity(i, j, k) << " v_11 " << v_11 ;
                              Cout <<  " v_21 " << v_21 <<  " v_22 " << v_22 <<  " v_12 " << v_12 << endl;
#endif

                            }
                          else     // pas de champ miroir forcage direct classique
                            {
                              delta_v = extrapolated_v - velocity(i,j,k);
                              velocity(i, j, k) = extrapolated_v ;
                              //Cout << " direction : " << direction << " i : " << i << " j : "<< j << " k : " << k << endl;
                              //Cout << " v(i,j,k) : " <<  velocity(i, j, k) << endl;
                            }

                          // Calcul de rho sur la face: moyenne des rho sur les elements voisins
                          double rho_fluide = (rho_fluide_field(i,j,k) + rho_fluide_field(i+ivoisin, j+jvoisin, k+kvoisin)) * 0.5; // rho sur la face est la moy des deux cotes
                          delta_qdm_cylindre += delta_v * volume_maille * rho_fluide;
                          masse_fluide += volume_maille * rho_fluide;
                          volume += volume_maille;
                        }
                    }
                }
            }

          integrale_force(itube, direction) = delta_qdm_cylindre;
          masse_fluide_cylindres(itube, direction) = masse_fluide;
          volume_cylindres(itube, direction) = volume;
        }
    }
  mp_sum_for_each_item(integrale_force);
  mp_sum_for_each_item(masse_fluide_cylindres);
  mp_sum_for_each_item(volume_cylindres);
}


//////////////////////////////
// Calcul des forces de pression par interpolation avec les valeurs de la pression a l interieur ou l exterieur du tube
//////////////////////////////
void Couplage_Tubes_IBC::calcul_F_pression(const IJK_Field_double& pressure, const IJK_Field_double& vx,
                                           DoubleTab& integrale_force_pression, DoubleTab& pression_teta,
                                           const Faisceau_Tubes& faisceau) const
{
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  ArrOfDouble nodes_pos[3]; //defini un tableau dynamique de 3 lignes pour la position des nodes
  ArrOfDouble elem_pos[3];  //defini un tableau dynamique de 3 lignes pour la position des elements
  {
    for (int dir = 0; dir < 3; dir++)
      {
        nodes_pos[dir] = geom.get_node_coordinates(dir); // node=noeud, les noeuds sont a l'intersection des mailles; on ne connait que la position des noeuds
        const int n = nodes_pos[dir].size_array()-1; // n= nombre de noeuds-1 = nombre d'elements dans la direction dir; les alements sont au centre des mailles
        elem_pos[dir].resize_array(n); // donne la taille du tableau contenant les elements
        for (int i = 0; i < n; i++)
          elem_pos[dir][i] = (nodes_pos[dir][i] + nodes_pos[dir][i+1]) * 0.5; // la position de l'element selon dir est moyenne de la positions des noeuds de part et d'autre
      }
  }
  int jmin = 0;
  int jmax = pressure.nj();
  //const int offset_i = splitting.get_offset_local(DIRECTION_I); // donne l'offset selon la direction i
  //const int offset_k = splitting.get_offset_local(DIRECTION_K);  // donne l'offset selon la direction k
  const int ntubes = faisceau.size();
  DoubleTab Pos_P;
  DoubleTab P_l_j;
  DoubleTab P_teta; // pression integre selon z pour le tube itube en l
  Pos_P.resize(ntubes , n_P_, 2); // tableau contenant pour chaque tube toutes les coordonnees x,z des point situes sur le cercle ou la pression est interpollees
  P_l_j.resize(ntubes , n_P_, jmax);
  P_teta.resize(ntubes, n_P_);
  const double teta = 2 * M_PI / n_P_;
  const double delta_x = geom.get_constant_delta(0);
  const double delta_z = geom.get_constant_delta(2);
  const double inv_delta_x = 1./delta_x;
  const double inv_delta_z = 1./delta_z;


  for (int itube = 0; itube < ntubes; itube++)
    {
      const Tube_base& tube = faisceau[itube].valeur();
      const double tube_r = tube.get_rayon();
      const Vecteur3 position_cylindre = tube.get_current_pos();
      double p_surface_cylindre_x = 0.;
      double p_surface_cylindre_z = 0.;
      const double dS = tube_r * teta * delta_z;

      for (int l =0; l < n_P_; l++)
        {
          P_teta(itube, l) = 0.;
          const double teta_l = l * teta;
          // coord = 0 correspond a x, coord = 1 correspond a z
          // x_P = x_cylindre + R * cos(teta_l); z_P = z_cylindre + R * sin(teta_l)
          Pos_P(itube, l, 0) = position_cylindre[0] + tube_r * cos(teta_l);
          Pos_P(itube, l, 1) = position_cylindre[2] + tube_r * sin(teta_l);
          const double x_P_Delta_x = Pos_P(itube , l, 0) * inv_delta_x;
          const double z_P_Delta_z = Pos_P(itube , l, 1) * inv_delta_z;
          const double resu_x = x_P_Delta_x - int(floor(x_P_Delta_x));
          const double resu_z = z_P_Delta_z - int(floor(z_P_Delta_z));
          // Attention ce sont les indices globaux !!!!!
          int i_min, i_max, k_min, k_max;
          if ( resu_x < 0.5)
            {
              i_min = int(floor(x_P_Delta_x))-1;
              i_max = i_min +1;
            }
          else
            {
              i_min = int(floor(x_P_Delta_x));
              i_max = i_min +1;
            }

          if ( resu_z < 0.5)
            {
              k_min = int(floor(z_P_Delta_z))-1;
              k_max = k_min +1;
            }
          else
            {
              k_min = int(floor(z_P_Delta_z));
              k_max = k_min +1;
            }
# if 0
          int i_min_local, i_max_local, k_min_local, k_max_local;
          i_min_local = i_min - 128;
          i_max_local = i_max - 128;
          k_min_local = k_min - 128;
          k_max_local = k_max - 128;


          // Test si on est bien a l'interieur de la zone locale sur ce processeur:
          {
            const int gh = pressure.ghost();
            const int ni = pressure.ni();
            const int nk = pressure.nk();
            Cout << " gh " << gh << " ni " <<  ni << " nk " << nk << endl;
            Cout << " l : " << l << " offset_i " << offset_i << " offset_k " << offset_k << endl;
            Cout << " i_min_l : " << i_min_local << " i_max_l : "<< i_max_local << " k_min_l : " << k_min_local << " k_max_l : " << k_max_local << " Pos_P(itube, l, 0) " <<  Pos_P(itube, l, 0) << " Pos_P(itube, l, 1) " << Pos_P(itube, l, 1) << endl;
            if (i_min_local < -gh || k_min_local < -gh || i_max_local >= ni + gh || k_max_local >= nk + gh)
              {
                Cerr << "Erreur, l'epaisseur ghost est insuffisante pour le calcul du champ mirroir ibc" << finl;
                Process::exit();
              }
          }

          //# endif
          const double z_11 = elem_pos[DIRECTION_K][k_min];
          const double x_11 = elem_pos[DIRECTION_I][i_min];


          const double x_12 = elem_pos[DIRECTION_I][i_min];

          const double z_22 = elem_pos[DIRECTION_K][k_max];
          const double x_22 = elem_pos[DIRECTION_I][i_max];
          const double z_21 = elem_pos[DIRECTION_K][k_min];
# endif
          const double x_21 = elem_pos[DIRECTION_I][i_max];
          const double z_12 = elem_pos[DIRECTION_K][k_max];


          const double alpha  = (z_12 - ( Pos_P(itube , l, 1) ) ) * inv_delta_z;
          const double gamma  = (x_21 - ( Pos_P(itube , l, 0) ) ) * inv_delta_x;
          // on fait l'interpolation bilineaire : mais comment on obtient la vitesse dans le repere global dc en i_min + offset_i par ex ????
          // j'ai l'impression que velocity(i,j,k) c'est dans le repere local pas ds le global
          // pour l'instant je le fais en monoproc il faudra dc mettre i_min + offset_i, j , k_min + offset_k
          // Valeur des vitesses aux quatre points du carre ds lequel est P
          for (int j = jmin; j < jmax; j++)
            {
              const double p_11 = pressure(i_min, j, k_min);
              const double p_21 = pressure(i_max, j, k_min);
              const double p_22 = pressure(i_max, j, k_max);
              const double p_12 = pressure(i_min, j, k_max);

              // calcule de la pression P par interpolation bilineaire pour le tube itube au point l,j :
              P_l_j(itube, l, j) = (1-alpha) * (1-gamma) * p_11 + alpha * gamma * p_22 + alpha * (1 - gamma) * p_12 + gamma * (1 - alpha) * p_21;
              p_surface_cylindre_x += P_l_j(itube, l, j) * cos(teta_l);
              p_surface_cylindre_z += P_l_j(itube, l, j) * sin(teta_l);
              P_teta(itube, l) += P_l_j(itube, l, j);
#if 0
              Cout << " l : " << l << " j : "<< j << endl;
              Cout << " i_min : " << i_min << " i_max : "<< i_max << " k_min : " << k_min << " k_max : " << k_max << " Pos_P(itube, l, 0) " <<  Pos_P(itube, l, 0) << " Pos_P(itube, l, 1) " << Pos_P(itube, l, 1) << endl;
              Cout << " alpha : " << alpha << " gamma : " << gamma << " x_11 " << x_11 << " z_11 " << z_11 << " x_21 " << x_21 ;
              Cout << " z_21 " << z_21 << " x_22 " << x_22 << " z_22 " << z_22 <<  " x_12 " << x_12 << " z_12 " << z_12 <<endl;
              Cout << " P_l_j(itube, l, j) : " << P_l_j(itube, l, j) << " p_11 " << p_11 <<  " p_21 " << p_21 <<  " p_22 " << p_22 <<  " p_12 " << p_12 << endl;
#endif

            }
          pression_teta(itube, l) = P_teta(itube, l) / jmax; // moy de la pression en l selon z
        }
      // Fx = - dS * [ Somme(j= a nj-1) Somme(l=0 a N-1) p(R,teta_l,j) * cos(teta_l) ]
      // Fz = - dS * [ Somme(j= a nj-1) Somme(l=0 a N-1) p(R,teta_l,j) * sin(teta_l) ]
      integrale_force_pression(itube, 0) = - dS * p_surface_cylindre_x;
      integrale_force_pression(itube, 1) = - dS * p_surface_cylindre_z;

    }


  mp_sum_for_each_item(integrale_force_pression);
  mp_sum_for_each_item(pression_teta);


}


//////////////////////////////
// Calcul des forces de pression en ne prenant les valeurs de la pression que si elles sont a l exterieur du tube
//////////////////////////////
void Couplage_Tubes_IBC::calcul_F_pression2(const IJK_Field_double& pressure, const IJK_Field_double& vx,
                                            DoubleTab& integrale_force_pression, DoubleTab& pression_teta,
                                            const Faisceau_Tubes& faisceau) const
{
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  ArrOfDouble nodes_pos[3]; //defini un tableau dynamique de 3 lignes pour la position des nodes
  ArrOfDouble elem_pos[3];  //defini un tableau dynamique de 3 lignes pour la position des elements
  {
    for (int dir = 0; dir < 3; dir++)
      {
        nodes_pos[dir] = geom.get_node_coordinates(dir); // node=noeud, les noeuds sont a l'intersection des mailles; on ne connait que la position des noeuds
        const int n = nodes_pos[dir].size_array()-1; // n= nombre de noeuds-1 = nombre d'elements dans la direction dir; les alements sont au centre des mailles
        elem_pos[dir].resize_array(n); // donne la taille du tableau contenant les elements
        for (int i = 0; i < n; i++)
          elem_pos[dir][i] = (nodes_pos[dir][i] + nodes_pos[dir][i+1]) * 0.5; // la position de l'element selon dir est moyenne de la positions des noeuds de part et d'autre
      }
  }
  int jmin = 0;
  int jmax = pressure.nj();
  const int offset_i = splitting.get_offset_local(DIRECTION_I); // donne l'offset selon la direction i
  const int offset_k = splitting.get_offset_local(DIRECTION_K);  // donne l'offset selon la direction k
  const int ntubes = faisceau.size();
  DoubleTab Pos_P;
  DoubleTab P_l_j;
  DoubleTab P_teta; // pression integre selon z pour le tube itube en l
  Pos_P.resize(ntubes , n_P_, 2); // tableau contenant pour chaque tube toutes les coordonnees x,z des point situes sur le cercle ou la pression est interpollees
  P_l_j.resize(ntubes , n_P_, jmax);
  P_teta.resize(ntubes, n_P_);
  const double teta = 2 * M_PI / n_P_;
  const double delta_x = geom.get_constant_delta(0);
  const double delta_y = geom.get_constant_delta(1);
  const double delta_z = geom.get_constant_delta(2);
  const double inv_delta_x = 1./delta_x;
  const double inv_delta_z = 1./delta_z;


  for (int itube = 0; itube < ntubes; itube++)
    {
      const Tube_base& tube = faisceau[itube].valeur();
      const double tube_r = tube.get_rayon();
      const Vecteur3 position_cylindre = tube.get_current_pos();
      double p_surface_cylindre_x = 0.;
      double p_surface_cylindre_z = 0.;
      const double dS = tube_r * teta * delta_y;
      // Cout << "teta : " << teta << " delta_z : " << delta_z << " delta_x : " << delta_x << " delta_y : " << delta_y << finl;

      for (int l =0; l < n_P_; l++)
        {
          P_teta(itube, l) = 0.;
          const double teta_l = l * teta;
          // coord = 0 correspond a x, coord = 1 correspond a z
          // x_P = x_cylindre + R * cos(teta_l); z_P = z_cylindre + R * sin(teta_l)
          Pos_P(itube, l, 0) = position_cylindre[0] + tube_r * cos(teta_l);
          Pos_P(itube, l, 1) = position_cylindre[2] + tube_r * sin(teta_l);
          const double x_P_Delta_x = Pos_P(itube , l, 0) * inv_delta_x;
          const double z_P_Delta_z = Pos_P(itube , l, 1) * inv_delta_z;
          // maille i,k contenant le point P
          // Attention ce sont les indices globaux !!!!!
          const int i_M = int(floor(x_P_Delta_x));
          const int k_M = int(floor(z_P_Delta_z));
          int i_p = 0;
          int k_p = 0;
          double ratio =0;
          ratio =0; // je fais ca ici car j'ai un bug de compilation "ratio defined but not used"
          ratio += ratio; // je rajoute ca car ca bug toujours  sur fluent
          // Retirer ration cree un bug de compilation :
          // Couplage_Tubes_IBC.cpp:2854:5: error: 'ratio' was not declared in this scope
          // Couplage_Tubes_IBC.cpp:2862:3: error: 'ratio' was not declared in this scope
          // commenter la part de code en #IF 0 ne change rien
          // a voir du cote des optimisations du compilateur
          // determination des coordonnes de l'element M situe au centre de la maille i,k

          const double x_M = elem_pos[DIRECTION_I][i_M];
          const double z_M = elem_pos[DIRECTION_K][k_M];
          const double d_OM2 = (x_M - position_cylindre[0]) * (x_M - position_cylindre[0]) + (z_M - position_cylindre[2]) * (z_M - position_cylindre[2]);
          const double R2 = tube_r * tube_r ;
          /*
          #if 0
            Cout << " offset_i : " << offset_i << " offset_k : " << offset_k << endl;
            Cout << " l : " << l << " x_l :" << Pos_P(itube, l, 0) << " z_L : " << Pos_P(itube, l, 1) << " x_M : " << x_M;
            Cout << " z_M : " << z_M << " i_M : "<< i_M << " k_M : "<< k_M << " d_OM2 : "<< d_OM2 << " R2 : " << R2 << endl;
            Cerr << " offset_i : " << offset_i << " offset_k : " << offset_k << endl;
            Cerr << " l : " << l << " x_l :" << Pos_P(itube, l, 0) << " z_L : " << Pos_P(itube, l, 1) << " x_M : " << x_M << " z_M : " << z_M ;
            Cerr << " i_M : "<< i_M << " k_M : "<< k_M << " d_OM2 : "<< d_OM2 << " R2 : " << R2 << endl;
            #endif
          */
          if ( d_OM2 < R2)
            {
              const double f = (tube_r + sqrt( delta_x * delta_x + delta_z * delta_z ))/tube_r; // on va chercher la valeur de la pression dans une maille situee en dehros du tube a une distance d= f*tube_r
              // calcul des coordonnees du point E, qui est situe sur la droite OM a la distance OM+racine(deltax*deltax + deltaz*deltaz).
              const double x_E =  position_cylindre[0] + f * ( Pos_P(itube, l, 0) -  position_cylindre[0] );
              const double z_E = position_cylindre[2] + f * ( Pos_P(itube, l, 1) -  position_cylindre[2] );
              const int i_E = int(floor(x_E * inv_delta_x));
              const int k_E = int(floor(z_E * inv_delta_z));
              // indice globaux de lelement ou on prend p
              i_p = i_E;
              k_p = k_E;
              ratio = 1. /f ;

              //Cout << " x_E : " << x_E << " z_E : " << z_E << " i_p : " << i_p << " k_p : " << k_p << " ratio : " << ratio << endl;
              //Cerr << " x_E : " << x_E << " z_E : " << z_E << " i_p : " << i_p << " k_p : " << k_p << " ratio : " << ratio << endl;
            }
          else
            {
              i_p = i_M;
              k_p = k_M;
              ratio = tube_r / sqrt(d_OM2);
              //Cout << " i_p : " << i_p << " k_p : " << k_p << " ratio : " << ratio << " sqrt_d_OM2 : "  << sqrt(d_OM2) << endl;
              //Cerr << " i_p : " << i_p << " k_p : " << k_p << " ratio : " << ratio << " sqrt_d_OM2 : "  << sqrt(d_OM2) << endl;
            }
          // Ajout BM:
          //Cout <<  " offset_i + pressure.ni() : " << offset_i + pressure.ni() << endl;
          //Cerr <<  " offset_i + pressure.ni() : " << offset_i + pressure.ni() << endl;
          if (i_p >= offset_i && i_p < offset_i + pressure.ni() && k_p >= offset_k && k_p < offset_k + pressure.nk())
            {

              for (int j = jmin; j < jmax; j++)
                {
                  // P_l_j(itube, l, j) = ratio * pressure_(i_p, j, k_p);
                  P_l_j(itube, l, j) = pressure(i_p-offset_i, j, k_p-offset_k);
                  //Cerr << " i_p-offset_i : " << i_p-offset_i << " i_p-offset_k : "<< k_p-offset_k << endl;
                  p_surface_cylindre_x += P_l_j(itube, l, j) * cos(teta_l);
                  p_surface_cylindre_z += P_l_j(itube, l, j) * sin(teta_l);
                  P_teta(itube, l) += P_l_j(itube, l, j);
#if 0
                  Cout << " l : " << l << " j : "<< j << " P_l_j(itube, l, j) : " << P_l_j(itube, l, j) << " pressure(i_M, j, k_M) : " ;
                  Cout << " p_surface_cylindre_x : " <<  p_surface_cylindre_x << " p_surface_cylindre_z : " << p_surface_cylindre_z << endl;
                  Cerr << " l : " << l << " j : "<< j << " P_l_j(itube, l, j) : " << P_l_j(itube, l, j) << " pressure(i_M, j, k_M) : " ;
                  Cerr<< " p_surface_cylindre_x : " <<  p_surface_cylindre_x << " p_surface_cylindre_z : " << p_surface_cylindre_z << endl;

#endif

                }
            }
          pression_teta(itube, l) = P_teta(itube, l) / jmax; // moy de la pression en l selon z
        }
      // Fx = - dS * [ Somme(j= a nj-1) Somme(l=0 a N-1) p(R,teta_l,j) * cos(teta_l) ]
      // Fz = - dS * [ Somme(j= a nj-1) Somme(l=0 a N-1) p(R,teta_l,j) * sin(teta_l) ]
      //Cout << "p_surface_cylindre_x : " << p_surface_cylindre_x << " p_surface_cylindre_z : " << p_surface_cylindre_z << " dS : " << dS << finl;
      integrale_force_pression(itube, 0) = - dS * p_surface_cylindre_x;
      integrale_force_pression(itube, 1) = - dS * p_surface_cylindre_z;
      //Cout << "integrale_force_pression(itube, 0) : " << integrale_force_pression(itube, 0) << " integrale_force_pression(itube, 1) : " << integrale_force_pression(itube, 0) << finl;

    }


  mp_sum_for_each_item(integrale_force_pression);
  mp_sum_for_each_item(pression_teta);


}



///////////////////////////////////////////////////
//
//////////////////////////////////////////////////

void Couplage_Tubes_IBC::update(double current_time,
                                double timestep)
{
  // a modifier plus tard pour les tubes libres:
  //Vecteur3 force_appliquee(0,0,0);
  Vecteur3 force_appliquee;
  Vecteur3 pression;
  Vecteur3 volume_cylindre;
  Vecteur3 masse_fluide_cylindres;
  Vecteur3 F_IBC;
  Vecteur3 F_IBC_post_proj;
  Vecteur3 m_a;

  // faisceau_[i] est un objet de type DERIV(Tube_base)
  // pour acceder a l'objet Tube_impose ou Tube_libre ou autre, il faut utiliser
  // faisceau_[i].valeur() qui est ici de type Tube_base.
  const int ntubes = faisceau_.size();
  for (int i = 0; i < ntubes; i++)
    {
      // calcul de la force appliquee
      const Tube_base& tube = faisceau_[i].valeur();  // faisceau_[itube].valeur() est de type Tube_base
      Vecteur3 terme_inertiel = tube.get_current_acceleration(); // pour l'instant terme_inertiel = acceleration
      Vecteur3 F_ibc_extrapole;
      for (int dir = 0; dir < 3; dir++)
        {
          terme_inertiel[dir] *=masse_fluide_cylindres_(i,dir); // maintenant : terme_inertiel[dir] = acceleration * m_F
          // F_ibc_extrapole = F_IBC(n) + [t(n+1) - t(n)] * [F_IBC(n) - F_IBC(n-1)] / [t(n) - t(n-1)] //Calcul d'ordre 1 de la force IBC au temps n+1
          F_ibc_extrapole[dir] = integrale_force_(i,dir) + timestep * d_integrale_force_(i,dir);
          // F fluide sur solide pour le newmark = - (F_IBC_extrapole - m_F*acceleration_cylindre)
          force_appliquee[dir]= -( F_ibc_extrapole[dir] - terme_inertiel[dir]);
          volume_cylindre[dir] = volume_cylindres_(i,dir);
          masse_fluide_cylindres[dir] = masse_fluide_cylindres_(i,dir);
          F_IBC[dir] = integrale_force_(i,dir);
          F_IBC_post_proj[dir] = integrale_force_post_proj_(i,dir);
          m_a[dir]   = terme_inertiel[dir];
        }

      pression[0] = integrale_force_pression_(i,0);
      pression[1] = integrale_force_pression_(i,1);
      pression[2] = 0.;
      // Affichage des caracteristiques des tubes
      const Vecteur3 acceleration = tube.get_current_acceleration();
      const Vecteur3 pos = tube.get_current_pos();
      const Vecteur3 v = tube.get_current_velocity();

      // Calcul de la hauteur du tube
      const IJK_Splitting& splitting = ref_splitting_.valeur();
      const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
      ArrOfDouble noeuds_y;
      noeuds_y = geom.get_node_coordinates(1); // rentre les coordonnees des noeuds selon y dans le tableau noeuds_y
      const int n_element_y = noeuds_y.size_array()-1; // nombre d'elements selon y = nombre de noeuds -1
      const double delta_y = geom.get_constant_delta(1);
      const double hauteur_cylindre = n_element_y * delta_y;
      //
      const double tube_r = tube.get_rayon();
      // Calcul de Cx, Cy et Cz : C=F/(0.5rho_F*D*U*U*h)
      Vecteur3 C = force_appliquee;
      const double Adim = 1./(0.5 * rho_fluide_pour_adim_ * vitesse_pour_adim_ * vitesse_pour_adim_ * 2 * tube_r * hauteur_cylindre);
      C *=Adim;
      Vecteur3 C_IBC = F_IBC;
      C_IBC *=Adim;

      // Affichage des carac des tubes
      Cout << "Caracteristiques du tube " << tube.nom_du_tube() << " de type " << tube.que_suis_je() << finl; // tube est de type nom_, que_suis_je() est une fonction TrioU qui permet de lire ce type
      Cout << "Pos : " << pos << " t= " << current_time << finl;
      Cout << "V : " << v << " t= " << current_time << finl;
      Cout << "A : " << acceleration << " t= " << current_time << finl;
      Cout << "C : " << C << " t= " << current_time << finl;
      Cout << "F : " << force_appliquee << " t= " << current_time << finl;
      Cout << "F_IBC : " << F_IBC << " t= " << current_time << finl;
      Cout << "F_IBC_post_proj : " << F_IBC_post_proj << " t= " << current_time << finl;
      Cout << "C_IBC : " << C_IBC << " t= " << current_time << finl;
      Cout << "F_P_non_adim : " << pression << " t= " << current_time << finl;
      pression *=Adim;
      Cout << "F_P : " << pression << " t= " << current_time << finl;
      Cout << "V_c : " << volume_cylindre << " t= " << current_time << finl;
      Cout << "M_F : " << masse_fluide_cylindres << " t= " << current_time << finl;
      Cout << "m_a : " << m_a << " t= " << current_time << finl;
      Cout << " " << finl;
      // Actualisation de la position du tube i
      faisceau_[i].valeur().update_vitesse_position(current_time, timestep, force_appliquee);

    }
}

#if 0
void Couplage_Tubes_IBC::update(int tstep,
                                double timestep,
                                const double vitesse_entree)
{
  if (tstep <= t_lache_ - 1)
    return;

  switch(mouvement_)
    {
    case MOUVEMENT_LIBRE:
      newmark_implicit_update(Vz_,
                              az_,
                              donnees_tubes_,
                              tstep,
                              position_initiale_cylindres_(Num_cylindre_mobile_,1),
                              timestep,
                              vitesse_entree);
      break;
    case MOUVEMENT_IMPOSE:
      {
        double pos = oscillation_cylindre(az_,
                                          Vz_,
                                          tstep,
                                          position_initiale_cylindres_(Num_cylindre_mobile_,1),
                                          timestep,
                                          masse_fluide_cylindres_(Num_cylindre_mobile_,2),
                                          vitesse_entree);
        donnees_tubes_(Num_cylindre_mobile_, 1) = pos;
        break;
      }
    default:
      Cerr << "Erreur dans Couplage_Tubes_IBC::update: mouvement " << mouvement_ << " non code" << finl;
      exit();
    }
}

// renvoie la position du cylindre
double Couplage_Tubes_IBC::oscillation_cylindre(double& az,
                                                double& vz,
                                                int tstep,
                                                double position_initiale_cylindre_oscillant,
                                                const double timestep,
                                                const double masse_fluide_dans_cylindre,
                                                const double vitesse_entree) const
{
  const int Num_cylindre = Num_cylindre_mobile_;
  double position_cylindre_oscillant;
  int t = tstep;
  double retard = t_lache_;
  double z_cylindre;
  const double A = A_;
  const double pulsation = pulsation_;
  const double pulsation_2 = pulsation_*pulsation_;
  const DoubleTab integrale_force = integrale_force_;
  double Fz;
  double Fie;
  Fie = - masse_fluide_dans_cylindre * az;
  z_cylindre = A*sin(pulsation * (t-retard)*timestep);
  double Az;
  Az = z_cylindre /  (2 * tube_r_);
  vz =  pulsation *A* cos(pulsation * (t-retard)*timestep);
  az = - pulsation_2 * z_cylindre;//- pulsation2 *A* sin(pulsation * (t-retard)*timestep);
  Fz = - (integrale_force(Num_cylindre,2) + Fie);
  position_cylindre_oscillant = position_initiale_cylindre_oscillant + z_cylindre;
  double Coeff_adim, Cx, Cz;
  Coeff_adim = 0.5 * rho_fluide_pour_adim_ * vitesse_entree * vitesse_entree * tube_r_ * 2 * h_cylindre_ ;
  Cx = - integrale_force(Num_cylindre,0) / Coeff_adim;
  Cz = Fz / Coeff_adim;
  Cout << "z_cylindre " << position_cylindre_oscillant << " v_cylindre " << vz << " a_cylindre " << az << " z_absolue " << z_cylindre << " Az " << Az<<finl;
  Cout << "F_inertie " <<- Fie << " F_IBC " << - integrale_force(Num_cylindre,2)  << " F_f_s " << Fz  << finl;
  Cout << "Cx " << Cx << " Cz " << Cz << finl;
  return position_cylindre_oscillant;
}

void Couplage_Tubes_IBC::newmark_implicit_update(double& vz,
                                                 double& az,
                                                 DoubleTab& donnees_tubes_,
                                                 int tstep,
                                                 double position_initiale_cylindre_oscillant,
                                                 const double timestep,
                                                 const double vitesse_entree) const
{
  const double beta = 0.25;
  const double gamma = 0.5;
  const int Num_cylindre = Num_cylindre_mobile_;
  const DoubleTab integrale_force = integrale_force_;
  const DoubleTab integrale_force_N_moins_1 = integrale_force_N_moins_1_;
  int t = tstep;
  const double A = A_;
  const double B= B_;
  const double  pulsation = pulsation_;
  const double pulsation2 = pulsation2_;
  double F_imposee;
  double retard = t_lache_;
  F_imposee = A*sin(pulsation * (t-retard)*timestep) + B*sin(pulsation2 * (t-retard)*timestep);
  Cout << "F_imposee  " << F_imposee  << finl;
  double z_N;
  ArrOfDouble a(8);
  a[0] = 1/( beta * timestep * timestep);
  a[1] = gamma /( beta * timestep);
  a[2] =  1 /( beta * timestep);
  a[3] = 1/(2 *beta)-1;
  a[4] = gamma/ beta -1;
  a[5] = 0.5* timestep * (gamma/ beta -2);
  a[6] = timestep * (1- gamma);
  a[7] = gamma * timestep;

  double m_s = rho_cylindre_* volume_cylindres_(0,2);
  double m_F = rho_fluide_pour_adim_ * volume_cylindres_(0,2);
  Cout << "Volume du cylindre " << volume_cylindres_(0,2)  << finl;
  double m_T = m_s - m_F;
  double m_R = m_s / m_F;
  Cout << "m_s " <<  m_s  << " m_F " <<  m_F  << " m_T " <<  m_T  << " m_R " <<  m_R  << finl;
  double c= c_;
  double k= k_;
  DoubleTab donnees_tubes = donnees_tubes_;
  Cout << "Coefficients newmark: " << a;
  z_N = (donnees_tubes(Num_cylindre,1) - position_initiale_cylindre_oscillant); // z du cylindre au temps n
  double Az_n;
  Az_n = z_N / (2 * tube_r_);
  Cout << "z_cylindre_n " << z_N  << " v_cylindre_n " << vz <<" a_cylindre_n " << az_ <<" Az_n " << Az_n <<finl;
  // calcul de la force du fluide sur le solide au temps n+1
  double F_IBC_extrapol, F_IBC_precedent, F_IBC_n;
  F_IBC_precedent = integrale_force_N_moins_1_(Num_cylindre,2);
  F_IBC_n = integrale_force(Num_cylindre,2);
  F_IBC_extrapol = 2 * integrale_force(Num_cylindre,2) - integrale_force_N_moins_1_(Num_cylindre,2); //Calcul d'ordre 1 de la force IBC au temps n+1
  Cout << "F_IBC_precedent " << F_IBC_precedent << " F_IBC_n " << F_IBC_n << " F_IBC_extrapol " << F_IBC_extrapol << finl;
  // F_IBC_extrapol =  integrale_force(Num_cylindre,2) // Calcul d'ordre 0
  double Fie, F_f_s;
  Fie = - rho_fluide_pour_adim_ * volume_cylindres_(0,2)*az; // au temps n
  F_f_s = - (F_IBC_n + Fie); // au temps n
  Cout << "F_inertie " << Fie <<" F_IBC_n_bis " << F_IBC_n << " F_f_s " << F_f_s  << finl; // au temps n
  double Coeff_adim, Cx, Cz;
  Coeff_adim = 0.5 * rho_fluide_pour_adim_ * vitesse_entree * vitesse_entree * tube_r_ * 2 * h_cylindre_ ;
  Cx = - integrale_force(Num_cylindre,0) / Coeff_adim;
  Cz = F_f_s / Coeff_adim;
  Cout << "Cx " << Cx <<" Cz " << Cz  << finl; // au temps n
  //
  double K_tilde, F_tilde, K_droite;
  K_tilde = a[0] * m_T + a[1] * c + k;
  K_droite = m_T * (a[0] * z_N + a[2] * vz + a[3] * az) + c * (a[1] * z_N + a[4] * vz + a[5] * az);
  const double para_libre = para_libre_;
  F_tilde =  K_droite - F_IBC_extrapol * para_libre  + F_imposee;
  Cout <<"K_tilde " << K_tilde << " K_droite " << K_droite  << " Ftilde " << F_tilde  << finl;
  // calcul de z(n+1), vz(n+1) et az(n+1)
  double z_N_plus_1;
  z_N_plus_1 = F_tilde / K_tilde ;
  // Calcul de l'energie macanique au temps n
  double Ec_n, Ep_n, Em_n, Ef;
  Ec_n = 0.5 * m_T * vz * vz;
  Ep_n = 0.5 * k * z_N*z_N;
  Em_n = Ec_n + Ep_n;
  Ef = -F_IBC_extrapol * vz * timestep;
  Cout << "Ec_n " <<  Ec_n  << " Ep_n " <<  Ep_n << " Em_n " <<  Em_n << " Ef " <<  Ef  << finl;
  double az_N;
  az_N = az; // on stocke l'accalaration au temps n ds az_N
  az = a[0] * (z_N_plus_1 - z_N) - a[2] *vz - a[3] *az; // az (n+1) = a_0 * (z_N_plus_1 - z_N) - a_2 vz(n) - a_3 az(n)
  vz = vz + a[6] * az_N + a[7] * az; // vz(n+1) = vz(n) + a_6 * az(n) + a_7 * az (n+1)
  // actualisation de la position du cylindre ds le tableau contenant les coordonnees
  donnees_tubes_(Num_cylindre,1) = z_N_plus_1 + position_initiale_cylindre_oscillant;
  Cout << "z_cylindre_n_1 " <<  donnees_tubes_(Num_cylindre,1) << " z_absolue " << z_N_plus_1  << " v_cylindre_n_1 " << vz <<" a_cylindre_n_1 " << az_ << finl;
  double Ec_n_1, Ep_n_1, Em_n_1, DeltaEm, eta;
  Ec_n_1 = 0.5 * m_T * vz * vz;
  Ep_n_1 = 0.5 * k * z_N_plus_1 * z_N_plus_1;
  Em_n_1 = Ec_n_1 + Ep_n_1;
  DeltaEm = Em_n_1 - Em_n;
  eta = (DeltaEm/Ef-1) *100;
  Cout << "Em " << Em_n  << " Em_(n+1) " << Em_n_1 << " DeltaEm " << DeltaEm << " Eta " <<eta << finl;
}
#endif
#if 0
void print()
{
  Cout << "Integrale des forces sur tous les cylindres\n" << finl;
  for (int i = 0; i < integrale_force_.dimension(0); i++)
    {
      Cout << "F_cylindre " << i << " position " << i%n_tubes_x << " " << i/n_tubes_z << " ";
      for (int j = 0; j < integrale_force_.dimension(1); j++)
        {
          Cout << -integrale_force_(i,j) << " ";
        }
      Cout << finl;
    }
  Cout << "\n" << finl;
  Cout << "Volume cylindres\n" << finl;
  for (int i = 0; i < integrale_force_.dimension(0); i++)
    {
      Cout << "V_tube " << i << " position " << i%n_tubes_x << " " << i/n_tubes_z << " ";
      for (int j = 0; j < integrale_force_.dimension(1); j++)
        {
          Cout << Volume_cylindre_(i,j) << " ";
        }
      Cout << finl;
    }

}
#endif


#if 0
// Premier essai: symetrique de la vitesse tangentielle, vitesse normale nulle
double g_etalement;
if (d2 > r_tube_min_carre )   // champ mirroir si d*d > (R-epsilon)*(R-epsilon)
  {
    if (d > tube_r)
      {
        g_etalement = 1.-(d - tube_r)/(h_maillage*2);

      }
    else if (d > tube_r - h_maillage)
      g_etalement = 1.;
    else
      g_etalement = (tube_r - h_maillage - d) / (tube_r - h_maillage - r_tube_min);

    // on calcule les coordonnees de P qui est le sym de M par rapport a l'interface : OP=beta*OM en notation vectorielle
    double beta;
    if (d > tube_r)
      beta = 1.;
    else
      beta = (2.*tube_r-d) / d;
    const double delta_x_P = beta * delta_x; // x_P - x_0
    const double delta_z_P = beta * delta_z; // z_P _ z_0

    coords(0,0) = position_cylindre[0] + delta_x_P;
    coords(0,1) = y_coord;
    coords(0,2) = position_cylindre[2] + delta_z_P;
    Vecteur3 vitesse_interpolee;
    ijk_interpolate(vx, coords, result); // la fonction remplit le tableau result
    vitesse_interpolee[0]=result[0];
    /* ijk_interpolate(vy, coords, result); */
    vitesse_interpolee[1]=0.;
    ijk_interpolate(vz, coords, result);
    vitesse_interpolee[2]=result[0];

    // Calcul du symetrique de la vitesse par rapport a la direction normale:
    Vecteur3 vitesse_forcee;
    // Retranche la vitesse du cylindre pour se placer dans le referentiel du cylindre
    for (int l=0; l<3; l++)
      vitesse_interpolee[l] -= vitesse_cylindre[l];
    // Produit scalaire vitesse.normale_au_cylindre  (la composante y de la normale est nulle
    Vecteur3 normale;
    normale[0] = delta_x / d;
    normale[1] = 0.;
    normale[2] = delta_z / d;
    double v_scalaire_n =
      vitesse_interpolee[0] *  normale[0]
      /* + vitesse_interpolee[1] *  normale[1] */
      + vitesse_interpolee[2] *  normale[2];
    // Symetrique

    for (int l=0; l<3; l++)
      vitesse_forcee[l] = vitesse_cylindre[l] + vitesse_interpolee[l] - normale[l] * v_scalaire_n;

    velocity(i, j, k) = g_etalement * vitesse_forcee[direction];
#endif
