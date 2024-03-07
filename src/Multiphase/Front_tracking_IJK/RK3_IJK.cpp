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
// File      : RK3_IJK.cpp
// Directory : $IJK_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////
#include <RK3_IJK.h>
#include <Param.h>
#include <Interprete_bloc.h>
#include <SFichier.h>
#include <stat_counters.h>
#include <IJK_Lata_writer.h>
#include <IJK_Navier_Stokes_tools.h>
#include <EFichier.h>
#include <communications.h>
#include <LecFicDiffuse_JDD.h>


Implemente_instanciable(IJK_problem_double, "IJK_problem_double", Interprete);
Entree& IJK_problem_double::interpreter(Entree& is)
{
  check_divergence_ = 0;
  Param param(que_suis_je());
  Nom ijk_splitting_name;
  compteur_post_instantanes_ = 0;
  dt_post_ = 100;
  dt_sauvegarde_ = 2000000000; // jamais
  current_time_ = 0.;
  nom_sauvegarde_ = nom_du_cas() + ".sauv";
  gravite_.resize_array(3);
  gravite_ = 0.;
  ghost_size_pour_ibc_ = 2; // par defaut, pareil que pour les operateurs convection etc...

  fichier_reprise_rho_ = "??";
  fichier_reprise_alpha_ = "??";
  fichier_reprise_vitesse_ = "??";
  timestep_reprise_rho_ = -2;
  timestep_reprise_alpha_ = -2;
  timestep_reprise_vitesse_ = -2;
  expression_vitesse_initiale_.dimensionner(3);


  timestep_facsec_ = -1; // disabled by default
  param.ajouter("ijk_splitting", &ijk_splitting_name, Param::REQUIRED);
  param.ajouter("timestep", &timestep_, Param::REQUIRED);
  param.ajouter("timestep_facsec", &timestep_facsec_);
  param.ajouter("nb_pas_dt_max", &nb_timesteps_, Param::REQUIRED);
  param.ajouter("pressure_gradient", &source_pressure_gradient_, Param::REQUIRED);
  param.ajouter("multigrid_solver", &poisson_solver_, Param::REQUIRED);
  param.ajouter_flag("check_divergence", &check_divergence_);
  param.ajouter("mu_liquide", &mu_liquide_, Param::REQUIRED);
  param.ajouter("mu_gaz", &mu_gaz_, Param::REQUIRED);
  param.ajouter("vitesse_entree", &vitesse_entree_, Param::REQUIRED);
  param.ajouter("rho_liquide", &rho_liquide_, Param::REQUIRED);
  param.ajouter("rho_gaz", &rho_gaz_, Param::REQUIRED);
  param.ajouter("couplage_tubes_ibc", &couplage_tubes_ibc_, Param::REQUIRED);
  param.ajouter("dt_post", &dt_post_);
  param.ajouter("champs_a_postraiter", &liste_post_instantanes_);
  param.ajouter("check_stop_file", &check_stop_file_);
  param.ajouter("dt_sauvegarde", &dt_sauvegarde_);
  param.ajouter("nom_sauvegarde", &nom_sauvegarde_);
  param.ajouter("nom_reprise", &nom_reprise_);
  param.ajouter("gravite", &gravite_);

  param.ajouter("expression_vx_init", &expression_vitesse_initiale_[0]);
  param.ajouter("expression_vy_init", &expression_vitesse_initiale_[1]);
  param.ajouter("expression_vz_init", &expression_vitesse_initiale_[2]);
  param.ajouter("expression_rho_init", &expression_rho_initiale_);
  param.ajouter("expression_alpha_init", &expression_alpha_initiale_);
  // ATTENTION les fichiers reprises sont des fichiers .lata ou sauv.lata
  // On peut reprendre uniquement la vitesse ou uniquement rho dans un fichier de post:
  param.ajouter("fichier_reprise_vitesse", &fichier_reprise_vitesse_);
  param.ajouter("fichier_reprise_rho", &fichier_reprise_rho_);
  param.ajouter("fichier_reprise_alpha", &fichier_reprise_alpha_);
  param.ajouter("timestep_reprise_vitesse", &timestep_reprise_vitesse_);
  param.ajouter("timestep_reprise_rho", &timestep_reprise_rho_);
  param.ajouter("timestep_reprise_alpha", &timestep_reprise_alpha_);

  param.ajouter("boundary_conditions", &boundary_conditions_);
  param.ajouter("ghost_size_pour_ibc", &ghost_size_pour_ibc_);

  param.ajouter("champ_miroir_pour_RHS", &champ_miroir_pour_RHS_, Param::REQUIRED);
  param.dictionnaire("Active", ACTIVE); // Attention le dictionnaire s'applique sur le dernier param ajoute !!!!
  param.dictionnaire("Desactive", DESACTIVE);
  param.ajouter("choix_schema_transport_rho", &choix_schema_transport_rho_, Param::REQUIRED);
  param.dictionnaire("despres_samuel", DESPRES_SAMUEL); // Attention le dictionnaire s'applique sur le dernier param ajoute !!!!
  param.dictionnaire("despres_dellacherie", DESPRES_DELLACHERIE);
  param.dictionnaire("amont", AMONT);


  param.lire_avec_accolades(is);
  if (gravite_.size_array() != 3)
    {
      Cerr << "Erreur: la gravite doit etre un vecteur de 3 composantes" << finl;
      exit();
    }
  splitting_ = ref_cast(IJK_Splitting, Interprete_bloc::objet_global(ijk_splitting_name));

  if (nom_reprise_ != "??")
    reprendre_probleme(nom_reprise_);

  run();
  return is;
}

Sortie& IJK_problem_double::printOn(Sortie& os) const
{
  return os;
}
Entree& IJK_problem_double::readOn(Entree& is)
{
  return is;
}

void IJK_problem_double::posttraiter_champs_instantanes(const char *lata_name, double current_time)
{
  const int latastep = compteur_post_instantanes_;
  dumplata_newtime(lata_name,current_time);
  if (liste_post_instantanes_.contient_("TOUS"))
    {
      liste_post_instantanes_.dimensionner_force(0);
      liste_post_instantanes_.add("VELOCITY");
      liste_post_instantanes_.add("VELOCITY_MIROIR");
      liste_post_instantanes_.add("PRESSURE");
      liste_post_instantanes_.add("RHO");
      liste_post_instantanes_.add("RHO_BIS");
      liste_post_instantanes_.add("ALPHA");
      liste_post_instantanes_.add("MU");
      liste_post_instantanes_.add("PRESSURE_RHS");
      liste_post_instantanes_.add("DV_DT");
      liste_post_instantanes_.add("DRHO_DT");
      //liste_post_instantanes_.add("LEVELSET");
      //liste_post_instantanes_.add("DLEVELSET_DT");

    }
  int n = liste_post_instantanes_.size();
  if (liste_post_instantanes_.contient_("VELOCITY"))
    n--,dumplata_vector(lata_name,"VELOCITY", velocity_[0], velocity_[1], velocity_[2], latastep);
  if (liste_post_instantanes_.contient_("VELOCITY_MIROIR"))
    n--,dumplata_vector(lata_name,"VELOCITY_MIROIR", velocity_tmp_[0], velocity_tmp_[1], velocity_tmp_[2], latastep);
  if (liste_post_instantanes_.contient_("PRESSURE"))
    n--,dumplata_scalar(lata_name,"PRESSURE", pressure_, latastep);
  if (liste_post_instantanes_.contient_("RHO"))
    n--,dumplata_scalar(lata_name,"RHO", rho_field_, latastep);
  if (liste_post_instantanes_.contient_("RHO_BIS"))
    n--,dumplata_scalar(lata_name,"RHO_BIS", rho_field_bis_, latastep);
  if (liste_post_instantanes_.contient_("ALPHA"))
    n--,dumplata_scalar(lata_name,"ALPHA", alpha_field_, latastep);
  if (liste_post_instantanes_.contient_("MU"))
    n--,dumplata_scalar(lata_name,"MU", molecular_mu_, latastep);
  if (liste_post_instantanes_.contient_("PRESSURE_RHS"))
    n--,dumplata_scalar(lata_name,"PRESSURE_RHS", pressure_rhs_, latastep);
  if (liste_post_instantanes_.contient_("DV_DT"))
    n--,dumplata_vector(lata_name,"DV_DT", d_velocity_[0], d_velocity_[1], d_velocity_[2], latastep);
  if (liste_post_instantanes_.contient_("DRHO_DT"))
    n--,dumplata_scalar(lata_name,"DRHO_DT", d_rho_field_, latastep);
  //if (liste_post_instantanes_.contient_("LEVELSET"))
  //  n--,dumplata_scalar(lata_name,"LEVELSET", fonction_levelset_, latastep);
  //if (liste_post_instantanes_.contient_("DLEVELSET_DT"))
  //  n--,dumplata_scalar(lata_name,"DLEVELSET_DT", d_fonction_levelset_, latastep);


  if (n>0)
    {
      Cerr << "Il y a des noms de champs a postraiter inconnus dans la liste de champs a postraiter" << finl;
      Process::exit();
    }
  compteur_post_instantanes_++;
}

void set_data(IJK_Field_double& f, double func(double, double, double))
{
  const IJK_Grid_Geometry& geom = f.get_splitting().get_grid_geometry();
  const int i_min = f.get_splitting().get_offset_local(DIRECTION_I);
  const int j_min = f.get_splitting().get_offset_local(DIRECTION_J);
  const int k_min = f.get_splitting().get_offset_local(DIRECTION_K);
  const int ni_tot = geom.get_nb_elem_tot(DIRECTION_I);
  const int nj_tot = geom.get_nb_elem_tot(DIRECTION_J);
  const int nk_tot = geom.get_nb_elem_tot(DIRECTION_K);
  for (int k = 0; k < f.nk(); k++)
    for (int j = 0; j < f.nj(); j++)
      for (int i = 0; i < f.ni(); i++)
        f(i,j,k)=func((double)(i+i_min)/ni_tot,
                      (double)(j+j_min)/nj_tot,
                      (double)(k+k_min)/nk_tot);
}

double funci_double(double x, double y, double z)
{
  //return 6000+sin(x*2*M_PI*13)*sin(y*2*M_PI*9)*sin(z*2*M_PI*7)*100;
  return 0;
}
double funcj_double(double x, double y, double z)
{
  //return sin(x*2*M_PI*6)*sin(y*2*M_PI*5)*sin(z*2*M_PI*11)*100;
  return 0.;
}
double funck_double(double x, double y, double z)
{
  //return sin(x*2*M_PI*8)*sin(y*2*M_PI*3)*sin(z*2*M_PI*12)*100;
  return 0.;
}
#define select(a,x,y,z) ((a==0)?(x):((a==1)?(y):(z)))


double function_rho_init(double x, double y, double z)
{
  if ( x<0.1 && sin(x*M_PI) * sin(z*2*M_PI*2) > 0)
    return 1000;
  else
    return 500;
}

double function_levelset_init(double x, double y, double z)
{
  //double levelset = sin(x*2*M_PI*2) * sin(z*2*M_PI*4);
  double levelset;
  if (x<0.2)
    levelset = sin(z*2*M_PI*8);
  else levelset =-1;
  return levelset;
}

//void calcul_integrale_pression_sur_ibc(0
void force_entry_velocity(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz, double v_imposed)
{
  const IJK_Splitting& splitting = vx.get_splitting();
  const int offset_i = splitting.get_offset_local(DIRECTION_I);
  if (offset_i > 0)
    return;
  {
    double imposed[3] = {0., 0., 0.};
    imposed[0] = v_imposed;
    for (int direction = 0; direction < 3; direction++)
      {
        IJK_Field_double& velocity = select(direction, vx, vy, vz);
        const int imin = 0;
        const int jmin = 0;
        const int kmin = 0;
        const int imax = 3;
        const int jmax = velocity.nj();
        const int kmax = velocity.nk();
        for (int k = kmin; k < kmax; k++)
          {
            for (int j = jmin; j < jmax; j++)
              {
                for (int i = imin; i < imax; i++)
                  {
                    velocity(i,j,k) = imposed[direction];
                  }
              }
          }
      }
  }
}

void IJK_problem_double::sauvegarder_probleme(const char *fichier_sauvegarde)
{
  Nom lata_name(fichier_sauvegarde);
  lata_name += ".lata";
  SFichier fichier;
  if (Process::je_suis_maitre())
    {
      fichier.ouvrir(fichier_sauvegarde);
      Cerr << "T= " << current_time_ << " Checkpointing dans le fichier " << fichier_sauvegarde << finl;
      fichier.precision(17);
      fichier.setf(std::ios_base::scientific);
      fichier << "{\n"
              << " tinit " << current_time_ << "\n"
              << " fichier_reprise_vitesse " << lata_name << "\n"
              << " fichier_reprise_rho " << lata_name << "\n";
      fichier << " fichier_reprise_alpha " << lata_name << "\n"
              << " timestep_reprise_vitesse 1\n"
              << " timestep_reprise_rho 1\n"
              << " timestep_reprise_alpha 1\n"
              << "}\n";
    }
  dumplata_header(lata_name, velocity_[0] /* on passe un champ pour ecrire la geometrie */);
  dumplata_newtime(lata_name,current_time_);
  dumplata_vector(lata_name,"VELOCITY", velocity_[0], velocity_[1], velocity_[2], 0);
  dumplata_scalar(lata_name,"RHO", rho_field_, 0);
  dumplata_scalar(lata_name,"ALPHA", alpha_field_, 0);
  //dumplata_scalar(lata_name,"LEVELSET", fonction_levelset_, 0);

  // Sauvegarde des donnees cylindres
  couplage_tubes_ibc_.sauvegarder_probleme(fichier);
}

void IJK_problem_double::sauvegarder_pression()
{
  Nom nom_fichier = "pression_teta";
  //nom_fichier +=current_time_;
  SFichier fichier (nom_fichier);
  //fichier << "T= " << current_time_ << finl;
  couplage_tubes_ibc_.sauvegarder_pression(fichier);
}

void IJK_problem_double::reprendre_probleme(const char *fichier_reprise)
{
  // Lecture par tous les processeurs, on retire les commentaires etc...
  LecFicDiffuse_JDD fichier(fichier_reprise);
  Param param(que_suis_je());
  param.ajouter("tinit", &current_time_);
  param.ajouter("fichier_reprise_vitesse", &fichier_reprise_vitesse_);
  param.ajouter("fichier_reprise_rho", &fichier_reprise_rho_);
  param.ajouter("fichier_reprise_alpha", &fichier_reprise_alpha_);
  param.ajouter("timestep_reprise_vitesse", &timestep_reprise_vitesse_);
  param.ajouter("timestep_reprise_rho", &timestep_reprise_rho_);
  param.ajouter("timestep_reprise_alpha", &timestep_reprise_alpha_);
  param.lire_avec_accolades(fichier);
  fichier_reprise_vitesse_=dirname(fichier_reprise)+fichier_reprise_vitesse_;
  couplage_tubes_ibc_.reprendre_probleme(fichier);
  // Appeler ensuite initialize() pour lire les fichiers lata etc...
}


void force_boundary_condition_v(IJK_Field_double& vx, IJK_Field_double& vy,  IJK_Field_double& vz)
{
  const int nj = vz.nj();
  const int ni = vz.ni();
  const int kmin = vz.get_splitting().get_offset_local(DIRECTION_K);
  const int nktot = vz.get_splitting().get_nb_items_global(IJK_Splitting::FACES_K, DIRECTION_K);
  if (kmin == 0)
    {
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          vz(i,j,0) = 0.;
    }
  if (kmin + vz.nk() == nktot)
    {
      const int k = vz.nk()-1;
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          vz(i,j,k) = 0.;
    }
}

double find_timestep(const IJK_Field_double& v, int direction)
{
  double max_v = 0.;
  const int ni = v.ni();
  const int nj = v.nj();
  const int nk = v.nk();
  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              max_v = std::max(max_v, fabs(v(i,j,k)));
            }
        }
    }
  max_v = Process::mp_max(max_v);
  const double delta = v.get_splitting().get_grid_geometry().get_constant_delta(direction);
  return delta / max_v;
}

// Methode appelee par run() une fois la memoire alouee pour les champs.
// Cette fonction remplit les valeurs initiales de vitesse et rho
void IJK_problem_double::initialise()
{
  Cout << que_suis_je() << "::initialise()" << finl;
  if ( fichier_reprise_rho_ == "??" )   // si on ne fait pas une reprise on initialise T et rho suivant une methode
    {
      Cout << "Initialisation rho, expression = " << expression_rho_initiale_ << finl;
      set_field_data(rho_field_, expression_rho_initiale_);
    }
  else
    {
      Cout << "Lecture rho initial dans fichier " << fichier_reprise_rho_ << " timestep= " << timestep_reprise_rho_ << finl;
      lire_dans_lata(fichier_reprise_rho_, timestep_reprise_rho_, "DOM", "RHO", rho_field_);
    }

  if ( fichier_reprise_alpha_ == "??" )   // si on ne fait pas une reprise on initialise T et alpha suivant une methode
    {
      Cout << "Initialisation alpha, expression = " << expression_alpha_initiale_ << finl;
      set_field_data(alpha_field_, expression_alpha_initiale_);
    }
  else
    {
      Cout << "Lecture alpha initial dans fichier " << fichier_reprise_alpha_ << " timestep= " << timestep_reprise_alpha_ << finl;
      lire_dans_lata(fichier_reprise_alpha_, timestep_reprise_alpha_, "DOM", "ALPHA", alpha_field_);
    }

  if (fichier_reprise_vitesse_ == "??")   // si on ne fait pas une reprise on initialise V
    {
      if (expression_vitesse_initiale_.size() != 3)
        {
          Cerr << "Erreur dans l'initialisation: la vitesse initiale doit etre fournie avec trois expressions" << finl;
          Process::exit();
        }
      Cout << "Initialisation vitesse \nvx = " << expression_vitesse_initiale_[0]
           << "\nvy = " <<  expression_vitesse_initiale_[1]
           << "\nvz = " <<  expression_vitesse_initiale_[2]  << finl;
      for (int i = 0; i < 3; i++)
        set_field_data(velocity_[i], expression_vitesse_initiale_[i]);
    }
  else
    {
      Cout << "Lecture vitesse initiale dans fichier " << fichier_reprise_vitesse_ << " timestep= " << timestep_reprise_vitesse_ << finl;
      lire_dans_lata(fichier_reprise_vitesse_, timestep_reprise_vitesse_, "DOM", "VELOCITY",
                     velocity_[0], velocity_[1], velocity_[2]); // fonction qui lit un champ a partir d'un lata .
    }

}

void IJK_problem_double::run()
{
  splitting_.get_local_mesh_delta(DIRECTION_K, 2 /* ghost cells */, delta_z_local_);
  Cerr << "IJK_problem_double::run()" << finl;
  allocate_velocity(velocity_, splitting_, ghost_size_pour_ibc_);
  allocate_velocity(velocity_tmp_, splitting_, ghost_size_pour_ibc_);
  allocate_velocity(d_velocity_, splitting_, 1);
  pressure_.allocate(splitting_, IJK_Splitting::ELEM, ghost_size_pour_ibc_);
  molecular_mu_.allocate(splitting_, IJK_Splitting::ELEM, 2);
  molecular_mu_.data() = 1.;
  Cerr << " Allocating 11 arrays, approx total size= "
       << molecular_mu_.data().size_array() * (int)sizeof(double) * 11 * 9.537E-07 << " MB per core" << finl;
  pressure_rhs_.allocate(splitting_, IJK_Splitting::ELEM, 1);
  rho_field_.allocate(splitting_, IJK_Splitting::ELEM, 2);

  rho_field_bis_.allocate(splitting_, IJK_Splitting::ELEM, 2);

  d_rho_field_.allocate(splitting_, IJK_Splitting::ELEM, 0);
  alpha_field_.allocate(splitting_, IJK_Splitting::ELEM, 2);
  //fonction_levelset_.allocate(splitting_, IJK_Splitting::ELEM, 2);
  //d_fonction_levelset_.allocate(splitting_, IJK_Splitting::ELEM, 0);

  velocity_diffusion_op_.typer("OpDiffIJK_double");
  velocity_diffusion_op_.set_bc(boundary_conditions_);
  velocity_diffusion_op_.initialize(splitting_);

  velocity_convection_op_.typer("OpConvIJKQuickSharp_double");
  velocity_convection_op_.initialize(splitting_);

  velocity_convection_op_.typer("OpConvIJKQuickSharp_double");
  scalar_convection_op_.initialize(splitting_);

  poisson_solver_.initialize(splitting_);

  rho_field_.data() = 1.;

  rho_field_bis_.data() = 1.;

  alpha_field_.data() = 1.;
  initialise();

  rho_field_.echange_espace_virtuel(2);
  rho_field_bis_.echange_espace_virtuel(2);
  alpha_field_.echange_espace_virtuel(2);
  recalculer_mu_de_rho(rho_field_, molecular_mu_, 2);

  pressure_projection(velocity_[0], velocity_[1],  velocity_[2], pressure_, 1.,
                      pressure_rhs_, check_divergence_, poisson_solver_);
  const double max_timestep = timestep_;

  couplage_tubes_ibc_.initialize(splitting_);

  // Preparer le fichier de postraitement et postraiter la condition initiale:
  const Nom lata_name = nom_du_cas() + Nom(".lata");
  dumplata_header(lata_name, velocity_[0] /* on passe un champ pour ecrire la geometrie */);
  posttraiter_champs_instantanes(lata_name, current_time_);

  int stop = 0;
  for (int tstep = 0; tstep < nb_timesteps_ && stop == 0; tstep++)
    {
      force_boundary_condition_v(velocity_[0], velocity_[1], velocity_[2]);


      //force_boundary_condition_rho(rho_field_, current_time_);

      if (timestep_facsec_ > 0.)
        {
          timestep_ = max_timestep;
          for (int dir = 0; dir < 3; dir++)
            timestep_ = std::min(timestep_, find_timestep(velocity_[dir], dir) * timestep_facsec_);
          Cout << "Automatic timestep: " << timestep_ << finl;
        }
      couplage_tubes_ibc_.update(current_time_,
                                 timestep_);

      euler_time_step();

      current_time_+= timestep_;
      // verification du fichier stop
      stop = 0;
      if (check_stop_file_ != "??")
        {
          if (je_suis_maitre())
            {
              EFichier f;
              stop = f.ouvrir(check_stop_file_);
              if (stop)
                {
                  // file exists, check if it contains 1:
                  f >> stop;
                }
            }
          envoyer_broadcast(stop, 0);
        }
      velocity_[1].data() = 0.;
      if (tstep == nb_timesteps_ - 1)
        stop = 1;
      if (tstep % dt_sauvegarde_ == dt_sauvegarde_-1 || stop)
        {
          sauvegarder_probleme(nom_sauvegarde_);
          sauvegarder_pression();
        }

      if (tstep % dt_post_ == dt_post_-1 || stop)
        {
          posttraiter_champs_instantanes(lata_name, current_time_);
        }
    }
}

void IJK_problem_double::recalculer_rho_de_levelset(const IJK_Field_double& levelset,
                                                    IJK_Field_double& rho,
                                                    int epaisseur_joint) const
{
  const int imax = rho.ni() + epaisseur_joint;
  const int jmax = rho.nj() + epaisseur_joint;
  const int kmax = rho.nk() + epaisseur_joint;

  for (int k = -epaisseur_joint; k < kmax ; k++)
    for (int j = -epaisseur_joint; j < jmax; j++)
      for (int i = -epaisseur_joint; i < imax; i++)
        {
          const double l = levelset(i, j, k);
          double r;
          if (l > 0)
            r = rho_liquide_;
          else
            r = rho_gaz_;
          rho(i,j,k) = r;
        }
}

void IJK_problem_double::recalculer_mu_de_rho(const IJK_Field_double& rho,
                                              IJK_Field_double& molecular_mu,
                                              int epaisseur_joint) const
{
  const int imax = rho.ni() + epaisseur_joint;
  const int jmax = rho.nj() + epaisseur_joint;
  const int kmax = rho.nk() + epaisseur_joint;

  if (rho_liquide_ == rho_gaz_)
    {
      if (mu_liquide_ != mu_gaz_)
        {
          Cerr << " Erreur dans recalculer_mu_de_rho: rho_l = rho_g mais mu_l != mu_g" << finl;
          exit();
        }
      molecular_mu.data() = mu_liquide_;
    }
  else
    {
      const double constante = (- mu_liquide_ * rho_gaz_ + mu_gaz_ * rho_liquide_) / (rho_liquide_ - rho_gaz_);
      const double facteur = (mu_liquide_ - mu_gaz_) / (rho_liquide_ - rho_gaz_);

      for (int k = -epaisseur_joint; k < kmax ; k++)
        for (int j = -epaisseur_joint; j < jmax; j++)
          for (int i = -epaisseur_joint; i < imax; i++)
            {
              const double r = rho(i, j, k);
              assert(r > 0.);
              molecular_mu(i,j,k) = r * facteur + constante;
            }
    }
}

void IJK_problem_double::euler_time_step()
{
  static Stat_Counter_Id projection_ = statistiques().new_counter(0, "projection");
  velocity_[0].echange_espace_virtuel(2);
  velocity_[1].echange_espace_virtuel(2);
  velocity_[2].echange_espace_virtuel(2);
  rho_field_.echange_espace_virtuel(2);
  rho_field_bis_.echange_espace_virtuel(2);
  alpha_field_.echange_espace_virtuel(2);

  IJK_Field_double rho_field_nplus1 = rho_field_;
  IJK_Field_double alpha_field_nplus1 = alpha_field_;

  // Convection de rho avec un schema Despres Lagoutiere
  // ATTENTION il faut faut utiliser le Despres apres avoir calculer la pression et le champ de vitesse, car le DESPRES ecrase rho_field_ au temps n par rho_field_ au temps n+1
  //scalar_convection_op_Despres_Lagoutiere_optimise_samuel(rho_field_nplus1, velocity_[0], velocity_[1], velocity_[2], d_rho_field_);


  switch(choix_schema_transport_rho_)
    {
    case DESPRES_SAMUEL:
      // Despres type samuel pour le transport de rho
      scalar_convection_op_Despres_Samuel(rho_field_nplus1, velocity_[0], velocity_[1], velocity_[2]);
      break;
    case DESPRES_DELLACHERIE:
      // Despres type Dellacherie pour le transport de rho
      scalar_convection_op_Despres_Dellacherie(rho_field_nplus1, velocity_[0], velocity_[1], velocity_[2]);
      break;
    case AMONT:
      // schema de transport decentre amont pour rho
      scalar_convection_op_Amont_rho(rho_field_nplus1, velocity_[0], velocity_[1], velocity_[2]);
      break;
    default:
      Cerr << "Erreur dans RK3_IJK::euler_time_step " << choix_schema_transport_rho_ << " non code" << finl;
      exit();
    }

  calcul_integrale_rho(rho_field_nplus1);
  //fonction_levelset_.echange_espace_virtuel(2);



  //  nouveau despres
  //scalar_convection_op_DESPRES_alpha_AMONT_rho(alpha_field_nplus1, rho_field_nplus1, velocity_[0], velocity_[1], velocity_[2]);
  //scalar_convection_op_DESPRES_alpha_AMONT_rho_TYPE_DELACHERIE(alpha_field_nplus1, rho_field_nplus1, velocity_[0], velocity_[1], velocity_[2]);

  //calculer_rho2(alpha_field_nplus1, rho_field_bis_);

  // Convection / diffusion qdm:
  // on a besoin d'une epaisseur de joint sur mu pour la diffusion
  //recalculer_rho_de_levelset(fonction_levelset_, rho_field_, 2);
  recalculer_mu_de_rho(rho_field_, molecular_mu_, 2);

# if 0
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Affichage vitesse initiale
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  IJK_Field_double& v_initiale = velocity_[0];
  Cout << " " << finl;
  Cout << "v_initiale : " << v_initiale(64, 0, 73) << " " << v_initiale(64, 0, 74) << " " << v_initiale(64, 0, 75) << " " << v_initiale(64, 0, 76) << finl;
  Cout << " " << finl;
#endif
  switch(champ_miroir_pour_RHS_)
    {
    case ACTIVE:
      // Calcul du champ miroir
      velocity_tmp_[0] = velocity_[0];
      velocity_tmp_[1] = velocity_[1];
      velocity_tmp_[2] = velocity_[2];
      velocity_tmp_[0].echange_espace_virtuel(ghost_size_pour_ibc_);
      velocity_tmp_[1].echange_espace_virtuel(ghost_size_pour_ibc_);
      velocity_tmp_[2].echange_espace_virtuel(ghost_size_pour_ibc_);
      couplage_tubes_ibc_.champ_miroir(velocity_tmp_[0], velocity_tmp_[1], velocity_tmp_[2], rho_field_, timestep_, pressure_);



      // Calcul d_velocity = v scalaire gradient(v)
      //(les compteurs statistiques sont deja presents a l'interieur des fonctions ajouter/calculer)
      //statistiques().begin_count(convection_counter_);
      velocity_convection_op_.calculer(velocity_tmp_[0], velocity_tmp_[1], velocity_tmp_[2],
                                       velocity_tmp_[0], velocity_tmp_[1], velocity_tmp_[2],
                                       d_velocity_[0], d_velocity_[1], d_velocity_[2]);
      //statistiques().end_count(convection_counter_);
      // Multiplication par rho aux faces (on va rediviser a la fin)
      calculer_rho_v(rho_field_, d_velocity_, d_velocity_);

      // Calcul d_velocity
      //(les compteurs statistiques sont deja presents a l'interieur des fonctions ajouter/calculer)
      //statistiques().begin_count(diffusion_counter_);
      velocity_diffusion_op_.set_nu(molecular_mu_);
      velocity_diffusion_op_.ajouter(velocity_tmp_[0], velocity_tmp_[1], velocity_tmp_[2],
                                     d_velocity_[0], d_velocity_[1], d_velocity_[2]);
      //statistiques().end_count(diffusion_counter_);
      statistiques().begin_count(source_counter_);


      break;
    case DESACTIVE:

      //couplage_tubes_ibc_.forcage_anticipe(velocity_[0], velocity_[1], velocity_[2], rho_field_, timestep_, pressure_);

      // Calcul d_velocity = v scalaire gradient(v)
      //(les compteurs statistiques sont deja presents a l'interieur des fonctions ajouter/calculer)
      //statistiques().begin_count(convection_counter_);
      velocity_convection_op_.calculer(velocity_[0], velocity_[1], velocity_[2],
                                       velocity_[0], velocity_[1], velocity_[2],
                                       d_velocity_[0], d_velocity_[1], d_velocity_[2]);
      //statistiques().end_count(convection_counter_);
      // Multiplication par rho aux faces (on va rediviser a la fin)
      calculer_rho_v(rho_field_, d_velocity_, d_velocity_);

      // Calcul d_velocity
      //(les compteurs statistiques sont deja presents a l'interieur des fonctions ajouter/calculer)
      //statistiques().begin_count(diffusion_counter_);
      velocity_diffusion_op_.set_nu(molecular_mu_);
      velocity_diffusion_op_.ajouter(velocity_[0], velocity_[1], velocity_[2],
                                     d_velocity_[0], d_velocity_[1], d_velocity_[2]);
      //statistiques().end_count(diffusion_counter_);
      statistiques().begin_count(source_counter_);

      break;
    default:
      Cerr << "Erreur dans RK3_IJK::euler_time_step " << champ_miroir_pour_RHS_ << " non code" << finl;
      exit();
    }


  // Convection rho:
  //scalar_convection_op_.calculer(rho_field_, velocity_[0], velocity_[1], velocity_[2], d_rho_field_);
  // Convection de rho avec un schema amont tire du code Triton
  //scalar_convection_op_Triton(rho_field_, velocity_[0], velocity_[1], velocity_[2], d_rho_field_);

  // Convection levelset:
  //scalar_convection_op_.calculer(fonction_levelset_, velocity_[0], velocity_[1], velocity_[2], d_fonction_levelset_);

  // diffusion creates non zero velocity at walls, restore zero

  // conditions a la paroi
  force_boundary_condition_v(d_velocity_[0], d_velocity_[1], d_velocity_[2]);
  for (int dir = 0; dir < 3; dir++)
    {
      const int kmax = d_velocity_[dir].nk();
      for (int k = 0; k < kmax; k++)
        {
          mass_solver_with_rho(d_velocity_[dir], rho_field_, delta_z_local_, k);

          terme_source_gravite(d_velocity_[dir], k, dir);
          /////////////////////////
          /// Calcul du champ predit
          /////////////////////////
          euler_explicit_update(d_velocity_[dir], velocity_[dir], k);
        }
    }
  //////////////////////////
  /// Forcage
  ///////////////////////////
  velocity_[0].echange_espace_virtuel(ghost_size_pour_ibc_);
  velocity_[1].echange_espace_virtuel(ghost_size_pour_ibc_);
  velocity_[2].echange_espace_virtuel(ghost_size_pour_ibc_);


# if 0
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Affichage vitesse predite
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  IJK_Field_double& v_predite = velocity_[0];
  Cout << " " << finl;
  Cout << "v_p : " << v_predite(64, 0, 54) << " " << v_predite(64, 0, 55) << " " << v_predite(64, 0, 76) << " " << v_predite(64, 0, 77) << finl;
  Cout << " " << finl;
#endif
  couplage_tubes_ibc_.force_ibc(velocity_[0], velocity_[1], velocity_[2], rho_field_, timestep_, pressure_, current_time_);
# if 0
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Affichage vitesse impose
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  IJK_Field_double& v_impose = velocity_[0];
  Cout << " " << finl;
  Cout << "v_s : " << v_impose(64, 0, 54) << " " << v_impose(64, 0, 55) << " " << v_impose(64, 0, 76) << " " << v_impose(64, 0, 77) << finl;
  Cout << " " << finl;
#endif
  // Conditions en entree
  force_entry_velocity(velocity_[0], velocity_[1], velocity_[2], vitesse_entree_);

  //////////////////////////
  /// Etape de projection
  /////////////////////////
  statistiques().end_count(source_counter_);
  //(le compteur statistique est deja present dans la fonction de projection)
  //statistiques().begin_count(projection_);
  pressure_projection_with_rho(rho_field_, velocity_[0], velocity_[1],  velocity_[2], pressure_, timestep_,
                               pressure_rhs_, check_divergence_, poisson_solver_);

  //statistiques().end_count(projection_);

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Calcul des forces post projection
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //  couplage_tubes_ibc_.calcul_force_post_projection(velocity_[0], velocity_[1], velocity_[2], rho_field_, timestep_, current_time_);
# if 0
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Affichage vitesse corrigee
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  IJK_Field_double& v_corrige = velocity_[0];
  Cout << " " << finl;
  Cout << "v_corrige : " << v_corrige(64, 0, 54) << " " << v_corrige(64, 0, 55) << " " << v_corrige(64, 0, 76) << " " << v_corrige(64, 0, 77) << finl;
  Cout << " " << finl;
# endif
# if 0
  // forcage apres l'etape de projection

  velocity_tmp_[0] = velocity_[0];
  velocity_tmp_[1] = velocity_[1];
  velocity_tmp_[2] = velocity_[2];
  velocity_tmp_[0].echange_espace_virtuel(ghost_size_pour_ibc_);
  velocity_tmp_[1].echange_espace_virtuel(ghost_size_pour_ibc_);
  velocity_tmp_[2].echange_espace_virtuel(ghost_size_pour_ibc_);
  couplage_tubes_ibc_.force_ibc(velocity_[0], velocity_[1], velocity_[2], rho_field_, timestep_, pressure_);

# endif

# if 0
  ///////////////////////////
  // dans velocity_ il y a le champ de vitesse au temps n+1, il est a div nulle
  // dans velocity_tmp il a le champ de vitesse au temps n+1 modifiee pour tenir compte des mailles qui vont entrer dans le solide. il n'est pas a div nulle
  velocity_tmp_[0] = velocity_[0];
  velocity_tmp_[1] = velocity_[1];
  velocity_tmp_[2] = velocity_[2];
  velocity_tmp_[0].echange_espace_virtuel(ghost_size_pour_ibc_);
  velocity_tmp_[1].echange_espace_virtuel(ghost_size_pour_ibc_);
  velocity_tmp_[2].echange_espace_virtuel(ghost_size_pour_ibc_);
  couplage_tubes_ibc_.forcage_anticipe(velocity_[0], velocity_[1], velocity_[2], rho_field_, timestep_, pressure_);
# endif

#if 0
  {
    const int kmax = rho_field_.nk();
    for (int k = 0; k < kmax; k++)
      {
        mass_solver_scalar(d_rho_field_, delta_z_local_, k);
        //euler_explicit_update(d_rho_field_, rho_field_, k);

        //mass_solver_scalar(d_fonction_levelset_, delta_z_local_, k);
        //euler_explicit_update(d_fonction_levelset_, fonction_levelset_, k);


#if 0
        // astuce temporaire en attendant de corriger le bug de l'operateur de convection scalaire sur les parois
        const int ni = rho_field_.ni();
        const int nj = rho_field_.nj();
        const double a = rho_liquide_;
        const double b = rho_gaz_;
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              rho_field_(i, j, k) = std::min(std::max(rho_field_(i, j, k), b), a);
            }
#endif
      }
  }
#endif

  rho_field_ = rho_field_nplus1;
  alpha_field_ = alpha_field_nplus1;

  if (Process::je_suis_maitre())
    {
      Cout << "Timings diff=" << statistiques().last_time(diffusion_counter_)
           << " conv=" << statistiques().last_time(convection_counter_)
           << " src=" << statistiques().last_time(source_counter_)
           << " proj=" << statistiques().last_time(projection_) << finl;
    }
}

// Hard coded constant pressure gradient in i direction, add contribution in m/s*volume of control volume
void IJK_problem_double::terme_source_gravite(IJK_Field_double& dv, int k_index, int dir) const
{
  const double constant = gravite_[dir];
  const int imax = dv.ni();
  const int jmax = dv.nj();
  for (int j = 0; j < jmax; j++)
    {
      for (int i = 0; i < imax; i++)
        {
          dv(i, j, k_index) += constant;
        }
    }
}

void IJK_problem_double::euler_explicit_update(const IJK_Field_double& dv, IJK_Field_double& v,
                                               const int k_layer) const
{
  const double delta_t = timestep_;
  const int imax = v.ni();
  const int jmax = v.nj();
  for (int j = 0; j < jmax; j++)
    {
      for (int i = 0; i < imax; i++)
        {
          double x = dv(i,j,k_layer);
          v(i,j,k_layer) += x * delta_t;
        }
    }
}


/////////////////////////////////////////
////////////////////////////////////////
// schema decentre amont TRITON pour rho

void IJK_problem_double::scalar_convection_op_Triton(const IJK_Field_double& rho_,
                                                     const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                                     IJK_Field_double& d_rho_)
{
  // il s'agit d'un schema decentre amont !!!!
  IJK_Field_double d_rho_field = d_rho_;
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  const double delta_x = geom.get_constant_delta(0);  // taille de la maille selon x
  const double delta_y = geom.get_constant_delta(1);
  const double delta_z = geom.get_constant_delta(2);
  // schema Despres Lagoutiere pour un fluide incompresible !!!!! si compressible il faut ajouter une terme supplementaire dans bilanx, bilany et bilanz
  // div (rho * u) = rho div u + u.Grad rho = u.Grad rho pour un fluide incompressible

  // convection selon x
  IJK_Field_double velocity;
  velocity= vx; // velocity est egale a la vitesse selon la direction 0,1 ou 2
  int imin = 0;
  int jmin = 0;
  int kmin = 0;
  int imax = rho_.ni(); // donne le nombre d'elements de rho selon x
  int jmax = rho_.nj();
  int kmax = rho_.nk();
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              double velocity_x_ijk = velocity(i, j, k);
              double velocity_x_suivante_ijk = velocity(i+1, j, k); // attention comment ca se passe a la frontiere a cause du i+1
              double tmp_rho_x_1;
              if ( velocity_x_suivante_ijk > 0 )
                {
                  tmp_rho_x_1 = rho_(i, j, k);
                }
              else tmp_rho_x_1 = rho_(i+1, j, k);// attention comment ca se passe a la frontiere a cause du i+1
              double tmp_rho_x_2;
              if ( velocity_x_ijk > 0 )
                {
                  tmp_rho_x_2 = rho_(i-1, j, k); // attention comment ca se passe a la frontiere a cause du i-1
                }
              else tmp_rho_x_2 = rho_(i, j, k);
              double bilan_x;
              bilan_x = (velocity_x_suivante_ijk * tmp_rho_x_1 - velocity_x_ijk * tmp_rho_x_2)* delta_y * delta_z; // si velocity(i, j, k)>0 en i et i+1 ca donne : f_x = rho(i,j,k) *velocity(i+1,j,k) - rho(i-1,j,k) *velocity(i,j,k), ce qui correspondant au bilan sur les volumes de controle ijk et i+1jk
              d_rho_field(i,j,k) =  bilan_x;
            }
        }
    }

// convection selon y
  velocity = vy ;// velocity est egale a la vitesse selon la direction 0,1 ou 2
  imax = rho_.ni(); // donne le nombre d'elements de rho selon x
  jmax = rho_.nj();
  kmax = rho_.nk();
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              double velocity_y_ijk = velocity(i, j, k);
              double velocity_y_suivante_ijk = velocity(i, j+1, k); // attention comment ca se passe a la frontiere a cause du i+1
              double tmp_rho_y_1;
              if ( velocity_y_suivante_ijk > 0 )
                {
                  tmp_rho_y_1 = rho_(i, j, k);
                }
              else tmp_rho_y_1 = rho_(i, j+1, k);// attention comment ca se passe a la frontiere a cause du i+1
              double tmp_rho_y_2;
              if ( velocity_y_ijk > 0 )
                {
                  tmp_rho_y_2 = rho_(i, j-1, k); // attention comment ca se passe a la frontiere a cause du i-1
                }
              else tmp_rho_y_2 = rho_(i, j, k);
              double bilan_y;
              bilan_y = (velocity_y_suivante_ijk * tmp_rho_y_1 - velocity_y_ijk * tmp_rho_y_2) * delta_x * delta_z; // si velocity(i, j, k)>0 en i et i+1 ca donne : f_x = rho(i,j,k) *velocity(i+1,j,k) - rho(i-1,j,k) *velocity(i,j,k), ce qui correspondant au bilan sur les volumes de controle ijk et i+1jk
              d_rho_field(i,j,k) +=  bilan_y;

            }
        }
    }

// convection selon z
  velocity =  vz; // velocity est egale a la vitesse selon la direction 0,1 ou 2
  imax = rho_.ni(); // donne le nombre d'elements de rho selon x
  jmax = rho_.nj();
  kmax = rho_.nk();
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              double velocity_z_ijk = velocity(i, j, k);
              double velocity_z_suivante_ijk = velocity(i, j, k+1); // attention comment ca se passe a la frontiere a cause du i+1
              double tmp_rho_z_1;
              if ( velocity_z_suivante_ijk > 0 )
                {
                  tmp_rho_z_1 = rho_(i, j, k);
                }
              else tmp_rho_z_1 = rho_(i, j, k+1);// attention comment ca se passe a la frontiere a cause du i+1
              double tmp_rho_z_2;
              if ( velocity_z_ijk > 0 )
                {
                  tmp_rho_z_2 = rho_(i, j, k-1); // attention comment ca se passe a la frontiere a cause du i-1
                }
              else tmp_rho_z_2 = rho_(i, j, k);
              double bilan_z;
              bilan_z = (velocity_z_suivante_ijk * tmp_rho_z_1 - velocity_z_ijk * tmp_rho_z_2) * delta_x * delta_y; // si velocity(i, j, k)>0 en i et i+1 ca donne : f_x = rho(i,j,k) *velocity(i+1,j,k) - rho(i-1,j,k) *velocity(i,j,k), ce qui correspondant au bilan sur les volumes de controle ijk et i+1jk
              d_rho_field(i,j,k) +=  bilan_z;
              d_rho_field_(i,j,k) = -d_rho_field(i,j,k);
            }
        }
    }
}



////////////////////////////////////////////////////
//////////////////////////////////////////////////////
// Despres type Samuel pour rho

inline double calcul_rho_amont_aval_samuel(const double i_cfl, const double rho_aval, const double rho_centre, const double rho_amont, const double v_amont, const double v_aval)
{
  // on raisonne en amont, centre et aval par rapport aux faces car c'est la qu'on fait les bilans sur rho
  double max_rho = std::max(rho_centre, rho_amont);
  double min_rho = std::min(rho_centre, rho_amont);
  double max_consistance = std::max(rho_aval, rho_centre);
  double min_consistance = std::min(rho_aval, rho_centre);
  double b_i = (rho_centre - max_rho) * i_cfl / v_aval  + max_rho * v_amont / v_aval + rho_centre * (1 - v_amont / v_aval);
  double B_i = (rho_centre - min_rho) * i_cfl / v_aval  + min_rho * v_amont / v_aval + rho_centre * (1 - v_amont / v_aval);
  double borne_inf = std::max(min_consistance, b_i);
  double borne_sup = std::min(max_consistance, B_i);

  // dans Despres lagouttiere les tests sont fais sur rho_aval par rapport a la face consideree
  // resu est la valeur de la masse volumique sur la face consideree
  // Reecriture de :
  //  if (rho_aval < b_i) {
  //  resu = b_i;
  //} else if (rho_aval > B_i) {
  //  resu = B_i;
  //} else {
  //  resu = rho_aval;
  //}


  return std::min(borne_sup, std::max(borne_inf, rho_aval));
}

#define IJK_DIR(x) i+((direction==DIRECTION_I)?x:0), j+((direction==DIRECTION_J)?x:0), k+((direction==DIRECTION_K)?x:0)
// ((direction==DIRECTION_I)?x:0) signifie : si la direction==DIRECTION_I alors x, sinon 0 : donc si la direction==DIRECTION_I on fait i+x sinon i+0
inline double calcul_rho_amont_aval2_samuel(const IJK_Field_double& velocity, const double cfl, const double inv_delta,
                                            const IJK_Field_double& rho, int i, int j, int k, int direction,
                                            double& conv_rho)
{
  // on a besoin de connaitre les valeurs de la vitesse sur la face i-3/2 : v_GG; i-0.5 : v_G et i+0.5 : v_D pour calculer la valeur de rho sur la face i-0.5 : rho_G
  // on a besoin de connaitre les valeurs de la vitesse sur la face i-0.5 : v_G; i+0.5 : v_D et i+3/2 : v_DD pour calculer la valeur de rho sur la face i+0.5 : rho_D
  const double v_GG = velocity(IJK_DIR(-1));
  const double v_G = velocity(IJK_DIR(0));
  const double v_D = velocity(IJK_DIR(1));
  const double v_DD = velocity(IJK_DIR(2));
  //
  double rho_face_G =0;
  double rho_face_D =0;
  double icfl;
  if ( v_G > 0)
    {
      if ( v_GG >0 )
        {
          icfl = 1. / cfl;
          rho_face_G = calcul_rho_amont_aval_samuel(icfl, rho(IJK_DIR(0)), rho(IJK_DIR(-1)), rho(IJK_DIR(-2)), velocity(IJK_DIR(-1)), velocity(IJK_DIR(0)));
        }
      else  rho_face_G = rho(IJK_DIR(-1));
    }
  else if ( v_G <0 )
    {
      if ( v_D <0   )
        {
          icfl = 1. / cfl;
          rho_face_G = calcul_rho_amont_aval_samuel(-icfl, rho(IJK_DIR(-1)), rho(IJK_DIR(0)), rho(IJK_DIR(1)), velocity(IJK_DIR(1)), velocity(IJK_DIR(0)));
        }
      else rho_face_G = rho(IJK_DIR(0));
    }
  else if (v_G == 0)
    {
      rho_face_G =0;
    }

  if ( v_D > 0)
    {
      if ( v_DD >0 )
        {
          icfl = 1. / cfl;
          rho_face_D = calcul_rho_amont_aval_samuel(icfl, rho(IJK_DIR(1)), rho(IJK_DIR(0)), rho(IJK_DIR(-1)), velocity(IJK_DIR(0)), velocity(IJK_DIR(1)));
        }
      else  rho_face_D = rho(IJK_DIR(0));
    }
  else if ( v_D <0 )
    {
      if ( v_DD <0   )
        {
          icfl = 1. / cfl;
          rho_face_D = calcul_rho_amont_aval_samuel(-icfl, rho(IJK_DIR(0)), rho(IJK_DIR(1)), rho(IJK_DIR(2)), velocity(IJK_DIR(2)), velocity(IJK_DIR(1)));
        }
      else rho_face_D = rho(IJK_DIR(1));
    }
  else if (v_D == 0)
    {
      rho_face_D =0;
    }

  return (v_D * rho_face_D - v_G * rho_face_G) * inv_delta - rho(i,j,k) * (v_D - v_G) * inv_delta; // bilan de rho sur la maille i,j,k ; inv_delta = surface / volume

}
void IJK_problem_double::scalar_convection_op_Despres_Samuel(IJK_Field_double& rho_,
                                                             const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz)
{
  const double dt = timestep_;
  IJK_Field_double rho_tmp;
  rho_tmp.allocate(splitting_, IJK_Splitting::ELEM, 2);
  rho_tmp.data() = 1.; // pourquoi on a besoin de ca???
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  const double delta_x = geom.get_constant_delta(0);  // taille de la maille selon x
  double conv_rho =0;
  int imin = 0;
  int jmin = 0;
  int kmin = 0;
  int imax = rho_.ni(); // donne le nombre d'elements de rho selon x il y en a n donc imax=n-1
  int jmax = rho_.nj();
  int kmax = rho_.nk();

  // convection selon x
  const double CFL_x = dt / delta_x;
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              conv_rho = calcul_rho_amont_aval2_samuel(vx, CFL_x, 1./delta_x, rho_, i, j, k, DIRECTION_I, conv_rho);
              rho_tmp(i,j,k) = rho_(i,j,k) - dt * conv_rho;
            }
        }
    }

  rho_ = rho_tmp;
  rho_.echange_espace_virtuel(2);

  const double delta_y = geom.get_constant_delta(1);
  const double CFL_y = dt / delta_y;
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              conv_rho = calcul_rho_amont_aval2_samuel(vy, CFL_y, 1./delta_y, rho_, i, j, k, DIRECTION_J, conv_rho);
              rho_tmp(i,j,k) = rho_(i,j,k) - dt * conv_rho;
            }
        }
    }

  rho_ = rho_tmp;
  rho_.echange_espace_virtuel(2);

  const double delta_z = geom.get_constant_delta(2);
  const double CFL_z = dt / delta_z;
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              conv_rho = calcul_rho_amont_aval2_samuel(vz, CFL_z, 1./delta_z, rho_, i, j, k, DIRECTION_K, conv_rho);
              rho_tmp(i,j,k) = rho_(i,j,k) - dt * conv_rho;
            }
        }
    }

  rho_ = rho_tmp;

}


///////////////////////////////////////////////
//////////////////////////////////////////////////
// Desppres type delacherie pour rho

inline double calcul_rho_amont_aval(const double i_v_cfl, const double rho_aval, const double rho_centre, const double rho_amont)
{
  // on raisonne en amont, centre et aval par rapport aux faces car c'est la qu'on fait les bilans sur rho
  double max_rho = std::max(rho_centre, rho_amont);
  double min_rho = std::min(rho_centre, rho_amont);
  double max_consistance = std::max(rho_aval, rho_centre);
  double min_consistance = std::min(rho_aval, rho_centre);
  double b_i = (rho_centre - max_rho) * i_v_cfl + max_rho;
  double B_i = (rho_centre - min_rho) * i_v_cfl + min_rho;
  double borne_inf = std::max(min_consistance, b_i);
  double borne_sup = std::min(max_consistance, B_i);
  // dans Despres lagouttiere les tests sont fais sur rho_aval par rapport a la face consideree
  // resu est la valeur de la masse volumique sur la face consideree
  // Reecriture de :
  //  if (rho_aval < b_i) {
  //  resu = b_i;
  //} else if (rho_aval > B_i) {
  //  resu = B_i;
  //} else {
  //  resu = rho_aval;
  //}

  return std::min(borne_sup, std::max(borne_inf, rho_aval));
}

#define IJK_DIR(x) i+((direction==DIRECTION_I)?x:0), j+((direction==DIRECTION_J)?x:0), k+((direction==DIRECTION_K)?x:0)
// ((direction==DIRECTION_I)?x:0) signifie : si la direction==DIRECTION_I alors x, sinon 0 : donc si la direction==DIRECTION_I on fait i+x sinon i+0
inline double calcul_rho_amont_aval2(const IJK_Field_double& velocity, const double cfl, const double inv_delta,
                                     const IJK_Field_double& rho,
                                     int i, int j, int k, int direction)
{
  const double v = (velocity(IJK_DIR(0)) + velocity(IJK_DIR(1))) * 0.5; // on prend la moyenne de la vitesse dans la maille i,j,k comme vitesse de convection : NON CONSERVATIF
  double rho_face_G;
  double rho_face_D;
  double ivcfl;
  // La face droite pour v>0 est equivalente a la face gauche pour v<0.  La face gauche pour v>0 est equivalente a la face droite pour v<0.
  if (v > 0)
    {
      ivcfl = 1. / (v * cfl);
      rho_face_G = calcul_rho_amont_aval(ivcfl, rho(IJK_DIR(0)), rho(IJK_DIR(-1)), rho(IJK_DIR(-2)));
      rho_face_D = calcul_rho_amont_aval(ivcfl, rho(IJK_DIR(1)), rho(IJK_DIR(0)),  rho(IJK_DIR(-1))); // idem que rho_face_G avec un indice decale en aval donc +1
    }
  else if (v < 0)
    {
      ivcfl = 1. / (v * cfl);
      rho_face_G = calcul_rho_amont_aval(-ivcfl, rho(IJK_DIR(-1)), rho(IJK_DIR(0)), rho(IJK_DIR(1))); // idem que rho_face_D avec un indice decale en aval donc -1
      rho_face_D = calcul_rho_amont_aval(-ivcfl, rho(IJK_DIR(0)), rho(IJK_DIR(1)), rho(IJK_DIR(2)));
    }
  else
    {
      rho_face_G = rho_face_D = 0.;
    }
  return v * (rho_face_D - rho_face_G) * inv_delta; // bilan de rho sur la maille i,j,k ; inv_delta = surface / volume

}
void IJK_problem_double::scalar_convection_op_Despres_Dellacherie(IJK_Field_double& rho_,
                                                                  const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz)
{
  const double dt = timestep_;
  IJK_Field_double rho_tmp;
  rho_tmp.allocate(splitting_, IJK_Splitting::ELEM, 2);
  rho_tmp.data() = 1.; // pourquoi on a besoin de ca???
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  const double delta_x = geom.get_constant_delta(0);  // taille de la maille selon x
  const double delta_y = geom.get_constant_delta(1);
  const double delta_z = geom.get_constant_delta(2);
  // schema Despres Lagoutiere pour un fluide incompresible !!!!!
  // div (rho * u) = rho div u + u.Grad rho = u.Grad rho pour un fluide incompressible
  int imin = 0;
  int jmin = 0;
  int kmin = 0;
  int imax = rho_.ni(); // donne le nombre d'elements de rho selon x il y en a n donc imax=n-1
  int jmax = rho_.nj();
  int kmax = rho_.nk();
  // convection selon x

  const double CFL_x = dt / delta_x;
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              double drho_dt = - calcul_rho_amont_aval2(vx, CFL_x, 1./delta_x, rho_, i, j, k, DIRECTION_I);
              rho_tmp(i,j,k) = rho_(i,j,k) + dt * drho_dt; // on met a jour rho avec la convection selon x
            }
        }
    }

  rho_ = rho_tmp;
  rho_.echange_espace_virtuel(2);

  const double CFL_y = dt / delta_y;
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              double drho_dt = - calcul_rho_amont_aval2(vy, CFL_y, 1./delta_y, rho_, i, j, k, DIRECTION_J);
              rho_tmp(i,j,k) = rho_(i,j,k) + dt * drho_dt; // on met a jour rho avec la convection selon y
            }
        }
    }

  rho_ = rho_tmp;
  rho_.echange_espace_virtuel(2);

  const double CFL_z = dt / delta_z;
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              double drho_dt = - calcul_rho_amont_aval2(vz, CFL_z, 1./delta_z, rho_, i, j, k, DIRECTION_K);
              rho_tmp(i,j,k) = rho_(i,j,k) + dt * drho_dt; // on met a jour rho avec la convection selon z donne donc rho au temps n+1 !
            }
        }
    }

  rho_ = rho_tmp;

}

/////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Despres type Samuel pour le scalaire alpha

inline double calcul_alpha_amont_aval(const double i_cfl, const double alpha_aval, const double alpha_centre, const double alpha_amont, const double rho_aval, const double rho_centre, const double rho_amont, const double v_amont, const double v_aval)
{
  Cout << "rho_aval : " << rho_aval << " rho_centre : " << rho_centre << " rho_amont : " << rho_amont << endl;
  Cout << "alpha_aval : " << alpha_aval << " alpha_centre : " << alpha_centre << " alpha_amont : " << alpha_amont << endl;
  // on raisonne en amont, centre et aval par rapport aux faces car c'est la qu'on fait les bilans sur alpha
  double max_alpha = std::max(alpha_centre, alpha_amont);
  double min_alpha = std::min(alpha_centre, alpha_amont);
  Cout << "min_alpha : " << min_alpha << " max_alpha : " << max_alpha << endl;
  // CONSISTANCE DU FLUX ALPHA
  double max_alpha_consistance = std::max(alpha_aval, alpha_centre);
  double min_alpha_consistance = std::min(alpha_aval, alpha_centre);

  //
  double inv_v_alpha = 1 / v_aval;
  double v_amont_v_valpha = inv_v_alpha * v_amont;
  double rho_L = 1000; //rho_O
  double rho_G = 500; // rho_1

  double max_rho = std::max(rho_centre, rho_amont);
  double min_rho = std::min(rho_centre, rho_amont);
  // CONSISTANCE DU FLUX RHO
  double max_rho_consistance = std::max(rho_aval, rho_centre);
  double min_rho_consistance = std::min(rho_aval, rho_centre);
  Cout << "min_alpha_consistance : " << min_alpha_consistance << " max_alpha_consistance : " << max_alpha_consistance << " min_rho_consistance : " << min_rho_consistance << " max_rho_consistance : " << max_rho_consistance << endl;
  //
  // intervalle de confiance sur alpha pour la consistance de rho
  double diff_rho = 1. / (rho_L - rho_G); //
  double d_i = (min_rho_consistance - rho_G) *diff_rho;
  double D_i = (max_rho_consistance - rho_G) *diff_rho;
  Cout << "d_i : " << d_i << " D_i : " << D_i << " diff_rho " << diff_rho << endl;
  // intervalle de confiance sur alpha pour la stabilite de rho
  double coeff_a_i = rho_G*max_alpha + (1-min_alpha)*rho_L;
  double coeff_A_i = rho_G*min_alpha + (1-max_alpha)*rho_L;
  double a_i_tmp = (rho_centre - max_rho) * i_cfl * inv_v_alpha + max_rho * v_amont_v_valpha + rho_centre * (1 - v_amont / v_aval);
  double A_i_tmp = (rho_centre - min_rho) * i_cfl * inv_v_alpha + min_rho * v_amont_v_valpha + rho_centre * (1 - v_amont / v_aval);
  double a_i = (a_i_tmp - rho_G) * diff_rho;
  double A_i = (A_i_tmp - rho_G) * diff_rho;
  Cout << "a_i : " << a_i << " A_i : " << A_i << " (rho_centre - min_rho) * i_cfl * inv_v_alpha : " << (rho_centre - min_rho) * i_cfl * inv_v_alpha << " min_rho * v_amont_v_valpha : " << min_rho * v_amont_v_valpha << " v_amont_v_valpha " << v_amont_v_valpha << endl;
  // intervalle de confiance sur alpha pour la stabilite de alpha
  double b_i = (alpha_centre - max_alpha) * i_cfl / v_aval  + max_alpha * v_amont / v_aval + alpha_centre * (1 - v_amont / v_aval);
  double B_i = (alpha_centre - min_alpha) * i_cfl / v_aval  + min_alpha * v_amont / v_aval + alpha_centre * (1 - v_amont / v_aval);
  Cout << "b_i : " << b_i << " B_i : " << B_i << " coeff_a_i : " << coeff_a_i << " coeff_A_i : " << coeff_A_i << " a_i_tmp : " << a_i_tmp << " a_i_tmp - rho_G  " << a_i_tmp - rho_G  << endl;

  // pour empeche D_i et A_i de devenir <0 ....
  if ( D_i <0)
    {
      D_i=0;
    }
  if ( A_i <0)
    {
      A_i=0;
    }

  // dans Despres lagouttiere les tests sont fais sur alpha_aval par rapport a la face consideree
  // resu est la valeur de alpha sur la face consideree
  // Reecriture de :
  //  if (alpha_aval < b_i) {
  //  resu = b_i;
  //} else if (alpha_aval > B_i) {
  //  resu = B_i;
  //} else {
  //  resu = alpha_aval;
  //}

  // determination de l'intervalle total :

  double borne_inf_2 = std::max(min_alpha_consistance, std::max(d_i, std::max(a_i, b_i)));
  double borne_sup_2 = std::min(max_alpha_consistance, std::min(D_i, std::min(A_i, B_i)));
  Cout << "borne_inf_2 = std::max(min_alpha_consistance, std::max(d_i, std::max(a_i, b_i))) : " << borne_inf_2 << " borne_sup_2 = std::min(max_alpha_consistance, std::min(D_i, std::min(A_i, B_i))) : " << borne_sup_2 << endl;
  // borne_inf_2 = borne_sup_2 =0;


  double borne_inf = std::max(min_alpha_consistance, b_i);
  double borne_sup = std::min(max_alpha_consistance, B_i);

  double resu = std::min(borne_sup_2, std::max(borne_inf_2, alpha_aval));
# if 0
  if ( resu<0)
    {
      resu=0;
    }
  if ( resu>1)
    {
      resu=1;
    }
#endif

  Cout << "borne_inf = std::max(min_alpha_consistance, b_i) : " << borne_inf << " borne_sup = std::min(max_alpha_consistance, B_i) : " << borne_sup << endl;
  Cout << " min(borne_sup_2, max(borne_inf_2, alpha_aval)) : " << std::min(borne_sup_2, std::max(borne_inf_2, alpha_aval)) <<" min(B_i, max(b_i, alpha_aval)) : " << std::min(B_i, std::max(b_i, alpha_aval)) << " resu : " << resu << endl;
  Cout << "                                                  " << endl;
  return resu;
  //return std::min(B_i, std::max(b_i, alpha_aval));



}

#define IJK_DIR(x) i+((direction==DIRECTION_I)?x:0), j+((direction==DIRECTION_J)?x:0), k+((direction==DIRECTION_K)?x:0)
// ((direction==DIRECTION_I)?x:0) signifie : si la direction==DIRECTION_I alors x, sinon 0 : donc si la direction==DIRECTION_I on fait i+x sinon i+0
inline void calcul_alpha_amont_aval2(const IJK_Field_double& velocity, const double cfl, const double inv_delta,
                                     const IJK_Field_double& alpha, const IJK_Field_double& rho,
                                     int i, int j, int k, int direction,
                                     double& conv_alpha, double& conv_rho)
{

  const double v_GG = velocity(IJK_DIR(-1));
  const double v_G = velocity(IJK_DIR(0));
  const double v_D = velocity(IJK_DIR(1));
  const double v_DD = velocity(IJK_DIR(2));
  Cout << "i :" << i << " j :" << j << " k :" << k << " v_GG :";
  Cout  << v_GG << " v_G :" << v_G <<  " v_D :" << v_D << " v_DD :" << v_DD << endl;
  double alpha_face_G = 0.;
  double alpha_face_D = 0.;
  double rho_face_G = 0.;
  double rho_face_D = 0.;
  double icfl;
  double rho_L = 1000;
  double rho_G = 500;
  // La face droite pour v>0 est equivalente a la face gauche pour v<0.  La face gauche pour v>0 est equivalente a la face droite pour v<0.
  // Decentrement amont pour rho
  // Despres pour alpha
  icfl = 1. / cfl;
  // calcul d'alpha en i+0.5 donc a droite
  Cout << "Calcul de alpha_face_D" << endl;
  if ( v_D >0  )
    {
      if (v_G >0)
        {
          alpha_face_D = calcul_alpha_amont_aval(icfl, alpha(IJK_DIR(1)), alpha(IJK_DIR(0)),  alpha(IJK_DIR(-1)), rho(IJK_DIR(1)), rho(IJK_DIR(0)),  rho(IJK_DIR(-1)), velocity(IJK_DIR(0 )), velocity(IJK_DIR(1))); // idem que alpha_face_G avec un indice decale en aval donc +1

          rho_face_D   = alpha_face_D * rho_L  + (1-alpha_face_D)* rho_G;
        }
      else if (v_G < 0)
        {
          alpha_face_D = alpha(IJK_DIR(0));
          rho_face_D   = alpha_face_D * rho_L  + (1-alpha_face_D)* rho_G;
        }
    }
  if ( v_D<0)
    {
      if ( v_DD<0)
        {
          alpha_face_D = calcul_alpha_amont_aval(-icfl, alpha(IJK_DIR(0)),  alpha(IJK_DIR(1)), alpha(IJK_DIR(2)), rho(IJK_DIR(0)),  rho(IJK_DIR(1)), rho(IJK_DIR(2)), velocity(IJK_DIR(2)), velocity(IJK_DIR(1)));
          rho_face_D   = alpha_face_D * rho_L  + (1-alpha_face_D)* rho_G;
        }
      else if ( v_DD>0)
        {
          alpha_face_D = alpha(IJK_DIR(1));
          rho_face_D   = alpha_face_D * rho_L  + (1-alpha_face_D)* rho_G;
        }

    }

  // calcul d'alpha en i-0.5 donc a gauche
  Cout << "Calcul de alpha_face_G" << endl;
  if ( v_G>0)
    {
      if (v_GG >0)
        {
          alpha_face_G = calcul_alpha_amont_aval(icfl, alpha(IJK_DIR(0)), alpha(IJK_DIR(-1)), alpha(IJK_DIR(-2)), rho(IJK_DIR(0)), rho(IJK_DIR(-1)), rho(IJK_DIR(-2)), velocity(IJK_DIR(-1)), velocity(IJK_DIR(0)));
          rho_face_G   = alpha_face_G * rho_L + (1-alpha_face_G) * rho_G;
        }
      else if (v_GG<0)
        {
          alpha_face_G = alpha(IJK_DIR(-1));
          rho_face_G   = alpha_face_G * rho_L + (1-alpha_face_G) * rho_G;
        }
    }

  if (v_G<0)
    {
      if ( v_D<0)
        {
          alpha_face_G = calcul_alpha_amont_aval(-icfl, alpha(IJK_DIR(-1)), alpha(IJK_DIR(0)), alpha(IJK_DIR(1)), rho(IJK_DIR(-1)), rho(IJK_DIR(0)), rho(IJK_DIR(1)), velocity(IJK_DIR(1)), velocity(IJK_DIR(0)));
          rho_face_G   = alpha_face_G * rho_L + (1-alpha_face_G) * rho_G;
        }

      else if ( v_D>0)
        {
          alpha_face_G = alpha(IJK_DIR(0));
          rho_face_G   = alpha_face_G * rho_L + (1-alpha_face_G) * rho_G;
        }

    }




# if 0
  if (  (v_G > 0) && (v_D > 0) )
    {
      icfl = 1. / cfl;
      alpha_face_G = calcul_alpha_amont_aval(icfl, alpha(IJK_DIR(0)), alpha(IJK_DIR(-1)), alpha(IJK_DIR(-2)), rho(IJK_DIR(0)), rho(IJK_DIR(-1)), rho(IJK_DIR(-2)), velocity(IJK_DIR(-1)), velocity(IJK_DIR(0)));
      alpha_face_D = calcul_alpha_amont_aval(icfl, alpha(IJK_DIR(1)), alpha(IJK_DIR(0)),  alpha(IJK_DIR(-1)), rho(IJK_DIR(1)), rho(IJK_DIR(0)),  rho(IJK_DIR(-1)), velocity(IJK_DIR(0 )), velocity(IJK_DIR(1))); // idem que alpha_face_G avec un indice decale en aval donc +1

      //rho_face_G   = alpha_face_G * rho(IJK_DIR(-1)) + (1-alpha_face_G) * rho(IJK_DIR(-1));
      //rho_face_D   = alpha_face_D * rho(IJK_DIR(0))  + (1-alpha_face_D) * rho(IJK_DIR(0));

      rho_face_G   = alpha_face_G * rho_L + (1-alpha_face_G) * rho_G;
      rho_face_D   = alpha_face_D * rho_L  + (1-alpha_face_D)* rho_G;
      //Cout << " rho_faceG : " << rho_face_G << " rho_faceD : " << rho_face_D << endl;

    }
  else if (  (v_G < 0) && (v_D < 0) )
    {
      icfl = 1. / cfl;
      alpha_face_G = calcul_alpha_amont_aval(-icfl, alpha(IJK_DIR(-1)), alpha(IJK_DIR(0)), alpha(IJK_DIR(1)), rho(IJK_DIR(-1)), rho(IJK_DIR(0)), rho(IJK_DIR(1)), velocity(IJK_DIR(1)), velocity(IJK_DIR(0))); // idem que alpha_face_D avec un indice decale en aval donc -1
      alpha_face_D = calcul_alpha_amont_aval(-icfl, alpha(IJK_DIR(0)),  alpha(IJK_DIR(1)), alpha(IJK_DIR(2)), rho(IJK_DIR(0)),  rho(IJK_DIR(1)), rho(IJK_DIR(2)), velocity(IJK_DIR(2)), velocity(IJK_DIR(1)));

      rho_face_G   = alpha_face_G * rho_L + (1-alpha_face_G) * rho_G;
      rho_face_D   = alpha_face_D * rho_L  + (1-alpha_face_D)* rho_G;

    }
  else if (  (v_G < 0) && (v_D > 0) )
    {
      alpha_face_G = alpha(IJK_DIR(0));
      alpha_face_D = alpha(IJK_DIR(0));

      rho_face_G   = alpha_face_G * rho_L + (1-alpha_face_G) * rho_G;
      rho_face_D   = alpha_face_D * rho_L  + (1-alpha_face_D)* rho_G;



    }
  else if (  (v_G > 0) && (v_D < 0) )
    {
      alpha_face_G = alpha(IJK_DIR(-1));
      alpha_face_D = alpha(IJK_DIR(1));

      rho_face_G   = alpha_face_G * rho_L + (1-alpha_face_G) * rho_G;
      rho_face_D   = alpha_face_D * rho_L  + (1-alpha_face_D)* rho_G;


    }
  else if (  (v_G == 0) && (v_D > 0) )
    {
      alpha_face_G = 0;
      alpha_face_D = alpha(IJK_DIR(0));

      rho_face_G = 0;
      rho_face_D   = alpha_face_D * rho_L  + (1-alpha_face_D)* rho_G;

    }
  else if (  (v_G == 0) && (v_D < 0) )
    {
      alpha_face_G = 0;
      alpha_face_D = alpha(IJK_DIR(1));

      rho_face_G = 0;
      rho_face_D   = alpha_face_D * rho_L  + (1-alpha_face_D)* rho_G;

    }
  else if (  (v_G > 0) && (v_D == 0) )
    {
      alpha_face_G = alpha(IJK_DIR(-1));
      alpha_face_D = 0;


      rho_face_G   = alpha_face_G * rho_L + (1-alpha_face_G) * rho_G;
      rho_face_D = 0;
    }
  else if (  (v_G < 0) && (v_D == 0) )
    {
      alpha_face_G = alpha(IJK_DIR(0));
      alpha_face_D = 0;


      rho_face_G   = alpha_face_G * rho_L + (1-alpha_face_G) * rho_G;
      rho_face_D = 0;
    }
  else
    {
      alpha_face_G = alpha_face_D = 0.;
      rho_face_G   = rho_face_D   = 0.;
    }
#endif

  conv_alpha = ((v_D * alpha_face_D - v_G * alpha_face_G) - alpha(i,j,k) * (v_D - v_G) ) * inv_delta;
  conv_rho   = ( v_D * rho_face_D   - v_G * rho_face_G) * inv_delta ;

}



void IJK_problem_double::scalar_convection_op_DESPRES_alpha_AMONT_rho(IJK_Field_double& alpha_, IJK_Field_double& rho_,
                                                                      const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz)
{
  const double dt = timestep_;
  IJK_Field_double alpha_tmp;
  alpha_tmp.allocate(splitting_, IJK_Splitting::ELEM, 2);
  alpha_tmp.data() = 1.; // pourquoi on a besoin de ca???
  IJK_Field_double rho_tmp;
  rho_tmp.allocate(splitting_, IJK_Splitting::ELEM, 2);
  rho_tmp.data() = 1.; // pourquoi on a besoin de ca???

  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  const double delta_x = geom.get_constant_delta(0);  // taille de la maille selon x

  double conv_alpha=0;
  double conv_rho=0;

  int imin = 0;
  int jmin = 0;
  int kmin = 0;
  int imax = alpha_.ni(); // rho_field_ et alpha_field_ ont la meme taille
  int jmax = alpha_.nj();
  int kmax = alpha_.nk();
  // on fait du splitting
  // convection selon x
  const double CFL_x = dt / delta_x;
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              calcul_alpha_amont_aval2(vx, CFL_x, 1./delta_x, alpha_, rho_, i, j, k, DIRECTION_I, conv_alpha, conv_rho);
              alpha_tmp(i,j,k) = alpha_(i,j,k) - dt * conv_alpha;
              rho_tmp(i,j,k) = rho_(i,j,k) - dt * conv_rho;
            }
        }
    }
  rho_ = rho_tmp;
  alpha_ = alpha_tmp;

  rho_.echange_espace_virtuel(2);
  alpha_.echange_espace_virtuel(2);


  const double delta_y = geom.get_constant_delta(1);
  const double CFL_y = dt / delta_y;
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              calcul_alpha_amont_aval2(vy, CFL_y, 1./delta_y, alpha_, rho_, i, j, k, DIRECTION_J, conv_alpha, conv_rho);
              alpha_tmp(i,j,k) = alpha_(i,j,k) - dt * conv_alpha;
              rho_tmp(i,j,k) = rho_(i,j,k) - dt * conv_rho;
            }
        }
    }
  rho_ = rho_tmp;
  alpha_ = alpha_tmp;
  rho_.echange_espace_virtuel(2);
  alpha_.echange_espace_virtuel(2);
  const double delta_z = geom.get_constant_delta(2);
  const double CFL_z = dt / delta_z;
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              calcul_alpha_amont_aval2(vz, CFL_z, 1./delta_z, alpha_, rho_, i, j, k, DIRECTION_K, conv_alpha, conv_rho);
              alpha_tmp(i,j,k) = alpha_(i,j,k) - dt * conv_alpha;
              rho_tmp(i,j,k) = rho_(i,j,k) - dt * conv_rho;
            }
        }
    }
  rho_ = rho_tmp;
  alpha_ = alpha_tmp;

# if 0


  for (int j = jmin; j < jmax; j++)
    {
      for (int i = imin; i < imax; i++)
        {
# if 0
          alpha_(i,j,kmax-1) = alpha_(i,j,kmax-2);
          alpha_(i,j,kmax) = alpha_(i,j,kmax-2);
          alpha_(i,j,kmax+1) = alpha_(i,j,kmax-2);
          alpha_(i,j,0) = alpha_(i,j,1);
          alpha_(i,j,-1) = alpha_(i,j,1);
          alpha_(i,j,-2) = alpha_(i,j,1);
#endif
          //alpha_(i,j,kmax-3) = 1;
          alpha_(i,j,kmax-2) = 1;
          rho_field_(i,j,kmax-2) = 1000;
          rho_field_(i,j,kmax-1) = 1000;
          rho_field_(i,j,kmax) = 1000;
          rho_field_(i,j,kmax+1) = 1000;
          alpha_(i,j,kmax-1) = 1;
          alpha_(i,j,kmax) = 1;
          alpha_(i,j,kmax+1) = 1;
          //alpha_(i,j,2) = 1;
          alpha_(i,j,1) = 1;
          alpha_(i,j,0) = 1;
          alpha_(i,j,-1) = 1;
          alpha_(i,j,-2) = 1;
          rho_field_(i,j,1) = 1000;
          rho_field_(i,j,0) = 1000;
          rho_field_(i,j,-1) = 1000;
          rho_field_(i,j,-2) = 1000;
        }
    }

  for (int k = kmin; k < kmax; k++)
    {
      for (int j =jmin; j < jmax; j++)
        {
          rho_field_(0,j,k) = 1000;
          rho_field_(imax-1,j,k) = 1000;
        }
    }
#endif

}

void IJK_problem_double::calculer_rho2(const IJK_Field_double& alpha_, IJK_Field_double& rho_field_bis2)

{


  int imin = 0;
  int jmin = 0;
  int kmin = 0;
  int imax = alpha_.ni(); // rho_field_ et alpha_ ont la meme taille
  int jmax = alpha_.nj();
  int kmax = alpha_.nk();
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              rho_field_bis2(i,j,k) = alpha_(i,j,k) * 1000 + (1 - alpha_(i,j,k)) * 500;
            }
        }
    }

}



/////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// convection de rho ac un schema decentre amont conservatif





// schema decentre amont pour rho
#define IJK_DIR(x) i+((direction==DIRECTION_I)?x:0), j+((direction==DIRECTION_J)?x:0), k+((direction==DIRECTION_K)?x:0)
// ((direction==DIRECTION_I)?x:0) signifie : si la direction==DIRECTION_I alors x, sinon 0 : donc si la direction==DIRECTION_I on fait i+x sinon i+0
inline double calcul_rho_amont(const IJK_Field_double& velocity, const double inv_delta,
                               const IJK_Field_double& rho,
                               int i, int j, int k, int direction)
{
  const double v_G = velocity(IJK_DIR(0));
  const double v_D = velocity(IJK_DIR(1));
  double rho_face_G =0;
  double rho_face_D =0;
  // decentrement amont
  if ( v_G > 0)
    {
      rho_face_G = rho(IJK_DIR(-1));
    }
  else if ( v_G <0 )
    {
      rho_face_G = rho(IJK_DIR(0));
    }
  else if (v_G == 0)
    {
      rho_face_G =0;
    }

  if ( v_D > 0)
    {
      rho_face_D = rho(IJK_DIR(0));
    }
  else if ( v_D <0 )
    {
      rho_face_D = rho(IJK_DIR(1));
    }
  else if (v_D == 0)
    {
      rho_face_D =0;
    }

  return (v_D * rho_face_D - v_G * rho_face_G) * inv_delta ; // bilan de rho sur la maille i,j,k ; inv_delta = surface / volume

}



void IJK_problem_double::scalar_convection_op_Amont_rho(IJK_Field_double& rho_,
                                                        const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz)
{
  const double dt = timestep_;
  IJK_Field_double rho_tmp;
  rho_tmp.allocate(splitting_, IJK_Splitting::ELEM, 2);
  rho_tmp.data() = 1.; // pourquoi on a besoin de ca???
  const IJK_Splitting& splitting = vx.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  const double delta_x = geom.get_constant_delta(0);  // taille de la maille selon x

  int imin = 0;
  int jmin = 0;
  int kmin = 0;
  int imax = rho_.ni();
  int jmax = rho_.nj();
  int kmax = rho_.nk();
  // on fait du splitting
  // convection selon x
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              double d_rho_dt = - calcul_rho_amont(vx, 1./delta_x, rho_, i, j, k, DIRECTION_I);
              rho_tmp(i,j,k) = rho_(i,j,k) + dt * d_rho_dt;
            }
        }
    }

  rho_ = rho_tmp;


  rho_.echange_espace_virtuel(2);
  const double delta_y = geom.get_constant_delta(1);
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              double d_rho_dt = - calcul_rho_amont(vy, 1./delta_y, rho_, i, j, k, DIRECTION_J);
              rho_tmp(i,j,k) = rho_(i,j,k) + dt * d_rho_dt;
            }
        }
    }

  rho_ = rho_tmp;
  rho_.echange_espace_virtuel(2);
  const double delta_z = geom.get_constant_delta(2);
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              double d_rho_dt = - calcul_rho_amont(vz, 1./delta_z, rho_, i, j, k, DIRECTION_K);
              rho_tmp(i,j,k) = rho_(i,j,k) + dt * d_rho_dt;
            }
        }
    }

  rho_ = rho_tmp;

}


void IJK_problem_double::calcul_integrale_rho( IJK_Field_double& rho_) const
{
  const IJK_Splitting& splitting = rho_.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  const double volume = geom.get_constant_delta(0) * geom.get_constant_delta(1) * geom.get_constant_delta(2);  // taille de la maille selon x
  double tmp_rho_domaine =0;
  int imin = 0;
  int jmin = 0;
  int kmin = 0;
  int imax = rho_.ni();
  int jmax = rho_.nj();
  int kmax = rho_.nk();
  double Vol_domaine = imax * jmax * kmax * volume;;
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              tmp_rho_domaine += rho_(i, j, k) * volume;
            }
        }
    }
  double rho_domaine = tmp_rho_domaine;
  mp_sum(rho_domaine);
  Cout << "Rho_domaine : " << rho_domaine << " t = " << current_time_ << " Vol_domaine " << Vol_domaine << endl;
}
