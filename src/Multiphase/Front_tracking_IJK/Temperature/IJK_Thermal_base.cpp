/****************************************************************************
* Copyright (c) 2023, CEA
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
// File      : IJK_Thermal_base.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Thermal_base.h>
#include <Param.h>
#include <IJK_Navier_Stokes_tools.h>
#include <stat_counters.h>
#include <IJK_FT.h>
#include <DebogIJK.h>
#include <Corrige_flux_FT.h>
#include <OpConvDiscIJKQuickScalar.h>


Implemente_base_sans_constructeur( IJK_Thermal_base, "IJK_Thermal_base", Objet_U ) ;


/********************************************
 * Methods inherited from Objet_U
 ********************************************/

IJK_Thermal_base::IJK_Thermal_base()
{
  rang_ = 0; // default value, used as an index for the list of thermal sub-problems
  calulate_grad_T_=0;
  calculate_local_energy_=0;

  /*
   * Physical parameters
   */
  cp_liquid_=1.;
  cp_vapour_=0.;
  lambda_liquid_=1.;
  lambda_vapour_=0.;
  single_phase_=1.;
  fo_ = 1.; // fourier number

  /*
   * Initialisation of the problem
   */
  expression_T_init_="??";
  fichier_reprise_temperature_= "??";
  timestep_reprise_temperature_=1;

  conv_temperature_negligible_=0;
  diff_temp_negligible_=0;
  type_temperature_convection_op_ = 1;  // Default value: 1 - Quick

  type_T_source_="??";
  wall_flux_=0;
  lambda_variable_=0; //terme source variable
  /*
   * Ceci est une initialisation des derivees des temperatures moyenne de chaque phase
   * Il n'est peut-etre pas pertinent de les mettre ici
   */
  dTv_ = 0.;
  dTl_ = 1.;
  Tl_ = 0.;
  Tv_ = 1.;
  Tv0_ = 1.;  // Serait-ce plutot Tref (une temperature de reference pour reconstruire le champ dim??)
  kl_ = -100000000000000.;
  kv_ = -200000000000000.;
  T0v_ = 1.;
  T0l_ = 0.;
  expression_T_ana_="??";
}

Sortie& IJK_Thermal_base::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  os<< "  {\n"
    << "    boundary_conditions {"  << "\n";
  /*
   * Boundary conditions (Periodicity or wall)
   */
  Nom bctype_kmin, bctype_kmax, bckmin, bckmax, valeur_kmin, valeur_kmax;
  if( boundary_conditions_.get_bctype_k_max()==boundary_conditions_.Paroi_Temperature_imposee)
    {
      bctype_kmax="Paroi_Temperature_imposee";
      bckmax = "temperature_imposee_kmax";
      valeur_kmax = boundary_conditions_.get_temperature_kmax();
    }
  else if( boundary_conditions_.get_bctype_k_max()==boundary_conditions_.Paroi_Flux_impose)
    {
      bctype_kmax="Paroi_Flux_impose";
      bckmax = "flux_impose_kmax";
      valeur_kmax = boundary_conditions_.get_flux_kmax();
    }
  else if( boundary_conditions_.get_bctype_k_max()==boundary_conditions_.Perio)
    {
      bctype_kmax="Perio";
      bctype_kmin="Perio";
      bckmax = " ";
      bckmin = " ";
      valeur_kmax = " ";
    }

  if( boundary_conditions_.get_bctype_k_min()==boundary_conditions_.Paroi_Temperature_imposee)
    {
      bctype_kmin="Paroi_Temperature_imposee";
      bckmin = "temperature_imposee_kmin";
      valeur_kmin = boundary_conditions_.get_temperature_kmin();
    }
  else if( boundary_conditions_.get_bctype_k_min()==boundary_conditions_.Paroi_Flux_impose)
    {
      bctype_kmin="Paroi_Flux_impose";
      bckmin = "flux_impose_kmin";
      valeur_kmin = boundary_conditions_.get_flux_kmin();
    }
  else if( boundary_conditions_.get_bctype_k_min()==boundary_conditions_.Perio)
    {
      bctype_kmin="Perio";
      bctype_kmin="Perio";
      bckmin = " ";
      bckmin = " ";
      valeur_kmin = " ";
    }
  os<< "      bctype_kmin" << " " << bctype_kmin << " \n";
  os<< "      bctype_kmax" << " " << bctype_kmax << " \n";

  os<< "      " << bckmin << " " << valeur_kmin << " \n";
  os<< "      " << bckmax << " " << valeur_kmax << " \n";
  os<< "    } \n" ;

  /*
   * Physical parameters
   */
  os<< "    lambda_liquid " <<  lambda_liquid_ << "\n";
  os<< "    lambda_vapour " <<  lambda_vapour_ << "\n";
  os<< "    cp_liquid " <<  cp_liquid_ << "\n";
  os<< "    cp_vapour " <<  cp_vapour_ << "\n";

  /*
   * Source term
   */
  os<< "    type_T_source " << type_T_source_ << "\n";
  if (type_T_source_=="SWARM")
    {
      os<< "      kl_source " <<  kl_ << "\n";
      os<< "      kv_source " <<  kv_ << "\n";
      os<< "      T0l_source " <<  T0l_ << "\n";
      os<< "      T0v_source " <<  T0v_ << "\n";
    }

  if( wall_flux_)
    os << "    wall_flux \n";

  /*
   * Resume calculation
   */
  os<< "    fichier_reprise_temperature" << " " << fichier_reprise_temperature_  << "\n";
  os<< "    timestep_reprise_temperature" << " " << timestep_reprise_temperature_ << "\n";

  /*
   * Analytical expression of temperature at t_initial
   */
  if ( expression_T_ana_!="??")
    os<< "    expression_T_ana" <<  " " << expression_T_ana_ << "\n";

  /*
   * Neglect an operator
   */
  if ( conv_temperature_negligible_)
    os<< "    conv_temperature_negligible \n ";

  if ( diff_temp_negligible_)
    os<< "    diff_temp_negligible \n";

  os<< "  }\n";
  return os;
}

// XD thermique listobj thermique -1 thermique_bloc 1 to add energy equation resolution if needed
// XD thermique_bloc interprete nul 1 not_set
Entree& IJK_Thermal_base::readOn( Entree& is )
{
  /*
   * Parse the datafile
   */
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades(is);

  Cout << "IJK_Thermal_base::readOn : Parameters summary. " << finl;
  printOn(Cout);
  return is;
}

void IJK_Thermal_base::set_param(Param& param)
{
  param.ajouter("fo", &fo_); // XD_ADD_P floattant not_set
  param.ajouter("cp_liquid", &cp_liquid_, Param::REQUIRED); // XD_ADD_P floattant Liquid specific heat at constant pressure
  param.ajouter("lambda_liquid", &lambda_liquid_, Param::REQUIRED); // XD_ADD_P floattant Liquid thermal conductivity

  param.ajouter("expression_T_init", &expression_T_init_); // XD_ADD_P chaine Expression of initial temperature (parser of x,y,z)

  param.ajouter("boundary_conditions", &boundary_conditions_, Param::REQUIRED); // XD_ADD_P bloc_lecture boundary conditions
  param.ajouter("type_T_source", &type_T_source_); // XD_ADD_P chaine(into=["dabiri","patch_dabiri","unweighted_dabiri"]) source term
  param.ajouter("expression_source_temperature", &expression_source_temperature_); // XD_ADD_P chaine source terms
  param.ajouter_flag("lambda_variable", &lambda_variable_);
  param.ajouter_flag("wall_flux", &wall_flux_); // XD_ADD_P rien not_set
  param.ajouter("kl_source", &kl_);
  param.ajouter("kv_source", &kv_);
  param.ajouter("T0l_source", &T0l_);
  param.ajouter("T0v_source", &T0v_);

  param.ajouter("fichier_reprise_temperature", &fichier_reprise_temperature_);
  param.ajouter("timestep_reprise_temperature", &timestep_reprise_temperature_);

  param.ajouter_flag("conv_temperature_negligible", &conv_temperature_negligible_); // XD_ADD_P rien neglect temperature convection
  param.ajouter("type_temperature_convection_op", &type_temperature_convection_op_); // XD_ADD_P chaine(into=["Amont","Quick","Centre2","Centre4"]) convection operator
  param.dictionnaire("Quick",1);
  param.dictionnaire("Centre2",2);
  param.ajouter_flag("diff_temp_negligible", &diff_temp_negligible_); // XD_ADD_P rien neglect temperature diffusion

  param.ajouter("expression_T_ana", &expression_T_ana_); // XD_ADD_P chaine Analytical expression T=f(x,y,z,t) for post-processing only

//  param.ajouter_flag("single_phase", &single_phase_);
  param.ajouter("calculate_local_energy", &calculate_local_energy_);
}


/********************************************
 * Public methods
 ********************************************/

int IJK_Thermal_base::initialize(const IJK_Splitting& splitting, const int idx)
{
  Cout << que_suis_je() << "::initialize()" << finl;
  rang_ = idx;
  int nalloc = 0;

  /*
   * Diffusion operator
   */
  if (single_phase_)
    {
      temperature_diffusion_op_.typer("OpDiffUniformIJKScalar_double");
      temperature_diffusion_op_.set_uniform_lambda(uniform_lambda_);
    }
  else
    {
      temperature_diffusion_op_.typer("OpDiffIJKScalar_double");
    }

  /*
   * Convection operator
   */
  switch(type_temperature_convection_op_)
    {
    case 1:
      temperature_convection_op_.typer("OpConvIJKQuickScalar_double");
      break;
    case 2:
      temperature_convection_op_.typer("OpConvCentre2IJKScalar_double");
      break;
    default:
      Cerr << "Undefined operator for the convection of the temperature. " << finl;
      Process::exit();
    }

  /*
   * Initialise the operators
   */
  temperature_convection_op_.initialize(splitting);
  temperature_diffusion_op_.initialize(splitting);

  /*
   * Fields
   */
  temperature_.allocate(splitting, IJK_Splitting::ELEM, 2);
  d_temperature_.allocate(splitting, IJK_Splitting::ELEM, 2);
  nalloc += 2;

  if (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("RHO_CP"))
    {
      rho_cp_.allocate(splitting, IJK_Splitting::ELEM, 2);
      nalloc += 1;
    }
  if (calculate_local_energy_)
    {
      rho_cp_T_.allocate(splitting, IJK_Splitting::ELEM, 2);
      nalloc += 1;
    }

  //	lambda_.allocate(splitting, IJK_Splitting::ELEM, 1);
  //  unit_field_.allocate(splitting, IJK_Splitting::ELEM, 2);
  //  unit_field_.data() = 1.;
  //  unit_field_.echange_espace_virtuel(unit_field_.ghost());

  /*
   * Storage for temperature gradient post-processing or method
   */
  if ((ref_ijk_ft_.non_nul()) && (!ref_ijk_ft_->get_disable_diphasique()))
    {
      Cout << "Allocating fields temperature_ft_ and storage" << finl;
      allocate_cell_vector(storage_, ref_ijk_ft_->get_splitting_ft(), 1);
      nalloc += 3;
    }

  /*
   * Dimensionless temperature field (thermostat)
   */
  if ((wall_flux_) || liste_post_instantanes_.contient_("SOURCE_TEMPERATURE")
      || liste_post_instantanes_.contient_("TEMPERATURE_PHYSIQUE_T")
      || liste_post_instantanes_.contient_("TEMPERATURE_ADIMENSIONNELLE_THETA")
      || (type_T_source_ != "??"))
    {
      Cout << "Allocating field for the thermal source term & co. " << finl;
      source_temperature_.allocate(splitting, IJK_Splitting::ELEM, 1);
      source_temperature_v_.allocate(splitting, IJK_Splitting::ELEM, 1);
      source_temperature_l_.allocate(splitting, IJK_Splitting::ELEM, 1);
      d_source_Tv_.allocate(splitting, IJK_Splitting::ELEM, 1);
      d_source_Tl_.allocate(splitting, IJK_Splitting::ELEM, 1);
      //d_source_T_.allocate(splitting, IJK_Splitting::ELEM, 1);
      temperature_physique_T_.allocate(splitting, IJK_Splitting::ELEM, 2);
      temperature_adimensionnelle_theta_.allocate(splitting, IJK_Splitting::ELEM, 2);
      nalloc += 7;
      // par defaut s'il n'y a pas de source renseignee, on utilise la source de Dabiri/Kawamura
      // cela veut dire que dans le cas des SWARMS il faut imperativement renseigner le nom de
      // la source
      if (type_T_source_ == "??")
        {
          Cerr << "Attention on demande des post-traitement sans avoir renseigner type_T_source" << finl;
          throw "Erreur post et type_T_source";
          // type_T_source_ = "dabiri";
        }
    }
  if (liste_post_instantanes_.contient_("TEMPERATURE_ADIM_BULLES"))
    {
      temperature_adim_bulles_.allocate(splitting, IJK_Splitting::ELEM, 2);
      nalloc += 1;
    }

  /*
   * RK3 sub-steps
   * Check that pointer is not null:
   */
  if (ref_ijk_ft_.non_nul() && ref_ijk_ft_->get_time_scheme()== ref_ijk_ft_->RK3_FT)
    {
      RK3_F_temperature_.allocate(splitting, IJK_Splitting::ELEM, 0);
      nalloc +=1;
    }

  /*
   * Resume the calculation
   */
  if (fichier_reprise_temperature_ == "??")   // si on ne fait pas une reprise on initialise V
    {
      if (expression_T_init_ != "??")
        {
          Cout << "Temperature initialization from expression \nTini = " << expression_T_init_ << finl;
          set_field_data(temperature_, expression_T_init_, ref_ijk_ft_->itfce().I(), 0.);
        }
      else
        {
          Cerr << "Please provide initial conditions for temperature, either by an expression or a field for restart. "
               << "You should consider using either fichier_reprise_temperature or expression_T_init keywords. " << finl;
          Process::exit();
        }
    }
  else
    {
      Cout << "Reading initial temperature field T"<< rang_<<" from file " << fichier_reprise_temperature_ << " timestep= " << timestep_reprise_temperature_ << finl;
      const Nom& geom_name = splitting.get_grid_geometry().le_nom();
      lire_dans_lata(fichier_reprise_temperature_, timestep_reprise_temperature_, geom_name, Nom("TEMPERATURE_")+Nom(idx),
                     temperature_); // fonction qui lit un champ a partir d'un lata .
      temperature_.echange_espace_virtuel(temperature_.ghost()); // It is essential to fill the EV because the first call to convection needs them.
    }

  /*
   * List of post-processed data
   */
  Cerr << " Initializing thermal fields dependant on the post-pro list : " ;
  if (liste_post_instantanes_.size())
    Cerr << liste_post_instantanes_;
  else
    Cerr << "empty";
  Cerr << finl;
  if (liste_post_instantanes_.size() && (liste_post_instantanes_.contient_("TEMPERATURE_ANA") || liste_post_instantanes_.contient_("ECART_T_ANA")))
    {
      temperature_ana_.allocate(splitting, IJK_Splitting::ELEM, 1);
      ecart_t_ana_.allocate(splitting, IJK_Splitting::ELEM, 1);
      nalloc +=2;
    }

  /* if ((liste_post_instantanes_.contient_("SOURCE_TEMPERATURE_ANA")) || (liste_post_instantanes_.contient_("ECART_SOURCE_TEMPERATURE_ANA")) )
     {
       source_temperature_ana_.allocate(splitting, IJK_Splitting::ELEM, 1);
       ecart_source_t_ana_.allocate(splitting, IJK_Splitting::ELEM, 1);
       nalloc +=2;
     } */
  if ((liste_post_instantanes_.size() && liste_post_instantanes_.contient_("GRAD_T"))
      || (ref_ijk_ft_.non_nul() && ref_ijk_ft_->t_debut_statistiques() <  1.e10 ))
    {
      allocate_velocity(grad_T_, splitting, 1);
      nalloc +=3 ;
    }

  Cout << "End of " << que_suis_je() << "::initialize()" << finl;
  return nalloc;
}

void IJK_Thermal_base::update_thermal_properties()
{
  if (single_phase_)
    {
      cp_vapour_=0.;
    }
  const double ene_ini = compute_global_energy();

  const IJK_Field_double& indic = ref_ijk_ft_->itfce().I();
  // Nombre de mailles du domaine NS :
  const int nx = indic.ni();
  const int ny = indic.nj();
  const int nz = indic.nk();
  const double rho_l = ref_ijk_ft_->get_rho_l();
  const double rho_v = ref_ijk_ft_->get_rho_v();
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        {
          double chi_l = indic(i,j,k);
          if (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("RHO_CP"))
            {
              rho_cp_(i,j,k) = rho_l*cp_liquid_*chi_l + rho_v*cp_vapour_*(1-chi_l);
            }
          if (calculate_local_energy_)
            {
              rho_cp_T_(i,j,k) = (rho_l*cp_liquid_*chi_l + rho_v*cp_vapour_*(1-chi_l))*temperature_(i,j,k);
            }
        }
  if (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("RHO_CP"))
    {
      rho_cp_.echange_espace_virtuel(rho_cp_.ghost());
    }
  if (calculate_local_energy_)
    {
      rho_cp_T_.echange_espace_virtuel(rho_cp_T_.ghost());
    }
  // Semble un endroit approprie pour calculer la variation d'energie due au transport de l'interface:
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T"<<rang_<<"-2-TransportIndic] time t=" << ref_ijk_ft_->get_current_time()
       << " " << ene_ini
       << " " << ene_post
       << " delta=" << ene_post-ene_ini << " [W.m-3]." << finl;
}

// Methode de calcul du pas de temps max base sur Fo pour l'equation de thermique.
// CFL value is not computed as it is the same as for the velocity equation.
// The calculation should be stable if Fo <= 1.0 (thanks to the 0.5 in the formula below).
double IJK_Thermal_base::compute_timestep(const double timestep,
                                          const double dxmin) const
{
  double alpha_max;
  double rho_l = ref_ijk_ft_->get_rho_l();
  double rho_v= ref_ijk_ft_->get_rho_v();
  if (single_phase_)
    {
      alpha_max = lambda_liquid_ / (rho_l * cp_liquid_);
    }
  else
    {
      alpha_max = std::max(lambda_liquid_ / (rho_l * cp_liquid_), lambda_vapour_ / (rho_v * cp_vapour_));
    }
  double dt_fo  = dxmin*dxmin/(alpha_max + 1.e-20) * fo_ * (1./6.); // Attention 0.125 vient du 3D. (1/6 au lieu de 1/8)
  if (diff_temp_negligible_) dt_fo = 1.e20;
  return dt_fo;
}

void IJK_Thermal_base::associer(const IJK_FT_double& ijk_ft)
{
  ref_ijk_ft_ = ijk_ft;
  liste_post_instantanes_ = ijk_ft.get_post().get_liste_post_instantanes();
}

void IJK_Thermal_base::euler_time_step(const double timestep)
{
  calculer_dT(ref_ijk_ft_->get_velocity());
  // Update the temperature :
  const int kmax = temperature_.nk();
  const double ene_ini = compute_global_energy();
  for (int k = 0; k < kmax; k++)
    {
      ref_ijk_ft_->euler_explicit_update(d_temperature_, temperature_, k);
    }
  temperature_.echange_espace_virtuel(temperature_.ghost());
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T"<<rang_<<"] time t=" << ref_ijk_ft_->get_current_time()
       << " " << ene_ini
       << " " << ene_post << " [W.m-3]." << finl;
  source_callback();
}

void IJK_Thermal_base::rk3_sub_step(const int rk_step, const double total_timestep,
                                    const double time)
{
  calculer_dT(ref_ijk_ft_->get_velocity());
  // Update the temperature :
  const int kmax = temperature_.nk();
  const double ene_ini = compute_global_energy();
  for (int k = 0; k < kmax; k++)
    {
      runge_kutta3_update(d_temperature_, RK3_F_temperature_, temperature_, rk_step, k, total_timestep);
    }
  temperature_.echange_espace_virtuel(temperature_.ghost());
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T"<<rang_<<"] time t=" << ref_ijk_ft_->get_current_time()
       << " " << ene_ini
       << " " << ene_post << " [W.m-3]. [step"<< rk_step << "]" << finl;
  source_callback();
}

void IJK_Thermal_base::sauvegarder_temperature(Nom& lata_name, int idx)
{
  fichier_reprise_temperature_ = lata_name;
  timestep_reprise_temperature_ = 1;
  dumplata_scalar(lata_name,Nom("TEMPERATURE_")+Nom(idx) , temperature_, 0 /*we store a 0 */);
}

/********************************************
 * Protected methods
 ********************************************/


// Mettre rk_step = -1 si schema temps different de rk3.
void IJK_Thermal_base::calculer_dT(const FixedVector<IJK_Field_double, 3>& velocity)
{
  const double current_time = ref_ijk_ft_->get_current_time();
  const double ene_ini = compute_global_energy(d_temperature_);

  /*
   * Correct the temperature field using either the ghost-fluid
   * approach or the laminar sub-resolution approach (and zero values for debug)
   */
  correct_temperature_for_eulerian_fluxes();

  compute_temperature_convection(velocity);

  const double ene_postConv = compute_global_energy(d_temperature_);
  add_temperature_diffusion();
  const double ene_postDiffu = compute_global_energy(d_temperature_);
  add_temperature_source();
  const double ene_postSource = compute_global_energy(d_temperature_);

  Cerr << "[Energy-Budget-T"<<rang_<<"-1-TimeResolution] time t=" << current_time
       << " " << ene_ini
       << " " << ene_postConv
       << " " << ene_postDiffu
       << " " << ene_postSource
       << " delta=" << ene_postSource-ene_ini << " [W.m-3]." << finl;

  const IJK_Field_double& T = temperature_;
  double Tmax = -1.e20;
  double Tmin = 1.e20;
  const int ni = T.ni();
  const int nj = T.nj();
  const int nk = T.nk();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          Tmax = std::max(Tmax, T(i,j,k));
          Tmin = std::min(Tmin, T(i,j,k));
        }
  Tmax = Process::mp_max(Tmax);
  Tmin = Process::mp_min(Tmin);
  Cerr <<"[Temperature-MinMax-" << rang_ <<"] t/Tmin/Tmax " << current_time << " "
       << Tmin << " " << Tmax
       << finl;
  return;
}

// Convect temperature field by velocity.
// The output is stored in d_temperature_ (it is a volume integral over the CV)
void IJK_Thermal_base::compute_temperature_convection(const FixedVector<IJK_Field_double, 3>& velocity)
{
  static Stat_Counter_Id cnt_conv_temp = statistiques().new_counter(1, "FT convection rho");
  statistiques().begin_count(cnt_conv_temp);
  if (conv_temperature_negligible_)
    {
      d_temperature_.data()=0;
    }
  else
    {
      temperature_convection_op_.calculer(temperature_, velocity[0], velocity[1], velocity[2], d_temperature_);
      const int ni = d_temperature_.ni();
      const int nj = d_temperature_.nj();
      const int nk = d_temperature_.nk();
      const IJK_Grid_Geometry& geom = d_temperature_.get_splitting().get_grid_geometry();
      const double dx = geom.get_constant_delta(DIRECTION_I);
      const double dy = geom.get_constant_delta(DIRECTION_J);
      const double dz = geom.get_constant_delta(DIRECTION_K);
      const double vol = dx*dy*dz;
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              d_temperature_(i,j,k) /= vol ;
            }
    }
  statistiques().end_count(cnt_conv_temp);
  DebogIJK::verifier("op_conv(rho)", d_temperature_);
  return;
}

void IJK_Thermal_base::add_temperature_diffusion()
{
  if (boundary_conditions_.get_bctype_k_min() == Boundary_Conditions_Thermique::Paroi_Temperature_imposee)
    {
      const double T_paroi_impose_kmin = boundary_conditions_.get_temperature_kmin();
      double lambda_de_t_paroi_kmin = lambda_liquid_;
      // calculer ici div(lambda*grad(T))*volume)
      calculer_flux_thermique_bord(temperature_, lambda_de_t_paroi_kmin,
                                   T_paroi_impose_kmin, boundary_flux_kmin_, 0 /* boundary kmin */);
    }
  else if (boundary_conditions_.get_bctype_k_min() == Boundary_Conditions_Thermique::Paroi_Flux_impose)
    {
      const double flux_paroi_impose_kmin = boundary_conditions_.get_flux_kmin();
      imposer_flux_thermique_bord(temperature_,
                                  flux_paroi_impose_kmin, boundary_flux_kmin_, 0 /* boundary kmin */);
    }
  else
    {
      // Perio... not needed. The flux on the "boundary" (in that case , it's not truely a boundary) will be computed as inside...
    }
  if (boundary_conditions_.get_bctype_k_max() == Boundary_Conditions_Thermique::Paroi_Temperature_imposee)
    {
      const double T_paroi_impose_kmax = boundary_conditions_.get_temperature_kmax();
      double lambda_de_t_paroi_kmax = lambda_liquid_;
      // calculer ici div(lambda*grad(T))*volume)
      calculer_flux_thermique_bord(temperature_, lambda_de_t_paroi_kmax,
                                   T_paroi_impose_kmax, boundary_flux_kmax_, 1 /* boundary kmax */);
    }
  else if (boundary_conditions_.get_bctype_k_max() == Boundary_Conditions_Thermique::Paroi_Flux_impose)
    {
      const double flux_paroi_impose_kmax = boundary_conditions_.get_flux_kmax();
      imposer_flux_thermique_bord(temperature_,
                                  flux_paroi_impose_kmax, boundary_flux_kmax_, 1 /* boundary kmax */);
      // Cerr << "not coded yet" << finl;
      // Process::exit();
    }
  else
    {
      // Perio... not needed. The flux on the "boundary" (in that case , it's not truely a boundary) will be computed as inside...
    }
  temperature_.echange_espace_virtuel(temperature_.ghost());
  DebogIJK::verifier("temp", temperature_);

  // Performance counters:
  static Stat_Counter_Id cnt_diff_temp = statistiques().new_counter(1, "FT diffusion temperature");
  statistiques().begin_count(cnt_diff_temp);

  if (diff_temp_negligible_) // si diffusion negligeable
    {
      div_coeff_grad_T_volume_.data()=0;
    }
  else
    {
      /*
       * Correct directly the fluxes in the operator
       */
      temperature_diffusion_op_.calculer(temperature_,
                                         div_coeff_grad_T_volume_,
                                         boundary_flux_kmin_, boundary_flux_kmax_);
    }

  /*
   * Correct the diffusive fluxes here or in the operator ?
   */
  //	correct_temperature_increment_diffusion();

  compute_diffusion_increment();

  statistiques().end_count(cnt_diff_temp);
  DebogIJK::verifier("div_coeff_grad_T_volume_", div_coeff_grad_T_volume_);

  if ((liste_post_instantanes_.contient_("GRAD_T") || (calulate_grad_T_))
      || (ref_ijk_ft_->t_debut_statistiques() <  1.e10 ))
    {
      calculer_gradient_temperature(temperature_, grad_T_);
    }
}


//////////////////////////////////////////

void IJK_Thermal_base::add_temperature_source() { ; }

//void IJK_Thermal_base::add_temperature_source()
//// DONE: ajouter un mot clef pour changer le type de source (constante, ponderee par la vitesse, ...)
//{
//  static Stat_Counter_Id cnt_source_temp = statistiques().new_counter(1, "FT source temperature");
//  statistiques().begin_count(cnt_source_temp);
//
//  // dans le cas ou les flux entrants et sortants sont identiques :
//  // DONE: changer cette condition non adaptee
//  if (type_T_source_!="??")
//    {
//      const IJK_Field_double&  vx = ref_ijk_ft_->get_velocity()[DIRECTION_I];
//      double rho_cp_u_moy = calculer_rho_cp_u_moyen(vx, cp_,ref_ijk_ft_->rho_field_);
//      // double rho_cp_moy = calculer_rho_cp_moyen(cp_, ref_ijk_ft_->rho_field_);
//      const IJK_Splitting& splitting = temperature_.get_splitting();
//      const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
//      const double dx =geom.get_constant_delta(DIRECTION_I);
//      const double dy = geom.get_constant_delta(DIRECTION_J);
//      const double dz = geom.get_constant_delta(DIRECTION_K);
//      const double vol = dx*dy*dz;
//      const double lz = geom.get_domain_length(DIRECTION_K) ;
//      const double h =lz/2.;
//      // double volume = 1.;
//      // for (int i = 0; i < 3; i++)
//      //   volume *= splitting.get_grid_geometry().get_constant_delta(i);
//      //IJK_Field_double& dT = d_temperature_(i,j,k);
//      const int nk = d_temperature_.nk();
//      const int ni = d_temperature_.ni();
//      const int nj = d_temperature_.nj();
//      // TODO: faire une methode calculer_rho_cp
//      //debut if source = ponderee par la vitesse
//      if (type_T_source_=="dabiri")
//        {
//
//          const double wall_flux = boundary_conditions_.get_flux_kmax();
//          const double qw = wall_flux;
//          const double dTm = 2*qw/rho_cp_u_moy ; // TODO: faux, il manque un facteur 1/h pour homogeneite
//
//          for (int k = 0; k < nk; k++)
//            for (int j = 0; j < nj; j++)
//              for (int i = 0; i < ni; i++)
//                {
//                  const double rho = ref_ijk_ft_->rho_field_(i,j,k);
//                  const double cp = cp_(i,j,k);
//                  const double u = (vx(i,j,k) +vx(i+1,j,k))/2;
//                  const double rho_cp_u = rho*cp*u;
//                  // const double lambda = lambda_(i,j,k);
//                  //const double f = rho_cp_u; // - div (lambda);
//                  if(lambda_variable_)
//                    {
//                      const double div_lambda = (lambda_(i+1,j,k)-lambda_(i-1,j,k))/(2*dx);
//                      source_temperature_(i,j,k) = (rho_cp_u-div_lambda)*dTm;
//                    }
//                  else
//                    {
//                      source_temperature_(i,j,k) = rho_cp_u*dTm;
//                    }
//                  // const double lambda = lambda_(i,j,k);
//                  const double Sc = qw/h ;//qw/lambda*h;
//                  const double source = (source_temperature_(i,j,k)-Sc)*vol;//-Sc)*volume; TODO: faux, que vient faire Sc en plus de source ?
//
//                  d_temperature_(i,j,k) += source;
//
//                }
//          //
//          calculer_temperature_physique_T(vx, dTm);
//          calculer_temperature_adimensionnelle_theta(vx, wall_flux);
//          calculer_Nusselt(vx);
//          return;
//        }
//      // TODO: remplacer dabiri par patch_dabiri apres verif
//      else if (type_T_source_=="patch_dabiri")
//        {
//
//          Cerr << "Type de source : patch_dabiri" << finl;
//          const double wall_flux = boundary_conditions_.get_flux_kmax();
//          const double qw = wall_flux;
//          const double dTm = -2*qw/(2*h*rho_cp_u_moy) ;
//
//          for (int k = 0; k < nk; k++)
//            for (int j = 0; j < nj; j++)
//              for (int i = 0; i < ni; i++)
//                {
//                  //const double chi = ref_ijk_ft_->itfce().I(i,j,k);
//                  //const double rho_cp_l = ref_ijk_ft_->get_rho_l() * cp_liquid_;
//                  //const double rho_cp_v = ref_ijk_ft_->get_rho_v() * cp_vapour_;
//                  const double u = (vx(i,j,k) +vx(i+1,j,k))/2;
//                  //const double rho_cp_u = (chi*rho_cp_l + (1-chi)*rho_cp_v)*u;
//                  if(lambda_variable_)
//                    {
//                      const double div_lambda = (lambda_(i+1,j,k)-lambda_(i-1,j,k))/(2*dx);
//                      source_temperature_(i,j,k) = (u-div_lambda)*dTm;
//                    }
//                  else
//                    {
//                      source_temperature_(i,j,k) = u*dTm;
//                    }
//                  const double source = source_temperature_(i,j,k);
//
//                  d_temperature_(i,j,k) += source;
//
//                }
//          //
//          calculer_temperature_physique_T(vx, dTm);
//          calculer_temperature_adimensionnelle_theta(vx, wall_flux);
//          calculer_Nusselt(vx);
//          return;
//        }
//      else if (type_T_source_=="unweighted_dabiri")
//        {
//          // Sans la ponderation de S par u (S=rhocpu.../<rhocpu>), cela permet d'avoir une source uniforme!
//          // Mais ce n'est plus le meme changement de variable, c'est quoi alors?
//          // TODO: Que faire de rho_cp en diphasique? Moyen ou local?
//          Cerr << "Type de source : unweighted_dabiri" << finl;
//          const double wall_flux = boundary_conditions_.get_flux_kmax();
//          const double qw = wall_flux;
//          const double liquid_fraction = calculer_v_moyen(ref_ijk_ft_->itfce().I());
//          const double rho_cp_l = ref_ijk_ft_->get_rho_l() * cp_liquid_;
//          const double rho_cp_v = ref_ijk_ft_->get_rho_v() * cp_vapour_;
//          const double rhocp_moy =  rho_cp_l*liquid_fraction + rho_cp_v*(1-liquid_fraction);
//          const double dTm = -2*qw/(2*h*rhocp_moy) ;
//          for (int k = 0; k < nk; k++)
//            for (int j = 0; j < nj; j++)
//              for (int i = 0; i < ni; i++)
//                {
//                  if(lambda_variable_)
//                    {
//                      Cerr << "Veut-on vraiment calculer une partie en lambda spatialement variable??" <<finl;
//                      Cerr << "Exit at IJK_Thermal_base::add_temperature_source" << finl;
//                      Process::exit();
//                      // const double div_lambda = (lambda_(i+1,j,k)-lambda_(i-1,j,k))/(2*dx);
//                      // source_temperature_(i,j,k) = (u-div_lambda)*dTm/<u>;
//                    }
//                  else
//                    {
//                      source_temperature_(i,j,k) = dTm;
//                    }
//                  const double source = source_temperature_(i,j,k);
//                  d_temperature_(i,j,k) += source;
//                }
//          //
//          calculer_temperature_physique_T(vx, dTm);
//          calculer_temperature_adimensionnelle_theta(vx, wall_flux);
//          calculer_Nusselt(vx);
//          return;
//        }
//      //debut source = SWARM
//      else if (type_T_source_=="SWARM")
//        {
//          //DONE: idem
//          double Tv=0;
//          double Tl=0;
//          double Vv = 0;
//          double Vl = 0;
//
//          //calcul de Tv et Tl
//          for (int k = 0; k < nk; k++)
//            for (int j = 0; j < nj; j++)
//              for (int i = 0; i < ni; i++)
//                {
//                  const double chi = ref_ijk_ft_->itfce().I(i, j, k);
//                  const double T = temperature_(i, j, k);
//                  Tv += T*(1.-chi);
//                  Vv += (1.-chi);
//                  Tl += T*chi;
//                  Vl += chi;
//                }
//          Tv = Process::mp_sum(Tv);
//          Tl = Process::mp_sum(Tl);
//          Vv = Process::mp_sum(Vv);
//          Vl = Process::mp_sum(Vl);
//          Tv /= Vv;
//          Tl /= Vl;
//          Cerr << "AY-test_source : " <<  Tv << finl;
//
//          //calcul de dT a partir de l'expression suivante : dT = k(T - Tm)/(rho*cp)
//          const double kv = kv_;
//          const double kl = kl_;
//          const double T0v = T0v_;
//          const double T0l = T0l_;
//          const double rho_cp_l = ref_ijk_ft_->get_rho_l() * cp_liquid_;
//          const double rho_cp_v = ref_ijk_ft_->get_rho_v() * cp_vapour_;
//          for (int k = 0; k < nk; k++)
//            for (int j = 0; j < nj; j++)
//              for (int i = 0; i < ni; i++)
//                {
//                  d_source_Tv_(i,j,k) = kv/rho_cp_v * (Tv - T0v);
//                  d_source_Tl_(i,j,k) = kl/rho_cp_l * (Tl - T0l);
//                }
//          Cerr << "AY-test_source1 : " <<  Tv << finl;
//          // TODO: remplacer euler_explicit_update par l'utilisation de timestep_ et utiliser d_source_Tv_ comme une constante
//          for (int k = 0; k < nk; k++)
//            {
//              ref_ijk_ft_->euler_explicit_update(d_source_Tv_, source_temperature_v_, k);
//              ref_ijk_ft_->euler_explicit_update(d_source_Tl_, source_temperature_l_, k);
//            }
//          Cerr << "AY-test_source2 : " <<  Tv << finl;
//          for (int k = 0; k < nk; k++)
//            for (int j = 0; j < nj; j++)
//              for (int i = 0; i < ni; i++)
//                {
//                  // chi vaut 1. dans le liquide et 0 dans les bulles
//                  const double chi = ref_ijk_ft_->itfce().I(i, j, k);
//                  source_temperature_(i,j,k) = (1.-chi)*source_temperature_v_(i,j,k) + chi*source_temperature_l_(i,j,k);
//                }
//          Cerr << "AY-test_source3 : " <<  Tv << finl;
//          for (int k = 0; k < nk; k++)
//            for (int j = 0; j < nj; j++)
//              for (int i = 0; i < ni; i++)
//                {
//                  d_temperature_(i,j,k) += source_temperature_(i,j,k)*vol;
//                }
//          Cerr << "AY-test_source 111 : " <<  source_temperature_(1,1,1) << finl;
//          Cerr << "AY-test_source vol : " <<  vol << finl;
//          Cerr << "source_temp " << " " << d_temperature_(1,1,1) << finl;
//          const double current_time = ref_ijk_ft_->get_current_time();
//
//          // GB : Ma comprehension est que ce ne sont pas des champs, mais un scalaire unique
//          const double Sl = source_temperature_l_(0,0,0);
//          const double Sv = source_temperature_v_(0,0,0);
//          const double dSl = d_source_Tl_(0,0,0);
//          const double dSv = d_source_Tv_(0,0,0);
//          Cerr <<"[ThermalInfo-" << rang_ <<"] t/Tl/Tv/Sl/Sv/dSldt/dSvdt " << current_time << " "
//               << Tl << " " << Tv << " "
//               << Sl << " " << Sv << " "
//               << dSl << " " << dSv
//               <<finl;
//          // si on utilise ca c'est pour remplir dans tous les cas d'utilisation d'une source le champs temperature_physique_T
//          calculer_temperature_physique_T_dummy();
//
//          return;
//        }
//    }
//  // dans ce cas le ce ne sont pas des flux thermiques identiques
//  else
//    {
//      Cerr << "no_source_for_temperature" << finl;
//
//      // calculer_source_temperature_ana();
//      //     Cerr << "source_temp" << " " << source_temperature_(24,24,0) << " " << source_temperature_(24,24,12) << " " <<  source_temperature_(24,24,23)  << finl;
//      //     Cerr << " debit_moy " << rho_cp_u_moy << " vx_max " << (vx(0,0,12)+vx(1,0,12))/2
//      //        << " vx_wall_0 " << (vx(0,0,0)+vx(1,0,0))/2 << " vx_wall_lz " << (vx(0,0,23)+vx(1,0,23))/2 << finl;
//      return;
//    }
//  //AYM Ce bloc n'est jamais atteint non ?
//  statistiques().end_count(cnt_source_temp);
//  Cerr << "source_temp" << " " <<  source_temperature_(1,1,1) << finl;
//  if (Process::je_suis_maitre())
//    {
//      const int nk = d_temperature_.nk();
//      Cerr << "source_temp mil" << " " <<  source_temperature_(1,1,int(nk/2)) << finl;
//    }
//}

void IJK_Thermal_base::source_callback()
{
  if (liste_post_instantanes_.contient_("TEMPERATURE_ADIM_BULLES"))
    {
      calculer_temperature_adim_bulles();
    }
}


/*
 * Aymeric: Renommer pour expliciter qu'il s'agit de la transformation
 * inverse de Kawamura avec le gradient de temperature moyenne
 */
void IJK_Thermal_base::calculer_temperature_physique_T(const IJK_Field_double&  vx, const double dTm)
{
  if (wall_flux_)
    {
      const IJK_Splitting& splitting = temperature_.get_splitting();

      //const IJK_Splitting::Localisation loc = field.get_localisation();
      const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
      double dx =geom.get_constant_delta(DIRECTION_I);
      double origin_x = geom.get_origin(DIRECTION_I) + (dx * 0.5) ;
      const int offset_i = splitting.get_offset_local(DIRECTION_I);

      const int nk = d_temperature_.nk();
      const int ni = d_temperature_.ni();
      const int nj = d_temperature_.nj();

      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  const double x = (i+ offset_i ) *dx + origin_x;    //MR A GB: ne doit-on pas soustraire par  -dx*0.5???
                  temperature_physique_T_(i,j,k) = (x*dTm)-temperature_(i,j,k);
                }
            }
        }
      temperature_physique_T_.echange_espace_virtuel(temperature_physique_T_.ghost());
      DebogIJK::verifier("temperature_physique_T", temperature_physique_T_);
      return;
    }
  else
    {
      Cerr << "no_source_for_temperature" << finl;
      return;
    }
}

void IJK_Thermal_base::calculer_temperature_adim_bulles()
{
  const int nk = temperature_.nk();
  const int ni = temperature_.ni();
  const int nj = temperature_.nj();

  // Calcul de moy1 = moy(chi*T+)/moy(chi) et moy2 = moy((1-chi)*T+) / moy(1-chi)
  double Tl = 0.;
  double Tv = 0.;
  double Vl = 0.;
  double Vv = 0.;
  const IJK_Field_double& T = temperature_;
  const IJK_Field_double& chi = ref_ijk_ft_->itfce().I(); // rappel : chi vaut 1. dans le liquide et 0 dans la vapeur

  // assuming uniform mesh : vol_cell=cste.
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          Tl += T(i,j,k)*chi(i, j, k);
          Tv += T(i,j,k)*(1.-chi(i, j, k));
          Vl += chi(i, j, k);
          Vv += (1.-chi(i, j, k));
        }
  Tl = Process::mp_sum(Tl);
  Tv = Process::mp_sum(Tv);
  Vv = Process::mp_sum(Vv);
  Vl = Process::mp_sum(Vl);
  Tl /= Vl;
  Tv /= Vv;

  // Calcul de Tl et Tv :

  // const double Tl = Tv0_ / (1 - Tl) * (1 - Tl / (1 + Tv)) / (1 + Tl * Tv / ((1 - Tl) * (1 + Tv)));
  // const double Tv = Tv0_ / (1 + Tv) * (1 + Tv / (1 - Tl)) / (1 + Tv * Tl / ((1 + Tv) * (1 - Tl)));
  const int ntot = temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_I)
                   *temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_J)
                   *temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const double time = ref_ijk_ft_->get_current_time();
  Cerr << "Tl_test : time= "<< time << " alpha_l= " << Vl/ntot<< "  TI=" <<Tl*Vl/ntot<< " Tl="<< Tl << finl;
  Cerr << "Tv_test : time= "<< time << " alpha_v= " << Vv/ntot<< " TIv=" <<Tv*Vv/ntot<< " Tv="<< Tv << finl;

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              temperature_adim_bulles_(i,j,k) = (temperature_(i,j,k) - Tv0_)/(Tl - Tv);
            }
        }
    }
//  temperature_adim_bulles_.echange_espace_virtuel(temperature_adim_bulles_.ghost());
//
//
//  TODO: ce qui suit ne devrait pas etre la, mais je le met ici temporairement avant de trouver une meilleure solution

  double E_tot = 0.;
  double E_liq_pure = 0., E_liq = 0;
  double E_vap_pure = 0., E_vap = 0;
  double E_mixt = 0.;
  calculer_energies(E_liq_pure, E_liq,
                    E_vap_pure, E_vap,
                    E_mixt, E_tot);

  /*
   * TODO: voir si on ne doit pas faire mieux, mais a priori les variations de Tl et Tv
   * sont lentes par rapport au reste donc ea devrait aller.
   * DONE: il y a manifestement un pb ici, car on ne peut pas avoir acces e Tv(n+1) encore,
   * donc il faut stocker Tv(n-1)
   */

  // Impression dans le fichier temperature_bulles.out
  if (Process::je_suis_maitre())
    {
      int reset = (!ref_ijk_ft_->get_reprise()) && (ref_ijk_ft_->get_tstep()==0);
      SFichier fic=Ouvrir_fichier(Nom("_source_temperature_")+Nom(rang_)+Nom("_bulles.out"),
                                  "tstep\ttime\tTl\tTv\tEtot\tElpu\tEl\tEvpu\tEv\tEm",
                                  reset);
      // la derivee_acceleration n'est connue que sur le maitre
      fic<< ref_ijk_ft_->get_tstep()<<" "<< ref_ijk_ft_->get_current_time() <<" "<< Tl << " " << Tv << " " << E_tot;
      fic<< " " << E_liq_pure <<  " " << E_liq;
      fic<< " " << E_vap_pure <<  " " << E_vap;
      fic<< " " << E_mixt;
      fic<<finl;
      fic.close();
    }
}

void IJK_Thermal_base::calculer_energies(double& E_liq_pure,
                                         double& E_liq,
                                         double& E_vap_pure,
                                         double& E_vap,
                                         double& E_mixt,
                                         double& E_tot)
{
  const int nk = temperature_.nk();
  const int ni = temperature_.ni();
  const int nj = temperature_.nj();
  const IJK_Field_double& chi = ref_ijk_ft_->itfce().I(); // rappel : chi vaut 1. dans le liquide et 0 dans la vapeur
  const double rhocpl = ref_ijk_ft_->get_rho_l()*cp_liquid_;
  const double rhocpv = ref_ijk_ft_->get_rho_v()*cp_vapour_;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double indic=chi(i, j, k);
          const double  rhocpm = indic*rhocpl+(1.-indic)*rhocpv;
          const double T = temperature_(i,j,k);
          E_tot += T * rhocpm;
          E_liq += indic * rhocpl * T;
          E_vap += (1.-indic) * rhocpv * T;
          if (std::fabs(indic)<1.e-8)
            {
              // vap pure
              E_vap_pure += rhocpv * T;
            }
          else if (std::fabs(1.-indic)<1.e-8)
            {
              // liq pure
              E_liq_pure += rhocpl * T;
            }
          else
            {
              // mixte :
              E_mixt += rhocpm * T;
            }
        }

  const int ntot = (temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_I)
                    * temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_J)
                    * temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K));
  E_vap_pure = Process::mp_sum(E_vap_pure)/ntot;
  E_liq_pure = Process::mp_sum(E_liq_pure)/ntot;
  E_vap = Process::mp_sum(E_vap)/ntot;
  E_liq = Process::mp_sum(E_liq)/ntot;
  E_tot = Process::mp_sum(E_tot)/ntot;
  E_mixt = Process::mp_sum(E_mixt)/ntot;
}

//void IJK_Thermal_base::calculer_temperature_physique_T_dummy()
//{
//  const int nk = d_temperature_.nk();
//  const int ni = d_temperature_.ni();
//  const int nj = d_temperature_.nj();
//
//  for (int k = 0; k < nk; k++)
//    {
//      for (int j = 0; j < nj; j++)
//        {
//          for (int i = 0; i < ni; i++)
//            {
//              temperature_physique_T_(i,j,k) = temperature_(i,j,k);
//            }
//        }
//    }
//  temperature_physique_T_.echange_espace_virtuel(temperature_physique_T_.ghost());
//  DebogIJK::verifier("temperature_physique_T", temperature_physique_T_);
//  return;
//}

//void IJK_Thermal_base::calculer_temperature_adimensionnelle_theta(const IJK_Field_double&  vx, const double wall_flux)
//{
//  if(wall_flux_)
//    {
//      const IJK_Splitting& splitting = temperature_.get_splitting();
//      const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
//      const double Lz = geom.get_domain_length(DIRECTION_K);
//      const double h = Lz/2.;
//      const double q_w = wall_flux;
//      const int kmin = temperature_.get_splitting().get_offset_local(DIRECTION_K);
//      const int kmax = splitting.get_grid_geometry().get_nb_elem_tot(DIRECTION_K);
//      const int nk = temperature_.nk();
//      const int ni = temperature_.ni();
//      const int nj = temperature_.nj();
//      double T_wall = 0;
//      T_wall = calculer_variable_wall(temperature_, cp_, ref_ijk_ft_->rho_field_, kmin, kmax);
//      /*   if (Process::je_suis_maitre())
//           {
//           T_wall = calculer_variable_wall(temperature_, cp_, ref_ijk_ft_->rho_field_, kmin, kmax);
//           Cerr << "calcul de T_wall sur maitre" << finl;
//           }
//           envoyer_broadcast(T_wall, 0); */
//      /*
//       * if(kmin+ nk == kmax) //|| (kmin ==0)
//       {
//       rank = Process::me();
//       T_wall = calculer_variable_wall(temperature_, cp_, ref_ijk_ft_->rho_field_, kmin, kmax);
//       envoyer_broadcast(T_wall, rank);
//       }
//       */
//      for (int k = 0; k < nk; k++)
//        {
//          //  const double T_mean = compute_spatial_mean(vx, temperature_, cp_, ref_ijk_ft_->rho_field_, kmin, nktot, k);
//          for (int j = 0; j < nj; j++)
//            {
//              for (int i = 0; i < ni; i++)
//                {
//                  double theta = T_wall - temperature_(i,j,k);
//                  // double theta = T_wall -T_mean  ;
//                  // temperature_adimensionnelle_theta_(i,j,k) = theta/theta_tau;
//                  const double lambda_l = lambda_liquid_;
//                  temperature_adimensionnelle_theta_(i,j,k) = theta*lambda_l/q_w/h;
//                  //    temperature_adimensionnelle_theta_(i,j,k) = -T_mean*lambda_l/q_w/h;
//                }
//            }
//        }
//      return;
//    }
//  else
//    {
//      Cerr << "no_source_for_temperature" << finl;
//      return;
//    }
//}


//void IJK_Thermal_base::calculer_Nusselt(const IJK_Field_double& vx)
//{
//  const double theta_adim_moy = calculer_temperature_adimensionnelle_theta_moy(vx, temperature_adimensionnelle_theta_,
//                                                                               cp_,ref_ijk_ft_->rho_field_);
//  double Nu = 0.;
//  if (std::fabs(theta_adim_moy)>1.e-10)
//    Nu = 2./theta_adim_moy;
//  const double rho_cp_u_moy = calculer_rho_cp_u_moyen(vx, cp_,ref_ijk_ft_->rho_field_);
//  // Impression dans le fichier source_temperature.out
//  if (Process::je_suis_maitre())
//    {
//      int reset = (!ref_ijk_ft_->get_reprise()) && (ref_ijk_ft_->get_tstep()_==0);
//      const Nom name = Nom("_temperature_")+Nom(rang_)+Nom(".out");
//      SFichier fic=Ouvrir_fichier(name,
//                                  "tstep\ttime\ttheta_adim_moy\tNu\trho_cp_u",
//                                  reset);
//
//      fic<< ref_ijk_ft_->get_tstep()_ <<" "<< ref_ijk_ft_->get_current_time() <<" "<< theta_adim_moy <<" "<< Nu << " " << rho_cp_u_moy  << finl ;
//      fic.close();
//    }
//}

void IJK_Thermal_base::set_field_T_ana()
{
  Cerr << "Setting analytical temperature "<< rang_ <<" field to "<< expression_T_ana_ << finl;
  set_field_data(temperature_ana_, expression_T_ana_, ref_ijk_ft_->get_current_time());
}

void IJK_Thermal_base::calculer_ecart_T_ana()
{
  if (liste_post_instantanes_.contient_("ECART_T_ANA"))
    {
      if (!liste_post_instantanes_.contient_("TEMPERATURE_ANA"))
        {
          set_field_data(temperature_ana_, expression_T_ana_, ref_ijk_ft_->get_current_time());
        }
      // do some work

      double ct = ref_ijk_ft_->get_current_time();
      Cerr << "GB: ERROR T FIELD " << ct;
      double err = 0.;
      set_field_data(temperature_ana_, expression_T_ana_, ct);
      const int ni = temperature_.ni();
      const int nj = temperature_.nj();
      const int nk = temperature_.nk();
      const int ntot=Process::mp_sum(ni*nj*nk);
      // La temperature est definie a une constante pres:
      //const double cst_temp = temperature_ana_(0,0,0) - curseur->temperature_(0,0,0);
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double val =   temperature_ana_(i,j,k) - temperature_(i,j,k); //- cst_temp;
              ecart_t_ana_(i,j,k) = val;
              err += val*val;
            }
      err=Process::mp_sum(err);
      err=sqrt(err/ntot);
      Cerr << " " << err ;
      if (!Process::je_suis_maitre())
        {
          Process::Journal() << "IJK_FT_Post::posttraiter_champs_instantanes : Champ ECART_T_ANA sur ce proc (ni,nj,nk,ntot):"
                             << " " << ni << " " << nj << " " << nk << " " << ntot << finl;
        }
      ecart_t_ana_.echange_espace_virtuel(ecart_t_ana_.ghost());
      Cerr << finl ;
      //  n++,dumplata_scalar(lata_name,"ECART_T_ANA", ecart_t_ana_, latastep);
    }
}

void IJK_Thermal_base::calculer_gradient_temperature(const IJK_Field_double& temperature, FixedVector<IJK_Field_double, 3>& grad_T)
{
  /*
   * Re-initialisation of the gradient vector
   */
  for (int dir = 0; dir < 3; dir++)
    {
      grad_T[dir].data() = 0.;
    }
//  add_gradient_temperature(temperature, 1. /*constant*/,  grad_T[0], grad_T[1], grad_T[2], boundary_conditions_, lambda_);
  for (int dir = 0; dir < 3; dir++)
    {
      grad_T[dir].echange_espace_virtuel(1);
    }
}

// Results are intensive (ie prop to area)
// Method fills storage_ so it changes the class
// Les interfaces connaissent le splitting_ft_ donc la correspondance doit etre appliquee au splitting ft pour convertir :
// convert_packed_to_ijk_cell.
// Donc il faut un champ de T etendu...

double IJK_Thermal_base::compute_global_energy(const IJK_Field_double& temperature)
{
  global_energy_ = 0.;
  const IJK_Field_double& indic = ref_ijk_ft_->itfce().I();
  const double rhocpl = get_rhocp_l();
  const double rhocpv = get_rhocp_v();
  const int nx = temperature.ni();
  const int ny = temperature.nj();
  const int nz = temperature.nk();
  // To be sure we're on a regular mesh
  assert(indic.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_K) >0);
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        {
          double chi_l = indic(i,j,k);
          global_energy_ += (rhocpl * chi_l + (1.- chi_l) * rhocpv) * temperature(i,j,k);
        }
  const int ntot = temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_I)
                   *temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_J)
                   *temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  global_energy_ = mp_sum(global_energy_)/(double)(ntot);
  return global_energy_;
}

/*
 * Getters and setters
 */

double IJK_Thermal_base::get_rhocp_l() const
{
  return  ref_ijk_ft_->get_rho_l() * cp_liquid_;
}

double IJK_Thermal_base::get_rhocp_v() const
{
  return  ref_ijk_ft_->get_rho_v() * cp_vapour_;
}

/*
 * Methods that do not belong to the class
 */

// From DNS_QC; Vectorize code later?
int IJK_Thermal_base::calculer_k_pour_bord(const IJK_Field_double& temperature, const bool bord_kmax)
{
  const int kmin = temperature.get_splitting().get_offset_local(DIRECTION_K);
  const int nktot = temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  int k;
  // calcul l'indice k de la couche de mailles voisine du bord. Si je n'ai pas de bord, on met k = -1
  if (!bord_kmax)
    {
      // on veut le bord "k_global = 0"
      if (kmin == 0)
        {
          // ce bord est chez moi... et il est en k=0
          k = 0;
        }
      else
        {
          // ce bord n'est pas chez moi
          k = -1;
        }
    }
  else
    {
      // on veut le bord kmax
      if (kmin + temperature.nk() == nktot)
        {
          // ce bord est chez moi... et il est en k= truc...
          k = temperature.nk() - 1;
        }
      else
        {
          k = -1;
        }
    }
  return k;
}

// From DNS_QC; Vectorize code later?
// valeur de retour: indice local du plan de temperature voisin utilise,
//  -1 si on n'a pas le bord sur ce processeur
// Calcule l'integrale sur chaque face du bord demande du flux de chaleur a travers la face
// positif si le flux va vers les k positifs.
int IJK_Thermal_base::calculer_flux_thermique_bord(const IJK_Field_double& temperature,
                                                   const double lambda_de_t_paroi,
                                                   const double T_paroi_impose,
                                                   IJK_Field_local_double& flux_bord,
                                                   const bool bord_kmax)
{
  const int kmin = temperature.get_splitting().get_offset_local(DIRECTION_K);
  int k = calculer_k_pour_bord(temperature, bord_kmax);
  if (k == -1)
    return k;

  // redimensionne flux_bord avec ni * nj:
  const int ni = temperature.ni(); // nombre d'element local sur ce processeur
  const int nj = temperature.nj();
  flux_bord.allocate(ni, nj, 1 /* 1 seule couche de valeurs en k */, 0 /* pas d'elements fantomes */);

  const IJK_Grid_Geometry& geometry = temperature.get_splitting().get_grid_geometry();
  const double delta_k = geometry.get_delta(DIRECTION_K)[k + kmin]; // k+kmin est l'indice global de la maille locale k
  double facteur = 2.0 / delta_k * geometry.get_constant_delta(DIRECTION_I) * geometry.get_constant_delta(DIRECTION_J);
  if (bord_kmax)
    facteur *= -1.; // he he... je vous laisse reflechir a ca :)
  // nan c'est pas simpa: la convention dans l'operateur de diffusion est
  //   d/dt = flux(i,j) - flux(i+1,j) + ... + flux(i,j,k) - flux(i,j,k+1)
  // donc si la paroi inferieure (k=0) est plus froide que le fluide, il faut que le flux stocke soit negatif.
  // et si la paroi inferieure (k=kmax) est plus froide que le fluide, il faut que le flux stocke soit positif.
  for (int j = 0; j < nj; j++)
    {
      for (int i = 0; i < ni; i++)
        {
          // Temperature de la maille voisine
          const double t = temperature(i,j,k);
          // le flux est positif s'il va vers les k croissants
          flux_bord(i,j,0) = (T_paroi_impose - t) * lambda_de_t_paroi * facteur;
        }
    }
  return k;
}
int IJK_Thermal_base::imposer_flux_thermique_bord(const IJK_Field_double& temperature,
                                                  const double flux_paroi_impose,
                                                  IJK_Field_local_double& flux_bord,
                                                  const bool bord_kmax)
{
  int k = calculer_k_pour_bord(temperature, bord_kmax);
  if (k == -1)
    return k;

  // redimensionne flux_bord avec ni * nj:
  const int ni = temperature.ni(); // nombre d'element local sur ce processeur
  const int nj = temperature.nj();
  flux_bord.allocate(ni, nj, 1 /* 1 seule couche de valeurs en k */, 0 /* pas d'elements fantomes */);
  //MR je multiplie le flux par la surface dxdy
  const IJK_Grid_Geometry& geometry = temperature.get_splitting().get_grid_geometry();

  double facteur = 1.* geometry.get_constant_delta(DIRECTION_I) * geometry.get_constant_delta(DIRECTION_J);
  if (bord_kmax)
    facteur *= -1.; // he he... je vous laisse reflechir a ca :)
  // nan c'est pas sympa: la convention dans l'operateur de diffusion est
  // d/dt = flux(i,j) - flux(i+1,j) + ... + flux(i,j,k) - flux(i,j,k+1)
  // donc si la paroi inferieure (k=0) est plus froide que le fluide, il faut que le flux stocke soit negatif.
  // et si la paroi superieure (k=kmax) est plus froide que le fluide, il faut que le flux stocke soit positif.
  for (int j = 0; j < nj; j++)
    {
      for (int i = 0; i < ni; i++)
        {
          // le flux est positif s'il va vers les k croissants
          flux_bord(i,j,0) = flux_paroi_impose * facteur;
        }
    }
  return k;
}
