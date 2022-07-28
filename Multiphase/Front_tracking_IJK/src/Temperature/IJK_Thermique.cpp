/****************************************************************************
* Copyright (c) 2019, CEA
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
// File      : IJK_Thermique.cpp
// Directory : $IJK_ROOT/src/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Thermique.h>
#include <Param.h>
#include <IJK_Navier_Stokes_tools.h>
#include <DebogIJK.h>
#include <stat_counters.h>
#include <IJK_FT.h>

Implemente_instanciable( IJK_Thermique, "IJK_Thermique", Objet_U ) ;

Sortie& IJK_Thermique::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  os<< "  {\n"
    << "    boundary_conditions {"  << "\n";
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

  os<< "    lambda_liquid " <<  lambda_liquid_ << "\n";
  os<< "    lambda_vapor " <<  lambda_vapor_ << "\n";
  os<< "    cp_liquid " <<  cp_liquid_ << "\n";
  os<< "    cp_vapor " <<  cp_vapor_ << "\n";
  os<< "    type_T_source " << type_T_source_ << "\n";
  if (type_T_source_=="SWARM")
    {
      os<< "      kl_source " <<  kl_ << "\n";
      os<< "      kv_source " <<  kv_ << "\n";
      os<< "      T0l_source " <<  T0l_ << "\n";
      os<< "      T0v_source " <<  T0v_ << "\n";
    }
  os<< "    fichier_reprise_temperature" << " " << fichier_reprise_temperature_  << "\n";
  os<< "    timestep_reprise_temperature" << " " << timestep_reprise_temperature_ << "\n";

  if (rho_cp_inv_)
    os<< "    rho_cp_inv \n";
  if (lambda_variable_)
    os<< "    lambda_variable \n";
  if (lambda_moy_arith_)
    os<< "    lambda_moy_arith \n";
  if (depracated_rho_cp_)
    os<< "    depracated_rho_cp \n";
  if (conserv_energy_global_)
    os<< "    conserv_energy_global \n";

  if ( expression_T_ana_!="??")
    os<< "    expression_T_ana" <<  " " << expression_T_ana_ << "\n";

  if( wall_flux_)
    os << "    wall_flux \n";

  if ( conv_temperature_negligible_)
    os<< "    conv_temperature_negligible \n ";

  if ( diff_temp_negligible_)
    os<< "    diff_temp_negligible \n";

  os<< "  }\n";
  return os;
}

Entree& IJK_Thermique::readOn( Entree& is )
{
  rang_ = 0; // default value
  cp_liquid_=1.;
  cp_vapor_=1.;
  lambda_liquid_=1.;
  lambda_vapor_=1.;
  fo_ = 1.;
  expression_T_init_="??";
  fichier_reprise_temperature_= "??";
  expression_T_ana_="??";
  type_T_source_="??";

// reprise_ = 0;
  timestep_reprise_temperature_=1;
  diff_temp_negligible_=0;
  conv_temperature_negligible_=0;
  //terme source variable
  lambda_variable_=0;

  type_temperature_convection_op_ = 3;  // Default value: 3 : Quick
  wall_flux_=0;
  lambda_moy_arith_=0;

// wall_flux_=0.;
//
//Ceci est une initialisation des derivees des temperatures moyenne de chaque phase
//Il n'est peut-etre pas pertinent de les mettre ici
  dTv_ = 0.;
  dTl_ = 1.;
  Tl_ = 0.;
  Tv_ = 1.;
  Tv0_ = 1.;  // Serait-ce plutot Tref (une temperature de reference pour reconstruire le champ dim??)
  kl_ = -100000000000000.;
  kv_ = -200000000000000.;
  T0v_ = 1.;
  T0l_ = 0.;

  type_temperature_convection_form_ = 1;  // Default value: 1 : non conservative

  Param param(que_suis_je());

  param.ajouter("cp_liquid", &cp_liquid_, Param::REQUIRED);
  param.ajouter("lambda_liquid", &lambda_liquid_, Param::REQUIRED);
  param.ajouter("cp_vapor", &cp_vapor_, Param::REQUIRED);
  param.ajouter("lambda_vapor", &lambda_vapor_, Param::REQUIRED);
  param.ajouter("fo", &fo_);
  param.ajouter("boundary_conditions", &boundary_conditions_, Param::REQUIRED);
  /*	  param.ajouter("type_convection_op", &type_convection_op_);
  	  param.dictionnaire("Quick",0);
  	  param.dictionnaire("Centre",1);
  	*/
  param.ajouter("kl_source", &kl_);
  param.ajouter("kv_source", &kv_);
  param.ajouter("T0l_source", &T0l_);
  param.ajouter("T0v_source", &T0v_);

  param.ajouter("expression_T_init", &expression_T_init_);
  param.ajouter("fichier_reprise_temperature", &fichier_reprise_temperature_);
  param.ajouter("timestep_reprise_temperature", &timestep_reprise_temperature_);
  param.ajouter_flag("conv_temperature_negligible", &conv_temperature_negligible_);
  param.ajouter("type_temperature_convection_op", &type_temperature_convection_op_);
  param.dictionnaire("Amont",1);
  param.dictionnaire("Quick",3);
  param.dictionnaire("Centre2",2);
  param.dictionnaire("Centre4",4);
  param.ajouter_flag("diff_temp_negligible", &diff_temp_negligible_);

  param.ajouter_flag("lambda_variable", &lambda_variable_);
  param.ajouter_flag("wall_flux", &wall_flux_);
  param.ajouter_flag("depracated_rho_cp", &depracated_rho_cp_);
  param.ajouter_flag("conserv_energy_global", &conserv_energy_global_);
  param.ajouter_flag("lambda_moy_arith", &lambda_moy_arith_);
  param.ajouter_flag("rho_cp_inv", &rho_cp_inv_);

  // Expression analytique de la temperature
  param.ajouter("expression_T_ana", &expression_T_ana_);
  param.ajouter("type_T_source", &type_T_source_);
  param.ajouter("expression_source_temperature", &expression_source_temperature_);

  param.ajouter("type_temperature_convection_form", &type_temperature_convection_form_);
  param.dictionnaire("non conservative",1);
  param.dictionnaire("conservative",2);

  param.lire_avec_accolades(is);
  Cout << "IJK_Thermique::readOn : Parameters summary. " << finl;
  printOn(Cout);

  return is;
}

int IJK_Thermique::initialize(const IJK_Splitting& splitting, const int idx)
{
  Cout << que_suis_je() << "::initialize()" << finl;
  rang_ = idx;
  int nalloc = 0;
  diffusion_temperature_op_.initialize(splitting);

  switch(type_temperature_convection_op_)
    {
    case 1:
      temperature_convection_op_amont_.initialize(splitting);
      break;
    case 2:
      temperature_convection_op_centre2_.initialize(splitting);
      break;
    case 3:
      temperature_convection_op_quick_.initialize(splitting);
      if (type_temperature_convection_form_ == 2)
        {
          rho_cp_convection_op_quick_.initialize(splitting);
        }
      break;
    case 4:
      temperature_convection_op_centre4_.initialize(splitting);
      break;
    default:
      Cerr << "Undefined operator for the convection of the temperature. " << finl;
      Process::exit();
    }

  temperature_.allocate(splitting, IJK_Splitting::ELEM, 2);
  cp_.allocate(splitting, IJK_Splitting::ELEM, 2);
  lambda_.allocate(splitting, IJK_Splitting::ELEM, 1);
  d_temperature_.allocate(splitting, IJK_Splitting::ELEM, 2);// 1
  div_lambda_grad_T_volume_.allocate(splitting, IJK_Splitting::ELEM, 0);;
  nalloc += 5;

  if ((ref_ijk_ft_.non_nul()) and (!ref_ijk_ft_->disable_diphasique_))
    {
      Cout << "Allocating fields temperature_ft_ and storage" << finl;
      allocate_cell_vector(storage_, ref_ijk_ft_->get_splitting_ft(), 1);
      temperature_ft_.allocate(ref_ijk_ft_->get_splitting_ft(), IJK_Splitting::ELEM, 1);
      nalloc += 4;
    }

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
          type_T_source_ = "dabiri";
        }
    }
  if (liste_post_instantanes_.contient_("TEMPERATURE_ADIM_BULLES"))
    {
      temperature_adim_bulles_.allocate(splitting, IJK_Splitting::ELEM, 2);
      nalloc += 1;
    }

  // Check that pointer is not null:
  if (ref_ijk_ft_.non_nul() && ref_ijk_ft_->get_time_scheme()== ref_ijk_ft_->RK3_FT)
    {
      RK3_F_temperature_.allocate(splitting, IJK_Splitting::ELEM, 0);
      nalloc +=1;
    }

  if (fichier_reprise_temperature_ == "??")   // si on ne fait pas une reprise on initialise V
    {
      if (expression_T_init_ != "??")
        {
          Cout << "Temperature initialization from expression \nTini = " << expression_T_init_ << finl;
          set_field_data(temperature_, expression_T_init_, ref_ijk_ft_->indicatrice_ns_, 0.);
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
    allocate_velocity(grad_T_, splitting, 1);
  nalloc +=3 ;

  if ((conserv_energy_global_) || (type_temperature_convection_form_==2))
    {
      T_rust_.allocate(splitting, IJK_Splitting::ELEM, 1);
      rho_cp_.allocate(splitting, IJK_Splitting::ELEM, 2);
      nalloc += 2;
    }

  // Compute initial energy :
  if (conserv_energy_global_)
    {
      E0_ = compute_global_energy(temperature_);
      d_T_rustine_.allocate(splitting, IJK_Splitting::ELEM, 1);
      nalloc += 1;
      if (ref_ijk_ft_.non_nul() && ref_ijk_ft_->get_time_scheme()== ref_ijk_ft_->RK3_FT)
        {
          RK3_F_rustine_.allocate(splitting, IJK_Splitting::ELEM, 0);
          nalloc +=1;
        }
      Cout << "Initial energy at time t=" << ref_ijk_ft_->get_current_time() << " is " << E0_ << " [W.m-3]." << finl;
      Cerr << "Initial energy at time t=" << ref_ijk_ft_->get_current_time() << " is " << E0_ << " [W.m-3]." << finl;
    }
  if (type_temperature_convection_form_==2)
    {
      rho_cp_T_.allocate(splitting, IJK_Splitting::ELEM, 2);
      div_rho_cp_T_.allocate(splitting, IJK_Splitting::ELEM, 0);
      nalloc += 2;
    }
  Cout << "End of " << que_suis_je() << "::initialize()" << finl;
  return nalloc;
}

void IJK_Thermique::associer(const IJK_FT_double& ijk_ft)
{
  ref_ijk_ft_ = ijk_ft;
  liste_post_instantanes_ = ijk_ft.post_.get_liste_post_instantanes();
}

void IJK_Thermique::sauvegarder_temperature(Nom& lata_name, int idx)
{
  fichier_reprise_temperature_ = lata_name;
  timestep_reprise_temperature_ = 1;
  dumplata_scalar(lata_name,Nom("TEMPERATURE_")+Nom(idx) , temperature_, 0 /*we store a 0 */);
}

void IJK_Thermique::update_thermal_properties()
{
  const double ene_ini = compute_global_energy();
  const IJK_Field_double& indic = ref_ijk_ft_->indicatrice_ns_;
  // Nombre de mailles du domaine NS :
  const int nx = indic.ni();
  const int ny = indic.nj();
  const int nz = indic.nk();
  const double rho_l = ref_ijk_ft_->rho_liquide_;
  const double rho_v = ref_ijk_ft_->rho_vapeur_;
  const bool geometric_mean = ((!lambda_moy_arith_) and (lambda_liquid_ > DMINFLOAT) and (lambda_vapor_ >DMINFLOAT));
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        {
          double chi_l = indic(i,j,k);
          cp_(i,j,k) = cp_liquid_ * chi_l + (1.- chi_l) * cp_vapor_;
          if (geometric_mean)
            lambda_(i,j,k) = lambda_liquid_ * lambda_vapor_ / ((1 - chi_l) * lambda_liquid_ + chi_l * lambda_vapor_);
          else
            lambda_(i,j,k) = lambda_liquid_  * chi_l + (1.- chi_l) * lambda_vapor_ ;
          if ((type_temperature_convection_form_ == 2) || (conserv_energy_global_))
            {
              rho_cp_(i,j,k) = rho_l*cp_liquid_*indic(i,j,k) + rho_v*cp_vapor_*(1-indic(i,j,k));
            }
          if (type_temperature_convection_form_ == 2)
            {
              rho_cp_T_(i,j,k) = (rho_l*cp_liquid_*indic(i,j,k) + rho_v*cp_vapor_*(1-indic(i,j,k)))*temperature_(i,j,k);
            }

        }
  cp_.echange_espace_virtuel(cp_.ghost());
  lambda_.echange_espace_virtuel(lambda_.ghost());
  if ((type_temperature_convection_form_ == 2) || (conserv_energy_global_))
    {
      rho_cp_.echange_espace_virtuel(rho_cp_.ghost());
    }
  if (type_temperature_convection_form_ == 2)
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
double IJK_Thermique::compute_timestep(const double timestep,
                                       const double rho_l, const double rho_v,
                                       const double dxmin) const
{
  // alpha = lambda/(rho*cp)
  const double alpha_max = std::max(lambda_liquid_/(rho_l*cp_liquid_), lambda_vapor_/(rho_v*cp_vapor_));
  const double dt_fo  = dxmin*dxmin/(alpha_max + 1.e-20) * fo_ * (1./6.); // Attention 0.125 vient du 3D. (1/6 au lieu de 1/8)
  return dt_fo;
}


// From DNS_QC; Vectorize code later?
int calculer_k_pour_bord(const IJK_Field_double& temperature, const bool bord_kmax)
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
int calculer_flux_thermique_bord(const IJK_Field_double& temperature,
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
int imposer_flux_thermique_bord(const IJK_Field_double& temperature,
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
  // nan c'est pas simpa: la convention dans l'operateur de diffusion est
  //   d/dt = flux(i,j) - flux(i+1,j) + ... + flux(i,j,k) - flux(i,j,k+1)
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



void IJK_Thermique::euler_time_step(const double timestep)
{
  calculer_dT(ref_ijk_ft_->velocity_);
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

void IJK_Thermique::rk3_sub_step(const int rk_step, const double total_timestep,
                                 const double fractionnal_timestep, const double time)
{
  calculer_dT(ref_ijk_ft_->velocity_);
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


// Mettre rk_step = -1 si schema temps different de rk3.
void IJK_Thermique::calculer_dT(const FixedVector<IJK_Field_double, 3>& velocity)
{
  const double current_time = ref_ijk_ft_->get_current_time();
  const double ene_ini = compute_global_energy(d_temperature_);
  switch (type_temperature_convection_form_)
    {
    case 1:
      compute_temperature_convection(velocity);
      break;
    case 2:
      compute_temperature_convection_conservative(velocity);
      break;
    default:
      Cerr << "Pb du switch de convection" << finl;
    }
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

  //Ici on stocke T_rustine pour plus tard
  if (conserv_energy_global_)
    {
      compute_T_rust(velocity);
    }
  return;
}

// Convect temperature field by velocity.
// The output is stored in d_temperature_ (it is a volume integral over the CV)
void IJK_Thermique::compute_temperature_convection(const FixedVector<IJK_Field_double, 3>& velocity)
{
  static Stat_Counter_Id cnt_conv_temp = statistiques().new_counter(1, "FT convection rho");
  statistiques().begin_count(cnt_conv_temp);
  if (conv_temperature_negligible_)
    {
      d_temperature_.data()=0;
    }
  else
    {
      switch (type_temperature_convection_op_)
        {
        case 1:
          temperature_convection_op_amont_.calculer(temperature_, temperature_, temperature_, velocity[0], velocity[1], velocity[2], d_temperature_, d_temperature_, d_temperature_);
          break;
        case 2:
          temperature_convection_op_centre2_.calculer(temperature_, velocity[0], velocity[1], velocity[2], d_temperature_);
          break;
        case 3:
          temperature_convection_op_quick_.calculer(temperature_, velocity[0], velocity[1], velocity[2], d_temperature_);
          break;
        case 4:
          temperature_convection_op_centre4_.calculer(temperature_,temperature_,temperature_, velocity[0], velocity[1], velocity[2], d_temperature_, d_temperature_, d_temperature_);
          break;

        default:
          Cerr << "Unknown convection operator for the temperature." << finl;
          Process::exit();
        }
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
              d_temperature_(i,j,k) /=vol ;
            }
    }
  statistiques().end_count(cnt_conv_temp);
  DebogIJK::verifier("op_conv(rho)", d_temperature_);
  return;
}

void IJK_Thermique::compute_temperature_convection_conservative(const FixedVector<IJK_Field_double, 3>& velocity)
{
  static Stat_Counter_Id cnt_conv_temp = statistiques().new_counter(1, "FT convection rho");
  statistiques().begin_count(cnt_conv_temp);
  if (conv_temperature_negligible_)
    {
      d_temperature_.data()=0;
    }
  else
    {
      switch (type_temperature_convection_op_)
        {
        case 1:
          Cerr << "Operateur amont non implemente" << finl;
          break;
        case 2:
          Cerr << "Operateur centre non implemente" << finl;
          break;
        case 3:
          rho_cp_convection_op_quick_.calculer(rho_cp_T_, velocity[0], ref_ijk_ft_->indicatrice_ns_, velocity[1], velocity[2], div_rho_cp_T_);
          break;
        case 4:
          Cerr << "Operateur centre4 non implemente" << finl;
          break;

        default:
          Cerr << "Unknown convection operator for the temperature." << finl;
          Process::exit();
        }
      const int ni = d_temperature_.ni();
      const int nj = d_temperature_.nj();
      const int nk = d_temperature_.nk();
      const IJK_Grid_Geometry& geom = d_temperature_.get_splitting().get_grid_geometry();
      const double dx = geom.get_constant_delta(DIRECTION_I);
      const double dy = geom.get_constant_delta(DIRECTION_J);
      const double dz = geom.get_constant_delta(DIRECTION_K);
      const double vol = dx*dy*dz;
      const double rho_l = ref_ijk_ft_->rho_liquide_;
      const double rho_v = ref_ijk_ft_->rho_vapeur_;
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double I = ref_ijk_ft_->indicatrice_ns_(i,j,k);
              // dT = 1/rho_cp_h * div( rho*cp*T*v )
              div_rho_cp_T_(i,j,k) /= vol;
              d_temperature_(i,j,k) = 1/((I*rho_l*cp_liquid_) + (1.-I)*rho_v*cp_vapor_)*div_rho_cp_T_(i,j,k);
            }
      // on met a jour T_rust_
      compute_T_rust(velocity);
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              // dT = 1/rho_cp_h * div( rho*cp*T*v ) - T/rho_cp_h *div( rho*cp*v)
              d_temperature_(i,j,k) -= T_rust_(i,j,k) ;
            }
      d_temperature_.echange_espace_virtuel(d_temperature_.ghost());
    }
  statistiques().end_count(cnt_conv_temp);
  DebogIJK::verifier("op_conv(rho)", d_temperature_);
  return;
}

void IJK_Thermique::add_temperature_diffusion()
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
  lambda_.echange_espace_virtuel(lambda_.ghost());
  DebogIJK::verifier("temp", temperature_);
  DebogIJK::verifier("lambda", lambda_);

  // Performance counters:
  static Stat_Counter_Id cnt_diff_temp = statistiques().new_counter(1, "FT diffusion temperature");
  statistiques().begin_count(cnt_diff_temp);
  if (diff_temp_negligible_) // si diffusion negligeable
    {
      div_lambda_grad_T_volume_.data()=0;
    }
  else
    {
      diffusion_temperature_op_.calculer(temperature_,
                                         lambda_,
                                         div_lambda_grad_T_volume_,
                                         boundary_flux_kmin_, boundary_flux_kmax_);
    }
  // Update d_temperature
  //d_temperature_ +=div_lambda_grad_T_volume_;

  const int ni = d_temperature_.ni();
  const int nj = d_temperature_.nj();
  const int nk = d_temperature_.nk();
  const IJK_Grid_Geometry& geom = d_temperature_.get_splitting().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  const double vol = dx*dy*dz;
  const double vol_inv = 1./vol;
  const double rhocp_l = ref_ijk_ft_->rho_liquide_*cp_liquid_;
  const double rhocp_v = ref_ijk_ft_->rho_vapeur_*cp_vapor_;
  const bool geometric_mean = ((rho_cp_inv_) and (rhocp_l > DMINFLOAT) and (rhocp_v >DMINFLOAT));
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          if (geometric_mean)
            {
              const double chi_l = ref_ijk_ft_->indicatrice_ns_(i,j,k);
              const double rhocpV_inv = (chi_l/rhocp_l  + (1 - chi_l)/rhocp_v) * vol_inv;
              const double ope = div_lambda_grad_T_volume_(i,j,k);
              const double resu = ope*rhocpV_inv;
              d_temperature_(i,j,k) +=resu ;
            }
          else
            {
              double rhocpV = 0;
              if (depracated_rho_cp_)
                {
                  const double rho = ref_ijk_ft_->rho_field_(i,j,k);
                  const double cp = cp_(i,j,k);
                  rhocpV = rho * cp *vol;
                }
              else
                {
                  const double chi_l = ref_ijk_ft_->indicatrice_ns_(i,j,k);
                  rhocpV = (rhocp_l * chi_l + rhocp_v * (1 - chi_l)) *  vol;
                }
              const double ope = div_lambda_grad_T_volume_(i,j,k);
              const double resu = ope/rhocpV;
              d_temperature_(i,j,k) +=resu ;
            }
        }

  statistiques().end_count(cnt_diff_temp);
  DebogIJK::verifier("div_lambda_grad_T_volume", div_lambda_grad_T_volume_);
// Cerr << "diff_temp" << " " <<  d_temperature_(1,1,1) << finl;
//TODO: deplacer ce bloc ailleurs...
  if (liste_post_instantanes_.contient_("GRAD_T")
      || (ref_ijk_ft_->t_debut_statistiques() <  1.e10 ))
    {
      calculer_gradient_temperature(temperature_,grad_T_);
    }
}

//////////////////////////////////////////

void IJK_Thermique::add_temperature_source()
// DONE: ajouter un mot clef pour changer le type de source (constante, ponderee par la vitesse, ...)
{
  static Stat_Counter_Id cnt_source_temp = statistiques().new_counter(1, "FT source temperature");
  statistiques().begin_count(cnt_source_temp);

  // dans le cas ou les flux entrants et sortants sont identiques :
  // DONE: changer cette condition non adaptee
  if (type_T_source_!="??")
    {
      const IJK_Field_double&  vx = ref_ijk_ft_->velocity_[DIRECTION_I];
      double rho_cp_u_moy = calculer_rho_cp_u_moyen(vx, cp_,ref_ijk_ft_->rho_field_);
      // double rho_cp_moy = calculer_rho_cp_moyen(cp_, ref_ijk_ft_->rho_field_);
      const IJK_Splitting& splitting = temperature_.get_splitting();
      const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
      const double dx =geom.get_constant_delta(DIRECTION_I);
      const double dy = geom.get_constant_delta(DIRECTION_J);
      const double dz = geom.get_constant_delta(DIRECTION_K);
      const double vol = dx*dy*dz;
      const double lz = geom.get_domain_length(DIRECTION_K) ;
      const double h =lz/2.;
      // double volume = 1.;
      // for (int i = 0; i < 3; i++)
      //   volume *= splitting.get_grid_geometry().get_constant_delta(i);
      //IJK_Field_double& dT = d_temperature_(i,j,k);
      const int nk = d_temperature_.nk();
      const int ni = d_temperature_.ni();
      const int nj = d_temperature_.nj();


      // TODO: faire une methode calculer_rho_cp
      //debut if source = ponderee par la vitesse
      if (type_T_source_=="dabiri")
        {

          const double wall_flux = boundary_conditions_.get_flux_kmax();
          const double qw = wall_flux;
          const double dTm = 2*qw/rho_cp_u_moy ; // TODO: faux, il manque un facteur 1/h pour homogeneite

          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  const double rho = ref_ijk_ft_->rho_field_(i,j,k);
                  const double cp = cp_(i,j,k);
                  const double u = (vx(i,j,k) +vx(i+1,j,k))/2;
                  const double rho_cp_u = rho*cp*u;
                  // const double lambda = lambda_(i,j,k);
                  //const double f = rho_cp_u; // - div (lambda);
                  if(lambda_variable_)
                    {
                      const double div_lambda = (lambda_(i+1,j,k)-lambda_(i-1,j,k))/(2*dx);
                      source_temperature_(i,j,k) = (rho_cp_u-div_lambda)*dTm;
                    }
                  else
                    {
                      source_temperature_(i,j,k) = rho_cp_u*dTm;
                    }
                  // const double lambda = lambda_(i,j,k);
                  const double Sc = qw/h ;//qw/lambda*h;
                  const double source = (source_temperature_(i,j,k)-Sc)*vol;//-Sc)*volume; TODO: faux, que vient faire Sc en plus de source ?

                  d_temperature_(i,j,k) += source;

                }
          //
          calculer_temperature_physique_T(vx, dTm);
          calculer_temperature_adimensionnelle_theta(vx, wall_flux);
          calculer_Nusselt(vx);
          return;
        }

      // TODO: remplacer dabiri par patch_dabiri apres verif
      else if (type_T_source_=="patch_dabiri")
        {

          Cerr << "Type de source : patch_dabiri" << finl;
          const double wall_flux = boundary_conditions_.get_flux_kmax();
          const double qw = wall_flux;
          const double dTm = -2*qw/(2*h*rho_cp_u_moy) ;

          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  //const double chi = ref_ijk_ft_->indicatrice_ns_(i,j,k);
                  //const double rho_cp_l = ref_ijk_ft_->rho_liquide_ * cp_liquid_;
                  //const double rho_cp_v = ref_ijk_ft_->rho_vapeur_ * cp_vapor_;
                  const double u = (vx(i,j,k) +vx(i+1,j,k))/2;
                  //const double rho_cp_u = (chi*rho_cp_l + (1-chi)*rho_cp_v)*u;
                  if(lambda_variable_)
                    {
                      const double div_lambda = (lambda_(i+1,j,k)-lambda_(i-1,j,k))/(2*dx);
                      source_temperature_(i,j,k) = (u-div_lambda)*dTm;
                    }
                  else
                    {
                      source_temperature_(i,j,k) = u*dTm;
                    }
                  const double source = source_temperature_(i,j,k);

                  d_temperature_(i,j,k) += source;

                }
          //
          calculer_temperature_physique_T(vx, dTm);
          calculer_temperature_adimensionnelle_theta(vx, wall_flux);
          calculer_Nusselt(vx);
          return;
        }
      else if (type_T_source_=="unweighted_dabiri")
        {
          // Sans la ponderation de S par u (S=rhocpu.../<rhocpu>), cela permet d'avoir une source uniforme!
          // Mais ce n'est plus le meme changement de variable, c'est quoi alors?
          // TODO: Que faire de rho_cp en diphasique? Moyen ou local?
          Cerr << "Type de source : unweighted_dabiri" << finl;
          const double wall_flux = boundary_conditions_.get_flux_kmax();
          const double qw = wall_flux;
          const double liquid_fraction = calculer_v_moyen(ref_ijk_ft_->indicatrice_ns_);
          const double rho_cp_l = ref_ijk_ft_->rho_liquide_ * cp_liquid_;
          const double rho_cp_v = ref_ijk_ft_->rho_vapeur_ * cp_vapor_;
          const double rhocp_moy =  rho_cp_l*liquid_fraction + rho_cp_v*(1-liquid_fraction);
          const double dTm = -2*qw/(2*h*rhocp_moy) ;
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  if(lambda_variable_)
                    {
                      Cerr << "Veut-on vraiment calculer une partie en lambda spatialement variable??" <<finl;
                      Cerr << "Exit at IJK_Thermique::add_temperature_source" << finl;
                      Process::exit();
                      // const double div_lambda = (lambda_(i+1,j,k)-lambda_(i-1,j,k))/(2*dx);
                      // source_temperature_(i,j,k) = (u-div_lambda)*dTm/<u>;
                    }
                  else
                    {
                      source_temperature_(i,j,k) = dTm;
                    }
                  const double source = source_temperature_(i,j,k);
                  d_temperature_(i,j,k) += source;
                }
          //
          calculer_temperature_physique_T(vx, dTm);
          calculer_temperature_adimensionnelle_theta(vx, wall_flux);
          calculer_Nusselt(vx);
          return;
        }
      //debut source = SWARM

      else if (type_T_source_=="SWARM")
        {
          //DONE: idem
          double Tv=0;
          double Tl=0;
          double Vv = 0;
          double Vl = 0;

          //calcul de Tv et Tl
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  const double chi = ref_ijk_ft_->indicatrice_ns_(i, j, k);
                  const double T = temperature_(i, j, k);
                  Tv += T*(1.-chi);
                  Vv += (1.-chi);
                  Tl += T*chi;
                  Vl += chi;
                }
          Tv = Process::mp_sum(Tv);
          Tl = Process::mp_sum(Tl);
          Vv = Process::mp_sum(Vv);
          Vl = Process::mp_sum(Vl);
          Tv /= Vv;
          Tl /= Vl;
          Cerr << "AY-test_source : " <<  Tv << finl;

          //calcul de dT a partir de l'expression suivante : dT = k(T - Tm)/(rho*cp)
          const double kv = kv_;
          const double kl = kl_;
          const double T0v = T0v_;
          const double T0l = T0l_;
          const double rho_cp_l = ref_ijk_ft_->rho_liquide_ * cp_liquid_;
          const double rho_cp_v = ref_ijk_ft_->rho_vapeur_ * cp_vapor_;
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  d_source_Tv_(i,j,k) = kv/rho_cp_v * (Tv - T0v);
                  d_source_Tl_(i,j,k) = kl/rho_cp_l * (Tl - T0l);
                }
          Cerr << "AY-test_source1 : " <<  Tv << finl;
          // TODO: remplacer euler_explicit_update par l'utilisation de timestep_ et utiliser d_source_Tv_ comme une constante
          for (int k = 0; k < nk; k++)
            {
              ref_ijk_ft_->euler_explicit_update(d_source_Tv_, source_temperature_v_, k);
              ref_ijk_ft_->euler_explicit_update(d_source_Tl_, source_temperature_l_, k);
            }
          Cerr << "AY-test_source2 : " <<  Tv << finl;
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  // chi vaut 1. dans le liquide et 0 dans les bulles
                  const double chi = ref_ijk_ft_->indicatrice_ns_(i, j, k);
                  source_temperature_(i,j,k) = (1.-chi)*source_temperature_v_(i,j,k) + chi*source_temperature_l_(i,j,k);
                }
          Cerr << "AY-test_source3 : " <<  Tv << finl;
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  d_temperature_(i,j,k) += source_temperature_(i,j,k)*vol;
                }
          Cerr << "AY-test_source 111 : " <<  source_temperature_(1,1,1) << finl;
          Cerr << "AY-test_source vol : " <<  vol << finl;
          Cerr << "source_temp " << " " << d_temperature_(1,1,1) << finl;
          const double current_time = ref_ijk_ft_->get_current_time();

          // GB : Ma comprehension est que ce ne sont pas des champs, mais un scalaire unique
          const double Sl = source_temperature_l_(0,0,0);
          const double Sv = source_temperature_v_(0,0,0);
          const double dSl = d_source_Tl_(0,0,0);
          const double dSv = d_source_Tv_(0,0,0);
          Cerr <<"[ThermalInfo-" << rang_ <<"] t/Tl/Tv/Sl/Sv/dSldt/dSvdt " << current_time << " "
               << Tl << " " << Tv << " "
               << Sl << " " << Sv << " "
               << dSl << " " << dSv
               <<finl;
          // si on utilise ca c'est pour remplir dans tous les cas d'utilisation d'une source le champs temperature_physique_T
          calculer_temperature_physique_T_dummy();

          return;
        }

    }
  // dans ce cas le ce ne sont pas des flux thermiques identiques
  else
    {
      Cerr << "no_source_for_temperature" << finl;

      // calculer_source_temperature_ana();
      //     Cerr << "source_temp" << " " << source_temperature_(24,24,0) << " " << source_temperature_(24,24,12) << " " <<  source_temperature_(24,24,23)  << finl;
      //     Cerr << " debit_moy " << rho_cp_u_moy << " vx_max " << (vx(0,0,12)+vx(1,0,12))/2
      //        << " vx_wall_0 " << (vx(0,0,0)+vx(1,0,0))/2 << " vx_wall_lz " << (vx(0,0,23)+vx(1,0,23))/2 << finl;
      return;
    }
  //AYM Ce bloc n'est jamais atteint non ?
  statistiques().end_count(cnt_source_temp);
  Cerr << "source_temp" << " " <<  source_temperature_(1,1,1) << finl;
  if (Process::je_suis_maitre())
    {
      const int nk = d_temperature_.nk();
      Cerr << "source_temp mil" << " " <<  source_temperature_(1,1,int(nk/2)) << finl;
    }
}

void IJK_Thermique::source_callback()
{
  if (liste_post_instantanes_.contient_("TEMPERATURE_ADIM_BULLES"))
    {
      calculer_temperature_adim_bulles();
    }
}


// renomer pour expliciter qu'il s'agit de la transformation
// inverse de Kawamura avec le gradient de temperature moyenne
void IJK_Thermique::calculer_temperature_physique_T(const IJK_Field_double&  vx, const double dTm)
{
  if(wall_flux_)
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

// TODO: supprimer cette methode qui n'est plus utile
// void IJK_Thermique::maj_Tl_Tv()
// {
//   const int nk = temperature_.nk();
//   const int ni = temperature_.ni();
//   const int nj = temperature_.nj();
//
//   // Calcul de moy1 = moy(chi*T+)/moy(chi) et moy2 = moy((1-chi)*T+) / moy(1-chi)
//   double moy1 = 0.;
//   double moy2 = 0.;
//   double Vl = 0.;
//   double Vv = 0.;
//
//   for (int k = 0; k < nk; k++)
//     for (int j = 0; j < nj; j++)
//       for (int i = 0; i < ni; i++)
//         {
//           const double chi = ref_ijk_ft_->indicatrice_ns_(i, j, k); // rappel : chi vaut 1. dans le liquide et 0 dans la vapeur
//           const double T = temperature_(i, j, k);
//           moy1 += T*chi;
//           moy2 += T*(1.-chi);
//           Vl += chi;
//           Vv += (1.-chi);
//         }
//   moy1 = Process::mp_sum(moy1);
//   moy2 = Process::mp_sum(moy2);
//   Vv = Process::mp_sum(Vv);
//   Vl = Process::mp_sum(Vl);
//   moy1 /= Vl;
//   moy2 /= Vv;
//
//   // Calcul de Tl et Tv :
//
//   const double Tl = Tv0_ / (1 - moy1) * (1 - moy1 / (1 + moy2)) / (1 + moy1 * moy2 / ((1 - moy1) * (1 + moy2)));
//   const double Tv = Tv0_ / (1 + moy2) * (1 + moy2 / (1 - moy1)) / (1 + moy2 * moy1 / ((1 + moy2) * (1 - moy1)));
//
//   // const double Tl = (Tl_ - Tv_) * moy1 + Tv0_;
//   // const double Tv = (Tl_ - Tv_) * moy2 + Tv0_;
//   Cerr << "Tl_test : " << Tl << endl;
//   Cerr << "Tv_test : " << Tv << endl;
//   Cerr << "Diviseur Tl_test : " << (1 + moy1 * moy2 / ((1 - moy1) * (1 + moy2))) << endl;
//   Cerr << "Diviseur Tv_test : " << (1 + moy2 * moy1 / ((1 + moy2) * (1 - moy1))) << endl;
//
//
//
//   // TODO: voir si on ne doit pas faire mieux, mais a priori les variations de Tl et Tv
//   // sont lentes par rapport au reste donc ea devrait aller.
//   // DONE: il y a manifestement un pb ici, car on ne peut pas avoir acces e Tv(n+1) encore,
//   // donc il faut stocker Tv(n-1)
//
//   const double dt = ref_ijk_ft_->timestep_;
//   dTv_ = (Tv - Tv_)/dt;
//   dTl_ = (Tl - Tl_)/dt;
//
//   Tv_ = Tv;
//   Tl_ = Tl;
// }

void IJK_Thermique::calculer_energies(double& E_liq_pure, double& E_lta, double& E_lth,
                                      double& E_vap_pure, double& E_vta, double& E_vth,
                                      double& E_mixt_arithm, double& E_mixt_harmo, double& E_tot, double& E_tot_h)
{
  const int nk = temperature_.nk();
  const int ni = temperature_.ni();
  const int nj = temperature_.nj();

  const IJK_Field_double& chi = ref_ijk_ft_->indicatrice_ns_; // rappel : chi vaut 1. dans le liquide et 0 dans la vapeur
  const double rhocpl = ref_ijk_ft_->rho_liquide_*cp_liquid_;
  const double rhocpv = ref_ijk_ft_->rho_vapeur_*cp_vapor_;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double indic=chi(i, j, k);
          const double rhocpm_harmo = 1./(indic/rhocpl+(1.-indic)/rhocpv);
          const double  rhocpm_arithm = indic*rhocpl+(1.-indic)*rhocpv;
          const double T = temperature_(i,j,k);
          const double Th = T*(rhocpm_harmo/rhocpm_arithm);
          E_tot += T*rhocpm_arithm;
          E_tot_h += Th*rhocpm_arithm;
          E_vta += (1.-indic)*rhocpv*T;
          E_lta += indic*rhocpl*T;
          E_vth += (1.-indic)*rhocpv*Th;
          E_lth += indic*rhocpl*Th;
          if (fabs(indic)<1.e-8)
            {
              // vap pure
              E_vap_pure +=rhocpv*T;
            }
          else if (fabs(1.-indic)<1.e-8)
            {
              // liq pure
              E_liq_pure +=rhocpl*T;
            }
          else
            {
              // mixte :
              E_mixt_arithm +=rhocpm_arithm*T;
              E_mixt_harmo +=rhocpm_harmo*T;
            }
        }

  const int ntot = temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_I)
                   *temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_J)
                   *temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  E_vap_pure = Process::mp_sum(E_vap_pure)/ntot;
  E_liq_pure = Process::mp_sum(E_liq_pure)/ntot;
  E_vta = Process::mp_sum(E_vta)/ntot;
  E_lta = Process::mp_sum(E_lta)/ntot;
  E_vth = Process::mp_sum(E_vth)/ntot;
  E_lth = Process::mp_sum(E_lth)/ntot;
  E_tot = Process::mp_sum(E_tot)/ntot;
  E_tot_h = Process::mp_sum(E_tot_h)/ntot;
  E_mixt_arithm = Process::mp_sum(E_mixt_arithm)/ntot;
  E_mixt_harmo = Process::mp_sum(E_mixt_harmo)/ntot;
}


void IJK_Thermique::calculer_temperature_adim_bulles()
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
  const IJK_Field_double& chi = ref_ijk_ft_->indicatrice_ns_; // rappel : chi vaut 1. dans le liquide et 0 dans la vapeur

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

  double E_tot_a = 0., E_tot_h = 0., E_liq_pure= 0., E_vap_pure= 0., E_mixt_arithm= 0., E_mixt_harmo= 0. ;
  double E_lta=0., E_lth=0., E_vta=0., E_vth =0.;
  calculer_energies(E_liq_pure, E_lta, E_lth,
                    E_vap_pure, E_vta, E_vth,
                    E_mixt_arithm, E_mixt_harmo, E_tot_a, E_tot_h);

  // TODO: voir si on ne doit pas faire mieux, mais a priori les variations de Tl et Tv
  // sont lentes par rapport au reste donc ea devrait aller.
  // DONE: il y a manifestement un pb ici, car on ne peut pas avoir acces e Tv(n+1) encore,
  // donc il faut stocker Tv(n-1)

  // Impression dans le fichier temperature_bulles.out
  if (Process::je_suis_maitre())
    {
      int reset = (!ref_ijk_ft_->reprise_) && (ref_ijk_ft_->tstep_==0);
      SFichier fic=Ouvrir_fichier(Nom("_source_temperature_")+Nom(rang_)+Nom("_bulles.out"),
                                  "tstep\ttime\tTl\tTv\tEtot_a\tEtot_h\tElpu\tElta\tElth\tEvpu\tEvta\tEvth\tEma\tEmh",
                                  reset);
      // la derivee_acceleration n'est connue que sur le maitre
      fic<< ref_ijk_ft_->tstep_<<" "<< ref_ijk_ft_->current_time_ <<" "<< Tl << " " << Tv << " " << E_tot_a << " " << E_tot_h ;
      fic<< " " << E_liq_pure <<  " " << E_lta <<  " " <<  E_lth;
      fic<< " " << E_vap_pure <<  " " << E_vta <<  " " <<  E_vth;
      fic<< " " << E_mixt_arithm <<  " " << E_mixt_harmo;
      fic<<finl;
      fic.close();
    }
}

void IJK_Thermique::calculer_temperature_physique_T_dummy()
{
  const int nk = d_temperature_.nk();
  const int ni = d_temperature_.ni();
  const int nj = d_temperature_.nj();

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              temperature_physique_T_(i,j,k) = temperature_(i,j,k);
            }
        }
    }
  temperature_physique_T_.echange_espace_virtuel(temperature_physique_T_.ghost());
  DebogIJK::verifier("temperature_physique_T", temperature_physique_T_);
  return;
}

void IJK_Thermique::calculer_temperature_adimensionnelle_theta(const IJK_Field_double&  vx, const double wall_flux)
{
  if(wall_flux_)
    {
      const IJK_Splitting& splitting = temperature_.get_splitting();
      const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
      const double Lz = geom.get_domain_length(DIRECTION_K);
      const double h = Lz/2.;
      const double q_w = wall_flux;
      const int kmin = temperature_.get_splitting().get_offset_local(DIRECTION_K);
      const int kmax = splitting.get_grid_geometry().get_nb_elem_tot(DIRECTION_K);
      const int nk = temperature_.nk();
      const int ni = temperature_.ni();
      const int nj = temperature_.nj();
      double T_wall = 0;
      T_wall = calculer_variable_wall(temperature_, cp_, ref_ijk_ft_->rho_field_, kmin, kmax);
      /*   if (Process::je_suis_maitre())
           {
           T_wall = calculer_variable_wall(temperature_, cp_, ref_ijk_ft_->rho_field_, kmin, kmax);
           Cerr << "calcul de T_wall sur maitre" << finl;
           }
           envoyer_broadcast(T_wall, 0); */
      /*
       * if(kmin+ nk == kmax) //|| (kmin ==0)
       {
       rank = Process::me();
       T_wall = calculer_variable_wall(temperature_, cp_, ref_ijk_ft_->rho_field_, kmin, kmax);
       envoyer_broadcast(T_wall, rank);
       }
       */
      for (int k = 0; k < nk; k++)
        {
          //  const double T_mean = compute_spatial_mean(vx, temperature_, cp_, ref_ijk_ft_->rho_field_, kmin, nktot, k);
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  double theta = T_wall - temperature_(i,j,k);
                  // double theta = T_wall -T_mean  ;
                  // temperature_adimensionnelle_theta_(i,j,k) = theta/theta_tau;
                  const double lambda_l = lambda_liquid_;
                  temperature_adimensionnelle_theta_(i,j,k) = theta*lambda_l/q_w/h;
                  //    temperature_adimensionnelle_theta_(i,j,k) = -T_mean*lambda_l/q_w/h;
                }
            }
        }
      return;
    }
  else
    {
      Cerr << "no_source_for_temperature" << finl;
      return;
    }

}


void IJK_Thermique::calculer_Nusselt(const IJK_Field_double& vx)
{
  const double theta_adim_moy = calculer_temperature_adimensionnelle_theta_moy(vx, temperature_adimensionnelle_theta_,
                                                                               cp_,ref_ijk_ft_->rho_field_);
  double Nu = 0.;
  if (fabs(theta_adim_moy)>1.e-10)
    Nu = 2./theta_adim_moy;
  const double rho_cp_u_moy = calculer_rho_cp_u_moyen(vx, cp_,ref_ijk_ft_->rho_field_);
  // Impression dans le fichier source_temperature.out
  if (Process::je_suis_maitre())
    {
      int reset = (!ref_ijk_ft_->reprise_) && (ref_ijk_ft_->tstep_==0);
      const Nom name = Nom("_temperature_")+Nom(rang_)+Nom(".out");
      SFichier fic=Ouvrir_fichier(name,
                                  "tstep\ttime\ttheta_adim_moy\tNu\trho_cp_u",
                                  reset);

      fic<< ref_ijk_ft_->tstep_ <<" "<< ref_ijk_ft_->current_time_ <<" "<< theta_adim_moy <<" "<< Nu << " " << rho_cp_u_moy  << finl ;
      fic.close();
    }
}

void IJK_Thermique::set_field_T_ana()
{
  Cerr << "Setting analytical temperature "<< rang_ <<" field to "<< expression_T_ana_ << finl;
  set_field_data(temperature_ana_, expression_T_ana_, ref_ijk_ft_->current_time_);
}

void IJK_Thermique::calculer_ecart_T_ana()
{
  if (liste_post_instantanes_.contient_("ECART_T_ANA"))
    {
      if (!liste_post_instantanes_.contient_("TEMPERATURE_ANA"))
        {
          set_field_data(temperature_ana_, expression_T_ana_, ref_ijk_ft_->current_time_);
        }
      // do some work

      double ct = ref_ijk_ft_->current_time_;
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
/*
   void IJK_Thermique::calculer_source_temperature_ana()
   {
   if (liste_post_instantanes_.contient_("ECART_SOURCE_TEMPERATURE_ANA"))
   {
   if (!liste_post_instantanes_.contient_("SOURCE_TEMPERATURE_ANA"))
   {
   set_field_data(source_temperature_ana_, expression_source_temperature_, ref_ijk_ft_->velocity_[0], ref_ijk_ft_->current_time_);
   }
// do some work

double ct = ref_ijk_ft_->current_time_;
Cerr << "MR: ERROR SOURCE T FIELD " << ct;
double err = 0.;
//set_field_data(source_temperature_ana_, curseur->expression_source_T_ana_, ct);
const int ni = source_temperature_.ni();
const int nj = source_temperature_.nj();
const int nk = source_temperature_.nk();
const int ntot=Process::mp_sum(ni*nj*nk);
// La temperature est definie a une constante pres:
//const double cst_temp = temperature_ana_(0,0,0) - curseur->temperature_(0,0,0);
for (int k = 0; k < nk; k++)
for (int j = 0; j < nj; j++)
for (int i = 0; i < ni; i++)
{
const double val =   source_temperature_ana_(i,j,k) - source_temperature_(i,j,k); //- cst_temp;
ecart_source_t_ana_(i,j,k) = val;
err += val*val;
}
err=Process::mp_sum(err);
err=sqrt(err/ntot);
Cerr << " " << err ;
if (!Process::je_suis_maitre())
{
Process::Journal() << "IJK_FT_Post::posttraiter_champs_instantanes : Champ ECART_SOURCE_TEMPERATURE_ANA sur ce proc (ni,nj,nk,ntot):"
<< " " << ni << " " << nj << " " << nk << " " << ntot << finl;
}
// ecart_t_ana_.echange_espace_virtuel(ecart_t_ana_.ghost());
Cerr << finl ;
//  n++,dumplata_scalar(lata_name,"ECART_SOURCE_TEMPERATURE_ANA", ecart_source_t_ana_, latastep);
}
} */

void IJK_Thermique::calculer_gradient_temperature(const IJK_Field_double& temperature, FixedVector<IJK_Field_double, 3>& grad_T)
{
  // Remise a zero :
  for (int dir = 0; dir < 3; dir++)
    {
      grad_T[dir].data() = 0.;
    }

  add_gradient_temperature(temperature, 1. /*constant*/ ,  grad_T[0], grad_T[1], grad_T[2],boundary_conditions_, lambda_);


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
void IJK_Thermique::compute_interfacial_temperature(ArrOfDouble& interfacial_temperature, ArrOfDouble& interfacial_phin_ai,
                                                    FixedVector<IJK_Field_double, 3> storage) const
{
  const IJK_Splitting& s =temperature_ft_.get_splitting();
  // Comme un cell_vector (ie cell-centered) mais en local :
  const int ni = temperature_ft_.ni();
  const int nj = temperature_ft_.nj();
  const int nk = temperature_ft_.nk();
  const IJK_Field_double& indic = ref_ijk_ft_->indicatrice_ft_;
  const IJK_Grid_Geometry& geom = s.get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  const double delta[3] = {dx, dy, dz};
  // "left" side of the "local" origin on the proc.
  double origin_x = geom.get_origin(DIRECTION_I)+s.get_offset_local(DIRECTION_I)*dx;
  double origin_y = geom.get_origin(DIRECTION_J)+s.get_offset_local(DIRECTION_J)*dy;
  double origin_z = geom.get_origin(DIRECTION_K)+s.get_offset_local(DIRECTION_K)*dz;
  double origin_local[3] = {origin_x,origin_y,origin_z};
  //int offsets[3] = {s.get_offset_local(DIRECTION_I), s.get_offset_local(DIRECTION_J), s.get_offset_local(DIRECTION_K)};
  const double lambda_l=DMINFLOAT+lambda_liquid_;
  const double lambda_v=DMINFLOAT+lambda_vapor_;
  for (int dir = 0; dir < 3; dir++)
    for (int k = 0; k < nk; k++)
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          {
            const double dTdx_twice = (temperature_ft_(i+1,j,k) - temperature_ft_(i-1,j,k))/dx;
            const double chi_x_l = indic(i,j,k)+0.5*(indic(i-1,j,k)+indic(i+1,j,k));
            const double chi_y_l = indic(i,j,k)+0.5*(indic(i,j-1,k)+indic(i,j+1,k));
            const double chi_z_l = indic(i,j,k)+0.5*(indic(i,j,k-1)+indic(i,j,k+1));
            const double lambda_eff_x_half =  1./(chi_x_l/(lambda_l) + (2-chi_x_l)/(lambda_v));
            const double qx = lambda_eff_x_half*dTdx_twice;
            const double qy = (temperature_ft_(i,j+1,k) - temperature_ft_(i,j-1,k))/(dy*(chi_y_l/lambda_l + (2-chi_y_l)/lambda_v));
            const double qz = (temperature_ft_(i,j,k+1) - temperature_ft_(i,j,k-1))/(dz*(chi_z_l/lambda_l + (2-chi_z_l)/lambda_v));
            storage[DIRECTION_I](i,j,k) = qx;
            storage[DIRECTION_J](i,j,k) = qy;
            storage[DIRECTION_K](i,j,k) = qz;
#if 0
// HACK TO CHECK that the flux is one
            const double x= origin_x+(i+0.5)*delta[0];
            const double y= origin_y+(j+0.5)*delta[1];
            const double z= origin_z+(k+0.5)*delta[2];
            const double r=sqrt(x*x+y*y+z*z);
            storage[DIRECTION_I](i,j,k) = x/r; // HACK!!!
            storage[DIRECTION_J](i,j,k) = y/r; // HACK!!!
            storage[DIRECTION_K](i,j,k) = z/r; // HACK!!!
#endif
          }
  storage[DIRECTION_I].echange_espace_virtuel(storage[DIRECTION_I].ghost());
  storage[DIRECTION_J].echange_espace_virtuel(storage[DIRECTION_J].ghost());
  storage[DIRECTION_K].echange_espace_virtuel(storage[DIRECTION_K].ghost());
  const Maillage_FT_IJK& maillage = ref_ijk_ft_->interfaces_.maillage_ft_ijk();
  Cerr << "ca passe" << finl;
  const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
  Cerr << "c'est passe" << finl;
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();
  const int nb_facettes=maillage.nb_facettes();
  const DoubleTab& sommets = maillage.sommets();
  const IntTab& facettes = maillage.facettes();

  // On prend T au centre de ijkb, et phin au centre de ijk... Mouais bof ou pas? En tout cas pas evident que ce soit mieux.
  Int3 ijkb;
  //
  // Computing interfacial flux at fa7 centre (weighted by ai)
  interfacial_phin_ai.resize_array(nb_facettes);
  interfacial_phin_ai = 0.;
  interfacial_temperature.resize_array(nb_facettes);
  interfacial_temperature = 0.;
  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    {
      // On rempli aussi les facettes_virtuelles
      //if (maillage.facette_virtuelle(fa7))
      //  continue;

      const double sf=surface_facettes[fa7];
      int index=intersections.index_facette()[fa7];
      while (index >= 0)
        {
          const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
          if (data.fraction_surface_intersection_< 1.e-6) // It's a relative value to the total surface.
            {
              index = data.index_element_suivant_;
              continue;
            }
          const int num_elem = data.numero_element_;
          // Attention, il faut bien que ce soit le splitting_ft car c'est comme cela que l'intersection a ete remplie :
          const Int3 ijk = s.convert_packed_to_ijk_cell(num_elem);
          const double surf = data.fraction_surface_intersection_ * sf;
          double phin = 0.;
          double dist=0.;
          ijkb[0]= 0.;
          ijkb[1]= 0.;
          ijkb[2]= 0.;
          for (int dir = 0; dir< 3; dir++)
            {
              const double ndir = normale_facettes(fa7,dir);
              // 0th order interpolation : ie we basically assumes that the flux computed at the cell centre is equal to the flux at the interface location.
              // Alternately, this means that flux is assumed constant within a cell.
              // It is not obvious whether an interpolation would be more accurate as it would rely on more temperatures in mixed cells :
              // void ijk_interpolate(const IJK_Field_double & field, const DoubleTab &coordinates, ArrOfDouble & result)
              //ijk_interpolate_skip_unknown_points(vitesse[direction], sommets, vinterp_component, 0.0 /* value for unknown points */);
              const double phidir = storage[dir](ijk[DIRECTION_I],ijk[DIRECTION_J],ijk[DIRECTION_K]);
              phin += ndir * phidir;
              //
              const int i0 = facettes(fa7,0);
              const int i1 = facettes(fa7,1);
              const int i2 = facettes(fa7,2);
              // Coordonnee du centre de gravite de la fraction de facette
              const double xG_frac = data.barycentre_[0]*sommets(i0,dir)+data.barycentre_[1]*sommets(i1,dir)+data.barycentre_[2]*sommets(i2,dir);
              const double x2 = (xG_frac - origin_local[dir])/delta[dir];
              ijkb[dir] = (int)(floor(x2)); // - offsets[dir];

              const double xc = (ijkb[dir]+0.5)*delta[dir];// centre de l'elem
              dist +=(xG_frac-xc)*ndir;
            }

          //
          // Value stored at cell center:
          const double Tc = temperature_ft_(ijkb[DIRECTION_I], ijkb[DIRECTION_J], ijkb[DIRECTION_K]);
          // Value stored at cell-center assumed equal to the value at cell-center... not perfect, but still...
          // Si dist>0, le centre de la cellule est dans la vapeur, donc le grad est du vapeur :
          const double Ti = (dist>0.) ? Tc + dist*phin/lambda_v  : Tc + dist*phin/lambda_l;
          interfacial_temperature[fa7] += Ti*surf;
          //
          interfacial_phin_ai[fa7] += phin*surf;
          index = data.index_element_suivant_;
        }
    }

  return;
}

double IJK_Thermique::compute_global_energy(const IJK_Field_double& temperature)
{
  global_energy_ = 0.;
  const IJK_Field_double& indic = ref_ijk_ft_->indicatrice_ns_;
  //const IJK_Grid_Geometry& geom = indic.get_splitting().get_grid_geometry();
  const double rhocpl = ref_ijk_ft_->rho_liquide_ *cp_liquid_;
  const double rhocpv = ref_ijk_ft_->rho_vapeur_ *cp_vapor_;
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
          //double cp = cp_(i,j,k) ;
          global_energy_ += (rhocpl * chi_l + (1.- chi_l) * rhocpv)*temperature(i,j,k);
        }
  const int ntot = temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_I)
                   *temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_J)
                   *temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  global_energy_ = mp_sum(global_energy_)/(double)(ntot);
  return global_energy_;
}

void IJK_Thermique::compute_T_rust(const FixedVector<IJK_Field_double, 3>& velocity)
{
  static Stat_Counter_Id cnt_conv_temp = statistiques().new_counter(1, "FT convection rho cp");
  statistiques().begin_count(cnt_conv_temp);
  // To be sure we're on a regular mesh
  assert(temperature_.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_K) >0);
  const double rho_l = ref_ijk_ft_->rho_liquide_;
  const double rho_v = ref_ijk_ft_->rho_vapeur_;

  // DONE: remplacer rho_cp par un champ rho_cp_ mis à jour dans update_thermal_properties. Necessaire pour que ça marche.
  //On calcule div(rho_cp*v) qu'on stocke dans T_rust
  switch (type_temperature_convection_op_)
    {
    case 1:
      temperature_convection_op_amont_.calculer(rho_cp_, rho_cp_, rho_cp_, velocity[0], velocity[1], velocity[2], T_rust_, T_rust_, T_rust_);
      break;
    case 2:
      temperature_convection_op_centre2_.calculer(rho_cp_, velocity[0], velocity[1], velocity[2], T_rust_);
      break;
    case 3:
      temperature_convection_op_quick_.calculer(rho_cp_, velocity[0], velocity[1], velocity[2], T_rust_);
      break;
    case 4:
      temperature_convection_op_centre4_.calculer(rho_cp_,rho_cp_,rho_cp_, velocity[0], velocity[1], velocity[2], T_rust_, T_rust_, T_rust_);
      break;

    default:
      Cerr << "Unknown convection operator for the temperature." << finl;
      Process::exit();
    }
  const int ni = T_rust_.ni();
  const int nj = T_rust_.nj();
  const int nk = T_rust_.nk();
  const IJK_Grid_Geometry& geom = T_rust_.get_splitting().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  const double vol = dx*dy*dz;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          T_rust_(i,j,k) /=vol ;
        }
  //On met a jour T_rust_ on le multipliant par T/rho_cp_hamro
  const int nx = temperature_.ni();
  const int ny = temperature_.nj();
  const int nz = temperature_.nk();
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        {
          const double chi_l = ref_ijk_ft_->indicatrice_ft_(i,j,k);
          T_rust_(i,j,k) *= temperature_(i,j,k)/(chi_l*rho_l*cp_liquid_ + (1-chi_l)*rho_v*cp_vapor_);
        }
  statistiques().end_count(cnt_conv_temp);
  return;
}

void IJK_Thermique::compute_dT_rustine(const double dE)
{
  const int ni = T_rust_.ni();
  const int nj = T_rust_.nj();
  const int nk = T_rust_.nk();
  const double rho_l = ref_ijk_ft_->rho_liquide_;
  const double rho_v = ref_ijk_ft_->rho_vapeur_;
  const IJK_Field_double indic = ref_ijk_ft_->indicatrice_ns_;
  double int_rhocpTrust = 0;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          int_rhocpTrust +=  (rho_l*cp_liquid_*indic(i,j,k) + rho_v*cp_vapor_*(1.-indic(i,j,k)))*T_rust_(i,j,k);
        }
  const int ntot = T_rust_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_I)
                   *T_rust_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_J)
                   *T_rust_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  int_rhocpTrust = mp_sum(int_rhocpTrust)/(double)(ntot);
  Cerr << "Le coeff de manque d'energie dE/int_rhocpTrust vaut : " << dE/int_rhocpTrust << finl;
  if (int_rhocpTrust)
    {
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              d_T_rustine_(i,j,k) = dE/ int_rhocpTrust * T_rust_(i,j,k);
            }
    }
}

void IJK_Thermique::euler_rustine_step(const double timestep, const double dE)
{
  compute_dT_rustine(dE);
  // Update the temperature :
  const int kmax = temperature_.nk();
  const double ene_ini = compute_global_energy();
  for (int k = 0; k < kmax; k++)
    {
      ref_ijk_ft_->euler_explicit_update(d_T_rustine_, temperature_, k);
    }
  temperature_.echange_espace_virtuel(temperature_.ghost());
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T"<<rang_<<" euler rustine] time t=" << ref_ijk_ft_->get_current_time()
       << " " << ene_ini
       << " " << ene_post << " [W.m-3]."
       << " dE "<< dE
       << finl;
  source_callback();
}

void IJK_Thermique::rk3_rustine_sub_step(const int rk_step, const double total_timestep,
                                         const double fractionnal_timestep, const double time, const double dE)
{
  compute_dT_rustine(dE);
  // Update the temperature :
  const int kmax = temperature_.nk();
  const double ene_ini = compute_global_energy();
  for (int k = 0; k < kmax; k++)
    {
      runge_kutta3_update(d_T_rustine_, RK3_F_rustine_, temperature_, rk_step, k, total_timestep);
    }
  temperature_.echange_espace_virtuel(temperature_.ghost());
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T"<<rang_<<"RK3 rustine step "<<rk_step<<"] time t=" << ref_ijk_ft_->get_current_time()
       << " " << ene_ini
       << " " << ene_post << " [W.m-3]. [step"<< rk_step << "]" << finl;
  source_callback();
}
#if 0
// Attention cette methode n'est utilisee nulle part !!!!!!!!
void IJK_Thermique::ecrire_reprise_thermique(SFichier& fichier )
{

  fichier << "   {\n";
  fichier << "     fichier_reprise_temperature " <<  fichier_reprise_temperature_ << "\n";
  fichier << "     timestep_reprise_temperature " <<  timestep_reprise_temperature_ << "\n";
  fichier << "     boundary_conditions " <<  boundary_conditions_ << "\n";
  fichier << "     lambda_liquid " <<  lambda_liquid_ << "\n";
  fichier << "     lambda_vapor " <<  lambda_vapor_ << "\n";
  fichier << "     cp_liquid " <<  cp_liquid_ << "\n";
  fichier << "     cp_vapor " <<  cp_vapor_ << "\n";
  fichier << "     type_T_source " << type_T_source_ << "\n";

  if ( expression_T_init_!="??")
    fichier << "     expression_T_init" <<  expression_T_init_ << "\n";

  if ( expression_T_ana_!="??")
    fichier << "     expression_T_ana" <<  expression_T_ana_ << "\n";

  if ( wall_flux_)
    fichier  << "     wall_flux \n";

  if ( conv_temperature_negligible_)
    fichier << "     conv_temperature_negligible \n ";

  if ( diff_temp_negligible_)
    fichier << "     diff_temp_negligible \n";

  fichier << "   }\n";
}
#endif
