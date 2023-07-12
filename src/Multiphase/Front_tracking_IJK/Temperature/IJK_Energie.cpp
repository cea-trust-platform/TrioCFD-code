/****************************************************************************
 * Copyright (c) 2019, CEA
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 *modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *this list of conditions and the following disclaimer in the documentation
 *and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *may be used to endorse or promote products derived from this software without
 *specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *****************************************************************************/
/////////////////////////////////////////////////////////////////////////////
//
// File      : IJK_Energie.cpp
// Directory : $IJK_ROOT/src/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <DebogIJK.h>
#include <IJK_Energie.h>
#include <IJK_FT.h>
#include <IJK_Navier_Stokes_tools.h>
#include <Param.h>
#include <stat_counters.h>

Implemente_instanciable(IJK_Energie, "IJK_Energie", Objet_U);

Sortie& IJK_Energie::printOn(Sortie& os) const
{
  Objet_U::printOn(os);
  os << "  {\n"
     << "    boundary_conditions {"
     << "\n";
  Nom bctype_kmin, bctype_kmax, bckmin, bckmax, valeur_kmin, valeur_kmax;
  if (boundary_conditions_.get_bctype_k_max() == boundary_conditions_.Perio)
    {
      bctype_kmax = "Perio";
      bctype_kmin = "Perio";
      bckmax = " ";
      bckmin = " ";
      valeur_kmax = " ";
    }
  else
    Cerr << "La seule condition aux limites validée pour la formulation en "
         "énergie est la condition Perio"
         << finl;

  if (boundary_conditions_.get_bctype_k_min() == boundary_conditions_.Perio)
    {
      bctype_kmin = "Perio";
      bctype_kmin = "Perio";
      bckmin = " ";
      bckmin = " ";
      valeur_kmin = " ";
    }
  else
    Cerr << "La seule condition aux limites validée pour la formulation en "
         "énergie est la condition Perio"
         << finl;

  os << "      bctype_kmin"
     << " " << bctype_kmin << " \n";
  os << "      bctype_kmax"
     << " " << bctype_kmax << " \n";

  os << "      " << bckmin << " " << valeur_kmin << " \n";
  os << "      " << bckmax << " " << valeur_kmax << " \n";
  os << "    } \n";

  os << "    lambda_liquid " << lambda_liquid_ << "\n";
  os << "    lambda_vapor " << lambda_vapor_ << "\n";
  os << "    cp_liquid " << cp_liquid_ << "\n";
  os << "    cp_vapor " << cp_vapor_ << "\n";

  os << "    fichier_reprise_temperature"
     << " " << fichier_reprise_temperature_ << "\n";
  os << "    timestep_reprise_temperature"
     << " " << timestep_reprise_temperature_ << "\n";

  if (expression_T_ana_ != "??")
    os << "    expression_T_ana"
       << " " << expression_T_ana_ << "\n";
  if (conv_temperature_negligible_)
    os << "    conv_temperature_negligible \n ";
  if (diff_temp_negligible_)
    os << "    diff_temp_negligible \n";
  os << "  }\n";
  return os;
}

Entree& IJK_Energie::readOn(Entree& is)
{
  rang_ = 0; // default value
  cp_liquid_ = 1.;
  cp_vapor_ = 1.;
  lambda_liquid_ = 1.;
  lambda_vapor_ = 1.;
  fo_ = 1.;
  expression_T_init_ = "??";
  fichier_reprise_temperature_ = "??";
  expression_T_ana_ = "??";

  // reprise_ = 0;
  timestep_reprise_temperature_ = 1;
  diff_temp_negligible_ = 0;
  conv_temperature_negligible_ = 0;
  // terme source variable

  // wall_flux_=0.;
  //
  // Ceci est une initialisation des derivees des temperatures moyenne de chaque
  // phase Il n'est peut-etre pas pertinent de les mettre ici

  Param param(que_suis_je());

  param.ajouter("cp_liquid", &cp_liquid_, Param::REQUIRED);
  param.ajouter("lambda_liquid", &lambda_liquid_, Param::REQUIRED);
  param.ajouter("cp_vapor", &cp_vapor_, Param::REQUIRED);
  param.ajouter("lambda_vapor", &lambda_vapor_, Param::REQUIRED);
  param.ajouter("fo", &fo_);
  param.ajouter("boundary_conditions", &boundary_conditions_, Param::REQUIRED);

  param.ajouter("expression_T_init", &expression_T_init_);
  param.ajouter("fichier_reprise_temperature", &fichier_reprise_temperature_);
  param.ajouter("timestep_reprise_temperature", &timestep_reprise_temperature_);
  param.ajouter_flag("conv_temperature_negligible",
                     &conv_temperature_negligible_);
  param.ajouter_flag("diff_temp_negligible", &diff_temp_negligible_);

  // Expression analytique de la temperature
  param.ajouter("expression_T_ana", &expression_T_ana_);

  param.lire_avec_accolades(is);
  Cout << "IJK_Energie::readOn : Parameters summary. " << finl;
  printOn(Cout);

  return is;
}

int IJK_Energie::initialize(const IJK_Splitting& splitting, const int idx)
{
  Cout << que_suis_je() << "::initialize()" << finl;
  rang_ = idx;

  cp_.initialize(cp_liquid_, cp_vapor_);
  lda_.initialize(lambda_liquid_, lambda_vapor_);
  rho_.initialize(ref_ijk_ft_->get_rho_l(), ref_ijk_ft_->get_rho_v());
  rho_cp_ = rho_ * cp_;

  int nalloc = 0;

  temperature_.allocate(splitting, IJK_Splitting::ELEM, 2);
  lambda_.allocate(splitting, IJK_Splitting::ELEM, 1);
  d_temperature_.allocate(splitting, IJK_Splitting::ELEM, 2); // 1
  div_lambda_grad_T_volume_.allocate(splitting, IJK_Splitting::ELEM, 0);
  D_rhocp_T_.allocate(splitting, IJK_Splitting::ELEM, 0);
  nalloc += 5;

  diffusion_temperature_op_.initialize(splitting);

  corrige_flux_.typer("Corrige_flux_FT_temperature_conv");

  corrige_flux_.set_physical_parameters(rho_cp_.liqu(), rho_cp_.vap(),
                                        lda_.liqu(), lda_.vap());
  corrige_flux_.initialize(
    ref_ijk_ft_->get_splitting_ns(),
    temperature_,
    ref_ijk_ft_->itfce(),
    ref_ijk_ft_);

  energy_convection_op_quick_interface_.set_corrige_flux(corrige_flux_);
  energy_convection_op_quick_interface_.initialize(splitting);


  if ((ref_ijk_ft_.non_nul()) and (!ref_ijk_ft_->disable_diphasique()))
    {
      temperature_ft_.allocate(ref_ijk_ft_->get_splitting_ft(),
                               IJK_Splitting::ELEM, 1);
      nalloc += 4;
    }

  if (fichier_reprise_temperature_ ==
      "??") // si on ne fait pas une reprise on initialise V
    {
      if (expression_T_init_ != "??")
        {
          Cout << "Temperature initialization from expression \nTini = "
               << expression_T_init_ << finl;
          set_field_data(temperature_, expression_T_init_,
                         ref_ijk_ft_->itfce().I(), 0.);
        }
      else
        {
          Cerr << "Please provide initial conditions for temperature, either by an "
               "expression or a field for restart. "
               << "You should consider using either fichier_reprise_temperature or "
               "expression_T_init keywords. "
               << finl;
          Process::exit();
        }
    }
  else
    {
      Cout << "Reading initial temperature field T" << rang_ << " from file "
           << fichier_reprise_temperature_
           << " timestep= " << timestep_reprise_temperature_ << finl;
      const Nom& geom_name = splitting.get_grid_geometry().le_nom();
      lire_dans_lata(
        fichier_reprise_temperature_, timestep_reprise_temperature_, geom_name,
        Nom("TEMPERATURE_") + Nom(idx),
        temperature_); // fonction qui lit un champ a partir d'un lata .
      temperature_.echange_espace_virtuel(
        temperature_.ghost()); // It is essential to fill the EV because the
      // first call to convection needs them.
    }

  Cerr << " Initializing thermal fields dependant on the post-pro list : ";
  if (liste_post_instantanes_.size())
    Cerr << liste_post_instantanes_;
  else
    Cerr << "empty";
  Cerr << finl;
  if (liste_post_instantanes_.size() &&
      (liste_post_instantanes_.contient_("TEMPERATURE_ANA") ||
       liste_post_instantanes_.contient_("ECART_T_ANA")))
    {
      temperature_ana_.allocate(splitting, IJK_Splitting::ELEM, 1);
      ecart_t_ana_.allocate(splitting, IJK_Splitting::ELEM, 1);
      nalloc += 2;
    }

  if ((liste_post_instantanes_.size() &&
       liste_post_instantanes_.contient_("GRAD_T")) ||
      (ref_ijk_ft_.non_nul() && ref_ijk_ft_->t_debut_statistiques() < 1.e10))
    allocate_velocity(grad_T_, splitting, 1);
  nalloc += 3;

  // rho_cp_.allocate(splitting, IJK_Splitting::ELEM, 2);
  // nalloc +=1;

  rho_cp_T_.allocate(splitting, IJK_Splitting::ELEM, 2);
  div_rho_cp_T_.allocate(splitting, IJK_Splitting::ELEM, 0);
  nalloc += 2;

  Cout << "End of " << que_suis_je() << "::initialize()" << finl;
  return nalloc;
}

void IJK_Energie::associer(const IJK_FT_double& ijk_ft)
{
  ref_ijk_ft_ = ijk_ft;
  liste_post_instantanes_ = ijk_ft.get_post().get_liste_post_instantanes();
}

void IJK_Energie::sauvegarder_temperature(Nom& lata_name, int idx)
{
  fichier_reprise_temperature_ = lata_name;
  timestep_reprise_temperature_ = 1;
  dumplata_scalar(lata_name, Nom("TEMPERATURE_") + Nom(idx), temperature_,
                  0 /*we store a 0 */);
}

void IJK_Energie::update_thermal_properties()
{
  const double ene_ini = compute_global_energy();
  lda_(ref_ijk_ft_->itfce().I(), lambda_);
  // Nombre de mailles du domaine NS :
  const int nx = lambda_.ni();
  const int ny = lambda_.nj();
  const int nz = lambda_.nk();
  double chi_l;
  for (int k = 0; k < nz; k++)
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
        {
          chi_l = ref_ijk_ft_->itfce().I(i, j, k);
          rho_cp_T_(i, j, k) = rho_cp_(chi_l) * temperature_(i, j, k);
        }
  rho_cp_T_.echange_espace_virtuel(rho_cp_T_.ghost());

  // Semble un endroit approprie pour calculer la variation d'energie due au
  // transport de l'interface:
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T" << rang_
       << "-2-TransportIndic] time t=" << ref_ijk_ft_->get_current_time() << " "
       << ene_ini << " " << ene_post << " delta=" << ene_post - ene_ini
       << " [W.m-3]." << finl;
}

// Methode de calcul du pas de temps max base sur Fo pour l'equation de
// thermique. CFL value is not computed as it is the same as for the velocity
// equation. The calculation should be stable if Fo <= 1.0 (thanks to the 0.5 in
// the formula below).
double IJK_Energie::compute_timestep(const double timestep,
                                     const double dxmin) const
{
  // alpha = lambda/(rho*cp)
  const double alpha_max = std::max(lambda_liquid_ / (rho_.liqu() * cp_liquid_),
                                    lambda_vapor_ / (rho_.vap() * cp_vapor_));
  // Attention 0.125 vient du 3D. (1/6 au lieu de 1/8)
  const double dt_fo = dxmin * dxmin / (alpha_max + 1.e-20) * fo_ * (1. / 6.);
  return dt_fo;
}

void IJK_Energie::euler_time_step(
  const FixedVector<IJK_Field_double, 3>& velocity)
{
  calculer_dT(velocity);
  // Update the temperature :
  const int kmax = temperature_.nk();
  const double ene_ini = compute_global_energy();
  Cerr << "d_temperature(0,0,0) = " << d_temperature_(0, 0, 0) << finl;
  for (int k = 0; k < kmax; k++)
    {
      ref_ijk_ft_->euler_explicit_update(d_temperature_, temperature_, k);
    }
  temperature_.echange_espace_virtuel(temperature_.ghost());
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T" << rang_
       << "] time t=" << ref_ijk_ft_->get_current_time() << " " << ene_ini
       << " " << ene_post << " [W.m-3]." << finl;
}

// Mettre rk_step = -1 si schema temps different de rk3.
void IJK_Energie::calculer_dT(
  const FixedVector<IJK_Field_double, 3>& velocity)
{
  const double current_time = ref_ijk_ft_->get_current_time();
  const double ene_ini = compute_global_energy(d_temperature_);
  compute_energy_convection(velocity);
  const double ene_postConv = compute_global_energy(d_temperature_);
  add_temperature_diffusion();
  const double ene_postDiffu = compute_global_energy(d_temperature_);
  const double ene_postSource = compute_global_energy(d_temperature_);
  add_temporal_rho_cp_term();
  divide_by_rho_cp_np1();

  Cerr << "[Energy-Budget-T" << rang_
       << "-1-TimeResolution] time t=" << current_time << " " << ene_ini << " "
       << ene_postConv << " " << ene_postDiffu << " " << ene_postSource
       << " delta=" << ene_postSource - ene_ini << " [W.m-3]." << finl;

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
          Tmax = std::max(Tmax, T(i, j, k));
          Tmin = std::min(Tmin, T(i, j, k));
        }
  Tmax = Process::mp_max(Tmax);
  Tmin = Process::mp_min(Tmin);
  Cerr << "[Temperature-MinMax-" << rang_ << "] t/Tmin/Tmax " << current_time
       << " " << Tmin << " " << Tmax << finl;

  return;
}

// Convect energy field by velocity.
// The output is stored in d_temperature_ (it is a volume integral over the CV)
void IJK_Energie::compute_energy_convection(
  const FixedVector<IJK_Field_double, 3>& velocity)
{
  static Stat_Counter_Id cnt_conv_temp =
    statistiques().new_counter(1, "FT convection rho");
  statistiques().begin_count(cnt_conv_temp);
  if (conv_temperature_negligible_)
    {
      d_temperature_.data() = 0;
    }
  else
    {
      const int ni = d_temperature_.ni();
      const int nj = d_temperature_.nj();
      const int nk = d_temperature_.nk();
      energy_convection_op_quick_interface_.calculer(
        temperature_, velocity[0], velocity[1], velocity[2], d_temperature_);
      const IJK_Grid_Geometry& geom =
        d_temperature_.get_splitting().get_grid_geometry();
      const double dx = geom.get_constant_delta(DIRECTION_I);
      const double dy = geom.get_constant_delta(DIRECTION_J);
      const double dz = geom.get_constant_delta(DIRECTION_K);
      const double vol = dx * dy * dz;
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              d_temperature_(i, j, k) /= vol;
            }
    }
  statistiques().end_count(cnt_conv_temp);
  DebogIJK::verifier("op_conv(rho)", d_temperature_);
  return;
}

// Cette méthode ajoute 1/dV int_{dS}{lambda_h grad_T dS} à d_temperature_
void IJK_Energie::add_temperature_diffusion()
{
  if (boundary_conditions_.get_bctype_k_min() ==
      Boundary_Conditions_Thermique::Paroi_Temperature_imposee)
    {
      const double T_paroi_impose_kmin =
        boundary_conditions_.get_temperature_kmin();
      double lambda_de_t_paroi_kmin = lambda_liquid_;
      // calculer ici div(lambda*grad(T))*volume)
      calculer_flux_thermique_bord(temperature_, lambda_de_t_paroi_kmin,
                                   T_paroi_impose_kmin, boundary_flux_kmin_,
                                   0 /* boundary kmin */);
    }
  else if (boundary_conditions_.get_bctype_k_min() ==
           Boundary_Conditions_Thermique::Paroi_Flux_impose)
    {
      const double flux_paroi_impose_kmin = boundary_conditions_.get_flux_kmin();
      imposer_flux_thermique_bord(temperature_, flux_paroi_impose_kmin,
                                  boundary_flux_kmin_, 0 /* boundary kmin */);
    }
  else
    {
      // Perio... not needed. The flux on the "boundary" (in that case , it's not
      // truely a boundary) will be computed as inside...
    }
  if (boundary_conditions_.get_bctype_k_max() ==
      Boundary_Conditions_Thermique::Paroi_Temperature_imposee)
    {
      const double T_paroi_impose_kmax =
        boundary_conditions_.get_temperature_kmax();
      double lambda_de_t_paroi_kmax = lambda_liquid_;
      // calculer ici div(lambda*grad(T))*volume)
      calculer_flux_thermique_bord(temperature_, lambda_de_t_paroi_kmax,
                                   T_paroi_impose_kmax, boundary_flux_kmax_,
                                   1 /* boundary kmax */);
    }
  else if (boundary_conditions_.get_bctype_k_max() ==
           Boundary_Conditions_Thermique::Paroi_Flux_impose)
    {
      const double flux_paroi_impose_kmax = boundary_conditions_.get_flux_kmax();
      imposer_flux_thermique_bord(temperature_, flux_paroi_impose_kmax,
                                  boundary_flux_kmax_, 1 /* boundary kmax */);
      // Cerr << "not coded yet" << finl;
      // Process::exit();
    }
  else
    {
      // Perio... not needed. The flux on the "boundary" (in that case , it's not
      // truely a boundary) will be computed as inside...
    }
  temperature_.echange_espace_virtuel(temperature_.ghost());
  lambda_.echange_espace_virtuel(lambda_.ghost());
  DebogIJK::verifier("temp", temperature_);
  DebogIJK::verifier("lambda", lambda_);

  // Performance counters:
  static Stat_Counter_Id cnt_diff_temp =
    statistiques().new_counter(1, "FT diffusion temperature");
  statistiques().begin_count(cnt_diff_temp);
  if (diff_temp_negligible_) // si diffusion negligeable
    {
      div_lambda_grad_T_volume_.data() = 0;
    }
  else
    {
      diffusion_temperature_op_.set_lambda(lambda_);
      diffusion_temperature_op_.calculer(
        temperature_, div_lambda_grad_T_volume_, boundary_flux_kmin_,
        boundary_flux_kmax_);
      // Update d_temperature
      // d_temperature_ +=div_lambda_grad_T_volume_;

      const int ni = d_temperature_.ni();
      const int nj = d_temperature_.nj();
      const int nk = d_temperature_.nk();
      const IJK_Grid_Geometry& geom =
        d_temperature_.get_splitting().get_grid_geometry();
      const double dx = geom.get_constant_delta(DIRECTION_I);
      const double dy = geom.get_constant_delta(DIRECTION_J);
      const double dz = geom.get_constant_delta(DIRECTION_K);
      const double vol = dx * dy * dz;
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double ope = div_lambda_grad_T_volume_(i, j, k);
              const double resu = ope / vol;
              d_temperature_(i, j, k) += resu;
            }
    }

  statistiques().end_count(cnt_diff_temp);
  DebogIJK::verifier("div_lambda_grad_T_volume", div_lambda_grad_T_volume_);
  // Cerr << "diff_temp" << " " <<  d_temperature_(1,1,1) << finl;
  // TODO: deplacer ce bloc ailleurs...
  if (liste_post_instantanes_.contient_("GRAD_T") ||
      (ref_ijk_ft_->t_debut_statistiques() < 1.e10))
    {
      calculer_gradient_temperature(temperature_, grad_T_);
    }
}

void IJK_Energie::add_temporal_rho_cp_term()
{
  if (boundary_conditions_.get_bctype_k_min() ==
      Boundary_Conditions_Thermique::Paroi_Temperature_imposee)
    {
      throw "La formulation en energie n'est pas implemente pour une temperature "
      "impose";
    }
  else if (boundary_conditions_.get_bctype_k_min() ==
           Boundary_Conditions_Thermique::Paroi_Flux_impose)
    {
      throw "La formulation en energie n'est pas implemente pour des flux impose";
    }
  else
    {
      // Perio... not needed. The flux on the "boundary" (in that case , it's not
      // truely a boundary) will be computed as inside...
    }
  temperature_.echange_espace_virtuel(temperature_.ghost());
  DebogIJK::verifier("temp", temperature_);

  // Performance counters:
  static Stat_Counter_Id cnt_rhocp_temp =
    statistiques().new_counter(1, "FT temporal rhocp temperature");
  statistiques().begin_count(cnt_rhocp_temp);
  if (conv_temperature_negligible_) // si diffusion negligeable
    {
      D_rhocp_T_.data() = 0;
    }
  else
    {
      // D_rhocp_T_ = - T^n * (rhocp^(n+1) - rhocp^(n))
      const int ni = D_rhocp_T_.ni();
      const int nj = D_rhocp_T_.nj();
      const int nk = D_rhocp_T_.nk();
      for (int i = 0; i < ni; i++)
        for (int j = 0; j < nj; j++)
          for (int k = 0; k < nk; k++)
            {
              const double i_np1 = ref_ijk_ft_->itfce().In(i, j, k);
              const double i_n = ref_ijk_ft_->itfce().I(i, j, k);
              D_rhocp_T_(i, j, k) =
                -temperature_(i, j, k) * (rho_cp_(i_np1) - rho_cp_(i_n));
            }
    }
  // Update d_temperature

  const int ni = d_temperature_.ni();
  const int nj = d_temperature_.nj();
  const int nk = d_temperature_.nk();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double resu = D_rhocp_T_(i, j, k);
          d_temperature_(i, j, k) += resu;
        }

  statistiques().end_count(cnt_rhocp_temp);
  DebogIJK::verifier("D_rhocp_T", D_rhocp_T_);
}

void IJK_Energie::divide_by_rho_cp_np1()
{
  if (boundary_conditions_.get_bctype_k_min() ==
      Boundary_Conditions_Thermique::Paroi_Temperature_imposee)
    {
      throw "La formulation en energie n'est pas implemente pour une temperature "
      "impose";
    }
  else if (boundary_conditions_.get_bctype_k_min() ==
           Boundary_Conditions_Thermique::Paroi_Flux_impose)
    {
      throw "La formulation en energie n'est pas implemente pour des flux impose";
    }
  else
    {
      // Perio... not needed. The flux on the "boundary" (in that case , it's not
      // truely a boundary) will be computed as inside...
    }
  temperature_.echange_espace_virtuel(temperature_.ghost());
  DebogIJK::verifier("temp", temperature_);

  // Performance counters:
  static Stat_Counter_Id cnt_div_rhocp_temp =
    statistiques().new_counter(1, "FT div rhocp np1 temperature");
  statistiques().begin_count(cnt_div_rhocp_temp);
  // Update d_temperature

  const int ni = d_temperature_.ni();
  const int nj = d_temperature_.nj();
  const int nk = d_temperature_.nk();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double rhocp_np1 =rho_cp_(ref_ijk_ft_->itfce().In(i,j,k));
          d_temperature_(i, j, k) /= rhocp_np1;
        }
  statistiques().end_count(cnt_div_rhocp_temp);
}

void IJK_Energie::calculer_energies(double& E_liq_pure, double& E_lta,
                                    double& E_lth, double& E_vap_pure,
                                    double& E_vta, double& E_vth,
                                    double& E_mixt_arithm, double& E_mixt_harmo,
                                    double& E_tot, double& E_tot_h)
{
  const int nk = temperature_.nk();
  const int ni = temperature_.ni();
  const int nj = temperature_.nj();

  const IJK_Field_double& chi =
    ref_ijk_ft_->itfce().I(); // rappel : chi vaut 1. dans le liquide et 0
  // dans la vapeur
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double indic = chi(i, j, k);
          const double rhocpm_harmo = rho_cp_.harmo(indic);
          const double rhocpm_arithm = rho_cp_(indic);
          const double T = temperature_(i, j, k);
          const double Th = T * (rhocpm_harmo / rhocpm_arithm);
          E_tot += T * rhocpm_arithm;
          E_tot_h += Th * rhocpm_arithm;
          E_vta += (1. - indic) * rho_cp_.vap() * T;
          E_lta += indic * rho_cp_.liqu() * T;
          E_vth += (1. - indic) * rho_cp_.vap() * Th;
          E_lth += indic * rho_cp_.liqu() * Th;
          if (std::abs(indic) < 1.e-12)
            {
              // vap pure
              E_vap_pure += rho_cp_.vap() * T;
            }
          else if (std::abs(1. - indic) < 1.e-12)
            {
              // liq pure
              E_liq_pure += rho_cp_.liqu() * T;
            }
          else
            {
              // mixte :
              E_mixt_arithm += rhocpm_arithm * T;
              E_mixt_harmo += rhocpm_harmo * T;
            }
        }

  const int ntot =
    temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM,
                                                     DIRECTION_I) *
    temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM,
                                                     DIRECTION_J) *
    temperature_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM,
                                                     DIRECTION_K);
  E_vap_pure = Process::mp_sum(E_vap_pure) / ntot;
  E_liq_pure = Process::mp_sum(E_liq_pure) / ntot;
  E_vta = Process::mp_sum(E_vta) / ntot;
  E_lta = Process::mp_sum(E_lta) / ntot;
  E_vth = Process::mp_sum(E_vth) / ntot;
  E_lth = Process::mp_sum(E_lth) / ntot;
  E_tot = Process::mp_sum(E_tot) / ntot;
  E_tot_h = Process::mp_sum(E_tot_h) / ntot;
  E_mixt_arithm = Process::mp_sum(E_mixt_arithm) / ntot;
  E_mixt_harmo = Process::mp_sum(E_mixt_harmo) / ntot;
}

void IJK_Energie::set_field_T_ana()
{
  Cerr << "Setting analytical temperature " << rang_ << " field to "
       << expression_T_ana_ << finl;
  set_field_data(temperature_ana_, expression_T_ana_,
                 ref_ijk_ft_->get_current_time());
}

void IJK_Energie::calculer_ecart_T_ana()
{
  if (liste_post_instantanes_.contient_("ECART_T_ANA"))
    {
      if (!liste_post_instantanes_.contient_("TEMPERATURE_ANA"))
        {
          set_field_data(temperature_ana_, expression_T_ana_,
                         ref_ijk_ft_->get_current_time());
        }
      // do some work

      double ct = ref_ijk_ft_->get_current_time();
      Cerr << "GB: ERROR T FIELD " << ct;
      double err = 0.;
      set_field_data(temperature_ana_, expression_T_ana_, ct);
      const int ni = temperature_.ni();
      const int nj = temperature_.nj();
      const int nk = temperature_.nk();
      const int ntot = Process::mp_sum(ni * nj * nk);
      // La temperature est definie a une constante pres:
      // const double cst_temp = temperature_ana_(0,0,0) -
      // curseur->temperature_(0,0,0);
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double val =
                temperature_ana_(i, j, k) - temperature_(i, j, k); //- cst_temp;
              ecart_t_ana_(i, j, k) = val;
              err += val * val;
            }
      err = Process::mp_sum(err);
      err = sqrt(err / ntot);
      Cerr << " " << err;
      if (!Process::je_suis_maitre())
        {
          Process::Journal() << "IJK_FT_Post::posttraiter_champs_instantanes : "
                             "Champ ECART_T_ANA sur ce proc (ni,nj,nk,ntot):"
                             << " " << ni << " " << nj << " " << nk << " " << ntot
                             << finl;
        }
      ecart_t_ana_.echange_espace_virtuel(ecart_t_ana_.ghost());
      Cerr << finl;
      //  n++,dumplata_scalar(lata_name,"ECART_T_ANA", ecart_t_ana_, latastep);
    }
}

void IJK_Energie::calculer_gradient_temperature(
  const IJK_Field_double& temperature,
  FixedVector<IJK_Field_double, 3>& grad_T)
{
  // Remise a zero :
  for (int dir = 0; dir < 3; dir++)
    {
      grad_T[dir].data() = 0.;
    }

  add_gradient_temperature(temperature, 1. /*constant*/, grad_T[0], grad_T[1],
                           grad_T[2], boundary_conditions_, lambda_);

  for (int dir = 0; dir < 3; dir++)
    {
      grad_T[dir].echange_espace_virtuel(1);
    }
}

// Results are intensive (ie prop to area)
// Les interfaces connaissent le splitting_ft_ donc la correspondance doit etre
// appliquee au splitting ft pour convertir : convert_packed_to_ijk_cell. Donc
// il faut un champ de T etendu...
void IJK_Energie::compute_interfacial_temperature2(
  ArrOfDouble& interfacial_temperature,
  ArrOfDouble& flux_normal_interp) const
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_geometry();
  const double dist =
    1.52 * std::pow(std::pow(geom.get_constant_delta(0), 2.) +
                    std::pow(geom.get_constant_delta(1), 2.) +
                    std::pow(geom.get_constant_delta(2), 2.),
                    0.5);
  const Maillage_FT_IJK& maillage = ref_ijk_ft_->itfce().maillage_ft_ijk();
  const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const int nb_facettes = maillage.nb_facettes();
  const IntTab& facettes = maillage.facettes();
  const DoubleTab& sommets = maillage.sommets();
  DoubleTab coord_facettes;
  coord_facettes.resize(nb_facettes, 3);
  coord_facettes = 0.;
  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      for (int som = 0; som < 3; som++)
        coord_facettes(fa7, dir) += sommets(facettes(fa7, som), dir);

  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      coord_facettes(fa7, dir) /= 3.;

  DoubleTab coo_liqu, coo_vap;
  ArrOfDouble temp_liqu, temp_vap;
  temp_liqu.set_smart_resize(1);
  temp_vap.set_smart_resize(1);
  coo_liqu.set_smart_resize(1);
  coo_vap.set_smart_resize(1);
//  Corrige_flux_FT_temperature_conv::calcul_temperature_flux_interface(
//    temperature_ft_, lambda_liquid_, lambda_vapor_, dist, coord_facettes,
//    normale_facettes, interfacial_temperature, flux_normal_interp, temp_liqu,
//    temp_vap, coo_liqu, coo_vap);
  corrige_flux_.calcul_temperature_flux_interface(temperature_ft_,
                                                  lambda_liquid_,
                                                  lambda_vapor_,
                                                  dist,
                                                  coord_facettes,
                                                  normale_facettes,
                                                  interfacial_temperature,
                                                  flux_normal_interp,
                                                  temp_liqu,
                                                  temp_vap,
                                                  coo_liqu,
                                                  coo_vap);
  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    {
      flux_normal_interp(fa7) *= surface_facettes(fa7);
      interfacial_temperature(fa7) *= surface_facettes(fa7);
    }
}

double IJK_Energie::compute_global_energy(const IJK_Field_double& temperature)
{
  global_energy_ = 0.;
  const IJK_Field_double& indic = ref_ijk_ft_->itfce().I();
  const int nx = temperature.ni();
  const int ny = temperature.nj();
  const int nz = temperature.nk();
  // To be sure we're on a regular mesh
  assert(indic.get_splitting().get_grid_geometry().get_constant_delta(
           DIRECTION_K) > 0);
  for (int k = 0; k < nz; k++)
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
        {
          double chi_l = indic(i, j, k);
          // double cp = cp_(i,j,k) ;
          global_energy_ += rho_cp_(chi_l) * temperature(i, j, k);
        }
  const int ntot =
    temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM,
                                                    DIRECTION_I) *
    temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM,
                                                    DIRECTION_J) *
    temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM,
                                                    DIRECTION_K);
  global_energy_ = mp_sum(global_energy_) / (double)(ntot);
  return global_energy_;
}

/*
 * Methods that do not belong to the class
 */

// From DNS_QC; Vectorize code later?
int IJK_Energie::calculer_k_pour_bord(const IJK_Field_double& temperature, const bool bord_kmax)
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
int IJK_Energie::calculer_flux_thermique_bord(const IJK_Field_double& temperature,
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
int IJK_Energie::imposer_flux_thermique_bord(const IJK_Field_double& temperature,
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

double IJK_Energie::get_rhocp_l() const { return rho_cp_.liqu(); }

double IJK_Energie::get_rhocp_v() const { return rho_cp_.vap(); }

double IJK_Energie::get_lda_l() const { return lambda_liquid_; }

double IJK_Energie::get_lda_v() const { return lambda_vapor_; }
