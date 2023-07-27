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
// File      : IJK_Thermal_Onefluid.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Thermal_Onefluid.h>
#include <IJK_FT.h>
#include <DebogIJK.h>
#include <IJK_Navier_Stokes_tools.h>

Implemente_instanciable_sans_constructeur( IJK_Thermal_Onefluid, "IJK_Thermal_Onefluid", IJK_Thermal_base ) ;

IJK_Thermal_Onefluid::IJK_Thermal_Onefluid()
{
  single_phase_=0;
  lambda_moy_arith_=0;
  type_temperature_convection_form_ = 0;  // Default value: 0 : non conservative
  conserv_energy_global_=0;
  rho_cp_moy_harmonic_=0;
  rho_cp_post_=0;
  E0_=0;
  deprecated_rho_cp_=0;
}

Sortie& IJK_Thermal_Onefluid::printOn( Sortie& os ) const
{
  IJK_Thermal_base::printOn( os );
  os<< "  {\n";

  os<< "    type_T_source " << type_T_source_ << "\n";

  if (rho_cp_moy_harmonic_)
    os<< "    rho_cp_moy_harmonic \n";
  if (lambda_variable_)
    os<< "    lambda_variable \n";
  if (lambda_moy_arith_)
    os<< "    lambda_moy_arith \n";
  if (deprecated_rho_cp_)
    os<< "    depracated_rho_cp \n";
  if (conserv_energy_global_)
    os<< "    conserv_energy_global \n";

  os<< "  \n}";
  return os;
}

Entree& IJK_Thermal_Onefluid::readOn( Entree& is )
{
  IJK_Thermal_base::readOn( is );
  return is;
}

void IJK_Thermal_Onefluid::set_param( Param& param )
{
  IJK_Thermal_base::set_param(param);
  param.ajouter("type_temperature_convection_form", &type_temperature_convection_form_);
  param.dictionnaire("non conservative",0);
  param.dictionnaire("conservative",1);
  param.ajouter("conserv_energy_global", &conserv_energy_global_);
  param.ajouter_flag("lambda_moy_arith_", &lambda_moy_arith_);
  param.ajouter_flag("rho_cp_moy_harmonic", &rho_cp_moy_harmonic_);
  param.ajouter_flag("deprecated_rho_cp", &deprecated_rho_cp_);
}

int IJK_Thermal_Onefluid::initialize(const IJK_Splitting& splitting, const int idx)
{
  Cout << que_suis_je() << "::initialize()" << finl;
  int nalloc = 0;
  nalloc = IJK_Thermal_base::initialize(splitting, idx);
  temperature_diffusion_op_.set_conductivity_coefficient(uniform_lambda_, lambda_, temperature_, temperature_, temperature_);
  lambda_.allocate(splitting, IJK_Splitting::ELEM, 1);
  nalloc += 2;

  if (rho_cp_moy_harmonic_)
    {
      rho_cp_inv_.allocate(splitting, IJK_Splitting::ELEM, 2);
      nalloc += 1;
    }
  else
    {
      if (deprecated_rho_cp_)
        {
          cp_.allocate(splitting, IJK_Splitting::ELEM, 2);
          nalloc += 1;
        }
      else
        {
          if (!rho_cp_post_)
            {
              rho_cp_.allocate(splitting, IJK_Splitting::ELEM, 2);
              nalloc += 1;
            }
        }
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
  if (type_temperature_convection_form_==1)
    {
      rho_cp_T_.allocate(splitting, IJK_Splitting::ELEM, 2);
      div_rho_cp_T_.allocate(splitting, IJK_Splitting::ELEM, 0);
      nalloc += 2;
    }
  Cout << "End of " << que_suis_je() << "::initialize()" << finl;
  return nalloc;
}

void IJK_Thermal_Onefluid::update_thermal_properties()
{
  IJK_Thermal_base::update_thermal_properties();
  const IJK_Field_double& indic = ref_ijk_ft_->itfce().I();
  // Nombre de mailles du domaine NS :
  const int nx = indic.ni();
  const int ny = indic.nj();
  const int nz = indic.nk();
  const double rho_l = ref_ijk_ft_->get_rho_l();
  const double rho_v = ref_ijk_ft_->get_rho_v();
  const bool harmonic_mean = ((!lambda_moy_arith_) and (lambda_liquid_ > DMINFLOAT) and (lambda_vapour_ >DMINFLOAT));
  const bool rho_cp_harmonic_mean = ((rho_cp_moy_harmonic_) and (rho_l*cp_liquid_ > DMINFLOAT) and (rho_v*cp_vapour_ >DMINFLOAT));
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        {
          double chi_l = indic(i,j,k);
          if (rho_cp_moy_harmonic_)
            if (rho_cp_harmonic_mean)
              rho_cp_inv_(i,j,k)= chi_l / (rho_l * cp_liquid_)  + (1 - chi_l) / (rho_v * cp_vapour_);
            else
              rho_cp_inv_(i,j,k) = rho_l * cp_liquid_ * chi_l + rho_v * cp_vapour_ * (1.- chi_l);
          else
            {
              if (deprecated_rho_cp_)
                cp_(i,j,k) = cp_liquid_ * chi_l + cp_vapour_ * (1.- chi_l);
              else if (!rho_cp_post_)
                rho_cp_(i,j,k) = rho_l * cp_liquid_ * chi_l + rho_v * cp_vapour_ * (1.- chi_l);
            }
          if (rho_cp_post_)
            rho_cp_(i,j,k) = rho_l * cp_liquid_ * chi_l + rho_v * cp_vapour_ * (1.- chi_l);
          if (harmonic_mean)
            lambda_(i,j,k) = lambda_liquid_ * lambda_vapour_ / ((1 - chi_l) * lambda_liquid_ + chi_l * lambda_vapour_);
          else
            lambda_(i,j,k) = lambda_liquid_  * chi_l + (1.- chi_l) * lambda_vapour_ ;
        }
  lambda_.echange_espace_virtuel(lambda_.ghost());
  if (rho_cp_post_)
    rho_cp_.echange_espace_virtuel(rho_cp_.ghost());
  if (rho_cp_moy_harmonic_)
    rho_cp_inv_.echange_espace_virtuel(rho_cp_inv_.ghost());
  else
    {
      if (deprecated_rho_cp_)
        cp_.echange_espace_virtuel(cp_.ghost());
      else
        {
          if (!rho_cp_post_)
            rho_cp_.echange_espace_virtuel(rho_cp_.ghost());
        }
    }
}

void IJK_Thermal_Onefluid::add_temperature_diffusion()
{
  lambda_.echange_espace_virtuel(lambda_.ghost());
  DebogIJK::verifier("lambda", lambda_);

  temperature_diffusion_op_.set_lambda(lambda_);
  IJK_Thermal_base::add_temperature_diffusion();
}

void IJK_Thermal_Onefluid::compute_diffusion_increment()
{
  // Update d_temperature
  const int ni = d_temperature_.ni();
  const int nj = d_temperature_.nj();
  const int nk = d_temperature_.nk();
  const double vol_inv = 1./vol_;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          if (rho_cp_moy_harmonic_)
            {
              const double rhocpV_inv = rho_cp_inv_(i,j,k) * vol_inv;
              const double ope = div_coeff_grad_T_volume_(i,j,k);
              const double resu = ope*rhocpV_inv;
              d_temperature_(i,j,k) +=resu ;
            }
          else
            {
              double rhocpV = 0;
              if (deprecated_rho_cp_)
                {
                  const double rho = ref_ijk_ft_->get_rho_field_ijk(i,j,k);
                  const double cp = cp_(i,j,k);
                  rhocpV = rho * cp * vol_;
                }
              else
                {
                  rhocpV = rho_cp_(i,j,k) * vol_;
                }
              const double ope = div_coeff_grad_T_volume_(i,j,k);
              const double resu = ope/rhocpV;
              d_temperature_(i,j,k) +=resu ;
            }
        }
}

double IJK_Thermal_Onefluid::compute_rho_cp_u_mean(const IJK_Field_double& vx)
{
  double rho_cp_u_mean = 0.;
  if (rho_cp_moy_harmonic_)
    {
      const double rho_l = ref_ijk_ft_->get_rho_l();
      const double rho_v = ref_ijk_ft_->get_rho_v();
      const bool rho_cp_harmonic_mean = ((rho_cp_moy_harmonic_) and (rho_l*cp_liquid_ > DMINFLOAT) and (rho_v*cp_vapour_ >DMINFLOAT));
      if (rho_cp_harmonic_mean)
        rho_cp_u_mean = calculer_rho_cp_u_moyen_inv(vx, rho_cp_inv_);
      else
        rho_cp_u_mean = calculer_rho_cp_u_moyen(vx, rho_cp_inv_);
    }
  else
    {
      if (deprecated_rho_cp_)
        {
          rho_cp_u_mean = calculer_rho_cp_u_moyen(vx, cp_, ref_ijk_ft_->get_rho_field());
        }
      else
        {
          rho_cp_u_mean = calculer_rho_cp_u_moyen(vx, rho_cp_);
        }
    }
  return rho_cp_u_mean;
}

double IJK_Thermal_Onefluid::get_rho_cp_ijk(int i, int j, int k) const
{
  double rho_cp = 0.;
  if (rho_cp_moy_harmonic_)
    {
      const double rho_l = ref_ijk_ft_->get_rho_l();
      const double rho_v = ref_ijk_ft_->get_rho_v();
      const bool rho_cp_harmonic_mean = ((rho_cp_moy_harmonic_) and (rho_l*cp_liquid_ > DMINFLOAT) and (rho_v*cp_vapour_ >DMINFLOAT));
      if (rho_cp_harmonic_mean)
        rho_cp = 1 / rho_cp_inv_(i,j,k);
      else
        rho_cp = rho_cp_inv_(i,j,k);
    }
  else
    {
      if (deprecated_rho_cp_)
        {
          rho_cp = cp_(i,j,k) * (ref_ijk_ft_->get_rho_field_ijk(i,j,k));
        }
      else
        {
          rho_cp = rho_cp_(i,j,k);
        }
    }
  return rho_cp;
}

double IJK_Thermal_Onefluid::get_rho_cp_u_ijk(const IJK_Field_double& vx, int i, int j, int k) const
{
  return get_rho_cp_ijk(i,j,k) * vx(i,j,k);
}

double IJK_Thermal_Onefluid::get_div_lambda_ijk(int i, int j, int k) const
{
  return (lambda_(i+1,j,k)-lambda_(i-1,j,k));
}

double IJK_Thermal_Onefluid::compute_temperature_dimensionless_theta_mean(const IJK_Field_double& vx)
{
  double theta_adim_moy = 0.;
  if (rho_cp_moy_harmonic_)
    {
      const double rho_l = ref_ijk_ft_->get_rho_l();
      const double rho_v = ref_ijk_ft_->get_rho_v();
      const bool rho_cp_harmonic_mean = ((rho_cp_moy_harmonic_) and (rho_l*cp_liquid_ > DMINFLOAT) and (rho_v*cp_vapour_ >DMINFLOAT));
      if (rho_cp_harmonic_mean)
        theta_adim_moy = calculer_temperature_adimensionnelle_theta_moy_inv(vx, temperature_adimensionnelle_theta_, rho_cp_inv_);
      else
        theta_adim_moy = calculer_temperature_adimensionnelle_theta_moy(vx, temperature_adimensionnelle_theta_, rho_cp_inv_); //TEST
    }
  else
    {
      if (deprecated_rho_cp_)
        {
          theta_adim_moy = calculer_temperature_adimensionnelle_theta_moy(vx, temperature_adimensionnelle_theta_, cp_, ref_ijk_ft_->get_rho_field());
        }
      else
        {
          theta_adim_moy = calculer_temperature_adimensionnelle_theta_moy(vx, temperature_adimensionnelle_theta_, rho_cp_);
        }
    }
  return theta_adim_moy;
}
