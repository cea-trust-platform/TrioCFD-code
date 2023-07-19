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

Implemente_instanciable( IJK_Thermal_Onefluid, "IJK_Thermal_Onefluid", IJK_Thermal_base ) ;

Sortie& IJK_Thermal_Onefluid::printOn( Sortie& os ) const
{
  IJK_Thermal_base::printOn( os );
  os<< "  {\n";

  os<< "    type_T_source " << type_T_source_ << "\n";

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

  os<< "  \n}";
  return os;
}

Entree& IJK_Thermal_Onefluid::readOn( Entree& is )
{
  IJK_Thermal_base::readOn( is );

  lambda_moy_arith_=0;
  type_temperature_convection_form_ = 0;  // Default value: 0 : non conservative
  conserv_energy_global_=0;

  Param param(que_suis_je());

  param.ajouter("type_temperature_convection_form", &type_temperature_convection_form_);
  param.dictionnaire("non conservative",0);
  param.dictionnaire("conservative",1);
  param.ajouter("conserv_energy_global", &conserv_energy_global_);

  param.lire_avec_accolades(is);
  Cout << "IJK_Thermal_Onefluid::readOn : Parameters summary. " << finl;
  IJK_Thermal_Onefluid::printOn(Cout);

  return is;
}

int IJK_Thermal_Onefluid::initialize(const IJK_Splitting& splitting, const int idx)
{
  Cout << que_suis_je() << "::initialize()" << finl;
  int nalloc = 0;
  nalloc = IJK_Thermal_base::initialize(splitting, idx);
  temperature_diffusion_op_.set_conductivity_coefficient(uniform_lambda_, lambda_, temperature_, temperature_, temperature_);
  cp_.allocate(splitting, IJK_Splitting::ELEM, 2);
  lambda_.allocate(splitting, IJK_Splitting::ELEM, 1);
  nalloc += 2;

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

void IJK_Thermal_Onefluid::update_thermal_properties()
{
  IJK_Thermal_base::update_thermal_properties();
  const IJK_Field_double& indic = ref_ijk_ft_->itfce().I();
  // Nombre de mailles du domaine NS :
  const int nx = indic.ni();
  const int ny = indic.nj();
  const int nz = indic.nk();
//  const double rho_l = ref_ijk_ft_->get_rho_l();
//  const double rho_v = ref_ijk_ft_->get_rho_v();
  const bool geometric_mean = ((!lambda_moy_arith_) and (lambda_liquid_ > DMINFLOAT) and (lambda_vapour_ >DMINFLOAT));
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        {
          double chi_l = indic(i,j,k);
          if (geometric_mean)
            lambda_(i,j,k) = lambda_liquid_ * lambda_vapour_ / ((1 - chi_l) * lambda_liquid_ + chi_l * lambda_vapour_);
          else
            lambda_(i,j,k) = lambda_liquid_  * chi_l + (1.- chi_l) * lambda_vapour_ ;
        }
  lambda_.echange_espace_virtuel(lambda_.ghost());
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
  const double rhocp_l = ref_ijk_ft_->get_rho_l()*cp_liquid_;
  const double rhocp_v = ref_ijk_ft_->get_rho_v()*cp_vapour_;
  const bool geometric_mean = ((rho_cp_inv_) and (rhocp_l > DMINFLOAT) and (rhocp_v >DMINFLOAT));
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          if (geometric_mean)
            {
              const double chi_l = ref_ijk_ft_->itfce().I(i,j,k);
              const double rhocpV_inv = (chi_l/rhocp_l  + (1 - chi_l)/rhocp_v) * vol_inv;
              const double ope = div_coeff_grad_T_volume_(i,j,k);
              const double resu = ope*rhocpV_inv;
              d_temperature_(i,j,k) +=resu ;
            }
          else
            {
              double rhocpV = 0;
              if (depracated_rho_cp_)
                {
                  const double rho = ref_ijk_ft_->get_rho_field_ijk(i,j,k);
                  const double cp = cp_(i,j,k);
                  rhocpV = rho * cp * vol;
                }
              else
                {
                  const double chi_l = ref_ijk_ft_->itfce().I(i,j,k);
                  rhocpV = (rhocp_l * chi_l + rhocp_v * (1 - chi_l)) *  vol;
                }
              const double ope = div_coeff_grad_T_volume_(i,j,k);
              const double resu = ope/rhocpV;
              d_temperature_(i,j,k) +=resu ;
            }
        }
}
