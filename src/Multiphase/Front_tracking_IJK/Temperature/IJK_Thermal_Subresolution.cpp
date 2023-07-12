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
// File      : IJK_Thermal_Subresolution.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Thermal_Subresolution.h>
#include <Param.h>
#include <IJK_Navier_Stokes_tools.h>
#include <DebogIJK.h>
#include <stat_counters.h>
#include <IJK_FT.h>
#include <Corrige_flux_FT.h>
#include <OpConvDiscIJKQuickScalar.h>

Implemente_instanciable_sans_constructeur( IJK_Thermal_Subresolution, "IJK_Thermal_Subresolution", IJK_Thermal_base ) ;

IJK_Thermal_Subresolution::IJK_Thermal_Subresolution()
{
  convective_flux_correction_ = 0;
  diffusion_flux_correction_ = 0;
  ghost_fluid_ = 0;
  override_vapour_mixed_values_ = 0;
}

Sortie& IJK_Thermal_Subresolution::printOn( Sortie& os ) const
{
  IJK_Thermal_base::printOn( os );
  return os;
}

Entree& IJK_Thermal_Subresolution::readOn( Entree& is )
{
  IJK_Thermal_base::readOn( is );
  Cout << "IJK_Thermal_Subresolution::readOn : Parameters summary. " << finl;
  printOn(Cout);
  return is;
}

void IJK_Thermal_Subresolution::set_param( Param& param )
{
  IJK_Thermal_base::set_param(param);
  param.ajouter_flag("convective_flux_correction", &convective_flux_correction_);
  param.ajouter_flag("diffusion_flux_correction", &diffusion_flux_correction_);
  param.ajouter_flag("ghost_fluid", &ghost_fluid_);
  param.ajouter_flag("override_vapour_mixed_values", &override_vapour_mixed_values_);
}

int IJK_Thermal_Subresolution::initialize(const IJK_Splitting& splitting, const int idx)
{
  Cout << que_suis_je() << "::initialize()" << finl;

  uniform_lambda_ = lambda_liquid_;
  uniform_alpha_ =	lambda_liquid_ / (ref_ijk_ft_->get_rho_l() * cp_liquid_);
  calulate_grad_T_ = 1;
  // TODO: Implement ghost fluid if necessary
  ghost_fluid_ = 0;

  int nalloc = 0;
  nalloc = IJK_Thermal_base::initialize(splitting, idx);

  /*TODO:
   * Change the operators to add fluxes corrections
   */
  if (diffusion_flux_correction_)
    {
      temperature_diffusion_op_.typer("OpDiffUniformIJKScalarCorrection_double");
    }
  if (convective_flux_correction_)
    {
      temperature_convection_op_.typer("OpConvIJKQuickScalar_double");
    }

  temperature_convection_op_.initialize(splitting);
  temperature_diffusion_op_.initialize(splitting);

  Cout << "End of " << que_suis_je() << "::initialize()" << finl;
  return nalloc;
}

void IJK_Thermal_Subresolution::update_thermal_properties()
{
  IJK_Thermal_base::update_thermal_properties();
}

void IJK_Thermal_Subresolution::compute_diffusion_increment()
{
  /*
   * Update d_temperature
   * d_temperature_ += div_lambda_grad_T_volume_;
   * It depends on the nature of the properties (one-fluid or single-fluid)
   */
  const int ni = d_temperature_.ni();
  const int nj = d_temperature_.nj();
  const int nk = d_temperature_.nk();
  const IJK_Grid_Geometry& geom = d_temperature_.get_splitting().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  const double vol = dx*dy*dz;
  const double rhocp_l = ref_ijk_ft_->get_rho_l() * cp_liquid_;
//  const double rhocp_v = ref_ijk_ft_->get_rho_v() * cp_vapour_;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          double rhocpVol;
          rhocpVol = rhocp_l * vol;
          const double ope = div_coeff_grad_T_volume_(i,j,k);
          const double resu = ope/rhocpVol;
          d_temperature_(i,j,k) += resu;
        }
}


void IJK_Thermal_Subresolution::correct_temperature_for_eulerian_fluxes()
{
  if (!ghost_fluid_)
    {
      if (override_vapour_mixed_values_)
        {
          const int ni = temperature_.ni();
          const int nj = temperature_.nj();
          const int nk = temperature_.nk();
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  const double indic = ref_ijk_ft_->itfce().I(i,j,k);
                  if (std::fabs(1.-indic)>1.e-8) // Mixed cells and pure vapour cells
                    { temperature_(i,j,k) = 0; }
                }
        }
// TODO: Temperature sub-resolution (Weak coupling)
//			else
//				{
//
//				}
    }
  else
    {

    }
}

/*
 * Correct directly the fluxes...
 */
//void IJK_Thermal_Subresolution::correct_temperature_increment_diffusion()
//{
//  if (!ghost_fluid_)
//    {
//      if (!diffusion_flux_correction_)
//        {
//          if (override_vapour_mixed_values_)
//            {
//              const int ni = div_coeff_grad_T_volume_.ni();
//              const int nj = div_coeff_grad_T_volume_.nj();
//              const int nk = div_coeff_grad_T_volume_.nk();
//              for (int k = 0; k < nk; k++)
//                for (int j = 0; j < nj; j++)
//                  for (int i = 0; i < ni; i++)
//                    {
//                      const double indic = ref_ijk_ft_->itfce().I(i,j,k);
//                      if (std::fabs(1.-indic)>1.e-8) // Mixed cells and pure vapour cells
//                        { div_coeff_grad_T_volume_(i,j,k) = 0; }
//                    }
//            }
//          else
//            {
//
//            }
//        }
//      // TODO: Temperature sub-resolution (Strong coupling)
//			//			else
//			//				{
//			//
//			//				}
//    }
//  else
//		{
//
//		}
//}

