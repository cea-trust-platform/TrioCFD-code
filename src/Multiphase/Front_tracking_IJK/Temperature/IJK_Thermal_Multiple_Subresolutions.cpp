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
// File      : IJK_Thermal_Multiple_Subresolutions.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Thermal_Multiple_Subresolutions.h>
#include <Param.h>
#include <IJK_Navier_Stokes_tools.h>
#include <DebogIJK.h>
#include <stat_counters.h>
#include <IJK_FT.h>
#include <Corrige_flux_FT.h>
#include <OpConvDiscIJKQuickScalar.h>

Implemente_instanciable( IJK_Thermal_Multiple_Subresolutions, "IJK_Thermal_Multiple_Subresolutions", IJK_Thermal_Subresolution ) ;

Sortie& IJK_Thermal_Multiple_Subresolutions::printOn( Sortie& os ) const
{
  IJK_Thermal_Subresolution::printOn( os );
  return os;
}

Entree& IJK_Thermal_Multiple_Subresolutions::readOn( Entree& is )
{
  IJK_Thermal_Subresolution::readOn( is );
  main_phase_ = 0;
  uniform_lambda_vap_=0;
  uniform_alpha_vap_=0;
  /*
   * Parse the datafile
   */
  Param param(que_suis_je());
  param.ajouter_flag("main_phase_", &main_phase_);
  param.lire_avec_accolades(is);

  /*
   * The main phase is called liquid because it is mainly used in bubbly flows
   * It will directly influence which phase uses the temperature_ field attribute.
   * Correction of temperature_ and temperature_vap_ may be different
   */
  if (main_phase_)
    {
      uniform_lambda_vap_ = lambda_vapour_;
      uniform_alpha_vap_ = lambda_vapour_ / (ref_ijk_ft_->get_rho_v() * cp_vapour_);
    }
  else
    {
      uniform_lambda_vap_ = lambda_liquid_;
      uniform_lambda_ = lambda_vapour_;
      uniform_alpha_vap_ = lambda_liquid_ / (ref_ijk_ft_->get_rho_l() * cp_liquid_);
      uniform_alpha_ = lambda_vapour_ / (ref_ijk_ft_->get_rho_v() * cp_vapour_);
    }
  return is;
}

int IJK_Thermal_Multiple_Subresolutions::initialize(const IJK_Splitting& splitting, const int idx)
{
  Cout << que_suis_je() << "::initialize()" << finl;

  int nalloc = 0;
  nalloc = IJK_Thermal_Subresolution::initialize(splitting, idx);

  temperature_vapour_.allocate(splitting, IJK_Splitting::ELEM, 2);
  nalloc += 1;

  Cout << "End of " << que_suis_je() << "::initialize()" << finl;
  return nalloc;
}

void IJK_Thermal_Multiple_Subresolutions::update_thermal_properties()
{
  IJK_Thermal_Subresolution::update_thermal_properties();
}

void IJK_Thermal_Multiple_Subresolutions::correct_temperature_vapour_for_eulerian_fluxes()
{
  if (!ghost_fluid_)
    {
      if (override_vapour_mixed_values_)
        {
          const int ni = temperature_vapour_.ni();
          const int nj = temperature_vapour_.nj();
          const int nk = temperature_vapour_.nk();
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  const double indic = ref_ijk_ft_->itfce().I(i,j,k);
                  if (std::fabs(1.-indic)>1.e-8) // Mixed cells and pure vapour cells
                    { temperature_vapour_(i,j,k) = 0; }
                }
        }
// TODO: Temperature sub-resolution
//			else
//				{
//
//				}
    }
}

void IJK_Thermal_Multiple_Subresolutions::correct_temperature_vapour_increment_diffusion()
{
  if (!ghost_fluid_)
    {
      if (!diffusion_flux_correction_)
        {
          if (override_vapour_mixed_values_)
            {
              const int ni = div_coeff_grad_T_vapour_volume_.ni();
              const int nj = div_coeff_grad_T_vapour_volume_.nj();
              const int nk = div_coeff_grad_T_vapour_volume_.nk();
              for (int k = 0; k < nk; k++)
                for (int j = 0; j < nj; j++)
                  for (int i = 0; i < ni; i++)
                    {
                      const double indic = ref_ijk_ft_->itfce().I(i,j,k);
                      if (std::fabs(indic)<(1.-1.e-8)) // Mixed cells and pure liquid cells
                        { div_coeff_grad_T_vapour_volume_(i,j,k) = 0; }
                    }
            }
          else
            {

            }
        }
      // TODO: Temperature sub-resolution (Strong coupling)
//			else
//				{
//
//				}
    }
}
