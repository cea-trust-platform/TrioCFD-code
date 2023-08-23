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
// File      : IJK_Thermal_Subresolution.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Thermal_Subresolution_included
#define IJK_Thermal_Subresolution_included

#include <IJK_Thermal_base.h>
#include <IJK_Field.h>
#include <Boundary_Conditions_Thermique.h>
#include <IJK_Splitting.h>
#include <Parser.h>
#include <IJK_Lata_writer.h>
#include <OpConvQuickIJKScalar.h>
#include <OpConvCentre2IJKScalar.h>
#include <Ouvrir_fichier.h>
#include <Corrige_flux_FT.h>
#include <TRUST_Ref.h>
#include <Operateur_IJK_elem_diff_base.h>
#include <OpConvAmontIJK.h>
#include <OpConvDiscQuickIJKScalar.h>
#include <OpConvCentre4IJK.h>
#include <IJK_One_Dimensional_Subproblems.h>


/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_Thermal_Subresolution
//
// <Description of class IJK_Thermal_Subresolution>
//
/////////////////////////////////////////////////////////////////////////////


class IJK_Thermal_Subresolution : public IJK_Thermal_base
{

  Declare_instanciable( IJK_Thermal_Subresolution ) ;

public :

  int initialize(const IJK_Splitting& splitting, const int idx) override;
  void update_thermal_properties() override;
  void set_param(Param& param) override;
  void compute_ghost_cell_numbers_for_subproblems(const IJK_Splitting& splitting, int ghost_init) override;

protected :

  void compute_diffusion_increment() override;
  void correct_temperature_for_eulerian_fluxes() override;
  void correct_temperature_for_visu() override;
  void compute_overall_probes_parameters() override;
  void compute_radial_convection_diffusion_operators(DoubleTab& radial_first_order_operator_raw,
                                                     DoubleTab& radial_second_order_operator_raw,
                                                     DoubleTab& radial_first_order_operator,
                                                     DoubleTab& radial_second_order_operator,
                                                     DoubleTab& radial_diffusion_matrix,
                                                     DoubleTab& radial_convection_matrix);
  void compute_first_order_operator_raw(DoubleTab& radial_first_order_operator);
  void compute_first_order_operator(DoubleTab& radial_first_order_operator, double dr);
  void compute_second_order_operator(DoubleTab& radial_second_order_operator, double dr);
  void compute_second_order_operator_raw(DoubleTab& radial_second_order_operator);
  void compute_radial_convection_operator(const DoubleTab& radial_first_order_operator,
                                          DoubleTab& radial_convection_matrix);
  void compute_radial_diffusion_operator(const DoubleTab& radial_second_order_operator,
                                         DoubleTab& radial_diffusion_matrix);
  void initialise_thermal_subproblems() override;
  void solve_thermal_subproblems() override;
  void apply_thermal_flux_correction() override;
  void clean_thermal_subproblems() override;
  /* compute_rho_cp_u_mean() May be clearly overridden later */
  double compute_rho_cp_u_mean(const IJK_Field_double& vx) override { return IJK_Thermal_base::compute_rho_cp_u_mean(vx); };

  int diffusive_flux_correction_;
  int convective_flux_correction_;

  int override_vapour_mixed_values_; // For debug purposes

  IJK_One_Dimensional_Subproblems thermal_local_subproblems_;
  int points_per_thermal_subproblem_;
  double coeff_distance_diagonal_ = 3.;
  double probe_length_;
  double dr_;
  DoubleVect radial_coordinates_;
  DoubleTab radial_first_order_operator_raw_;
  DoubleTab radial_second_order_operator_raw_;
  DoubleTab radial_first_order_operator_;
  DoubleTab radial_second_order_operator_;
  DoubleTab radial_diffusion_matrix_;
  DoubleTab radial_convection_matrix_;

};

#endif /* IJK_Thermal_Subresolution_included */
