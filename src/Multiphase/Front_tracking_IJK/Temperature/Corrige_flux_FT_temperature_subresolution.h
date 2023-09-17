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
// File      : Corrige_flux_FT_temperature_subresolution.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Corrige_flux_FT_temperature_subresolution_included
#define Corrige_flux_FT_temperature_subresolution_included

#include <Corrige_flux_FT_base.h>
#include <IJK_One_Dimensional_Subproblems.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Corrige_flux_FT_temperature_subresolution
//
// <Description of class Corrige_flux_FT_temperature_subresolution>
//
/////////////////////////////////////////////////////////////////////////////
#define FACES_DIR {0, 0, 1, 1, 2, 2}

class Corrige_flux_FT_temperature_subresolution : public Corrige_flux_FT_base
{

  Declare_instanciable( Corrige_flux_FT_temperature_subresolution ) ;

public :

  void initialize_with_subproblems(const IJK_Splitting& splitting,
                                   const IJK_Field_double& field,
                                   const IJK_Interfaces& interfaces,
                                   const IJK_FT_double& ijk_ft,
                                   Intersection_Interface_ijk_face& intersection_ijk_face,
                                   Intersection_Interface_ijk_cell& intersection_ijk_cell,
                                   const IJK_One_Dimensional_Subproblems& thermal_subproblems) override;

  void associate_thermal_problems(const IJK_One_Dimensional_Subproblems& thermal_subproblems);

  void set_convection_negligible(const int& convection_negligible) override { convection_negligible_ = convection_negligible; };
  void set_diffusion_negligible(const int& diffusion_negligible) override { diffusion_negligible_ = diffusion_negligible; };
  /*
   * On va calculer sur la grille IJ du layer k_layer tous les flux a proximite de
   * l'interface. On remplace les flux donnes en entree par ces flux la.
   */
  void corrige_flux_faceIJ(IJK_Field_local_double *const flux,
                           const int k_layer, const int dir) override { ; };

  void corrige_flux_conv_faceIJ(IJK_Field_local_double *const flux,
                                const int k_layer, const int dir) { corrige_flux_faceIJ(flux, k_layer, dir); };

  void corrige_flux_diff_faceIJ(IJK_Field_local_double *const flux,
                                const int k_layer, const int dir) override { ; };

  void calcul_temperature_flux_interface(const IJK_Field_double& temperature, const double ldal, const double ldav,
                                         const double dist, const DoubleTab& positions, const DoubleTab& normale,
                                         ArrOfDouble& temperature_interp, ArrOfDouble& flux_normal_interp,
                                         ArrOfDouble& temp_liqu, ArrOfDouble& temp_vap, DoubleTab& coo_liqu,
                                         DoubleTab& coo_vap) const override { ; };

  void update_intersections() override;
  void update() override;
  void associate_indices_and_check_subproblems_consistency();
  void compute_temperature_cell_centre(IJK_Field_double& temperature) const override;
  void set_zero_temperature_increment(IJK_Field_double& d_temperature) const override;
  void compute_temperature_face_centre() override;
  void compute_thermal_fluxes_face_centre() override;
  void clean() override;
  void check_pure_fluxes_duplicates(const DoubleVect& fluxes, DoubleVect& fluxes_unique, IntVect& pure_face_unique, const int known_unique);
protected :
  DoubleVect dist_;
  IntVect pure_face_unique_;
  DoubleVect convective_fluxes_;
  DoubleVect convective_fluxes_unique_;
  DoubleVect diffusive_fluxes_;
  DoubleVect diffusive_fluxes_unique_;
  const IJK_One_Dimensional_Subproblems * thermal_subproblems_;
  bool has_checked_consistency_;
  ArrOfInt ijk_intersections_subproblems_indices_;
  int convection_negligible_ = 0;
  int diffusion_negligible_ = 0;
  int debug_=0;
};

#endif /* Corrige_flux_FT_temperature_subresolution_included */
