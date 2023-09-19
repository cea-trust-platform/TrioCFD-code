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
// File      : Corrige_flux_FT.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Corrige_flux_FT_included
#define Corrige_flux_FT_included

#include <Corrige_flux_FT_base.h>
#include <IJK_One_Dimensional_Subproblems.h>
#include <TRUST_Deriv.h>
/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Corrige_flux_FT
//
// <Description of class Corrige_flux_FT>
//
/////////////////////////////////////////////////////////////////////////////

class Corrige_flux_FT : public DERIV(Corrige_flux_FT_base)
{

  Declare_instanciable( Corrige_flux_FT ) ;

public :

  inline void initialize(const IJK_Splitting& splitting,
                         const IJK_Field_double& field,
                         const IJK_Interfaces& interfaces,
                         const IJK_FT_double& ijk_ft,
                         Intersection_Interface_ijk_face& intersection_ijk_face,
                         Intersection_Interface_ijk_cell& intersection_ijk_cell);
  inline void initialize_with_subproblems(const IJK_Splitting& splitting,
                                          const IJK_Field_double& field,
                                          const IJK_Interfaces& interfaces,
                                          const IJK_FT_double& ijk_ft,
                                          Intersection_Interface_ijk_face& intersection_ijk_face,
                                          Intersection_Interface_ijk_cell& intersection_ijk_cell,
                                          const IJK_One_Dimensional_Subproblems& thermal_local_subproblems);
  inline void set_physical_parameters(const double rhocpl,
                                      const double rhocpv,
                                      const double ldal,
                                      const double ldav);
  inline void corrige_flux_faceIJ(IJK_Field_local_double *const flux,
                                  const int k_layer,
                                  const int dir);
  inline void corrige_flux_diff_faceIJ(IJK_Field_local_double *const flux,
                                       const int k_layer, const int dir);
  inline void update_intersections();
  inline void update();
  inline void calcul_temperature_flux_interface(const IJK_Field_double& temperature, const double ldal, const double ldav,
                                                const double dist, const DoubleTab& positions, const DoubleTab& normale,
                                                ArrOfDouble& temperature_interp, ArrOfDouble& flux_normal_interp,
                                                ArrOfDouble& temp_liqu, ArrOfDouble& temp_vap, DoubleTab& coo_liqu,
                                                DoubleTab& coo_vap) const;
  inline void compute_temperature_cell_centre(IJK_Field_double& temperature) const;
  inline void set_zero_temperature_increment(IJK_Field_double& d_temperature) const;
  inline void compute_temperature_face_centre();
  inline void compute_thermal_fluxes_face_centre();
  inline void set_convection_negligible(const int& convection_negligible);
  inline void set_diffusion_negligible(const int& diffusion_negligible);
  inline void clean();
  inline void compute_ijk_pure_faces_indices();
  inline void sort_ijk_intersections_subproblems_indices_by_k_layers();
};

inline void Corrige_flux_FT::initialize(const IJK_Splitting& splitting,
                                        const IJK_Field_double& field,
                                        const IJK_Interfaces& interfaces,
                                        const IJK_FT_double& ijk_ft,
                                        Intersection_Interface_ijk_face& intersection_ijk_face,
                                        Intersection_Interface_ijk_cell& intersection_ijk_cell)
{
  valeur().initialize(splitting, field, interfaces, ijk_ft, intersection_ijk_face, intersection_ijk_cell);
}

inline void Corrige_flux_FT::initialize_with_subproblems(const IJK_Splitting& splitting,
                                                         const IJK_Field_double& field,
                                                         const IJK_Interfaces& interfaces,
                                                         const IJK_FT_double& ijk_ft,
                                                         Intersection_Interface_ijk_face& intersection_ijk_face,
                                                         Intersection_Interface_ijk_cell& intersection_ijk_cell,
                                                         const IJK_One_Dimensional_Subproblems& thermal_local_subproblems)
{
  valeur().initialize_with_subproblems(splitting, field, interfaces, ijk_ft, intersection_ijk_face, intersection_ijk_cell, thermal_local_subproblems);
}

inline void Corrige_flux_FT::set_physical_parameters(const double rhocpl,
                                                     const double rhocpv,
                                                     const double ldal,
                                                     const double ldav)
{
  valeur().set_physical_parameters(rhocpl, rhocpv, ldal, ldav);
}

inline void Corrige_flux_FT::corrige_flux_faceIJ(IJK_Field_local_double *const flux,
                                                 const int k_layer, const int dir)
{
  valeur().corrige_flux_faceIJ(flux, k_layer, dir);
}

inline void Corrige_flux_FT::update_intersections()
{
  valeur().update_intersections();
}

inline void Corrige_flux_FT::update()
{
  valeur().update();
}

inline void Corrige_flux_FT::calcul_temperature_flux_interface(const IJK_Field_double& temperature,
                                                               const double ldal,
                                                               const double ldav,
                                                               const double dist,
                                                               const DoubleTab& positions,
                                                               const DoubleTab& normale,
                                                               ArrOfDouble& temperature_interp,
                                                               ArrOfDouble& flux_normal_interp,
                                                               ArrOfDouble& temp_liqu,
                                                               ArrOfDouble& temp_vap,
                                                               DoubleTab& coo_liqu,
                                                               DoubleTab& coo_vap) const
{
  valeur().calcul_temperature_flux_interface(temperature,
                                             ldal,
                                             ldav,
                                             dist,
                                             positions,
                                             normale,
                                             temperature_interp,
                                             flux_normal_interp,
                                             temp_liqu,
                                             temp_vap,
                                             coo_liqu,
                                             coo_vap);
}

inline void Corrige_flux_FT::compute_temperature_cell_centre(IJK_Field_double& temperature) const
{
  valeur().compute_temperature_cell_centre(temperature);
}

inline void Corrige_flux_FT::set_zero_temperature_increment(IJK_Field_double& d_temperature) const
{
  valeur().set_zero_temperature_increment(d_temperature);
}

inline void Corrige_flux_FT::compute_temperature_face_centre()
{
  valeur().compute_temperature_face_centre();
}

inline void Corrige_flux_FT::compute_thermal_fluxes_face_centre()
{
  valeur().compute_thermal_fluxes_face_centre();
}

inline void Corrige_flux_FT::corrige_flux_diff_faceIJ(IJK_Field_local_double *const flux,
                                                      const int k_layer, const int dir)
{
  valeur().corrige_flux_diff_faceIJ(flux, k_layer, dir);
}

inline void Corrige_flux_FT::set_convection_negligible(const int& convection_negligible)
{
  valeur().set_convection_negligible(convection_negligible);
}

inline void Corrige_flux_FT::set_diffusion_negligible(const int& diffusion_negligible)
{
  valeur().set_diffusion_negligible(diffusion_negligible);
}

inline void Corrige_flux_FT::clean()
{
  valeur().clean();
}

inline void Corrige_flux_FT::compute_ijk_pure_faces_indices()
{
  valeur().compute_ijk_pure_faces_indices();
}

inline void Corrige_flux_FT::sort_ijk_intersections_subproblems_indices_by_k_layers()
{
  valeur().sort_ijk_intersections_subproblems_indices_by_k_layers();
}

#endif /* Corrige_flux_FT_included */
