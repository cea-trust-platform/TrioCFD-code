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
// File      : IJK_Ghost_Fluid_tools.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Field.h>
#include <IJK_Interfaces.h>
#include <Maillage_FT_IJK.h>

#define INVALID_TEST -1.e30
//Be coherent with LOCAL EPS of Intersection Interface IJK
#define LIQUID_INDICATOR_TEST 1.-1.e-12
#define VAPOUR_INDICATOR_TEST 1.e-12
#define NEIGHBOURS_I {-1, 1, 0, 0, 0, 0}
#define NEIGHBOURS_J {0, 0, -1, 1, 0, 0}
#define NEIGHBOURS_K {0, 0, 0, 0, -1, 1}


/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_Ghost_Fluid_tools
//
// <Description of class IJK_Ghost_Fluid_tools>
//
/////////////////////////////////////////////////////////////////////////////

void compute_eulerian_normal_distance_facet_barycentre_field(const IJK_Interfaces& interface,
                                                             IJK_Field_double& distance,
                                                             FixedVector<IJK_Field_double, 3>& normal_vect,
                                                             FixedVector<IJK_Field_double, 3>& facets_barycentre,
                                                             const int& n_iter);

void compute_eulerian_curvature_field_from_distance_field(const IJK_Field_double& distance,
                                                          IJK_Field_double& curvature,
                                                          const IJK_Field_local_double& boundary_flux_kmin,
                                                          const IJK_Field_local_double& boundary_flux_kmax);

void compute_eulerian_curvature_field_from_normal_vector_field(const FixedVector<IJK_Field_double, 3>& normal_vect,
                                                               IJK_Field_double& curvature);

void compute_eulerian_curvature_field_from_interface(const FixedVector<IJK_Field_double, 3>& normal_vect,
                                                     const IJK_Interfaces& interfaces,
                                                     IJK_Field_double& interfacial_area,
                                                     IJK_Field_double& curvature,
                                                     const int& n_iter,
                                                     int igroup);

void compute_eulerian_normal_temperature_gradient_interface(const IJK_Field_double& distance,
                                                            const IJK_Field_double& indicator,
                                                            const IJK_Field_double& interfacial_area,
                                                            const IJK_Field_double& curvature,
                                                            const	IJK_Field_double& temperature,
                                                            IJK_Field_double& grad_T_interface);

void propagate_eulerian_normal_temperature_gradient_interface(const IJK_Interfaces& interfaces,
                                                              const IJK_Field_double& distance,
                                                              IJK_Field_double& grad_T_interface,
                                                              const int stencil_width);

void compute_eulerian_extended_temperature(const IJK_Field_double& indicator,
                                           const IJK_Field_double& distance,
                                           const IJK_Field_double& curvature,
                                           IJK_Field_double& grad_T_interface,
                                           IJK_Field_double& temperature);
