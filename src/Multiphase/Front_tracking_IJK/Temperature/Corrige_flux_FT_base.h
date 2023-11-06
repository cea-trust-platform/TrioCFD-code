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
// File      : Corrige_flux_FT.h
// Directory : $IJK_ROOT/src/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Corrige_flux_FT_base_included
#define Corrige_flux_FT_base_included

#include <MonofluidVar.h>
#include <Boundary_Conditions_Thermique.h>
#include <IJK_Splitting.h>
#include <Objet_U.h>
#include <Parser.h>
#include <IJK_Interfaces.h>
#include <IJK_Lata_writer.h>
#include <Intersection_Interface_ijk.h>
#include <Ouvrir_fichier.h>
#include <TRUST_Ref.h>
#include <ParcoursIJKDir.h>
#include <IJK_One_Dimensional_Subproblems.h>

class IJK_FT_double;

/*! @brief : class Corrige_flux_FT
 * API pour modifier un champ de flux à partir de donnees à l'interface. Cette
 * classe est abstraite, elle est vouée à être héritée pour être adaptée aux
 * différents cas / conditions aux limites. Par ex sur la température pour
 * imposer la continuité de lda grad T à l'interface. Mais ça pourrait être
 * pareil pour la vitesse pour imposer la continuité tangentielle et normale
 * dans le cas sans changement de phase, ou pour imposer la bonne différence de
 * vitesse normale dans le cas du changement de phase.
 */

class Corrige_flux_FT_base : public Objet_U
{
  Declare_base( Corrige_flux_FT_base ) ;
public:
  virtual void initialize(const IJK_Splitting& splitting,
                          const IJK_Field_double& field,
                          const IJK_Interfaces& interfaces,
                          const IJK_FT_double& ijk_ft,
                          Intersection_Interface_ijk_face& intersection_ijk_face,
                          Intersection_Interface_ijk_cell& intersection_ijk_cell);

  virtual void initialize_with_subproblems(const IJK_Splitting& splitting,
                                           const IJK_Field_double& field,
                                           const IJK_Interfaces& interfaces,
                                           const IJK_FT_double& ijk_ft,
                                           Intersection_Interface_ijk_face& intersection_ijk_face,
                                           Intersection_Interface_ijk_cell& intersection_ijk_cell,
                                           const IJK_One_Dimensional_Subproblems& thermal_local_subproblems);

  virtual void set_fluxes_feedback_params(const int discrete_integral, const int levels) { ; };

  void set_physical_parameters(const double rhocpl,
                               const double rhocpv,
                               const double ldal,
                               const double ldav);

  virtual void set_convection_diffusion_correction(const int& convective_flux_correction, const int& diffusive_flux_correction) { ; };
  /*
   * On va calculer sur la grille IJ du layer k_layer tous les flux a proximite de
   * l'interface. On remplace les flux donnes en entree par ces flux la.
   */
  virtual void corrige_flux_faceIJ(IJK_Field_local_double *const flux,
                                   const int k_layer, const int dir)=0;

  virtual void corrige_flux_diff_faceIJ(IJK_Field_local_double *const flux,
                                        const int k_layer, const int dir) { ; };

  virtual void correct_flux_spherical(Simd_double& a,
                                      Simd_double& b,
                                      const int& i,
                                      const int& j,
                                      const int& k_layer,
                                      const int dir) { ; };

  virtual void update_intersections() { ; };
  virtual void update()=0;

  virtual void calcul_temperature_flux_interface(const IJK_Field_double& temperature, const double ldal, const double ldav,
                                                 const double dist, const DoubleTab& positions, const DoubleTab& normale,
                                                 ArrOfDouble& temperature_interp, ArrOfDouble& flux_normal_interp,
                                                 ArrOfDouble& temp_liqu, ArrOfDouble& temp_vap, DoubleTab& coo_liqu,
                                                 DoubleTab& coo_vap) const = 0 ;


  virtual void compute_temperature_cell_centre(IJK_Field_double& temperature) const { ; };
  virtual void set_zero_temperature_increment(IJK_Field_double& d_temperature) const { ; };

  virtual void compute_thermal_convective_fluxes() { ; };
  virtual void compute_thermal_diffusive_fluxes() { ; };

  virtual void set_convection_negligible(const int& convection_negligible) { ; };
  virtual void set_diffusion_negligible(const int& diffusion_negligible) { ; };
  virtual void clean() { ; };
  virtual void compute_ijk_pure_faces_indices() { ; };
  virtual void sort_ijk_intersections_subproblems_indices_by_k_layers() { ; };

  virtual void set_debug(const int& debug) { ; };
  virtual void set_distance_cell_faces_from_lrs(const int& distance_cell_faces_from_lrs) { ; };
  virtual void set_correction_cell_neighbours(const int& correct_temperature_cell_neighbours, const int& neighbours_colinearity_weighting) { ; };
  virtual void set_cell_faces_neighbours_corrected_bool(FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_bool) { ; };
  virtual void set_eulerian_normal_vectors_ns_normed(FixedVector<IJK_Field_double, 3>& eulerian_normal_vectors_ns_normed) { ; };

  virtual void set_correction_cell_faces_neighbours(const int& find_cell_neighbours_for_fluxes_spherical_correction,
                                                    const int& use_cell_neighbours_for_fluxes_spherical_correction,
                                                    const int& find_reachable_fluxes) { ; };
  virtual void initialise_cell_neighbours_indices_to_correct() { ; };
  virtual void compute_cell_neighbours_faces_indices_for_spherical_correction(const int& n_iter_distance) { ; };
  virtual void compute_cell_neighbours_faces_indices_to_correct(FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_bool,
                                                                FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_convective,
                                                                FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_diffusive,
                                                                FixedVector<IJK_Field_double, 3>& neighbours_weighting_colinearity,
                                                                const int& compute_fluxes_values) { ; };
  virtual void compute_temperature_cell_centre_neighbours(IJK_Field_double& temperature_neighbours,
                                                          IJK_Field_int& neighbours_weighting,
                                                          IJK_Field_double& neighbours_weighting_colinearity) const { ; };
  virtual void replace_temperature_cell_centre_neighbours(IJK_Field_double& temperature,
                                                          IJK_Field_double& temperature_neighbours,
                                                          IJK_Field_int& neighbours_weighting,
                                                          IJK_Field_double& neighbours_weighting_colinearity) const { ; };

  virtual void store_cell_faces_corrected(FixedVector<IJK_Field_int,3>& cell_faces_corrected_bool,
                                          FixedVector<IJK_Field_double,3>& cell_faces_corrected_convective,
                                          FixedVector<IJK_Field_double,3>& cell_faces_corrected_diffusive) { ; };
  virtual void clear_vectors() { ; };
  virtual void compute_min_max_ijk_reachable_fluxes(const FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_all_bool,
                                                    FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_min_max_bool) { ; };

protected:
  const IJK_Interfaces *interfaces_;
  const IJK_Field_double *field_;
  const IJK_Splitting *splitting_;
  REF(IJK_FT_double) ref_ijk_ft_;


  double rhocp_l_, rhocp_v_;
  double lda_l_, lda_v_;

  Intersection_Interface_ijk_face * intersection_ijk_face_;
  Intersection_Interface_ijk_cell * intersection_ijk_cell_;
  /*
   * TODO: mettre ces méthodes dans une petite classe pour parcourir
   * les trois directions.
   */
  bool test_if_stencil_inclut_bout_interface_liquide() const ;
  bool test_if_stencil_inclut_bout_interface_vapeur() const ;
  ParcoursIJKDir parcours_;
};

#endif /* Corrige_flux_FT_base_included */
