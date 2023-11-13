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
#define NEIGHBOURS_FACES_I {0, 1, 0, 0, 0, 0}
#define NEIGHBOURS_FACES_J {0, 0, 0, 1, 0, 0}
#define NEIGHBOURS_FACES_K {0, 0, 0, 0, 0, 1}
/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Corrige_flux_FT_temperature_subresolution
//
// <Description of class Corrige_flux_FT_temperature_subresolution>
//
/////////////////////////////////////////////////////////////////////////////
#define FACES_DIR {0, 0, 1, 1, 2, 2}
// TODO: Be careful to operators ! (Left - Right) values
#define FLUX_SIGN_DIFF {-1, -1, -1, -1, -1, -1}
#define FLUX_SIGN_CONV {1, 1, 1, 1, 1, 1}
#define FLUX_SIGN {FLUX_SIGN_CONV, FLUX_SIGN_DIFF}


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
  void set_convection_diffusion_correction(const int& convective_flux_correction, const int& diffusive_flux_correction) override
  {
    convective_flux_correction_ = convective_flux_correction;
    diffusive_flux_correction_ = diffusive_flux_correction;
  }
  void set_convection_negligible(const int& convection_negligible) override { convection_negligible_ = convection_negligible; };
  void set_diffusion_negligible(const int& diffusion_negligible) override { diffusion_negligible_ = diffusion_negligible; };
  void set_fluxes_feedback_params(const int discrete_integral, const int levels) override { discrete_integral_ = discrete_integral; levels_ = levels; };
  void set_distance_cell_faces_from_lrs(const int& distance_cell_faces_from_lrs) override { distance_cell_faces_from_lrs_=distance_cell_faces_from_lrs; };
  void set_correction_cell_neighbours(const int& correct_temperature_cell_neighbours,
                                      const int& neighbours_colinearity_weighting) override
  {
    find_temperature_cell_neighbours_ = correct_temperature_cell_neighbours;
    neighbours_colinearity_weighting_ = neighbours_colinearity_weighting;
  }

  void set_cell_faces_neighbours_corrected_bool(FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_bool) override
  {
    cell_faces_neighbours_corrected_bool_ = &cell_faces_neighbours_corrected_bool;
  }

  void set_eulerian_normal_vectors_ns_normed(FixedVector<IJK_Field_double, 3>& eulerian_normal_vectors_ns_normed) override
  {
    eulerian_normal_vectors_ns_normed_ = &eulerian_normal_vectors_ns_normed;
  }

  void set_correction_cell_faces_neighbours(const int& find_cell_neighbours_for_fluxes_spherical_correction,
                                            const int& use_cell_neighbours_for_fluxes_spherical_correction,
                                            const int& find_reachable_fluxes,
                                            const int& use_reachable_fluxes) override
  {
    find_cell_neighbours_for_fluxes_spherical_correction_ = find_cell_neighbours_for_fluxes_spherical_correction;
    use_cell_neighbours_for_fluxes_spherical_correction_ = use_cell_neighbours_for_fluxes_spherical_correction;
    find_reachable_fluxes_ = find_reachable_fluxes;
    use_reachable_fluxes_ = use_reachable_fluxes;
  }
  void set_debug(const int& debug) override { debug_ = debug; };
  /*
   * On va calculer sur la grille IJ du layer k_layer tous les flux a proximite de
   * l'interface. On remplace les flux donnes en entree par ces flux la.
   */
  void corrige_flux_faceIJ_any_flux(IJK_Field_local_double *const flux,
                                    std::vector<ArrOfInt>& index_face_i_flux_x_sorted,
                                    std::vector<ArrOfInt>& index_face_j_flux_x_sorted,
                                    std::vector<ArrOfInt>& index_face_i_flux_y_sorted,
                                    std::vector<ArrOfInt>& index_face_j_flux_y_sorted,
                                    std::vector<ArrOfInt>& index_face_i_flux_z_sorted,
                                    std::vector<ArrOfInt>& index_face_j_flux_z_sorted,
                                    std::vector<ArrOfDouble>& subgrid_fluxes_x,
                                    std::vector<ArrOfDouble>& subgrid_fluxes_y,
                                    std::vector<ArrOfDouble>& subgrid_fluxes_z,
                                    const int k_layer,
                                    const int dir);

  void corrige_flux_faceIJ(IJK_Field_local_double *const flux,
                           const int k_layer, const int dir) override;

  void corrige_flux_conv_faceIJ(IJK_Field_local_double *const flux,
                                const int k_layer, const int dir) { corrige_flux_faceIJ(flux, k_layer, dir); };

  void corrige_flux_diff_faceIJ(IJK_Field_local_double *const flux,
                                const int k_layer, const int dir) override;

  void correct_flux_spherical(Simd_double& a,
                              Simd_double& b,
                              const int& i,
                              const int& j,
                              const int& k_layer,
                              const int dir) override;

  void calcul_temperature_flux_interface(const IJK_Field_double& temperature, const double ldal, const double ldav,
                                         const double dist, const DoubleTab& positions, const DoubleTab& normale,
                                         ArrOfDouble& temperature_interp, ArrOfDouble& flux_normal_interp,
                                         ArrOfDouble& temp_liqu, ArrOfDouble& temp_vap, DoubleTab& coo_liqu,
                                         DoubleTab& coo_vap) const override { ; };

  void update_intersections() override;
  void update() override;
  void associate_indices_and_check_subproblems_consistency();
  void compute_temperature_cell_centre(IJK_Field_double& temperature) const override;
  void compute_temperature_cell_centre_neighbours(IJK_Field_double& temperature_neighbours,
                                                  IJK_Field_int& neighbours_weighting,
                                                  IJK_Field_double& neighbours_weighting_colinearity) const override;
  void replace_temperature_cell_centre_neighbours(IJK_Field_double& temperature,
                                                  IJK_Field_double& temperature_neighbours,
                                                  IJK_Field_int& neighbours_weighting,
                                                  IJK_Field_double& neighbours_weighting_colinearity) const override;
  void initialise_cell_neighbours_indices_to_correct() override;
  void initialise_any_cell_neighbours_indices_to_correct(std::vector<ArrOfInt>& index_face_i_flux_x_faces_sorted,
                                                         std::vector<ArrOfInt>& index_face_j_flux_x_faces_sorted,
                                                         std::vector<ArrOfInt>& index_face_i_flux_y_faces_sorted,
                                                         std::vector<ArrOfInt>& index_face_j_flux_y_faces_sorted,
                                                         std::vector<ArrOfInt>& index_face_i_flux_z_faces_sorted,
                                                         std::vector<ArrOfInt>& index_face_j_flux_z_faces_sorted);
  void initialise_any_cell_neighbours_indices_to_correct_with_flux(std::vector<ArrOfInt>& index_face_i_flux_x_faces_sorted,
                                                                   std::vector<ArrOfInt>& index_face_j_flux_x_faces_sorted,
                                                                   std::vector<ArrOfInt>& index_face_i_flux_y_faces_sorted,
                                                                   std::vector<ArrOfInt>& index_face_j_flux_y_faces_sorted,
                                                                   std::vector<ArrOfInt>& index_face_i_flux_z_faces_sorted,
                                                                   std::vector<ArrOfInt>& index_face_j_flux_z_faces_sorted,
                                                                   std::vector<ArrOfDouble>& flux_x,
                                                                   std::vector<ArrOfDouble>& flux_y,
                                                                   std::vector<ArrOfDouble>& flux_z,
                                                                   const bool& ini_index);
  void compute_cell_neighbours_faces_indices_for_spherical_correction(const int& n_iter_distance) override;
  void compute_cell_neighbours_faces_indices_to_correct(FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_bool,
                                                        FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_convective,
                                                        FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_diffusive,
                                                        FixedVector<IJK_Field_double, 3>& neighbours_weighting_colinearity) override;
  void complete_neighbours_and_weighting_colinearity(FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_bool,
                                                     FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_convective,
                                                     FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_diffusive,
                                                     FixedVector<IJK_Field_double, 3>& neighbours_weighting_colinearity,
                                                     const int& compute_fluxes_values);
  void compute_cell_neighbours_fluxes_to_correct(FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_convective,
                                                 FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_corrected_diffusive,
                                                 const int& subproblem_index,
                                                 const int& index_i, const int& index_j, const int& index_k,
                                                 const int& index_i_perio, const int& index_j_perio, const int& index_k_perio,
                                                 const double& dist,
                                                 const int& dir,
                                                 const double& colinearity,
                                                 const int& compute_fluxes_values);
  void compute_cell_neighbours_convective_fluxes_to_correct(double& convective_flux,
                                                            const int& subproblem_index,
                                                            const int& index_i, const int& index_j, const int& index_k,
                                                            const double& dist,
                                                            const int& dir,
                                                            const double& colinearity);
  void compute_cell_neighbours_thermal_convective_fluxes_face_centre(double& convective_flux,
                                                                     const int& subproblem_index,
                                                                     const int& index_i, const int& index_j, const int& index_k,
                                                                     const double& dist,
                                                                     const int& dir,
                                                                     const double& colinearity);
  void compute_cell_neighbours_thermal_convective_fluxes_face_centre_discrete_integral(double& convective_flux,
                                                                                       const int& subproblem_index,
                                                                                       const int& index_i, const int& index_j, const int& index_k,
                                                                                       const double& dist,
                                                                                       const int& dir,
                                                                                       const double& colinearity);
  void compute_cell_neighbours_diffusive_fluxes_to_correct(double& diffusive_flux,
                                                           const int& subproblem_index,
                                                           const int& index_i,
                                                           const int& index_j,
                                                           const int& index_k,
                                                           const double& dist,
                                                           const int& dir,
                                                           const double& colinearity);
  void compute_cell_neighbours_thermal_diffusive_fluxes_face_centre(double& diffusive_flux,
                                                                    const int& subproblem_index,
                                                                    const int& index_i, const int& index_j, const int& index_k,
                                                                    const double& dist,
                                                                    const int& dir,
                                                                    const double& colinearity);
  void compute_cell_neighbours_thermal_diffusive_fluxes_face_centre_discrete_integral(double& diffusive_flux,
                                                                                      const int& subproblem_index,
                                                                                      const int& index_i, const int& index_j, const int& index_k,
                                                                                      const double& dist,
                                                                                      const int& dir,
                                                                                      const double& colinearity);
  void compute_cell_neighbours_thermal_fluxes_face_centre(double& flux,
                                                          const int fluxes_type,
                                                          const int& subproblem_index,
                                                          const int& index_i, const int& index_j, const int& index_k,
                                                          const double& dist,
                                                          const int& dir,
                                                          const double& colinearity);
  void compute_cell_neighbours_thermal_fluxes_face_centre_discrete_integral(double& flux,
                                                                            const int fluxes_type,
                                                                            const int& subproblem_index,
                                                                            const int& index_i, const int& index_j, const int& index_k,
                                                                            const double& dist,
                                                                            const int& dir,
                                                                            const double& colinearity);
  void replace_cell_neighbours_thermal_convective_diffusive_fluxes_faces(const FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_min_max_bool,
                                                                         const FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_fluxes_corrected,
                                                                         const int& fluxes_type) override;
  void replace_cell_neighbours_thermal_fluxes_faces(const FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_min_max_bool,
                                                    const FixedVector<IJK_Field_double, 3>& cell_faces_neighbours_fluxes_corrected,
                                                    std::vector<ArrOfDouble>& flux_x,
                                                    std::vector<ArrOfDouble>& flux_y,
                                                    std::vector<ArrOfDouble>& flux_z);

  void set_zero_temperature_increment(IJK_Field_double& d_temperature) const override;
  void compute_thermal_convective_fluxes() override;
  void compute_thermal_diffusive_fluxes() override;
  void compute_thermal_convective_fluxes_face_centre();
  void compute_thermal_diffusive_fluxes_face_centre();
  void compute_thermal_fluxes_face_centre(DoubleVect& fluxes, const int fluxes_type);
  double compute_thermal_flux_face_centre(const int fluxes_type, const int& index_subproblem, const double& dist, const int& dir);
  void compute_thermal_convective_fluxes_face_centre_discrete_integral();
  void compute_thermal_diffusive_fluxes_face_centre_discrete_integral();
  void compute_thermal_fluxes_face_centre_discrete_integral(DoubleVect& fluxes, const int fluxes_type);
  DoubleVect compute_thermal_flux_face_centre_discrete_integral(const int fluxes_type, const int& index_subproblem, const double& dist, const int& dir);
  void get_discrete_surface_at_level(const int& dir, const int& level);
  void clean() override;
  void compute_ijk_pure_faces_indices() override;
  void sort_ijk_intersections_subproblems_indices_by_k_layers() override;
  void sort_ijk_intersections_subproblems_indices_fluxes_by_k_layers(std::vector<ArrOfInt>& index_face_i_flux_x,
                                                                     std::vector<ArrOfInt>& index_face_j_flux_x,
                                                                     std::vector<ArrOfInt>& index_face_i_flux_y,
                                                                     std::vector<ArrOfInt>& index_face_j_flux_y,
                                                                     std::vector<ArrOfInt>& index_face_i_flux_z,
                                                                     std::vector<ArrOfInt>& index_face_j_flux_z,
                                                                     std::vector<ArrOfDouble>& flux_x,
                                                                     std::vector<ArrOfDouble>& flux_y,
                                                                     std::vector<ArrOfDouble>& flux_z,
                                                                     const DoubleVect& fluxes_subgrid,
                                                                     const int ini_index);
  void store_cell_faces_corrected(FixedVector<IJK_Field_int,3>& cell_faces_corrected_bool,
                                  FixedVector<IJK_Field_double,3>& cell_faces_corrected_convective,
                                  FixedVector<IJK_Field_double,3>& cell_faces_corrected_diffusive) override;
  void store_any_cell_faces_corrected(FixedVector<IJK_Field_int,3>& cell_faces_corrected_bool,
                                      FixedVector<IJK_Field_double,3>& cell_faces_corrected,
                                      const DoubleVect& fluxes,
                                      const int counter);
  void check_pure_fluxes_duplicates(const DoubleVect& fluxes, DoubleVect& fluxes_unique, IntVect& pure_face_unique, const int known_unique);
  void clear_vectors() override;
  void compute_min_max_ijk_reachable_fluxes(const FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_all_bool,
                                            const IJK_Field_int& neighbours_temperature_to_correct,
                                            FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_min_max_bool,
                                            const int& max_flux_per_dir,
                                            const int& check_cell_center_neighbour,
                                            const int& remove_external_neighbour_values,
                                            IJK_Field_int& neighbours_temperature_to_correct_trimmed) override;
  void remove_min_max_ijk_reachable_fluxes_discontinuous(const FixedVector<IJK_Field_int, 3>& cell_faces_neighbours_corrected_all_bool,
                                                         FixedVector<IJK_Field_local_int, 3>& cell_faces_neighbours_corrected_min_max_bool);
protected :
  enum fluxes_type_ { convection, diffusion };
  DoubleVect dist_;
  IntVect pure_face_unique_;
  DoubleVect convective_fluxes_;
  DoubleVect convective_fluxes_unique_;
  DoubleVect diffusive_fluxes_;
  DoubleVect diffusive_fluxes_unique_;
  const IJK_One_Dimensional_Subproblems * thermal_subproblems_;
  bool has_checked_consistency_;
  ArrOfInt ijk_intersections_subproblems_indices_;
  /*
   * To be used in the operator
   * Store (i, j, flux) at a given row (k)
   * convective_flux_x_sorted_.pushback(double_tab_ij_flux)
   */
  std::vector<ArrOfInt> index_face_i_flux_x_sorted_;
  std::vector<ArrOfInt> index_face_j_flux_x_sorted_;
  std::vector<ArrOfInt> index_face_i_flux_y_sorted_;
  std::vector<ArrOfInt> index_face_j_flux_y_sorted_;
  std::vector<ArrOfInt> index_face_i_flux_z_sorted_;
  std::vector<ArrOfInt> index_face_j_flux_z_sorted_;
  std::vector<ArrOfDouble> convective_flux_x_sorted_;
  std::vector<ArrOfDouble> convective_flux_y_sorted_;
  std::vector<ArrOfDouble> convective_flux_z_sorted_;
  std::vector<ArrOfDouble> diffusive_flux_x_sorted_;
  std::vector<ArrOfDouble> diffusive_flux_y_sorted_;
  std::vector<ArrOfDouble> diffusive_flux_z_sorted_;
  FixedVector<IntVect,4> ijk_faces_to_correct_;

  std::vector<ArrOfInt> index_face_i_flux_x_neighbours_diag_faces_sorted_;
  std::vector<ArrOfInt> index_face_j_flux_x_neighbours_diag_faces_sorted_;
  std::vector<ArrOfInt> index_face_i_flux_y_neighbours_diag_faces_sorted_;
  std::vector<ArrOfInt> index_face_j_flux_y_neighbours_diag_faces_sorted_;
  std::vector<ArrOfInt> index_face_i_flux_z_neighbours_diag_faces_sorted_;
  std::vector<ArrOfInt> index_face_j_flux_z_neighbours_diag_faces_sorted_;

  std::vector<ArrOfInt> index_face_i_flux_x_neighbours_all_faces_sorted_;
  std::vector<ArrOfInt> index_face_j_flux_x_neighbours_all_faces_sorted_;
  std::vector<ArrOfInt> index_face_i_flux_y_neighbours_all_faces_sorted_;
  std::vector<ArrOfInt> index_face_j_flux_y_neighbours_all_faces_sorted_;
  std::vector<ArrOfInt> index_face_i_flux_z_neighbours_all_faces_sorted_;
  std::vector<ArrOfInt> index_face_j_flux_z_neighbours_all_faces_sorted_;

  std::vector<ArrOfInt> index_face_i_flux_x_neighbours_min_max_faces_sorted_;
  std::vector<ArrOfInt> index_face_j_flux_x_neighbours_min_max_faces_sorted_;
  std::vector<ArrOfInt> index_face_i_flux_y_neighbours_min_max_faces_sorted_;
  std::vector<ArrOfInt> index_face_j_flux_y_neighbours_min_max_faces_sorted_;
  std::vector<ArrOfInt> index_face_i_flux_z_neighbours_min_max_faces_sorted_;
  std::vector<ArrOfInt> index_face_j_flux_z_neighbours_min_max_faces_sorted_;

  std::vector<ArrOfDouble> convective_flux_x_min_max_faces_sorted_;
  std::vector<ArrOfDouble> convective_flux_y_min_max_faces_sorted_;
  std::vector<ArrOfDouble> convective_flux_z_min_max_faces_sorted_;
  std::vector<ArrOfDouble> diffusive_flux_x_min_max_faces_sorted_;
  std::vector<ArrOfDouble> diffusive_flux_y_min_max_faces_sorted_;
  std::vector<ArrOfDouble> diffusive_flux_z_min_max_faces_sorted_;

  int convection_negligible_;
  int diffusion_negligible_;
  int debug_;
  int levels_;
  int discrete_integral_;
  bool flux_init_;

  int distance_cell_faces_from_lrs_;
  int find_temperature_cell_neighbours_;
  int find_cell_neighbours_for_fluxes_spherical_correction_;
  int use_cell_neighbours_for_fluxes_spherical_correction_;
  int neighbours_colinearity_weighting_;

  int find_reachable_fluxes_;
  int use_reachable_fluxes_;
  FixedVector<IJK_Field_int, 3>* cell_faces_neighbours_corrected_bool_;
  FixedVector<IJK_Field_double, 3>* eulerian_normal_vectors_ns_normed_;

  int convective_flux_correction_;
  int diffusive_flux_correction_;

};

#endif /* Corrige_flux_FT_temperature_subresolution_included */
