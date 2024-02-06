/****************************************************************************
* Copyright (c) 2024, CEA
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
// File      : IJK_One_Dimensional_Subproblem_Geometry.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_One_Dimensional_Subproblem_Geometry_included
#define IJK_One_Dimensional_Subproblem_Geometry_included

#include <Objet_U.h>
#include <IJK_Field.h>
#include <IJK_Interfaces.h>
#include <FixedVector.h>
#include <Vecteur3.h>

#define INVALID_TEMPERATURE 1e10
#define INVALID_FIELD 1e10
#define INVALID_VELOCITY 1e-12
#define INVALID_INTERP 1.e20
#define INVALID_INTERP_TEST 1.e19
#define INVALID_VELOCITY_CFL 1e-20
#define INVALID_SOURCE_TERM 1e-20
#define NEIGHBOURS_FIRST_DIR {-1., -1., 1., 1.}
#define NEIGHBOURS_SECOND_DIR {-1., 1., -1., 1.}
#define NEIGHBOURS_I {-1, 1, 0, 0, 0, 0}
#define NEIGHBOURS_J {0, 0, -1, 1, 0, 0}
#define NEIGHBOURS_K {0, 0, 0, 0, -1, 1}
#define NEIGHBOURS_FACES_I {0, 1, 0, 0, 0, 0}
#define NEIGHBOURS_FACES_J {0, 0, 0, 1, 0, 0}
#define NEIGHBOURS_FACES_K {0, 0, 0, 0, 0, 1}
#define LIQUID_INDICATOR_TEST 1.-1.e-12
#define VAPOUR_INDICATOR_TEST 1.e-12
#define FACES_DIR {0, 0, 1, 1, 2, 2}
#define FLUX_SIGN_DIFF {-1, -1, -1, -1, -1, -1}
#define FLUX_SIGN_CONV {1, 1, 1, 1, 1, 1}

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_One_Dimensional_Subproblem_Geometry
//
// <Description of class IJK_One_Dimensional_Subproblem_Geometry>
//
/////////////////////////////////////////////////////////////////////////////

class IJK_FT_double;
class IJK_One_Dimensional_Subproblem_Geometry : public Objet_U
{

  Declare_instanciable( IJK_One_Dimensional_Subproblem_Geometry ) ;
  friend class IJK_One_Dimensional_Subproblem;

public :

  void compute_interface_basis_vectors();
  void compute_pure_spherical_basis_vectors();
  void compute_local_discretisation();

  void interpolate_indicator_on_probes();
  double find_cell_related_indicator_on_probes(const int& last_index);

  void compute_distance_cell_centre();
  void compute_distance_faces_centres();

  void compute_distance_cell_centres_neighbours();
  double compute_cell_weighting(const double& dx_contrib,
                                const double& dy_contrib,
                                const double& dz_contrib);


  void compute_distance_last_cell_faces_neighbours();
  double compute_cell_faces_weighting(const double& dx_contrib,
                                      const double& dy_contrib,
                                      const double& dz_contrib,
                                      const int& dir);
  Vecteur3 compute_relative_vector_cell_faces(const double& dx_contrib,
                                              const double& dy_contrib,
                                              const double& dz_contrib);
  double compute_colinearity(const double& dx_contrib,
                             const double& dy_contrib,
                             const double& dz_contrib);
  double compute_colinearity_cell_faces(const double& dx_contrib,
                                        const double& dy_contrib,
                                        const double& dz_contrib,
                                        const int& dir);
  double compute_distance_cell_faces(const double& dx_contrib,
                                     const double& dy_contrib,
                                     const double& dz_contrib);

  int get_dxyz_increment_max();
  int get_dxyz_over_two_increment_max();
  void get_maximum_remaining_distance(int& dx_remaining,
                                      int& dy_remaining,
                                      int& dz_remaining);

protected :

  void compute_vertex_position(const int& vertex_number,
                               const int& face_dir,
                               Vecteur3& bary_vertex,
                               double& distance_vertex_centre,
                               double& tangential_distance_vertex_centre,
                               Vecteur3& tangential_distance_vector_vertex_centre);

  REF(IJK_FT_double) ref_ijk_ft_;
  const IJK_Interfaces * interfaces_;
  int debug_ = 0;

  int * points_per_thermal_subproblem_;
  double probe_length_ = 0.;
  DoubleTab coordinates_cartesian_compo_;

  int disable_probe_collision_;
  int enable_resize_probe_collision_;
  int disable_find_cell_centre_probe_tip_;
  int resize_probe_collision_index_;
  DoubleVect indicator_interp_;

  Vecteur3 facet_barycentre_;
  Vecteur3 normal_vector_compo_;

  int index_i_ = 0, index_j_ = 0, index_k_ = 0;

  double cell_centre_distance_ = 0;
  double cell_centre_tangential_distance_ = 0.;
  Vecteur3 tangential_distance_vector_;
  FixedVector<bool,6> pure_liquid_neighbours_;
  FixedVector<double,6> face_centres_distance_;
  FixedVector<double,6> face_centres_tangential_distance_;
  FixedVector<Vecteur3,6> face_tangential_distance_vector_;
  FixedVector<FixedVector<double,4>,6> vertices_centres_distance_;
  FixedVector<FixedVector<double,4>,6> vertices_centres_tangential_distance_;
  FixedVector<FixedVector<Vecteur3,4>,6> vertices_tangential_distance_vector_;
  double modified_probe_length_from_vertices_ = 0.;
  bool has_computed_cell_centre_distance_=false;
  bool has_computed_cell_faces_distance_=false;
  int correct_fluxes_ = 0;

  /*
   * Identify neighbours cells centre for temperature correction
   */
  int correct_temperature_cell_neighbours_ = 0;
  int correct_neighbours_rank_ = 1;
  int neighbours_corrected_rank_ = 1;
  int neighbours_weighting_= 0;
  int neighbours_colinearity_weighting_ = 0;
  int neighbours_distance_weighting_ = 0;
  int neighbours_colinearity_distance_weighting_ = 0;
  FixedVector<int,3> pure_neighbours_corrected_sign_;
  std::vector<std::vector<std::vector<bool>>> pure_neighbours_to_correct_;
  std::vector<std::vector<std::vector<double>>> pure_neighbours_corrected_distance_;
  std::vector<std::vector<std::vector<double>>> pure_neighbours_corrected_colinearity_;
  int dxyz_increment_bool_ = 0;
  int dxyz_over_two_increment_bool_ = 0;

  int find_cell_neighbours_for_fluxes_spherical_correction_ = 0;

  /*
   * Identify neighbours faces centres for flux correction
   */
  int compute_reachable_fluxes_ = 0;
  int neighbours_last_faces_weighting_ = 0;
  int neighbours_last_faces_colinearity_weighting_ = 0;
  int neighbours_last_faces_colinearity_face_weighting_ = 0.;
  int neighbours_last_faces_distance_weighting_ = 0.;
  int neighbours_last_faces_distance_colinearity_weighting_ = 0.;
  int neighbours_last_faces_distance_colinearity_face_weighting_ = 0.;
  int neighbours_face_corrected_rank_ = 1;
  std::vector<std::vector<std::vector<std::vector<bool>>>> pure_neighbours_last_faces_to_correct_;
  std::vector<std::vector<std::vector<std::vector<double>>>> pure_neighbours_last_faces_corrected_distance_;
  std::vector<std::vector<std::vector<std::vector<double>>>> pure_neighbours_last_faces_corrected_colinearity_;

};

#endif /* IJK_One_Dimensional_Subproblem_Geometry_included */
