//TRUST_NO_INDENT
/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
// File      : Grid_Level_Data.h
// Directory : $IJK_ROOT/src/IJK/solveur_mg
//
/////////////////////////////////////////////////////////////////////////////
//
// WARNING: DO NOT EDIT THIS FILE! Only edit the template file Grid_Level_Data.h.P
//
#ifndef Grid_Level_Data_H
#define Grid_Level_Data_H

#include <IJK_Field.h>
#include <IJK_Splitting.h>


// Data for each grid
class Grid_Level_Data_double : public Objet_U
{
  Declare_instanciable(Grid_Level_Data_double);
 public:
  void initialize(const IJK_Splitting &, int ghost, int additional_k_layers);

  int get_ghost_size() const { return ghost_size_; }
  const IJK_Splitting & get_splitting() const { return grid_splitting_; }
  const IJK_Grid_Geometry & get_grid_geometry() const { return grid_splitting_.get_grid_geometry(); }
  // Compute the ijk_faces_coefficients from ijk_rho_
  void compute_faces_coefficients_from_rho();
  void compute_faces_coefficients_from_inv_rho();
  // Returns the reference to the rho_field (to fill the data)
  IJK_Field_double & get_update_rho() { return ijk_rho_; }
  const IJK_Field_double & get_rho() const { return ijk_rho_; }
  IJK_Field_double & get_update_x() { return ijk_x_; }
  IJK_Field_double & get_update_rhs() { return ijk_rhs_; }
  IJK_Field_double & get_update_residue() { return ijk_residue_; }
  const IJK_Field_double & get_faces_coefficients() const { return ijk_faces_coefficients_; }

 protected:
  void compute_faces_coefficients_from_rho_cst_i_cst_j_cst_k();
  void compute_faces_coefficients_from_inv_rho_cst_i_cst_j_cst_k();
  void compute_faces_coefficients_from_rho_cst_i_cst_j_var_k();
  void compute_faces_coefficients_from_inv_rho_cst_i_cst_j_var_k();
  //IJK_Grid_Geometry grid_geometry_;
  IJK_Splitting grid_splitting_;
  int ghost_size_;
  bool perio_k_;

  VECT(ArrOfDouble_with_ghost) local_delta_xyz_;
  bool uniform_[3];

  // Poisson coefficient on each cell on this level
  IJK_Field_double ijk_rho_;
  // matrix coefficients at faces for each level (computed from ijk_rho_levels_)
  // These coefficients contain the wall boundary conditions.
  // The field has 4 components per element
  // coeff_faces(i,j,k,0) contains data for faces of normal x
  // coeff_faces(i,j,k,1) contains data for faces of normal y
  // coeff_faces(i,j,k,2) contains data for faces of normal z
  // coeff_faces(i,j,k,3) contains the sum of the coefficients at all faces of element (i,j,k)
  IJK_Field_double ijk_faces_coefficients_;

  // Temporary storage for multigrid algorithm:
  // Storage for the unknown:
  IJK_Field_double ijk_x_;
  // Storage for the right-hand side:
  IJK_Field_double ijk_rhs_;
  // Storage for the residue:
  IJK_Field_double ijk_residue_;
};

Declare_vect(Grid_Level_Data_double);

// Data for each grid
class Grid_Level_Data_float : public Objet_U
{
  Declare_instanciable(Grid_Level_Data_float);
 public:
  void initialize(const IJK_Splitting &, int ghost, int additional_k_layers);

  int get_ghost_size() const { return ghost_size_; }
  const IJK_Splitting & get_splitting() const { return grid_splitting_; }
  const IJK_Grid_Geometry & get_grid_geometry() const { return grid_splitting_.get_grid_geometry(); }
  // Compute the ijk_faces_coefficients from ijk_rho_
  void compute_faces_coefficients_from_rho();
  void compute_faces_coefficients_from_inv_rho();
  void compute_faces_coefficients_from_double_coeffs(const Grid_Level_Data_double &);
  // Returns the reference to the rho_field (to fill the data)
  IJK_Field_float & get_update_rho() { return ijk_rho_; }
  const IJK_Field_float & get_rho() const { return ijk_rho_; }
  IJK_Field_float & get_update_x() { return ijk_x_; }
  IJK_Field_float & get_update_rhs() { return ijk_rhs_; }
  IJK_Field_float & get_update_residue() { return ijk_residue_; }
  const IJK_Field_float & get_faces_coefficients() const { return ijk_faces_coefficients_; }

 protected:
  void compute_faces_coefficients_from_rho_cst_i_cst_j_cst_k();
  void compute_faces_coefficients_from_inv_rho_cst_i_cst_j_cst_k();
  void compute_faces_coefficients_from_rho_cst_i_cst_j_var_k();
  void compute_faces_coefficients_from_inv_rho_cst_i_cst_j_var_k();
  //IJK_Grid_Geometry grid_geometry_;
  IJK_Splitting grid_splitting_;
  int ghost_size_;
  bool perio_k_;

  VECT(ArrOfDouble_with_ghost) local_delta_xyz_;
  bool uniform_[3];

  // Poisson coefficient on each cell on this level
  IJK_Field_float ijk_rho_;
  // matrix coefficients at faces for each level (computed from ijk_rho_levels_)
  // These coefficients contain the wall boundary conditions.
  // The field has 4 components per element
  // coeff_faces(i,j,k,0) contains data for faces of normal x
  // coeff_faces(i,j,k,1) contains data for faces of normal y
  // coeff_faces(i,j,k,2) contains data for faces of normal z
  // coeff_faces(i,j,k,3) contains the sum of the coefficients at all faces of element (i,j,k)
  IJK_Field_float ijk_faces_coefficients_;

  // Temporary storage for multigrid algorithm:
  // Storage for the unknown:
  IJK_Field_float ijk_x_;
  // Storage for the right-hand side:
  IJK_Field_float ijk_rhs_;
  // Storage for the residue:
  IJK_Field_float ijk_residue_;
};

Declare_vect(Grid_Level_Data_float);
#endif
