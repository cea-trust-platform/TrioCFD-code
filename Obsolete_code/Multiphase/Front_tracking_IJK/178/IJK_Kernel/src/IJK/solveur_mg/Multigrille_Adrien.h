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
// File      : Multigrille_Adrien.h
// Directory : $IJK_ROOT/src/IJK/solveur_mg
//
/////////////////////////////////////////////////////////////////////////////

//
// WARNING: DO NOT EDIT THIS FILE! Only edit the template file Multigrille_Adrien.h.P
//
#ifndef Multigrille_Adrien_H
#define Multigrille_Adrien_H
#include <Multigrille_base.h>
#include <Equation_base.h>
#include <ArrOfInt.h>
#include <Grid_Level_Data.h>
#include <Coarsen_Operator_base.h>
class IJK_Splitting;
class Multigrille_Adrien : public Multigrille_base
{
  Declare_instanciable_sans_constructeur(Multigrille_Adrien);
public:
  Multigrille_Adrien();
  int completer(const Equation_base& eq);
  void initialize(const IJK_Splitting&);
  int nb_grid_levels() const
  {
    return coarsen_operators_.size() + 1;
  }

  void set_rho(const DoubleVect& rho);
  void set_rho(const IJK_Field_double& rho);
  void set_rho(const IJK_Field_float& rho);
  void reset_rho();
  void set_inv_rho(const IJK_Field_float& inv_rho);
  void prepare_secmem(IJK_Field_float& x) const;
  void dump_lata(const Nom& field, const IJK_Field_float& data, int tstep) const;
  void set_inv_rho_float(const IJK_Field_float& inv_rho, bool set_coarse_matrix_flag, bool use_coeffs_from_double);
  void set_inv_rho_float(const IJK_Field_double& inv_rho, bool set_coarse_matrix_flag, bool use_coeffs_from_double);
  void set_inv_rho(const IJK_Field_double& inv_rho);
  void prepare_secmem(IJK_Field_double& x) const;
  void dump_lata(const Nom& field, const IJK_Field_double& data, int tstep) const;
  void set_inv_rho_double(const IJK_Field_float& inv_rho, bool set_coarse_matrix_flag, bool use_coeffs_from_double);
  void set_inv_rho_double(const IJK_Field_double& inv_rho, bool set_coarse_matrix_flag, bool use_coeffs_from_double);
protected:
  void ajouter_param(Param& param);
  int needed_kshift_for_jacobi(int level) const;
  void completer_float(const IJK_Splitting&);
  void alloc_field(IJK_Field_float& x, int level, bool with_additional_k_layers = false) const;
  void jacobi_residu(IJK_Field_float& x,
                     const IJK_Field_float *secmem, /* if null pointer, take secmem = 0 (to compute A*x) */
                     const int grid_level,
                     const int n_jacobi,
                     IJK_Field_float *residu) const;

  void coarsen(const IJK_Field_float& fine, IJK_Field_float& coarse, int fine_level) const;
  void interpolate_sub_shiftk(const IJK_Field_float& coarse, IJK_Field_float& fine, int fine_level) const;
  void setup_coarse_grid(const IJK_Field_float& fine, IJK_Field_float& coarse) const;
  void compute_coefficients(IJK_Field_float& coeffs_faces,
                            const int grid_level) const;
  void set_rho_float(const IJK_Field_float& rho, bool set_coarse_matrix, bool use_coeffs_from_double);
  void set_rho_float(const IJK_Field_double& rho, bool set_coarse_matrix, bool use_coeffs_from_double);
  virtual IJK_Field_float& get_storage_float(StorageId, int level);
  void completer_double(const IJK_Splitting&);
  void alloc_field(IJK_Field_double& x, int level, bool with_additional_k_layers = false) const;
  void jacobi_residu(IJK_Field_double& x,
                     const IJK_Field_double *secmem, /* if null pointer, take secmem = 0 (to compute A*x) */
                     const int grid_level,
                     const int n_jacobi,
                     IJK_Field_double *residu) const;

  void coarsen(const IJK_Field_double& fine, IJK_Field_double& coarse, int fine_level) const;
  void interpolate_sub_shiftk(const IJK_Field_double& coarse, IJK_Field_double& fine, int fine_level) const;
  void setup_coarse_grid(const IJK_Field_double& fine, IJK_Field_double& coarse) const;
  void compute_coefficients(IJK_Field_double& coeffs_faces,
                            const int grid_level) const;
  void set_rho_double(const IJK_Field_float& rho, bool set_coarse_matrix, bool use_coeffs_from_double);
  void set_rho_double(const IJK_Field_double& rho, bool set_coarse_matrix, bool use_coeffs_from_double);
  virtual IJK_Field_double& get_storage_double(StorageId, int level);
  void completer_double_for_residue(const IJK_Splitting& splitting);

private:

  // number of isotropic coarsening steps
  int nb_isotropic_coarsening_;
  // Thickness of the ghost zone: determines the maximum number of sweeps of the Jacobi
  // smoother that are performed in the same pass (1, 2, 4 or 6 can be efficient depending
  // on the architecture and the problem)
  int ghost_size_;

  // Coarsening operators (read in the .data file)
  VECT(DERIV(Coarsen_Operator_base)) coarsen_operators_;

  // Data for each grid (number of grids = number of coarsening operators + 1)
  VECT(Grid_Level_Data_float) grids_data_float_;
  VECT(Grid_Level_Data_double) grids_data_double_;
};


#endif
