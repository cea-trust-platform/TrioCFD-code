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
// File      : Coarsen_Operator_Uniform.cpp
// Directory : $IJK_ROOT/src/IJK/solveur_mg
//
/////////////////////////////////////////////////////////////////////////////
//
// WARNING: DO NOT EDIT THIS FILE! Only edit the template file Coarsen_Operator_Uniform.cpp.P
//
#include <Coarsen_Operator_Uniform.h>
#include <IJK_Grid_Geometry.h>
#include <stat_counters.h> 

Implemente_instanciable_sans_constructeur(Coarsen_Operator_Uniform, "Coarsen_Operator_Uniform", Coarsen_Operator_base);

Coarsen_Operator_Uniform::Coarsen_Operator_Uniform()
{
  coarsen_factors_.resize_array(3);
  // default: uniform refinement by factor 2
  coarsen_factors_ = 2;
}

// Reads coarsening parameters in the .data file. See ajouter_param()
Entree & Coarsen_Operator_Uniform::readOn(Entree & is) 
{ 
  Coarsen_Operator_base::readOn(is); 
  // Check if valid parameters
  bool ok = false;
  if (coarsen_factors_[0] == 1 && coarsen_factors_[1] == 1 && coarsen_factors_[2] == 2) ok = true;
  if (coarsen_factors_[0] == 1 && coarsen_factors_[1] == 2 && coarsen_factors_[2] == 2) ok = true;
  if (coarsen_factors_[0] == 2 && coarsen_factors_[1] == 1 && coarsen_factors_[2] == 2) ok = true;
  if (coarsen_factors_[0] == 2 && coarsen_factors_[1] == 2 && coarsen_factors_[2] == 2) ok = true;
  if (coarsen_factors_[0] == 1 && coarsen_factors_[1] == 1 && coarsen_factors_[2] == 3) ok = true;
  if (coarsen_factors_[0] == 1 && coarsen_factors_[1] == 3 && coarsen_factors_[2] == 2) ok = true;
  if (!ok) {
    Cerr << "Error in Coarsen_Operator_Uniform::readOn: invalid combination of coarsening factors" << finl;
    Process::exit();
  }
  return is;
}

Sortie & Coarsen_Operator_Uniform::printOn(Sortie & os) const { 
  return Coarsen_Operator_base::printOn(os); 
}

void Coarsen_Operator_Uniform::ajouter_param(Param & param)
{
  param.ajouter("coarsen_i", &coarsen_factors_[0]);
  param.ajouter("coarsen_j", &coarsen_factors_[1]);
  param.ajouter("coarsen_k", &coarsen_factors_[2]);
  Coarsen_Operator_base::ajouter_param(param);
}

void Coarsen_Operator_Uniform::initialize_grid_data(const Grid_Level_Data_float & fine, Grid_Level_Data_float & coarse,
						    int additional_k_layers)
{
  const IJK_Grid_Geometry & src_grid_geom = fine.get_grid_geometry();
  VECT(ArrOfDouble) coarse_delta(3);
  ArrOfInt nlocal(3);

  for (int dir = 0; dir < 3; dir++) {
    const int coarsen_factor = coarsen_factors_[dir];
    const ArrOfDouble & fine_delta = src_grid_geom.get_delta(dir);
    
    const int src_n = fine_delta.size_array();
    if (src_n % coarsen_factor != 0) {
      Cerr << "Coarsen_Operator_Uniform::initialize_grid_data: source grid has " << src_n 
	   << " elements in direction " << dir 
	   << " and cannot be refined by factor " << coarsen_factor << finl;
      Process::exit();
    }
    nlocal[dir] = fine.get_rho().nb_elem_local(dir);
    if (nlocal[dir] % coarsen_factor != 0) {
      Cerr << "Coarsen_Operator_Uniform::initialize_grid_data: local source grid processor has " << nlocal[dir]
	   << " elements in direction " << dir 
	   << " and cannot be refined by factor " << coarsen_factor << finl;
      Process::exit();
    }
    nlocal[dir] /= coarsen_factor;
    const int n = src_n / coarsen_factor;
    ArrOfDouble & delta = coarse_delta[dir];
    delta.resize_array(n);
    
    int src_index = 0;
    delta.resize_array(n);
    for (int dest_index = 0; dest_index < n; dest_index++) {
      double d = 0;
      for (int j = 0; j < coarsen_factor; j++) {
	d += fine_delta[src_index];
	src_index++;
      }
      delta[dest_index] = d;
    }
  }

  IJK_Grid_Geometry grid_geom;
  grid_geom.initialize_origin_deltas(src_grid_geom.get_origin(0),
				     src_grid_geom.get_origin(1),
				     src_grid_geom.get_origin(2),
				     coarse_delta[0],
				     coarse_delta[1],
				     coarse_delta[2],
				     src_grid_geom.get_periodic_flag(0),
				     src_grid_geom.get_periodic_flag(1),
				     src_grid_geom.get_periodic_flag(2));

  IJK_Splitting coarse_splitting;
  // Same processor mapping as fine mesh
  IntTab processor_mapping;
  fine.get_splitting().get_processor_mapping(processor_mapping);
  // Splitting is identical, divide ncells by the coarsening factor
  VECT(ArrOfInt) slice_sizes(3);
  for (int dir = 0; dir < 3; dir++) {
    fine.get_splitting().get_slice_size(dir, IJK_Splitting::ELEM, slice_sizes[dir]);
    const int n = slice_sizes[dir].size_array();
    for (int i = 0; i < n; i++) 
      slice_sizes[dir][i] /= coarsen_factors_[dir];
  }  
  coarse_splitting.initialize(grid_geom, slice_sizes[0], slice_sizes[1], slice_sizes[2], processor_mapping);
  const int ghost_zone_size = fine.get_ghost_size();
  coarse.initialize(coarse_splitting, ghost_zone_size, additional_k_layers);
}


void Coarsen_Operator_Uniform::coarsen(const IJK_Field_float & fine, IJK_Field_float & coarse, int compute_weighted_average) const
{
  static Stat_Counter_Id coarsen_counter_ = statistiques().new_counter(2, "multigrille : uniform coarsen ");
  statistiques().begin_count(coarsen_counter_);

  const int ni2 = coarse.ni();
  const int nj2 = coarse.nj();
  const int nk2 = coarse.nk();

  float coef;
  if (compute_weighted_average)
    coef = 1. / ((float) coarsen_factors_[0] * coarsen_factors_[1] * coarsen_factors_[2]);
  else
    coef = 1.;

 
 if (coarsen_factors_[2] == 2) {
    if (coarsen_factors_[1] == 2) {
      if (coarsen_factors_[0] == 2) {
          const int Kstart = 0;
  const int Kend = nk2;
  const int deltaK = 1;
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*2;
	int j = J*2;
	int k = K*2;
	float sum = 0.;
	for (int ii = 0; ii < 2; ii++)
	  for (int jj = 0; jj < 2; jj++)
	    for (int kk = 0; kk < 2; kk++)
	      sum += fine(i + ii, j + jj, k + kk);
	coarse(I,J,K) = sum * coef;
      }
    }
  }
;
      } else {
	  const int Kstart = 0;
  const int Kend = nk2;
  const int deltaK = 1;
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*1;
	int j = J*2;
	int k = K*2;
	float sum = 0.;
	for (int ii = 0; ii < 1; ii++)
	  for (int jj = 0; jj < 2; jj++)
	    for (int kk = 0; kk < 2; kk++)
	      sum += fine(i + ii, j + jj, k + kk);
	coarse(I,J,K) = sum * coef;
      }
    }
  }
;
      }
    } else {
      if (coarsen_factors_[0] == 2) {
          const int Kstart = 0;
  const int Kend = nk2;
  const int deltaK = 1;
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*2;
	int j = J*1;
	int k = K*2;
	float sum = 0.;
	for (int ii = 0; ii < 2; ii++)
	  for (int jj = 0; jj < 1; jj++)
	    for (int kk = 0; kk < 2; kk++)
	      sum += fine(i + ii, j + jj, k + kk);
	coarse(I,J,K) = sum * coef;
      }
    }
  }
      } else {
	  const int Kstart = 0;
  const int Kend = nk2;
  const int deltaK = 1;
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*1;
	int j = J*1;
	int k = K*2;
	float sum = 0.;
	for (int ii = 0; ii < 1; ii++)
	  for (int jj = 0; jj < 1; jj++)
	    for (int kk = 0; kk < 2; kk++)
	      sum += fine(i + ii, j + jj, k + kk);
	coarse(I,J,K) = sum * coef;
      }
    }
  }
;
      }
    }
 } else if (coarsen_factors_[2] == 3
	    && coarsen_factors_[1] == 1
	    && coarsen_factors_[0] == 1) {
          const int Kstart = 0;
  const int Kend = nk2;
  const int deltaK = 1;
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*1;
	int j = J*1;
	int k = K*3;
	float sum = 0.;
	for (int ii = 0; ii < 1; ii++)
	  for (int jj = 0; jj < 1; jj++)
	    for (int kk = 0; kk < 3; kk++)
	      sum += fine(i + ii, j + jj, k + kk);
	coarse(I,J,K) = sum * coef;
      }
    }
  }
;
 } else if (coarsen_factors_[2] == 2
            && coarsen_factors_[1] == 3
            && coarsen_factors_[0] == 1) {
          const int Kstart = 0;
  const int Kend = nk2;
  const int deltaK = 1;
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*1;
	int j = J*3;
	int k = K*2;
	float sum = 0.;
	for (int ii = 0; ii < 1; ii++)
	  for (int jj = 0; jj < 3; jj++)
	    for (int kk = 0; kk < 2; kk++)
	      sum += fine(i + ii, j + jj, k + kk);
	coarse(I,J,K) = sum * coef;
      }
    }
  }
;
 } else {
   
    Process::exit();
  }
  statistiques().end_count(coarsen_counter_);

}

void Coarsen_Operator_Uniform::interpolate_sub_shiftk(const IJK_Field_float & coarse, IJK_Field_float & fine, const int kshift) const
{
  static Stat_Counter_Id interpolate_counter_ = statistiques().new_counter(2, "multigrille : interpolate (uniform)");
  statistiques().begin_count(interpolate_counter_);

#if 0
  const int ni = fine.ni();
  const int nj = fine.nj();
  const int nk = fine.nk();

  const int ni2 = coarse.ni();
  const int nj2 = coarse.nj();
  const int nk2 = coarse.nk();
  
  int Kstart, Kend, deltaK;
  if (kshift <= 0) {
    Kstart = 0;
    Kend = nk2 - 1;
    deltaK = 1;
  } else {
    Kstart = nk2 - 1;
    Kend = 0;
    deltaK = -1;
  }
  if (ni2 * 2 == ni) {
    for (int K = Kstart; K != Kend; K += deltaK) {
      for (int J = 0; J < nj2; J++) {
	for (int I = 0; I < ni2; I++) {
	  int i = I*2;
	  int j = J*2;
	  int k = K*2;
	  const double val = coarse(I,J,K);
	  for (int dk = 0; dk < 2; dk++)
	    for (int dj = 0; dj < 2; dj++)
	      for (int di = 0; di < 2; di++)
		//fine.get_in_allocated_area(i+di,j+dj,k+dk + kshift) = fine(i+di,j+dj,k+dk) - val;
		fine(i+di,j+dj,k+dk) -= val;
	}
      }
    }
  } else {
    assert(0);
    Process::exit();
  }
  fine.shift_k_origin(kshift);
  flop_count += ni*nj*nk;
#else
  const int ni2 = coarse.ni();
  const int nj2 = coarse.nj();
  const int nk2 = coarse.nk();
 
 if (coarsen_factors_[2] == 2) {
    if (coarsen_factors_[1] == 2) {
      if (coarsen_factors_[0] == 2) {
        if (kshift == 1) {
  Cerr << "error Coarsen_Operator_Uniform: kshift=1 will not work" << finl;
  Process::exit();
 }
  int Kstart, Kend, deltaK;
  if (kshift <= 0) {
    Kstart = 0;
    Kend = nk2;
    deltaK = 1;
  } else {
    Kstart = nk2 - 1;
    Kend = -1;
    deltaK = -1;
  }
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*2;
	int j = J*2;
	int k = K*2;
	float val = coarse(I,J,K);
	for (int ii = 0; ii < 2; ii++)
	  for (int jj = 0; jj < 2; jj++)
	    for (int kk = 0; kk < 2; kk++) // this loop must be reversed if kshift=1 !
	      fine.get_in_allocated_area(i + ii, j + jj, k + kk + kshift) = fine(i + ii, j + jj, k + kk) - val;	
      }
    }
  }
;
      } else {
	if (kshift == 1) {
  Cerr << "error Coarsen_Operator_Uniform: kshift=1 will not work" << finl;
  Process::exit();
 }
  int Kstart, Kend, deltaK;
  if (kshift <= 0) {
    Kstart = 0;
    Kend = nk2;
    deltaK = 1;
  } else {
    Kstart = nk2 - 1;
    Kend = -1;
    deltaK = -1;
  }
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*1;
	int j = J*2;
	int k = K*2;
	float val = coarse(I,J,K);
	for (int ii = 0; ii < 1; ii++)
	  for (int jj = 0; jj < 2; jj++)
	    for (int kk = 0; kk < 2; kk++) // this loop must be reversed if kshift=1 !
	      fine.get_in_allocated_area(i + ii, j + jj, k + kk + kshift) = fine(i + ii, j + jj, k + kk) - val;	
      }
    }
  }
;
      }
    } else {
      if (coarsen_factors_[0] == 2) {
        if (kshift == 1) {
  Cerr << "error Coarsen_Operator_Uniform: kshift=1 will not work" << finl;
  Process::exit();
 }
  int Kstart, Kend, deltaK;
  if (kshift <= 0) {
    Kstart = 0;
    Kend = nk2;
    deltaK = 1;
  } else {
    Kstart = nk2 - 1;
    Kend = -1;
    deltaK = -1;
  }
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*2;
	int j = J*1;
	int k = K*2;
	float val = coarse(I,J,K);
	for (int ii = 0; ii < 2; ii++)
	  for (int jj = 0; jj < 1; jj++)
	    for (int kk = 0; kk < 2; kk++) // this loop must be reversed if kshift=1 !
	      fine.get_in_allocated_area(i + ii, j + jj, k + kk + kshift) = fine(i + ii, j + jj, k + kk) - val;	
      }
    }
  }
      } else {
	if (kshift == 1) {
  Cerr << "error Coarsen_Operator_Uniform: kshift=1 will not work" << finl;
  Process::exit();
 }
  int Kstart, Kend, deltaK;
  if (kshift <= 0) {
    Kstart = 0;
    Kend = nk2;
    deltaK = 1;
  } else {
    Kstart = nk2 - 1;
    Kend = -1;
    deltaK = -1;
  }
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*1;
	int j = J*1;
	int k = K*2;
	float val = coarse(I,J,K);
	for (int ii = 0; ii < 1; ii++)
	  for (int jj = 0; jj < 1; jj++)
	    for (int kk = 0; kk < 2; kk++) // this loop must be reversed if kshift=1 !
	      fine.get_in_allocated_area(i + ii, j + jj, k + kk + kshift) = fine(i + ii, j + jj, k + kk) - val;	
      }
    }
  }
;
      }
    }
 } else if (coarsen_factors_[2] == 3
	    && coarsen_factors_[1] == 1
	    && coarsen_factors_[0] == 1) {
        if (kshift == 1) {
  Cerr << "error Coarsen_Operator_Uniform: kshift=1 will not work" << finl;
  Process::exit();
 }
  int Kstart, Kend, deltaK;
  if (kshift <= 0) {
    Kstart = 0;
    Kend = nk2;
    deltaK = 1;
  } else {
    Kstart = nk2 - 1;
    Kend = -1;
    deltaK = -1;
  }
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*1;
	int j = J*1;
	int k = K*3;
	float val = coarse(I,J,K);
	for (int ii = 0; ii < 1; ii++)
	  for (int jj = 0; jj < 1; jj++)
	    for (int kk = 0; kk < 3; kk++) // this loop must be reversed if kshift=1 !
	      fine.get_in_allocated_area(i + ii, j + jj, k + kk + kshift) = fine(i + ii, j + jj, k + kk) - val;	
      }
    }
  }
;
 } else if (coarsen_factors_[2] == 2
            && coarsen_factors_[1] == 3
            && coarsen_factors_[0] == 1) {
        if (kshift == 1) {
  Cerr << "error Coarsen_Operator_Uniform: kshift=1 will not work" << finl;
  Process::exit();
 }
  int Kstart, Kend, deltaK;
  if (kshift <= 0) {
    Kstart = 0;
    Kend = nk2;
    deltaK = 1;
  } else {
    Kstart = nk2 - 1;
    Kend = -1;
    deltaK = -1;
  }
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*1;
	int j = J*3;
	int k = K*2;
	float val = coarse(I,J,K);
	for (int ii = 0; ii < 1; ii++)
	  for (int jj = 0; jj < 3; jj++)
	    for (int kk = 0; kk < 2; kk++) // this loop must be reversed if kshift=1 !
	      fine.get_in_allocated_area(i + ii, j + jj, k + kk + kshift) = fine(i + ii, j + jj, k + kk) - val;	
      }
    }
  }
;
 } else {
   
    Process::exit();
  }
  fine.shift_k_origin(kshift);
#endif
  statistiques().end_count(interpolate_counter_);

}
void Coarsen_Operator_Uniform::initialize_grid_data(const Grid_Level_Data_double & fine, Grid_Level_Data_double & coarse,
						    int additional_k_layers)
{
  const IJK_Grid_Geometry & src_grid_geom = fine.get_grid_geometry();
  VECT(ArrOfDouble) coarse_delta(3);
  ArrOfInt nlocal(3);

  for (int dir = 0; dir < 3; dir++) {
    const int coarsen_factor = coarsen_factors_[dir];
    const ArrOfDouble & fine_delta = src_grid_geom.get_delta(dir);
    
    const int src_n = fine_delta.size_array();
    if (src_n % coarsen_factor != 0) {
      Cerr << "Coarsen_Operator_Uniform::initialize_grid_data: source grid has " << src_n 
	   << " elements in direction " << dir 
	   << " and cannot be refined by factor " << coarsen_factor << finl;
      Process::exit();
    }
    nlocal[dir] = fine.get_rho().nb_elem_local(dir);
    if (nlocal[dir] % coarsen_factor != 0) {
      Cerr << "Coarsen_Operator_Uniform::initialize_grid_data: local source grid processor has " << nlocal[dir]
	   << " elements in direction " << dir 
	   << " and cannot be refined by factor " << coarsen_factor << finl;
      Process::exit();
    }
    nlocal[dir] /= coarsen_factor;
    const int n = src_n / coarsen_factor;
    ArrOfDouble & delta = coarse_delta[dir];
    delta.resize_array(n);
    
    int src_index = 0;
    delta.resize_array(n);
    for (int dest_index = 0; dest_index < n; dest_index++) {
      double d = 0;
      for (int j = 0; j < coarsen_factor; j++) {
	d += fine_delta[src_index];
	src_index++;
      }
      delta[dest_index] = d;
    }
  }

  IJK_Grid_Geometry grid_geom;
  grid_geom.initialize_origin_deltas(src_grid_geom.get_origin(0),
				     src_grid_geom.get_origin(1),
				     src_grid_geom.get_origin(2),
				     coarse_delta[0],
				     coarse_delta[1],
				     coarse_delta[2],
				     src_grid_geom.get_periodic_flag(0),
				     src_grid_geom.get_periodic_flag(1),
				     src_grid_geom.get_periodic_flag(2));

  IJK_Splitting coarse_splitting;
  // Same processor mapping as fine mesh
  IntTab processor_mapping;
  fine.get_splitting().get_processor_mapping(processor_mapping);
  // Splitting is identical, divide ncells by the coarsening factor
  VECT(ArrOfInt) slice_sizes(3);
  for (int dir = 0; dir < 3; dir++) {
    fine.get_splitting().get_slice_size(dir, IJK_Splitting::ELEM, slice_sizes[dir]);
    const int n = slice_sizes[dir].size_array();
    for (int i = 0; i < n; i++) 
      slice_sizes[dir][i] /= coarsen_factors_[dir];
  }  
  coarse_splitting.initialize(grid_geom, slice_sizes[0], slice_sizes[1], slice_sizes[2], processor_mapping);
  const int ghost_zone_size = fine.get_ghost_size();
  coarse.initialize(coarse_splitting, ghost_zone_size, additional_k_layers);
}


void Coarsen_Operator_Uniform::coarsen(const IJK_Field_double & fine, IJK_Field_double & coarse, int compute_weighted_average) const
{
  static Stat_Counter_Id coarsen_counter_ = statistiques().new_counter(2, "multigrille : uniform coarsen ");
  statistiques().begin_count(coarsen_counter_);

  const int ni2 = coarse.ni();
  const int nj2 = coarse.nj();
  const int nk2 = coarse.nk();

  double coef;
  if (compute_weighted_average)
    coef = 1. / ((double) coarsen_factors_[0] * coarsen_factors_[1] * coarsen_factors_[2]);
  else
    coef = 1.;

 
 if (coarsen_factors_[2] == 2) {
    if (coarsen_factors_[1] == 2) {
      if (coarsen_factors_[0] == 2) {
          const int Kstart = 0;
  const int Kend = nk2;
  const int deltaK = 1;
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*2;
	int j = J*2;
	int k = K*2;
	double sum = 0.;
	for (int ii = 0; ii < 2; ii++)
	  for (int jj = 0; jj < 2; jj++)
	    for (int kk = 0; kk < 2; kk++)
	      sum += fine(i + ii, j + jj, k + kk);
	coarse(I,J,K) = sum * coef;
      }
    }
  }
;
      } else {
	  const int Kstart = 0;
  const int Kend = nk2;
  const int deltaK = 1;
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*1;
	int j = J*2;
	int k = K*2;
	double sum = 0.;
	for (int ii = 0; ii < 1; ii++)
	  for (int jj = 0; jj < 2; jj++)
	    for (int kk = 0; kk < 2; kk++)
	      sum += fine(i + ii, j + jj, k + kk);
	coarse(I,J,K) = sum * coef;
      }
    }
  }
;
      }
    } else {
      if (coarsen_factors_[0] == 2) {
          const int Kstart = 0;
  const int Kend = nk2;
  const int deltaK = 1;
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*2;
	int j = J*1;
	int k = K*2;
	double sum = 0.;
	for (int ii = 0; ii < 2; ii++)
	  for (int jj = 0; jj < 1; jj++)
	    for (int kk = 0; kk < 2; kk++)
	      sum += fine(i + ii, j + jj, k + kk);
	coarse(I,J,K) = sum * coef;
      }
    }
  }
      } else {
	  const int Kstart = 0;
  const int Kend = nk2;
  const int deltaK = 1;
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*1;
	int j = J*1;
	int k = K*2;
	double sum = 0.;
	for (int ii = 0; ii < 1; ii++)
	  for (int jj = 0; jj < 1; jj++)
	    for (int kk = 0; kk < 2; kk++)
	      sum += fine(i + ii, j + jj, k + kk);
	coarse(I,J,K) = sum * coef;
      }
    }
  }
;
      }
    }
 } else if (coarsen_factors_[2] == 3
	    && coarsen_factors_[1] == 1
	    && coarsen_factors_[0] == 1) {
          const int Kstart = 0;
  const int Kend = nk2;
  const int deltaK = 1;
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*1;
	int j = J*1;
	int k = K*3;
	double sum = 0.;
	for (int ii = 0; ii < 1; ii++)
	  for (int jj = 0; jj < 1; jj++)
	    for (int kk = 0; kk < 3; kk++)
	      sum += fine(i + ii, j + jj, k + kk);
	coarse(I,J,K) = sum * coef;
      }
    }
  }
;
 } else if (coarsen_factors_[2] == 2
            && coarsen_factors_[1] == 3
            && coarsen_factors_[0] == 1) {
          const int Kstart = 0;
  const int Kend = nk2;
  const int deltaK = 1;
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*1;
	int j = J*3;
	int k = K*2;
	double sum = 0.;
	for (int ii = 0; ii < 1; ii++)
	  for (int jj = 0; jj < 3; jj++)
	    for (int kk = 0; kk < 2; kk++)
	      sum += fine(i + ii, j + jj, k + kk);
	coarse(I,J,K) = sum * coef;
      }
    }
  }
;
 } else {
   
    Process::exit();
  }
  statistiques().end_count(coarsen_counter_);

}

void Coarsen_Operator_Uniform::interpolate_sub_shiftk(const IJK_Field_double & coarse, IJK_Field_double & fine, const int kshift) const
{
  static Stat_Counter_Id interpolate_counter_ = statistiques().new_counter(2, "multigrille : interpolate (uniform)");
  statistiques().begin_count(interpolate_counter_);

#if 0
  const int ni = fine.ni();
  const int nj = fine.nj();
  const int nk = fine.nk();

  const int ni2 = coarse.ni();
  const int nj2 = coarse.nj();
  const int nk2 = coarse.nk();
  
  int Kstart, Kend, deltaK;
  if (kshift <= 0) {
    Kstart = 0;
    Kend = nk2 - 1;
    deltaK = 1;
  } else {
    Kstart = nk2 - 1;
    Kend = 0;
    deltaK = -1;
  }
  if (ni2 * 2 == ni) {
    for (int K = Kstart; K != Kend; K += deltaK) {
      for (int J = 0; J < nj2; J++) {
	for (int I = 0; I < ni2; I++) {
	  int i = I*2;
	  int j = J*2;
	  int k = K*2;
	  const double val = coarse(I,J,K);
	  for (int dk = 0; dk < 2; dk++)
	    for (int dj = 0; dj < 2; dj++)
	      for (int di = 0; di < 2; di++)
		//fine.get_in_allocated_area(i+di,j+dj,k+dk + kshift) = fine(i+di,j+dj,k+dk) - val;
		fine(i+di,j+dj,k+dk) -= val;
	}
      }
    }
  } else {
    assert(0);
    Process::exit();
  }
  fine.shift_k_origin(kshift);
  flop_count += ni*nj*nk;
#else
  const int ni2 = coarse.ni();
  const int nj2 = coarse.nj();
  const int nk2 = coarse.nk();
 
 if (coarsen_factors_[2] == 2) {
    if (coarsen_factors_[1] == 2) {
      if (coarsen_factors_[0] == 2) {
        if (kshift == 1) {
  Cerr << "error Coarsen_Operator_Uniform: kshift=1 will not work" << finl;
  Process::exit();
 }
  int Kstart, Kend, deltaK;
  if (kshift <= 0) {
    Kstart = 0;
    Kend = nk2;
    deltaK = 1;
  } else {
    Kstart = nk2 - 1;
    Kend = -1;
    deltaK = -1;
  }
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*2;
	int j = J*2;
	int k = K*2;
	double val = coarse(I,J,K);
	for (int ii = 0; ii < 2; ii++)
	  for (int jj = 0; jj < 2; jj++)
	    for (int kk = 0; kk < 2; kk++) // this loop must be reversed if kshift=1 !
	      fine.get_in_allocated_area(i + ii, j + jj, k + kk + kshift) = fine(i + ii, j + jj, k + kk) - val;	
      }
    }
  }
;
      } else {
	if (kshift == 1) {
  Cerr << "error Coarsen_Operator_Uniform: kshift=1 will not work" << finl;
  Process::exit();
 }
  int Kstart, Kend, deltaK;
  if (kshift <= 0) {
    Kstart = 0;
    Kend = nk2;
    deltaK = 1;
  } else {
    Kstart = nk2 - 1;
    Kend = -1;
    deltaK = -1;
  }
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*1;
	int j = J*2;
	int k = K*2;
	double val = coarse(I,J,K);
	for (int ii = 0; ii < 1; ii++)
	  for (int jj = 0; jj < 2; jj++)
	    for (int kk = 0; kk < 2; kk++) // this loop must be reversed if kshift=1 !
	      fine.get_in_allocated_area(i + ii, j + jj, k + kk + kshift) = fine(i + ii, j + jj, k + kk) - val;	
      }
    }
  }
;
      }
    } else {
      if (coarsen_factors_[0] == 2) {
        if (kshift == 1) {
  Cerr << "error Coarsen_Operator_Uniform: kshift=1 will not work" << finl;
  Process::exit();
 }
  int Kstart, Kend, deltaK;
  if (kshift <= 0) {
    Kstart = 0;
    Kend = nk2;
    deltaK = 1;
  } else {
    Kstart = nk2 - 1;
    Kend = -1;
    deltaK = -1;
  }
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*2;
	int j = J*1;
	int k = K*2;
	double val = coarse(I,J,K);
	for (int ii = 0; ii < 2; ii++)
	  for (int jj = 0; jj < 1; jj++)
	    for (int kk = 0; kk < 2; kk++) // this loop must be reversed if kshift=1 !
	      fine.get_in_allocated_area(i + ii, j + jj, k + kk + kshift) = fine(i + ii, j + jj, k + kk) - val;	
      }
    }
  }
      } else {
	if (kshift == 1) {
  Cerr << "error Coarsen_Operator_Uniform: kshift=1 will not work" << finl;
  Process::exit();
 }
  int Kstart, Kend, deltaK;
  if (kshift <= 0) {
    Kstart = 0;
    Kend = nk2;
    deltaK = 1;
  } else {
    Kstart = nk2 - 1;
    Kend = -1;
    deltaK = -1;
  }
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*1;
	int j = J*1;
	int k = K*2;
	double val = coarse(I,J,K);
	for (int ii = 0; ii < 1; ii++)
	  for (int jj = 0; jj < 1; jj++)
	    for (int kk = 0; kk < 2; kk++) // this loop must be reversed if kshift=1 !
	      fine.get_in_allocated_area(i + ii, j + jj, k + kk + kshift) = fine(i + ii, j + jj, k + kk) - val;	
      }
    }
  }
;
      }
    }
 } else if (coarsen_factors_[2] == 3
	    && coarsen_factors_[1] == 1
	    && coarsen_factors_[0] == 1) {
        if (kshift == 1) {
  Cerr << "error Coarsen_Operator_Uniform: kshift=1 will not work" << finl;
  Process::exit();
 }
  int Kstart, Kend, deltaK;
  if (kshift <= 0) {
    Kstart = 0;
    Kend = nk2;
    deltaK = 1;
  } else {
    Kstart = nk2 - 1;
    Kend = -1;
    deltaK = -1;
  }
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*1;
	int j = J*1;
	int k = K*3;
	double val = coarse(I,J,K);
	for (int ii = 0; ii < 1; ii++)
	  for (int jj = 0; jj < 1; jj++)
	    for (int kk = 0; kk < 3; kk++) // this loop must be reversed if kshift=1 !
	      fine.get_in_allocated_area(i + ii, j + jj, k + kk + kshift) = fine(i + ii, j + jj, k + kk) - val;	
      }
    }
  }
;
 } else if (coarsen_factors_[2] == 2
            && coarsen_factors_[1] == 3
            && coarsen_factors_[0] == 1) {
        if (kshift == 1) {
  Cerr << "error Coarsen_Operator_Uniform: kshift=1 will not work" << finl;
  Process::exit();
 }
  int Kstart, Kend, deltaK;
  if (kshift <= 0) {
    Kstart = 0;
    Kend = nk2;
    deltaK = 1;
  } else {
    Kstart = nk2 - 1;
    Kend = -1;
    deltaK = -1;
  }
  for (int K = Kstart; K != Kend; K += deltaK) {
    for (int J = 0; J < nj2; J++) {
      for (int I = 0; I < ni2; I++) {
	int i = I*1;
	int j = J*3;
	int k = K*2;
	double val = coarse(I,J,K);
	for (int ii = 0; ii < 1; ii++)
	  for (int jj = 0; jj < 3; jj++)
	    for (int kk = 0; kk < 2; kk++) // this loop must be reversed if kshift=1 !
	      fine.get_in_allocated_area(i + ii, j + jj, k + kk + kshift) = fine(i + ii, j + jj, k + kk) - val;	
      }
    }
  }
;
 } else {
   
    Process::exit();
  }
  fine.shift_k_origin(kshift);
#endif
  statistiques().end_count(interpolate_counter_);

}

