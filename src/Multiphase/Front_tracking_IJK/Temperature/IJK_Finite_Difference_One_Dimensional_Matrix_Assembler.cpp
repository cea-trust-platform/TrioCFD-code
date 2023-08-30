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
// File      : IJK_Finite_Difference_One_Dimensional_Matrix_Assembler.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Finite_Difference_One_Dimensional_Matrix_Assembler.h>

Implemente_instanciable_sans_constructeur( IJK_Finite_Difference_One_Dimensional_Matrix_Assembler, "IJK_Finite_Difference_One_Dimensional_Matrix_Assembler", Objet_U ) ;

IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::IJK_Finite_Difference_One_Dimensional_Matrix_Assembler()
{
  precision_order_ = 2;
  reduce_side_precision_ = 0;
  core_matrix_type_ = centred;
  // Forward
  set_operators_indices(first_order_derivative_forward_vector_, first_order_derivative_forward_, 0);
  // Centred
  set_operators_indices(first_order_derivative_centred_vector_, first_order_derivative_centred_, 1);
  // Backward
  set_operators_indices(first_order_derivative_backward_vector_, first_order_derivative_backward_, 2);

  // Forward
  set_operators_indices(second_order_derivative_forward_vector_, second_order_derivative_forward_, 0);
  // Centred
  set_operators_indices(second_order_derivative_centred_vector_, second_order_derivative_centred_, 1);
  // Backward
  set_operators_indices(second_order_derivative_backward_vector_, second_order_derivative_backward_, 2);
}

Sortie& IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::readOn( Entree& is )
{
  /*
   * Parse the datafile
   */
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades(is);
  Cout << "IJK_Thermal_base::readOn : Parameters summary. " << finl;
  printOn(Cout);
  return is;
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::set_param(Param& param)
{
  param.ajouter("precision_order", &precision_order_);
  param.ajouter_flag("reduce_side_precision", &reduce_side_precision_);
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::set_operators_size(const std::vector<std::vector<double>>& fd_operator_vector,
                                                                                FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE>& fd_operator)
{
  int i,j;
  for(i=0; i<MAX_ORDER_DERIVATIVE; i++)
    {
      const int fd_operator_size = (int) fd_operator_vector[i].size();
      fd_operator[i][0] = DoubleVect(fd_operator_size);
      fd_operator[i][1] = DoubleVect(fd_operator_size);
      for (j=0; j<fd_operator_size; j++)
        fd_operator[i][1](j) = fd_operator_vector[i][j];
    }
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::set_operators_indices(const std::vector<std::vector<double>>& fd_operator_vector,
                                                                                   FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE>& fd_operator,
                                                                                   const int& fd_coefficient_type)
{
  int i,j;
  set_operators_size(fd_operator_vector, fd_operator);
  switch (fd_coefficient_type)
    {
    case forward:
      for(i=0; i<MAX_ORDER_DERIVATIVE; i++)
        for (j=0; j < fd_operator[i][0].size(); j++)
          fd_operator[i][0](j) = j;
      break;
    case centred:
      for(i=0; i<MAX_ORDER_DERIVATIVE; i++)
        {
          const int fd_operator_size = fd_operator[i][0].size();
          const int offset = (int) (fd_operator_size/2);
          for (j=0; j < fd_operator_size; j++)
            fd_operator[i][0](j) = j-offset;
        }
      break;
    case backward:
      for(i=0; i<MAX_ORDER_DERIVATIVE; i++)
        for (j=0; j < fd_operator[i][0].size(); j++)
          fd_operator[i][0](j) = -j;
      break;
    default:
      for(i=0; i<MAX_ORDER_DERIVATIVE; i++)
        {
          const int fd_operator_size = fd_operator[i][0].size();
          const int offset = (int) (fd_operator_size/2);
          for (j=0; j < fd_operator_size; j++)
            fd_operator[i][0](j) = j-offset;
        }
      break;
    }
}

int IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::non_zero_stencil_values(const FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE>& fd_operator)
{
  int non_zero_elem = 0;
  const int stencil_width = fd_operator[precision_order_-1][1].size();
  for (int i=0; i<stencil_width; i++)
    if (fd_operator[precision_order_-1][1][i] != 0)
      non_zero_elem++;
  return non_zero_elem;
}

/*
 *
 */
int IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::build(Matrice& matrix, const int& nb_elem, const int& derivative_order)
{
  int non_zero_elem = 0;
  int stencil_forward = 0;
  int stencil_centred = 0;
  int stencil_backward = 0;
  FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE> * forward_derivative = nullptr;
  FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE> * centred_derivative = nullptr;
  FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE> * backward_derivative = nullptr;
  switch(derivative_order)
    {
    case first:
      stencil_forward = non_zero_stencil_values(first_order_derivative_forward_);
      stencil_centred = non_zero_stencil_values(first_order_derivative_centred_);
      stencil_backward = non_zero_stencil_values(first_order_derivative_backward_);
      forward_derivative = &first_order_derivative_forward_;
      centred_derivative = &first_order_derivative_centred_;
      backward_derivative = &first_order_derivative_backward_;
      break;
    case second:
      stencil_forward = non_zero_stencil_values(second_order_derivative_forward_);
      stencil_centred = non_zero_stencil_values(second_order_derivative_centred_);
      stencil_backward = non_zero_stencil_values(second_order_derivative_backward_);
      forward_derivative = &second_order_derivative_forward_;
      centred_derivative = &second_order_derivative_centred_;
      backward_derivative = &second_order_derivative_backward_;
      break;
    default:
      stencil_forward = non_zero_stencil_values(first_order_derivative_forward_);
      stencil_centred = non_zero_stencil_values(first_order_derivative_centred_);
      stencil_backward = non_zero_stencil_values(first_order_derivative_backward_);
      forward_derivative = &first_order_derivative_forward_;
      centred_derivative = &first_order_derivative_centred_;
      backward_derivative = &first_order_derivative_backward_;
      break;
    }
  const int core_lines = (nb_elem - 2);
  non_zero_elem = stencil_forward + stencil_backward + stencil_centred * (nb_elem - 2);
  // Cast the finite difference matrix
  matrix.typer("Matrice_Bloc");
  Matrice_Bloc& block_matrix =ref_cast(Matrice_Bloc, matrix.valeur());
  block_matrix.dimensionner(1,1);
  block_matrix.get_bloc(0,0).typer("Matrice_Morse");
  Matrice_Morse& sparse_matrix  = ref_cast(Matrice_Morse, block_matrix.get_bloc(0,0).valeur());
  sparse_matrix.dimensionner(nb_elem, nb_elem, non_zero_elem);

  ArrOfDouble& matrix_values = sparse_matrix.get_set_coeff();
  ArrOfInt& non_zero_coeff_per_line = sparse_matrix.get_set_tab1();
  ArrOfInt& matrix_column_indices = sparse_matrix.get_set_tab2();
  // Fortran start at one
  int non_zero_values_counter = 0;

  /*
   * Ini
   */
  {
    const int forward_derivative_size = (*forward_derivative)[precision_order_-1][0].size();
    for (int i=0; i < forward_derivative_size; i++)
      {
        const int indices = (int) (*forward_derivative)[precision_order_-1][0](i);
        const double fd_coeff = (*forward_derivative)[precision_order_-1][1](i);
        if (fd_coeff != 0)
          {
            matrix_column_indices[non_zero_values_counter] = indices + FORTRAN_INDEX_INI;
            matrix_values[non_zero_values_counter] = fd_coeff;
            non_zero_values_counter++;
          }
      }
    non_zero_coeff_per_line[1] = non_zero_values_counter + FORTRAN_INDEX_INI;
  }
  /*
   * Core
   */
  {
    const int centred_derivative_size = (*centred_derivative)[precision_order_-1][0].size();
    for (int j=1; j<(core_lines+1); j++)
      {
        for (int i=0; i < centred_derivative_size; i++)
          {
            const int indices = (int) (*centred_derivative)[precision_order_-1][0](i);
            const double fd_coeff = (*centred_derivative)[precision_order_-1][1](i);
            if (fd_coeff != 0)
              {
                matrix_column_indices[non_zero_values_counter] = indices + j + FORTRAN_INDEX_INI;
                matrix_values[non_zero_values_counter] = fd_coeff;
                non_zero_values_counter++;
              }
          }
        non_zero_coeff_per_line[j + 1] = non_zero_values_counter + FORTRAN_INDEX_INI;
      }
  }

  /*
   * End
   */
  {
    const int backward_derivative_size = (*backward_derivative)[precision_order_-1][0].size();
    const int last_column = (nb_elem - 1);
    for (int i=0; i < backward_derivative_size; i++)
      {
        const int indices = (int) (*backward_derivative)[precision_order_-1][0](i);
        const double fd_coeff = (*backward_derivative)[precision_order_-1][1](i);
        if (fd_coeff != 0)
          {
            matrix_column_indices[non_zero_values_counter] = last_column + indices + FORTRAN_INDEX_INI;
            matrix_values[non_zero_values_counter] = fd_coeff;
            non_zero_values_counter++;
          }
      }
    non_zero_coeff_per_line[last_column + 1] = non_zero_values_counter + FORTRAN_INDEX_INI;
  }

  /*
   * FIXME : Do three blocks ?
   */
  /*
  block_matrix.dimensionner(3,1);
  block_matrix.get_bloc(0,0) = typer("Matrice_Morse");
  block_matrix.get_bloc(1,0) = typer("Matrice_Morse");
  block_matrix.get_bloc(1,0) = typer("Matrice_Morse_Sym");
  block_matrix.get_bloc(2,0) = typer("Matrice_Morse");
  Matrice_Morse& ini_sparse_matrix  = ref_cast(Matrice_Morse, block_matrix.get_bloc(0,0).valeur());
  Matrice_Morse& end_sparse_matrix  = ref_cast(Matrice_Morse, block_matrix.get_bloc(2,0).valeur());
  Matrice_Morse& core_sparse_matrix;
  Matrice_Morse_Sym& core_sparse_matrix_sym;
  switch(core_matrix_type_)
  	{
  	case centred:
  //			core_sparse_matrix  = ref_cast(Matrice_Morse_Sym, block_matrix.get_bloc(1,0).valeur());
  		core_sparse_matrix  = ref_cast(Matrice_Morse, block_matrix.get_bloc(1,0).valeur());
  		break;
  	case forward:
  		core_sparse_matrix  = ref_cast(Matrice_Morse, block_matrix.get_bloc(1,0).valeur());
  		break;
  	case backward:
  		core_sparse_matrix  = ref_cast(Matrice_Morse, block_matrix.get_bloc(1,0).valeur());
  		break;
  	default:
  //			core_sparse_matrix  = ref_cast(Matrice_Morse_Sym, block_matrix.get_bloc(1,0).valeur());
  		core_sparse_matrix  = ref_cast(Matrice_Morse, block_matrix.get_bloc(1,0).valeur());
  		break;
  	}
  ini_sparse_matrix.dimensionner(1, nb_elem, stencil_forward);
  end_sparse_matrix.dimensionner(nb_elem-2, nb_elem, stencil_backward);
  core_sparse_matrix.dimensionner(1, nb_elem, stencil_centred);
   */

  block_matrix.assert_check_block_matrix_structure();
  return non_zero_elem;
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::sum_matrices_subproblems(Matrice& matrix_A, Matrice& matrix_B)
{
  Matrice_Bloc& block_matrix_A =ref_cast(Matrice_Bloc, matrix_A.valeur());
  Matrice_Bloc& block_matrix_B =ref_cast(Matrice_Bloc, matrix_B.valeur());
  const int dim = block_matrix_A.dim(0);
  for (int i=0; i<dim ; i++)
    {
      Matrice_Morse& sparse_matrix_A  = ref_cast(Matrice_Morse, block_matrix_A.get_bloc(i,i).valeur());
      Matrice_Morse& sparse_matrix_B  = ref_cast(Matrice_Morse, block_matrix_B.get_bloc(i,i).valeur());
      sparse_matrix_A += sparse_matrix_B;
    }
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::sum_matrices(Matrice& matrix_A, Matrice& matrix_B)
{
  Matrice_Bloc& block_matrix_A =ref_cast(Matrice_Bloc, matrix_A.valeur());
  // block_matrix_A.get_bloc(0,0).typer("Matrice_Morse");
  Matrice_Morse& sparse_matrix_A  = ref_cast(Matrice_Morse, block_matrix_A.get_bloc(0,0).valeur());

  Matrice_Bloc& block_matrix_B =ref_cast(Matrice_Bloc, matrix_B.valeur());
  // block_matrix_A.get_bloc(0,0).typer("Matrice_Morse");
  Matrice_Morse& sparse_matrix_B  = ref_cast(Matrice_Morse, block_matrix_B.get_bloc(0,0).valeur());

  sparse_matrix_A += sparse_matrix_B;
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::initialise_matrix_subproblems(Matrice& matrix_subproblems,
                                                                                           Matrice& fd_operator,
                                                                                           const int& nb_subproblems)
{
  matrix_subproblems.typer("Matrice_Bloc");
  Matrice_Bloc& block_matrix_subproblems =ref_cast(Matrice_Bloc, matrix_subproblems.valeur());
  block_matrix_subproblems.dimensionner(nb_subproblems,nb_subproblems);

  Matrice_Bloc& fd_operator_bloc =ref_cast(Matrice_Bloc, fd_operator.valeur());
  Matrice_Morse& fd_operator_morse =ref_cast(Matrice_Morse, fd_operator_bloc.get_bloc(0,0).valeur());

  // const int nb_sub_problems = block_matrix_subproblems.dim(0);
  for (int i=0; i<nb_subproblems; i++)
    {
      block_matrix_subproblems.get_bloc(i,i).typer("Matrice_Morse");
      Matrice_Morse& sparse_matrix  = ref_cast(Matrice_Morse, block_matrix_subproblems.get_bloc(i,i).valeur());
      sparse_matrix = Matrice_Morse(fd_operator_morse);
    }
}
void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::scale_matrix_by_vector(Matrice& matrix,
                                                                                    const DoubleVect& vector,
                                                                                    const int& boundary_conditions)
{
  /*
   * Don't scale the first and last row (where B.Cs are applied)
   */
  DoubleVect vector_tmp = vector;
  if (boundary_conditions)
    {
      vector_tmp[0] = 1.;
      vector_tmp[vector.size() - 1] = 1.;
    }
  /*
   * TODO FIXME: I want a matrix in the end
   */
  // matrix *= vector_tmp;
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::scale_matrix_subproblem_by_vector(Matrice * matrix,
                                                                                               const DoubleVect& vector,
                                                                                               const int& subproblem_index,
                                                                                               const int& boundary_conditions)
{
  Matrice_Bloc& block_matrix =ref_cast(Matrice_Bloc, (*matrix).valeur());
  Matrice& sub_matrix = block_matrix.get_bloc(subproblem_index, subproblem_index);
  scale_matrix_by_vector(sub_matrix, vector, boundary_conditions);
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::scale_matrix_subproblem_by_vector(Matrice& matrix,
                                                                                               DoubleVect& vector,
                                                                                               const int& subproblem_index,
                                                                                               const int& boundary_conditions)
{
  Matrice_Bloc& block_matrix =ref_cast(Matrice_Bloc, matrix.valeur());
  Matrice& sub_matrix = block_matrix.get_bloc(subproblem_index, subproblem_index);
  scale_matrix_by_vector(sub_matrix, vector, boundary_conditions);
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::impose_boundary_conditions(Matrice& modified_matrix,
                                                                                        DoubleVect& modified_rhs,
                                                                                        const int& ini_boundary_conditions,
                                                                                        const double& interfacial_value,
                                                                                        const int& end_boundary_conditions,
                                                                                        const double& end_value)
{
  Matrice matrix = Matrice(modified_matrix);
  DoubleVect rhs = modified_rhs;
  Matrice_Bloc& block_matrix =ref_cast(Matrice_Bloc, matrix.valeur());
  Matrice_Morse& sparse_matrix = ref_cast(Matrice_Morse, block_matrix.get_bloc(0,0).valeur());
  ArrOfDouble& matrix_values = sparse_matrix.get_set_coeff();
  ArrOfInt& non_zero_coeff_per_line = sparse_matrix.get_set_tab1();
  const int nb_lines = sparse_matrix.nb_lignes();
  const int nb_coeff = sparse_matrix.nb_coeff();
  {
    /*
     * Consider a constant interfacial temperature value at saturation
     * A jump condition could be implemented later by equating
     * the fluxes crossing the interface
     */
    switch(ini_boundary_conditions)
      {
      case dirichlet:
        {
          const int non_zero_elem_ini = non_zero_coeff_per_line[1] - FORTRAN_INDEX_INI;
          matrix_values[0] = 1.;
          for (int i=1; i<non_zero_elem_ini; i++)
            matrix_values[i] = 0.;
          rhs[0] = interfacial_value;
        }
        break;
      case neumann:
        rhs[0] = interfacial_value;
        break;
      case flux_jump:
        Cerr << "Flux_jump condition necessitates a sub-problem in each phase : not implemented yet !" << finl;
        break;
      default:
        {
          const int non_zero_elem_ini = non_zero_coeff_per_line[1] - FORTRAN_INDEX_INI;
          matrix_values[0] = 1.;
          for (int i=1; i<non_zero_elem_ini; i++)
            matrix_values[i] = 0.;
          rhs[0] = 0.;
        }
        break;
      }
  }
  /*
   * The probe's end will never be coupled to a sub-resolution in the other phase
   */
  {
    switch(end_boundary_conditions)
      {
      case dirichlet:
        {
          const int non_zero_elem_end = non_zero_coeff_per_line[nb_lines] - FORTRAN_INDEX_INI;
          matrix_values[nb_coeff - 1] = 1.;
          for (int i=nb_lines - 1; i > (nb_lines-non_zero_elem_end); i--)
            matrix_values[i] = 0.;
          rhs[0] = interfacial_value;
        }
        break;
      case neumann:
        break;
      default:
        {
          const int non_zero_elem_end = non_zero_coeff_per_line[nb_lines] - FORTRAN_INDEX_INI;
          matrix_values[nb_coeff - 1] = 1.;
          for (int i=nb_lines - 1; i > (nb_lines-non_zero_elem_end); i--)
            matrix_values[i] = 0.;
          rhs[0] = -1.;
        }
        break;
      }
  }
  sparse_matrix.compacte(remove);
  if ((ini_boundary_conditions == dirichlet) || (end_boundary_conditions == dirichlet))
    modify_rhs_for_bc(matrix,
                      modified_matrix,
                      rhs,
                      modified_rhs,
                      ini_boundary_conditions,
                      end_boundary_conditions);
}

/*
 * Modifcations of the rhs to apply BCs to the resolved 1D field
 */
void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::modify_rhs_for_bc(const Matrice& matrix,
                                                                               Matrice& modified_matrix,
                                                                               const DoubleVect& rhs,
                                                                               DoubleVect& modified_rhs,
                                                                               const int& ini_boundary_conditions,
                                                                               const int& end_boundary_conditions)
{
  /*
   * Copy the matrix into modified matrix
   */
  modified_matrix = Matrice(matrix);
  Matrice_Bloc& block_matrix =ref_cast(Matrice_Bloc, modified_matrix.valeur());
  Matrice_Morse& sparse_matrix = ref_cast(Matrice_Morse, block_matrix.get_bloc(0,0).valeur());
  ArrOfDouble& matrix_values = sparse_matrix.get_set_coeff();
  ArrOfInt& non_zero_coeff_per_line = sparse_matrix.get_set_tab1();
  ArrOfInt& matrix_column_indices = sparse_matrix.get_set_tab2();
  const int nb_lines = sparse_matrix.nb_lignes();
  const int nb_column = sparse_matrix.nb_colonnes();

  /*
   * Get the modified rhs
   */
  modified_rhs = rhs;
  if (ini_boundary_conditions != dirichlet)
    modified_rhs[0] = 0.;
  if (end_boundary_conditions != dirichlet)
    modified_rhs[modified_rhs.size() - 1] = 0.;
  DoubleVect rhs_tmp = modified_rhs;
  modified_rhs = matrix * modified_rhs;
  modified_rhs -= rhs_tmp;
  modified_rhs *= -1.;
  modified_rhs += rhs;

  /*
   * Remove useless coefficients
   */
  for (int j=2; j<nb_lines; j++)
    {
      const int non_zero_elem_core_series = non_zero_coeff_per_line[j] - FORTRAN_INDEX_INI;
      const int non_zero_elem_core = non_zero_coeff_per_line[j] - non_zero_coeff_per_line[j-1];
      const int index_sparse = non_zero_elem_core_series - non_zero_elem_core;
      const int column_index = matrix_column_indices[index_sparse] - FORTRAN_INDEX_INI;
      if (column_index == 0)
        matrix_values[index_sparse] = 0.;
      if ((nb_column - 1) == 0)
        matrix_values[index_sparse] = 0.;
    }

  sparse_matrix.compacte(remove);
}


