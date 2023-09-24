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

// static IntVect arg_sort(DoubleVect matrix_column_indices)
static std::vector<int> arg_sort(DoubleVect matrix_column_indices)
{
  const int n = matrix_column_indices.size();
  // IntVect indices(n);
  std::vector<int> indices(n);
  for (int j=0; j<n; j++)
    indices[j]=j;
  std::sort(indices.begin(), indices.end(), [&matrix_column_indices](int i, int j) {return matrix_column_indices[i] < matrix_column_indices[j];});
  return indices;
}

static void sort_stencil(Matrice_Morse& sparse_matrix)
{
  DoubleVect& matrix_values = sparse_matrix.get_set_coeff();
  IntVect& non_zero_coeff_per_line = sparse_matrix.get_set_tab1();
  IntVect& matrix_column_indices = sparse_matrix.get_set_tab2();
  for (int i = 0; i + 1 < non_zero_coeff_per_line.size(); i++) //indice de ligne
    {
      const int index_ini = non_zero_coeff_per_line[i] - FORTRAN_INDEX_INI;
      const int index_end = non_zero_coeff_per_line[i+1] - FORTRAN_INDEX_INI;
      int n = index_end - index_ini;
      DoubleVect matrix_column_indices_local(n);
      DoubleVect matrix_values_local(n);
      for (int j=0; j<n; j++)
        {
          matrix_column_indices_local(j) = matrix_column_indices(j + index_ini);
          matrix_values_local(j) = matrix_values(j + index_ini);
        }
      std::vector<int> indices;
      indices = arg_sort(matrix_column_indices_local);
      for (int j=0; j<n; j++)
        if (j != indices[j])
          matrix_values[j + index_ini] = matrix_values_local[indices[j]];
    }
}

IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::IJK_Finite_Difference_One_Dimensional_Matrix_Assembler()
{
  precision_order_ = 2;
  reduce_side_precision_ = 0;
  core_matrix_type_ = centred;
  equation_type_ = advection_diffusion;
  nb_elem_ = 0;
  non_zero_elem_ = 0;
  stencil_forward_max_ = 0;
  stencil_centred_max_ = 0;
  stencil_backward_max_ = 0;
  computed_stencil_=false;
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

int IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::get_max_operators(const FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE>& fd_operator_conv,
                                                                              const FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE>& fd_operator_diff,
                                                                              int (&stencil_left_offset) [2], int (&stencil_right_offset) [2])
{
  std::vector<int> indices_op;
  const int stencil_width_conv = fd_operator_conv[precision_order_-1][0].size();
  const int stencil_width_diff = fd_operator_diff[precision_order_-1][0].size();
  int min_conv = 0;
  int min_diff = 0;
  int max_conv = 0;
  int max_diff = 0;
  int i;
  for (i=0; i<stencil_width_conv; i++)
    {
      const int index = (int) fd_operator_conv[precision_order_-1][0](i);
      min_conv = std::min(min_conv, (int) fd_operator_conv[precision_order_-1][0](i));
      max_conv = std::max(max_conv, (int) fd_operator_conv[precision_order_-1][0](i));
      if (std::find(indices_op.begin(), indices_op.end(), index) == indices_op.end())
        indices_op.push_back(index);
    }
  for (i=0; i<stencil_width_diff; i++)
    {
      const int index = (int) fd_operator_diff[precision_order_-1][0](i);
      min_diff = std::min(min_diff, (int) fd_operator_diff[precision_order_-1][0](i));
      max_diff = std::max(max_diff, (int) fd_operator_diff[precision_order_-1][0](i));
      if (std::find(indices_op.begin(), indices_op.end(), index) == indices_op.end())
        indices_op.push_back(index);
    }
  const int left_offset = (int) fabs(min_diff-min_conv);
  const int right_offset = (int) fabs(max_diff-max_conv);
  if (min_conv < min_diff)
    {
      stencil_left_offset[0] = 0;
      stencil_left_offset[1] = left_offset;
    }
  else
    {
      stencil_left_offset[0] = left_offset;
      stencil_left_offset[1] = 0;
    }
  if (max_conv < max_diff)
    {
      stencil_right_offset[0] = right_offset;
      stencil_right_offset[1] = 0;
    }
  else
    {
      stencil_right_offset[0] = 0;
      stencil_right_offset[1] = right_offset;
    }
  int non_zero_elem = (int) indices_op.size();
  return non_zero_elem;
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::get_max_stencil(const int& nb_elem,
                                                                             int& non_zero_elem_max,
                                                                             int& stencil_forward_max,
                                                                             int& stencil_centred_max,
                                                                             int& stencil_backward_max)
{
  int stencil_forward_size = 0;
  int stencil_centred_size = 0;
  int stencil_backward_size = 0;
  switch (equation_type_)
    {
    case advection_diffusion:
      stencil_forward_size = get_max_operators(first_order_derivative_forward_, second_order_derivative_forward_, forward_left_offset_, forward_right_offset_);
      stencil_centred_size = get_max_operators(first_order_derivative_centred_, second_order_derivative_centred_, centred_left_offset_, centred_right_offset_);
      stencil_backward_size = get_max_operators(first_order_derivative_backward_, second_order_derivative_backward_, backward_left_offset_, backward_right_offset_);
      break;
    case linear_diffusion:
      stencil_forward_size = second_order_derivative_forward_[precision_order_-1][1].size();
      stencil_centred_size = second_order_derivative_centred_[precision_order_-1][1].size();
      stencil_backward_size = second_order_derivative_backward_[precision_order_-1][1].size();
      break;
    default:
      stencil_forward_size = get_max_operators(first_order_derivative_forward_, second_order_derivative_forward_, forward_left_offset_, forward_right_offset_);
      stencil_centred_size = get_max_operators(first_order_derivative_centred_, second_order_derivative_centred_, centred_left_offset_, centred_right_offset_);
      stencil_backward_size = get_max_operators(first_order_derivative_backward_, second_order_derivative_backward_, backward_left_offset_, backward_right_offset_);
      break;
    }
  stencil_forward_max = stencil_forward_size;
  stencil_centred_max = stencil_centred_size;
  stencil_backward_max = stencil_backward_size;
  const int core_lines = (nb_elem - 2);
  non_zero_elem_max = stencil_forward_max + stencil_backward_max + stencil_centred_max * core_lines;
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::compute_stencil_properties(const int nb_elem)
{
  nb_elem_ = nb_elem;
  if (!computed_stencil_)
    {
      get_max_stencil(nb_elem_, non_zero_elem_, stencil_forward_max_, stencil_centred_max_, stencil_backward_max_);
      computed_stencil_=true;
    }
}

int IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::build(Matrice& matrix, const int& nb_elem, const int& derivative_order)
{
  if (known_pattern_)
    return build_with_known_pattern(matrix, nb_elem, derivative_order);
  else
    return build_with_unknown_pattern(matrix, nb_elem, derivative_order);
}

int IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::build_with_known_pattern(Matrice& matrix, const int& nb_elem, const int& derivative_order)
{
  FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE> * forward_derivative = nullptr;
  FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE> * centred_derivative = nullptr;
  FixedVector<FixedVector<DoubleVect,2>,MAX_ORDER_DERIVATIVE> * backward_derivative = nullptr;
  switch(derivative_order)
    {
    case first:
      forward_derivative = &first_order_derivative_forward_;
      centred_derivative = &first_order_derivative_centred_;
      backward_derivative = &first_order_derivative_backward_;
      break;
    case second:
      forward_derivative = &second_order_derivative_forward_;
      centred_derivative = &second_order_derivative_centred_;
      backward_derivative = &second_order_derivative_backward_;
      break;
    default:
      forward_derivative = &first_order_derivative_forward_;
      centred_derivative = &first_order_derivative_centred_;
      backward_derivative = &first_order_derivative_backward_;
      break;
    }
  compute_stencil_properties(nb_elem);
  // Cast the finite difference matrix
  matrix.typer("Matrice_Morse");
  Matrice_Morse& sparse_matrix  = ref_cast(Matrice_Morse, matrix.valeur());
  sparse_matrix.dimensionner(nb_elem, nb_elem, non_zero_elem_);

  ArrOfDouble& matrix_values = sparse_matrix.get_set_coeff();
  ArrOfInt& non_zero_coeff_per_line = sparse_matrix.get_set_tab1();
  ArrOfInt& matrix_column_indices = sparse_matrix.get_set_tab2();

  Cerr << matrix_values[0] << finl;
  Cerr << non_zero_coeff_per_line[0] << finl;
  Cerr << matrix_column_indices[0] << finl;

  Cerr << (*forward_derivative)[precision_order_-1][0](0) << finl;
  Cerr << (*centred_derivative)	[precision_order_-1][0](0) << finl;
  Cerr << (*backward_derivative)[precision_order_-1][0](0) << finl;

//
//  // Fortran start at one
//  int non_zero_values_counter = 0;
//
//  /*
//   * TODO: Re-write the matrix filling routines
//   */
  const int core_lines = (nb_elem - 2);
  int non_zero_values_counter = 0;

  /*
   * Ini
   */
  {
    const int forward_derivative_size = (*forward_derivative)[precision_order_-1][0].size();
    const int forward_left_offset = forward_left_offset_[derivative_order];
    const int forward_right_offset = forward_right_offset_[derivative_order];
    for (int i=0; i < forward_left_offset; i++)
      {
        matrix_column_indices[non_zero_values_counter] = i + FORTRAN_INDEX_INI;
        non_zero_values_counter++;
      }
    for (int i=0; i < forward_derivative_size; i++)
      {
        const int indices = (int) (*forward_derivative)[precision_order_-1][0](i);
        const double fd_coeff = (*forward_derivative)[precision_order_-1][1](i);
        matrix_column_indices[non_zero_values_counter] = indices + forward_left_offset + FORTRAN_INDEX_INI;
        matrix_values[non_zero_values_counter] = fd_coeff;
        non_zero_values_counter++;
      }
    for (int i=0; i < forward_right_offset; i++)
      {
        matrix_column_indices[non_zero_values_counter] = i + (forward_left_offset + forward_derivative_size) + FORTRAN_INDEX_INI;
        non_zero_values_counter++;
      }
    non_zero_coeff_per_line[1] = non_zero_values_counter + FORTRAN_INDEX_INI;
  }
  /*
   * Core
   */
  {
    const int centred_derivative_size = (*centred_derivative)[precision_order_-1][0].size();
    const int centred_left_offset = centred_left_offset_[derivative_order];
    const int centred_right_offset = centred_right_offset_[derivative_order];
    for (int j=1; j<(core_lines+1); j++)
      {
        for (int i=0; i < centred_left_offset; i++)
          {
            matrix_column_indices[non_zero_values_counter] = i + j + FORTRAN_INDEX_INI;
            non_zero_values_counter++;
          }
        for (int i=0; i < centred_derivative_size; i++)
          {
            const int indices = (int) (*centred_derivative)[precision_order_-1][0](i);
            const double fd_coeff = (*centred_derivative)[precision_order_-1][1](i);
            matrix_column_indices[non_zero_values_counter] = indices + j + centred_left_offset + FORTRAN_INDEX_INI;
            matrix_values[non_zero_values_counter] = fd_coeff;
            non_zero_values_counter++;
          }
        for (int i=0; i < centred_right_offset; i++)
          {
            matrix_column_indices[non_zero_values_counter] = i + (centred_left_offset + centred_derivative_size) + FORTRAN_INDEX_INI;
            non_zero_values_counter++;
          }
        non_zero_coeff_per_line[j + 1] = non_zero_values_counter + FORTRAN_INDEX_INI;
      }
  }

  /*
   * End
   */
  {
    const int backward_derivative_size = (*backward_derivative)[precision_order_-1][0].size();
    int backward_left_offset = backward_left_offset_[derivative_order];
    int backward_right_offset = backward_right_offset_[derivative_order];
    int end_counter = non_zero_values_counter - 1 + backward_derivative_size + backward_left_offset + backward_right_offset;
    const int last_column = (nb_elem - 1);
    for (int i=0; i < backward_right_offset; i++)
      {
        matrix_column_indices[end_counter] = last_column -i + FORTRAN_INDEX_INI;
        end_counter--;
        non_zero_values_counter++;
      }
    for (int i=0; i < backward_derivative_size; i++)
      {
        const int indices = (int) (*backward_derivative)[precision_order_-1][0](i);
        const double fd_coeff = (*backward_derivative)[precision_order_-1][1](i);
        matrix_column_indices[end_counter] = last_column + indices + FORTRAN_INDEX_INI;
        matrix_values[end_counter] = fd_coeff;
        end_counter--;
        non_zero_values_counter++;

      }
    for (int i=0; i < backward_left_offset; i++)
      {
        matrix_column_indices[end_counter] = last_column - i - (backward_right_offset + backward_derivative_size) + FORTRAN_INDEX_INI;
        end_counter--;
        non_zero_values_counter++;
      }
    non_zero_coeff_per_line[last_column + 1] = non_zero_values_counter + FORTRAN_INDEX_INI;
  }

  return non_zero_elem_;
}

/*
 *
 */
int IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::build_with_unknown_pattern(Matrice& matrix, const int& nb_elem, const int& derivative_order)
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
  matrix.typer("Matrice_Morse");
  Matrice_Morse& sparse_matrix  = ref_cast(Matrice_Morse, matrix.valeur());
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
    int end_counter = non_zero_values_counter - 1 + stencil_backward;
    const int last_column = (nb_elem - 1);
    for (int i=0; i < backward_derivative_size; i++)
      {
        const int indices = (int) (*backward_derivative)[precision_order_-1][0](i);
        const double fd_coeff = (*backward_derivative)[precision_order_-1][1](i);
        if (fd_coeff != 0)
          {
            matrix_column_indices[end_counter] = last_column + indices + FORTRAN_INDEX_INI;
            matrix_values[end_counter] = fd_coeff;
            end_counter--;
            non_zero_values_counter++;
          }
      }
    non_zero_coeff_per_line[last_column + 1] = non_zero_values_counter + FORTRAN_INDEX_INI;
  }

  sparse_matrix.compacte(remove);

  return non_zero_elem;
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::sum_matrices_subproblems(Matrice& matrix_A, Matrice& matrix_B)
{
  Matrice_Bloc& block_matrix_A =ref_cast(Matrice_Bloc, matrix_A.valeur());
  Matrice_Bloc& block_matrix_B =ref_cast(Matrice_Bloc, matrix_B.valeur());
  const int dim = block_matrix_A.dim(0);
  for (int i=0; i<dim ; i++)
    {
      Matrice_Morse& sparse_matrix_A  = ref_cast(Matrice_Morse,  block_matrix_A.get_bloc(i,i).valeur());
      Matrice_Morse& sparse_matrix_B  = ref_cast(Matrice_Morse, block_matrix_B.get_bloc(i,i).valeur());
      sparse_matrix_A += sparse_matrix_B;

      // Don't forget to call my own sort_stencil() if matrices have been compacted !
      if (!known_pattern_)
        {
          sort_stencil(sparse_matrix_A);
          sparse_matrix_A.sort_stencil();
        }
      // Cerr << "Convection and Diffusion matrices have been assembled. "<< finl;
    }
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::sum_matrices(Matrice& matrix_A, Matrice& matrix_B)
{
  Matrice_Morse& sparse_matrix_A  = ref_cast(Matrice_Morse, matrix_A.valeur());
  Matrice_Morse& sparse_matrix_B  = ref_cast(Matrice_Morse, matrix_B.valeur());
  sparse_matrix_A += sparse_matrix_B;
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::initialise_matrix_subproblems(Matrice& matrix_subproblems,
                                                                                           Matrice& fd_operator,
                                                                                           const int& nb_subproblems)
{
  matrix_subproblems.typer("Matrice_Bloc");
  Matrice_Bloc& block_matrix_subproblems =ref_cast(Matrice_Bloc, matrix_subproblems.valeur());
  block_matrix_subproblems.dimensionner(nb_subproblems,nb_subproblems);

  Matrice_Morse& fd_operator_sparse = ref_cast(Matrice_Morse, fd_operator.valeur());

  for (int i=0; i<nb_subproblems; i++)
    {
      for (int j=0; j<nb_subproblems; j++)
        if (j != i)
          {
            // Same structure for zeros
            block_matrix_subproblems.get_bloc(i,j).typer("Matrice_Morse");
            Matrice_Morse& sparse_matrix_zeros  = ref_cast(Matrice_Morse, block_matrix_subproblems.get_bloc(i,j).valeur());
            sparse_matrix_zeros.dimensionner(fd_operator_sparse.nb_lignes(), fd_operator_sparse.nb_colonnes(), fd_operator_sparse.nb_coeff());
            sparse_matrix_zeros.compacte(remove);
          }
      block_matrix_subproblems.get_bloc(i,i).typer("Matrice_Morse");
      Matrice_Morse& sparse_matrix  = ref_cast(Matrice_Morse, block_matrix_subproblems.get_bloc(i,i).valeur());
      sparse_matrix = Matrice_Morse(fd_operator_sparse);
    }
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::reinitialise_matrix_subproblem(Matrice& matrix_subproblems,
                                                                                            Matrice& fd_operator,
                                                                                            const int& nb_subproblems)
{
  Matrice_Bloc& block_matrix_subproblems =ref_cast(Matrice_Bloc, matrix_subproblems.valeur());

  Matrice_Morse& sparse_matrix  = ref_cast(Matrice_Morse, block_matrix_subproblems.get_bloc(nb_subproblems,nb_subproblems).valeur());
  Matrice_Morse& fd_operator_morse =ref_cast(Matrice_Morse, fd_operator.valeur());

  sparse_matrix = Matrice_Morse(fd_operator_morse);
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::scale_matrix_by_vector(Matrice& matrix,
                                                                                    const DoubleVect& vector,
                                                                                    const int& boundary_conditions)
{
  /*
   * Don't scale the first and last row (where B.Cs are applied)
   */
  Matrice_Morse& sparse_matrix = ref_cast(Matrice_Morse, matrix.valeur());
  DoubleVect vector_tmp = vector;
  if (boundary_conditions)
    {
      vector_tmp[0] = 1.;
      vector_tmp[vector.size() - 1] = 1.;
    }
  sparse_matrix *= vector_tmp;
  if (!known_pattern_)
    sparse_matrix.compacte(remove);
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

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::impose_boundary_conditions_subproblem(Matrice * matrix,
                                                                                                   DoubleVect * global_rhs,
                                                                                                   DoubleVect& local_rhs,
                                                                                                   const int& ini_boundary_conditions,
                                                                                                   const double& interfacial_value,
                                                                                                   const int& end_boundary_conditions,
                                                                                                   const double& end_value,
                                                                                                   const int& subproblem_index,
                                                                                                   const double& dr_inv)
{
  Matrice_Bloc& block_matrix =ref_cast(Matrice_Bloc, (*matrix).valeur());
  Matrice& sub_matrix = block_matrix.get_bloc(subproblem_index, subproblem_index);
  impose_boundary_conditions(sub_matrix,
                             local_rhs,
                             ini_boundary_conditions,
                             interfacial_value,
                             end_boundary_conditions,
                             end_value,
                             dr_inv);
  /*
   * Fill global RHS
   */
  const int global_rhs_size = (*global_rhs).size();
  const int local_rhs_size = local_rhs.size();
  const int index_start = global_rhs_size - local_rhs_size;
  for (int i=0; i<local_rhs.size(); i++)
    (*global_rhs)[i + index_start] = local_rhs[i];
}


void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::impose_boundary_conditions(Matrice& modified_matrix,
                                                                                        DoubleVect& modified_rhs,
                                                                                        const int& ini_boundary_conditions,
                                                                                        const double& interfacial_value,
                                                                                        const int& end_boundary_conditions,
                                                                                        const double& end_value,
                                                                                        const double& dr_inv)
{
  Matrice_Morse& sparse_matrix = ref_cast(Matrice_Morse, modified_matrix.valeur());
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
          modified_rhs[0] = interfacial_value;
        }
        break;
      case neumann:
        {
          // Refill with first order finite difference coefficients
          const int stencil_forward = non_zero_stencil_values(first_order_derivative_forward_);
          const int non_zero_elem_ini = non_zero_coeff_per_line[1] - FORTRAN_INDEX_INI;
          int counter_fd_coeff = 0;
          for (int i=0; i < non_zero_elem_ini; i++)
            {
              if (i < stencil_forward)
                {
                  const double fd_coeff = first_order_derivative_forward_[precision_order_-1][1](counter_fd_coeff);
                  matrix_values[i] = fd_coeff * dr_inv;
                  counter_fd_coeff++;
                }
              else
                matrix_values[i] = 0.;
            }
          modified_rhs[0] = interfacial_value;
        }
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
          modified_rhs[0] = 0.;
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
          const int non_zero_elem_end = sparse_matrix.nb_vois(nb_lines - 1);
          const int last_index = nb_coeff - 1;
          matrix_values[last_index] = 1.;
          for (int i=last_index - 1; i > ( last_index - non_zero_elem_end); i--)
            matrix_values[i] = 0.;
          modified_rhs[modified_rhs.size() - 1] = end_value;
        }
        break;
      case neumann:
        {
          // Refill with first order finite difference coefficients
          const int non_zero_elem_end = sparse_matrix.nb_vois(nb_lines - 1);
          const int stencil_backward = non_zero_stencil_values(first_order_derivative_backward_);
          int counter_fd_coeff = 0;
          for (int i=0; i < non_zero_elem_end; i++)
            {
              const int index = i + nb_coeff - non_zero_elem_end;
              if (index + stencil_backward >= nb_coeff)
                {
                  const double fd_coeff = first_order_derivative_backward_[precision_order_-1][1](counter_fd_coeff);
                  matrix_values[index] = fd_coeff * dr_inv;
                  counter_fd_coeff++;
                }
              else
                matrix_values[index] = 0.;
            }
          modified_rhs[modified_rhs.size() - 1] = end_value;
        }
        break;
      default:
        {
          const int non_zero_elem_end = sparse_matrix.nb_vois(nb_lines - 1);
          const int last_index = nb_coeff - 1;
          matrix_values[last_index] = 1.;
          for (int i = last_index - 1; i > (last_index - non_zero_elem_end); i--)
            matrix_values[i] = 0.;
          modified_rhs[modified_rhs.size() - 1] = -1.;
        }
        break;
      }
  }
  // TODO: Maybe remove this one
  if (!known_pattern_)
    sparse_matrix.compacte(remove);

  // sparse_matrix.sort_stencil();

  if ((ini_boundary_conditions != neumann && ini_boundary_conditions != flux_jump)
      || (end_boundary_conditions != neumann))
    modify_rhs_for_bc(modified_matrix,
                      modified_rhs,
                      ini_boundary_conditions,
                      end_boundary_conditions);
}

/*
 * Modifcations of the rhs to apply BCs to the resolved 1D field
 */
void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::modify_rhs_for_bc(Matrice& modified_matrix,
                                                                               DoubleVect& modified_rhs,
                                                                               const int& ini_boundary_conditions,
                                                                               const int& end_boundary_conditions)
{
  /*
   * Copy the matrix into modified matrix
   */
  Matrice_Morse& sparse_matrix = ref_cast(Matrice_Morse, modified_matrix.valeur());
  ArrOfDouble& matrix_values = sparse_matrix.get_set_coeff();

  ArrOfInt& non_zero_coeff_per_line = sparse_matrix.get_set_tab1();
  ArrOfInt& matrix_column_indices = sparse_matrix.get_set_tab2();
  const int nb_lines = sparse_matrix.nb_lignes();
  const int nb_column = sparse_matrix.nb_colonnes();

  const DoubleVect rhs = modified_rhs;
  /*
   * Get the modified rhs
   */
  if (ini_boundary_conditions == neumann && ini_boundary_conditions == flux_jump)
    modified_rhs[0] = 0.;
  if (end_boundary_conditions == neumann)
    modified_rhs[modified_rhs.size() - 1] = 0.;

  /*
   * Build B_BCs = A * X_BC
   * and B_modif = B - B_BC
   */
  // modified_rhs = sparse_matrix * rhs;
  sparse_matrix.multvect(rhs, modified_rhs);


  modified_rhs -= rhs;
  modified_rhs *= -1.;
  modified_rhs += rhs;

  /*
   * Remove useless coefficients
   * and build A_modif
   * A_modif X = B_modif
   */
  for (int j=2; j<nb_lines; j++)
    {
      const int non_zero_elem_core_series = non_zero_coeff_per_line[j] - FORTRAN_INDEX_INI;
      const int non_zero_elem_core = non_zero_coeff_per_line[j] - non_zero_coeff_per_line[j-1];
      for (int i = 0; i < non_zero_elem_core; i++)
        {
          const int index_sparse = non_zero_elem_core_series - non_zero_elem_core + i;
          const int column_index = matrix_column_indices[index_sparse] - FORTRAN_INDEX_INI;
          if (column_index == 0)
            matrix_values[index_sparse] = 0.;
          if (column_index == (nb_column - 1))
            matrix_values[index_sparse] = 0.;
        }
    }

  // TODO: Maybe remove this one
  if (!known_pattern_)
    sparse_matrix.compacte(remove);


// sparse_matrix.sort_stencil();
//  for (int i=0; i<modified_rhs.size(); i++)
//    if (fabs(modified_rhs[i]) > 1.e9)
//      Cerr << "Check the modified RHS" << modified_rhs[i] << finl;
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::add_source_terms(DoubleVect * thermal_subproblems_rhs_assembly, const DoubleVect& rhs_assembly)
{

}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::compute_operator(const Matrice * fd_operator, const DoubleVect& solution, DoubleVect& res)
{
  const Matrice_Morse& sparse_matrix = ref_cast(Matrice_Morse, (*fd_operator).valeur());
  sparse_matrix.multvect(solution, res);
  // res = (*fd_operator) * solution;
}

