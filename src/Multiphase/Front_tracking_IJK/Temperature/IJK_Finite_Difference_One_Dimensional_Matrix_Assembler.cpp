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

  //Identity
  for(int i=0; i<MAX_ORDER_DERIVATIVE; i++)
    {
      identity_coefficient_[i][0] = DoubleVect(1);
      identity_coefficient_[i][0](0) = 0.;
      identity_coefficient_[i][1] = DoubleVect(1);
      identity_coefficient_[i][1](0) = 1.;
    }
}

Sortie& IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::printOn( Sortie& os ) const
{
  // Objet_U::printOn( os );
  Nom front_space = "    ";
  Nom end_space = " ";
  Nom escape = "\n";
  os << escape;
  os << front_space << "{" << escape;
  if (reduce_side_precision_)
    os << front_space << " reduce_side_precision" << escape;
  os << front_space << " precision_order" << end_space << precision_order_ << escape;
  os << front_space << "}" << escape;
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
                                                                              int (&stencil_left_offset) [3], int (&stencil_right_offset) [3])
{
  std::vector<int> indices_op;
  const int stencil_width_conv = fd_operator_conv[precision_order_-1][0].size();
  const int stencil_width_diff = fd_operator_diff[precision_order_-1][0].size();
  int min_conv = 0;
  int min_diff = 0;
  int max_conv = 0;
  int max_diff = 0;
  int min_identity = 0;
  int max_identity = 0;
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
  const int min_conv_diff_identity = std::min(std::min(min_conv, min_diff), min_identity);
  const int max_conv_diff_identity = std::max(std::max(max_conv, max_diff), max_identity);

  const int left_offset_conv = (int) fabs(min_conv-min_conv_diff_identity);
  const int right_offset_conv = (int) fabs(max_conv-max_conv_diff_identity);
  const int left_offset_diff = (int) fabs(min_diff-min_conv_diff_identity);
  const int right_offset_diff = (int) fabs(max_diff-max_conv_diff_identity);
  const int left_offset_identity = (int) fabs(min_identity-min_conv_diff_identity);
  const int right_offset_identity = (int) fabs(max_identity-max_conv_diff_identity);

  stencil_left_offset[0] = left_offset_conv;
  stencil_left_offset[1] = left_offset_diff;
  stencil_left_offset[2] = left_offset_identity;
  stencil_right_offset[0] = right_offset_conv;
  stencil_right_offset[1] = right_offset_diff;
  stencil_right_offset[2] = right_offset_identity;

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
    case identity:
      forward_derivative = &identity_coefficient_;
      centred_derivative = &identity_coefficient_;
      backward_derivative = &identity_coefficient_;
      break;
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

  int derivative_offset = derivative_order;
  if (derivative_order==identity)
    derivative_offset = 2;
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
    const int forward_left_offset = forward_left_offset_[derivative_offset];
    const int forward_right_offset = forward_right_offset_[derivative_offset];
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
    const int centred_left_offset = centred_left_offset_[derivative_offset];
    const int centred_right_offset = centred_right_offset_[derivative_offset];
    for (int j=1; j<(core_lines+1); j++)
      {
        for (int i=0; i < centred_left_offset; i++)
          {
            // matrix_column_indices[non_zero_values_counter] = i + j + FORTRAN_INDEX_INI; // Add MG 12/02/24
            matrix_column_indices[non_zero_values_counter] = i + j + FORTRAN_INDEX_INI - centred_left_offset; // Add MG 12/02/24
            non_zero_values_counter++;
          }
        for (int i=0; i < centred_derivative_size; i++)
          {
            const int indices = (int) (*centred_derivative)[precision_order_-1][0](i);
            const double fd_coeff = (*centred_derivative)[precision_order_-1][1](i);
            // matrix_column_indices[non_zero_values_counter] = indices + j + centred_left_offset + FORTRAN_INDEX_INI;
            matrix_column_indices[non_zero_values_counter] = indices + j + FORTRAN_INDEX_INI;
            matrix_values[non_zero_values_counter] = fd_coeff;
            non_zero_values_counter++;
          }
        for (int i=0; i < centred_right_offset; i++)
          {
            // matrix_column_indices[non_zero_values_counter] = i + (centred_left_offset + centred_derivative_size) + FORTRAN_INDEX_INI;
            matrix_column_indices[non_zero_values_counter] = i + j + centred_derivative_size + FORTRAN_INDEX_INI;
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
    int backward_left_offset = backward_left_offset_[derivative_offset];
    int backward_right_offset = backward_right_offset_[derivative_offset];
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
    case identity:
      forward_derivative = &identity_coefficient_;
      centred_derivative = &identity_coefficient_;
      backward_derivative = &identity_coefficient_;
      break;
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

  sparse_matrix.compacte(REMOVE);

  return non_zero_elem;
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::sum_any_matrices_subproblems(Matrice& matrix_A, Matrice& matrix_B,
                                                                                          const int& use_sparse_matrix,
                                                                                          const int& debug)
{
  if (use_sparse_matrix)
    sum_sparse_matrices_subproblems(matrix_A, matrix_B, debug);
  else
    sum_matrices_subproblems(matrix_A, matrix_B);

}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::sum_sparse_matrices_subproblems(Matrice& matrix_A, Matrice& matrix_B, const int& debug)
{
  Matrice_Morse& sparse_matrix_A  = ref_cast(Matrice_Morse, matrix_A.valeur());
  Matrice_Morse& sparse_matrix_B  = ref_cast(Matrice_Morse, matrix_B.valeur());
  sparse_matrix_A += sparse_matrix_B;

  if (debug)
    {
      const IntVect& tab1 = sparse_matrix_A.get_tab1();
      Cerr << "tab1[tab1.size_array() - 1]" << tab1[tab1.size_array() - 1] << finl;
    }

  // Don't forget to call my own sort_stencil() if matrices have been compacted !
  if (!known_pattern_)
    {
      sort_stencil(sparse_matrix_A);
      sparse_matrix_A.sort_stencil();
    }
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

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::pre_initialise_matrix_subproblems(Matrice& matrix_subproblems,
                                                                                               Matrice& fd_operator,
                                                                                               const int& max_subproblems_predicted)
{
  matrix_subproblems.typer("Matrice_Bloc");
  Matrice_Bloc& block_matrix_subproblems =ref_cast(Matrice_Bloc, matrix_subproblems.valeur());
  block_matrix_subproblems.dimensionner(max_subproblems_predicted,max_subproblems_predicted);

  Matrice_Morse& fd_operator_sparse = ref_cast(Matrice_Morse, fd_operator.valeur());

  for (int i=0; i<max_subproblems_predicted; i++)
    {
      for (int j=0; j<max_subproblems_predicted; j++)
        if (j != i)
          {
            // Same structure for zeros
            block_matrix_subproblems.get_bloc(i,j).typer("Matrice_Morse");
            Matrice_Morse& sparse_matrix_zeros  = ref_cast(Matrice_Morse, block_matrix_subproblems.get_bloc(i,j).valeur());
            sparse_matrix_zeros.dimensionner(fd_operator_sparse.nb_lignes(), fd_operator_sparse.nb_colonnes(), fd_operator_sparse.nb_coeff());
            sparse_matrix_zeros.compacte(REMOVE);
          }
      block_matrix_subproblems.get_bloc(i,i).typer("Matrice_Morse");
      Matrice_Morse& sparse_matrix  = ref_cast(Matrice_Morse, block_matrix_subproblems.get_bloc(i,i).valeur());
      sparse_matrix.dimensionner(fd_operator_sparse.nb_lignes(), fd_operator_sparse.nb_colonnes(), fd_operator_sparse.nb_coeff());
    }
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::pre_initialise_sparse_matrix_subproblems(Matrice& matrix_subproblems,
                                                                                                      Matrice& fd_operator,
                                                                                                      const int& max_subproblems)
{
  Matrice_Morse& fd_operator_sparse = ref_cast(Matrice_Morse, fd_operator.valeur());
  const int nb_rows = (fd_operator_sparse.nb_lignes()) * max_subproblems;
  const int nb_column = (fd_operator_sparse.nb_colonnes()) * max_subproblems;
  const int nb_coeff = (fd_operator_sparse.nb_coeff()) * max_subproblems;

  matrix_subproblems.typer("Matrice_Morse");
  Matrice_Morse& sparse_matrix_subproblems  = ref_cast(Matrice_Morse, matrix_subproblems.valeur());
  sparse_matrix_subproblems.dimensionner(nb_rows, nb_column, nb_coeff);
}



void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::complete_empty_matrices_initialisation(Matrice& matrix_subproblems,
                                                                                                    Matrice& fd_operator,
                                                                                                    const int& empty_problem_start_index,
                                                                                                    const int& empty_problem_end_index)
{
  Matrice_Bloc& block_matrix_subproblems =ref_cast(Matrice_Bloc, matrix_subproblems.valeur());
  Matrice_Morse& fd_operator_morse =ref_cast(Matrice_Morse, fd_operator.valeur());
  const DoubleVect& coeff = fd_operator_morse.get_coeff();
  const IntVect& tab1 = fd_operator_morse.get_tab1();
  const IntVect& tab2 = fd_operator_morse.get_tab2();
  for (int i=empty_problem_start_index; i<empty_problem_end_index; i++)
    {
      Matrice_Morse& sparse_matrix  = ref_cast(Matrice_Morse, block_matrix_subproblems.get_bloc(i,i).valeur());
      sparse_matrix.get_set_coeff() = coeff;
      sparse_matrix.get_set_tab1() = tab1;
      sparse_matrix.get_set_tab2() = tab2;
    }
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::initialise_sparse_matrix_subproblems(Matrice& matrix_subproblems,
                                                                                                  Matrice& fd_operator,
                                                                                                  const int& nb_subproblems,
                                                                                                  const int& first_time_step_varying_probes,
                                                                                                  FixedVector<ArrOfInt, 6>& first_indices_sparse_matrix,
                                                                                                  const int& first_initialisation,
                                                                                                  int& initialise_sparse_indices)
{

  if (first_initialisation)
    {
      matrix_subproblems.typer("Matrice_Morse");
    }
  if (initialise_sparse_indices)
    {
      for (int l=0; l<6; l++)
        first_indices_sparse_matrix[l].reset();
    }

  Matrice_Morse& sparse_matrix_subproblems =ref_cast(Matrice_Morse, matrix_subproblems.valeur());
  Matrice_Morse& fd_operator_sparse = ref_cast(Matrice_Morse, fd_operator.valeur());

  const IntVect& tab1_operator = fd_operator_sparse.get_tab1();
  const IntVect& tab2_operator = fd_operator_sparse.get_tab2();
  const DoubleVect& coeff_operator = fd_operator_sparse.get_coeff();

  const int nb_rows_operator = fd_operator_sparse.nb_lignes();
  const int nb_columns_operator = fd_operator_sparse.nb_colonnes();

  const int nb_rows = nb_rows_operator * nb_subproblems;
  const int nb_columns = nb_columns_operator * nb_subproblems;

  const int nb_val_non_zero_per_operator = fd_operator_sparse.nb_coeff();
  const int nb_val_non_zero_sparse = nb_subproblems * nb_val_non_zero_per_operator;

  IntVect& tab1 = sparse_matrix_subproblems.get_set_tab1();
  IntVect& tab2 = sparse_matrix_subproblems.get_set_tab2();
  DoubleVect& coeff = sparse_matrix_subproblems.get_set_coeff();

  int j;
  /*
   * dimensionner resize each IntVect tab1 and tab2 and DoubleVect coeff
   */
  sparse_matrix_subproblems.dimensionner(nb_rows, nb_columns, nb_val_non_zero_sparse);
  for (int i=0; i<nb_subproblems; i++)
    {
      const int first_index = i * nb_val_non_zero_per_operator;
      const int column_index = i * nb_columns_operator;
      const int row_index = i * nb_rows_operator;
      if (initialise_sparse_indices)
        {
          first_indices_sparse_matrix[0].append_array(first_index);
          first_indices_sparse_matrix[1].append_array(nb_val_non_zero_per_operator);
          first_indices_sparse_matrix[2].append_array(column_index);
          first_indices_sparse_matrix[3].append_array(nb_columns_operator);
          first_indices_sparse_matrix[4].append_array(row_index);
          first_indices_sparse_matrix[5].append_array(nb_rows_operator);
        }
      for (j=0; j<nb_val_non_zero_per_operator; j++)
        {
          tab2(j + first_index) = tab2_operator(j) + column_index;
          coeff(j + first_index) = coeff_operator(j);
        }
      for (j=1; j<=nb_rows_operator; j++)
        tab1(j + row_index) = tab1_operator(j) + first_index;
    }
  initialise_sparse_indices = 0;
//if (!first_time_step_varying_probes)
//	{
//		Matrice_Morse& sparse_matrix  = ref_cast(Matrice_Morse, block_matrix_subproblems.get_bloc(i,i).valeur());
//		sparse_matrix = Matrice_Morse(fd_operator_sparse);
//	}
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::initialise_matrix_subproblems(Matrice& matrix_subproblems,
                                                                                           Matrice& fd_operator,
                                                                                           const int& nb_subproblems,
                                                                                           const int& first_time_step_varying_probes)
{
  /*
   * Work only for constant number of points on probes
   * TODO: Adapt for varying numbers of points
   */
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
            sparse_matrix_zeros.compacte(REMOVE);
          }
      block_matrix_subproblems.get_bloc(i,i).typer("Matrice_Morse");
      if (!first_time_step_varying_probes)
        {
          Matrice_Morse& sparse_matrix  = ref_cast(Matrice_Morse, block_matrix_subproblems.get_bloc(i,i).valeur());
          sparse_matrix = Matrice_Morse(fd_operator_sparse);
        }
    }
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::reinitialise_any_matrix_subproblem(Matrice * matrix_subproblems,
                                                                                                const Matrice * fd_operator,
                                                                                                const int& nb_subproblems,
                                                                                                const int& use_sparse_matrix,
                                                                                                FixedVector<ArrOfInt,6> * first_indices_sparse_matrix,
                                                                                                const int& first_initialisation,
                                                                                                const int& keep_global_probes_discretisation)
{
  if (use_sparse_matrix)
    reinitialise_sparse_matrix_subproblem(matrix_subproblems, fd_operator, nb_subproblems, first_indices_sparse_matrix, first_initialisation, keep_global_probes_discretisation);
  else
    reinitialise_matrix_subproblem(matrix_subproblems, fd_operator, nb_subproblems, keep_global_probes_discretisation);
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::reinitialise_sparse_matrix_subproblem(Matrice * matrix_subproblems,
                                                                                                   const Matrice * fd_operator,
                                                                                                   const int& nb_subproblems,
                                                                                                   FixedVector<ArrOfInt,6> * first_indices_sparse_matrix,
                                                                                                   const int& first_initialisation,
                                                                                                   const int& keep_global_probes_discretisation)
{
  Matrice_Morse& sparse_matrix_subproblems =ref_cast(Matrice_Morse, (*matrix_subproblems).valeur());
  IntVect& tab1 = sparse_matrix_subproblems.get_set_tab1();
  IntVect& tab2 = sparse_matrix_subproblems.get_set_tab2();
  DoubleVect& coeff = sparse_matrix_subproblems.get_set_coeff();

  const Matrice_Morse& fd_operator_morse =ref_cast(Matrice_Morse, (*fd_operator).valeur());
  const IntVect& tab1_operator = fd_operator_morse.get_tab1();
  const IntVect& tab2_operator = fd_operator_morse.get_tab2();
  const DoubleVect& coeff_operator = fd_operator_morse.get_coeff();

  int nb_rows_operator = fd_operator_morse.nb_lignes();
  int nb_columns_operator = fd_operator_morse.nb_colonnes();
  int nb_coeff_operator = fd_operator_morse.nb_coeff();

  int nb_rows = nb_rows_operator;
  int nb_columns = nb_columns_operator;
  int nb_coeff = nb_coeff_operator;

  int nb_rows_ini = sparse_matrix_subproblems.nb_lignes();
  int nb_columns_ini = sparse_matrix_subproblems.nb_colonnes();
  int nb_coeff_ini = sparse_matrix_subproblems.nb_coeff();

  /*
   * Re-write the entire matrix (can be usefull if probe discretisation changes)
   */
  if (!nb_subproblems)
    {
      if (!keep_global_probes_discretisation)
        {
          nb_rows_ini = 0;
          nb_columns_ini = 0;
          nb_coeff_ini = 0;
        }
      else
        {
          nb_rows_ini = (*first_indices_sparse_matrix)[4][nb_subproblems];
          nb_columns_ini = (*first_indices_sparse_matrix)[2][nb_subproblems];
          nb_coeff_ini = (*first_indices_sparse_matrix)[0][nb_subproblems];
        }
    }
  else
    {
      if (!keep_global_probes_discretisation)
        {
          nb_rows += nb_rows_ini;
          nb_columns += nb_columns_ini;
          nb_coeff += nb_coeff_ini;
        }
      else
        {
          nb_rows_ini = (*first_indices_sparse_matrix)[4][nb_subproblems];
          nb_columns_ini = (*first_indices_sparse_matrix)[2][nb_subproblems];
          nb_coeff_ini = (*first_indices_sparse_matrix)[0][nb_subproblems];
          nb_rows = nb_rows_ini;
          nb_columns = nb_columns_ini;
          nb_coeff = nb_coeff_ini;
        }
    }

  if (!keep_global_probes_discretisation)
    sparse_matrix_subproblems.dimensionner(nb_rows, nb_columns, nb_coeff);
  assert(nb_rows>=nb_rows_ini && nb_columns>=nb_columns_ini && nb_rows>=nb_rows_ini);

  if (first_initialisation && !keep_global_probes_discretisation)
    {
      (*first_indices_sparse_matrix)[0][nb_subproblems] = nb_coeff_ini;
      (*first_indices_sparse_matrix)[1][nb_subproblems] = nb_coeff_operator;
      (*first_indices_sparse_matrix)[2][nb_subproblems] = nb_columns_ini;
      (*first_indices_sparse_matrix)[3][nb_subproblems] = nb_columns_operator;
      (*first_indices_sparse_matrix)[4][nb_subproblems] = nb_rows_ini;
      (*first_indices_sparse_matrix)[5][nb_subproblems] = nb_rows_operator;
      // (*first_indices_sparse_matrix)[4][nb_subproblems] = nb_coeff_ini;
    }

  const int coeff_offset = nb_coeff_ini;
  const int column_offset = nb_columns_ini;
  const int row_offset = nb_coeff_ini;

  int i;
  if (!keep_global_probes_discretisation)
    {
      for (i=0; i<nb_coeff_operator; i++)
        {
          tab2(i + coeff_offset) = tab2_operator(i) + column_offset;
          coeff(i + coeff_offset) = coeff_operator(i);
        }
      for (i=1; i<=nb_rows_operator; i++)
        tab1(i + row_offset) = tab1_operator(i) + coeff_offset;
    }
  else
    for (i=0; i<nb_coeff_operator; i++)
      coeff(i + coeff_offset) = coeff_operator(i);
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::reinitialise_matrix_subproblem(Matrice * matrix_subproblems,
                                                                                            const Matrice * fd_operator,
                                                                                            const int& nb_subproblems,
                                                                                            const int& keep_global_probes_discretisation)
{
  Matrice_Bloc& block_matrix_subproblems =ref_cast(Matrice_Bloc, (*matrix_subproblems).valeur());
  Matrice_Morse& sparse_matrix  = ref_cast(Matrice_Morse, block_matrix_subproblems.get_bloc(nb_subproblems,nb_subproblems).valeur());
  const Matrice_Morse& fd_operator_morse =ref_cast(Matrice_Morse, (*fd_operator).valeur());
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
    sparse_matrix.compacte(REMOVE);
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::scale_matrix_subproblem_by_vector(Matrice * matrix,
                                                                                               const DoubleVect& vector,
                                                                                               const int& subproblem_index,
                                                                                               const int& boundary_conditions,
                                                                                               const int& use_sparse_matrix,
                                                                                               const FixedVector<ArrOfInt,6> * first_indices_sparse_matrix)
{
  if (use_sparse_matrix)
    {
      Matrice local_sub_matrix;
      make_operation_on_sub_matrix_sparse(local_sub_matrix,
                                          matrix,
                                          subproblem_index,
                                          first_indices_sparse_matrix);
      scale_matrix_by_vector(local_sub_matrix, vector, boundary_conditions);
      recombined_local_sub_matrix_with_matrix(local_sub_matrix,
                                              matrix,
                                              subproblem_index,
                                              first_indices_sparse_matrix);
    }
  else
    {
      Matrice_Bloc& block_matrix = ref_cast(Matrice_Bloc, (*matrix).valeur());
      Matrice& sub_matrix = block_matrix.get_bloc(subproblem_index, subproblem_index);
      scale_matrix_by_vector(sub_matrix, vector, boundary_conditions);
    }
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::make_operation_on_sub_matrix_sparse(Matrice& local_sub_matrix,
                                                                                                 Matrice * matrix,
                                                                                                 const int& subproblem_index,
                                                                                                 const FixedVector<ArrOfInt,6> * first_indices_sparse_matrix)
{
  local_sub_matrix.typer("Matrice_Morse");
  Matrice_Morse& local_sparse_sub_matrix = ref_cast(Matrice_Morse, local_sub_matrix.valeur());

  const int nb_coeff = (*first_indices_sparse_matrix)[1][subproblem_index];
  const int nb_column = (*first_indices_sparse_matrix)[3][subproblem_index];
  const int nb_rows = (*first_indices_sparse_matrix)[5][subproblem_index];
  local_sparse_sub_matrix.dimensionner(nb_rows, nb_column, nb_coeff);

  IntVect& tab1_local = local_sparse_sub_matrix.get_set_tab1();
  IntVect& tab2_local = local_sparse_sub_matrix.get_set_tab2();
  DoubleVect& coeff_local = local_sparse_sub_matrix.get_set_coeff();

  Matrice_Morse& sparse_matrix = ref_cast(Matrice_Morse, (*matrix).valeur());
  const IntVect& tab1 = sparse_matrix.get_tab1();
  const IntVect& tab2 = sparse_matrix.get_tab2();
  const DoubleVect& coeff = sparse_matrix.get_coeff();

  const int coeff_offset = (*first_indices_sparse_matrix)[0][subproblem_index];
  const int column_offset = (*first_indices_sparse_matrix)[2][subproblem_index];
  const int row_offset = (*first_indices_sparse_matrix)[4][subproblem_index];
  int i;
  for (i=0; i<nb_coeff; i++)
    {
      tab2_local(i) = tab2(i + coeff_offset) - column_offset;
      coeff_local(i) = coeff(i + coeff_offset);
    }
  tab1_local(0) = FORTRAN_INDEX_INI;
  for (i=1; i<=nb_rows; i++)
    tab1_local(i) = tab1(i + row_offset) - coeff_offset;
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::recombined_local_sub_matrix_with_matrix(Matrice& local_sub_matrix,
                                                                                                     Matrice * matrix,
                                                                                                     const int& subproblem_index,
                                                                                                     const FixedVector<ArrOfInt,6> * first_indices_sparse_matrix)
{
  Matrice_Morse& local_sparse_sub_matrix = ref_cast(Matrice_Morse, local_sub_matrix.valeur());
  const IntVect& tab1_local = local_sparse_sub_matrix.get_tab1();
  const IntVect& tab2_local = local_sparse_sub_matrix.get_tab2();
  const DoubleVect& coeff_local = local_sparse_sub_matrix.get_coeff();

  Matrice_Morse& sparse_matrix = ref_cast(Matrice_Morse, (*matrix).valeur());
  IntVect& tab1 = sparse_matrix.get_set_tab1();
  IntVect& tab2 = sparse_matrix.get_set_tab2();
  DoubleVect& coeff = sparse_matrix.get_set_coeff();

  const int nb_coeff = (*first_indices_sparse_matrix)[1][subproblem_index];
  const int nb_rows = (*first_indices_sparse_matrix)[5][subproblem_index];

  const int coeff_offset = (*first_indices_sparse_matrix)[0][subproblem_index];
  const int column_offset = (*first_indices_sparse_matrix)[2][subproblem_index];
  const int row_offset = (*first_indices_sparse_matrix)[4][subproblem_index];

  int i;
  for (i=0; i<nb_coeff; i++)
    {
      tab2(i + coeff_offset) = tab2_local(i) + column_offset;
      coeff(i + coeff_offset) = coeff_local(i);
    }
  for (i=1; i<=nb_rows; i++)
    tab1(i + row_offset) = tab1_local(i) + coeff_offset;
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::impose_boundary_conditions_subproblem(Matrice * matrix,
                                                                                                   DoubleVect * global_rhs,
                                                                                                   DoubleVect& local_rhs,
                                                                                                   const int& ini_boundary_conditions,
                                                                                                   const double& interfacial_value,
                                                                                                   const int& end_boundary_conditions,
                                                                                                   const double& end_value,
                                                                                                   const int& subproblem_index,
                                                                                                   const double& dr_inv,
                                                                                                   const int& first_time_step_temporal,
                                                                                                   const int& first_time_step_explicit,
                                                                                                   const DoubleVect& temperature_ini_temporal_schemes,
                                                                                                   const int& start_index,
                                                                                                   const FixedVector<ArrOfInt, 6> * first_indices_sparse_matrix,
                                                                                                   const int& use_sparse_matrix)
{
  if (use_sparse_matrix)
    {
      Matrice local_sub_matrix;
      make_operation_on_sub_matrix_sparse(local_sub_matrix,
                                          matrix,
                                          subproblem_index,
                                          first_indices_sparse_matrix);
      impose_boundary_conditions(local_sub_matrix,
                                 local_rhs,
                                 ini_boundary_conditions,
                                 interfacial_value,
                                 end_boundary_conditions,
                                 end_value,
                                 dr_inv,
                                 first_time_step_temporal,
                                 first_time_step_explicit,
                                 temperature_ini_temporal_schemes);
      recombined_local_sub_matrix_with_matrix(local_sub_matrix,
                                              matrix,
                                              subproblem_index,
                                              first_indices_sparse_matrix);
    }
  else
    {
      Matrice_Bloc& block_matrix =ref_cast(Matrice_Bloc, (*matrix).valeur());
      Matrice& sub_matrix = block_matrix.get_bloc(subproblem_index, subproblem_index);
      impose_boundary_conditions(sub_matrix,
                                 local_rhs,
                                 ini_boundary_conditions,
                                 interfacial_value,
                                 end_boundary_conditions,
                                 end_value,
                                 dr_inv,
                                 first_time_step_temporal,
                                 first_time_step_explicit,
                                 temperature_ini_temporal_schemes);
    }
  /*
   * Fill global RHS
   */
  for (int i=0; i<local_rhs.size(); i++)
    (*global_rhs)[i + start_index] = local_rhs[i];
}


void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::impose_boundary_conditions(Matrice& modified_matrix,
                                                                                        DoubleVect& modified_rhs,
                                                                                        const int& ini_boundary_conditions,
                                                                                        const double& interfacial_value,
                                                                                        const int& end_boundary_conditions,
                                                                                        const double& end_value,
                                                                                        const double& dr_inv,
                                                                                        const int& first_time_step_temporal,
                                                                                        const int& first_time_step_explicit,
                                                                                        const DoubleVect& temperature_ini_temporal_schemes)
{
  Matrice_Morse& sparse_matrix = ref_cast(Matrice_Morse, modified_matrix.valeur());
  ArrOfDouble& matrix_values = sparse_matrix.get_set_coeff();
  ArrOfInt& non_zero_coeff_per_line = sparse_matrix.get_set_tab1();
  const int nb_lines = sparse_matrix.nb_lignes();
  const int nb_coeff = sparse_matrix.nb_coeff();
  double interfacial_value_rhs;
  int ini_boundary_conditions_static_temporal = ini_boundary_conditions;
  if (first_time_step_temporal)
    {
      interfacial_value_rhs = temperature_ini_temporal_schemes[0];
      ini_boundary_conditions_static_temporal = dirichlet;
    }
  else
    interfacial_value_rhs = interfacial_value;
  {
    /*
     * Consider a constant interfacial temperature value at saturation
     * A jump condition could be implemented later by equating
     * the fluxes crossing the interface
     */
    switch(ini_boundary_conditions_static_temporal)
      {
      case dirichlet:
        {
          const int non_zero_elem_ini = non_zero_coeff_per_line[1] - FORTRAN_INDEX_INI;
          matrix_values[0] = 1.;
          for (int i=1; i<non_zero_elem_ini; i++)
            matrix_values[i] = 0.;
          modified_rhs[0] = interfacial_value_rhs;
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
      case implicit:
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
  if (!first_time_step_temporal)
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
            const int stencil_backward = first_order_derivative_backward_[precision_order_-1][1].size_array();
            // const int stencil_backward = non_zero_stencil_values(first_order_derivative_backward_);
            int counter_fd_coeff = stencil_backward - 1;
            for (int i=0; i < non_zero_elem_end; i++)
              {
                const int index = i + nb_coeff - non_zero_elem_end;
                if (index + stencil_backward >= nb_coeff)
                  {
                    const double fd_coeff = first_order_derivative_backward_[precision_order_-1][1](counter_fd_coeff);
                    matrix_values[index] = fd_coeff * dr_inv;
                    counter_fd_coeff--;
                  }
                else
                  matrix_values[index] = 0.;
              }
            modified_rhs[modified_rhs.size() - 1] = end_value;
          }
          break;
        case implicit:
          {
            // Refill with first order finite difference coefficients
            const int non_zero_elem_end = sparse_matrix.nb_vois(nb_lines - 1);
            const int stencil_backward = first_order_derivative_backward_[precision_order_-1][1].size_array();
            // const int stencil_backward = non_zero_stencil_values(first_order_derivative_backward_);
            int counter_fd_coeff = stencil_backward - 1;
            for (int i=0; i < non_zero_elem_end; i++)
              {
                const int index = i + nb_coeff - non_zero_elem_end;
                if (index + stencil_backward >= nb_coeff)
                  {
                    const double fd_coeff = first_order_derivative_backward_[precision_order_-1][1](counter_fd_coeff);
                    matrix_values[index] = fd_coeff * dr_inv;
                    counter_fd_coeff--;
                  }
                else
                  matrix_values[index] = 0.;
              }
            const int stencil_centred = first_order_derivative_centred_[precision_order_-1][1].size_array();
            // const int stencil_centred = non_zero_stencil_values(first_order_derivative_centred_);
            counter_fd_coeff = 0;
            for (int i=0; i < non_zero_elem_end; i++)
              {
                const int index = i + nb_coeff - non_zero_elem_end;
                if (index + stencil_centred >= nb_coeff)
                  {
                    const double fd_coeff = first_order_derivative_centred_[precision_order_-1][1](counter_fd_coeff);
                    matrix_values[index] += - (fd_coeff * dr_inv);
                    counter_fd_coeff++;
                  }
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
  else
    Cerr << "First iteration is done with Euler implicit or explicit" << finl;
  // TODO: Maybe remove this one
  if (!known_pattern_)
    sparse_matrix.compacte(REMOVE);

  // sparse_matrix.sort_stencil();
  if (!(first_time_step_temporal && first_time_step_explicit))
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

  const DoubleVect rhs_ini = modified_rhs;
  DoubleVect rhs = modified_rhs;

  const int ini_boundary_bool = ((ini_boundary_conditions==default_bc) || (ini_boundary_conditions==dirichlet));
  const int end_boundary_bool = ((end_boundary_conditions==default_bc) || (end_boundary_conditions==dirichlet));
  /*
   * Get the modified rhs
   */
  if (ini_boundary_conditions == neumann && ini_boundary_conditions == flux_jump)
    rhs[0] = 0.;
  if (end_boundary_conditions == neumann)
    rhs[rhs.size() - 1] = 0.;

  /*
   * Build B_BCs = A * X_BC
   * and B_modif = B - B_BC
   */
  // modified_rhs = sparse_matrix * rhs;
  sparse_matrix.multvect(rhs, modified_rhs);


  modified_rhs -= rhs;
  modified_rhs *= -1.;
  modified_rhs += rhs;

  if (ini_boundary_conditions == neumann && ini_boundary_conditions == flux_jump)
    modified_rhs[0] = rhs_ini[0];
  if (end_boundary_conditions == neumann)
    modified_rhs[modified_rhs.size() - 1] = rhs_ini[rhs_ini.size() - 1];

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
          if (column_index == 0 && ini_boundary_bool)
            matrix_values[index_sparse] = 0.;
          if (column_index == (nb_column - 1) && end_boundary_bool)
            matrix_values[index_sparse] = 0.;
        }
    }

  // TODO: Maybe remove this one
  if (!known_pattern_)
    sparse_matrix.compacte(REMOVE);


// sparse_matrix.sort_stencil();
//  for (int i=0; i<modified_rhs.size(); i++)
//    if (fabs(modified_rhs[i]) > 1.e9)
//      Cerr << "Check the modified RHS" << modified_rhs[i] << finl;
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::add_source_terms(DoubleVect * thermal_subproblems_rhs_assembly,
                                                                              DoubleVect& rhs_assembly,
                                                                              const DoubleVect& source_terms,
                                                                              const int& index_start,
                                                                              const int& boundary_condition_interface,
                                                                              const int& boundary_condition_end)
{
  /*
   * Fill global RHS
   */
  const int local_rhs_size = source_terms.size();
  int index_ini = 0;
  int index_end = local_rhs_size;
  switch(boundary_condition_interface)
    {
    case dirichlet:
      index_ini++;
      break;
    case neumann:
      index_ini++;
      break;
    case flux_jump:
      break;
    case implicit:
      break;
    default:
      index_ini++;
      break;
    }
  switch(boundary_condition_end)
    {
    case dirichlet:
      index_end--;
      break;
    case neumann:
      index_end--;
      break;
    case flux_jump:
      break;
    case implicit:
      // index_end--;
      break;
    default:
      index_end--;
      break;
    }
  for (int i=index_ini; i<index_end; i++)
    {
      (*thermal_subproblems_rhs_assembly)[i + index_start] += source_terms[i];
      rhs_assembly[i] += source_terms[i];
    }
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::compute_operator(const Matrice * fd_operator, const DoubleVect& solution, DoubleVect& res)
{
  const Matrice_Morse& sparse_matrix = ref_cast(Matrice_Morse, (*fd_operator).valeur());
  sparse_matrix.multvect(solution, res);
  // res = (*fd_operator) * solution;
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::apply_euler_time_step(Matrice * convection_matrix,
                                                                                   Matrice * diffusion_matrix,
                                                                                   const int& subproblem_index,
                                                                                   const double& local_time_step,
                                                                                   const double& alpha)
{
  Matrice_Bloc& convection_block_matrix =ref_cast(Matrice_Bloc, (*convection_matrix).valeur());
  Matrice& convection_sub_matrix = convection_block_matrix.get_bloc(subproblem_index, subproblem_index);
  convection_sub_matrix *= (local_time_step * alpha);

  Matrice_Bloc& diffusion_block_matrix =ref_cast(Matrice_Bloc, (*diffusion_matrix).valeur());
  Matrice& diffusion_sub_matrix = diffusion_block_matrix.get_bloc(subproblem_index, subproblem_index);
  diffusion_sub_matrix *= (local_time_step * alpha);
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::correct_sign_temporal_schemes_subproblems(Matrice * convection_matrix,
                                                                                                       Matrice * diffusion_matrix,
                                                                                                       const int& subproblem_index,
                                                                                                       const double& local_time_step,
                                                                                                       const double& alpha)
{
  Matrice_Bloc& convection_block_matrix =ref_cast(Matrice_Bloc, (*convection_matrix).valeur());
  Matrice& convection_sub_matrix = convection_block_matrix.get_bloc(subproblem_index, subproblem_index);
  convection_sub_matrix *= (- local_time_step * alpha);

  Matrice_Bloc& diffusion_block_matrix =ref_cast(Matrice_Bloc, (*diffusion_matrix).valeur());
  Matrice& diffusion_sub_matrix = diffusion_block_matrix.get_bloc(subproblem_index, subproblem_index);
  diffusion_sub_matrix *= (- local_time_step * alpha);
}

void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::reduce_solver_matrix(const Matrice_Morse& sparse_matrix_solver,
                                                                                  Matrice_Morse& sparse_matrix_solver_reduced,
                                                                                  const int& nb_points,
                                                                                  const int& pre_initialise_thermal_subproblems_list)
{
  if (pre_initialise_thermal_subproblems_list)
    {
      const IntVect& tab1 = sparse_matrix_solver.get_tab1();
      const IntVect& tab2 = sparse_matrix_solver.get_tab2();
      const DoubleVect& coeff = sparse_matrix_solver.get_coeff();
      sparse_matrix_solver_reduced.dimensionner(nb_points, nb_points, tab1[nb_points] - FORTRAN_INDEX_INI);
      IntVect& tab1_reduced = sparse_matrix_solver_reduced.get_set_tab1();
      IntVect& tab2_reduced = sparse_matrix_solver_reduced.get_set_tab2();
      DoubleVect& coeff_reduced = sparse_matrix_solver_reduced.get_set_coeff();
      const int nb_coeff = coeff_reduced.size();
      const int nb_rows = tab1_reduced.size();

//      tab1_reduced = tab1;
//      tab2_reduced = tab2;
//      coeff_reduced = coeff;
//      tab1_reduced.resize(nb_rows, RESIZE_OPTIONS::COPY_NOINIT);
//      tab2_reduced.resize(nb_coeff, RESIZE_OPTIONS::COPY_NOINIT);
//      coeff_reduced.resize(nb_coeff, RESIZE_OPTIONS::COPY_NOINIT);

      int i;
      for (i=0; i<nb_rows; i++)
        tab1_reduced(i) = tab1(i);
      for (i=0; i<nb_coeff; i++)
        {
          tab2_reduced(i) = tab2(i);
          coeff_reduced(i) = coeff(i);
        }
    }
  else
    sparse_matrix_solver_reduced = sparse_matrix_solver;
}


