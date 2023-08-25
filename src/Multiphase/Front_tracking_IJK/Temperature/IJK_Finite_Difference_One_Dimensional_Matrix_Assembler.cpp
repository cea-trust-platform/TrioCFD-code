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
  const int stencil_width = fd_operator[precision_order_][1].size();
  for (int i=0; i<stencil_width; i++)
    if (fd_operator[precision_order_][1][i] != 0)
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
  switch(derivative_order)
    {
    case first:
      stencil_forward = non_zero_stencil_values(first_order_derivative_forward_);
      stencil_centred = non_zero_stencil_values(first_order_derivative_centred_);
      stencil_backward = non_zero_stencil_values(first_order_derivative_backward_);
      break;
    case second:
      stencil_forward = non_zero_stencil_values(second_order_derivative_forward_);
      stencil_centred = non_zero_stencil_values(second_order_derivative_centred_);
      stencil_backward = non_zero_stencil_values(second_order_derivative_backward_);
      break;
    default:
      stencil_forward = non_zero_stencil_values(first_order_derivative_forward_);
      stencil_centred = non_zero_stencil_values(first_order_derivative_centred_);
      stencil_backward = non_zero_stencil_values(first_order_derivative_backward_);
      break;
    }
  non_zero_elem = stencil_forward + stencil_backward + stencil_centred * (nb_elem - 2);
  // Cast the finite difference matrix
  matrix.typer("Matrice_Bloc");
  Matrice_Bloc& block_matrix =ref_cast(Matrice_Bloc, matrix.valeur());
  block_matrix.dimensionner(1,1);
  Matrice_Morse& sparse_matrix  = ref_cast(Matrice_Morse, block_matrix.get_bloc(0,0).valeur());
  sparse_matrix.dimensionner(nb_elem, nb_elem, non_zero_elem);


  /*
   * FIXME : Do three blocks ?
   */
  /*
  block_matrix.dimensionner(3,1);
  Matrice_Morse& ini_sparse_matrix  = ref_cast(Matrice_Morse, block_matrix.get_bloc(0,0).valeur());
  Matrice_Morse& end_sparse_matrix  = ref_cast(Matrice_Morse, block_matrix.get_bloc(2,0).valeur());
  Matrice_Morse_Sym& core_sparse_matrix;
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

/*
 * Modifcations of the rhs to apply BCs to the resolved 1D field
 * Dirichlet - Dirichlet
 * Dirichlet - Neumann
 * Neumann - Neumann
 */
void IJK_Finite_Difference_One_Dimensional_Matrix_Assembler::modify_rhs_for_bc(DoubleVect& rhs)
{

}
