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
// File      : Operateur_IJK_elem_diff.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#include <Operateur_IJK_elem_diff.h>

Implemente_instanciable_sans_constructeur( Operateur_IJK_elem_diff, "Operateur_IJK_elem_diff", DERIV(OpDiffIJKScalarGeneric_double) ) ;

Operateur_IJK_elem_diff::Operateur_IJK_elem_diff()
{
  diffusion_op_words_ = Motcles(6);
  {
    diffusion_op_words_[0] = "standard";
    diffusion_op_words_[1] = "uniform";
    diffusion_op_words_[2] = "anisotropic";
    diffusion_op_words_[3] = "vectorial";
    diffusion_op_words_[4] = "vectorialanisotropic";
    diffusion_op_words_[5] = "structural";
  }
  prefix_ = Nom("OpDiff");
  suffix_ = Nom("IJKScalar_double");
}

Sortie& Operateur_IJK_elem_diff::printOn( Sortie& os ) const
{
  DERIV(OpDiffIJKScalarGeneric_double)::printOn( os );
  return os;
}

Entree& Operateur_IJK_elem_diff::readOn( Entree& is )
{
  Cerr << "Read and Cast Operateur_IJK_elem_diff :" << finl;
  Motcle word;
  is >> word;
  Nom type = "";
  type += prefix_;
  int diffusion_rank = diffusion_op_words_.search(word);
  switch(diffusion_rank)
    {
    case 0 :
      {
        type += "";
        break;
      }
    case 1 :
      {
        type += "Uniform";
        break;
      }
    case 2 :
      {
        type += "Anisotropic";
        break;
      }
    case 3 :
      {
        type += "Vectorial";
        break;
      }
    case 4 :
      {
        type += "VectorialAnisotropic";
        break;
      }
    case 5 :
      {
        type += "StructuralOnly";
        break;
      }
    default :
      {
        Cerr << "ERROR : Diffusion operators that are already implemented are:" << finl;
        Cerr << diffusion_op_words_ << finl;
        abort();
      }
    }
  type += suffix_;
  typer(type);
  is >> valeur();
  return is;
}

