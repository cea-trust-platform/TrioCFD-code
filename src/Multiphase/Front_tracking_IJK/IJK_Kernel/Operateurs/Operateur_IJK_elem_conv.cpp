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
// File      : Operateur_IJK_elem_conv.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#include <Operateur_IJK_elem_conv.h>

Implemente_instanciable_sans_constructeur( Operateur_IJK_elem_conv, "Operateur_IJK_elem_conv", DERIV(OpConvIJKElemCommon_double) ) ;

Operateur_IJK_elem_conv::Operateur_IJK_elem_conv()
{
  convection_op_words_ = Motcles(4);
  {
    convection_op_words_[0] = "centre2";
    convection_op_words_[1] = "quick";
    convection_op_words_[2] = "discquick";
    convection_op_words_[3] = "quickinterface";
  }
  prefix_ = Nom("OpConv");
  suffix_ = Nom("IJKScalar_double");
}

Sortie& Operateur_IJK_elem_conv::printOn( Sortie& os ) const
{
  DERIV(OpConvIJKElemCommon_double)::printOn( os );
  return os;
}

Entree& Operateur_IJK_elem_conv::readOn( Entree& is )
{
  Cerr << "Read and Cast Operateur_IJK_elem_conv :" << finl;
  Motcle word;
  is >> word;
  Nom type = "";
  type += prefix_;
  convection_rank_ = convection_op_words_.search(word);
  switch(convection_rank_)
    {
    case 0 :
      {
        type += "Centre2";
        break;
      }
    case 1 :
      {
        type += "Quick";
        break;
      }
    case 2 :
      {
        type += "DiscQuick";
        break;
      }
    case 3 :
      {
        type += "QuickInterface";
        break;
      }
    default :
      {
        Cerr << "ERROR : Scalar convection operators that are already implemented are:" << finl;
        Cerr << convection_op_words_ << finl;
        abort();
      }
    }
  type += suffix_;
  typer(type);
  is >> valeur();
  return is;
}
