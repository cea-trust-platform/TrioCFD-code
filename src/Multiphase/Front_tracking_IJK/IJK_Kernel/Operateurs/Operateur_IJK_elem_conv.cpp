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
#include <Param.h>

Implemente_instanciable_sans_constructeur( Operateur_IJK_elem_conv, "Operateur_IJK_elem_conv", DERIV(Operateur_IJK_elem_conv_base_double) );

Operateur_IJK_elem_conv::Operateur_IJK_elem_conv()
{
  convection_op_words_ = Motcles(5);
  {
    convection_op_words_[0] = "centre";
    convection_op_words_[1] = "centre2";
    convection_op_words_[2] = "quick";
    convection_op_words_[3] = "discquick";
    convection_op_words_[4] = "quickinterface";
  }
  prefix_ = Nom("OpConv");
  suffix_ = Nom("IJKScalar_double");
  convection_op_ = "";
  convection_op_option_ = "";
  convection_rank_ = 0;
  is_cast_ = false;
}

void Operateur_IJK_elem_conv::reset_operator()
{
  convection_op_ = "";
  convection_op_option_ = "";
  convection_rank_ = 0;
  is_cast_ = false;
}

Sortie& Operateur_IJK_elem_conv::printOn(Sortie& os) const
{
  // DERIV(Operateur_IJK_elem_conv_base_double)::printOn( os );
  os << convection_op_words_[convection_rank_];
  return os;
}

Entree& Operateur_IJK_elem_conv::readOn(Entree& is)
{
  typer_convection_op(is);
  return is;
}

void Operateur_IJK_elem_conv::set_param(Param& param)
{
  param.ajouter_non_std("velocity_convection_form", (this));
}

int Operateur_IJK_elem_conv::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  return 1;
}

Entree& Operateur_IJK_elem_conv::typer_convection_op(Entree& is)
{
  if (is_cast_)
    reset_operator();
  Cerr << "Read and Cast Operateur_IJK_elem_conv :" << finl;
  Motcle word;
  is >> word;
  Motcle type(get_convection_op_type(word));
  typer(type);
  Cerr << "Operateur_IJK_elem_conv cast to " << type << finl;
  is >> valeur();
  is_cast_ = true;
  return is;
}

void Operateur_IJK_elem_conv::typer_convection_op(const char *convection_op)
{
  if (is_cast_)
    reset_operator();
  Cerr << "Read and Cast Operateur_IJK_elem_conv :" << finl;
  Motcle word(convection_op);
  Motcle type(get_convection_op_type(word));
  typer(type);
  Cerr << "Operateur_IJK_elem_conv cast to " << type << finl;
  is_cast_ = true;
}

Nom Operateur_IJK_elem_conv::get_convection_op_type(Motcle word)
{
  Nom type(prefix_);
  convection_rank_ = convection_op_words_.search(word);
  switch(convection_rank_)
    {
    case 0:
      {
        convection_op_ += "Centre2";
        break;
      }
    case 1:
      {
        convection_op_ += "Centre2";
        break;
      }
    case 2:
      {
        convection_op_ += "Quick";
        break;
      }
    case 3:
      {
        convection_op_ += "DiscQuick";
        break;
      }
    case 4:
      {
        convection_op_ += "QuickInterface";
        break;
      }
    default:
      {
        Cerr << "ERROR : Scalar convection operators that are already implemented are:" << finl;
        Cerr << convection_op_words_ << finl;
        abort();
      }
    }
  type += convection_op_;
  type += suffix_;
  return type;
}
