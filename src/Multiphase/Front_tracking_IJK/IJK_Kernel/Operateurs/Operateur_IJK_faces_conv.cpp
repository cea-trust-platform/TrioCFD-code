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
// File      : Deriv_Operateur_IJK_faces_conv.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#include <Operateur_IJK_faces_conv.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur( Operateur_IJK_faces_conv, "Operateur_IJK_faces_conv", DERIV(Operateur_IJK_faces_conv_base_double) ) ;

Operateur_IJK_faces_conv::Operateur_IJK_faces_conv()
{
  convection_op_words_ = Motcles(9);
  {
    convection_op_words_[0] = "amont";
    convection_op_words_[1] = "Amont";
    convection_op_words_[2] = "centre";
    convection_op_words_[3] = "Centre";
    convection_op_words_[4] = "centre4";
    convection_op_words_[5] = "Centre4";
    convection_op_words_[6] = "quick";
    convection_op_words_[7] = "Quick";
    convection_op_words_[8] = "quicksharp";
  }
  /* non_conservative_simple : rho div(u u)
   * non_conservative_rhou     : div(rho u u) - u div(rho u)
   * conservative              : div(rho u u)
   */
  convection_op_ = "";
  convection_op_options_ = Motcles(3);
  {
    convection_op_options_[0] = "non_conservative_simple";
    convection_op_options_[1] = "non_conservative_rhou";
    convection_op_options_[2] = "conservative";
  }

  prefix_ = Nom("OpConv");
  suffix_ = Nom("IJKScalar_double");
}

Sortie& Operateur_IJK_faces_conv::printOn( Sortie& os ) const
{
  DERIV(Operateur_IJK_faces_conv_base_double)::printOn( os );
  return os;
}

Entree& Operateur_IJK_faces_conv::readOn( Entree& is )
{
  typer_convection_op(is);
  Param param(que_suis_je());
  set_param(param);
  param.lire_sans_accolade(is);
  return is;
}

void Operateur_IJK_faces_conv::set_param(Param& param)
{
  param.ajouter_non_std("velocity_convection_op", (this));
}

int Operateur_IJK_faces_conv::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="velocity_convection_op")
    {
      Motcle motlu;
      is >> motlu;
      if (Process::je_suis_maitre())
        Cerr << mot << " " << motlu << finl;
      int option_rank = convection_op_options_.search(motlu);
      if (option_rank < convection_op_options_.size() && option_rank !=-1)
        convection_option_ = convection_op_options_[option_rank];
      else
        convection_option_ = convection_op_options_[0];
    }
  return 1;
}

Entree& Operateur_IJK_faces_conv::typer_convection_op(Entree& is) // (const char * convection_op)
{
  Cerr << "Read and Cast Operateur_IJK_faces_conv :" << finl;
  Motcle convection_op;
  is >> convection_op;
  Motcle word(convection_op);
  Nom type(get_convection_op_type(word));
  typer(type);
  return is;
}

void Operateur_IJK_faces_conv::typer_convection_op(const char * convection_op) // (const char * convection_op)
{
  Cerr << "Read and Cast Operateur_IJK_faces_conv :" << finl;
  Motcle word(convection_op);
  Nom type(get_convection_op_type(word));
  typer(type);
}

Nom Operateur_IJK_faces_conv::get_convection_op_type( Motcle word ) // (const char * convection_op)
{
  Motcle convection_op;
  Motcle convection_key(convection_op);
  convection_rank_ = convection_op_words_.search(convection_key);
  Nom type(prefix_);
  switch(convection_rank_)
    {
    case 0 :
      {
        convection_op_ += "Centre4";
        break;
      }
    case 1 :
      {
        convection_op_ += "Quick";
        break;
      }
    case 2 :
      {
        convection_op_ += "DiscQuick";
        break;
      }
    case 3 :
      {
        convection_op_ += "Centre4";
        break;
      }
    case 4 :
      {
        convection_op_ += "Quick";
        break;
      }
    case 5 :
      {
        convection_op_ += "DiscQuick";
        break;
      }
    default :
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

