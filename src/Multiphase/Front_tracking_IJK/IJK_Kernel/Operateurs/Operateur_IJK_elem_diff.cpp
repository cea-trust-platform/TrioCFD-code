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
#include <Param.h>

Implemente_instanciable_sans_constructeur( Operateur_IJK_elem_diff, "Operateur_IJK_elem_diff", DERIV( Operateur_IJK_elem_diff_base_double ) );

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
  diffusion_op_ = "";
  diffusion_op_options_ = "";
  is_cast_ = false;
  diffusion_rank_ = 0;
}

void Operateur_IJK_elem_diff::reset_operator()
{
  diffusion_op_ = "";
  diffusion_op_options_ = "";
  is_cast_ = false;
  diffusion_rank_ = 0;
}

Sortie& Operateur_IJK_elem_diff::printOn(Sortie& os) const
{
  // DERIV(Operateur_IJK_elem_diff_base_double)::printOn( os );
  os << diffusion_op_words_[diffusion_rank_];
  return os;
}

Entree& Operateur_IJK_elem_diff::readOn(Entree& is)
{
  typer_diffusion_op(is);
//  Param param(que_suis_je());
//  set_param(param);
//  param.lire_sans_accolade(is);
  return is;
}

void Operateur_IJK_elem_diff::set_param(Param& param)
{
  param.ajouter_non_std("velocity_diffusion_form", (this));
}

int Operateur_IJK_elem_diff::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  return 1;
}

Entree& Operateur_IJK_elem_diff::typer_diffusion_op(Entree& is)
{
  if (is_cast_)
    reset_operator();
  Cerr << "Read and Cast Operateur_IJK_elem_diff :" << finl;
  Motcle word;
  is >> word;
  Nom type(get_diffusion_op_type(word));
  typer(type);
  Cerr << "Operateur_IJK_elem_diff cast to " << type << finl;
  is >> valeur();
  is_cast_ = true;
  return is;
}

void Operateur_IJK_elem_diff::typer_diffusion_op(const char *diffusion_op) // (const char * convection_op)
{
  if (is_cast_)
    reset_operator();
  Cerr << "Read and Cast Operateur_IJK_elem_diff :" << finl;
  Motcle word(diffusion_op);
  Nom type(get_diffusion_op_type(word));
  typer(type);
  Cerr << "Operateur_IJK_elem_diff cast to " << type << finl;
  is_cast_ = true;
}

Nom Operateur_IJK_elem_diff::get_diffusion_op_type(Motcle word)
{
  Nom type(prefix_);
  diffusion_rank_ = diffusion_op_words_.search(word);
  switch(diffusion_rank_)
    {
    case 0:
      {
        diffusion_op_ += "";
        break;
      }
    case 1:
      {
        diffusion_op_ += "Uniform";
        break;
      }
    case 2:
      {
        diffusion_op_ += "Anisotropic";
        break;
      }
    case 3:
      {
        diffusion_op_ += "Vectorial";
        break;
      }
    case 4:
      {
        diffusion_op_ += "VectorialAnisotropic";
        break;
      }
    case 5:
      {
        diffusion_op_ += "StructuralOnly";
        break;
      }
    default:
      {
        Cerr << "ERROR : Diffusion operators (elem) that are already implemented are:" << finl;
        Cerr << diffusion_op_words_ << finl;
        abort();
      }
    }
  type += diffusion_op_;
  type += suffix_;
  typer(type);
  return type;
}

void Operateur_IJK_elem_diff::set_conductivity_coefficient(const double& uniform_lambda, const IJK_Field_local_double& lambda, IJK_Field_local_double& coeff_field_x,
                                                           IJK_Field_local_double& coeff_field_y, IJK_Field_local_double& coeff_field_z)
{
  switch(diffusion_rank_)
    {
    case 0:
      {
        // Standard
        valeur().set_lambda(lambda);
        break;
      }
    case 1:
      {
        // Uniform
        valeur().set_uniform_lambda(uniform_lambda);
        break;
      }
    case 2:
      {
        // Anisotropic
        valeur().set_coeff_x_y_z(coeff_field_x, coeff_field_y, coeff_field_z);
        break;
      }
    case 3:
      {
        // Vectorial
        valeur().set_coeff_x_y_z(coeff_field_x, coeff_field_y, coeff_field_z);
        break;
      }
    case 4:
      {
        // VectorialAnisotropic
        valeur().set_coeff_x_y_z(coeff_field_x, coeff_field_y, coeff_field_z);
        break;
      }
    case 5:
      {
        // StructuralOnly
        valeur().set_coeff_x_y_z(coeff_field_x, coeff_field_y, coeff_field_z);
        break;
      }
    default:
      {
        Cerr << "ERROR : The conductivity can not be set properly" << finl;
        abort();
      }
    }
}
