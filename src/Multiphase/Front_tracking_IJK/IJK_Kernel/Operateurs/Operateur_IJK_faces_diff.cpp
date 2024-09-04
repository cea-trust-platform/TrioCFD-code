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
// File      : Deriv_Operateur_IJK_faces_diff.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#include <Operateur_IJK_faces_diff.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur( Operateur_IJK_faces_diff, "Operateur_IJK_faces_diff", OWN_PTR( Operateur_IJK_faces_diff_base_double ) );

Operateur_IJK_faces_diff::Operateur_IJK_faces_diff()
{
  /* simple_arithmetic : div(mu grad(u))
   * full_arithmetic   : Tenseur des contraintes complet : div[mu (grad(u)+grad^T(u))]
   *                     mu : moyenne arithmetique
   * full_adaptative   : Tenseur des contraintes complet : div[mu (grad(u)+grad^T(u))]
   * 										 mu : switch from arithmetic to geometric mean depending on the direction (Not available yet)
   */
  diffusion_op_words_ = Motcles(17);
  {
    diffusion_op_words_[0] = "standard"; //
    diffusion_op_words_[1] = "turb"; //
    diffusion_op_words_[2] = "anisotropic"; //
    diffusion_op_words_[3] = "laminar_transpose"; //
    diffusion_op_words_[4] = "tensorial_zero_wall";
    diffusion_op_words_[5] = "tensorial_anisotropic_zero_wall";
    diffusion_op_words_[6] = "laminar_transpose_anisotropic"; //
    diffusion_op_words_[7] = "laminar_transpose_tensorial_zero_wall"; //
    diffusion_op_words_[8] = "laminar_transpose_tensorial_anisotropic_zero_wall"; //
    diffusion_op_words_[9] = "laminar_transpose_divergence"; //
    diffusion_op_words_[10] = "laminar_transpose_divergence_anisotropic"; //
    diffusion_op_words_[11] = "laminar_transpose_divergence_tensorial_anisotropic_zero_wall"; //
    diffusion_op_words_[12] = "laminar_transpose_divergence_tensorial_zero_wall"; //
    diffusion_op_words_[13] = "structural_zero_wall"; //
    diffusion_op_words_[14] = "simple_arithmetic";
    diffusion_op_words_[15] = "full_arithmetic";
    diffusion_op_words_[16] = "full_adaptative";
  }
  diffusion_op_ = "";
  diffusion_op_options_ = Motcles(3);
  {
    diffusion_op_options_[0] = "simple_arithmetic";
    diffusion_op_options_[1] = "full_arithmetic";
    diffusion_op_options_[2] = "full_adaptative";
  }
  diffusion_option_ = "";
  prefix_ = Nom("OpDiff");
  suffix_ = Nom("IJK_double");
  is_cast_ = false;
  diffusion_rank_ = 0;
  diffusion_option_rank_ = 0;
}

Sortie& Operateur_IJK_faces_diff::printOn(Sortie& os) const
{
  // OWN_PTR(Operateur_IJK_faces_diff_base_double)::printOn( os );
  os << diffusion_op_words_[diffusion_rank_] << "_" << diffusion_option_ << "\n";
  return os;
}

Entree& Operateur_IJK_faces_diff::readOn(Entree& is)
{

  typer_diffusion_op(is);
  return is;
}

void Operateur_IJK_faces_diff::set_param(Param& param)
{
  param.ajouter_non_std("velocity_diffusion_form", (this));
}

int Operateur_IJK_faces_diff::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "velocity_diffusion_form")
    {
      Motcle motlu;
      is >> motlu;
      if (Process::je_suis_maitre())
        Cerr << mot << " " << motlu << finl;
      int option_rank = diffusion_op_options_.search(motlu);
      if (option_rank < diffusion_op_options_.size() && option_rank != -1)
        diffusion_option_ = diffusion_op_options_[option_rank];
      else
        diffusion_option_ = diffusion_op_options_[0];
    }
  return 1;
}

Entree& Operateur_IJK_faces_diff::typer_diffusion_op(Entree& is)
{
  Cerr << "Read and Cast Operateur_IJK_faces_diff :" << finl;
  Motcle word;
  is >> word;
  Nom type(get_diffusion_op_type(word));
  typer(type);
  is >> valeur();
  is_cast_ = true;
  return is;
}

void Operateur_IJK_faces_diff::typer_diffusion_op(const char *diffusion_op)
{
  Cerr << "Read and Cast Operateur_IJK_faces_diff :" << finl;
  Motcle word(diffusion_op);
  Motcle type(get_diffusion_op_type(word));
  typer(type);
  is_cast_ = true;
}

Nom Operateur_IJK_faces_diff::get_diffusion_op_type(Motcle word)
{
  diffusion_option_rank_ = diffusion_op_options_.search(word);
  if (diffusion_option_rank_ == -1)
    diffusion_option_rank_ = 0;
  diffusion_option_ = diffusion_op_options_[diffusion_option_rank_];

  diffusion_rank_ = diffusion_op_words_.search(word);
  Nom type(prefix_);
  // TODO: Use enum instead such as in IJK_FT:   enum TimeScheme { EULER_EXPLICITE, RK3_FT }; ??
  switch(diffusion_rank_)
    {
    case 0:
      {
        break;
      }
    case 1:
      {
        diffusion_op_ += "Turb";
        break;
      }
    case 2:
      {
        diffusion_op_ += "Anisotropic";
        break;
      }
    case 3:
      {
        diffusion_op_ += "StdWithLaminarTranspose";
        break;
      }
    case 4:
      {
        diffusion_op_ += "TensorialZeroatwall";
        break;
      }
    case 5:
      {
        diffusion_op_ += "TensorialAnisotropicZeroatwall";
        break;
      }
    case 6:
      {
        diffusion_op_ += "StdWithLaminarTransposeAnisotropic";
        break;
      }
    case 7:
      {
        diffusion_op_ += "StdWithLaminarTransposeTensorialZeroatwall";
        break;
      }
    case 8:
      {
        diffusion_op_ += "StdWithLaminarTransposeTensorialAnisotropicZeroatwall";
        break;
      }
    case 9:
      {
        diffusion_op_ += "StdWithLaminarTransposeAndDivergence";
        break;
      }
    case 10:
      {
        diffusion_op_ += "StdWithLaminarTransposeAndDivergenceAnisotropic";
        break;
      }
    case 11:
      {
        diffusion_op_ += "StdWithLaminarTransposeAndDivergenceTensorialAnisotropicZeroatwall";
        break;
      }
    case 12:
      {
        diffusion_op_ += "StdWithLaminarTransposeAndDivergenceTensorialZeroatwall";
        break;
      }
    case 13:
      {
        diffusion_op_ += "StructuralOnlyZeroatwall";
        break;
      }
    case 14:
      {
        break;
      }
    case 15:
      {
        diffusion_op_ += "StdWithLaminarTranspose";
        break;
      }
    case 16:
      {
        // TODO: Full adaptive (viscosity with direction dependancy !)
        Cerr << "Unknown velocity diffusion operator! " << finl;
        abort();
        //        Process::exit();
        break;
      }
    default:
      {
        Cerr << "ERROR : Diffusion operators (faces) that are already implemented are:" << finl;
        Cerr << diffusion_op_words_ << finl;
        abort();
      }
    }
  type += diffusion_op_;
  type += suffix_;
  return type;
}
