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
// File      : Deriv_Operateur_IJK_faces_diff.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////
#ifndef Operateur_IJK_faces_diff_included
#define Operateur_IJK_faces_diff_included

#include <TRUST_Deriv.h>
#include <Operateur_IJK_faces_diff_base.h>

class Operateur_IJK_faces_diff: public DERIV( Operateur_IJK_faces_diff_base_double )
{
  Declare_instanciable( Operateur_IJK_faces_diff );
public:
  inline void initialize(const IJK_Splitting &splitting, const int harmonic_nu);
  Entree& typer_diffusion_op(Entree &is);
  void typer_diffusion_op(const char *diffusion_op);
  int lire_motcle_non_standard(const Motcle &mot, Entree &is) override;
  void set_param(Param &param);
  Nom get_diffusion_op_type(Motcle word);
  int get_diffusion_op_option_rank() { return diffusion_option_rank_; }
  Nom get_diffusion_op_option() { return diffusion_option_; }
  Nom get_diffusion_op() { return diffusion_op_; }
  int harmonic_nu_;

protected:
  Motcles diffusion_op_words_;
  Motcles diffusion_op_options_;
  Nom prefix_;
  Nom suffix_;
  int diffusion_rank_;
  Nom diffusion_op_;
  Nom diffusion_option_;
  int diffusion_option_rank_;
  bool is_cast_;
};

inline void Operateur_IJK_faces_diff::initialize(const IJK_Splitting& splitting, const int harmonic_nu = 0)
{
  if (!is_cast_)
    typer_diffusion_op("standard");
  // diffusion_option_ = diffusion_op_options_[0];
  valeur().initialize(splitting);
  harmonic_nu_ = harmonic_nu;
}

#endif /* Operateur_IJK_faces_diff_included */
