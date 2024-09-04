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
// File      : Operateur_IJK_elem_diff.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////
#ifndef Operateur_IJK_elem_diff_included
#define Operateur_IJK_elem_diff_included

#include <TRUST_Deriv.h>
#include <Operateur_IJK_elem_diff_base.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Operateur_IJK_elem_diff
//
// <Description of class Operateur_IJK_elem_diff>
//
/////////////////////////////////////////////////////////////////////////////

class Operateur_IJK_elem_diff: public OWN_PTR( Operateur_IJK_elem_diff_base_double )
{

  Declare_instanciable( Operateur_IJK_elem_diff );

public:
  inline void initialize(const IJK_Splitting &splitting);

  void reset_operator();
  void typer_diffusion_op(const char *diffusion_op);
  Entree& typer_diffusion_op(Entree &is);
  int lire_motcle_non_standard(const Motcle &mot, Entree &is) override;
  void set_param(Param &param);
  Nom get_diffusion_op_type(Motcle word);

  void set_conductivity_coefficient(const double &uniform_lambda, const IJK_Field_local_double &lambda, IJK_Field_local_double &coeff_field_x, IJK_Field_local_double &coeff_field_y,
                                    IJK_Field_local_double &coeff_field_z);

protected:
  Motcles diffusion_op_words_;
  Nom prefix_;
  Nom suffix_;
  int diffusion_rank_;
  Nom diffusion_op_;
  Nom diffusion_op_options_;
  bool is_cast_;
};

inline void Operateur_IJK_elem_diff::initialize(const IJK_Splitting& splitting)
{
  if (!is_cast_)
    typer_diffusion_op("standard");
  valeur().initialize(splitting);
}

#endif /* Operateur_IJK_elem_diff_included */
