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
// File      : Deriv_Operateur_IJK_faces_conv.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Operateur_IJK_faces_conv_included
#define Operateur_IJK_faces_conv_included

#include <TRUST_Deriv.h>
#include <Operateur_IJK_faces_conv_base.h>
#include <OpConvAmontIJK.h>

class Operateur_IJK_faces_conv : public DERIV( Operateur_IJK_faces_conv_base_double )
{
  Declare_instanciable( Operateur_IJK_faces_conv ) ;

public:
  inline void initialize(const IJK_Splitting& splitting);

  int lire_motcle_non_standard(const Motcle& mot, Entree& is) override;
  void search_convection_option_type(Entree& is);
  void set_convection_option_type(Motcle word);
  void set_param(Param& param);
  void typer_convection_op(const char * convection_op);
  Entree& typer_convection_op(Entree& is);
  Nom get_convection_op_type( Motcle word );

  inline int get_convection_op_option_rank() { return convection_option_rank_; };
  inline Nom get_convection_op_option() { return convection_option_; };
  inline Nom get_convection_op() { return convection_op_; };

protected:
  Motcles convection_op_words_;
  Motcles convection_op_options_;
  Nom prefix_;
  Nom suffix_;
  Nom convection_op_;
  Nom convection_option_;
  int convection_rank_;
  int convection_option_rank_;
  bool is_cast_;
};

inline void Operateur_IJK_faces_conv::initialize(const IJK_Splitting& splitting)
{
  if (!is_cast_)
    {
      typer_convection_op("quick");
      // convection_option_ = convection_op_options_[0];
    }
  valeur().initialize(splitting);
}

#endif /* Operateur_IJK_faces_conv_included */
