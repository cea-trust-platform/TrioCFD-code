/****************************************************************************
* Copyright (c) 2019, CEA
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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Op_Diff_K_Eps_Bas_Re_VDF_base.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Operateurs/Op_Dift
//
//////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_K_Eps_Bas_Re_VDF_base.h>

Implemente_base(Op_Diff_K_Eps_Bas_Re_VDF_base,"Op_Diff_K_Eps_Bas_Re_VDF_base",Op_Diff_K_Eps_Bas_Re_base);

Sortie& Op_Diff_K_Eps_Bas_Re_VDF_base::printOn(Sortie& s ) const { return s << que_suis_je() ; }
Entree& Op_Diff_K_Eps_Bas_Re_VDF_base::readOn(Entree& s ) { return s ; }

void Op_Diff_K_Eps_Bas_Re_VDF_base::completer()
{
  Operateur_base::completer();
  iter->completer_();
  iter->associer_champ_convecte_ou_inc(equation().inconnue(), nullptr);
  iter->set_name_champ_inco(equation().inconnue().le_nom().getString());
  iter->set_convective_op_pb_type(false /* diff op */, 0 /* pas pb_multiphase */);
}
