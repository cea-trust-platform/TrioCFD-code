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
// File:        Source_Transport_K_Eps_Bas_Reynolds_VEF_Face.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Source_Transport_K_Eps_Bas_Reynolds_VEF_Face_included
#define Source_Transport_K_Eps_Bas_Reynolds_VEF_Face_included

#include <Source_Transport_VEF_Face_base.h>
#include <TRUST_Ref.h>

class Transport_K_Eps_Bas_Reynolds;

class Source_Transport_K_Eps_Bas_Reynolds_VEF_Face : public Source_Transport_VEF_Face_base
{
  Declare_instanciable_sans_constructeur(Source_Transport_K_Eps_Bas_Reynolds_VEF_Face);
public :
  Source_Transport_K_Eps_Bas_Reynolds_VEF_Face(double cte1 = C11__, double cte2 = C21__) : Source_Transport_VEF_Face_base(cte1,cte2) { }
  DoubleTab& ajouter(DoubleTab& ) const override;
  void associer_pb(const Probleme_base& ) override;

protected :
  REF(Transport_K_Eps_Bas_Reynolds) eqn_keps_bas_re;
};

#endif /* Source_Transport_K_Eps_Bas_Reynolds_VEF_Face_included */