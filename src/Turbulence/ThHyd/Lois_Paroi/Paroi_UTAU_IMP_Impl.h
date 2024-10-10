/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
// File:        Paroi_UTAU_IMP_Impl.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Lois_Paroi
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Paroi_UTAU_IMP_Impl_included
#define Paroi_UTAU_IMP_Impl_included

#include <TRUSTTabs_forward.h>
#include <TRUST_Deriv.h>
#include <Parser_U.h>

class Champ_Don_base;

/*! @brief CLASS: Paroi_UTAU_IMP_Impl
 *
 * .SECTION  voir aussi
 *  Turbulence_paroi_base
 *
 */
class Paroi_UTAU_IMP_Impl
{
public:
  Entree& lire_donnees(Entree& s);
  double calculer_utau(const DoubleVect& pos, double norm_v, double d_visco);

protected:
  OWN_PTR(Champ_Don_base) u_star;
  OWN_PTR(Champ_Don_base) diam_hydr;
  Parser_U lambda_c;
  int u_star_ok;
  int lambda_c_ok;
};

#endif
