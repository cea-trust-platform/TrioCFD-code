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
// File:        MuLambda_TBLE_base.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Lois_Paroi/TBLE/MuLambda
//
//////////////////////////////////////////////////////////////////////////////

#ifndef MuLambda_TBLE_base_included
#define MuLambda_TBLE_base_included

#include <Objet_U.h>
#include <TRUST_Ref.h>

class Eq_couch_lim;

class Milieu_base;
/*! @brief Classe MuLambda_TBLE_base Classe abstraite calculant Mu et Lambda sur le maillage TBLE
 *
 *
 *
 */
class MuLambda_TBLE_base : public Objet_U
{
  Declare_base(MuLambda_TBLE_base);

public :

  virtual void initialiser(const Milieu_base&) = 0 ;

  // "ind" represente le numero d'une maille TBLE sur l'instance eq_couch_lim referencee dans cette classe
  virtual double getNu(OBS_PTR(Eq_couch_lim) eq_T, int ind) = 0 ;
  virtual double getAlpha(OBS_PTR(Eq_couch_lim) eq_T, int ind) = 0 ;


protected :

};


#endif
