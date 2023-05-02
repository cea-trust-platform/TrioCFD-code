/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Source_Con_Phase_field_base.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Multiphase/Phase_field/src
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Source_Con_Phase_field_base_included
#define Source_Con_Phase_field_base_included


#include <Source_base.h>
#include <TRUSTTab.h>
#include <Champ_Don.h>
#include <TRUST_Ref.h>

class Probleme_base;

/*! @brief class Source_Con_Phase_field_base
 *
 *
 *
 * @sa Source_base
 */
class Source_Con_Phase_field_base : public Source_base
{

  Declare_base(Source_Con_Phase_field_base);

public:
  virtual void premier_demi_dt()=0;
  virtual const DoubleVect& get_u_carre()=0;
  virtual void calculer_champ_fonc_c(const double t, Champ_Don& champ_fonc_c, const DoubleTab& val_c) const = 0;
};

#endif
